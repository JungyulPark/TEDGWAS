# =============================================================================
# TED-TRAP Step 3b (PATCHED) — FinnGen Coloc at TSHR locus
# Fix: Comprehensive NA filtering + explicit SNP ordering before coloc.abf
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(coloc)
    library(ieugwasr)
})

token <- Sys.getenv("OPENGWAS_JWT")
# options(ieugwasr_api = "https://api.opengwas.io/") # Disabled for stability

for (d in c("TrackA_MR/results", "TrackA_MR/logs")) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

sink("TrackA_MR/logs/03b_finngen_coloc_patched.log", split = TRUE)
cat("=== Step 3b (Patched): FinnGen R12 Coloc at TSHR locus ===\n\n")

# --- Config ---
FG_FILE <- "c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz"
FG_N <- 500348L
FG_CASES <- 858L
TSHR_CHR <- 14
TSHR_START <- 80921085L
TSHR_END <- 82081779L

# =============================================================================
# 1. Extract FinnGen TSHR locus
# =============================================================================

cat("--- Reading FinnGen GRAVES_OPHT, TSHR locus ---\n")
fg <- fread(FG_FILE,
    select = c(
        "#chrom", "pos", "ref", "alt", "rsids",
        "beta", "sebeta", "pval", "af_alt"
    )
)
setnames(fg, "#chrom", "chrom")

# Filter locus
fg_loc <- fg[chrom == TSHR_CHR & pos >= TSHR_START & pos <= TSHR_END]
rm(fg)
gc(verbose = FALSE)

cat(sprintf("FinnGen TSHR locus: %d SNPs (raw)\n", nrow(fg_loc)))

# Clean rsid
fg_loc[, rsid := sub(",.*", "", rsids)]

# STRICT cleaning — remove any row with NA in ANY required field
fg_loc <- fg_loc[
    !is.na(rsid) & grepl("^rs", rsid) &
        !is.na(beta) & !is.na(sebeta) & sebeta > 0 &
        !is.na(pos) & !is.na(af_alt) &
        af_alt > 0 & af_alt < 1
]
# Dedupe
fg_loc <- fg_loc[!duplicated(rsid)]
cat(sprintf("FinnGen after strict cleaning: %d SNPs\n", nrow(fg_loc)))

# =============================================================================
# 2. Fetch eQTL TSHR locus from OpenGWAS
# =============================================================================

cat("\n--- Fetching eQTL TSHR from OpenGWAS ---\n")
eqtl <- suppressMessages(ieugwasr::associations(
    variants = sprintf("%d:%d-%d", TSHR_CHR, TSHR_START, TSHR_END),
    id = "eqtl-a-ENSG00000165409",
    opengwas_jwt = token
))
cat(sprintf("eQTL raw: %d SNPs\n", nrow(eqtl)))

# STRICT cleaning — match FinnGen rigor
eqtl <- as.data.table(eqtl)
eqtl <- eqtl[
    !is.na(rsid) & grepl("^rs", rsid) &
        !is.na(beta) & !is.na(se) & se > 0 &
        !is.na(position) & !is.na(eaf) &
        eaf > 0 & eaf < 1
]
eqtl <- eqtl[!duplicated(rsid)]
cat(sprintf("eQTL after strict cleaning: %d SNPs\n", nrow(eqtl)))

# =============================================================================
# 3. Find common SNPs + ORDER THEM IDENTICALLY
# =============================================================================

common_snps <- intersect(eqtl$rsid, fg_loc$rsid)
cat(sprintf("\nCommon SNPs: %d\n", length(common_snps)))

if (length(common_snps) < 100) {
    sink()
    stop(sprintf("Too few common SNPs: %d", length(common_snps)))
}

# CRITICAL FIX: explicit ordering
common_snps <- sort(common_snps)
eqtl_m <- eqtl[match(common_snps, eqtl$rsid)]
fg_m <- fg_loc[match(common_snps, fg_loc$rsid)]

# Verify no NA after matching
stopifnot(all(!is.na(eqtl_m$beta)), all(!is.na(fg_m$beta)))
stopifnot(nrow(eqtl_m) == nrow(fg_m))
stopifnot(all(eqtl_m$rsid == fg_m$rsid))

cat(sprintf("After ordering + verification: %d SNPs in both datasets\n", nrow(eqtl_m)))

# =============================================================================
# 4. Build coloc inputs — defensive MAF + no NAs
# =============================================================================

# eQTL MAF
maf_eqtl <- pmin(eqtl_m$eaf, 1 - eqtl_m$eaf)
maf_eqtl[maf_eqtl <= 0 | maf_eqtl >= 0.5 | is.na(maf_eqtl)] <- 0.3

# FinnGen MAF
maf_fg <- pmin(fg_m$af_alt, 1 - fg_m$af_alt)
maf_fg[maf_fg <= 0 | maf_fg >= 0.5 | is.na(maf_fg)] <- 0.3

eqtl_input <- list(
    snp      = unname(as.character(common_snps)),
    beta     = as.numeric(eqtl_m$beta),
    varbeta  = as.numeric(eqtl_m$se)^2,
    MAF      = as.numeric(maf_eqtl),
    position = as.integer(eqtl_m$position),
    N        = 31684L,
    type     = "quant",
    sdY      = 1
)

fg_input <- list(
    snp      = unname(as.character(common_snps)),
    beta     = as.numeric(fg_m$beta),
    varbeta  = as.numeric(fg_m$sebeta)^2,
    MAF      = as.numeric(maf_fg),
    position = as.integer(fg_m$pos),
    N        = FG_N,
    type     = "cc",
    s        = FG_CASES / FG_N
)

# Validate before coloc
cat("\n--- Validating inputs ---\n")
cat(sprintf(
    "  eQTL: %d SNPs, type=quant, any NA: %s\n",
    length(eqtl_input$snp),
    any(is.na(c(eqtl_input$beta, eqtl_input$varbeta, eqtl_input$MAF)))
))
cat(sprintf(
    "  FG:   %d SNPs, type=cc, s=%.5f, any NA: %s\n",
    length(fg_input$snp), fg_input$s,
    any(is.na(c(fg_input$beta, fg_input$varbeta, fg_input$MAF)))
))

coloc::check_dataset(eqtl_input)
coloc::check_dataset(fg_input, req = "s")

# =============================================================================
# 5. Run coloc
# =============================================================================

cat("\n--- Running coloc.abf ---\n")
result <- coloc.abf(dataset1 = eqtl_input, dataset2 = fg_input)

pp <- result$summary
cat(sprintf("\n📊 RESULTS (eQTLGen TSHR × FinnGen R12 GRAVES_OPHT):\n"))
cat(sprintf("  PP.H0 = %.4f\n", pp["PP.H0.abf"]))
cat(sprintf("  PP.H1 = %.4f\n", pp["PP.H1.abf"]))
cat(sprintf("  PP.H2 = %.4f\n", pp["PP.H2.abf"]))
cat(sprintf("  PP.H3 = %.4f\n", pp["PP.H3.abf"]))
cat(sprintf(
    "  PP.H4 = %.4f  %s\n", pp["PP.H4.abf"],
    ifelse(pp["PP.H4.abf"] > 0.95, "⭐⭐⭐ Very Strong",
        ifelse(pp["PP.H4.abf"] > 0.80, "⭐⭐ Strong",
            ifelse(pp["PP.H4.abf"] > 0.50, "⭐ Suggestive", "❌ None")
        )
    )
))

top_idx <- which.max(result$results$SNP.PP.H4)
cat(sprintf(
    "\n  Top shared SNP: %s (SNP-PP.H4 = %.3f)\n",
    result$results$snp[top_idx],
    result$results$SNP.PP.H4[top_idx]
))

# Save
out_tab <- data.frame(
    Outcome = "FinnGen_R12_GRAVES_OPHT",
    N_SNPs = length(common_snps),
    PP_H0 = round(pp["PP.H0.abf"], 4),
    PP_H1 = round(pp["PP.H1.abf"], 4),
    PP_H2 = round(pp["PP.H2.abf"], 4),
    PP_H3 = round(pp["PP.H3.abf"], 4),
    PP_H4 = round(pp["PP.H4.abf"], 4),
    Top_SNP = result$results$snp[top_idx],
    Top_SNP_PP = round(result$results$SNP.PP.H4[top_idx], 3),
    Interpretation = ifelse(pp["PP.H4.abf"] > 0.95, "Very Strong",
        ifelse(pp["PP.H4.abf"] > 0.80, "Strong",
            ifelse(pp["PP.H4.abf"] > 0.50, "Suggestive", "None")
        )
    ),
    stringsAsFactors = FALSE
)
write.csv(out_tab, "TrackA_MR/results/03b_coloc_FinnGen_GRAVES_OPHT.csv", row.names = FALSE)
saveRDS(result, "TrackA_MR/results/03b_coloc_FinnGen_full.rds")

cat("\n💾 Saved: 03b_coloc_FinnGen_GRAVES_OPHT.csv\n")
sink()
cat("\n✅ Step 3b patched & complete.\n")
