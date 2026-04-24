# =============================================================================
# TED-TRAP Step 3b + Step 4 (v3) — FinnGen R12 local analysis
#
# Step 3b: Colocalization at TSHR locus using FinnGen R12 GRAVES_OPHT (TED)
#           → Sensitivity outcome for Step 3 (OpenGWAS FinnGen failed earlier)
#
# Step 4 (v3): MVMR with valid TED outcome
#           Exposure 1: TSHR cis-eQTL (eQTLGen via OpenGWAS)
#           Exposure 2: GD liability — FinnGen R12 E4_GRAVES_STRICT (local)
#           Outcome:    TED — FinnGen R12 GRAVES_OPHT (local)
# =============================================================================

# --- Packages ---
required_pkgs <- c("data.table", "coloc", "ieugwasr", "dplyr")
for (p in required_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "http://cran.us.r-project.org")
}
suppressPackageStartupMessages({
    library(data.table)
    library(coloc)
    library(ieugwasr)
    library(dplyr)
})

token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(token) < 100) stop("OPENGWAS_JWT not set.")

# --- Output dirs ---
for (d in c("TrackA_MR/results", "TrackA_MR/logs")) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

log_file <- "TrackA_MR/logs/03b_04v3_finngen_local.log"
sink(log_file, split = TRUE)
cat("=== Step 3b + 4v3: FinnGen R12 local analysis ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# =============================================================================
# GLOBAL CONFIG
# =============================================================================

FG_DIR <- "c:/ProjectTEDGWAS" # <--- FIXED DIRECTORY

FG_FILES <- list(
    TED_main       = "finngen_R12_GRAVES_OPHT.gz",
    GD_exposure    = "finngen_R12_E4_GRAVES_STRICT.gz",
    TED_strict     = "finngen_R12_E4_GRAVES_OPHT_STRICT.gz"
)

FG_META <- list(
    TED_main    = list(N = 500348, cases = 858, controls = 499490, label = "GRAVES_OPHT"),
    GD_exposure = list(N = 500348, cases = 3962, controls = 496386, label = "E4_GRAVES_STRICT"),
    TED_strict  = list(N = 500348, cases = 753, controls = 499595, label = "E4_GRAVES_OPHT_STRICT")
)

TSHR_LOCUS <- list(chr = 14, start = 80921085, end = 82081779) # GRCh37/hg19

cat("Configuration:\n")
cat(sprintf("  FG_DIR: %s\n", FG_DIR))
for (nm in names(FG_FILES)) {
    cat(sprintf("  %-12s: %s (cases=%d)\n", nm, FG_FILES[[nm]], FG_META[[nm]]$cases))
}
cat(sprintf(
    "  TSHR locus (GRCh37): chr%d:%d-%d\n\n",
    TSHR_LOCUS$chr, TSHR_LOCUS$start, TSHR_LOCUS$end
))

# =============================================================================
# HELPER 1 — Fast locus extraction from FinnGen .gz file
# =============================================================================

extract_finngen_locus <- function(file_path, chr, start, end, label) {
    cat(sprintf("\n  Extracting %s at chr%d:%d-%d\n", label, chr, start, end))

    dt <- fread(file_path,
        select = c(
            "#chrom", "pos", "ref", "alt", "rsids",
            "beta", "sebeta", "pval", "af_alt"
        )
    )
    setnames(dt, "#chrom", "chrom")

    dt_locus <- dt[chrom == chr & pos >= start & pos <= end]
    cat(sprintf("    Total SNPs in file: %d\n", nrow(dt)))
    cat(sprintf("    SNPs in locus: %d\n", nrow(dt_locus)))

    dt_locus[, rsid := sub(",.*", "", rsids)]
    dt_locus <- dt_locus[!is.na(rsid) & rsid != "" & grepl("^rs", rsid)]

    dt_locus <- dt_locus[!is.na(beta) & !is.na(sebeta) & sebeta > 0]
    dt_locus <- dt_locus[!duplicated(rsid)]

    cat(sprintf("    After cleaning: %d SNPs\n", nrow(dt_locus)))

    rm(dt)
    gc(verbose = FALSE)
    dt_locus
}

# =============================================================================
# HELPER 2 — Fast genome-wide SNP lookup (for MVMR)
# =============================================================================

extract_finngen_snps <- function(file_path, snp_list, label) {
    cat(sprintf("\n  Looking up %d SNPs in %s\n", length(snp_list), label))

    dt <- fread(file_path,
        select = c(
            "#chrom", "pos", "ref", "alt", "rsids",
            "beta", "sebeta", "pval", "af_alt"
        )
    )
    setnames(dt, "#chrom", "chrom")
    dt[, rsid := sub(",.*", "", rsids)]

    hits <- dt[rsid %in% snp_list]
    hits <- hits[!is.na(beta) & !is.na(sebeta) & sebeta > 0]
    hits <- hits[!duplicated(rsid)]

    cat(sprintf(
        "    Matched: %d / %d SNPs (%.1f%%)\n",
        nrow(hits), length(snp_list),
        100 * nrow(hits) / length(snp_list)
    ))

    rm(dt)
    gc(verbose = FALSE)
    hits
}

# ==============================================================================
# PART 1 — STEP 3b: Colocalization at TSHR locus
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PART 1: Colocalization at TSHR locus (FinnGen R12 GRAVES_OPHT)\n")
cat(rep("=", 70), "\n", sep = "")

fg_ted_locus <- extract_finngen_locus(
    file.path(FG_DIR, FG_FILES$TED_main),
    TSHR_LOCUS$chr, TSHR_LOCUS$start, TSHR_LOCUS$end,
    FG_META$TED_main$label
)

cat("\n  Fetching eQTL TSHR from OpenGWAS...\n")
eqtl_locus <- tryCatch(
    ieugwasr::associations(
        variants = sprintf("%d:%d-%d", TSHR_LOCUS$chr, TSHR_LOCUS$start, TSHR_LOCUS$end),
        id = "eqtl-a-ENSG00000165409",
        opengwas_jwt = token
    ),
    error = function(e) {
        cat(sprintf("    [ERROR] %s\n", e$message))
        NULL
    }
)

if (is.null(eqtl_locus) || nrow(eqtl_locus) == 0) {
    cat("  ⚠️ eQTL fetch failed — skipping coloc\n")
    coloc_result <- NULL
} else {
    cat(sprintf("    eQTL SNPs: %d\n", nrow(eqtl_locus)))

    eqtl_clean <- eqtl_locus[!is.na(eqtl_locus$beta) & !is.na(eqtl_locus$se) & eqtl_locus$se > 0, ]
    eqtl_clean <- eqtl_clean[!duplicated(eqtl_clean$rsid), ]

    common_snps <- intersect(eqtl_clean$rsid, fg_ted_locus$rsid)
    cat(sprintf("\n  Common SNPs: %d\n", length(common_snps)))

    if (length(common_snps) < 50) {
        cat("  ⚠️ Too few common SNPs; skipping coloc\n")
        coloc_result <- NULL
    } else {
        # CRITICAL FIX: Cast to character arrays and proper data matching parameters
        eqtl_sub <- eqtl_clean[match(common_snps, eqtl_clean$rsid), ]
        fg_sub <- fg_ted_locus[match(common_snps, fg_ted_locus$rsid)]

        eqtl_input <- list(
            snp      = as.character(eqtl_sub$rsid),
            beta     = as.numeric(eqtl_sub$beta),
            varbeta  = as.numeric(eqtl_sub$se)^2,
            MAF      = pmin(as.numeric(eqtl_sub$eaf), 1 - as.numeric(eqtl_sub$eaf)),
            position = as.integer(eqtl_sub$position),
            N        = 31684L,
            type     = "quant"
        )
        eqtl_input$MAF[is.na(eqtl_input$MAF) | eqtl_input$MAF <= 0 | eqtl_input$MAF >= 0.5] <- 0.3

        maf_fg <- pmin(as.numeric(fg_sub$af_alt), 1 - as.numeric(fg_sub$af_alt))
        maf_fg[is.na(maf_fg) | maf_fg <= 0 | maf_fg >= 0.5] <- 0.3

        fg_input <- list(
            snp      = as.character(fg_sub$rsid),
            beta     = as.numeric(fg_sub$beta),
            varbeta  = as.numeric(fg_sub$sebeta)^2,
            MAF      = as.numeric(maf_fg),
            position = as.integer(fg_sub$pos),
            N        = as.integer(FG_META$TED_main$N),
            type     = "cc",
            s        = FG_META$TED_main$cases / FG_META$TED_main$N
        )

        cat(sprintf(
            "  FinnGen s = %.5f (valid: %s)\n",
            fg_input$s, fg_input$s > 0 && fg_input$s < 1
        ))

        cat("\n  Running coloc.abf...\n")
        coloc_result <- tryCatch(
            coloc.abf(dataset1 = eqtl_input, dataset2 = fg_input),
            error = function(e) {
                cat(sprintf("  [ERROR] %s\n", e$message))
                NULL
            }
        )

        if (!is.null(coloc_result)) {
            pp <- coloc_result$summary
            cat(sprintf("\n  📊 RESULTS:\n"))
            cat(sprintf("    PP.H0 = %.4f\n", pp["PP.H0.abf"]))
            cat(sprintf("    PP.H1 = %.4f\n", pp["PP.H1.abf"]))
            cat(sprintf("    PP.H2 = %.4f\n", pp["PP.H2.abf"]))
            cat(sprintf("    PP.H3 = %.4f\n", pp["PP.H3.abf"]))
            cat(sprintf(
                "    PP.H4 = %.4f %s\n", pp["PP.H4.abf"],
                ifelse(pp["PP.H4.abf"] > 0.95, "⭐⭐⭐ Very Strong",
                    ifelse(pp["PP.H4.abf"] > 0.80, "⭐⭐ Strong",
                        ifelse(pp["PP.H4.abf"] > 0.50, "⭐ Suggestive", "❌ None")
                    )
                )
            ))

            top_idx <- which.max(coloc_result$results$SNP.PP.H4)
            cat(sprintf(
                "    Top shared SNP: %s (PP.H4 = %.3f)\n",
                coloc_result$results$snp[top_idx],
                coloc_result$results$SNP.PP.H4[top_idx]
            ))

            coloc_finngen_table <- data.frame(
                Outcome = "FinnGen_R12_GRAVES_OPHT",
                N_SNPs = length(common_snps),
                PP_H0 = round(pp["PP.H0.abf"], 4),
                PP_H1 = round(pp["PP.H1.abf"], 4),
                PP_H2 = round(pp["PP.H2.abf"], 4),
                PP_H3 = round(pp["PP.H3.abf"], 4),
                PP_H4 = round(pp["PP.H4.abf"], 4),
                Top_SNP = coloc_result$results$snp[top_idx],
                stringsAsFactors = FALSE
            )
            write.csv(coloc_finngen_table, "TrackA_MR/results/03b_coloc_FinnGen_GRAVES_OPHT.csv", row.names = FALSE)
            saveRDS(coloc_result, "TrackA_MR/results/03b_coloc_FinnGen_full.rds")
            cat("  💾 Saved: 03b_coloc_FinnGen_GRAVES_OPHT.csv\n")
        }
    }
}

# ==============================================================================
# PART 2 — STEP 4 v3: MVMR with FinnGen R12 outcomes
# ==============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("PART 2: MVMR with FinnGen R12 TED outcome\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("=== 4.1 Extracting instruments ===\n")
cat("\nExposure 1: TSHR cis-eQTL (OpenGWAS)\n")
exp1_tshr <- tryCatch(
    ieugwasr::tophits(
        id = "eqtl-a-ENSG00000165409", pval = 5e-8, clump = 1,
        r2 = 0.001, kb = 10000, opengwas_jwt = token
    ),
    error = function(e) NULL
)
if (is.null(exp1_tshr) || nrow(exp1_tshr) < 1) {
    exp1_tshr <- ieugwasr::tophits(
        id = "eqtl-a-ENSG00000165409", pval = 5e-6, clump = 1,
        r2 = 0.001, kb = 10000, opengwas_jwt = token
    )
}
cat(sprintf("  TSHR IVs: %d\n", nrow(exp1_tshr)))

cat("\nExposure 2: GD liability (FinnGen R12 E4_GRAVES_STRICT local)\n")
cat("  Reading full file and extracting P<5e-8 SNPs...\n")

gd_full <- fread(file.path(FG_DIR, FG_FILES$GD_exposure),
    select = c(
        "#chrom", "pos", "ref", "alt", "rsids",
        "beta", "sebeta", "pval", "af_alt"
    )
)
setnames(gd_full, "#chrom", "chrom")
gd_full[, rsid := sub(",.*", "", rsids)]
gd_full <- gd_full[!is.na(beta) & !is.na(sebeta) & sebeta > 0 & grepl("^rs", rsid)]

gd_sig <- gd_full[pval < 5e-8]
cat(sprintf("  GD significant SNPs (P<5e-8): %d\n", nrow(gd_sig)))

cat("  Clumping (r²<0.001, 10Mb)...\n")
clump_input <- data.frame(
    rsid = gd_sig$rsid,
    pval = gd_sig$pval,
    id = "finngen_GD"
)

gd_clumped <- tryCatch(
    ieugwasr::ld_clump(
        dat = clump_input, clump_r2 = 0.001, clump_kb = 10000,
        pop = "EUR", opengwas_jwt = token
    ),
    error = function(e) {
        cat(sprintf("  [CLUMP ERROR] %s\n", e$message))
        NULL
    }
)

if (is.null(gd_clumped)) {
    cat("  Falling back to simple position-based pruning (>1Mb apart)\n")
    gd_sig <- gd_sig[order(pval)]
    kept <- integer(0)
    for (i in seq_len(nrow(gd_sig))) {
        if (length(kept) == 0) {
            kept <- c(kept, i)
            next
        }
        current_chr <- gd_sig$chrom[i]
        current_pos <- gd_sig$pos[i]
        collisions <- sapply(kept, function(k) {
            gd_sig$chrom[k] == current_chr && abs(gd_sig$pos[k] - current_pos) < 1e6
        })
        if (!any(collisions)) kept <- c(kept, i)
    }
    gd_clumped_rsids <- gd_sig$rsid[kept]
    cat(sprintf("  Position-pruned: %d SNPs\n", length(gd_clumped_rsids)))
} else {
    gd_clumped_rsids <- gd_clumped$rsid
    cat(sprintf("  LD-clumped: %d SNPs\n", length(gd_clumped_rsids)))
}

exp2_gd <- gd_full[rsid %in% gd_clumped_rsids]
rm(gd_full, gd_sig)
gc(verbose = FALSE)

cat("\n=== 4.2 Retrieving TSHR eQTL effects for all IVs ===\n")
all_iv_rsids <- unique(c(exp1_tshr$rsid, exp2_gd$rsid))
cat(sprintf("  Total unique IVs: %d\n", length(all_iv_rsids)))

tshr_all <- tryCatch(
    ieugwasr::associations(
        variants = all_iv_rsids,
        id = "eqtl-a-ENSG00000165409",
        opengwas_jwt = token
    ),
    error = function(e) {
        cat(sprintf("  [ERROR] %s\n", e$message))
        NULL
    }
)
if (is.null(tshr_all)) {
    sink()
    stop("TSHR associations failed")
}
cat(sprintf("  TSHR effects retrieved: %d\n", nrow(tshr_all)))

cat("\n=== 4.3 Retrieving GD effects for all IVs (local) ===\n")
gd_all <- extract_finngen_snps(
    file.path(FG_DIR, FG_FILES$GD_exposure),
    all_iv_rsids, FG_META$GD_exposure$label
)

cat("\n=== 4.4 Retrieving TED outcome effects for all IVs (local) ===\n")
ted_all <- extract_finngen_snps(
    file.path(FG_DIR, FG_FILES$TED_main),
    all_iv_rsids, FG_META$TED_main$label
)

cat("\n=== 4.5 Harmonizing alleles ===\n")
final_snps <- Reduce(intersect, list(tshr_all$rsid, gd_all$rsid, ted_all$rsid))
cat(sprintf("  SNPs in all 3 datasets: %d\n", length(final_snps)))

t_sub <- tshr_all[match(final_snps, tshr_all$rsid), ]
g_sub <- gd_all[match(final_snps, gd_all$rsid)] # data.table index matching subset
o_sub <- ted_all[match(final_snps, ted_all$rsid)]

mvmr_df <- data.table(
    SNP    = final_snps,
    beta_T = t_sub$beta,   se_T = t_sub$se,
    ea_T   = t_sub$ea,     oa_T = t_sub$nea,
    beta_G = g_sub$beta,   se_G = g_sub$sebeta,
    ea_G   = g_sub$alt,    oa_G = g_sub$ref,
    beta_Y = o_sub$beta,   se_Y = o_sub$sebeta,
    ea_Y   = o_sub$alt,    oa_Y = o_sub$ref
)

flip_g <- mvmr_df$ea_T != mvmr_df$ea_G & mvmr_df$ea_T == mvmr_df$oa_G
mvmr_df[flip_g, beta_G := -beta_G]

flip_y <- mvmr_df$ea_T != mvmr_df$ea_Y & mvmr_df$ea_T == mvmr_df$oa_Y
mvmr_df[flip_y, beta_Y := -beta_Y]

keep <- (mvmr_df$ea_T == mvmr_df$ea_G | mvmr_df$ea_T == mvmr_df$oa_G) &
    (mvmr_df$ea_T == mvmr_df$ea_Y | mvmr_df$ea_T == mvmr_df$oa_Y)
mvmr_df <- mvmr_df[keep]
mvmr_df <- mvmr_df[complete.cases(mvmr_df[, .(beta_T, beta_G, beta_Y, se_T, se_G, se_Y)])]

cat(sprintf("  Final harmonized SNPs: %d\n", nrow(mvmr_df)))
cat(sprintf("  Flipped GD: %d; Flipped TED: %d\n", sum(flip_g), sum(flip_y)))

cat("\n=== 4.6 MVMR-IVW ===\n\n")
wt <- 1 / (mvmr_df$se_Y^2)
fit_mvmr <- lm(beta_Y ~ 0 + beta_T + beta_G, data = mvmr_df, weights = wt)
sum_mvmr <- summary(fit_mvmr)
print(sum_mvmr)

tshr_b <- coef(sum_mvmr)[1, "Estimate"]
tshr_se <- coef(sum_mvmr)[1, "Std. Error"]
tshr_p <- coef(sum_mvmr)[1, "Pr(>|t|)"]
gd_b <- coef(sum_mvmr)[2, "Estimate"]
gd_se <- coef(sum_mvmr)[2, "Std. Error"]
gd_p <- coef(sum_mvmr)[2, "Pr(>|t|)"]

cat("\n=== 4.7 Conditional F-statistics ===\n")
mod_T <- lm(beta_T ~ 0 + beta_G, data = mvmr_df, weights = 1 / mvmr_df$se_T^2)
cond_F_T <- mean((residuals(mod_T) / mvmr_df$se_T)^2)
mod_G <- lm(beta_G ~ 0 + beta_T, data = mvmr_df, weights = 1 / mvmr_df$se_G^2)
cond_F_G <- mean((residuals(mod_G) / mvmr_df$se_G)^2)

cat(sprintf("  Cond. F (TSHR | GD): %.2f %s\n", cond_F_T, ifelse(cond_F_T > 10, "✅", "⚠️")))
cat(sprintf("  Cond. F (GD | TSHR): %.2f %s\n", cond_F_G, ifelse(cond_F_G > 10, "✅", "⚠️")))

cat("\n=== 4.8 Univariable MR (TSHR → TED) reference ===\n")
tshr_only <- mvmr_df[SNP %in% exp1_tshr$rsid]
cat(sprintf("  TSHR-only IVs used: %d\n", nrow(tshr_only)))

if (nrow(tshr_only) >= 1) {
    if (nrow(tshr_only) == 1) {
        uv_b <- tshr_only$beta_Y / tshr_only$beta_T
        uv_se <- tshr_only$se_Y / abs(tshr_only$beta_T)
        uv_p <- 2 * pnorm(-abs(uv_b / uv_se))
        cat(sprintf("  Wald ratio: β=%+.4f  SE=%.4f  P=%.3e\n", uv_b, uv_se, uv_p))
    } else {
        uv_fit <- lm(beta_Y ~ 0 + beta_T, data = tshr_only, weights = 1 / tshr_only$se_Y^2)
        uv_b <- coef(uv_fit)[1]
        uv_se <- summary(uv_fit)$coefficients[1, "Std. Error"]
        uv_p <- summary(uv_fit)$coefficients[1, "Pr(>|t|)"]
        cat(sprintf("  IVW: β=%+.4f  SE=%.4f  P=%.3e\n", uv_b, uv_se, uv_p))
    }
} else {
    uv_b <- NA
    uv_se <- NA
    uv_p <- NA
}

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("SE1 VERDICT — TED-SPECIFIC SIGNAL (FinnGen R12 TED outcome)\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("MVMR estimates:\n")
cat(sprintf(
    "  TSHR       β_adj = %+.4f  SE = %.4f  P = %.3e  %s\n",
    tshr_b, tshr_se, tshr_p, ifelse(tshr_p < 0.001, "***", ifelse(tshr_p < 0.05, "*", ""))
))
cat(sprintf(
    "  GD liab    β     = %+.4f  SE = %.4f  P = %.3e  %s\n",
    gd_b, gd_se, gd_p, ifelse(gd_p < 0.001, "***", ifelse(gd_p < 0.05, "*", ""))
))

cat("\nUnivariable TSHR → TED reference:\n")
cat(sprintf("  TSHR       β_uv  = %+.4f  SE = %.4f  P = %.3e\n", uv_b, uv_se, uv_p))

if (!is.na(uv_b) && uv_b != 0) {
    attn <- abs(tshr_b / uv_b)
    cat(sprintf("\nAttenuation: |β_adj / β_uv| = %.3f\n", attn))
    cat("\n--- SE1 Interpretation ---\n")
    if (tshr_p < 0.05 && attn > 0.7) {
        cat("✅ SE1 ACHIEVED: TSHR retains TED-direct effect independent of GD liability.\n   Narrative: TSHR → TED is direct, not solely via GD.\n")
    } else if (tshr_p < 0.05 && attn > 0.3) {
        cat("🟡 SE1 PARTIAL: TSHR has both mediated and direct TED effects.\n")
    } else if (attn < 0.3 && gd_p < 0.05) {
        cat("⚠️  TSHR → TED effect largely mediated by GD liability.\n   Narrative: TSHR → GD → TED (sequential pathway).\n")
    } else {
        cat("❓ Inconclusive (check power / IV strength)\n")
    }
}

mvmr_summary <- data.frame(Analysis = c("MVMR_TSHR_adj", "MVMR_GD_liability", "UVMR_TSHR_only"), Beta = c(tshr_b, gd_b, uv_b), SE = c(tshr_se, gd_se, uv_se), P = c(tshr_p, gd_p, uv_p), N_SNPs = c(nrow(mvmr_df), nrow(mvmr_df), nrow(tshr_only)), Cond_F = c(cond_F_T, cond_F_G, NA), Outcome = "FinnGen_R12_GRAVES_OPHT", stringsAsFactors = FALSE)
write.csv(mvmr_summary, "TrackA_MR/results/04v3_mvmr_finngen_summary.csv", row.names = FALSE)
saveRDS(list(mvmr_df = mvmr_df, fit = fit_mvmr, cond_F_T = cond_F_T, cond_F_G = cond_F_G, uv_b = uv_b, uv_p = uv_p), "TrackA_MR/results/04v3_mvmr_finngen_full.rds")

cat("\n💾 Saved:\n  TrackA_MR/results/03b_coloc_FinnGen_GRAVES_OPHT.csv\n  TrackA_MR/results/04v3_mvmr_finngen_summary.csv\n")
sink()
cat("\n✅ Step 3b + Step 4v3 complete.\n")
