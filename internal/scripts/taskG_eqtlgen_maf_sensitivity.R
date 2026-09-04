#!/usr/bin/env Rscript
# =============================================================================
# Task 6 — eQTLGen allele-frequency sensitivity analysis
# =============================================================================
# Exposure effect sizes throughout this study were reconstructed as
#
#     se   = 1 / sqrt(2 * eaf * (1 - eaf) * (N + Z^2))      (Zhu 2016)
#     beta = Z * se
#
# using 1000 Genomes European allele frequencies (n = 503) for `eaf`. eQTLGen
# distributes its own allele-frequency file (n = 26,609), which is the better
# source. This script re-derives beta/se from that file and re-runs the analyses
# that depend on them, so the manuscript can state whether the choice matters.
#
# What can and cannot move, worked out in advance so the run only has to confirm
# magnitudes:
#
#   * beta and se scale by the SAME factor sqrt(eaf_old(1-eaf_old) /
#     eaf_new(1-eaf_new)), so Z = beta/se is EXACTLY preserved. Single-instrument
#     Wald-ratio estimates therefore cannot change direction or P value. That
#     covers TSHR in all three outcomes and CTLA4 in discovery.
#   * IVW weights are 1/SE^2 and the scaling is per-SNP, so multi-instrument
#     pooled estimates CAN move. IGF1R (4 instruments) is the case that matters.
#   * coloc.abf here consumes beta/varbeta, so every posterior depends on eaf.
#   * Table S8 is computed from the observed standard errors and moves with them.
#
# Run from the repository root on the machine holding the licence-bound data:
#     Rscript internal/scripts/taskG_eqtlgen_maf_sensitivity.R
#
# Writes:
#   internal/EQTLGEN_MAF_SENSITIVITY.md                 the report
#   TrackA_MR/v5_upgrade/08_maf_sensitivity/*.csv       per-analysis comparisons
#
# NOTE: written in the remote container without access to the data, so it is
# untested against real inputs. It asserts its assumptions loudly rather than
# proceeding on a guess — if an assertion fires, fix the path or the column name
# rather than removing the assertion.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(TwoSampleMR)
})

# ---------------------------- configuration ---------------------------------
WORK    <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
EQTLGEN <- "c:/ProjectTEDGWAS/TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"

# eQTLGen's own allele frequencies. Download from
#   https://www.eqtlgen.org/cis-eqtls.html  ->  "Allele frequency file"
# and point this at it. Expected columns include SNP, AlleleA/AlleleB (or
# AssessedAllele/OtherAllele) and an allele frequency for the assessed allele.
EQTLGEN_AF <- "c:/ProjectTEDGWAS/TrackA_MR/data/2018-07-18_SNP_AF_for_AlleleB_combination_60_never_combined_where_AlleleB_upper.txt.gz"

FRQ_EUR  <- file.path(WORK, "04_druggable_mr/data/g1000_eur_freq.frq")
MANIFEST <- file.path(WORK, "04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v2_verified.csv")
OUTDIR   <- file.path(WORK, "08_maf_sensitivity")
REPORT   <- "internal/EQTLGEN_MAF_SENSITIVITY.md"

OUTCOMES <- list(
  BBJ_Graves       = list(path = file.path(WORK, "04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz"),
                          cases = 2809L, N = 175465L, lmm = FALSE),
  UKB_hyperthyroid = list(path = file.path(WORK, "04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz"),
                          cases = 3731L, N = 484598L, lmm = TRUE),
  FinnGen_GO       = list(path = "c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz",
                          cases = 858L,  N = 500348L, lmm = FALSE)
)

BACKBONE   <- c("TSHR", "IGF1R", "CTLA4")
CANDIDATES <- c("TNFSF14", "IFNGR1", "MAPKAPK5", "HSD3B7", "VKORC1", "PRSS36")
GENES      <- c(BACKBONE, CANDIDATES)
WINDOW     <- 1000000L
PRIORS     <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(MANIFEST), file.exists(FRQ_EUR), file.exists(EQTLGEN))
if (!file.exists(EQTLGEN_AF))
  stop("eQTLGen allele-frequency file not found at:\n  ", EQTLGEN_AF,
       "\nDownload it from https://www.eqtlgen.org/cis-eqtls.html and set EQTLGEN_AF.")

cat("== Task 6: eQTLGen allele-frequency sensitivity ==\n")
cat("R:", R.version.string, " coloc:", as.character(packageVersion("coloc")), "\n")
cat("run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ------------------------- allele-frequency sources --------------------------
frq_eur <- fread(FRQ_EUR)                       # CHR SNP A1 A2 MAF NCHROBS
af_gen  <- fread(EQTLGEN_AF)
cat("eQTLGen AF file columns:", paste(names(af_gen), collapse = ", "), "\n")

# Normalise the eQTLGen AF column names. The distributed file has varied across
# releases, so match on what is actually present instead of assuming.
snp_col <- intersect(c("SNP", "rsid", "RSID"), names(af_gen))[1]
af_col  <- intersect(c("AlleleB_all", "AlleleB_freq", "AF", "allele_freq",
                       "AssessedAllele_freq"), names(af_gen))[1]
b_col   <- intersect(c("AlleleB", "AssessedAllele", "EffectAllele"), names(af_gen))[1]
a_col   <- intersect(c("AlleleA", "OtherAllele"), names(af_gen))[1]
if (any(is.na(c(snp_col, af_col, b_col, a_col))))
  stop("Could not identify SNP / allele / frequency columns in the eQTLGen AF file.\n",
       "Present: ", paste(names(af_gen), collapse = ", "),
       "\nSet snp_col / af_col / b_col / a_col by hand and re-run.")
setnames(af_gen, c(snp_col, af_col, b_col, a_col), c("SNP", "afB", "alleleB", "alleleA"))
af_gen <- af_gen[, .(SNP, alleleA, alleleB, afB)]
cat(sprintf("eQTLGen AF: %s rows, using SNP=%s allele=%s freq=%s\n\n",
            format(nrow(af_gen), big.mark = ","), snp_col, b_col, af_col))

# ------------------------------ instruments ----------------------------------
inst <- fread(MANIFEST)
stopifnot(nrow(inst) == 6135L)                  # v2 verified manifest
stopifnot(all(inst$pvalue < 5e-8))
cat(sprintf("instruments: %d rows, %d genes\n", nrow(inst), uniqueN(inst$gene_symbol)))

#' Effect-allele frequency under each source, aligned to the eQTLGen effect allele.
eaf_both <- function(dt) {
  m <- merge(dt, frq_eur[, .(SNP, fA1 = A1, fA2 = A2, MAF)],
             by.x = "snp", by.y = "SNP", all.x = TRUE)
  m[, eaf_1kg := fifelse(effect_allele == fA1, MAF,
                  fifelse(effect_allele == fA2, 1 - MAF, NA_real_))]
  m <- merge(m, af_gen, by.x = "snp", by.y = "SNP", all.x = TRUE)
  m[, eaf_gen := fifelse(effect_allele == alleleB, afB,
                  fifelse(effect_allele == alleleA, 1 - afB, NA_real_))]
  m
}

ei <- eaf_both(inst)
n_miss_1kg <- sum(is.na(ei$eaf_1kg)); n_miss_gen <- sum(is.na(ei$eaf_gen))
both <- ei[!is.na(eaf_1kg) & !is.na(eaf_gen)]
both[, d := eaf_gen - eaf_1kg]
both[, scale := sqrt((eaf_1kg * (1 - eaf_1kg)) / (eaf_gen * (1 - eaf_gen)))]

cat(sprintf("\nallele-frequency comparison (n = %d instruments with both):\n", nrow(both)))
cat(sprintf("  missing eaf : 1000G %d, eQTLGen %d\n", n_miss_1kg, n_miss_gen))
cat(sprintf("  |difference|: median %.4f, 90th pct %.4f, max %.4f\n",
            median(abs(both$d)), quantile(abs(both$d), .9), max(abs(both$d))))
cat(sprintf("  beta/se scaling factor: median %.4f, range %.4f-%.4f\n",
            median(both$scale), min(both$scale), max(both$scale)))
fwrite(both[, .(gene_symbol, snp, effect_allele, zscore, pvalue, n_samples,
                eaf_1kg, eaf_gen, diff = d, scale)],
       file.path(OUTDIR, "TaskG_01_af_comparison.csv"))

# The invariance claim in Limitations, checked rather than asserted.
both[, `:=`(se_1kg = 1 / sqrt(2 * eaf_1kg * (1 - eaf_1kg) * (n_samples + zscore^2)),
            se_gen = 1 / sqrt(2 * eaf_gen * (1 - eaf_gen) * (n_samples + zscore^2)))]
both[, `:=`(b_1kg = zscore * se_1kg, b_gen = zscore * se_gen)]
zdiff <- max(abs(both$b_1kg / both$se_1kg - both$b_gen / both$se_gen))
cat(sprintf("  max |Z_1kg - Z_eQTLGen| = %.3e  (should be ~0: Z is invariant)\n\n", zdiff))
stopifnot(zdiff < 1e-8)

# ------------------------------- outcomes ------------------------------------
load_outcome <- function(oc) {
  d <- fread(oc$path)
  setnames(d, tolower(gsub("^#", "", names(d))))
  rs <- intersect(c("rsid", "rsids", "variant_id", "snp"), names(d))[1]
  bc <- intersect(c("beta"), names(d))[1]
  sc <- intersect(c("standard_error", "sebeta", "se"), names(d))[1]
  ec <- intersect(c("effect_allele", "alt"), names(d))[1]
  oc_ <- intersect(c("other_allele", "ref"), names(d))[1]
  if (any(is.na(c(rs, bc, sc, ec, oc_))))
    stop("outcome columns not identified in ", oc$path, ": ",
         paste(names(d), collapse = ", "))
  out <- d[, .(SNP = get(rs), beta = as.numeric(get(bc)), se = as.numeric(get(sc)),
               ea = toupper(get(ec)), oa = toupper(get(oc_)))]
  out <- out[!is.na(beta) & !is.na(se) & se > 0 & SNP != "" & !is.na(SNP)]
  if (isTRUE(oc$lmm)) {                         # BOLT-LMM -> log-odds
    mu <- oc$cases / oc$N; k <- mu * (1 - mu)
    out[, `:=`(beta = beta / k, se = se / k)]   # Z and P preserved
  }
  unique(out, by = "SNP")
}
cat("loading outcomes...\n")
OUT_DT <- lapply(OUTCOMES, load_outcome)
for (nm in names(OUT_DT)) cat(sprintf("  %-18s %s variants\n", nm, format(nrow(OUT_DT[[nm]]), big.mark = ",")))

# --------------------------------- MR ----------------------------------------
mr_one <- function(gene, oc_name, which_eaf) {
  e <- both[gene_symbol == gene]
  if (!nrow(e)) return(NULL)
  e[, `:=`(beta.exposure = if (which_eaf == "1kg") b_1kg else b_gen,
           se.exposure   = if (which_eaf == "1kg") se_1kg else se_gen)]
  o <- OUT_DT[[oc_name]]
  m <- merge(e[, .(SNP = snp, effect_allele.exposure = effect_allele,
                   other_allele.exposure = other_allele,
                   beta.exposure, se.exposure)],
             o[, .(SNP, beta.outcome = beta, se.outcome = se,
                   effect_allele.outcome = ea, other_allele.outcome = oa)], by = "SNP")
  if (!nrow(m)) return(NULL)
  swap <- m$effect_allele.exposure == m$other_allele.outcome
  m[swap, beta.outcome := -beta.outcome]
  m[, `:=`(id.exposure = gene, id.outcome = oc_name,
           exposure = gene, outcome = oc_name, mr_keep = TRUE)]
  res <- suppressMessages(mr(as.data.frame(m),
           method_list = if (nrow(m) == 1) "mr_wald_ratio" else "mr_ivw"))
  data.table(gene = gene, outcome = oc_name, af = which_eaf, nsnp = nrow(m),
             beta = res$b[1], se = res$se[1], pval = res$pval[1],
             or = exp(res$b[1]))
}

cat("\nrunning MR under both allele-frequency sources...\n")
mr_res <- rbindlist(lapply(GENES, function(g)
  rbindlist(lapply(names(OUTCOMES), function(o)
    rbindlist(lapply(c("1kg", "gen"), function(a) mr_one(g, o, a)))))), fill = TRUE)
fwrite(mr_res, file.path(OUTDIR, "TaskG_02_mr_both_af.csv"))

mr_wide <- dcast(mr_res, gene + outcome + nsnp ~ af,
                 value.var = c("beta", "se", "pval", "or"))
mr_wide[, `:=`(d_beta = beta_gen - beta_1kg, d_or = or_gen - or_1kg,
               logp_ratio = -log10(pval_gen) + log10(pval_1kg))]
fwrite(mr_wide, file.path(OUTDIR, "TaskG_03_mr_comparison.csv"))

# ---------------------------- colocalization ---------------------------------
cat("loading eQTLGen (several minutes)...\n")
eq <- fread(EQTLGEN, select = c("Pvalue", "SNP", "SNPChr", "SNPPos", "GeneSymbol",
                                "NrSamples", "Zscore", "AssessedAllele", "OtherAllele"))
GENE_POS <- unique(inst[, .(gene_symbol, chr, pos_hg19)])[
  , .(chr = chr[1], lo = min(pos_hg19) - WINDOW, hi = max(pos_hg19) + WINDOW),
  by = gene_symbol]

coloc_one <- function(gene, oc_name, which_eaf) {
  gp <- GENE_POS[gene_symbol == gene]
  if (!nrow(gp)) return(NULL)
  e <- eq[GeneSymbol == gene & SNPChr == gp$chr & SNPPos >= gp$lo & SNPPos <= gp$hi]
  if (nrow(e) < 50) return(NULL)
  e <- eaf_both(e[, .(snp = SNP, effect_allele = AssessedAllele,
                      other_allele = OtherAllele, zscore = Zscore,
                      n_samples = NrSamples, pvalue = Pvalue)])
  e <- e[!is.na(eaf_1kg) & !is.na(eaf_gen)]
  eaf <- if (which_eaf == "1kg") e$eaf_1kg else e$eaf_gen
  e[, se_e := 1 / sqrt(2 * eaf * (1 - eaf) * (n_samples + zscore^2))]
  e[, beta_e := zscore * se_e]

  o <- OUT_DT[[oc_name]]
  m <- merge(e[, .(SNP = snp, ea_e = effect_allele, oa_e = other_allele, beta_e, se_e)],
             o, by = "SNP")
  if (nrow(m) < 50) return(NULL)
  swap <- m$ea_e == m$oa
  m[swap, beta := -beta]
  oc <- OUTCOMES[[oc_name]]
  r <- tryCatch(coloc.abf(
        dataset1 = list(snp = m$SNP, beta = m$beta, varbeta = m$se^2,
                        type = "cc", s = oc$cases / oc$N, N = oc$N),
        dataset2 = list(snp = m$SNP, beta = m$beta_e, varbeta = m$se_e^2,
                        type = "quant", N = 31684L),
        p1 = PRIORS$p1, p2 = PRIORS$p2, p12 = PRIORS$p12),
      error = function(e) NULL)
  if (is.null(r)) return(NULL)
  s <- r$summary
  data.table(gene = gene, outcome = oc_name, af = which_eaf, n_overlap = nrow(m),
             PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]],
             PP.H2 = s[["PP.H2.abf"]], PP.H3 = s[["PP.H3.abf"]],
             PP.H4 = s[["PP.H4.abf"]],
             top_snp = r$results$snp[which.max(r$results$SNP.PP.H4)])
}

cat("running colocalization under both allele-frequency sources...\n")
co_res <- rbindlist(lapply(GENES, function(g)
  rbindlist(lapply(names(OUTCOMES), function(o)
    rbindlist(lapply(c("1kg", "gen"), function(a) coloc_one(g, o, a)))))), fill = TRUE)
fwrite(co_res, file.path(OUTDIR, "TaskG_04_coloc_both_af.csv"))

co_wide <- dcast(co_res, gene + outcome ~ af,
                 value.var = c("PP.H2", "PP.H3", "PP.H4", "top_snp", "n_overlap"))
co_wide[, d_H4 := PP.H4_gen - PP.H4_1kg]
co_wide[, crosses := (PP.H4_1kg >= 0.8) != (PP.H4_gen >= 0.8)]
co_wide[, top_changed := top_snp_1kg != top_snp_gen]
fwrite(co_wide, file.path(OUTDIR, "TaskG_05_coloc_comparison.csv"))

# -------------------------------- report -------------------------------------
max_dH4 <- max(abs(co_wide$d_H4), na.rm = TRUE)
n_cross <- sum(co_wide$crosses, na.rm = TRUE)
n_topch <- sum(co_wide$top_changed, na.rm = TRUE)
mr_sig  <- mr_wide[, sum((pval_1kg < 0.05) != (pval_gen < 0.05), na.rm = TRUE)]

lines <- c(
  "# eQTLGen allele-frequency sensitivity analysis",
  "",
  sprintf("_Run %s. Script: `internal/scripts/taskG_eqtlgen_maf_sensitivity.R`._", Sys.Date()),
  "",
  "Exposure effect sizes are reconstructed as `se = 1/sqrt(2·eaf·(1-eaf)·(N + Z²))`,",
  "`beta = Z·se`. This compares 1000 Genomes European frequencies (n = 503), used",
  "throughout the manuscript, against eQTLGen's own frequencies (n = 26,609).",
  "",
  "## Allele frequencies",
  "",
  sprintf("- instruments with both frequencies: **%d** of %d", nrow(both), nrow(inst)),
  sprintf("- missing: 1000G %d, eQTLGen %d", n_miss_1kg, n_miss_gen),
  sprintf("- |difference|: median %.4f, 90th percentile %.4f, max %.4f",
          median(abs(both$d)), quantile(abs(both$d), .9), max(abs(both$d))),
  sprintf("- beta/se scaling factor: median %.4f (range %.4f–%.4f)",
          median(both$scale), min(both$scale), max(both$scale)),
  sprintf("- max |ΔZ| = %.3e — Z is invariant, as expected", zdiff),
  "",
  "## Mendelian randomization",
  "",
  sprintf("- estimates whose nominal significance changes at P = 0.05: **%d**", mr_sig),
  "- single-instrument loci cannot move (Z invariant); any change is confined to",
  "  multi-instrument loci, where inverse-variance weights depend on the frequency.",
  "",
  "## Colocalization",
  "",
  sprintf("- maximum |ΔPP.H4| across %d gene × outcome pairs: **%.4f**", nrow(co_wide), max_dH4),
  sprintf("- pairs crossing the PP.H4 = 0.80 threshold: **%d**", n_cross),
  sprintf("- pairs whose top-posterior variant changes: **%d**", n_topch),
  "",
  "## Verdict",
  "",
  if (n_cross == 0 && mr_sig == 0)
    paste0("No estimate crosses a reporting threshold and no colocalization ",
           "conclusion changes. The manuscript's caveat can be shortened to a ",
           "statement that the sensitivity analysis was performed and the ",
           "conclusions were unchanged.")
  else
    paste0("**Conclusions move.** ", n_cross, " colocalization pair(s) cross the ",
           "0.80 threshold and ", mr_sig, " MR estimate(s) change nominal ",
           "significance. Re-report the affected numbers from the eQTLGen-frequency ",
           "run, which is the better primary analysis."),
  "",
  "## Files",
  "",
  "- `TaskG_01_af_comparison.csv` — per-instrument frequencies and scaling factor",
  "- `TaskG_02_mr_both_af.csv`, `TaskG_03_mr_comparison.csv`",
  "- `TaskG_04_coloc_both_af.csv`, `TaskG_05_coloc_comparison.csv`"
)
writeLines(lines, REPORT)

cat(sprintf("\n== summary ==\n  max |dPP.H4| = %.4f\n  threshold crossings = %d\n"
            , max_dH4, n_cross))
cat(sprintf("  MR significance changes = %d\n  top-SNP changes = %d\n", mr_sig, n_topch))
cat("\nwrote", REPORT, "and 5 CSVs in", OUTDIR, "\n")
