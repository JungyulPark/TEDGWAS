# ═══════════════════════════════════════════════════════════════
# IGF1R cis-MR on BBJ Graves Disease — VERIFICATION SCRIPT
# Purpose: Determine which set of P-values is the true ground truth
#
# DISCREPANCY DETECTED:
#   Method            Figure 2     Supp Table S1
#   IVW               0.089        0.089          ✓ match
#   MR-Egger          0.221        0.812          ✗ MISMATCH
#   Weighted median   0.103        0.168          ✗ MISMATCH  
#   Weighted mode     0.621        (missing)      ?
#
# This script re-runs the analysis from scratch with full reproducibility.
# Output will determine the correct values to use in submission.
#
# Author: Claude (TED-TRAP project)
# Date: 2026-04-29
# ═══════════════════════════════════════════════════════════════

# ── Required packages ─────────────────────────────────────────
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(dplyr)
})

# ── Reproducibility ───────────────────────────────────────────
set.seed(42)

# ── Configuration ─────────────────────────────────────────────
EXPOSURE_GENE      <- "IGF1R"
EXPOSURE_ENSEMBL   <- "ENSG00000140443"
EXPOSURE_ID        <- paste0("eqtl-a-", EXPOSURE_ENSEMBL)
OUTCOME_ID         <- "ebi-a-GCST90018627"   # BBJ Graves disease
OUT_DIR            <- "C:/ProjectTEDGWAS/TrackA_MR/verification_IGF1R_BBJ"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n══════════════════════════════════════════════════════════\n")
cat("  IGF1R cis-MR on BBJ Graves — Re-run Verification\n")
cat("══════════════════════════════════════════════════════════\n")
cat("Exposure:  ", EXPOSURE_ID, "(", EXPOSURE_GENE, ")\n")
cat("Outcome:   ", OUTCOME_ID, "(BBJ Graves disease)\n")
cat("Output dir:", OUT_DIR, "\n\n")

# ── Step 1: Extract exposure instruments ─────────────────────
cat("─── Step 1: Extracting cis-eQTL instruments ───\n")
exposure <- tryCatch({
  extract_instruments(
    outcomes = EXPOSURE_ID,
    p1       = 5e-8,
    clump    = TRUE,
    r2       = 0.001,
    kb       = 10000
  )
}, error = function(e) {
  cat("Error at p<5e-8:", conditionMessage(e), "\n")
  NULL
})

if (is.null(exposure) || nrow(exposure) == 0) {
  cat("Retrying with p<5e-6 (eQTL conventional threshold)...\n")
  exposure <- extract_instruments(
    outcomes = EXPOSURE_ID,
    p1       = 5e-6,
    clump    = TRUE,
    r2       = 0.001,
    kb       = 10000
  )
}

n_iv_pre <- nrow(exposure)
cat("IVs extracted before harmonization:", n_iv_pre, "\n")
cat("IV F-stat range:",
    sprintf("%.1f – %.1f", 
            min((exposure$beta.exposure / exposure$se.exposure)^2),
            max((exposure$beta.exposure / exposure$se.exposure)^2)),
    "\n\n")

write.csv(exposure,
          file.path(OUT_DIR, "01_IGF1R_instruments_eQTLGen.csv"),
          row.names = FALSE)

# ── Step 2: Extract outcome data ─────────────────────────────
cat("─── Step 2: Extracting BBJ Graves outcome ───\n")
outcome <- extract_outcome_data(
  snps     = exposure$SNP,
  outcomes = OUTCOME_ID
)
cat("Outcome SNPs found:", nrow(outcome), "/ requested:", n_iv_pre, "\n\n")

# ── Step 3: Harmonize ────────────────────────────────────────
cat("─── Step 3: Harmonizing exposure-outcome ───\n")
harmonised <- harmonise_data(exposure, outcome, action = 2)
harmonised_kept <- harmonised[harmonised$mr_keep == TRUE, ]
cat("SNPs kept (mr_keep == TRUE):", nrow(harmonised_kept), "\n")
cat("SNPs dropped (palindromic/etc):",
    nrow(harmonised) - nrow(harmonised_kept), "\n\n")

write.csv(harmonised,
          file.path(OUT_DIR, "02_IGF1R_BBJ_harmonized_full.csv"),
          row.names = FALSE)
write.csv(harmonised_kept,
          file.path(OUT_DIR, "03_IGF1R_BBJ_harmonized_kept.csv"),
          row.names = FALSE)

# ── Step 4: Run ALL 5 MR estimators ──────────────────────────
cat("─── Step 4: Running MR with all 5 estimators ───\n")
mr_results <- mr(
  harmonised_kept,
  method_list = c("mr_wald_ratio",
                  "mr_ivw",
                  "mr_egger_regression",
                  "mr_weighted_median",
                  "mr_weighted_mode")
)

cat("\n┌──────────────────────────────────────────────────────────┐\n")
cat("│  RE-RUN RESULTS — IGF1R cis-MR on BBJ Graves            │\n")
cat("└──────────────────────────────────────────────────────────┘\n")
print(mr_results[, c("method", "nsnp", "b", "se", "pval")],
      row.names = FALSE)

write.csv(mr_results,
          file.path(OUT_DIR, "04_IGF1R_BBJ_MR_all5_estimators.csv"),
          row.names = FALSE)

# ── Step 5: Sensitivity statistics ───────────────────────────
cat("\n─── Step 5: Sensitivity tests ───\n")
het <- mr_heterogeneity(harmonised_kept)
plt <- mr_pleiotropy_test(harmonised_kept)

cat("\nHeterogeneity (Cochran Q):\n")
print(het, row.names = FALSE)
cat("\nMR-Egger intercept (directional pleiotropy):\n")
print(plt, row.names = FALSE)

write.csv(het, file.path(OUT_DIR, "05_heterogeneity.csv"), row.names = FALSE)
write.csv(plt, file.path(OUT_DIR, "06_egger_intercept.csv"), row.names = FALSE)

# ── Step 6: SIDE-BY-SIDE COMPARISON ──────────────────────────
cat("\n══════════════════════════════════════════════════════════\n")
cat("  GROUND-TRUTH DETERMINATION\n")
cat("══════════════════════════════════════════════════════════\n")

get_p <- function(method_name) {
  v <- mr_results$pval[mr_results$method == method_name]
  if (length(v) == 0) NA else v
}

ivw_p   <- get_p("Inverse variance weighted")
egger_p <- get_p("MR Egger")
wm_p    <- get_p("Weighted median")
wmod_p  <- get_p("Weighted mode")

format_p <- function(p) {
  if (is.na(p)) return("    NA  ")
  if (p < 0.001) return(sprintf("%.2e", p))
  return(sprintf("%.4f", p))
}

cat("\nMethod              Re-run      Figure 2     Supp S1\n")
cat("───────────────────────────────────────────────────────\n")
cat(sprintf("IVW                 %s    0.089        0.089\n",   format_p(ivw_p)))
cat(sprintf("MR-Egger            %s    0.221        0.812\n",   format_p(egger_p)))
cat(sprintf("Weighted median     %s    0.103        0.168\n",   format_p(wm_p)))
cat(sprintf("Weighted mode       %s    0.621        (missing)\n", format_p(wmod_p)))
cat("───────────────────────────────────────────────────────\n\n")

# Verdict logic
verdict_egger <- if (is.na(egger_p)) "Cannot compute" else {
  d_fig2 <- abs(egger_p - 0.221)
  d_supp <- abs(egger_p - 0.812)
  if (d_fig2 < d_supp) "Figure 2 matches" else "Supp S1 matches"
}
verdict_wm <- if (is.na(wm_p)) "Cannot compute" else {
  d_fig2 <- abs(wm_p - 0.103)
  d_supp <- abs(wm_p - 0.168)
  if (d_fig2 < d_supp) "Figure 2 matches" else "Supp S1 matches"
}

cat("VERDICT (closest match):\n")
cat("  MR-Egger         →", verdict_egger, "\n")
cat("  Weighted median  →", verdict_wm, "\n\n")

# ── Step 7: Save session info ────────────────────────────────
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt"))

cat("══════════════════════════════════════════════════════════\n")
cat("  DONE. Send Claude these files:\n")
cat("══════════════════════════════════════════════════════════\n")
cat("  1. 04_IGF1R_BBJ_MR_all5_estimators.csv  (MAIN result)\n")
cat("  2. 05_heterogeneity.csv\n")
cat("  3. 06_egger_intercept.csv\n")
cat("  4. The full console output above (copy-paste)\n")
cat("══════════════════════════════════════════════════════════\n")
