# ═══════════════════════════════════════════════════════════════
# IGF1R cis-MR on UKB Hyperthyroidism — VERIFICATION SCRIPT
# Purpose: Sanity check whether UKB MR-Egger and Weighted median 
#          values in Supp Table S1 are correct (after BBJ values 
#          were found to be wrong).
#
# CURRENT SUPP TABLE S1 VALUES (need verification):
#   Method           β        SE      OR (95% CI)         P
#   IVW              +0.152   0.091   1.16 (0.97–1.39)    0.094
#   MR-Egger         +0.089   0.211   1.09 (0.72–1.65)    0.681
#   Weighted median  +0.165   0.107   1.18 (0.96–1.45)    0.123
#
# IMPORTANT NOTE on UKB scaling:
#   UKB Hyperthyroidism is from Neale-lab style LMM. Raw β/SE will be
#   in LMM scale (small numbers, e.g., 0.001). The values in Supp Table
#   are RESCALED to log-odds: β_logodds = β_LMM / (μ × (1−μ))
#   where μ = case prevalence = 3731 / 484598 ≈ 0.0077.
#   Rescaling factor ≈ 1 / 0.00764 ≈ 130.9.
#   
#   ⚠️  P-VALUES are scale-invariant — these are the primary comparison.
#   β/SE comparison requires applying the rescaling factor.
#
# Author: Claude (TED-TRAP project)
# Date: 2026-04-29
# ═══════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(dplyr)
})
set.seed(42)

# ── Configuration ─────────────────────────────────────────────
EXPOSURE_GENE      <- "IGF1R"
EXPOSURE_ENSEMBL   <- "ENSG00000140443"
EXPOSURE_ID        <- paste0("eqtl-a-", EXPOSURE_ENSEMBL)
OUTCOME_ID         <- "ebi-a-GCST90038636"   # ← UKB Hyperthyroidism
OUT_DIR            <- "C:/ProjectTEDGWAS/TrackA_MR/verification_IGF1R_UKB"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# UKB rescaling factor (LMM → log-odds)
UKB_N_CASES        <- 3731
UKB_N_TOTAL        <- 484598
UKB_PREVALENCE     <- UKB_N_CASES / UKB_N_TOTAL
RESCALE_FACTOR     <- 1 / (UKB_PREVALENCE * (1 - UKB_PREVALENCE))

cat("\n══════════════════════════════════════════════════════════\n")
cat("  IGF1R cis-MR on UKB Hyperthyroidism — Sanity Check\n")
cat("══════════════════════════════════════════════════════════\n")
cat("Exposure:        ", EXPOSURE_ID, "(", EXPOSURE_GENE, ")\n")
cat("Outcome:         ", OUTCOME_ID, "(UKB Hyperthyroid)\n")
cat("UKB prevalence:  ", sprintf("%.4f", UKB_PREVALENCE), "\n")
cat("Rescale factor:  ", sprintf("%.2f", RESCALE_FACTOR), "(LMM → log-odds)\n")
cat("Output dir:      ", OUT_DIR, "\n\n")

# ── Step 1: Extract instruments ──────────────────────────────
cat("─── Step 1: Extracting cis-eQTL instruments ───\n")
exposure <- tryCatch({
  extract_instruments(outcomes = EXPOSURE_ID, p1 = 5e-8, clump = TRUE,
                      r2 = 0.001, kb = 10000)
}, error = function(e) NULL)

if (is.null(exposure) || nrow(exposure) == 0) {
  cat("Retrying with p<5e-6...\n")
  exposure <- extract_instruments(outcomes = EXPOSURE_ID, p1 = 5e-6,
                                  clump = TRUE, r2 = 0.001, kb = 10000)
}
cat("IVs extracted:", nrow(exposure), "\n\n")
write.csv(exposure, file.path(OUT_DIR, "01_IGF1R_instruments.csv"), row.names = FALSE)

# ── Step 2: Outcome ──────────────────────────────────────────
cat("─── Step 2: Extracting UKB hyperthyroid outcome ───\n")
outcome <- extract_outcome_data(snps = exposure$SNP, outcomes = OUTCOME_ID)
cat("Outcome SNPs found:", nrow(outcome), "\n\n")

# ── Step 3: Harmonize ────────────────────────────────────────
cat("─── Step 3: Harmonizing ───\n")
harmonised <- harmonise_data(exposure, outcome, action = 2)
harmonised_kept <- harmonised[harmonised$mr_keep == TRUE, ]
cat("SNPs kept:", nrow(harmonised_kept), "/", nrow(harmonised), "\n\n")
write.csv(harmonised_kept, file.path(OUT_DIR, "02_IGF1R_UKB_harmonized.csv"), row.names = FALSE)

# ── Step 4: Run all 5 MR estimators ──────────────────────────
cat("─── Step 4: Running all 5 MR estimators (RAW LMM SCALE) ───\n")
mr_results <- mr(harmonised_kept,
                 method_list = c("mr_wald_ratio", "mr_ivw",
                                 "mr_egger_regression",
                                 "mr_weighted_median",
                                 "mr_weighted_mode"))

# Add rescaled (log-odds) columns
mr_results$b_rescaled  <- mr_results$b  * RESCALE_FACTOR
mr_results$se_rescaled <- mr_results$se * RESCALE_FACTOR
mr_results$or_rescaled <- exp(mr_results$b_rescaled)
mr_results$or_lo       <- exp(mr_results$b_rescaled - 1.96 * mr_results$se_rescaled)
mr_results$or_hi       <- exp(mr_results$b_rescaled + 1.96 * mr_results$se_rescaled)

cat("\n┌──────────────────────────────────────────────────────────┐\n")
cat("│  RE-RUN RESULTS — IGF1R cis-MR on UKB Hyperthyroid       │\n")
cat("└──────────────────────────────────────────────────────────┘\n\n")

cat("RAW LMM scale (R direct output):\n")
print(mr_results[, c("method", "nsnp", "b", "se", "pval")], row.names = FALSE)

cat("\nRESCALED to log-odds (× ", sprintf("%.2f", RESCALE_FACTOR), "):\n", sep = "")
disp <- mr_results[, c("method", "nsnp", "b_rescaled", "se_rescaled", "pval")]
disp$or_ci <- sprintf("%.2f (%.2f-%.2f)", mr_results$or_rescaled,
                                          mr_results$or_lo,
                                          mr_results$or_hi)
print(disp[, c("method", "nsnp", "b_rescaled", "se_rescaled", "or_ci", "pval")],
      row.names = FALSE)

write.csv(mr_results, file.path(OUT_DIR, "03_IGF1R_UKB_MR_all5_estimators.csv"),
          row.names = FALSE)

# ── Step 5: Sensitivity ──────────────────────────────────────
het <- mr_heterogeneity(harmonised_kept)
plt <- mr_pleiotropy_test(harmonised_kept)
cat("\n─── Sensitivity ───\n")
cat("Heterogeneity:\n"); print(het, row.names = FALSE)
cat("Egger intercept:\n"); print(plt, row.names = FALSE)
write.csv(het, file.path(OUT_DIR, "04_heterogeneity.csv"), row.names = FALSE)
write.csv(plt, file.path(OUT_DIR, "05_egger_intercept.csv"), row.names = FALSE)

# ── Step 6: COMPARISON vs Supp Table S1 (current values) ─────
cat("\n══════════════════════════════════════════════════════════\n")
cat("  SANITY CHECK — Re-run vs current Supp Table S1 (UKB)\n")
cat("══════════════════════════════════════════════════════════\n\n")

get_p <- function(m) {
  v <- mr_results$pval[mr_results$method == m]
  if (length(v) == 0) NA else v
}
get_b_rs <- function(m) {
  v <- mr_results$b_rescaled[mr_results$method == m]
  if (length(v) == 0) NA else v
}
get_se_rs <- function(m) {
  v <- mr_results$se_rescaled[mr_results$method == m]
  if (length(v) == 0) NA else v
}

ivw_p   <- get_p("Inverse variance weighted")
ivw_b   <- get_b_rs("Inverse variance weighted")
ivw_se  <- get_se_rs("Inverse variance weighted")
egger_p <- get_p("MR Egger")
egger_b <- get_b_rs("MR Egger")
egger_se<- get_se_rs("MR Egger")
wm_p    <- get_p("Weighted median")
wm_b    <- get_b_rs("Weighted median")
wm_se   <- get_se_rs("Weighted median")

fp <- function(p) if (is.na(p)) "    NA  " else sprintf("%.4f", p)
fbs <- function(b) if (is.na(b)) "    NA " else sprintf("%+.3f", b)

cat("PRIMARY COMPARISON — P-values (scale-invariant):\n")
cat("───────────────────────────────────────────────────\n")
cat("Method              Re-run      Current Supp S1\n")
cat("───────────────────────────────────────────────────\n")
cat(sprintf("IVW                 %s      0.094\n", fp(ivw_p)))
cat(sprintf("MR-Egger            %s      0.681\n", fp(egger_p)))
cat(sprintf("Weighted median     %s      0.123\n", fp(wm_p)))
cat("───────────────────────────────────────────────────\n\n")

cat("SECONDARY COMPARISON — β rescaled to log-odds:\n")
cat("───────────────────────────────────────────────────\n")
cat("Method              Re-run β    Current Supp β\n")
cat("───────────────────────────────────────────────────\n")
cat(sprintf("IVW                 %s     +0.152\n", fbs(ivw_b)))
cat(sprintf("MR-Egger            %s     +0.089\n", fbs(egger_b)))
cat(sprintf("Weighted median     %s     +0.165\n", fbs(wm_b)))
cat("───────────────────────────────────────────────────\n\n")

# Verdict
verdict <- function(re, supp, label, tol_p = 0.02, tol_b = 0.05) {
  if (is.na(re)) return(paste(label, ": cannot compute"))
  diff <- abs(re - supp)
  if (diff < tol_p) return(paste(label, ": ✅ MATCH (Δ =", sprintf("%.4f", diff), ")"))
  return(paste(label, ": ❌ MISMATCH (Δ =", sprintf("%.4f", diff), ")"))
}

cat("VERDICT:\n")
cat("  IVW P              ", verdict(ivw_p,   0.094, "IVW          "), "\n")
cat("  MR-Egger P         ", verdict(egger_p, 0.681, "MR-Egger     "), "\n")
cat("  Weighted median P  ", verdict(wm_p,    0.123, "WM           "), "\n\n")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("══════════════════════════════════════════════════════════\n")
cat("  DONE. Send Claude:\n")
cat("══════════════════════════════════════════════════════════\n")
cat("  1. 03_IGF1R_UKB_MR_all5_estimators.csv (MAIN result)\n")
cat("  2. 04_heterogeneity.csv\n")
cat("  3. 05_egger_intercept.csv\n")
cat("  4. Console output (especially VERDICT section)\n")
cat("══════════════════════════════════════════════════════════\n")
