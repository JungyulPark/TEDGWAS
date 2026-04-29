library(TwoSampleMR)
library(data.table)
library(ieugwasr)

set.seed(20260429)

# ------------------------------------------------------------
# Step 1: Load FinnGen R12 GRAVES_STRICT as exposure
# ------------------------------------------------------------
exposure_path <- "c:/ProjectTEDGWAS/finngen_R12_E4_GRAVES_STRICT.gz"
gd_strict <- fread(exposure_path)

cat("Exposure file dimensions:", dim(gd_strict), "\n")
cat("First 3 rows:\n")
print(head(gd_strict, 3))

# Filter for genome-wide significance
gd_strict_sig <- gd_strict[pval < 5e-8]
cat("\nSNPs at P<5e-8 in GRAVES_STRICT:", nrow(gd_strict_sig), "\n")

# Format for TwoSampleMR
exposure_dat <- format_data(
  as.data.frame(gd_strict_sig),
  type = "exposure",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval",
  chr_col = "#chrom",
  pos_col = "pos"
)
exposure_dat$exposure <- "GD_liability_FinnGen_R12_strict"

# ------------------------------------------------------------
# Step 2: LD clumping (same params as MVMR)
# ------------------------------------------------------------
exposure_clumped <- clump_data(
  exposure_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 5e-8,
  pop = "EUR"
)
cat("\nIndependent instruments after clumping:", nrow(exposure_clumped), "\n")

fwrite(exposure_clumped, "c:/ProjectTEDGWAS/TrackA_MR/results/Phase3_GD_liability_instruments.csv")

# ------------------------------------------------------------
# Step 3: Load outcome (FinnGen R12 GRAVES_OPHT) and extract instruments
# ------------------------------------------------------------
outcome_path <- "c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz"
gd_opht <- fread(outcome_path)

outcome_dat_raw <- gd_opht[rsids %in% exposure_clumped$SNP]
cat("\nOverlapping SNPs in outcome:", nrow(outcome_dat_raw), "/", nrow(exposure_clumped), "\n")

outcome_dat <- format_data(
  as.data.frame(outcome_dat_raw),
  type = "outcome",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval",
  chr_col = "#chrom",
  pos_col = "pos"
)
outcome_dat$outcome <- "GravesOphthalmopathy_FinnGen_R12"

# ------------------------------------------------------------
# Step 4: Harmonize
# ------------------------------------------------------------
harmonised <- harmonise_data(
  exposure_dat = exposure_clumped,
  outcome_dat = outcome_dat,
  action = 2  # allele frequency-based palindrome resolution
)
cat("\nHarmonised SNPs (after palindrome filtering):", sum(harmonised$mr_keep), "\n")

fwrite(harmonised, "c:/ProjectTEDGWAS/TrackA_MR/results/Phase3_GD_liability_harmonised.csv")

# ------------------------------------------------------------
# Step 5: Run MR with all five estimators
# ------------------------------------------------------------
mr_results <- mr(harmonised, method_list = c(
  "mr_ivw",
  "mr_weighted_median", 
  "mr_egger_regression",
  "mr_simple_mode",
  "mr_weighted_mode"
))
cat("\n=== MR RESULTS ===\n")
print(mr_results)

# ------------------------------------------------------------
# Step 6: Sensitivity analyses
# ------------------------------------------------------------
het_test    <- mr_heterogeneity(harmonised)
pleio_test  <- mr_pleiotropy_test(harmonised)
steiger     <- directionality_test(harmonised)

cat("\n=== HETEROGENEITY (Cochran Q) ===\n")
print(het_test)

cat("\n=== PLEIOTROPY (MR-Egger intercept) ===\n")
print(pleio_test)

cat("\n=== STEIGER DIRECTIONALITY ===\n")
print(steiger)

# ------------------------------------------------------------
# Check chromosome representation
# ------------------------------------------------------------
cat("\n=== CHROMOSOME DISTRIBUTION OF INSTRUMENTS ===\n")
if("chr.exposure" %in% colnames(harmonised)) {
  print(table(harmonised$chr.exposure))
} else {
  cat("chr.exposure column not found\n")
}

# ------------------------------------------------------------
# Step 7: Summary output (single-row format for table integration)
# ------------------------------------------------------------
ivw <- mr_results[mr_results$method == "Inverse variance weighted", ]
wm  <- mr_results[mr_results$method == "Weighted median", ]
egg <- mr_results[mr_results$method == "MR Egger", ]

summary_row <- data.frame(
  exposure          = "GD_liability_FinnGen_R12_strict",
  outcome           = "GravesOphthalmopathy_FinnGen_R12",
  n_iv              = nrow(harmonised[harmonised$mr_keep, ]),
  ivw_beta          = ivw$b,
  ivw_se            = ivw$se,
  ivw_pval          = ivw$pval,
  wm_beta           = wm$b,
  wm_se             = wm$se,
  wm_pval           = wm$pval,
  egger_beta        = egg$b,
  egger_se          = egg$se,
  egger_pval        = egg$pval,
  egger_intercept   = pleio_test$egger_intercept,
  egger_intercept_p = pleio_test$pval,
  cochran_q         = het_test$Q[het_test$method == "Inverse variance weighted"],
  cochran_q_pval    = het_test$Q_pval[het_test$method == "Inverse variance weighted"],
  steiger_correct   = if(is.null(steiger) || nrow(steiger)==0) NA else steiger$correct_causal_direction,
  steiger_pval      = if(is.null(steiger) || nrow(steiger)==0) NA else steiger$steiger_pval
)

cat("\n=== SUMMARY ROW (FOR TABLE 2C) ===\n")
print(summary_row)

fwrite(summary_row, "c:/ProjectTEDGWAS/TrackA_MR/results/Phase3_GD_liability_univariable_MR_summary.csv")

# ------------------------------------------------------------
# Step 8: SNP-level forest plot data (for Supp if needed)
# ------------------------------------------------------------
single <- mr_singlesnp(harmonised)
fwrite(single, "c:/ProjectTEDGWAS/TrackA_MR/results/Phase3_GD_liability_singlesnp.csv")

cat("\n✅ Phase 3 univariable MR complete. Files saved to TrackA_MR/results/:\n")
