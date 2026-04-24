#==============================================================================
# TED-TRAP Upgrade — Phase 1 Script 07
# MRlap: Sample-overlap correction between eQTLGen and outcome GWAS
#==============================================================================
# Purpose:
#   Resolve MAJOR issue #9: Sample overlap between eQTLGen (which contains
#   UK Biobank participants) and the Pan-UKB hyperthyroidism GWAS can bias
#   standard IVW MR.
#
#   MRlap (Mounier & Kutalik 2023, Bioinformatics 39:btad038) computes an
#   overlap-corrected MR estimate via LDSC-based phenotypic/genetic correlation.
#
# Required inputs:
#   - eQTLGen full summary stats (exposure)
#   - Outcome GWAS full summary stats (the harder to obtain)
#   - LDSC reference: eur_w_ld_chr (from Bulik-Sullivan lab)
#
# Output:
#   TrackA_MR/results/07_mrlap_corrected.csv
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(MRlap)
library(data.table)

log_file <- "TrackA_MR/logs/07_mrlap.log"
sink(log_file, split = TRUE)
cat("=== MRlap — Sample Overlap Correction ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- Required: LDSC reference files ---
ld_path  <- "TrackA_MR/data/ldsc_ref/eur_w_ld_chr/"
hm3_file <- "TrackA_MR/data/ldsc_ref/w_hm3.snplist"

if (!dir.exists(ld_path)) {
  cat("❌ LDSC reference missing. Download from:\n")
  cat("   https://ibg.colorado.edu/mtag/ (eur_w_ld_chr.tar.bz2)\n")
  cat("   Extract to: TrackA_MR/data/ldsc_ref/eur_w_ld_chr/\n")
  stop("Download LDSC reference first.")
}

# --- Summary statistics paths ---
# Format required by MRlap:
#   columns: SNP, A1, A2, N, Z (or beta/SE + N)
exposure_ss <- "TrackA_MR/data/eqtlgen/eqtlgen_TSHR_fullss.tsv"   # prepared separately
outcome_ss  <- "TrackA_MR/data/ukb_hyperthyroid/ebi-a-GCST90038636_fullss.tsv"

if (!file.exists(exposure_ss) || !file.exists(outcome_ss)) {
  cat("❌ Summary statistics files missing:\n")
  cat(sprintf("   Exposure: %s\n", exposure_ss))
  cat(sprintf("   Outcome:  %s\n", outcome_ss))
  cat("\nTo prepare:\n")
  cat("  1. Download eQTLGen cis-eQTL full SS for TSHR region\n")
  cat("  2. Download target GWAS full SS (ebi-a-GCST90038636)\n")
  cat("  3. Harmonize to MRlap columns: SNP, A1, A2, N, Z (or BETA/SE)\n")
  stop("Prepare full summary statistics to proceed.")
}

# --- Run MRlap ---
cat("Running MRlap for TSHR → Hyperthyroidism...\n")
result <- tryCatch({
  MRlap(
    exposure = exposure_ss,
    exposure_name = "TSHR_eQTLGen",
    outcome = outcome_ss,
    outcome_name = "Hyperthyroidism_ebiaGCST90038636",
    ld = ld_path,
    hm3 = hm3_file,
    MR_threshold = 5e-8,
    MR_pruning_dist = 500,
    MR_pruning_LD = 0.001,
    verbose = TRUE
  )
}, error = function(e) {
  cat(sprintf("[ERROR] MRlap: %s\n", e$message))
  NULL
})

if (!is.null(result)) {
  cat("\n--- MRlap Results ---\n")
  cat(sprintf("  Raw IVW α (uncorrected): %.4f (SE=%.4f, p=%.2e)\n",
              result$MRcorrection$observed_effect,
              result$MRcorrection$observed_effect_se,
              result$MRcorrection$observed_effect_p))
  cat(sprintf("  Corrected α (overlap-adjusted): %.4f (SE=%.4f, p=%.2e)\n",
              result$MRcorrection$corrected_effect,
              result$MRcorrection$corrected_effect_se,
              result$MRcorrection$corrected_effect_p))
  cat(sprintf("  Sample overlap fraction estimate: %.3f\n",
              result$GeneticArchitecture$sample_overlap_estimate %||% NA))
  cat(sprintf("  LDSC genetic correlation (rg): %.3f\n",
              result$GeneticArchitecture$rg %||% NA))

  # Save
  out_df <- data.frame(
    exposure        = "TSHR_eQTLGen",
    outcome         = "Hyperthyroidism",
    raw_alpha       = result$MRcorrection$observed_effect,
    raw_alpha_se    = result$MRcorrection$observed_effect_se,
    raw_alpha_p     = result$MRcorrection$observed_effect_p,
    corrected_alpha = result$MRcorrection$corrected_effect,
    corrected_alpha_se = result$MRcorrection$corrected_effect_se,
    corrected_alpha_p  = result$MRcorrection$corrected_effect_p,
    stringsAsFactors = FALSE
  )
  fwrite(out_df, "TrackA_MR/results/07_mrlap_corrected.csv")

  cat("\n✅ Saved: TrackA_MR/results/07_mrlap_corrected.csv\n")
}

sink()
