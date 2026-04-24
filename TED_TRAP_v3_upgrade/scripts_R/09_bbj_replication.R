#==============================================================================
# TED-TRAP Upgrade — Phase 1 Script 09
# BioBank Japan (BBJ) Graves disease replication (ancestry check)
#==============================================================================
# Purpose:
#   Replicate the TSHR MR finding in a non-European cohort to rule out
#   ancestry-specific artifacts (Sakaue et al. 2021, Nat Genet 53:1415).
#
#   BBJ has published a Graves disease GWAS (N ≈ 180K Japanese, ~1.3K cases)
#   which provides an orthogonal population.
#
# Source:
#   Sakaue et al. (2021). A cross-population atlas of genetic associations
#   for 220 human phenotypes. Nat Genet 53:1415-1424.
#   BBJ portal: https://pheweb.jp/
#   Downloads: https://humandbs.dbcls.jp/en/hum0197-v3
#
# Output:
#   TrackA_MR/results/09_bbj_replication.csv
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(TwoSampleMR)
library(data.table)
library(dplyr)

log_file <- "TrackA_MR/logs/09_bbj.log"
sink(log_file, split = TRUE)
cat("=== BBJ Graves Replication ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- BBJ GWAS file (user downloads) ---
bbj_graves_file <- "TrackA_MR/data/bbj/BBJ_Graves_disease_summary.tsv.gz"

# Alternative: some BBJ phenotypes are on OpenGWAS with bbj-a-* IDs
# ieugwasr::gwasinfo() to search for "Graves"
bbj_opengwas_id <- "bbj-a-123"   # placeholder; replace with correct ID

# Check if file exists, else try OpenGWAS
use_opengwas <- FALSE
if (file.exists(bbj_graves_file)) {
  cat("Using local BBJ file.\n")
  bbj <- fread(bbj_graves_file)
} else {
  cat("Local BBJ file not found; attempting OpenGWAS...\n")
  use_opengwas <- TRUE
  # List BBJ-related Graves datasets
  library(ieugwasr)
  bbj_list <- tryCatch(
    gwasinfo() %>% filter(grepl("Graves", trait, ignore.case = TRUE) &
                           grepl("Japan", population, ignore.case = TRUE)),
    error = function(e) NULL
  )
  cat("Available BBJ-Graves GWAS on OpenGWAS:\n")
  if (!is.null(bbj_list)) print(bbj_list[, c("id", "trait", "sample_size")])

  if (is.null(bbj_list) || nrow(bbj_list) == 0) {
    cat("\n❌ No BBJ Graves data found on OpenGWAS.\n")
    cat("   Download manually from:\n")
    cat("     https://humandbs.dbcls.jp/en/hum0197-v3\n")
    cat("     or https://jenger.riken.jp/en/\n")
    stop("BBJ data not available.")
  }

  bbj_opengwas_id <- bbj_list$id[1]
}

# --- Load TSHR instruments from previous extraction ---
iv_file <- "TrackA_MR/data/instruments/TSHR_instruments_permissive.tsv"
if (!file.exists(iv_file)) {
  iv_file <- "TrackA_MR/data/instruments/TSHR_instruments_strict.tsv"
}
if (!file.exists(iv_file)) {
  stop("TSHR instruments not found. Run 02_eqtlgen_exposure_extraction_all_genes.R first.")
}

exposure_data <- fread(iv_file)

# --- Extract outcome data ---
if (use_opengwas) {
  outcome_data <- extract_outcome_data(snps = exposure_data$SNP,
                                        outcomes = bbj_opengwas_id)
} else {
  # Local file: manual harmonisation
  bbj_sub <- bbj[bbj$SNP %in% exposure_data$SNP, ]
  outcome_data <- format_data(
    dat = bbj_sub,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "AF",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    pval_col = "P",
    samplesize_col = "N"
  )
  outcome_data$id.outcome <- "BBJ_Graves"
  outcome_data$outcome    <- "BBJ Graves disease"
}

cat(sprintf("TSHR IVs available in BBJ: %d / %d\n",
            nrow(outcome_data), nrow(exposure_data)))

# --- Harmonize and run MR ---
h <- harmonise_data(exposure_data, outcome_data, action = 2)
h <- h[h$mr_keep, ]
cat(sprintf("SNPs post-harmonisation: %d\n", nrow(h)))

if (nrow(h) >= 2) {
  mr_res <- mr(h, method_list = "mr_ivw")
  cat("\n--- BBJ Replication MR (TSHR) ---\n")
  print(mr_res)

  fwrite(mr_res, "TrackA_MR/results/09_bbj_replication.csv")
  cat("\n✅ Saved: TrackA_MR/results/09_bbj_replication.csv\n")

  # Interpretation
  pri <- fread("TrackA_MR/results/03_mr_v3_full.csv")
  pri_tshr <- pri[pri$gene == "TSHR" & pri$outcome_role == "Primary" &
                    pri$method == "Inverse variance weighted", ]
  if (nrow(pri_tshr) > 0) {
    direction_match <- sign(pri_tshr$b[1]) == sign(mr_res$b[1])
    cat("\n--- Concordance with European primary MR ---\n")
    cat(sprintf("  Primary (European) β: %.4f\n", pri_tshr$b[1]))
    cat(sprintf("  BBJ (Japanese) β:     %.4f\n", mr_res$b[1]))
    cat(sprintf("  Direction match: %s\n", direction_match))
  }
} else {
  cat("❌ Insufficient SNPs for MR after harmonisation.\n")
}

sink()
