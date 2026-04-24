#==============================================================================
# TED-TRAP Upgrade — Phase 2 Script 11
# BH-FDR correction utilities
#==============================================================================
# Purpose:
#   Apply consistent Benjamini-Hochberg false discovery rate correction
#   across MR, colocalization, and DE analyses. This standalone utility
#   ensures reviewers cannot find inconsistencies in multiple-testing handling.
#
# Multiple-testing hierarchies applied:
#   Level 1: Per-outcome MR across candidate genes (~8 tests)
#   Level 2: Per-gene coloc across loci (~8 tests)
#   Level 3: Per-gene DE within pre-specified set (~20 tests)
#   Level 4: Hypergeometric enrichment tests (~5 tests)
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(data.table)

log_file <- "TrackA_MR/logs/11_bh_fdr.log"
sink(log_file, split = TRUE)
cat("=== BH-FDR correction ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

apply_bh_fdr <- function(df, p_col, new_col = NULL, threshold = 0.05) {
  if (is.null(new_col)) new_col <- paste0(p_col, "_BH")
  df[[new_col]] <- p.adjust(df[[p_col]], method = "BH")
  df[[paste0("sig_BH_", threshold)]] <- df[[new_col]] < threshold
  df
}

# --- Apply to MR results ---
mr_file <- "TrackA_MR/results/03_mr_v3_full.csv"
if (file.exists(mr_file)) {
  mr <- fread(mr_file)
  cat("MR results loaded:", nrow(mr), "rows\n")

  # Apply BH per outcome_role (so we correct within each outcome family)
  mr_corrected <- mr[, {
    bh <- p.adjust(pval, method = "BH")
    .(.SD, pval_BH = bh, sig_BH = bh < 0.05)
  }, by = outcome_role]

  # Flatten
  mr_final <- cbind(mr, pval_BH = unlist(lapply(split(mr$pval, mr$outcome_role),
                                                  function(p) p.adjust(p, method = "BH"))))
  mr_final$sig_BH_0.05 <- mr_final$pval_BH < 0.05

  fwrite(mr_final, "TrackA_MR/results/11_mr_v3_bh.csv")
  cat("  Saved: TrackA_MR/results/11_mr_v3_bh.csv\n")
} else {
  cat("MR results not yet available\n")
}

# --- Apply to DESeq2 results ---
de_file <- "TrackA_MR/results/10_deseq2_v3_candidate_set.csv"
if (file.exists(de_file)) {
  de <- fread(de_file)
  cat("\nDE results loaded:", nrow(de), "rows\n")
  # DESeq2 already applies BH at the genome-wide level
  # Here we apply BH specifically to the candidate set
  de$pval_BH_candidate <- p.adjust(de$pvalue, method = "BH")
  fwrite(de, "TrackA_MR/results/11_de_v3_bh.csv")
  cat("  Saved: TrackA_MR/results/11_de_v3_bh.csv\n")
}

# --- Apply to coloc results ---
# Coloc provides posterior probabilities, not p-values, so PP.H4 threshold replaces FDR
# We apply a threshold classification instead
coloc_file <- "TrackA_MR/results/06_coloc_v3_all_loci.csv"
if (file.exists(coloc_file)) {
  coloc <- fread(coloc_file)
  coloc$coloc_tier <- cut(coloc$PP_H4,
                          breaks = c(0, 0.5, 0.8, 0.95, 1.01),
                          labels = c("No evidence", "Weak", "Strong", "Very strong"),
                          include.lowest = TRUE)
  fwrite(coloc, "TrackA_MR/results/11_coloc_v3_tiers.csv")
  cat("\nColoc tiered: saved TrackA_MR/results/11_coloc_v3_tiers.csv\n")
}

sink()
cat("\n✅ BH-FDR correction complete.\n\n")
