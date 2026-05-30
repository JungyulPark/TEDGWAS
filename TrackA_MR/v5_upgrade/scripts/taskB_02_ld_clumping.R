#!/usr/bin/env Rscript
# =============================================================================
# Task B.2 — LD Clumping Sensitivity for TSHR Locus
# =============================================================================
# Purpose : Assess number of independent SNPs at TSHR locus across multiple
#           P-value × r² × ancestry thresholds. Sensitivity check before
#           committing to GCTA-COJO formal analysis.
# Output  : 02_tshr_finemap/results/TaskB_02_ld_clumping_sensitivity_v1.csv
# Input   : 02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv
# Spec    : WEEK1_TASK_SPEC_v1.1.md, Section 4, Step B.2
#
# LD REFERENCE (verified, see TaskB_00_ld_reference_metadata_v1.csv):
#   1000G_phase3_EUR_hg19 (n=503)
#   1000G_phase3_EAS_hg19 (n=504)
#   Source: 1000 Genomes Project Consortium 2015 Nature 526:68-74, PMID 26432245
#
# IMPORTANT: ld_clump() from ieugwasr uses OpenGWAS API by default, which has
# rate limits. For local execution, set `bfile` and `plink_bin` arguments to
# use local 1000G bfile + local PLINK installation — much faster and rate-free.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
})

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
INPUT_FILE <- file.path(WORK_DIR, "02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv")
OUTPUT_FILE <- file.path(WORK_DIR, "02_tshr_finemap/results/TaskB_02_ld_clumping_sensitivity_v1.csv")
LOG_FILE <- file.path(WORK_DIR, "00_meta/analytical_log.md")

# LD reference bfile prefixes (Gemini to confirm paths)
LD_REFS <- list(
  EUR = list(
    id = "1000G_phase3_EUR_hg19",
    bfile = "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/EUR_chr14"
  ),
  EAS = list(
    id = "1000G_phase3_EAS_hg19",
    bfile = "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3/EAS_chr14"
  )
)

# PLINK binary path (Gemini to confirm)
PLINK_BIN <- "c:/ProjectTEDGWAS/tools/plink.exe"

# Sensitivity grid (per spec v1.1)
P_THRESHOLDS <- c(5e-8, 5e-6)
R2_THRESHOLDS <- c(0.001, 0.01, 0.1)
WINDOW_KB <- 10000

# -----------------------------------------------------------------------------
# Read locus extraction
# -----------------------------------------------------------------------------
locus <- fread(INPUT_FILE)
cat("Read locus extraction:", nrow(locus), "variants\n")

# ieugwasr::ld_clump requires columns: rsid, pval (and optionally id)
clump_input <- locus[, .(rsid = snp, pval = pvalue)]

# -----------------------------------------------------------------------------
# Run clumping across grid
# -----------------------------------------------------------------------------
results <- list()
row_idx <- 1L

for (pthresh in P_THRESHOLDS) {
  for (r2 in R2_THRESHOLDS) {
    for (ancestry in names(LD_REFS)) {

      ref <- LD_REFS[[ancestry]]
      cat(sprintf("\nClumping: P<%g, r2<%g, %s\n", pthresh, r2, ancestry))

      # Pre-filter to P threshold (clumping doesn't always do this internally)
      input_filtered <- clump_input[pval < pthresh]

      if (nrow(input_filtered) == 0) {
        cat("  No SNPs pass P threshold.\n")
        results[[row_idx]] <- data.frame(
          pthresh = pthresh, r2_thresh = r2,
          ld_reference_id = ref$id,
          n_independent_snps = 0L,
          lead_snps = NA_character_,
          includes_rs179252 = FALSE,
          notes = "No SNPs pass P threshold"
        )
        row_idx <- row_idx + 1L
        next
      }

      # Try local bfile first; fall back to API if not available
      clumped <- tryCatch({
        ld_clump(
          dat = input_filtered,
          clump_r2 = r2,
          clump_p = pthresh,
          clump_kb = WINDOW_KB,
          bfile = ref$bfile,
          plink_bin = PLINK_BIN
        )
      }, error = function(e) {
        cat("  Local PLINK failed; trying OpenGWAS API:", conditionMessage(e), "\n")
        pop_code <- if (ancestry == "EUR") "EUR" else "EAS"
        tryCatch({
          ld_clump(dat = input_filtered, clump_r2 = r2, clump_p = pthresh,
                   clump_kb = WINDOW_KB, pop = pop_code)
        }, error = function(e2) {
          cat("  API clumping also failed:", conditionMessage(e2), "\n")
          NULL
        })
      })

      if (is.null(clumped) || nrow(clumped) == 0) {
        n_indep <- NA_integer_
        leads_str <- NA_character_
        has_rs179252 <- NA
        note <- "Clumping failed or empty"
      } else {
        n_indep <- nrow(clumped)
        leads_str <- paste(clumped$rsid, collapse = ",")
        has_rs179252 <- "rs179252" %in% clumped$rsid
        note <- ""
      }

      results[[row_idx]] <- data.frame(
        pthresh = pthresh, r2_thresh = r2,
        ld_reference_id = ref$id,
        n_independent_snps = n_indep,
        lead_snps = leads_str,
        includes_rs179252 = has_rs179252,
        notes = note
      )
      row_idx <- row_idx + 1L
    }
  }
}

results_df <- do.call(rbind, results)

# Sanity check: 12 rows expected (2 P × 3 r² × 2 ancestries)
stopifnot(nrow(results_df) == 12)

# -----------------------------------------------------------------------------
# Write output
# -----------------------------------------------------------------------------
fwrite(results_df, OUTPUT_FILE)
cat("\nWrote:", OUTPUT_FILE, "\n")
print(results_df)

# Log
log_entry <- sprintf(
  "## %s — Task B.2 run\n- Script: taskB_02_ld_clumping.R\n- Grid: 2 P × 3 r2 × 2 ancestries = 12 cells\n- Output: %s\n- rs179252 retained in any clump: %s\n\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  basename(OUTPUT_FILE),
  any(results_df$includes_rs179252, na.rm = TRUE)
)
cat(log_entry, file = LOG_FILE, append = TRUE)
