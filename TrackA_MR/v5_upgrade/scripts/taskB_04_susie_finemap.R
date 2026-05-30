#!/usr/bin/env Rscript
# =============================================================================
# Task B.4 — SuSiE Fine-mapping at TSHR Locus
# =============================================================================
# Purpose : Bayesian fine-mapping at TSHR locus using SuSiE-RSS (summary stats
#           + LD matrix). Identifies credible sets of variants.
# Output  : 02_tshr_finemap/results/TaskB_04_tshr_susie_credible_sets_v1.tsv
# Input   : TaskB_01 locus extraction + 1000G EUR/EAS bfiles for LD matrix
# Spec    : WEEK1_TASK_SPEC_v1.1.md, Section 4, Step B.4
#
# METHOD REFERENCES (verified):
#   Wang et al. 2020 J R Stat Soc B (SuSiE original)
#   Zou et al. 2022 PLOS Genet (susie_rss for summary stats)
#   susieR package: https://stephenslab.github.io/susieR/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
  library(Matrix)
})

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
INPUT_FILE <- file.path(WORK_DIR, "02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv")
OUTPUT_FILE <- file.path(WORK_DIR, "02_tshr_finemap/results/TaskB_04_tshr_susie_credible_sets_v1.tsv")
LOG_FILE <- file.path(WORK_DIR, "00_meta/analytical_log.md")

# LD reference bfile prefixes
BFILE_EUR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/EUR_chr14"
BFILE_EAS <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3/EAS_chr14"
PLINK_BIN <- "c:/ProjectTEDGWAS/tools/plink.exe"

# SuSiE parameters (per spec v1.1)
SUSIE_L <- 10            # max number of credible sets to discover
SUSIE_COVERAGE <- 0.95   # posterior coverage for each CS
EQTLGEN_N <- 31684L      # eQTLGen sample size

# -----------------------------------------------------------------------------
# Read locus
# -----------------------------------------------------------------------------
locus <- fread(INPUT_FILE)
cat("Read", nrow(locus), "variants for SuSiE\n")

# Filter to variants with valid z-scores and reasonable signal density
locus <- locus[!is.na(zscore) & !is.na(pos_hg19)]
cat("After filtering:", nrow(locus), "variants\n")

# -----------------------------------------------------------------------------
# Helper: compute LD matrix from local PLINK bfile
# -----------------------------------------------------------------------------
# Strategy: write SNP list, run PLINK --r square, read result.
# This avoids OpenGWAS API rate limits and gives full control over reference.

compute_ld_matrix <- function(snps, bfile_prefix, tmp_prefix) {
  snp_file <- paste0(tmp_prefix, ".snps")
  writeLines(snps, snp_file)

  cmd <- sprintf(
    '"%s" --bfile "%s" --extract "%s" --r square --keep-allele-order --out "%s"',
    PLINK_BIN, bfile_prefix, snp_file, tmp_prefix
  )
  cat("Running:", cmd, "\n")
  system(cmd)

  ld_file <- paste0(tmp_prefix, ".ld")
  bim_file <- paste0(tmp_prefix, ".bim")  # might not be created with --r; check
  # PLINK output for --r square is space-separated correlation matrix
  if (!file.exists(ld_file)) {
    stop("PLINK LD matrix not generated: ", ld_file)
  }

  ld <- as.matrix(fread(ld_file))
  # The SNPs in the output are the *intersection* of input list and bfile;
  # need to read the .bim or .nosex to know the order.
  # Use the temp prefix .nosex or .extract list — PLINK preserves bfile order.
  # Easier: re-run with --write-snplist and read that.
  snplist_file <- paste0(tmp_prefix, ".snplist")
  if (!file.exists(snplist_file)) {
    # Fall back: run PLINK with --write-snplist
    system(sprintf('"%s" --bfile "%s" --extract "%s" --write-snplist --out "%s"',
                   PLINK_BIN, bfile_prefix, snp_file, tmp_prefix))
  }
  retained_snps <- readLines(snplist_file)

  rownames(ld) <- retained_snps
  colnames(ld) <- retained_snps

  list(ld = ld, snps = retained_snps)
}

# -----------------------------------------------------------------------------
# Run SuSiE for one ancestry
# -----------------------------------------------------------------------------
run_susie_for_ancestry <- function(locus, bfile_prefix, ld_id, tmp_prefix) {
  cat("\n=== SuSiE with", ld_id, "===\n")

  ld_result <- compute_ld_matrix(locus$snp, bfile_prefix, tmp_prefix)
  ld_matrix <- ld_result$ld
  retained <- ld_result$snps

  # Subset locus to SNPs present in LD reference
  sub <- locus[snp %in% retained]
  sub <- sub[match(retained, snp)]  # reorder to match LD matrix
  cat("Variants retained after LD intersection:", nrow(sub), "\n")

  if (nrow(sub) < 10) {
    warning("Too few variants for stable SuSiE; skipping.")
    return(NULL)
  }

  # Run SuSiE-RSS
  fit <- tryCatch({
    susie_rss(
      z = sub$zscore,
      R = ld_matrix,
      n = EQTLGEN_N,
      L = SUSIE_L,
      coverage = SUSIE_COVERAGE
    )
  }, error = function(e) {
    cat("SuSiE failed:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(fit)) return(NULL)

  # Extract credible sets
  cs_list <- fit$sets$cs
  if (is.null(cs_list) || length(cs_list) == 0) {
    cat("No credible sets identified.\n")
    return(data.frame(
      ld_reference_id = ld_id,
      cs_id = NA_integer_,
      n_snps_in_cs = NA_integer_,
      cs_purity = NA_real_,
      top_snp = NA_character_,
      top_snp_pip = NA_real_,
      top_snp_pos_hg19 = NA_integer_,
      top_snp_p = NA_real_,
      ld_top_with_rs179252_eur = NA_real_,
      ld_top_with_rs179252_eas = NA_real_,
      coverage = NA_real_,
      converged = fit$converged
    ))
  }

  pips <- susie_get_pip(fit)
  out_rows <- list()
  for (i in seq_along(cs_list)) {
    cs_idx <- cs_list[[i]]
    top_idx <- cs_idx[which.max(pips[cs_idx])]
    top_snp <- sub$snp[top_idx]

    # cs purity from SuSiE output
    purity <- if (!is.null(fit$sets$purity)) fit$sets$purity[i, "min.abs.corr"] else NA_real_

    out_rows[[i]] <- data.frame(
      ld_reference_id = ld_id,
      cs_id = i,
      n_snps_in_cs = length(cs_idx),
      cs_purity = purity,
      top_snp = top_snp,
      top_snp_pip = pips[top_idx],
      top_snp_pos_hg19 = sub$pos_hg19[top_idx],
      top_snp_p = sub$pvalue[top_idx],
      ld_top_with_rs179252_eur = NA_real_,   # filled in below
      ld_top_with_rs179252_eas = NA_real_,
      coverage = SUSIE_COVERAGE,
      converged = fit$converged
    )
  }

  do.call(rbind, out_rows)
}

# -----------------------------------------------------------------------------
# Execute for EUR and EAS
# -----------------------------------------------------------------------------
tmp_dir <- file.path(WORK_DIR, "02_tshr_finemap/data/susie_tmp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

results_eur <- run_susie_for_ancestry(locus, BFILE_EUR, "1000G_phase3_EUR_hg19",
                                       file.path(tmp_dir, "eur"))
results_eas <- run_susie_for_ancestry(locus, BFILE_EAS, "1000G_phase3_EAS_hg19",
                                       file.path(tmp_dir, "eas"))

combined <- rbind(results_eur, results_eas)

# -----------------------------------------------------------------------------
# Fill in LD with rs179252 from local bfiles
# -----------------------------------------------------------------------------
# Compute LD locally via PLINK rather than relying on online OpenGWAS APIs
get_r2_with_rs179252 <- function(target_snp, pop) {
  if (is.na(target_snp)) return(NA_real_)
  if (target_snp == "rs179252") return(1.0)
  
  bfile <- if (pop == "EUR") {
    "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/EUR_chr14"
  } else {
    "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3/EAS_chr14"
  }
  
  temp_extract <- tempfile(fileext = ".txt")
  writeLines(c("rs179252", target_snp), temp_extract)
  
  temp_out <- tempfile()
  cmd <- sprintf('"%s" --bfile "%s" --r2 --ld-snp rs179252 --extract "%s" --ld-window-r2 0 --ld-window 99999 --ld-window-kb 99999 --out "%s"',
                 PLINK_BIN, bfile, temp_extract, temp_out)
  
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  ld_file <- paste0(temp_out, ".ld")
  r2_val <- NA_real_
  if (file.exists(ld_file)) {
    ld_data <- fread(ld_file)
    target_row <- ld_data[SNP_B == target_snp]
    if (nrow(target_row) > 0) {
      r2_val <- target_row$R2[1]
    }
    unlink(c(ld_file, paste0(temp_out, ".log")))
  }
  unlink(temp_extract)
  return(r2_val)
}

combined$ld_top_with_rs179252_eur <- vapply(combined$top_snp,
                                            get_r2_with_rs179252,
                                            numeric(1), pop = "EUR")
combined$ld_top_with_rs179252_eas <- vapply(combined$top_snp,
                                            get_r2_with_rs179252,
                                            numeric(1), pop = "EAS")


# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
fwrite(combined, OUTPUT_FILE, sep = "\t")
cat("\nWrote:", OUTPUT_FILE, "\n")
print(combined)

log_entry <- sprintf(
  "## %s — Task B.4 run\n- Script: taskB_04_susie_finemap.R\n- Method: SuSiE-RSS (Wang 2020 / Zou 2022)\n- L=%d, coverage=%.2f, N=%d\n- Output: %s\n- Total credible sets (EUR+EAS): %d\n\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  SUSIE_L, SUSIE_COVERAGE, EQTLGEN_N,
  basename(OUTPUT_FILE),
  sum(!is.na(combined$cs_id))
)
cat(log_entry, file = LOG_FILE, append = TRUE)
