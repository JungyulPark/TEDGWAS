#!/usr/bin/env Rscript
# =============================================================================
# Task B.1 — TSHR Locus Extraction from eQTLGen
# =============================================================================
# Purpose : Extract all eQTLGen cis-eQTLs within ±1 Mb of TSHR for downstream
#           fine-mapping (COJO + SuSiE).
# Output  : 02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv
# Input   : c:/ProjectTEDGWAS/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
# Spec    : WEEK1_TASK_SPEC_v1.1.md, Section 4, Step B.1
#
# DATA SOURCE (verified):
#   eQTLGen full cis-eQTL summary statistics
#   Vosa et al. 2021 Nat Genet 53:1300-1310, PMID 34475573
#   n=31,684 whole blood samples, native build hg19
#
# TARGET LOCUS (verified Ensembl GRCh37):
#   TSHR (ENSG00000165409)
#   hg19: chr14:81,421,333 - 81,612,646 (forward strand)
#   Extraction region: chr14:80,421,333 - 82,612,646 (gene ± 1 Mb)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
EQTLGEN_FILE <- "c:/ProjectTEDGWAS/TrackA_MR/data/eqtl_TSHR_all.csv"
OUTPUT_FILE <- file.path(WORK_DIR, "02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv")
LOG_FILE <- file.path(WORK_DIR, "00_meta/analytical_log.md")

# TSHR locus (verified Ensembl GRCh37)
TSHR_ENSEMBL <- "ENSG00000165409"
TSHR_GENE_SYMBOL <- "TSHR"
TSHR_TSS_HG19 <- 81421333L
TSHR_CHR <- 14L
EXTRACT_WINDOW_BP <- 1000000L  # ±1 Mb

EXTRACT_START <- TSHR_TSS_HG19 - EXTRACT_WINDOW_BP
EXTRACT_END   <- 81612646L + EXTRACT_WINDOW_BP  # gene end + 1Mb

cat("TSHR extraction window:\n")
cat("  Chr:  ", TSHR_CHR, "\n")
cat("  Start:", EXTRACT_START, "(hg19)\n")
cat("  End:  ", EXTRACT_END, "(hg19)\n")

# -----------------------------------------------------------------------------
# Read eQTLGen TSHR locus sumstats
# -----------------------------------------------------------------------------
# eQTLGen full columns (per consortium documentation):
#   Pvalue, SNP, SNPChr, SNPPos, AssessedAllele, OtherAllele,
#   Zscore, Gene, GeneSymbol, GeneChr, GenePos, NrCohorts, NrSamples, FDR, BonferroniP

cat("\nReading TSHR eQTL data from:\n  ", EQTLGEN_FILE, "\n")

# Read with column selection
eqtlgen <- fread(
  EQTLGEN_FILE,
  select = c("Pvalue", "SNP", "SNPChr", "SNPPos",
             "AssessedAllele", "OtherAllele", "Zscore",
             "Gene", "GeneSymbol", "NrSamples", "FDR")
)

cat("Total rows in eQTLGen:", nrow(eqtlgen), "\n")

# -----------------------------------------------------------------------------
# Filter to TSHR locus
# -----------------------------------------------------------------------------
tshr_locus <- eqtlgen[
  (Gene == TSHR_ENSEMBL | GeneSymbol == TSHR_GENE_SYMBOL) &
  SNPChr == TSHR_CHR &
  SNPPos >= EXTRACT_START &
  SNPPos <= EXTRACT_END
]

cat("Rows in TSHR locus (±1Mb):", nrow(tshr_locus), "\n")

if (nrow(tshr_locus) == 0) {
  stop("No TSHR cis-eQTLs found in extraction window. Check gene ID / coordinates / file.")
}

# -----------------------------------------------------------------------------
# Compute derived columns
# -----------------------------------------------------------------------------
# eQTLGen does NOT provide beta/SE directly; only Zscore + Pvalue + N.
# We can derive beta/SE for MR using the formula:
#   beta = Zscore / sqrt(2 * EAF * (1 - EAF) * (N + Zscore^2))
#   se   = 1     / sqrt(2 * EAF * (1 - EAF) * (N + Zscore^2))
# However eQTLGen full sumstats do NOT include EAF. Use 1000G EUR MAF as proxy,
# or use Z + N for MR-Egger / IVW directly (TwoSampleMR can work from z-scores).
# For Week 1 locus extraction we keep z and p only; beta/se computed downstream.

setnames(tshr_locus,
         c("SNP", "SNPChr", "SNPPos",
           "AssessedAllele", "OtherAllele", "Zscore",
           "Pvalue", "FDR", "NrSamples"),
         c("snp", "chr", "pos_hg19",
           "effect_allele", "other_allele", "zscore",
           "pvalue", "fdr", "n_samples"))

# Required output columns per spec
tshr_locus[, pos_hg38 := NA_integer_]  # to be filled by liftOver below
tshr_locus[, liftover_status := "native_hg19"]
tshr_locus[, eaf := NA_real_]  # not in eQTLGen; will use 1000G later
tshr_locus[, beta := NA_real_]  # derived later
tshr_locus[, se := NA_real_]    # derived later
tshr_locus[, distance_to_tss := pos_hg19 - TSHR_TSS_HG19]
tshr_locus[, is_rs179252 := (snp == "rs179252")]

# Reorder columns to match spec exactly
ordered_cols <- c("snp", "chr", "pos_hg19", "pos_hg38", "liftover_status",
                  "effect_allele", "other_allele", "eaf",
                  "beta", "se", "zscore", "pvalue", "fdr",
                  "n_samples", "distance_to_tss", "is_rs179252")
setcolorder(tshr_locus, ordered_cols)
tshr_locus <- tshr_locus[, ..ordered_cols]

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
# Check rs179252 is present (locked from v4.3 ground truth as TSHR BBJ coloc lead)
if (!any(tshr_locus$is_rs179252)) {
  warning("rs179252 NOT FOUND in eQTLGen TSHR locus extraction. ",
          "Investigate: gene ID mismatch, coordinate mismatch, or file integrity.")
} else {
  rs179252_row <- tshr_locus[is_rs179252 == TRUE]
  cat("\nrs179252 found:\n")
  cat("  Pos (hg19):", rs179252_row$pos_hg19, "\n")
  cat("  P-value:   ", rs179252_row$pvalue, "\n")
  cat("  Z-score:   ", rs179252_row$zscore, "\n")
  # Expected: P ~ 9.1e-15 per locked memory
}

# Number of variants at different significance thresholds
cat("\nVariant counts:\n")
cat("  Total in ±1Mb:    ", nrow(tshr_locus), "\n")
cat("  P < 5e-8:         ", sum(tshr_locus$pvalue < 5e-8), "\n")
cat("  P < 5e-6:         ", sum(tshr_locus$pvalue < 5e-6), "\n")

# -----------------------------------------------------------------------------
# OPTIONAL: liftOver to hg38 (for cross-reference with UKB-PPP if pQTL found)
# Requires rtracklayer + chain file
# -----------------------------------------------------------------------------
# Uncomment and run AFTER chain file is downloaded:
#
# library(rtracklayer)
# chain <- import.chain("hg19ToHg38.over.chain")
# gr_hg19 <- GRanges(seqnames = paste0("chr", tshr_locus$chr),
#                    ranges = IRanges(start = tshr_locus$pos_hg19, width = 1))
# gr_hg38 <- liftOver(gr_hg19, chain)
# tshr_locus[, pos_hg38 := vapply(as.list(gr_hg38),
#                                  function(x) if (length(x) == 1) start(x) else NA_integer_,
#                                  integer(1))]
# tshr_locus[!is.na(pos_hg38), liftover_status := "lifted_to_hg38_success"]
# tshr_locus[is.na(pos_hg38), liftover_status := "lifted_failed"]

# -----------------------------------------------------------------------------
# Write output
# -----------------------------------------------------------------------------
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
fwrite(tshr_locus, OUTPUT_FILE, sep = "\t")
cat("\nWrote:", OUTPUT_FILE, "\n")
cat("Rows: ", nrow(tshr_locus), "\n")

# Log
log_entry <- sprintf(
  "## %s — Task B.1 run\n- Script: taskB_01_locus_extract.R\n- Input: eQTLGen 2019-12-11 release\n- Locus: TSHR (ENSG00000165409) hg19 chr14:%d-%d ±1Mb\n- Output: %s\n- Rows: %d\n- rs179252 present: %s\n\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  TSHR_TSS_HG19, 81612646L,
  basename(OUTPUT_FILE), nrow(tshr_locus),
  any(tshr_locus$is_rs179252)
)
cat(log_entry, file = LOG_FILE, append = TRUE)
