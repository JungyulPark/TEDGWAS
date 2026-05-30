#!/usr/bin/env Rscript
# =============================================================================
# Task B.3 Helper — Convert eQTLGen TSHR locus to COJO .ma format
# =============================================================================
# Usage  : Rscript taskB_03_helper_make_ma.R <input.tsv> <output.ma>
# Purpose: eQTLGen provides Z-score + N but not beta/SE/EAF directly.
#          COJO requires .ma format: SNP A1 A2 freq b se p N
#          We derive beta and SE using 1000G EUR MAF as EAF proxy.
#
# Formula (per Zhu et al. 2016 SMR derivation):
#   se   = 1 / sqrt(2 * eaf * (1 - eaf) * (N + Z^2))
#   beta = Z * se
#
# CAVEAT: Using EUR MAF as proxy is standard for European-focused MR but is
# imperfect for cross-ancestry analysis. For EAS COJO, ideally use BBJ/EAS
# allele frequencies if available.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript taskB_03_helper_make_ma.R <input.tsv> <output.ma>")
}

INPUT_FILE <- args[1]
OUTPUT_FILE <- args[2]

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Read locus extraction
# -----------------------------------------------------------------------------
locus <- fread(INPUT_FILE)
cat("Read", nrow(locus), "variants from", INPUT_FILE, "\n")

# -----------------------------------------------------------------------------
# Get EAF from 1000G EUR
# -----------------------------------------------------------------------------
# OPTION A — Read pre-computed allele frequencies from 1000G EUR .frq file
#            (generated separately with: plink --bfile EUR_chr14 --freq --out EUR_chr14_freq)
EUR_FRQ_FILE <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/EUR_chr14_freq.frq"

if (file.exists(EUR_FRQ_FILE)) {
  frq <- fread(EUR_FRQ_FILE)
  # PLINK .frq columns: CHR SNP A1 A2 MAF NCHROBS
  # MAF is for A1; need to flip if A1 != effect_allele
  setnames(frq, c("SNP", "A1", "A2", "MAF"), c("snp", "frq_a1", "frq_a2", "maf"))

  merged <- merge(locus, frq[, .(snp, frq_a1, frq_a2, maf)], by = "snp", all.x = TRUE)

  # Align EAF to effect_allele
  merged[, eaf := fifelse(effect_allele == frq_a1, maf,
                  fifelse(effect_allele == frq_a2, 1 - maf, NA_real_))]
  matched <- !is.na(merged$eaf)
  cat("Matched", sum(matched), "/", nrow(merged), "SNPs to 1000G EUR allele frequencies\n")
} else {
  warning("1000G EUR .frq file not found at: ", EUR_FRQ_FILE,
          "\nUsing EAF=0.30 as conservative default (may bias beta/SE).")
  merged <- locus
  merged[, eaf := 0.30]  # conservative default
}

# Drop SNPs without valid EAF
merged <- merged[!is.na(eaf) & eaf > 0 & eaf < 1]

# -----------------------------------------------------------------------------
# Derive beta and SE
# -----------------------------------------------------------------------------
# se = 1 / sqrt(2 * eaf * (1-eaf) * (N + Z^2))
# beta = Z * se
merged[, se := 1 / sqrt(2 * eaf * (1 - eaf) * (n_samples + zscore^2))]
merged[, beta := zscore * se]

# -----------------------------------------------------------------------------
# Write .ma file (whitespace-delimited)
# -----------------------------------------------------------------------------
ma <- merged[, .(SNP = snp,
                 A1 = effect_allele,
                 A2 = other_allele,
                 freq = eaf,
                 b = beta,
                 se = se,
                 p = pvalue,
                 N = n_samples)]

fwrite(ma, OUTPUT_FILE, sep = " ", quote = FALSE)
cat("Wrote .ma file:", OUTPUT_FILE, "(", nrow(ma), "rows)\n")
