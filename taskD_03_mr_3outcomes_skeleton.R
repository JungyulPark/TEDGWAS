#!/usr/bin/env Rscript
# =============================================================================
# Task D.3 (SKELETON) — Druggable-Genome-Wide MR across 3 outcomes
# =============================================================================
# Reads SNP-level instruments from D.2a and runs MR against:
#   1. BBJ Graves disease       (discovery,  EAS, 2,809 cases)
#   2. UKB hyperthyroidism      (replication, EUR, 3,731 cases)
#   3. FinnGen R12 Graves ophth (TED-specific sensitivity, 858 cases, local .gz)
#
# Output: 04_druggable_mr/results/TaskD_03_mr_results_3outcomes_v1.csv
#
# IMPORTANT — this is a SKELETON. Before running:
#   (a) Confirm the THREE outcome GWAS file paths / IEU IDs below.
#   (b) Bonferroni denominator is READ from D.2 output (never hardcoded/simulated).
#   (c) Cross-ancestry note: eQTLGen instruments are EUR; BBJ outcome is EAS.
#       Harmonize on allele; flag palindromic SNPs; document as limitation.
#   (d) Multi-IV genes: IVW + MR-Egger + WM + Wmode. Single-IV: Wald only.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)   # Hemani 2018 eLife 7:e34408 (PMID 29846171)
})

WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
SNP_FILE <- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv")
GENE_FILE<- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02b_eqtl_instrument_gene_summary_v1.csv")
DENOM_FILE<-file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02_denominator_note.txt")
OUT_FILE <- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_03_mr_results_3outcomes_v1.csv")
LOG_FILE <- file.path(WORK_DIR, "00_meta/analytical_log.md")

# --- Outcome sources (CONFIRM THESE before running) -------------------------
# Sakaue 2021 Nat Genet 53:1415-1424 (PMID 34594039) for BBJ + UKB.
# FinnGen R12 GO local .gz.
OUTCOMES <- list(
  BBJ_Graves        = list(type="file_or_ieu", path_or_id="<<CONFIRM: BBJ Graves disease sumstats>>", ancestry="EAS", n_case=2809),
  UKB_hyperthyroid  = list(type="file_or_ieu", path_or_id="<<CONFIRM: UKB hyperthyroidism sumstats>>", ancestry="EUR", n_case=3731),
  FinnGen_GO        = list(type="file",        path_or_id="c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz", ancestry="EUR", n_case=858)
)

# -----------------------------------------------------------------------------
# Read instruments + denominator (denominator NEVER hardcoded)
# -----------------------------------------------------------------------------
snp_level  <- fread(SNP_FILE)
gene_level <- fread(GENE_FILE)

denom_lines <- readLines(DENOM_FILE)
denom_val   <- as.integer(sub(".*= ", "", grep("^Bonferroni denominator = ", denom_lines, value=TRUE)))
if (is.na(denom_val) || denom_val < 1) stop("Could not read Bonferroni denominator from D.2 output. STOP.")
BONF_THRESHOLD <- 0.05 / denom_val
cat(sprintf("Bonferroni denominator (from D.2): %d -> threshold %.3e\n", denom_val, BONF_THRESHOLD))

genes_to_test <- gene_level[has_valid_instrument == TRUE, gene_symbol]
cat("Genes to test:", length(genes_to_test), "\n")

# -----------------------------------------------------------------------------
# Derive beta/se for exposure from Z + N (eQTLGen provides Z, P, N)
#   se_hat = 1 / sqrt(2*p*(1-p)*(N + Z^2));  beta_hat = Z * se_hat
#   where p = effect allele frequency. Use 1000G EUR MAF as proxy (consistent
#   with Week 1 .ma construction). MAF must be merged in from EUR_chr14_freq or
#   a genome-wide EUR .frq. <<CONFIRM eaf source merge>>
# -----------------------------------------------------------------------------
# PLACEHOLDER: merge EAF then compute. Gemini to wire EAF source.
# snp_level[, se := 1/sqrt(2*eaf*(1-eaf)*(n_samples + zscore^2))]
# snp_level[, beta := zscore * se]

# -----------------------------------------------------------------------------
# MR loop (skeleton) — per gene x per outcome
# -----------------------------------------------------------------------------
run_one_gene_outcome <- function(gene_sym, exposure_dt, outcome_dt) {
  # harmonise_data() expects standard columns; build exposure & outcome frames.
  # Single-IV -> mr() returns Wald ratio automatically.
  # Multi-IV  -> request IVW, Egger, WM, Wmode.
  # Return a tidy data.table row per method.
  # ... (Gemini fills harmonisation; structure below) ...
  NULL
}

results <- list()
for (g in genes_to_test) {
  exp_g <- snp_level[gene_symbol == g]
  n_iv  <- nrow(exp_g)
  for (oc_name in names(OUTCOMES)) {
    # res <- run_one_gene_outcome(g, exp_g, outcome_loaded[[oc_name]])
    # results[[length(results)+1]] <- res
  }
}

# out <- rbindlist(results)
# fwrite(out, OUT_FILE)

# -----------------------------------------------------------------------------
# Output schema (per spec Section 5, D.3):
#   gene_symbol, ensembl_gene_id, outcome, method, n_iv, beta, se, pvalue,
#   or, or_ci_lower, or_ci_upper, egger_intercept, egger_intercept_p,
#   cochran_q, cochran_q_p, steiger_direction, steiger_p, notes
# -----------------------------------------------------------------------------

cat("D.3 skeleton ready. Gemini to: (1) confirm 3 outcome paths, (2) wire EAF source\n",
    "for beta/se, (3) implement harmonise + mr per gene/outcome, (4) write output.\n")
cat("Bonferroni threshold to apply in D.4:", BONF_THRESHOLD, "\n")
