#!/usr/bin/env Rscript
# =============================================================================
# Task B.5 — Integration of COJO + SuSiE Results into Final TSHR Fine-map Summary
# =============================================================================
# Purpose : Reconcile COJO (frequentist conditional) and SuSiE (Bayesian
#           fine-mapping) results to produce final list of independent signals
#           at TSHR locus.
# Output  : 02_tshr_finemap/results/TaskB_05_tshr_finemap_summary_v1.csv
# Inputs  : TaskB_03_tshr_cojo_independent_signals_v1.tsv
#           TaskB_04_tshr_susie_credible_sets_v1.tsv
#           TaskB_01_tshr_locus_extraction_v1.tsv (for marginal P lookup)
# Spec    : WEEK1_TASK_SPEC_v1.1.md, Section 4, Step B.5
#
# DECISION CRITERIA:
#   - "usable_as_independent_iv" = TRUE iff:
#       (1) Conditional P < 5e-6 (COJO), AND
#       (2) r² with rs179252 < 0.1 in BOTH EUR and EAS, AND
#       (3) Signal appears in SuSiE credible set with PIP > 0.5
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
RES_DIR <- file.path(WORK_DIR, "02_tshr_finemap/results")

LOCUS_FILE <- file.path(RES_DIR, "TaskB_01_tshr_locus_extraction_v1.tsv")
COJO_FILE  <- file.path(RES_DIR, "TaskB_03_tshr_cojo_independent_signals_v1.tsv")
SUSIE_FILE <- file.path(RES_DIR, "TaskB_04_tshr_susie_credible_sets_v1.tsv")
OUTPUT_FILE <- file.path(RES_DIR, "TaskB_05_tshr_finemap_summary_v1.csv")
LOG_FILE <- file.path(WORK_DIR, "00_meta/analytical_log.md")

# -----------------------------------------------------------------------------
# Load all upstream
# -----------------------------------------------------------------------------
locus <- fread(LOCUS_FILE)
cojo  <- fread(COJO_FILE)
susie <- fread(SUSIE_FILE)

cat("Loaded:\n")
cat("  Locus variants:", nrow(locus), "\n")
cat("  COJO signals:  ", nrow(cojo), "\n")
cat("  SuSiE CS:      ", sum(!is.na(susie$cs_id)), "\n")

# -----------------------------------------------------------------------------
# Pivot COJO and SuSiE to wide-by-ancestry format
# -----------------------------------------------------------------------------
cojo_eur <- cojo[ld_reference_id == "1000G_phase3_EUR_hg19"]
cojo_eas <- cojo[ld_reference_id == "1000G_phase3_EAS_hg19"]
susie_eur <- susie[ld_reference_id == "1000G_phase3_EUR_hg19" & !is.na(cs_id)]
susie_eas <- susie[ld_reference_id == "1000G_phase3_EAS_hg19" & !is.na(cs_id)]

# Build signal_id list — union of all SNPs flagged as independent in any analysis
all_signals <- unique(c(
  cojo_eur$lead_snp,
  cojo_eas$lead_snp,
  susie_eur$top_snp,
  susie_eas$top_snp
))
all_signals <- all_signals[!is.na(all_signals)]

cat("Total candidate independent signals across methods:", length(all_signals), "\n")

# -----------------------------------------------------------------------------
# Build integration table, one row per candidate signal
# -----------------------------------------------------------------------------
build_signal_row <- function(snp, idx) {
  c_eur <- cojo_eur[lead_snp == snp]
  c_eas <- cojo_eas[lead_snp == snp]
  s_eur <- susie_eur[top_snp == snp]
  s_eas <- susie_eas[top_snp == snp]
  marg  <- locus[snp == ..snp]

  data.frame(
    signal_id = sprintf("TSHR_S%d", idx),
    lead_snp_cojo_eur  = if (nrow(c_eur) > 0) c_eur$lead_snp[1] else NA_character_,
    lead_snp_cojo_eas  = if (nrow(c_eas) > 0) c_eas$lead_snp[1] else NA_character_,
    lead_snp_susie_eur = if (nrow(s_eur) > 0) s_eur$top_snp[1] else NA_character_,
    lead_snp_susie_eas = if (nrow(s_eas) > 0) s_eas$top_snp[1] else NA_character_,
    concordant_across_methods =
      length(unique(na.omit(c(c_eur$lead_snp[1], s_eur$top_snp[1])))) == 1 ||
      length(unique(na.omit(c(c_eas$lead_snp[1], s_eas$top_snp[1])))) == 1,
    concordant_across_ancestries =
      !is.na(c_eur$lead_snp[1]) && !is.na(c_eas$lead_snp[1]) &&
      c_eur$lead_snp[1] == c_eas$lead_snp[1],
    chr = 14L,
    pos_hg19 = if (nrow(marg) > 0) marg$pos_hg19[1] else NA_integer_,
    pos_hg38 = if (nrow(marg) > 0) marg$pos_hg38[1] else NA_integer_,
    marginal_p   = if (nrow(marg) > 0) marg$pvalue[1] else NA_real_,
    conditional_p_eur = if (nrow(c_eur) > 0) c_eur$conditional_p[1] else NA_real_,
    conditional_p_eas = if (nrow(c_eas) > 0) c_eas$conditional_p[1] else NA_real_,
    susie_pip_eur     = if (nrow(s_eur) > 0) s_eur$top_snp_pip[1] else NA_real_,
    susie_pip_eas     = if (nrow(s_eas) > 0) s_eas$top_snp_pip[1] else NA_real_,
    ld_with_rs179252_eur = if (nrow(c_eur) > 0) c_eur$ld_with_rs179252_eur[1] else NA_real_,
    ld_with_rs179252_eas = if (nrow(c_eas) > 0) c_eas$ld_with_rs179252_eas[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

rows <- lapply(seq_along(all_signals), function(i) {
  build_signal_row(all_signals[i], i)
})
combined <- do.call(rbind, rows)

# -----------------------------------------------------------------------------
# Decision: usable_as_independent_iv
# -----------------------------------------------------------------------------
combined$usable_as_independent_iv <- with(combined, {
  cond_p_ok <- pmin(conditional_p_eur, conditional_p_eas, na.rm = TRUE) < 5e-6
  ld_low <- pmax(ld_with_rs179252_eur, ld_with_rs179252_eas, na.rm = TRUE) < 0.1
  pip_high <- pmax(susie_pip_eur, susie_pip_eas, na.rm = TRUE) > 0.5
  # rs179252 itself is by definition the original signal — not "independent" from itself
  is_rs179252 <- (lead_snp_cojo_eur == "rs179252" |
                  lead_snp_susie_eur == "rs179252")
  is_rs179252[is.na(is_rs179252)] <- FALSE

  ifelse(is_rs179252, FALSE, cond_p_ok & ld_low & pip_high)
})

combined$notes <- with(combined, ifelse(
  is.na(usable_as_independent_iv), "Insufficient data",
  ifelse(usable_as_independent_iv, "Usable as independent IV for multi-IV MR",
         "Not independent (may be collinear with rs179252 or weak signal)")
))

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
fwrite(combined, OUTPUT_FILE)

cat("\nWrote:", OUTPUT_FILE, "(", nrow(combined), "signals)\n")
cat("\nSignals usable as independent IVs:",
    sum(combined$usable_as_independent_iv, na.rm = TRUE), "\n")

print(combined[, c("signal_id", "lead_snp_cojo_eur", "lead_snp_cojo_eas",
                   "conditional_p_eur", "conditional_p_eas",
                   "ld_with_rs179252_eur", "ld_with_rs179252_eas",
                   "usable_as_independent_iv")])

# -----------------------------------------------------------------------------
# Log + decision flag for Task C
# -----------------------------------------------------------------------------
n_usable <- sum(combined$usable_as_independent_iv, na.rm = TRUE)

decision_msg <- if (n_usable >= 1) {
  sprintf("RESOLVED: %d independent signal(s) beyond rs179252 found. Multi-IV MR feasible.", n_usable)
} else {
  "UNRESOLVED: No independent signals beyond rs179252. TSHR remains single-IV. Reframe as 'locus-level prioritization' not 'robust multi-IV MR'."
}

cat("\nTSHR multi-IV decision:", decision_msg, "\n")

log_entry <- sprintf(
  "## %s — Task B.5 integration\n- Output: %s\n- Total candidate signals: %d\n- Usable as independent IVs: %d\n- Decision: %s\n\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  basename(OUTPUT_FILE), nrow(combined), n_usable, decision_msg
)
cat(log_entry, file = LOG_FILE, append = TRUE)
