#==============================================================================
# TED-TRAP Upgrade — Phase 1 Script 08
# Multivariable MR: TED effect independent of underlying Graves disease
#==============================================================================
# Purpose:
#   Resolve MAJOR issue #8: Graves disease (GD) ≠ TED; only 20-40% of GD
#   patients develop clinically significant TED. Primary MR on GD genetics
#   may reflect GD liability rather than TED-specific causation.
#
#   This script runs MVMR (Sanderson 2019) with:
#     - Exposure 1: TSHR cis-eQTL expression
#     - Exposure 2: GD genetic liability (separate instrument set)
#     - Outcome: FinnGen R12 Graves ophthalmopathy
#
#   The goal is to isolate the DIRECT effect of TSHR expression on TED,
#   conditional on (adjusting for) GD liability.
#
# Reference: Sanderson E, Davey Smith G, Windmeijer F, Bowden J (2019).
#            An examination of multivariable Mendelian randomization in the
#            single-sample and two-sample summary data settings.
#            Int J Epidemiol 48:713–727.
#
# Output:
#   TrackA_MR/results/08_mvmr_ted_given_gd.csv
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(TwoSampleMR)
library(MVMR)
library(MendelianRandomization)
library(data.table)

log_file <- "TrackA_MR/logs/08_mvmr.log"
sink(log_file, split = TRUE)
cat("=== Multivariable MR: TSHR → TED | GD ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- Exposures ---
# Primary: TSHR cis-eQTL from eQTLGen
tshr_id  <- "eqtl-a-ENSG00000165409"
# Second: GD liability genetic instrument set (top GWAS SNPs from Graves GWAS)
gd_id    <- "ebi-a-GCST90018627"

# Outcome: FinnGen R12 Graves ophthalmopathy (TED-proximal)
outcome_id <- "finngen_R12_E4_GRAVES_OPHT"

# --- Extract joint instrument set ---
# For MVMR, instruments need to be chosen for BOTH exposures, clumped,
# and the harmonised data brought together.
cat("Step 1: Extract TSHR instruments\n")
tshr_ivs <- extract_instruments(outcomes = tshr_id, p1 = 5e-6,
                                  clump = TRUE, r2 = 0.001, kb = 10000)
cat(sprintf("  TSHR: %d IVs\n", nrow(tshr_ivs)))

cat("Step 2: Extract GD instruments\n")
gd_ivs <- extract_instruments(outcomes = gd_id, p1 = 5e-8,
                                clump = TRUE, r2 = 0.001, kb = 10000)
cat(sprintf("  GD: %d IVs\n", nrow(gd_ivs)))

# --- Joint SNP list ---
all_ivs <- unique(c(tshr_ivs$SNP, gd_ivs$SNP))
cat(sprintf("  Combined IV panel: %d unique SNPs\n\n", length(all_ivs)))

# --- Re-extract both exposures on the joint SNP list ---
tshr_full <- extract_outcome_data(snps = all_ivs, outcomes = tshr_id)
gd_full   <- extract_outcome_data(snps = all_ivs, outcomes = gd_id)
out_full  <- extract_outcome_data(snps = all_ivs, outcomes = outcome_id)

# Keep only SNPs present in all three
common <- Reduce(intersect, list(tshr_full$SNP, gd_full$SNP, out_full$SNP))
cat(sprintf("SNPs available across all 3 datasets: %d\n", length(common)))

tshr_full <- tshr_full[tshr_full$SNP %in% common, ]
gd_full   <- gd_full[gd_full$SNP %in% common, ]
out_full  <- out_full[out_full$SNP %in% common, ]

# Align order
tshr_full <- tshr_full[match(common, tshr_full$SNP), ]
gd_full   <- gd_full[match(common, gd_full$SNP), ]
out_full  <- out_full[match(common, out_full$SNP), ]

# --- Format for MVMR ---
# Required format: BXGs (matrix of beta exposures), seBXGs, BYG, seBYG
BXGs <- cbind(tshr_full$beta.outcome, gd_full$beta.outcome)
seBXGs <- cbind(tshr_full$se.outcome, gd_full$se.outcome)
colnames(BXGs) <- c("TSHR_exp", "GD_liability")
BYG <- out_full$beta.outcome
seBYG <- out_full$se.outcome

# --- Conditional F-statistic (Sanderson & Windmeijer 2016) ---
mvmr_input <- format_mvmr(
  BXGs  = BXGs,
  BYG   = BYG,
  seBXGs = seBXGs,
  seBYG  = seBYG,
  RSID   = common
)

cat("Computing conditional F-statistics...\n")
strength_mvmr <- strength_mvmr(r_input = mvmr_input, gencov = 0)
cat("Conditional F-stats per exposure:\n")
print(strength_mvmr)

# --- IVW MVMR estimate ---
mvmr_ivw <- ivw_mvmr(r_input = mvmr_input)
cat("\n--- MVMR IVW Estimates ---\n")
print(mvmr_ivw)

# --- Heterogeneity (Q_A) ---
pleio_test <- tryCatch(pleiotropy_mvmr(r_input = mvmr_input, gencov = 0),
                       error = function(e) NULL)
if (!is.null(pleio_test)) {
  cat("\n--- MVMR heterogeneity (Q_A) test ---\n")
  print(pleio_test)
}

# --- Save ---
out_df <- data.frame(
  exposure = c("TSHR_eQTL", "GD_liability"),
  effect_on_TED    = mvmr_ivw[, "Estimate"],
  se               = mvmr_ivw[, "Std. Error"],
  t_value          = mvmr_ivw[, "t value"],
  p_value          = mvmr_ivw[, "Pr(>|t|)"],
  conditional_F    = as.numeric(strength_mvmr),
  stringsAsFactors = FALSE
)
fwrite(out_df, "TrackA_MR/results/08_mvmr_ted_given_gd.csv")

cat("\n\n=== INTERPRETATION ===\n")
cat("If TSHR's effect on TED persists AFTER adjusting for GD liability\n")
cat("(i.e., p_TSHR remains significant in the MVMR model), this provides\n")
cat("evidence that TSHR expression affects TED through pathways beyond\n")
cat("mere GD susceptibility — supporting tissue-level TSHR biology.\n\n")
cat("If TSHR's effect attenuates to null, then the primary MR signal is\n")
cat("likely mediated through GD development itself, NOT through a TED-\n")
cat("specific mechanism. This would be an important honest limitation.\n")

sink()
cat("\n✅ Phase 1 Script 08 complete.\n\n")
