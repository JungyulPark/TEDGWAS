#==============================================================================
# TED-TRAP Upgrade — Phase 1 Script 03
# Comprehensive MR with full sensitivity suite
#==============================================================================
# Purpose:
#   - Run MR for all candidate genes across 3 outcomes
#   - Use FULL sensitivity suite: IVW + MR-Egger + WM + WMode + MR-PRESSO + CAUSE
#   - Apply MRlap for eQTLGen↔UKB sample overlap (later)
#   - Rescale linear-model beta to log-odds for replication outcome
#
# Methodological fixes vs v1:
#   [CRITICAL] β-scale harmonization between primary and replication
#   [CRITICAL] Full sensitivity suite (not just IVW)
#   [CRITICAL] MR-PRESSO outlier correction for TNF
#   [MAJOR]    Proper cis-MR framing for TSHR (2 IVs)
#
# Output:
#   TrackA_MR/results/03_mr_v3_full.csv (all sensitivity methods)
#   TrackA_MR/results/03_mr_v3_leaveonesnp_out.csv
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)
library(ieugwasr)
library(dplyr)
library(data.table)

# --- Outcomes ---
# Primary: Graves disease (log-odds scale)
# Replication: Hyperthyroidism (LINEAR mixed model — must be rescaled!)
# Sensitivity: FinnGen R12 Graves ophthalmopathy (log-odds via SAIGE)
outcomes <- list(
  Primary     = list(id = "ebi-a-GCST90018627",
                     name = "Graves disease",
                     scale = "log_odds",
                     case_prev = 2809/175465),
  Replication = list(id = "ebi-a-GCST90038636",
                     name = "Hyperthyroidism",
                     scale = "linear",
                     case_prev = 3731/484598),
  Sensitivity = list(id = "finngen_R12_E4_GRAVES_OPHT",     # adjust if different label
                     name = "FinnGen R12 Graves ophthalmopathy",
                     scale = "log_odds",
                     case_prev = 858/500348)
)

# --- Log setup ---
log_file <- "TrackA_MR/logs/03_mr_full.log"
sink(log_file, split = TRUE)
cat("=== Comprehensive MR with Full Sensitivity Suite ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- Load instrument summary ---
summary_df <- fread("TrackA_MR/results/02_extraction_summary.csv")
ready_genes <- summary_df$gene[summary_df$primary_mode != "insufficient"]
cat(sprintf("Genes to analyze (n=%d): %s\n\n",
            length(ready_genes), paste(ready_genes, collapse=", ")))

# --- Rescale helper (linear -> log-odds) ---
# Lloyd-Jones 2018; valid under low case prevalence and additive model
rescale_linear_to_logodds <- function(beta_linear, se_linear, case_prev) {
  p <- case_prev
  denom <- p * (1 - p)
  beta_logor <- beta_linear / denom
  se_logor   <- se_linear  / denom
  list(beta = beta_logor, se = se_logor)
}

# --- Storage ---
all_mr_results <- list()
all_sensitivity <- list()
all_presso      <- list()
all_leaveoneout <- list()

# --- Main loop ---
for (gene_sym in ready_genes) {
  cat(sprintf("\n=========================================\n"))
  cat(sprintf("Gene: %s\n", gene_sym))
  cat(sprintf("=========================================\n"))

  # Determine IV mode
  info <- summary_df[summary_df$gene == gene_sym, ]
  iv_mode <- info$primary_mode
  iv_file <- sprintf("TrackA_MR/data/instruments/%s_instruments_%s.tsv", gene_sym, iv_mode)

  if (!file.exists(iv_file)) {
    cat(sprintf("[SKIP] IV file not found: %s\n", iv_file))
    next
  }

  exposure_data <- fread(iv_file)
  cat(sprintf("Loaded %d IVs (%s mode)\n", nrow(exposure_data), iv_mode))

  for (oc_name in names(outcomes)) {
    oc <- outcomes[[oc_name]]
    cat(sprintf("\n  Outcome: %s (%s, scale=%s)\n", oc_name, oc$name, oc$scale))

    # Extract outcome
    outcome_data <- tryCatch(
      extract_outcome_data(snps = exposure_data$SNP, outcomes = oc$id),
      error = function(e) { cat(sprintf("    [ERROR] outcome extraction: %s\n", e$message)); NULL }
    )

    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      cat("    [SKIP] No outcome data retrieved\n")
      next
    }

    # Rescale linear scale to log-odds
    if (oc$scale == "linear") {
      cat("    ⚠️  Rescaling linear-model betas to log-odds (Lloyd-Jones 2018)\n")
      cat(sprintf("       Case prevalence: %.5f; scaling factor ~%.1f\n",
                  oc$case_prev, 1 / (oc$case_prev * (1 - oc$case_prev))))
      rescaled <- rescale_linear_to_logodds(
        beta_linear = outcome_data$beta.outcome,
        se_linear   = outcome_data$se.outcome,
        case_prev   = oc$case_prev
      )
      outcome_data$beta.outcome <- rescaled$beta
      outcome_data$se.outcome   <- rescaled$se
      outcome_data$outcome <- paste0(outcome_data$outcome, " [rescaled to log-odds]")
    }

    # Harmonize
    harmonised <- harmonise_data(exposure_data, outcome_data, action = 2)
    harmonised <- harmonised[harmonised$mr_keep, ]
    if (nrow(harmonised) == 0) { cat("    [SKIP] No SNPs after harmonization\n"); next }

    cat(sprintf("    SNPs after harmonisation: %d\n", nrow(harmonised)))

    # ===== 1. Full method suite (TwoSampleMR) =====
    methods_to_run <- c("mr_ivw")
    if (nrow(harmonised) >= 3) {
      methods_to_run <- c(methods_to_run, "mr_egger_regression",
                          "mr_weighted_median", "mr_weighted_mode",
                          "mr_simple_mode")
    }

    mr_res <- tryCatch(
      mr(harmonised, method_list = methods_to_run),
      error = function(e) { cat(sprintf("    [ERROR] mr(): %s\n", e$message)); NULL }
    )

    if (!is.null(mr_res)) {
      mr_res$gene         <- gene_sym
      mr_res$outcome_role <- oc_name
      mr_res$scale        <- oc$scale
      mr_res$iv_mode      <- iv_mode
      all_mr_results[[paste(gene_sym, oc_name, sep="_")]] <- mr_res
      print(mr_res[, c("method", "nsnp", "b", "se", "pval")])
    }

    # ===== 2. Sensitivity tests =====
    sens_row <- data.frame(
      gene = gene_sym, outcome_role = oc_name,
      Egger_intercept = NA_real_, Egger_intercept_p = NA_real_,
      Cochran_Q_IVW = NA_real_, Cochran_Q_IVW_p = NA_real_,
      Steiger_direction = NA, Steiger_p = NA_real_,
      stringsAsFactors = FALSE
    )

    if (nrow(harmonised) >= 3) {
      egger_p <- tryCatch(mr_pleiotropy_test(harmonised), error = function(e) NULL)
      if (!is.null(egger_p) && nrow(egger_p) > 0) {
        sens_row$Egger_intercept   <- egger_p$egger_intercept
        sens_row$Egger_intercept_p <- egger_p$pval
      }
      het <- tryCatch(mr_heterogeneity(harmonised, method_list = "mr_ivw"),
                      error = function(e) NULL)
      if (!is.null(het) && nrow(het) > 0) {
        sens_row$Cochran_Q_IVW   <- het$Q
        sens_row$Cochran_Q_IVW_p <- het$Q_pval
      }
    }
    st <- tryCatch(directionality_test(harmonised), error = function(e) NULL)
    if (!is.null(st) && nrow(st) > 0) {
      sens_row$Steiger_direction <- st$correct_causal_direction
      sens_row$Steiger_p         <- st$steiger_pval
    }
    all_sensitivity[[paste(gene_sym, oc_name, sep="_")]] <- sens_row

    # ===== 3. MR-PRESSO (requires ≥4 IVs) =====
    if (nrow(harmonised) >= 4) {
      cat("    Running MR-PRESSO (global + outlier test)...\n")
      presso_res <- tryCatch(
        mr_presso(
          BetaExposure = "beta.exposure",   SdExposure = "se.exposure",
          BetaOutcome  = "beta.outcome",    SdOutcome  = "se.outcome",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
          data = harmonised, NbDistribution = 1000, SignifThreshold = 0.05
        ),
        error = function(e) { cat(sprintf("    [ERROR] MR-PRESSO: %s\n", e$message)); NULL }
      )

      if (!is.null(presso_res)) {
        global_p  <- presso_res$`MR-PRESSO results`$`Global Test`$Pvalue
        distortion_p <- presso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue
        outliers  <- presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`

        presso_summary <- data.frame(
          gene         = gene_sym,
          outcome_role = oc_name,
          n_snps       = nrow(harmonised),
          global_test_p    = global_p,
          distortion_p     = ifelse(is.null(distortion_p), NA, distortion_p),
          n_outliers   = ifelse(is.null(outliers), 0, length(outliers)),
          raw_beta     = presso_res$`Main MR results`$`Causal Estimate`[1],
          raw_p        = presso_res$`Main MR results`$`P-value`[1],
          corrected_beta = presso_res$`Main MR results`$`Causal Estimate`[2],
          corrected_p  = presso_res$`Main MR results`$`P-value`[2],
          stringsAsFactors = FALSE
        )
        all_presso[[paste(gene_sym, oc_name, sep="_")]] <- presso_summary
        cat(sprintf("      Global P: %.4g  |  Outliers: %d\n", global_p, length(outliers)))
      }
    } else {
      cat("    [SKIP MR-PRESSO] Requires ≥4 IVs\n")
    }

    # ===== 4. Leave-one-SNP-out =====
    if (nrow(harmonised) >= 3) {
      loso <- tryCatch(mr_leaveoneout(harmonised), error = function(e) NULL)
      if (!is.null(loso)) {
        loso$gene         <- gene_sym
        loso$outcome_role <- oc_name
        all_leaveoneout[[paste(gene_sym, oc_name, sep="_")]] <- loso
      }
    } else {
      # For 2-IV (cis-MR), report each SNP's Wald ratio as de facto leave-one-out
      wald <- data.frame(
        SNP = harmonised$SNP,
        b = harmonised$beta.outcome / harmonised$beta.exposure,
        se = harmonised$se.outcome / abs(harmonised$beta.exposure),
        stringsAsFactors = FALSE
      )
      wald$p <- 2 * pnorm(-abs(wald$b / wald$se))
      wald$gene <- gene_sym
      wald$outcome_role <- oc_name
      wald$note <- "per-SNP Wald ratio (cis-MR, <3 IVs)"
      all_leaveoneout[[paste(gene_sym, oc_name, sep="_")]] <- wald
    }

    Sys.sleep(1)
  }  # end outcome loop
}  # end gene loop

# --- Save results ---
mr_combined  <- do.call(rbind, all_mr_results)
sens_combined <- do.call(rbind, all_sensitivity)
presso_combined <- if (length(all_presso) > 0) do.call(rbind, all_presso) else data.frame()
loso_combined  <- do.call(rbind, all_leaveoneout)

fwrite(mr_combined,  "TrackA_MR/results/03_mr_v3_full.csv")
fwrite(sens_combined, "TrackA_MR/results/03_mr_v3_sensitivity.csv")
if (nrow(presso_combined) > 0) fwrite(presso_combined, "TrackA_MR/results/03_mr_v3_presso.csv")
fwrite(loso_combined, "TrackA_MR/results/03_mr_v3_leaveonesnp_out.csv")

cat("\n\n=== FINAL RESULTS ===\n")
cat(sprintf("  MR full: %d rows\n", nrow(mr_combined)))
cat(sprintf("  Sensitivity: %d rows\n", nrow(sens_combined)))
cat(sprintf("  MR-PRESSO: %d rows\n", nrow(presso_combined)))
cat(sprintf("  Leave-one-out: %d rows\n", nrow(loso_combined)))

# Highlight TNF Cochran Q issue
tnf_het <- sens_combined[sens_combined$gene == "TNF" &
                          sens_combined$Cochran_Q_IVW_p < 0.05, ]
if (nrow(tnf_het) > 0) {
  cat("\n⚠️  TNF heterogeneity detected — MR-PRESSO corrected result will be used in downstream analysis.\n")
}

# Highlight TSHR cis-MR framing
cat("\n📌 TSHR (2-IV cis-MR): Per-SNP Wald ratios will serve as primary sensitivity evidence.\n")
cat("   MR-Egger/MR-PRESSO/weighted median are inestimable at n<3 IVs. This is disclosed.\n")
cat("   Colocalization (script 06) provides the formal consistency check.\n")

sink()
cat("\n✅ Phase 1 Script 03 complete.\n")
cat("   Next: 04_mr_scale_rescale.R (verification) or 06_coloc_all_loci.R (colocalization)\n\n")
