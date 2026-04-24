# =============================================================================
# TED-TRAP Step 4 — Multivariable MR (MVMR)
# Purpose: Isolate TED-specific causal effect of TSHR independent of GD liability
# SE1 test: Is TSHR signal preserved after adjusting for Graves disease?
# =============================================================================

required_pkgs <- c("TwoSampleMR", "MVMR", "ieugwasr")
for (p in required_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
        if (p == "MVMR") {
            if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "http://cran.us.r-project.org")
            Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
            remotes::install_github("WSpiller/MVMR")
        } else {
            install.packages(p, repos = "http://cran.us.r-project.org")
        }
    }
}
library(TwoSampleMR)
library(MVMR)
library(ieugwasr)

token <- Sys.getenv("OPENGWAS_JWT")
# options(ieugwasr_api = "https://api.opengwas.io/")

if (!dir.exists("TrackA_MR/results")) dir.create("TrackA_MR/results", recursive = TRUE)
if (!dir.exists("TrackA_MR/logs")) dir.create("TrackA_MR/logs", recursive = TRUE)
log_file <- "TrackA_MR/logs/04_mvmr.log"
sink(log_file, split = TRUE)

cat("=== Step 4: MVMR — TED-specific signal isolation ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

EXPOSURE_1_TSHR <- "eqtl-a-ENSG00000165409"
EXPOSURE_2_GD <- "ebi-a-GCST90018627"
OUTCOME_TED <- "finn-b-E4_GRAVES_OPHT"
OUTCOME_TED_ALTS <- c("finn-b-E4_GRAVES_OPHT", "finn-b-E4_GRAVES_STRICT", "finn-b-E4_THYROID")

cat("Design:\n")
cat(sprintf("  Exposure 1: TSHR cis-eQTL (%s)\n", EXPOSURE_1_TSHR))
cat(sprintf("  Exposure 2: GD liability (%s)\n", EXPOSURE_2_GD))
cat(sprintf("  Outcome:    TED/GO (%s)\n\n", OUTCOME_TED))

cat("=== 4.1 Extracting instruments ===\n\n")

cat("Exposure 1: TSHR cis-eQTL\n")
exp1_tshr <- tryCatch(
    extract_instruments(outcomes = EXPOSURE_1_TSHR, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000),
    error = function(e) NULL
)

if (is.null(exp1_tshr) || nrow(exp1_tshr) < 1) {
    cat("  No IVs at P<5e-8. Falling back to P<5e-6...\n")
    exp1_tshr <- extract_instruments(outcomes = EXPOSURE_1_TSHR, p1 = 5e-6, clump = TRUE, r2 = 0.001, kb = 10000)
}
if (!is.null(exp1_tshr)) cat(sprintf("  TSHR IVs: %d SNPs\n", nrow(exp1_tshr)))

cat("\nExposure 2: Graves disease liability\n")
exp2_gd <- extract_instruments(outcomes = EXPOSURE_2_GD, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
if (!is.null(exp2_gd)) cat(sprintf("  GD IVs at P<5e-8: %d SNPs\n", nrow(exp2_gd)))

all_snps <- unique(c(if (!is.null(exp1_tshr)) exp1_tshr$SNP, if (!is.null(exp2_gd)) exp2_gd$SNP))
cat(sprintf("\nCombined unique IVs: %d SNPs\n", length(all_snps)))

cat("\n=== 4.2 Harmonizing IV effects across both exposures ===\n")
tshr_effects <- extract_outcome_data(snps = all_snps, outcomes = EXPOSURE_1_TSHR)
cat(sprintf("  TSHR effects retrieved: %d SNPs\n", nrow(tshr_effects)))

gd_effects <- extract_outcome_data(snps = all_snps, outcomes = EXPOSURE_2_GD)
cat(sprintf("  GD effects retrieved: %d SNPs\n", nrow(gd_effects)))

common_snps <- intersect(tshr_effects$SNP, gd_effects$SNP)
cat(sprintf("  SNPs in both exposures: %d\n", length(common_snps)))

cat("\n=== 4.3 Retrieving outcome data (TED/FinnGen) ===\n")
outcome_data <- NULL
OUTCOME_USED <- NULL
for (oc_try in OUTCOME_TED_ALTS) {
    cat(sprintf("  Trying %s ...\n", oc_try))
    outcome_data <- tryCatch(extract_outcome_data(snps = common_snps, outcomes = oc_try), error = function(e) {
        cat(sprintf("    [ERROR] %s\n", e$message))
        NULL
    })
    if (!is.null(outcome_data) && nrow(outcome_data) > 5) {
        cat(sprintf("    ✅ %s returned %d SNPs\n", oc_try, nrow(outcome_data)))
        OUTCOME_USED <- oc_try
        break
    }
}

if (is.null(outcome_data) || nrow(outcome_data) < 5) {
    cat("\n⚠️  FinnGen TED outcome fetch failed. Trying alternative:\n")
    outcome_data <- tryCatch(extract_outcome_data(snps = common_snps, outcomes = "ebi-a-GCST90038636"), error = function(e) NULL)
    OUTCOME_USED <- "ebi-a-GCST90038636"
    if (!is.null(outcome_data)) cat(sprintf("  Using UKB Hyperthyroidism instead: %d SNPs\n", nrow(outcome_data)))
}

cat(sprintf("\nFinal outcome used: %s\n", OUTCOME_USED))

cat("\n=== 4.4 Building MVMR input matrix ===\n")
final_snps <- Reduce(intersect, list(tshr_effects$SNP, gd_effects$SNP, outcome_data$SNP))
cat(sprintf("SNPs in all 3 datasets: %d\n", length(final_snps)))

tshr_sub <- tshr_effects[match(final_snps, tshr_effects$SNP), ]
gd_sub <- gd_effects[match(final_snps, gd_effects$SNP), ]
out_sub <- outcome_data[match(final_snps, outcome_data$SNP), ]

mvmr_input <- data.frame(
    SNP = final_snps,
    beta_TSHR = tshr_sub$beta.outcome, se_TSHR = tshr_sub$se.outcome, ea_TSHR = tshr_sub$effect_allele.outcome, oa_TSHR = tshr_sub$other_allele.outcome,
    beta_GD = gd_sub$beta.outcome, se_GD = gd_sub$se.outcome, ea_GD = gd_sub$effect_allele.outcome, oa_GD = gd_sub$other_allele.outcome,
    beta_OUT = out_sub$beta.outcome, se_OUT = out_sub$se.outcome, ea_OUT = out_sub$effect_allele.outcome, oa_OUT = out_sub$other_allele.outcome,
    stringsAsFactors = FALSE
)

flip_gd <- mvmr_input$ea_TSHR != mvmr_input$ea_GD & mvmr_input$ea_TSHR == mvmr_input$oa_GD
mvmr_input$beta_GD[flip_gd] <- -mvmr_input$beta_GD[flip_gd]

flip_out <- mvmr_input$ea_TSHR != mvmr_input$ea_OUT & mvmr_input$ea_TSHR == mvmr_input$oa_OUT
mvmr_input$beta_OUT[flip_out] <- -mvmr_input$beta_OUT[flip_out]

mvmr_input <- mvmr_input[!is.na(mvmr_input$beta_GD) & !is.na(mvmr_input$beta_OUT), ]
cat(sprintf("Final harmonized SNPs: %d\n", nrow(mvmr_input)))

cat("\n=== 4.5 Running MVMR (MVMR package) ===\n\n")
if (nrow(mvmr_input) < 2) stop("MVMR requires at least 2 SNPs. Aborting.")

F_data <- format_mvmr(BXGs = cbind(mvmr_input$beta_TSHR, mvmr_input$beta_GD), seBXGs = cbind(mvmr_input$se_TSHR, mvmr_input$se_GD), BYG = mvmr_input$beta_OUT, seBYG = mvmr_input$se_OUT, RSID = mvmr_input$SNP)

cat("--- Conditional F-statistics (instrument strength) ---\n  (F > 10 indicates strong instruments)\n\n")
Fstats <- strength_mvmr(r_input = F_data, gencov = 0)
print(Fstats)

cat("\n--- MVMR-IVW causal estimates ---\n")
res_mvmr <- ivw_mvmr(r_input = F_data)
print(res_mvmr)

tshr_b <- res_mvmr[1, "Estimate"]
tshr_se <- res_mvmr[1, "Std. Error"]
tshr_p <- res_mvmr[1, "Pr(>|t|)"]
gd_b <- res_mvmr[2, "Estimate"]
gd_se <- res_mvmr[2, "Std. Error"]
gd_p <- res_mvmr[2, "Pr(>|t|)"]

cat("\n--- Heterogeneity (Cochran's Q_A) ---\n")
q_het <- tryCatch(pleiotropy_mvmr(r_input = F_data, gencov = 0), error = function(e) {
    cat(sprintf("Q calc error: %s\n", e$message))
    NULL
})
if (!is.null(q_het)) print(q_het)

cat("\n\n=== 4.6 Comparison: Univariable vs Multivariable ===\n\n")
uvmr_dat <- data.frame(SNP = mvmr_input$SNP, beta.exposure = mvmr_input$beta_TSHR, se.exposure = mvmr_input$se_TSHR, beta.outcome = mvmr_input$beta_OUT, se.outcome = mvmr_input$se_OUT, effect_allele.exposure = mvmr_input$ea_TSHR, other_allele.exposure = mvmr_input$oa_TSHR, mr_keep = TRUE, exposure = "TSHR", outcome = "TED", id.exposure = "TSHR", id.outcome = "TED")
uvmr_res <- mr(uvmr_dat, method_list = c("mr_ivw", "mr_wald_ratio"))
cat("Univariable MR (TSHR → TED):\n")
print(uvmr_res[, c("method", "nsnp", "b", "se", "pval")])

cat("\n\n", rep("=", 70), "\nSE1 — TED-SPECIFIC SIGNAL (MVMR)\n", rep("=", 70), "\n\n", sep = "")
cat("Primary MVMR estimate (TSHR adjusted for GD liability):\n")
cat(sprintf("  TSHR β_adj = %+.4f (SE=%.4f)  P=%.3e\n", tshr_b, tshr_se, tshr_p))
cat(sprintf("  GD liability β = %+.4f (SE=%.4f)  P=%.3e\n", gd_b, gd_se, gd_p))

cat("\nUnivariable reference (TSHR → TED):\n")
uvmr_beta <- uvmr_res$b[1]
uvmr_p <- uvmr_res$pval[1]
cat(sprintf("  TSHR β_uv  = %+.4f  P=%.3e\n", uvmr_beta, uvmr_p))

attenuation <- abs(tshr_b / uvmr_beta)
cat(sprintf("\nAttenuation ratio: β_adj / β_uv = %.3f\n", tshr_b / uvmr_beta))

cat("\n--- SE1 Verdict ---\n")
if (tshr_p < 0.05 && attenuation > 0.7) {
    cat("✅ SE1 ACHIEVED: TSHR retains TED-specific causal effect after GD adjustment.\n   Narrative: 'TSHR directly causes TED independent of Graves disease liability.'\n")
} else if (tshr_p < 0.05 && attenuation > 0.3) {
    cat("🟡 SE1 PARTIAL: TSHR effect attenuated but retained.\n   Narrative: 'TSHR has both GD-mediated and TED-direct effects.'\n")
} else if (attenuation < 0.3) {
    cat("⚠️  SE1 NOT ACHIEVED: TSHR effect largely mediated via GD liability.\n   Narrative: 'TSHR → GD → TED; direct TED effect is modest.'\n")
} else {
    cat("❓ SE1 INCONCLUSIVE: Low power (n SNPs small/non-sig).\n")
}

mvmr_summary <- data.frame(Exposure = c("TSHR", "GD_liability"), Beta_adjusted = c(tshr_b, gd_b), SE_adjusted = c(tshr_se, gd_se), P_adjusted = c(tshr_p, gd_p), N_SNPs = nrow(mvmr_input), Outcome = OUTCOME_USED, stringsAsFactors = FALSE)
write.csv(mvmr_summary, "TrackA_MR/results/04_mvmr_summary.csv", row.names = FALSE)
write.csv(as.data.frame(Fstats), "TrackA_MR/results/04_mvmr_Fstats.csv", row.names = FALSE)
saveRDS(list(mvmr = res_mvmr, uvmr = uvmr_res, Fstats = Fstats, q = q_het, input = mvmr_input), "TrackA_MR/results/04_mvmr_full.rds")

cat("\n💾 Saved:\n  TrackA_MR/results/04_mvmr_summary.csv\n  TrackA_MR/results/04_mvmr_Fstats.csv\n  TrackA_MR/results/04_mvmr_full.rds\n")
sink()
cat("\n✅ Step 4 complete.\n")
