# =============================================================================
# TED-TRAP Step 5 — Cross-Phenotype MR for Off-Target Effects
#
# Hypothesis: IGF1R genetic upregulation has effects on teprotumumab's known
# adverse event phenotypes (hearing loss, hyperglycemia), providing genetic
# support for IGF-1R blockade's off-target biology.
#
# Design:
#   Exposure: IGF1R cis-eQTL (same 10 IVs from Step 2)
#   Outcomes:
#     Primary off-target:   Hearing difficulty / aid use (UKB, Wells 2019)
#     Secondary off-target: Fasting glucose (MAGIC)
#     Tertiary off-target:  HbA1c (MAGIC)
#     Negative controls:    Height (ieu-a-89); Skin color (ukb-a-*)
# =============================================================================

suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(ieugwasr)
    library(data.table)
})

token <- Sys.getenv("OPENGWAS_JWT")
# options(ieugwasr_api = "https://api.opengwas.io/")

for (d in c("TrackA_MR/results", "TrackA_MR/logs")) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

log_file <- "TrackA_MR/logs/05_cross_phenotype_mr.log"
sink(log_file, split = TRUE)

cat("=== Step 5: Cross-Phenotype MR for IGF1R Off-Target ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# =============================================================================
# CONFIG
# =============================================================================

EXPOSURE_IGF1R <- "eqtl-a-ENSG00000140443"

# Outcomes grouped by role
off_target_outcomes <- list(
    # Primary off-targets (SE2)
    Hearing_Difficulty = list(
        id = "ukb-b-12002",
        label = "Self-reported hearing difficulty (UKB)",
        role = "primary_offtarget"
    ),
    Hearing_AidUse = list(
        id = "ukb-b-14135",
        label = "Hearing aid use (UKB)",
        role = "primary_offtarget"
    ),
    # Secondary off-targets — metabolic
    Fasting_Glucose = list(
        id = "ieu-b-113",
        label = "Fasting glucose (MAGIC)",
        role = "secondary_offtarget"
    ),
    HbA1c = list(
        id = "ieu-b-103",
        label = "HbA1c (MAGIC)",
        role = "secondary_offtarget"
    ),
    Type2_Diabetes = list(
        id = "ebi-a-GCST006867",
        label = "Type 2 diabetes (Xue 2018)",
        role = "secondary_offtarget"
    ),

    # Negative controls
    Height = list(
        id = "ieu-a-89",
        label = "Height (negative control)",
        role = "negative_control"
    ),
    Skin_Color = list(
        id = "ukb-b-19560",
        label = "Skin color (negative control)",
        role = "negative_control"
    )
)

cat("Design:\n")
cat(sprintf("  Exposure: IGF1R cis-eQTL (%s)\n\n", EXPOSURE_IGF1R))
cat("  Outcomes:\n")
for (nm in names(off_target_outcomes)) {
    oc <- off_target_outcomes[[nm]]
    cat(sprintf("    [%-20s] %s (%s)\n", oc$role, nm, oc$id))
}

# =============================================================================
# STEP 5.1 — Extract IGF1R IVs (reuse Step 2 logic)
# =============================================================================

cat("\n\n=== 5.1 Extracting IGF1R IVs ===\n")

igf_ivs <- ieugwasr::tophits(
    id = EXPOSURE_IGF1R, pval = 5e-8, clump = 1,
    r2 = 0.001, kb = 10000, opengwas_jwt = token
)

cat(sprintf("IGF1R IVs (strict, P<5e-8): %d SNPs\n", nrow(igf_ivs)))

# Convert to TwoSampleMR format
igf_exposure <- data.frame(
    SNP = igf_ivs$rsid,
    beta.exposure = igf_ivs$beta,
    se.exposure = igf_ivs$se,
    effect_allele.exposure = igf_ivs$ea,
    other_allele.exposure = igf_ivs$nea,
    eaf.exposure = igf_ivs$eaf,
    pval.exposure = igf_ivs$p,
    exposure = "IGF1R_expression",
    id.exposure = EXPOSURE_IGF1R,
    stringsAsFactors = FALSE
)
igf_exposure$F_stat <- igf_exposure$beta.exposure^2 / igf_exposure$se.exposure^2
igf_exposure <- igf_exposure[igf_exposure$F_stat >= 10, ]
cat(sprintf("After F≥10 filter: %d SNPs\n", nrow(igf_exposure)))
cat(sprintf("Mean F-statistic: %.1f\n\n", mean(igf_exposure$F_stat)))

# =============================================================================
# STEP 5.2 — MR across all outcomes
# =============================================================================

cat("\n=== 5.2 Running MR across outcomes ===\n")

all_mr_results <- list()
all_sens_results <- list()

for (oc_name in names(off_target_outcomes)) {
    oc <- off_target_outcomes[[oc_name]]
    cat(sprintf("\n%s\n", paste(rep("-", 60), collapse = "")))
    cat(sprintf("%s → %s [%s]\n", oc_name, oc$label, oc$role))
    cat(paste(rep("-", 60), collapse = ""), "\n")

    # Extract outcome
    out_data <- tryCatch(
        extract_outcome_data(snps = igf_exposure$SNP, outcomes = oc$id),
        error = function(e) {
            cat(sprintf("  [ERROR] %s\n", e$message))
            NULL
        }
    )

    if (is.null(out_data) || nrow(out_data) == 0) {
        cat("  [SKIP] No outcome data\n")
        next
    }

    # Harmonise
    dat <- tryCatch(
        harmonise_data(igf_exposure, out_data, action = 2),
        error = function(e) NULL
    )
    if (is.null(dat)) {
        cat("  [SKIP] Harmonise failed\n")
        next
    }
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0) {
        cat("  [SKIP] No SNPs after harmonise\n")
        next
    }

    cat(sprintf("  SNPs post-harmonise: %d\n", nrow(dat)))

    # Run MR
    methods <- if (nrow(dat) >= 3) {
        c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
    } else if (nrow(dat) == 2) {
        c("mr_ivw", "mr_wald_ratio")
    } else {
        "mr_wald_ratio"
    }

    mr_res <- tryCatch(mr(dat, method_list = methods), error = function(e) NULL)

    if (!is.null(mr_res)) {
        mr_res$outcome_name <- oc_name
        mr_res$outcome_role <- oc$role
        all_mr_results[[oc_name]] <- mr_res

        for (j in seq_len(nrow(mr_res))) {
            cat(sprintf(
                "    [%-25s] n=%d  β=%+.4f (SE=%.4f)  P=%.3e  %s\n",
                mr_res$method[j], mr_res$nsnp[j],
                mr_res$b[j], mr_res$se[j], mr_res$pval[j],
                ifelse(mr_res$pval[j] < 0.001, "***",
                    ifelse(mr_res$pval[j] < 0.05, "*", "")
                )
            ))
        }
    }

    # Sensitivity
    if (nrow(dat) >= 3) {
        pt <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
        ht <- tryCatch(mr_heterogeneity(dat, method_list = "mr_ivw"), error = function(e) NULL)

        sens_row <- data.frame(
            outcome_name = oc_name,
            outcome_role = oc$role,
            n_iv = nrow(dat),
            Egger_intercept = if (!is.null(pt)) pt$egger_intercept[1] else NA,
            Egger_int_p = if (!is.null(pt)) pt$pval[1] else NA,
            Cochran_Q = if (!is.null(ht)) ht$Q[1] else NA,
            Cochran_Q_p = if (!is.null(ht)) ht$Q_pval[1] else NA,
            stringsAsFactors = FALSE
        )
        all_sens_results[[oc_name]] <- sens_row
    }

    Sys.sleep(1)
}

# =============================================================================
# STEP 5.3 — Consolidate
# =============================================================================

cat("\n\n=== 5.3 Consolidating results ===\n")

mr_all <- if (length(all_mr_results) > 0) do.call(rbind, all_mr_results) else NULL
sens_all <- if (length(all_sens_results) > 0) do.call(rbind, all_sens_results) else NULL

if (!is.null(mr_all)) {
    write.csv(mr_all, "TrackA_MR/results/05_cross_phenotype_mr.csv", row.names = FALSE)
}
if (!is.null(sens_all)) {
    write.csv(sens_all, "TrackA_MR/results/05_cross_phenotype_sensitivity.csv", row.names = FALSE)
}

# =============================================================================
# STEP 5.4 — Summary table
# =============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("SE2 CROSS-PHENOTYPE MR RESULTS (IGF1R → Off-Target Phenotypes)\n")
cat(rep("=", 70), "\n\n", sep = "")

if (!is.null(mr_all)) {
    # Focus on IVW (and Wald if <3 IVs)
    primary_method_rows <- mr_all[mr_all$method %in% c("Inverse variance weighted", "Wald ratio"), ]

    summary_table <- data.frame(
        Outcome = primary_method_rows$outcome_name,
        Role = primary_method_rows$outcome_role,
        N_IV = primary_method_rows$nsnp,
        Method = primary_method_rows$method,
        Beta = round(primary_method_rows$b, 4),
        SE = round(primary_method_rows$se, 4),
        P = signif(primary_method_rows$pval, 3),
        Signif = ifelse(primary_method_rows$pval < 0.001, "***",
            ifelse(primary_method_rows$pval < 0.05, "*", "")
        ),
        stringsAsFactors = FALSE
    )

    # Order by role then P
    summary_table <- summary_table[order(summary_table$Role, summary_table$P), ]
    print(summary_table, row.names = FALSE)
}

# =============================================================================
# STEP 5.5 — SE2 Verdict
# =============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("SE2 VERDICT — IGF1R Off-Target Cross-Phenotype Evidence\n")
cat(rep("=", 70), "\n\n", sep = "")

if (!is.null(mr_all)) {
    # Check hearing outcomes
    hearing_sig <- any(
        mr_all$method == "Inverse variance weighted" &
            grepl("Hearing", mr_all$outcome_name) &
            mr_all$pval < 0.05
    )

    # Check glucose outcomes
    glucose_sig <- any(
        mr_all$method == "Inverse variance weighted" &
            (grepl("Glucose", mr_all$outcome_name) | grepl("HbA1c", mr_all$outcome_name) |
                grepl("Diabetes", mr_all$outcome_name)) &
            mr_all$pval < 0.05
    )

    # Check negative controls
    neg_control_sig <- any(
        mr_all$method == "Inverse variance weighted" &
            (grepl("Height", mr_all$outcome_name) | grepl("Skin", mr_all$outcome_name)) &
            mr_all$pval < 0.05
    )

    cat(
        "Primary off-target (hearing):     ",
        ifelse(hearing_sig, "✅ Significant signal", "⚪ Null"), "\n"
    )
    cat(
        "Secondary off-target (metabolic): ",
        ifelse(glucose_sig, "✅ Significant signal", "⚪ Null"), "\n"
    )
    cat(
        "Negative controls (specificity):  ",
        ifelse(!neg_control_sig, "✅ Null (good specificity)", "⚠️  Non-null (confounding risk)"), "\n"
    )

    cat("\n--- SE2 Overall Interpretation ---\n")
    if ((hearing_sig || glucose_sig) && !neg_control_sig) {
        cat("✅ SE2 ACHIEVED: IGF1R has specific genetic effects on\n")
        cat("   teprotumumab's known off-target phenotypes, but not on\n")
        cat("   non-specific negative controls.\n")
        cat("\n   Clinical implication: IGF1R genetic signal supports a\n")
        cat("   mechanistic basis for observed adverse events, which would be\n")
        cat("   expected to be enhanced by therapeutic IGF-1R blockade.\n")
    } else if (hearing_sig || glucose_sig) {
        cat("🟡 SE2 PARTIAL: Off-target signal detected, but negative controls\n")
        cat("   also show non-null effects — interpret with caution for specificity.\n")
    } else {
        cat("⚪ SE2 NULL: No significant genetic effect of IGF1R on off-target\n")
        cat("   phenotypes at population level. Does NOT preclude therapeutic\n")
        cat("   off-target effects, which reflect acute pharmacologic blockade.\n")
    }
}

sink()
cat("\n✅ Step 5 complete.\n")
cat("   Results: TrackA_MR/results/05_cross_phenotype_mr.csv\n")
cat("   Sensitivity: TrackA_MR/results/05_cross_phenotype_sensitivity.csv\n\n")
