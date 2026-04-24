# =============================================================================
# TED-TRAP Step 2 — Core MR Execution
# Primary: Graves disease (ebi-a-GCST90018627)
# Replication: Hyperthyroidism (ebi-a-GCST90038636)
# Sensitivity: FinnGen R12 Graves ophthalmopathy
# =============================================================================

library(TwoSampleMR)
library(ieugwasr)

# --- Setup ---
token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(token) < 100) stop("Token not loaded.")


# --- Gene config from Step 1b verdict ---
gene_config <- list(
    TSHR  = list(id = "eqtl-a-ENSG00000165409", mode = "wald", pval = 5e-8),
    IGF1R = list(id = "eqtl-a-ENSG00000140443", mode = "full", pval = 5e-8),
    TNF   = list(id = "eqtl-a-ENSG00000232810", mode = "no_presso", pval = 5e-8),
    PPARG = list(id = "eqtl-a-ENSG00000132170", mode = "no_presso", pval = 5e-6),
    ARRB1 = list(id = "eqtl-a-ENSG00000137486", mode = "full", pval = 5e-8),
    IRS1  = list(id = "eqtl-a-ENSG00000169047", mode = "no_presso", pval = 5e-8),
    AKT1  = list(id = "eqtl-a-ENSG00000142208", mode = "no_presso", pval = 5e-8),
    CTLA4 = list(id = "eqtl-a-ENSG00000163599", mode = "no_presso", pval = 5e-8)
)

# --- Outcomes ---
outcomes <- list(
    Primary     = list(id = "ebi-a-GCST90018627", label = "Graves disease"),
    Replication = list(id = "ebi-a-GCST90038636", label = "Hyperthyroidism"),
    Sensitivity = list(id = "finn-b-E4_GRAVES_OPHT", label = "FinnGen Graves ophthalmopathy")
)

# --- Storage ---
all_mr <- list()
all_sens <- list()
all_wald <- list()
all_presso <- list()

# --- Setup directories ---
for (d in c(
    "TrackA_MR/data/instruments", "TrackA_MR/results",
    "TrackA_MR/logs", "TrackA_MR/figures"
)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

log_file <- "TrackA_MR/logs/02_core_mr.log"
sink(log_file, split = TRUE)

cat("=== Step 2 — Core MR Execution ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# =============================================================================
# MAIN LOOP
# =============================================================================

for (g in names(gene_config)) {
    cfg <- gene_config[[g]]
    cat(sprintf(
        "\n%s\nGENE: %s (mode: %s, P<%g)\n%s\n",
        paste(rep("=", 60), collapse = ""),
        g, cfg$mode, cfg$pval,
        paste(rep("=", 60), collapse = "")
    ))

    # --- Extract exposure IVs ---
    cat(sprintf("→ Extracting IVs from %s\n", cfg$id))
    exp_data <- tryCatch(
        {
            raw <- ieugwasr::tophits(id = cfg$id, pval = cfg$pval, clump = TRUE, r2 = 0.001, kb = 10000, opengwas_jwt = token)
            if (is.null(raw) || nrow(raw) == 0) {
                NULL
            } else {
                suppressMessages(TwoSampleMR::format_data(raw, type = "exposure", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "p", effect_allele_col = "ea", other_allele_col = "nea", eaf_col = "eaf", id_col = "id"))
            }
        },
        error = function(e) {
            cat(sprintf("  [ERROR] %s\n", e$message))
            NULL
        }
    )

    if (is.null(exp_data) || nrow(exp_data) == 0) {
        cat("  [SKIP] No IVs retrieved\n")
        next
    }

    # F-statistic
    exp_data$F_stat <- (exp_data$beta.exposure)^2 / (exp_data$se.exposure)^2
    exp_data <- exp_data[exp_data$F_stat >= 10, ]
    cat(sprintf(
        "  %d IVs retained (F≥10, min F = %.1f)\n",
        nrow(exp_data), min(exp_data$F_stat)
    ))

    # Save instruments
    iv_file <- sprintf("TrackA_MR/data/instruments/%s_instruments.tsv", g)
    write.table(exp_data, iv_file, sep = "\t", row.names = FALSE, quote = FALSE)

    # --- Loop across outcomes ---
    for (oc_name in names(outcomes)) {
        oc <- outcomes[[oc_name]]
        cat(sprintf("\n  Outcome: %s [%s]\n", oc$label, oc_name))

        # Extract outcome data
        out_data <- tryCatch(
            {
                raw_out <- ieugwasr::associations(variants = exp_data$SNP, id = oc$id, opengwas_jwt = token)
                if (is.null(raw_out) || nrow(raw_out) == 0) {
                    NULL
                } else {
                    suppressMessages(TwoSampleMR::format_data(raw_out, type = "outcome", snp_col = "rsid", beta_col = "beta", se_col = "se", pval_col = "p", effect_allele_col = "ea", other_allele_col = "nea", eaf_col = "eaf", id_col = "id"))
                }
            },
            error = function(e) {
                cat(sprintf("    [ERROR outcome] %s\n", e$message))
                NULL
            }
        )

        if (is.null(out_data) || nrow(out_data) == 0) {
            cat("    [SKIP] No outcome data\n")
            next
        }

        # Harmonise
        dat <- tryCatch(suppressMessages(harmonise_data(exp_data, out_data, action = 2)),
            error = function(e) NULL
        )
        if (is.null(dat)) {
            cat("    [SKIP] Harmonise failed\n")
            next
        }
        dat <- dat[dat$mr_keep, ]
        if (nrow(dat) == 0) {
            cat("    [SKIP] No valid SNPs after harmonise\n")
            next
        }

        cat(sprintf("    SNPs post-harmonise: %d\n", nrow(dat)))

        # ========== MR methods by mode ==========

        # Determine methods to run
        methods_list <- if (nrow(dat) == 1) {
            "mr_wald_ratio"
        } else if (nrow(dat) == 2) {
            c("mr_ivw")
        } else if (nrow(dat) >= 3) {
            if (cfg$mode == "full" || nrow(dat) >= 4) {
                c(
                    "mr_ivw", "mr_egger_regression", "mr_weighted_median",
                    "mr_weighted_mode", "mr_simple_mode"
                )
            } else {
                c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
            }
        } else {
            c("mr_ivw")
        }

        # Run MR
        mr_res <- tryCatch(suppressMessages(mr(dat, method_list = methods_list)),
            error = function(e) {
                cat(sprintf("    [ERROR mr()] %s\n", e$message))
                NULL
            }
        )

        if (!is.null(mr_res) && nrow(mr_res) > 0) {
            mr_res$gene <- g
            mr_res$outcome_role <- oc_name
            mr_res$n_iv <- nrow(dat)
            all_mr[[paste(g, oc_name, sep = "_")]] <- mr_res

            # Compact print
            for (j in seq_len(nrow(mr_res))) {
                cat(sprintf(
                    "      [%-25s] n=%d  β=%+.4f (SE=%.4f)  P=%.3e\n",
                    mr_res$method[j], mr_res$nsnp[j],
                    mr_res$b[j], mr_res$se[j], mr_res$pval[j]
                ))
            }
        }

        # ========== Per-SNP Wald ratios (for TSHR and all <3 IV cases) ==========
        if (nrow(dat) <= 3) {
            wald_snp <- data.frame(
                gene = g, outcome_role = oc_name, SNP = dat$SNP,
                beta_exp = dat$beta.exposure, se_exp = dat$se.exposure,
                beta_out = dat$beta.outcome, se_out = dat$se.outcome,
                wald_beta = dat$beta.outcome / dat$beta.exposure,
                stringsAsFactors = FALSE
            )
            wald_snp$wald_se <- wald_snp$se_out / abs(wald_snp$beta_exp)
            wald_snp$wald_p <- 2 * pnorm(-abs(wald_snp$wald_beta / wald_snp$wald_se))
            all_wald[[paste(g, oc_name, sep = "_")]] <- wald_snp
        }

        # ========== Sensitivity: heterogeneity, pleiotropy, Steiger ==========
        sens_row <- data.frame(
            gene = g, outcome_role = oc_name, n_iv = nrow(dat),
            Egger_intercept = NA_real_, Egger_intercept_se = NA_real_, Egger_intercept_p = NA_real_,
            Cochran_Q = NA_real_, Cochran_Q_p = NA_real_,
            Steiger_correct = NA, Steiger_p = NA_real_,
            stringsAsFactors = FALSE
        )

        if (nrow(dat) >= 3) {
            pt <- tryCatch(suppressMessages(mr_pleiotropy_test(dat)), error = function(e) NULL)
            if (!is.null(pt) && nrow(pt) > 0) {
                sens_row$Egger_intercept <- pt$egger_intercept
                sens_row$Egger_intercept_se <- pt$se
                sens_row$Egger_intercept_p <- pt$pval
            }
            ht <- tryCatch(suppressMessages(mr_heterogeneity(dat, method_list = "mr_ivw")), error = function(e) NULL)
            if (!is.null(ht) && nrow(ht) > 0) {
                sens_row$Cochran_Q <- ht$Q
                sens_row$Cochran_Q_p <- ht$Q_pval
            }
        }
        st <- tryCatch(suppressMessages(directionality_test(dat)), error = function(e) NULL)
        if (!is.null(st) && nrow(st) > 0) {
            sens_row$Steiger_correct <- st$correct_causal_direction
            sens_row$Steiger_p <- st$steiger_pval
        }
        all_sens[[paste(g, oc_name, sep = "_")]] <- sens_row

        # ========== MR-PRESSO (only if cfg$mode == "full" AND n >= 4) ==========
        if (cfg$mode == "full" && nrow(dat) >= 4) {
            if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
                cat("    [INFO] MR-PRESSO not installed; skipping\n")
            } else {
                cat("    Running MR-PRESSO...\n")
                presso <- tryCatch(
                    suppressMessages(MRPRESSO::mr_presso(
                        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                        SdOutcome = "se.outcome", SdExposure = "se.exposure",
                        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                        data = dat, NbDistribution = 1000, SignifThreshold = 0.05
                    )),
                    error = function(e) {
                        cat(sprintf("      [ERROR PRESSO] %s\n", e$message))
                        NULL
                    }
                )
                if (!is.null(presso)) {
                    res <- presso$`Main MR results`
                    gp <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
                    out_idx <- presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
                    all_presso[[paste(g, oc_name, sep = "_")]] <- data.frame(
                        gene = g, outcome_role = oc_name, n_iv = nrow(dat),
                        global_test_p = gp,
                        n_outliers = if (is.null(out_idx)) 0L else length(out_idx),
                        raw_beta = res$`Causal Estimate`[1], raw_p = res$`P-value`[1],
                        corrected_beta = res$`Causal Estimate`[2], corrected_p = res$`P-value`[2],
                        stringsAsFactors = FALSE
                    )
                    cat(sprintf(
                        "      Global P=%.3e | Outliers: %d | Raw β=%+.4f | Corrected β=%+.4f\n",
                        gp, if (is.null(out_idx)) 0 else length(out_idx),
                        res$`Causal Estimate`[1], res$`Causal Estimate`[2]
                    ))
                }
            }
        }

        Sys.sleep(1) # rate-limit kindness
    } # end outcome loop
} # end gene loop

# =============================================================================
# CONSOLIDATE & SAVE
# =============================================================================

cat("\n\n=== Consolidating results ===\n")

mr_all <- if (length(all_mr) > 0) do.call(rbind, all_mr) else NULL
sens_all <- if (length(all_sens) > 0) do.call(rbind, all_sens) else NULL
wald_all <- if (length(all_wald) > 0) do.call(rbind, all_wald) else NULL
presso_all <- if (length(all_presso) > 0) do.call(rbind, all_presso) else NULL

if (!is.null(mr_all)) write.csv(mr_all, "TrackA_MR/results/02_mr_main.csv", row.names = FALSE)
if (!is.null(sens_all)) write.csv(sens_all, "TrackA_MR/results/02_mr_sensitivity.csv", row.names = FALSE)
if (!is.null(wald_all)) write.csv(wald_all, "TrackA_MR/results/02_mr_wald_perSNP.csv", row.names = FALSE)
if (!is.null(presso_all)) write.csv(presso_all, "TrackA_MR/results/02_mr_presso.csv", row.names = FALSE)

cat(sprintf("  MR main: %d rows\n", if (!is.null(mr_all)) nrow(mr_all) else 0))
cat(sprintf("  Sensitivity: %d rows\n", if (!is.null(sens_all)) nrow(sens_all) else 0))
cat(sprintf("  Wald per-SNP: %d rows\n", if (!is.null(wald_all)) nrow(wald_all) else 0))
cat(sprintf("  MR-PRESSO: %d rows\n", if (!is.null(presso_all)) nrow(presso_all) else 0))

# =============================================================================
# PRIMARY ENDPOINT — TSHR vs IGF1R for Graves disease
# =============================================================================

cat("\n\n", rep("=", 60), "\n", sep = "")
cat("PRIMARY ENDPOINT RESULTS (Graves disease)\n")
cat(rep("=", 60), "\n", sep = "")

if (!is.null(mr_all)) {
    primary <- mr_all[mr_all$outcome_role == "Primary" &
        mr_all$gene %in% c("TSHR", "IGF1R"), ]
    if (nrow(primary) > 0) {
        primary_subset <- primary[, c("gene", "method", "nsnp", "b", "se", "pval")]
        primary_subset$b <- round(primary_subset$b, 4)
        primary_subset$se <- round(primary_subset$se, 4)
        primary_subset$pval <- signif(primary_subset$pval, 3)
        print(primary_subset, row.names = FALSE)
    }
}

sink()

cat("\n✅ Step 2 Complete.\n")
cat("   Logs: TrackA_MR/logs/02_core_mr.log\n")
cat("   Results: TrackA_MR/results/02_mr_*.csv\n\n")
