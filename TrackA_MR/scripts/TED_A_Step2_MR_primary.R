library(TwoSampleMR)

# === 설정 ===
genes <- list(
    IGF1R = "eqtl-a-ENSG00000140443",
    TSHR  = "eqtl-a-ENSG00000165409",
    IL6   = "eqtl-a-ENSG00000136244",
    TGFB1 = "eqtl-a-ENSG00000105329",
    HAS2  = "eqtl-a-ENSG00000153446",
    PPARG = "eqtl-a-ENSG00000132170",
    ARRB1 = "eqtl-a-ENSG00000137486",
    IGF1  = "eqtl-a-ENSG00000017427",
    IRS1  = "eqtl-a-ENSG00000169047",
    AKT1  = "eqtl-a-ENSG00000142208",
    TNF   = "eqtl-a-ENSG00000232810",
    CTLA4 = "eqtl-a-ENSG00000163599"
)

outcome_id <- "ieu-a-1098"

# === MR 루프 ===
all_results <- list()

for (gene_name in names(genes)) {
    cat(sprintf("\n========== %s ==========\n", gene_name))

    # 1) IV 추출
    exposure <- tryCatch(
        extract_instruments(
            outcomes = genes[[gene_name]],
            p1 = 5e-6,
            clump = TRUE,
            r2 = 0.001,
            kb = 10000
        ),
        error = function(e) {
            cat("  ERROR:", e$message, "\n")
            NULL
        }
    )

    if (is.null(exposure) || nrow(exposure) == 0) {
        cat("  ❌ Insufficient IVs — SKIP\n")
        next
    }
    cat(sprintf("  IVs: %d\n", nrow(exposure)))

    # 2) Outcome extraction
    outcome <- tryCatch(
        extract_outcome_data(
            snps = exposure$SNP,
            outcomes = outcome_id
        ),
        error = function(e) {
            cat("  ERROR Outcome Extraction:", e$message, "\n")
            NULL
        }
    )

    if (is.null(outcome) || nrow(outcome) == 0) {
        cat("  ❌ No outcome data — SKIP\n")
        next
    }

    # 3) Harmonize
    harmonised <- harmonise_data(exposure, outcome)
    cat(sprintf("  Harmonised SNPs: %d\n", nrow(harmonised[harmonised$mr_keep, ])))

    if (sum(harmonised$mr_keep) == 0) {
        cat("  ❌ No harmonised SNPs — SKIP\n")
        next
    }

    # 4) MR
    mr_result <- tryCatch(mr(harmonised, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")), error = function(e) NULL)

    # 5) Steiger directionality
    steiger <- tryCatch(
        directionality_test(harmonised),
        error = function(e) {
            cat("  Steiger error:", e$message, "\n")
            NULL
        }
    )

    # 6) Heterogeneity & Pleiotropy
    het <- tryCatch(mr_heterogeneity(harmonised), error = function(e) NULL)
    pleio <- tryCatch(mr_pleiotropy_test(harmonised), error = function(e) NULL)

    # 8) Store results
    all_results[[gene_name]] <- list(mr = mr_result, steiger = steiger, het = het, pleio = pleio, n_iv = nrow(exposure), harmonised = harmonised)

    if (!is.null(mr_result)) {
        ivw <- mr_result[mr_result$method == "Inverse variance weighted", ]
        if (nrow(ivw) > 0) {
            cat(sprintf("  IVW: β=%.3f, SE=%.3f, P=%.2e\n", ivw$b, ivw$se, ivw$pval))
        }
    }
    Sys.sleep(3) # API rate limit
}

saveRDS(all_results, "c:/ProjectTEDGWAS/TrackA_MR/results/MR_all_results.rds")
cat("Finished MR loops.\n")
