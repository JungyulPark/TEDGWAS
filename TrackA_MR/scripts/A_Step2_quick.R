library(TwoSampleMR)

genes <- list(
    IGF1R = "eqtl-a-ENSG00000140443",
    TSHR  = "eqtl-a-ENSG00000165409",
    CTLA4 = "eqtl-a-ENSG00000163599",
    IL6   = "eqtl-a-ENSG00000136244"
)

run_mr <- function(outcome_id, label) {
    cat(sprintf("\n\n======== RUNNING %s (%s) ========\n", label, outcome_id))
    for (gene_name in c("IGF1R", "TSHR")) {
        cat(sprintf("\n--- %s ---\n", gene_name))
        exposure <- tryCatch(extract_instruments(outcomes = genes[[gene_name]], p1 = 5e-6, clump = TRUE, r2 = 0.001, kb = 10000), error = function(e) NULL)
        if (is.null(exposure) || nrow(exposure) == 0) {
            cat("  ❌ Insufficient IVs\n")
            next
        }
        cat("  IVs extracted: ", nrow(exposure), "\n")
        Sys.sleep(1)

        outcome <- tryCatch(extract_outcome_data(snps = exposure$SNP, outcomes = outcome_id), error = function(e) NULL)
        if (is.null(outcome) || nrow(outcome) == 0) {
            cat("  ❌ No outcome data\n")
            next
        }
        cat("  Outcomes extracted: ", nrow(outcome), "\n")

        harmonised <- harmonise_data(exposure, outcome)
        if (sum(harmonised$mr_keep) == 0) {
            cat("  ❌ No harmonised SNPs\n")
            next
        }

        mr_result <- tryCatch(mr(harmonised, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")), error = function(e) NULL)

        if (!is.null(mr_result)) {
            ivw <- mr_result[mr_result$method == "Inverse variance weighted", ]
            if (nrow(ivw) > 0) {
                cat(sprintf("  [%s] IVW: β=%.3f, SE=%.3f, P=%.2e\n", gene_name, ivw$b, ivw$se, ivw$pval))
            }
        }
        Sys.sleep(2)
    }
}

run_mr("ebi-a-GCST90018627", "PRIMARY (Graves)")
run_mr("ebi-a-GCST90038636", "REPLICATION (Hyperthyroidism)")
