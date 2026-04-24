# TED-TRAP Phase A-1: MR Re-execution for TSHR and IGF1R
# Purpose: Reproduce TSHR and IGF1R MR results across all 3 GWAS outcomes
# Expected values to reproduce:
#   TSHR Primary: beta=-1.394, P=6.79e-17, nSNP=2
#   TSHR Replication: beta=-0.012, P=3.88e-29, nSNP=5
#   IGF1R Primary: beta=0.217, P=5.27e-2, nSNP=11
#   IGF1R Replication: beta=0.001, P=4.12e-2, nSNP=11

library(TwoSampleMR)

# === Exposure data (eQTLGen cis-eQTLs) ===
genes <- list(
    TSHR = "eqtl-a-ENSG00000165409",
    IGF1R = "eqtl-a-ENSG00000140443"
)

# === Outcome data ===
outcomes <- list(
    Primary = "ebi-a-GCST90018627", # Graves disease
    Replication = "ebi-a-GCST90038636", # Hyperthyroidism
    FinnGen = "DiZ67P" # FinnGen Graves ophthalmopathy (if loaded)
)

results_all <- data.frame()
sensitivity_all <- data.frame()

for (gene_name in names(genes)) {
    cat("\n========================================\n")
    cat("Processing:", gene_name, "\n")
    cat("========================================\n")

    # Extract instruments
    exp <- extract_instruments(outcomes = genes[[gene_name]])

    if (is.null(exp) || nrow(exp) == 0) {
        cat("No instruments found for", gene_name, "\n")
        next
    }
    cat("Instruments found:", nrow(exp), "\n")

    for (out_name in names(outcomes)) {
        cat("\n--- ", gene_name, " vs ", out_name, " (", outcomes[[out_name]], ") ---\n", sep = "")

        # Extract outcome data
        out_data <- extract_outcome_data(
            snps = exp$SNP,
            outcomes = outcomes[[out_name]]
        )

        if (is.null(out_data) || nrow(out_data) == 0) {
            cat("No outcome data found\n")
            next
        }

        # Harmonize
        dat <- harmonise_data(exp, out_data)
        dat <- dat[dat$mr_keep, ]
        cat("Harmonized SNPs:", nrow(dat), "\n")

        if (nrow(dat) == 0) {
            cat("No harmonized SNPs\n")
            next
        }

        # Run MR (all methods)
        mr_res <- mr(dat, method_list = c(
            "mr_ivw", "mr_egger_regression",
            "mr_weighted_median", "mr_weighted_mode"
        ))

        cat("\nMR Results:\n")
        print(mr_res[, c("method", "nsnp", "b", "se", "pval")])

        # Add gene and outcome info
        mr_res$gene <- gene_name
        mr_res$outcome_label <- out_name
        results_all <- rbind(results_all, mr_res)

        # Sensitivity analyses
        sens_row <- data.frame(
            gene = gene_name,
            outcome = out_name,
            nsnp = nrow(dat),
            stringsAsFactors = FALSE
        )

        # Egger intercept
        if (nrow(dat) >= 3) {
            egger <- mr_pleiotropy_test(dat)
            sens_row$egger_intercept <- egger$egger_intercept
            sens_row$egger_intercept_se <- egger$se
            sens_row$egger_intercept_pval <- egger$pval
            cat("Egger intercept P:", egger$pval, "\n")

            # Cochran Q
            het <- mr_heterogeneity(dat)
            ivw_het <- het[het$method == "Inverse variance weighted", ]
            if (nrow(ivw_het) > 0) {
                sens_row$cochran_Q <- ivw_het$Q
                sens_row$cochran_Q_pval <- ivw_het$Q_pval
                cat("Cochran Q P:", ivw_het$Q_pval, "\n")
            }
        } else {
            cat("nSNP <3, cannot run pleiotropy/heterogeneity\n")
            sens_row$egger_intercept <- NA
            sens_row$egger_intercept_se <- NA
            sens_row$egger_intercept_pval <- NA
            sens_row$cochran_Q <- NA
            sens_row$cochran_Q_pval <- NA
        }

        # Steiger directionality
        steiger <- directionality_test(dat)
        sens_row$steiger_correct <- steiger$correct_causal_direction
        sens_row$steiger_pval <- steiger$steiger_pval
        cat(
            "Steiger direction:", steiger$correct_causal_direction,
            "P:", steiger$steiger_pval, "\n"
        )

        sensitivity_all <- rbind(sensitivity_all, sens_row)
    }
}

# Save results
outdir <- "c:/ProjectTEDGWAS/TrackA_MR/results"
write.csv(results_all, file.path(outdir, "MR_TSHR_IGF1R_reexecution.csv"), row.names = FALSE)
write.csv(sensitivity_all, file.path(outdir, "MR_TSHR_IGF1R_sensitivity.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("Results saved to:\n")
cat("  MR_TSHR_IGF1R_reexecution.csv\n")
cat("  MR_TSHR_IGF1R_sensitivity.csv\n")
cat("========================================\n")

# Print summary for verification
cat("\n========== VERIFICATION SUMMARY ==========\n")
ivw <- results_all[results_all$method == "Inverse variance weighted", ]
for (i in 1:nrow(ivw)) {
    cat(sprintf(
        "%s vs %s: beta=%.6f, SE=%.6f, P=%.2e, nSNP=%d\n",
        ivw$gene[i], ivw$outcome_label[i],
        ivw$b[i], ivw$se[i], ivw$pval[i], ivw$nsnp[i]
    ))
}
