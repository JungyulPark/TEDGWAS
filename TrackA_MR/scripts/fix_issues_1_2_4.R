library(TwoSampleMR)

cat("=== TSHR IV 정확한 수 ===\n")
tshr_exp <- tryCatch(extract_instruments("eqtl-a-ENSG00000165409", p1 = 5e-6, clump = TRUE), error = function(e) NULL)
if (!is.null(tshr_exp)) {
    cat("Exposure SNPs (pre-harmonize):", nrow(tshr_exp), "\n")
    cat("SNPs:", paste(tshr_exp$SNP, collapse = ", "), "\n")
    out_dat <- tryCatch(extract_outcome_data(tshr_exp$SNP, "ebi-a-GCST90018627"), error = function(e) NULL)
    if (!is.null(out_dat)) {
        h <- harmonise_data(tshr_exp, out_dat)
        cat("After harmonize (mr_keep=TRUE):", sum(h$mr_keep), "\n")
    }
}

cat("\n=== TNF Sensitivity ===\n")
tnf_exp <- tryCatch(extract_instruments("eqtl-a-ENSG00000232810", p1 = 5e-6, clump = TRUE), error = function(e) NULL)
out_rep <- tryCatch(extract_outcome_data(tnf_exp$SNP, "ebi-a-GCST90038636"), error = function(e) NULL)
if (!is.null(out_rep)) {
    h_tnf <- harmonise_data(tnf_exp, out_rep)
    if (sum(h_tnf$mr_keep) >= 3) {
        het <- mr_heterogeneity(h_tnf)
        print(het[, c("method", "Q", "Q_df", "Q_pval")])
        pleio <- mr_pleiotropy_test(h_tnf)
        cat("Egger intercept P:", pleio$pval, "\n")
    }
}

cat("\n=== LD r2 via OpenGWAS API ===\n")
# Use OpenGWAS LD endpoint
token <- Sys.getenv("OPENGWAS_JWT")
snps <- c("rs3783947", "rs179252")
r <- tryCatch(
    {
        res <- ieugwasr::ld_matrix(snps, pop = "EUR", r2 = TRUE, with_alleles = FALSE)
        print(res)
    },
    error = function(e) cat("LD API error:", conditionMessage(e), "\n")
)
