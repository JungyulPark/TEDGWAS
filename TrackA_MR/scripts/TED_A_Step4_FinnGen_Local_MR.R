library(TwoSampleMR)
library(data.table)

genes <- list(TSHR = "eqtl-a-ENSG00000165409", IGF1R = "eqtl-a-ENSG00000140443", TNF = "eqtl-a-ENSG00000232810")

cat("=== 1. Fetching Exposures ===\n")
exp_dat <- data.frame()
for (g in names(genes)) {
    e <- tryCatch(extract_instruments(genes[[g]], p1 = 5e-6, clump = TRUE), error = function(e) NULL)
    if (!is.null(e) && nrow(e) > 0) {
        e$gene <- g
        exp_dat <- rbind(exp_dat, e)
    }
}
snps_to_keep <- unique(exp_dat$SNP)
cat("Total unique IV SNPs to query:", length(snps_to_keep), "\n")

cat("=== 2. Loading Local FinnGen GRAVES_OPHT ===\n")
finn <- fread("c:/ProjectTEDGWAS/TrackA_MR/data/finngen_R12_GRAVES_OPHT.gz")

finn_sub <- finn[rsids %in% snps_to_keep]
cat("Found", nrow(finn_sub), "SNPs in FinnGen outcome.\n")

if (nrow(finn_sub) > 0) {
    out_dat <- format_data(
        as.data.frame(finn_sub),
        type = "outcome",
        snp_col = "rsids",
        beta_col = "beta",
        se_col = "sebeta",
        eaf_col = "af_alt",
        effect_allele_col = "alt",
        other_allele_col = "ref",
        pval_col = "pval"
    )

    cat("=== 3. Running MR ===\n")
    res_df <- data.frame()
    for (g in unique(exp_dat$gene)) {
        e_sub <- exp_dat[exp_dat$gene == g, ]
        h <- harmonise_data(e_sub, out_dat)
        if (sum(h$mr_keep) > 0) {
            mr_res <- mr(h, method_list = c("mr_ivw"))
            if (nrow(mr_res) > 0) {
                mr_res$n_iv <- sum(h$mr_keep)
                mr_res$gene <- g
                res_df <- rbind(res_df, mr_res)
            }
        }
    }
    print(res_df)
    write.csv(res_df, "c:/ProjectTEDGWAS/TrackA_MR/results/MR_FinnGen_Local_Summary.csv", row.names = FALSE)
} else {
    cat("No requested SNPs found in FinnGen data.\n")
}
