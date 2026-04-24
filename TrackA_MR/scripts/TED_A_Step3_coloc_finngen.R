library(coloc)
library(data.table)

cat("=== Coloc: TSHR eQTLs vs FinnGen R12 Graves Disease ===\n")

cat("Loading extracted FinnGen region data...\n")
finn_sub <- tryCatch(fread("c:/ProjectTEDGWAS/TrackA_MR/data/finngen_TSHR_region.csv"), error = function(e) NULL)

if (is.null(finn_sub) || nrow(finn_sub) == 0) {
    cat("Failed to load FinnGen data.\n")
} else {
    cat("FinnGen SNPs in region: ", nrow(finn_sub), "\n")

    cat("Loading eQTLGen full data for TSHR...\n")
    eqtl_all <- fread("c:/ProjectTEDGWAS/TrackA_MR/data/eqtl_TSHR_all.csv")

    finn_rsid_col <- grep("rsid", names(finn_sub), value = TRUE, ignore.case = TRUE)[1]
    finn_pval_col <- grep("pval|p_value", names(finn_sub), value = TRUE, ignore.case = TRUE)[1]
    finn_af_col <- grep("af_alt|maf|eaf", names(finn_sub), value = TRUE, ignore.case = TRUE)[1]

    cat(sprintf("Using columns: rsid='%s', pval='%s', af='%s'\n", finn_rsid_col, finn_pval_col, finn_af_col))

    shared_snps <- intersect(finn_sub[[finn_rsid_col]], eqtl_all$SNP)
    cat("Shared SNPs: ", length(shared_snps), "\n")

    if (length(shared_snps) > 0) {
        gwas_sub <- finn_sub[finn_sub[[finn_rsid_col]] %in% shared_snps, ]
        eqtl_sub <- eqtl_all[eqtl_all$SNP %in% shared_snps, ]

        gwas_sub <- gwas_sub[!duplicated(gwas_sub[[finn_rsid_col]]), ]
        eqtl_sub <- eqtl_sub[!duplicated(eqtl_sub$SNP), ]

        merged <- merge(gwas_sub, eqtl_sub, by.x = finn_rsid_col, by.y = "SNP")

        maf <- merged[[finn_af_col]]
        maf <- ifelse(maf < 0.5, maf, 1 - maf)
        p_gwas <- merged[[finn_pval_col]]
        p_eqtl <- merged$Pvalue

        valid <- !is.na(maf) & !is.na(p_gwas) & !is.na(p_eqtl) & maf > 0 & maf < 0.5
        merged <- merged[valid, ]
        maf <- maf[valid]
        cat("Valid SNPs for coloc: ", nrow(merged), "\n")

        if (nrow(merged) > 0) {
            dataset1 <- list(pvalues = merged$Pvalue, N = 31684, type = "quant", snp = merged[[finn_rsid_col]], MAF = maf)
            dataset2 <- list(pvalues = merged[[finn_pval_col]], type = "cc", s = 0.016, N = 175465, snp = merged[[finn_rsid_col]], MAF = maf)

            res <- coloc.abf(dataset1, dataset2)
            cat("\n=== Coloc Summary ===\n")
            print(res$summary)
            saveRDS(res, "c:/ProjectTEDGWAS/TrackA_MR/results/coloc_finngen.rds")

            h4 <- res$summary["PP.H4.abf"]
            h3 <- res$summary["PP.H3.abf"]
            cat(sprintf("\nPP.H4 (Shared causal variant): %.3f\n", h4))
            cat(sprintf("PP.H3 (Distinct causal variants): %.3f\n", h3))
        }
    }
}
