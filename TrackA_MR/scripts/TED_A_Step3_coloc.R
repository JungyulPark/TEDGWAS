library(TwoSampleMR)
library(ieugwasr)
library(coloc)
library(data.table)

cat("\n=== 1. TSHR IV SNPs ===\n")
iv <- extract_instruments("eqtl-a-ENSG00000165409", p1 = 5e-6)
print(iv[, c("SNP", "chr.exposure", "pos.exposure", "pval.exposure")])

cat("\n=== 2. Colocalization ===\n")
cat("Fetching GWAS region data...\n")
gwas_reg <- tryCatch(associations(id = "ebi-a-GCST90018627", chr = "14", start = 80500000, end = 81700000), error = function(e) NULL)

if (is.null(gwas_reg) || nrow(gwas_reg) == 0) {
    cat("Failed to fetch GWAS region.\n")
} else {
    cat("Loading eQTLGen full data for TSHR...\n")
    eqtl_all <- fread("c:/ProjectTEDGWAS/TrackA_MR/data/eqtl_TSHR_all.csv")

    shared_snps <- intersect(gwas_reg$rsid, eqtl_all$SNP)
    cat("Shared SNPs: ", length(shared_snps), "\n")

    if (length(shared_snps) > 0) {
        gwas_sub <- gwas_reg[gwas_reg$rsid %in% shared_snps, ]
        eqtl_sub <- eqtl_all[eqtl_all$SNP %in% shared_snps, ]

        gwas_sub <- gwas_sub[!duplicated(gwas_sub$rsid), ]
        eqtl_sub <- eqtl_sub[!duplicated(eqtl_sub$SNP), ]

        merged <- merge(gwas_sub, eqtl_sub, by.x = "rsid", by.y = "SNP")

        # Calculate MAF from GWAS EAF
        maf <- ifelse(merged$eaf < 0.5, merged$eaf, 1 - merged$eaf)

        # Handle NAs in MAF
        valid <- !is.na(maf) & !is.na(merged$Pvalue) & !is.na(merged$p) & maf > 0 & maf < 0.5
        merged <- merged[valid, ]
        maf <- maf[valid]
        cat("Valid SNPs for coloc: ", nrow(merged), "\n")

        dataset1 <- list(pvalues = merged$Pvalue, N = 31684, type = "quant", snp = merged$rsid, MAF = maf)
        dataset2 <- list(pvalues = merged$p, N = 175465, type = "cc", s = 0.016, snp = merged$rsid, MAF = maf)

        res <- coloc.abf(dataset1, dataset2)
        cat("\n=== Coloc Summary ===\n")
        print(res$summary)
    }
}
