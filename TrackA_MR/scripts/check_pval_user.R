library(data.table)
gwas_data <- fread("c:/ProjectTEDGWAS/TrackA_MR/data/finngen_TSHR_region.csv")

# Rename columns safely for subset
names(gwas_data)[1] <- "CHR"
names(gwas_data)[2] <- "POS"

gwas_region <- subset(gwas_data, CHR == 14 & POS > 80500000 & POS < 81700000)

cat("min(gwas_region$pval)  => ", min(gwas_region$pval, na.rm = TRUE), "\n")
cat("nrow(gwas_region)      => ", nrow(gwas_region), "\n")

cat("\nTop 5 SNPs by p-value:\n")
print(head(gwas_region[order(pval), c("CHR", "POS", "rsids", "pval", "beta")], 5))
