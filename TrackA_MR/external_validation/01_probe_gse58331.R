library(GEOquery)
library(Biobase)

setwd("c:/ProjectTEDGWAS/TrackA_MR/external_validation/")
gse <- getGEO("GSE58331", GSEMatrix = TRUE, AnnotGPL = TRUE)
eset <- gse[[1]]

pdata <- pData(eset)
cat("=== colnames(pdata) ===\n")
print(colnames(pdata))

cat("\n=== first 3 rows of specific columns ===\n")
cols_of_interest <- c("title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")
print(head(pdata[, intersect(colnames(pdata), cols_of_interest)], 5))

cat("\n=== Unique values in characteristics ===\n")
if("characteristics_ch1" %in% colnames(pdata)) print(table(pdata$characteristics_ch1))
if("characteristics_ch1.1" %in% colnames(pdata)) print(table(pdata$characteristics_ch1.1))
if("characteristics_ch1.2" %in% colnames(pdata)) print(table(pdata$characteristics_ch1.2))
if("source_name_ch1" %in% colnames(pdata)) print(table(pdata$source_name_ch1))
