library(GEOquery)
gset <- getGEO("GSE58331", GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
eset <- gset[[1]]
pd <- pData(eset)
sink("gse58331_pdata_summary.txt")
cat("=== Column names ===\n")
print(colnames(pd))
cat("\n=== Sample summary table ===\n")
# Print characteristics columns
for(c in colnames(pd)){
  if(grepl("characteristic|source|title|disease|diagnos|group|status|tissue", c, ignore.case=TRUE)){
    cat("\nColumn:", c, "\n")
    print(table(pd[[c]], useNA="ifany"))
  }
}
sink()
cat("Done!\n")
