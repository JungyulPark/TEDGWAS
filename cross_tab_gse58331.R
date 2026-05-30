library(GEOquery)
gset <- getGEO("GSE58331", GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
eset <- gset[[1]]
pd <- pData(eset)
sink("gse58331_cross_tab.txt")
print(table(pd[["diagnosis:ch1"]], pd[["tissue:ch1"]]))
sink()
cat("Done cross tab!\n")
