library(GEOquery)
library(limma)
gset <- getGEO("GSE105149", GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
eset <- gset[[1]]
pd <- pData(eset)
pd_sub <- pd[pd[["diagnosis:ch1"]] %in% c("Normal", "TED"), ]
pd_sub$group <- ifelse(pd_sub[["diagnosis:ch1"]] == "TED", "TED", "Control")
expr_sub <- exprs(eset)[, rownames(pd_sub)]
group_factor <- factor(pd_sub$group, levels = c("Control", "TED"))
design <- model.matrix(~ group_factor)
colnames(design) <- c("Intercept", "TED_vs_Control")
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "TED_vs_Control", number = Inf, adjust.method = "BH")
results$probe_id <- rownames(results)
fdata <- fData(eset)
sym_col <- grep("Gene.?symbol|^Symbol$|GENE_SYMBOL", colnames(fdata), ignore.case = TRUE)[1]
results_merged <- merge(results, fdata[, c("ID", colnames(fdata)[sym_col])], by.x="probe_id", by.y="ID")
colnames(results_merged)[ncol(results_merged)] <- "gene"

sink("tshr_probes_gse105149.txt")
cat("=== TSHR probes ===\n")
print(results_merged[grepl("TSHR", results_merged$gene), ])
cat("\n=== IGF1R probes ===\n")
print(results_merged[grepl("IGF1R", results_merged$gene), ])
cat("\n=== CTLA4 probes ===\n")
print(results_merged[grepl("CTLA4", results_merged$gene), ])
sink()
cat("Done!\n")
