library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)

eset <- readRDS("C:/ProjectTEDGWAS/TrackA_MR/GEO_validation/GSE105149_eset.rds")
pd <- pData(eset)

# 4 TED vs 7 normal
pd_sub <- pd[pd[["diagnosis:ch1"]] %in% c("Normal", "TED"), ]
pd_sub$group <- ifelse(pd_sub[["diagnosis:ch1"]] == "TED", "TED", "Control")

expr_sub <- exprs(eset)[, rownames(pd_sub)]
if (max(expr_sub, na.rm = TRUE) > 50) {
  expr_sub <- log2(expr_sub + 1)
}

group_factor <- factor(pd_sub$group, levels = c("Control", "TED"))
design <- model.matrix(~ group_factor)
colnames(design) <- c("Intercept", "TED_vs_Control")

fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "TED_vs_Control", number = Inf, adjust.method = "BH")
res$probe_id <- rownames(res)

fdata <- fData(eset)
sym_col <- grep("Gene.?symbol|^Symbol$|GENE_SYMBOL", colnames(fdata), ignore.case = TRUE)[1]
probe_to_gene <- fdata[, c("ID", colnames(fdata)[sym_col]), drop = FALSE]
colnames(probe_to_gene) <- c("probe_id", "gene")
probe_to_gene$gene <- sapply(strsplit(as.character(probe_to_gene$gene), "///"), `[`, 1)
probe_to_gene$gene <- trimws(probe_to_gene$gene)

res_merged <- merge(res, probe_to_gene, by = "probe_id")

target_genes <- c("TSHR", "IGF1R", "INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1")
res_targets <- res_merged %>%
  filter(toupper(gene) %in% target_genes) %>%
  select(gene, probe_id, logFC, P.Value) %>%
  arrange(gene, P.Value)

print(head(res_targets, 100))
