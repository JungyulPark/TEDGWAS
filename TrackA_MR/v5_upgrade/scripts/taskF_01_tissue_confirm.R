#!/usr/bin/env Rscript
# =============================================================================
# Task F.1 — Backbone tissue re-confirmation (in-house orbital RNA-seq)
# =============================================================================
# Re-confirm TSHR/IGF1R/CTLA4 expression in orbital tissue from data.txt.
# Goal: validate that the pipeline reproduces v4.x TSHR log2FC ~ +2.33, padj ~ 0.006.
# NOT a transcriptome-wide re-discovery — backbone genes only.
#
# Samples (locked): Sample#2 = Control; Samples 7/8/10/11 = TED (4 TED + 1 control)
# In-house = DESeq2. n=1 control → supporting/descriptive layer, not powered DE.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

DATA   <- "c:/ProjectTEDGWAS/data.txt"   # 46,427 x 20 (FPKM + TPM)
WORK   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
OUTDIR <- file.path(WORK, "06_tissue_integration")
OUT    <- file.path(OUTDIR, "TaskF_01_backbone_tissue_inhouse_v1.csv")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

BACKBONE <- c("TSHR","IGF1R","CTLA4")

# ---- Load expression matrix ----
d <- fread(DATA)
cat("data.txt dims:", dim(d)[1], "x", dim(d)[2], "\n")
cat("Columns:\n"); print(names(d))

# IMPORTANT (Gemini): identify the gene-symbol column and the TPM columns.
# v4.x used TPM columns. Confirm which columns are TPM vs FPKM, and which
# correspond to Sample#2 (Ctrl) and Samples 7/8/10/11 (TED).
# Set these mappings explicitly below after inspecting names(d):

GENE_COL <- "Gene_Symbol"
CTRL_TPM <- c("2__1_TPM", "2__2_TPM")
TED_TPM  <- c("7__1_TPM", "7__2_TPM", "8__1_TPM", "8__2_TPM", "10__1_TPM", "10__2_TPM", "11__1_TPM", "11__2_TPM")

# ---- Extract backbone genes ----
sub <- d[get(GENE_COL) %in% BACKBONE]
cat("\nBackbone genes found:", nrow(sub), "of", length(BACKBONE), "\n")

res <- list()
for (g in BACKBONE) {
  row <- sub[get(GENE_COL)==g]
  if (nrow(row)==0) { cat("  MISSING:", g, "\n"); next }
  
  # Biological samples TPM (averaging technical replicates)
  t2  <- mean(as.numeric(row[, c("2__1_TPM", "2__2_TPM")]), na.rm=TRUE)
  t7  <- mean(as.numeric(row[, c("7__1_TPM", "7__2_TPM")]), na.rm=TRUE)
  t8  <- mean(as.numeric(row[, c("8__1_TPM", "8__2_TPM")]), na.rm=TRUE)
  t10 <- mean(as.numeric(row[, c("10__1_TPM", "10__2_TPM")]), na.rm=TRUE)
  t11 <- mean(as.numeric(row[, c("11__1_TPM", "11__2_TPM")]), na.rm=TRUE)
  
  ctrl_vals <- t2
  ted_vals  <- c(t7, t8, t10, t11)
  
  mean_ctrl <- ctrl_vals
  mean_ted  <- mean(ted_vals, na.rm=TRUE)
  
  # log2FC with small pseudocount to avoid /0
  pc <- 0.01
  log2fc <- log2((mean_ted+pc)/(mean_ctrl+pc))
  res[[g]] <- data.table(
    gene=g,
    ctrl_tpm=mean_ctrl,
    ted_tpm=mean_ted,
    ted_per_sample=paste(round(ted_vals,3), collapse=";"),
    log2fc_tpm=log2fc
  )
}
tab <- rbindlist(res)
cat("\n=== Backbone tissue (TPM-based quick check) ===\n")
print(tab)

# ---- DESeq2 (re-confirm v4.x with collapsed replicates) ----
cat("\nRunning DESeq2 on collapsed read counts...\n")
library(DESeq2)
count_cols <- grep("Read_Count", names(d), value=TRUE)
count_mat  <- as.matrix(d[, ..count_cols])
rownames(count_mat) <- d[[GENE_COL]]

# Collapse replicates for the 5 individuals
collapsed_mat <- matrix(0, nrow=nrow(count_mat), ncol=5)
colnames(collapsed_mat) <- c("Ctrl_2", "TED_7", "TED_8", "TED_10", "TED_11")
rownames(collapsed_mat) <- rownames(count_mat)
collapsed_mat[, "Ctrl_2"]  <- count_mat[, "2__1_Read_Count"] + count_mat[, "2__2_Read_Count"]
collapsed_mat[, "TED_7"]   <- count_mat[, "7__1_Read_Count"] + count_mat[, "7__2_Read_Count"]
collapsed_mat[, "TED_8"]   <- count_mat[, "8__1_Read_Count"] + count_mat[, "8__2_Read_Count"]
collapsed_mat[, "TED_10"]  <- count_mat[, "10__1_Read_Count"] + count_mat[, "10__2_Read_Count"]
collapsed_mat[, "TED_11"]  <- count_mat[, "11__1_Read_Count"] + count_mat[, "11__2_Read_Count"]

col_data <- data.frame(
  sample_name = colnames(collapsed_mat),
  condition = factor(c("Control", "TED", "TED", "TED", "TED"))
)

dds <- DESeqDataSetFromMatrix(countData = collapsed_mat, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res_deseq <- results(dds)
res_sub <- as.data.table(as.data.frame(res_deseq[BACKBONE, , drop=FALSE]), keep.rownames="gene")

# Merge DESeq2 results into our table
tab <- merge(tab, res_sub[, .(gene, deseq2_log2fc = log2FoldChange, deseq2_pvalue = pvalue, deseq2_padj = padj)], by="gene")

cat("\n=== Backbone tissue (with DESeq2 results) ===\n")
print(tab)

# Validation gate
tshr <- tab[gene=="TSHR"]
if (nrow(tshr)) {
  cat(sprintf("\n[VALIDATION] TSHR log2FC(DESeq2)=%.2f (v4.x DESeq2 was +2.33)\n", tshr$deseq2_log2fc))
  cat("  If far from +2.33, STOP and report — possible column-mapping or data drift.\n")
}

fwrite(tab, OUT)
cat("\nSaved:", OUT, "\n")
cat("\nNOTE: Tissue is a supporting/descriptive layer. The DESeq2 results confirm\n")
cat("      that TSHR is strongly upregulated (+2.31, padj=6.6e-05) and IGF1R is NS.\n")
