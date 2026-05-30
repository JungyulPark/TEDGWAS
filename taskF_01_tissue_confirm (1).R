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

GENE_COL <- "gene_symbol"      # <<< CONFIRM actual column name
CTRL_TPM <- c("S2_TPM")        # <<< CONFIRM: Sample#2 control TPM col
TED_TPM  <- c("S7_TPM","S8_TPM","S10_TPM","S11_TPM")  # <<< CONFIRM TED TPM cols

# ---- Extract backbone genes ----
sub <- d[get(GENE_COL) %in% BACKBONE]
cat("\nBackbone genes found:", nrow(sub), "of", length(BACKBONE), "\n")

res <- list()
for (g in BACKBONE) {
  row <- sub[get(GENE_COL)==g]
  if (nrow(row)==0) { cat("  MISSING:", g, "\n"); next }
  ctrl_vals <- as.numeric(row[, ..CTRL_TPM])
  ted_vals  <- as.numeric(row[, ..TED_TPM])
  mean_ctrl <- mean(ctrl_vals, na.rm=TRUE)
  mean_ted  <- mean(ted_vals, na.rm=TRUE)
  # log2FC with small pseudocount to avoid /0 (TSHR ctrl ~0.10)
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

# ---- DESeq2 (re-confirm v4.x) ----
# If raw counts available, run DESeq2 for proper padj. If only TPM/FPKM,
# report TPM fold-change and note that v4.x padj came from the count-based
# DESeq2 model. Gemini: use the SAME count input v4.x used to reproduce
# TSHR padj~0.006. If counts not available here, flag and use v4.x stored padj.

# Validation gate
tshr <- tab[gene=="TSHR"]
if (nrow(tshr)) {
  cat(sprintf("\n[VALIDATION] TSHR log2FC(TPM)=%.2f (v4.x DESeq2 was +2.33)\n", tshr$log2fc_tpm))
  cat("  If far from +2.33, STOP and report — possible column-mapping or data drift.\n")
}

fwrite(tab, OUT)
cat("\nSaved:", OUT, "\n")
cat("\nNOTE: attach v4.x DESeq2 padj (TSHR 0.006, IGF1R 0.710) in integration table;\n")
cat("      n=1 control limits per-gene DE inference — tissue is supporting layer.\n")
