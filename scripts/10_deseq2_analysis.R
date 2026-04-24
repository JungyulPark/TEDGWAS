#==============================================================================
# TED-TRAP Upgrade — Phase 2 Script 10
# RNA-seq Differential Expression with DESeq2 (replaces Wilcoxon)
#==============================================================================
# Purpose:
#   - Resolve CRITICAL issue #2: Wilcoxon rank-sum test is mathematically
#     incoherent for n=4 vs n=1 (exact null has only 5 permutations;
#     minimum two-sided p ≈ 0.4).
#   - Use DESeq2 with empirical Bayes dispersion shrinkage across all genes,
#     which tolerates n=1 per group through across-gene information sharing.
#   - Apply BH-FDR correction within a pre-specified candidate gene set
#     (MR-triangulated genes).
#   - Frame RNA-seq as orthogonal tissue-localization, not discovery DE.
#
# Required input:
#   - Gene-level count matrix (rows = genes, cols = samples)
#   - Sample metadata (sample ID, condition [TED/Control])
#
# Output:
#   TrackA_MR/results/10_deseq2_v3_results.csv
#   TrackA_MR/results/10_deseq2_v3_candidate_set.csv  (BH-FDR applied)
#   TrackA_MR/figures/10_ma_plot.pdf
#   TrackA_MR/figures/10_pca.pdf
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(DESeq2)
library(apeglm)
library(data.table)
library(ggplot2)
library(dplyr)

log_file <- "TrackA_MR/logs/10_deseq2.log"
sink(log_file, split = TRUE)
cat("=== DESeq2 Differential Expression Analysis ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- File paths ---
count_file <- "TrackA_MR/data/rnaseq/counts_gene.tsv"
meta_file  <- "TrackA_MR/data/rnaseq/sample_metadata.tsv"

if (!file.exists(count_file)) {
  cat("❌ Count matrix not found at:\n")
  cat(sprintf("   %s\n\n", count_file))
  cat("EXPECTED FORMAT:\n")
  cat("  Tab-separated file with:\n")
  cat("  - First column: gene_id (Ensembl ID or HGNC symbol)\n")
  cat("  - Remaining columns: raw count per sample (integer counts, NOT TPM)\n\n")
  cat("If only TPM is available, DESeq2 requires counts. Options:\n")
  cat("  1. Re-run upstream quantification (STAR + featureCounts, or salmon + tximport)\n")
  cat("  2. Use length-scaled TPM via tximport\n")
  cat("  3. Apply limma-voom on log2(TPM+1) as a reasonable approximation (disclosed)\n\n")
  stop("Provide count matrix to proceed.")
}

if (!file.exists(meta_file)) {
  cat("❌ Metadata not found at:\n")
  cat(sprintf("   %s\n\n", meta_file))
  cat("EXPECTED FORMAT:\n")
  cat("  sample_id<TAB>condition\n")
  cat("  Ctrl_01<TAB>Control\n")
  cat("  TED_01<TAB>TED\n")
  cat("  ... etc.\n")
  stop("Provide metadata to proceed.")
}

# --- Load data ---
counts <- fread(count_file, data.table = FALSE)
meta   <- fread(meta_file, data.table = FALSE)

# Prepare count matrix (rows = genes, cols = samples)
rownames(counts) <- counts[[1]]
counts <- as.matrix(counts[, -1])
storage.mode(counts) <- "integer"

cat(sprintf("Loaded: %d genes × %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("Samples: %s\n", paste(colnames(counts), collapse=", ")))

# Align metadata
rownames(meta) <- meta$sample_id
meta <- meta[colnames(counts), ]
meta$condition <- factor(meta$condition, levels = c("Control", "TED"))
cat(sprintf("Condition distribution: Control=%d, TED=%d\n\n",
            sum(meta$condition == "Control"),
            sum(meta$condition == "TED")))

# --- Build DESeq2 dataset ---
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

# --- Pre-filter low-count genes (improves DESeq2 stability) ---
keep <- rowSums(counts(dds) >= 10) >= 2
cat(sprintf("Pre-filter: keeping %d / %d genes with ≥10 counts in ≥2 samples\n",
            sum(keep), nrow(dds)))
dds <- dds[keep, ]

# --- Run DESeq2 ---
cat("Running DESeq2 (this may take 1-3 min for small datasets)...\n")
dds <- DESeq(dds)

# --- Extract results with apeglm shrinkage (recommended 2024-2025) ---
# apeglm shrinks effect sizes via empirical Bayes, improving small-n estimation
res <- lfcShrink(dds, coef = "condition_TED_vs_Control", type = "apeglm")

# Convert to data.frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

# Sort by p-value
res_df <- res_df[order(res_df$pvalue), ]

# Save full results
fwrite(res_df, "TrackA_MR/results/10_deseq2_v3_results.csv")
cat(sprintf("✅ Full DESeq2 results: %d genes\n", nrow(res_df)))

# --- Candidate gene set testing (MR-triangulated) ---
cat("\n--- Candidate Gene Set Analysis (pre-specified) ---\n")
candidate_genes <- c(
  # Three-way intersection (15 genes from Phase 1 pathway)
  "IGF1R", "IGF1", "IL1B", "CXCL8",              # IGF1R-only (4)
  "TSHR", "ADIPOQ", "FABP4",                     # TSHR-only (3)
  "HAS1", "HAS2", "HAS3", "TNF", "IL6",
  "ARRB1", "PPARG", "CEBPA",                     # Shared (8)
  # Off-target candidates
  "INSR", "IRS1", "IRS2", "AKT2", "PIK3CD"
)

# Subset
cand_res <- res_df[res_df$gene %in% candidate_genes, ]
cat(sprintf("Candidate genes found in DE results: %d / %d\n",
            nrow(cand_res), length(candidate_genes)))

# Apply BH-FDR to the restricted candidate set (independent of genome-wide FDR)
cand_res$padj_candidate_BH <- p.adjust(cand_res$pvalue, method = "BH")

# Add zone annotation
zone_map <- data.frame(
  gene = candidate_genes,
  zone = c(rep("IGF1R-only", 4), rep("TSHR-only", 3),
           rep("Shared", 8), rep("Off-target", 5)),
  stringsAsFactors = FALSE
)
cand_res <- merge(cand_res, zone_map, by = "gene", all.x = TRUE)
cand_res <- cand_res[order(cand_res$pvalue), ]

fwrite(cand_res, "TrackA_MR/results/10_deseq2_v3_candidate_set.csv")
cat("✅ Candidate-set results (BH corrected):\n")
print(cand_res[, c("gene", "zone", "log2FoldChange", "pvalue", "padj_candidate_BH")],
      row.names = FALSE)

# --- QC plots ---
cat("\n--- Generating QC plots ---\n")

# MA plot
pdf("TrackA_MR/figures/10_ma_plot.pdf", width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), alpha = 0.05, colSig = "red",
       main = "DESeq2 MA plot (shrunken LFC)")
abline(h = 0, col = "blue", lty = 2)
dev.off()

# PCA
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel() +
  xlab(paste0("PC1 (", percent_var[1], "%)")) +
  ylab(paste0("PC2 (", percent_var[2], "%)")) +
  theme_bw(base_size = 12) +
  ggtitle("Sample-level PCA (VST-transformed)")
ggsave("TrackA_MR/figures/10_pca.pdf", p_pca, width = 7, height = 5)

# Volcano plot of candidate genes
cand_res$sig <- ifelse(cand_res$padj_candidate_BH < 0.05, "BH q<0.05", "n.s.")
p_volcano <- ggplot(cand_res, aes(x = log2FoldChange, y = -log10(pvalue),
                                   color = sig, label = gene)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3) +
  scale_color_manual(values = c("BH q<0.05" = "red", "n.s." = "grey50")) +
  theme_bw(base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  xlab("log2(TED / Control)") + ylab("-log10(p)") +
  ggtitle("Candidate genes: DESeq2 with apeglm shrinkage (BH within candidate set)")
ggsave("TrackA_MR/figures/10_volcano_candidate.pdf", p_volcano, width = 8, height = 6)

cat("✅ QC plots saved in TrackA_MR/figures/\n")

# --- Summary ---
cat("\n\n=== METHOD SUMMARY (for Methods section) ===\n")
cat("Differential expression analysis was performed on gene-level integer counts\n")
cat("using DESeq2 (Love et al. 2014, Genome Biol 15:550), which applies empirical\n")
cat("Bayes dispersion shrinkage across all genes — allowing reliable inference even\n")
cat("with unbalanced small sample sizes. Effect sizes were moderated via apeglm\n")
cat("(Zhu et al. 2019, Bioinformatics 35:2084). We pre-specified a candidate gene\n")
cat("set comprising the 15 three-way-intersection genes (IGF-1R pathway ∩ TSHR\n")
cat("pathway ∩ TED disease gene set) plus 5 insulin-signaling off-target genes\n")
cat("nominated by pathway analysis; multiple testing within this candidate set\n")
cat("was controlled using Benjamini-Hochberg FDR at q < 0.05. We explicitly frame\n")
cat("this RNA-seq analysis as orthogonal tissue-level confirmation of MR-nominated\n")
cat("signals, not as genome-wide discovery; the limited sample size (4 inactive TED\n")
cat("vs 1 control) precludes powered genome-wide DE (Schurch et al. 2016, RNA\n")
cat("22:839). Direction-of-effect concordance, not p-values alone, anchors\n")
cat("interpretation for genes failing the BH threshold.\n")

sink()
cat("\n✅ Phase 2 Script 10 complete.\n\n")
