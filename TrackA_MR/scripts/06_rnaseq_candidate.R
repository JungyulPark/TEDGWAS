# =============================================================================
# TED-TRAP Step 6 — Pre-Specified Candidate Gene Tissue Validation
#
# Design: Hypothesis-driven DESeq2 on 20 pre-specified genes from MR/coloc
# evidence. NOT a genome-wide DE — orthogonal to any prior descriptive analysis.
#
# Primary aim: SE3 — Orbital insulin signaling cassette (INSR, IRS2, FOXO1,
# PIK3R1, PDPK1) directionality in TED vs Control.
#
# Statistical method: DESeq2 + apeglm LFC shrinkage + CAMERA gene-set test
# Multiple testing: BH-FDR within the pre-specified 20-gene set only
# =============================================================================

# --- Packages ---
required_cran <- c("data.table", "ggplot2", "dplyr")
required_bioc <- c("DESeq2", "apeglm", "limma")

for (p in required_cran) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (p in required_bioc) {
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(apeglm)
    library(limma)
    library(ggplot2)
    library(dplyr)
})

for (d in c("TrackA_MR/results", "TrackA_MR/figures", "TrackA_MR/logs")) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
sink("TrackA_MR/logs/06_rnaseq_candidate.log", split = TRUE)
cat("=== Step 6: Pre-specified Candidate Gene Tissue Validation ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# =============================================================================
# 1. Load and prepare count matrix
# =============================================================================

DATA_PATH <- "c:/ProjectTEDGWAS/data.txt"

cat("--- Loading RNA-seq data ---\n")
dat <- fread(DATA_PATH)
cat(sprintf("Loaded: %d rows × %d cols\n", nrow(dat), ncol(dat)))
cat("First 5 columns:\n")
print(colnames(dat)[1:min(10, ncol(dat))])

# Identify count columns (raw counts, not TPM/FPKM)
count_cols <- grep("_count$|_Count$|_COUNT$|__[0-9]+_count", colnames(dat), value = TRUE)
if (length(count_cols) == 0) {
    # Try alternative patterns
    all_cols <- colnames(dat)
    # Look for numeric sample columns that are NOT TPM/FPKM
    numeric_cols <- all_cols[sapply(dat, is.numeric)]
    tpm_cols <- grep("TPM|FPKM|tpm|fpkm", numeric_cols, value = TRUE)
    count_cols <- setdiff(numeric_cols, tpm_cols)
    # Keep only sample-like columns (typical pattern: __1_, __2_ suffix)
    count_cols <- grep("__[12]", count_cols, value = TRUE)
    count_cols <- count_cols[!grepl("TPM|FPKM|tpm|fpkm", count_cols)]
}

cat(sprintf("\nDetected count columns (%d):\n", length(count_cols)))
print(count_cols)

if (length(count_cols) != 10) {
    cat("\n[WARN] Expected 10 count columns (5 samples × 2 tech reps).\n")
    cat("Showing all column names for manual review:\n")
    print(colnames(dat))
    sink()
    stop("Column pattern mismatch — adjust count_cols manually")
}

# Gene ID column
gene_col <- if ("Gene_Symbol" %in% colnames(dat)) {
    "Gene_Symbol"
} else if ("gene_symbol" %in% colnames(dat)) {
    "gene_symbol"
} else if ("gene" %in% colnames(dat)) {
    "gene"
} else {
    colnames(dat)[1]
}
cat(sprintf("\nGene column: %s\n", gene_col))

# Extract counts as matrix
counts_mat <- as.matrix(dat[, ..count_cols])
rownames(counts_mat) <- dat[[gene_col]]

# Integer conversion (DESeq2 requires integer counts)
counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"

# Remove genes with all-zero counts
keep_nonzero <- rowSums(counts_mat) > 0
counts_mat <- counts_mat[keep_nonzero, ]
cat(sprintf(
    "\nCount matrix: %d genes × %d libraries\n",
    nrow(counts_mat), ncol(counts_mat)
))

# =============================================================================
# 2. Collapse technical replicates → 5 samples
# =============================================================================

cat("\n--- Collapsing technical replicates ---\n")

# Sample ID pattern: "2__1_..." "2__2_..." → patient #2
extract_patient_id <- function(col) {
    sub("__[12].*$", "", col)
}
patient_ids <- sapply(count_cols, extract_patient_id)
unique_patients <- unique(patient_ids)
cat(sprintf("Unique patients: %d\n", length(unique_patients)))
cat("Patient IDs:\n")
print(unique_patients)

# Sum technical replicates (standard DESeq2 approach)
counts_collapsed <- sapply(unique_patients, function(pid) {
    cols <- count_cols[patient_ids == pid]
    rowSums(counts_mat[, cols, drop = FALSE])
})
colnames(counts_collapsed) <- unique_patients
storage.mode(counts_collapsed) <- "integer"
cat(sprintf(
    "After collapsing: %d genes × %d samples\n",
    nrow(counts_collapsed), ncol(counts_collapsed)
))

# =============================================================================
# 3. Sample metadata — Control identified via TSHR TPM (lowest = C)
# =============================================================================

cat("\n--- Sample assignment (from TSHR TPM) ---\n")
# Based on prior analysis: Sample #2 is Control (TSHR TPM 0.10), rest are TED

control_id <- "2" # From TSHR TPM analysis
colData <- data.frame(
    patient = colnames(counts_collapsed),
    condition = ifelse(colnames(counts_collapsed) == control_id, "Control", "TED"),
    stringsAsFactors = FALSE
)
colData$condition <- factor(colData$condition, levels = c("Control", "TED"))
rownames(colData) <- colData$patient

print(colData)
cat(sprintf("\nControl samples: %d\n", sum(colData$condition == "Control")))
cat(sprintf("TED samples: %d\n", sum(colData$condition == "TED")))

# =============================================================================
# 4. Pre-specified 20-gene candidate set
# =============================================================================

candidate_genes <- list(
    # MR Primary targets
    MR_Primary = c("TSHR", "IGF1R", "IGF1"),
    # MR Secondary candidates
    MR_Secondary = c("ARRB1", "PPARG", "IRS1", "AKT1", "TNF", "CTLA4"),
    # Insulin signaling cassette (SE3 core)
    Insulin_Cassette = c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1"),
    # TED-specific biology
    TED_Biology = c("HAS1", "HAS2", "HAS3", "ADIPOQ", "FABP4", "CEBPA")
)

all_candidates <- unique(unlist(candidate_genes))
cat(sprintf("\nPre-specified candidate set: %d genes in 4 groups\n", length(all_candidates)))

# Which candidates are present in data?
found_genes <- all_candidates[all_candidates %in% rownames(counts_collapsed)]
missing_genes <- setdiff(all_candidates, found_genes)
cat(sprintf("Present in data: %d (%s)\n", length(found_genes), paste(found_genes, collapse = ", ")))
if (length(missing_genes) > 0) {
    cat(sprintf("Missing: %s\n", paste(missing_genes, collapse = ", ")))
}

# =============================================================================
# 5. DESeq2 analysis — FULL matrix (for normalization stability)
# =============================================================================

cat("\n--- DESeq2 normalization & testing ---\n")

dds <- DESeqDataSetFromMatrix(
    countData = counts_collapsed,
    colData = colData,
    design = ~condition
)
# Pre-filter low count genes (standard)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("After low-count filter: %d genes\n", nrow(dds)))

# Run DESeq2
dds <- DESeq(dds, quiet = TRUE)
cat("DESeq2 model fit complete.\n")

# Raw results (for reference)
res_raw <- results(dds, contrast = c("condition", "TED", "Control"))

# apeglm LFC shrinkage (Zhu 2018) — recommended for small-N
res_shrink <- tryCatch(
    lfcShrink(dds, coef = "condition_TED_vs_Control", type = "apeglm"),
    error = function(e) {
        cat(sprintf("  [WARN] apeglm failed: %s\n  Falling back to ashr\n", e$message))
        lfcShrink(dds, coef = "condition_TED_vs_Control", type = "ashr")
    }
)
cat("LFC shrinkage complete.\n")

# =============================================================================
# 6. Extract candidate-only results + apply BH-FDR within candidate set
# =============================================================================

cat("\n--- Candidate-specific BH-FDR adjustment ---\n")

res_candidates <- as.data.frame(res_shrink[found_genes, ])
res_candidates$gene <- rownames(res_candidates)

# Assign group
res_candidates$group <- NA_character_
for (grp_name in names(candidate_genes)) {
    res_candidates$group[res_candidates$gene %in% candidate_genes[[grp_name]]] <- grp_name
}

# Candidate-only BH-FDR (critical: NOT genome-wide padj)
res_candidates$padj_candidate <- p.adjust(res_candidates$pvalue, method = "BH")

# Reorder
res_candidates <- res_candidates[, c(
    "gene", "group", "baseMean",
    "log2FoldChange", "lfcSE",
    "pvalue", "padj_candidate"
)]
res_candidates <- res_candidates[order(res_candidates$pvalue), ]

cat("\n--- Candidate-set DESeq2 results (BH-FDR within 20 genes) ---\n")
print_df <- res_candidates
print_df$log2FoldChange <- round(print_df$log2FoldChange, 3)
print_df$lfcSE <- round(print_df$lfcSE, 3)
print_df$pvalue <- signif(print_df$pvalue, 3)
print_df$padj_candidate <- signif(print_df$padj_candidate, 3)
print_df$baseMean <- round(print_df$baseMean, 1)
print(print_df, row.names = FALSE)

write.csv(res_candidates, "TrackA_MR/results/06_candidate_deseq2.csv", row.names = FALSE)

# =============================================================================
# 7. CAMERA gene-set test for Insulin Cassette (SE3 core)
# =============================================================================

cat("\n--- CAMERA gene-set test: Insulin Cassette ---\n")

# Use variance-stabilized counts for limma framework
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Build design matrix (for limma-style gene-set test)
design <- model.matrix(~condition, data = colData)

# Insulin cassette gene indices
insulin_cassette <- candidate_genes$Insulin_Cassette
insulin_cassette_in <- intersect(insulin_cassette, rownames(vsd_mat))
cat(sprintf(
    "Insulin cassette genes tested: %d (%s)\n",
    length(insulin_cassette_in), paste(insulin_cassette_in, collapse = ", ")
))

if (length(insulin_cassette_in) >= 3) {
    # CAMERA requires a list of index vectors
    cassette_idx <- list(Insulin_Cassette = match(insulin_cassette_in, rownames(vsd_mat)))

    camera_res <- camera(vsd_mat, cassette_idx,
        design = design,
        contrast = "conditionTED",
        inter.gene.cor = 0.01
    )
    cat("\nCAMERA result:\n")
    print(camera_res)

    write.csv(camera_res, "TrackA_MR/results/06_camera_insulin_cassette.csv")
} else {
    cat("  [SKIP] Insufficient cassette genes for CAMERA\n")
    camera_res <- NULL
}

# =============================================================================
# 8. SE3 Verdict — Direction concordance & statistical significance
# =============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("SE3 VERDICT — Orbital Insulin Cassette (INSR + co-effectors)\n")
cat(rep("=", 70), "\n\n", sep = "")

# Criterion 1: ≥ 4/5 insulin cassette genes log2FC > 0
cassette_results <- res_candidates[res_candidates$group == "Insulin_Cassette", ]
cat("Insulin cassette gene directionality (log2FC):\n")
for (i in seq_len(nrow(cassette_results))) {
    g <- cassette_results$gene[i]
    lfc <- cassette_results$log2FoldChange[i]
    p <- cassette_results$pvalue[i]
    padj <- cassette_results$padj_candidate[i]
    mark <- ifelse(lfc > 0, "⬆️", "⬇️")
    sig <- ifelse(padj < 0.05, " ***", ifelse(p < 0.05, " *", ""))
    cat(sprintf(
        "  %-7s %s log2FC=%+.3f (p=%.3g, padj=%.3g)%s\n",
        g, mark, lfc, p, padj, sig
    ))
}

n_up <- sum(cassette_results$log2FoldChange > 0)
n_total <- nrow(cassette_results)
cat(sprintf("\nUpregulated: %d / %d\n", n_up, n_total))

# Criterion 2: CAMERA P < 0.05
cat(sprintf(
    "CAMERA gene-set test P: %s\n",
    if (!is.null(camera_res)) signif(camera_res$PValue, 3) else "N/A"
))

cat("\n--- SE3 Interpretation ---\n")
if (n_up >= 4 && !is.null(camera_res) && camera_res$PValue < 0.05) {
    cat("✅ SE3 ACHIEVED: Insulin signaling cassette coordinately upregulated in TED orbital tissue.\n")
    cat("   Narrative: Orbital tissue-level insulin signaling activation provides\n")
    cat("   mechanistic basis for teprotumumab hyperglycemia independent of\n")
    cat("   germline genetic susceptibility.\n")
} else if (n_up >= 4) {
    cat("🟡 SE3 PARTIAL: Directional concordance confirmed, CAMERA test marginal.\n")
} else if (n_up == n_total && n_total >= 3) {
    cat("🟡 SE3 SUGGESTIVE: All cassette members trend upward.\n")
} else {
    cat("⚪ SE3 NOT CLEAR: Mixed or weak directional signal.\n")
}

# Highlight INSR specifically
insr_row <- res_candidates[res_candidates$gene == "INSR", ]
if (nrow(insr_row) > 0) {
    cat(sprintf(
        "\n🔑 INSR highlight: log2FC = %+.3f (p = %.3g, padj = %.3g)\n",
        insr_row$log2FoldChange, insr_row$pvalue, insr_row$padj_candidate
    ))
}

# =============================================================================
# 9. MR–RNA-seq direction concordance table (for Triangulation)
# =============================================================================

cat("\n\n--- Direction Concordance with MR findings ---\n")

# Expected direction from Step 2 (GD outcome, log-odds scale)
mr_dirs <- data.frame(
    gene = c("TSHR", "IGF1R", "ARRB1", "PPARG", "IRS1", "AKT1", "TNF", "CTLA4"),
    MR_beta = c(-1.408, 0.200, -0.177, -0.097, -0.138, 0.081, NA, -0.943),
    stringsAsFactors = FALSE
)

tri_df <- merge(res_candidates, mr_dirs, by = "gene", all.x = TRUE)
tri_df$RNAseq_dir <- ifelse(tri_df$log2FoldChange > 0, "+", "-")
tri_df$MR_dir <- ifelse(is.na(tri_df$MR_beta), "NA",
    ifelse(tri_df$MR_beta > 0, "+", "-")
)
tri_df$concordance <- ifelse(tri_df$MR_dir == "NA", "—",
    ifelse(tri_df$MR_dir == tri_df$RNAseq_dir, "✅", "✖")
)

tri_subset <- tri_df[!is.na(tri_df$MR_beta), c(
    "gene", "group", "MR_beta", "MR_dir",
    "log2FoldChange", "RNAseq_dir", "concordance"
)]
tri_subset$log2FoldChange <- round(tri_subset$log2FoldChange, 3)
tri_subset$MR_beta <- round(tri_subset$MR_beta, 3)
cat("\nTriangulation table (MR × RNA-seq direction):\n")
print(tri_subset, row.names = FALSE)

write.csv(tri_df, "TrackA_MR/results/06_triangulation_MR_RNAseq.csv", row.names = FALSE)

# =============================================================================
# 10. Candidate-only forest plot (novel figure, different from prior paper)
# =============================================================================

cat("\n--- Generating forest plot ---\n")

plot_df <- res_candidates[!is.na(res_candidates$log2FoldChange), ]
plot_df$gene <- factor(plot_df$gene, levels = rev(plot_df$gene[order(plot_df$group, plot_df$log2FoldChange)]))
plot_df$group <- factor(plot_df$group,
    levels = c(
        "MR_Primary", "MR_Secondary",
        "Insulin_Cassette", "TED_Biology"
    )
)
plot_df$significance <- ifelse(plot_df$padj_candidate < 0.05, "padj<0.05",
    ifelse(plot_df$pvalue < 0.05, "p<0.05 (nominal)",
        "n.s."
    )
)

p_forest <- ggplot(plot_df, aes(x = log2FoldChange, y = gene, color = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(
        aes(
            xmin = log2FoldChange - 1.96 * lfcSE,
            xmax = log2FoldChange + 1.96 * lfcSE
        ),
        height = 0.3
    ) +
    geom_point(aes(size = -log10(pvalue), shape = significance)) +
    scale_shape_manual(values = c("padj<0.05" = 16, "p<0.05 (nominal)" = 17, "n.s." = 1)) +
    scale_color_brewer(palette = "Set1") +
    labs(
        title = "Pre-specified candidate genes: TED vs Control",
        subtitle = "DESeq2 with apeglm LFC shrinkage | BH-FDR within 20-gene set",
        x = "log2 Fold Change (TED / Control)",
        y = NULL, color = "Group", shape = "Significance", size = "-log10(p)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "right",
        panel.grid.minor = element_blank()
    )

ggsave("TrackA_MR/figures/06_candidate_forest.pdf", p_forest,
    width = 9, height = 7
)
ggsave("TrackA_MR/figures/06_candidate_forest.png", p_forest,
    width = 9, height = 7, dpi = 300
)
cat("  Saved: 06_candidate_forest.pdf / .png\n")

# =============================================================================
# SAVE SUMMARY
# =============================================================================

sink()
cat("\n✅ Step 6 complete.\n")
cat("   Results: TrackA_MR/results/06_candidate_deseq2.csv\n")
cat("            TrackA_MR/results/06_camera_insulin_cassette.csv\n")
cat("            TrackA_MR/results/06_triangulation_MR_RNAseq.csv\n")
cat("   Figure:  TrackA_MR/figures/06_candidate_forest.pdf\n\n")
