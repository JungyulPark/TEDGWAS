# === Setup ===
library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)
library(stringr)

setwd("c:/ProjectTEDGWAS/TrackA_MR/external_validation/")

# === Step 1: Download GSE58331 ===
gse <- getGEO("GSE58331", GSEMatrix = TRUE, AnnotGPL = TRUE)
eset <- gse[[1]]

cat("Total samples in GSE58331:", ncol(eset), "\n")
cat("Total probes:", nrow(eset), "\n")

# === Step 2: Extract metadata + filter to orbital adipose only ===
pdata <- pData(eset)

# Filter to orbital adipose tissue only (tissue: Anterior Orbit)
adipose_samples <- pdata[grepl("Anterior Orbit", pdata$characteristics_ch1.2, ignore.case = TRUE), ]
cat("Anterior Orbit (adipose) samples:", nrow(adipose_samples), "\n")

# Filter to TED + control only (exclude sarcoidosis, GPA, NSOI)
adipose_TED_or_control <- adipose_samples[
  grepl("TED|Normal", adipose_samples$characteristics_ch1, ignore.case = TRUE), 
]

# Group assignment
adipose_TED_or_control$group <- ifelse(
  grepl("TED", adipose_TED_or_control$characteristics_ch1, ignore.case = TRUE),
  "TED", "Control"
)

cat("\nSample distribution by group:\n")
print(table(adipose_TED_or_control$group))
# Expected: TED 23, Control 20

# === Step 3: Subset expression matrix ===
selected_GSMs <- rownames(adipose_TED_or_control)
expr <- exprs(eset)[, selected_GSMs]

# === Step 4: Map probes to gene symbols ===
fdata <- fData(eset)
gene_col <- "Gene symbol"

# Annotate probes
probe_to_gene <- fdata[, c("ID", gene_col), drop = FALSE]
colnames(probe_to_gene) <- c("probe_id", "gene")

# Some symbols might be separated by "///"
# We will just split them and take the first one or keep as is. Usually AnnotGPL has clean symbols.
probe_to_gene$gene <- sapply(strsplit(as.character(probe_to_gene$gene), "///"), `[`, 1)
probe_to_gene$gene <- trimws(probe_to_gene$gene)

probe_to_gene <- probe_to_gene[probe_to_gene$gene != "" & !is.na(probe_to_gene$gene), ]

# === Step 5: 20-gene panel filter ===
target_genes <- c(
  "TSHR", "IGF1R", "IGF1",
  "ARRB1", "PPARG", "IRS1", "AKT1", "TNF", "CTLA4",
  "INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1",
  "HAS1", "HAS2", "HAS3", "ADIPOQ", "FABP4", "CEBPA"
)

# Get probes for our 20 genes
panel_probes <- probe_to_gene[probe_to_gene$gene %in% target_genes, ]
cat("\nProbes mapping to our 20 genes:", nrow(panel_probes), "\n")
cat("Genes detected:", length(unique(panel_probes$gene)), "/ 20\n")

# === Step 6: limma DE analysis (전체 expression matrix에서) ===
group_factor <- factor(adipose_TED_or_control$group, levels = c("Control", "TED"))
design <- model.matrix(~ group_factor)
colnames(design) <- c("Intercept", "TED_vs_Control")

fit <- lmFit(expr, design)
fit <- eBayes(fit)
results_full <- topTable(fit, coef = "TED_vs_Control", number = Inf, adjust.method = "BH")
results_full$probe_id <- rownames(results_full)

# Merge with gene symbols
results_full <- merge(results_full, probe_to_gene, by = "probe_id")

# === Step 7: Filter to 20-gene panel + summarize per gene ===
panel_results <- results_full[results_full$gene %in% target_genes, ]

# Take probe with smallest P
panel_summary <- panel_results %>%
  group_by(gene) %>%
  slice(which.min(P.Value)) %>%
  ungroup() %>%
  select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
  arrange(P.Value)

# === Step 8: Save outputs ===
write.csv(panel_summary, "GSE58331_orbital_adipose_20gene_panel_results.csv", row.names = FALSE)
write.csv(results_full, "GSE58331_orbital_adipose_full_DE_results.csv", row.names = FALSE)

# === Summary report ===
cat("\n=== GSE58331 Orbital Adipose Cross-Validation Summary ===\n")
cat("Total TED samples:", sum(group_factor == "TED"), "\n")
cat("Total Control samples:", sum(group_factor == "Control"), "\n")
cat("\n20-gene panel results:\n")
print(panel_summary, n = 20)

cat("\n=== Direction concordance check ===\n")
in_house_directions <- c(
  TSHR = "UP", INSR = "UP", IRS2 = "UP", FOXO1 = "UP", 
  PIK3R1 = "UP", PDPK1 = "UP"
)

for (g in names(in_house_directions)) {
  gene_row <- panel_summary[panel_summary$gene == g, ]
  if (nrow(gene_row) == 0) {
    cat(sprintf("%-10s: NOT DETECTED in GSE58331\n", g))
  } else {
    gse_dir <- ifelse(gene_row$logFC > 0, "UP", "DOWN")
    concordant <- ifelse(in_house_directions[g] == gse_dir, "✓ concordant", "✗ discordant")
    cat(sprintf("%-10s: in-house=%s, GSE58331=%s (logFC=%+.3f, P=%.4f) %s\n",
                g, in_house_directions[g], gse_dir, gene_row$logFC, gene_row$P.Value, concordant))
  }
}
