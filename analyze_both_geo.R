library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)

# Function to run limma on a GSE dataset
run_limma_gse <- function(gse_id, selected_diagnoses, group_mapping) {
  gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
  eset <- gset[[1]]
  pd <- pData(eset)
  
  # Check and filter to Anterior Orbit for GSE58331 (tissue-specific)
  if (gse_id == "GSE58331") {
    pd <- pd[grepl("Anterior Orbit", pd$characteristics_ch1.2, ignore.case = TRUE), ]
  }
  
  # Filter to selected diagnoses
  pd_sub <- pd[pd[["diagnosis:ch1"]] %in% selected_diagnoses, ]
  pd_sub$group <- group_mapping[pd_sub[["diagnosis:ch1"]]]
  
  cat(sprintf("\n--- %s Group Distribution ---\n", gse_id))
  print(table(pd_sub$group))
  
  # Extract expression data
  expr_sub <- exprs(eset)[, rownames(pd_sub)]
  
  # Setup design matrix
  group_factor <- factor(pd_sub$group, levels = c("Control", "TED"))
  design <- model.matrix(~ group_factor)
  colnames(design) <- c("Intercept", "TED_vs_Control")
  
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "TED_vs_Control", number = Inf, adjust.method = "BH")
  results$probe_id <- rownames(results)
  
  # Map probes
  fdata <- fData(eset)
  sym_col <- grep("Gene.?symbol|^Symbol$|GENE_SYMBOL", colnames(fdata), ignore.case = TRUE)[1]
  
  probe_to_gene <- fdata[, c("ID", colnames(fdata)[sym_col]), drop = FALSE]
  colnames(probe_to_gene) <- c("probe_id", "gene")
  probe_to_gene$gene <- sapply(strsplit(as.character(probe_to_gene$gene), "///"), `[`, 1)
  probe_to_gene$gene <- trimws(probe_to_gene$gene)
  probe_to_gene <- probe_to_gene[probe_to_gene$gene != "" & !is.na(probe_to_gene$gene), ]
  
  results_merged <- merge(results, probe_to_gene, by = "probe_id")
  
  # Target genes
  target_genes <- c("TSHR", "IGF1R", "CTLA4", "HAS1")
  
  # Filter and take probe with minimum P-value for each gene
  summary_df <- results_merged %>%
    filter(gene %in% target_genes) %>%
    group_by(gene) %>%
    slice(which.min(P.Value)) %>%
    ungroup() %>%
    select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, probe_id) %>%
    mutate(GSE = gse_id)
  
  return(list(summary = summary_df, full = results_merged))
}

# 1. Analyze GSE58331
gse58331_mapping <- c("Normal" = "Control", "TED" = "TED")
res58331 <- run_limma_gse("GSE58331", c("Normal", "TED"), gse58331_mapping)

# 2. Analyze GSE105149
gse105149_mapping <- c("Normal" = "Control", "TED" = "TED")
res105149 <- run_limma_gse("GSE105149", c("Normal", "TED"), gse105149_mapping)

# Combine and print summary
combined_summary <- bind_rows(res58331$summary, res105149$summary)

cat("\n=== COMBINED SUMMARY FOR TARGET GENES ===\n")
print(combined_summary)

# Write to CSV
write.csv(combined_summary, "GSE_external_replication_summary.csv", row.names = FALSE)
cat("\nResults written to GSE_external_replication_summary.csv\n")
