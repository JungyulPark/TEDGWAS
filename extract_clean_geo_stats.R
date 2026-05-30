library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)

# We will run the analysis for both GSE58331 and GSE105149 using the same clean mapping
# and print log2FC and P-value for the target genes of interest.

run_limma_clean <- function(gse_id, selected_diagnoses, group_mapping, target_genes) {
  # Load local RDS if possible to save time, else get from GEO
  eset_path <- sprintf("C:/ProjectTEDGWAS/TrackA_MR/GEO_validation/%s_eset.rds", gse_id)
  if (file.exists(eset_path)) {
    eset <- readRDS(eset_path)
    pd <- pData(eset)
  } else {
    gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
    eset <- gset[[1]]
    pd <- pData(eset)
  }
  
  if (gse_id == "GSE58331") {
    # filter to Anterior Orbit for tissue specificity
    pd <- pd[grepl("Anterior Orbit", pd$characteristics_ch1.2, ignore.case = TRUE), ]
  }
  
  # Filter to selected diagnoses in characteristics column or diagnosis column
  # In GSE58331, diagnosis column is diagnosis:ch1
  # In GSE105149, diagnosis column is diagnosis:ch1
  diag_col <- grep("diagnosis", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(diag_col)) {
    diag_col <- "characteristics_ch1"
  }
  
  pd_sub <- pd[pd[[diag_col]] %in% selected_diagnoses, ]
  pd_sub$group <- group_mapping[pd_sub[[diag_col]]]
  
  cat(sprintf("\n--- %s Group Distribution ---\n", gse_id))
  print(table(pd_sub$group))
  
  expr_sub <- exprs(eset)[, rownames(pd_sub)]
  
  # check log scale
  if (max(expr_sub, na.rm = TRUE) > 50) {
    expr_sub <- log2(expr_sub + 1)
  }
  
  group_factor <- factor(pd_sub$group, levels = c("Control", "TED"))
  design <- model.matrix(~ group_factor)
  colnames(design) <- c("Intercept", "TED_vs_Control")
  
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "TED_vs_Control", number = Inf, adjust.method = "BH")
  results$probe_id <- rownames(results)
  
  fdata <- fData(eset)
  sym_col <- grep("Gene.?symbol|^Symbol$|GENE_SYMBOL", colnames(fdata), ignore.case = TRUE)[1]
  
  probe_to_gene <- fdata[, c("ID", colnames(fdata)[sym_col]), drop = FALSE]
  colnames(probe_to_gene) <- c("probe_id", "gene")
  probe_to_gene$gene <- sapply(strsplit(as.character(probe_to_gene$gene), "///"), `[`, 1)
  probe_to_gene$gene <- trimws(probe_to_gene$gene)
  probe_to_gene <- probe_to_gene[probe_to_gene$gene != "" & !is.na(probe_to_gene$gene), ]
  
  results_merged <- merge(results, probe_to_gene, by = "probe_id")
  
  summary_df <- results_merged %>%
    filter(toupper(gene) %in% toupper(target_genes)) %>%
    group_by(gene) %>%
    slice(which.min(P.Value)) %>%
    ungroup() %>%
    select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, probe_id) %>%
    mutate(GSE = gse_id)
  
  return(list(summary = summary_df, full = results_merged))
}

target_genes <- c("TSHR", "IGF1R", "INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1", "CTLA4", "HAS1")

res58331 <- run_limma_clean("GSE58331", c("Normal", "TED"), c("Normal" = "Control", "TED" = "TED"), target_genes)
res105149 <- run_limma_clean("GSE105149", c("Normal", "TED"), c("Normal" = "Control", "TED" = "TED"), target_genes)

combined <- bind_rows(res58331$summary, res105149$summary)
print(combined)
