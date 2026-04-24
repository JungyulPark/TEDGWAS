#==============================================================================
# TED-TRAP Upgrade — Phase 2 Script 12
# Gene-set testing with CAMERA/ROAST (candidate set enrichment)
#==============================================================================
# Purpose:
#   Resolve MAJOR gene-set testing rigor. For the RNA-seq data, we use
#   CAMERA (Wu & Smyth 2012) and/or ROAST (Wu et al. 2010), which properly
#   account for inter-gene correlation within a gene set — an important
#   correction that per-gene BH alone does not provide.
#
#   These methods test whether the MR-triangulated candidate gene set shows
#   differential expression as a group, which is the scientifically
#   correct question for a candidate-gene confirmation design.
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(limma)
library(edgeR)
library(data.table)

log_file <- "TrackA_MR/logs/12_gene_set_testing.log"
sink(log_file, split = TRUE)
cat("=== Gene-Set Testing: CAMERA/ROAST ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- Load count matrix ---
count_file <- "TrackA_MR/data/rnaseq/counts_gene.tsv"
meta_file  <- "TrackA_MR/data/rnaseq/sample_metadata.tsv"

if (!file.exists(count_file) || !file.exists(meta_file)) {
  cat("❌ Count matrix or metadata missing. See script 10 for expected format.\n")
  stop()
}

counts <- fread(count_file, data.table = FALSE)
rownames(counts) <- counts[[1]]
counts <- as.matrix(counts[, -1])

meta <- fread(meta_file, data.table = FALSE)
rownames(meta) <- meta$sample_id
meta <- meta[colnames(counts), ]

cat(sprintf("Loaded %d genes × %d samples\n", nrow(counts), ncol(counts)))

# --- edgeR-style normalization and voom transformation ---
dge <- DGEList(counts = counts, group = meta$condition)
keep <- filterByExpr(dge, group = meta$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
cat(sprintf("After filterByExpr: %d genes\n", nrow(dge)))

design <- model.matrix(~condition, data = meta)
v <- voom(dge, design, plot = FALSE)

# --- Define gene sets ---
# Primary: MR-triangulated candidate set
gene_sets <- list(
  "MR_triangulated" = c(
    "IGF1R", "IGF1", "IL1B", "CXCL8",
    "TSHR", "ADIPOQ", "FABP4",
    "HAS1", "HAS2", "HAS3", "TNF", "IL6",
    "ARRB1", "PPARG", "CEBPA"
  ),
  "IGF1R_pathway" = c(
    "IGF1R","IGF1","IGF2","IRS1","IRS2","SHC1","GRB2","SOS1",
    "PIK3CA","PIK3CB","PIK3R1","AKT1","AKT2","MTOR","RPS6KB1",
    "KRAS","BRAF","MAP2K1","MAP2K2","MAPK1","MAPK3",
    "JAK1","JAK2","STAT3","STAT5A",
    "ARRB1","ARRB2",
    "HAS1","HAS2","HAS3","PPARG","CEBPA",
    "IL6","CXCL8","TNF","IL1B"
  ),
  "TSHR_pathway" = c(
    "TSHR","GNAS","GNAQ","GNA11",
    "ADCY1","ADCY2","ADCY3","ADCY4","ADCY5",
    "ADCY6","ADCY7","ADCY8","ADCY9",
    "PRKACA","PRKACB","PRKAR1A","CREB1","CREB3",
    "ARRB1","ARRB2","GRK2","GRK5",
    "FOXO1","FOXO3",
    "HAS1","HAS2","HAS3","IL6","TNF",
    "PPARG","CEBPA","FABP4","ADIPOQ"
  ),
  "Insulin_signaling_KEGG_hsa04910" = c(
    "INSR", "IRS1", "IRS2", "PIK3CA", "PIK3R1", "AKT1", "AKT2",
    "PDPK1", "FOXO1", "FOXO3", "GSK3B", "MTOR", "RPS6KB1",
    "PYGL", "PCK1", "G6PC", "MAPK1", "MAPK3"
  )
)

# Convert gene symbols to row indices
gene_set_indices <- lapply(gene_sets, function(gs) {
  which(rownames(v) %in% gs)
})

cat("\nGene set sizes (genes found in data):\n")
for (nm in names(gene_set_indices)) {
  cat(sprintf("  %s: %d / %d\n", nm, length(gene_set_indices[[nm]]),
              length(gene_sets[[nm]])))
}

# --- CAMERA gene-set test ---
cat("\n--- CAMERA (Wu & Smyth 2012) ---\n")
camera_res <- camera(y = v, index = gene_set_indices, design = design,
                     contrast = "conditionTED")
print(camera_res)

# --- ROAST gene-set test ---
cat("\n--- ROAST (Wu et al. 2010) ---\n")
roast_res <- mroast(y = v, index = gene_set_indices, design = design,
                    contrast = "conditionTED", nrot = 9999)
print(roast_res)

# --- Save ---
camera_df <- data.frame(gene_set = rownames(camera_res), camera_res)
roast_df  <- data.frame(gene_set = rownames(roast_res),  roast_res)
fwrite(camera_df, "TrackA_MR/results/12_camera_results.csv")
fwrite(roast_df,  "TrackA_MR/results/12_roast_results.csv")

cat("\n✅ CAMERA/ROAST complete.\n")
cat("   Interpretation: Look for Direction=Up for 'Insulin_signaling_KEGG_hsa04910'\n")
cat("   in TED vs Control — this directly confirms/refutes the insulin-pathway\n")
cat("   off-target hypothesis at the gene-set level.\n")

sink()
