#!/usr/bin/env Rscript
# =============================================================================
# Task D.1 — Druggable Gene Master List
# =============================================================================
# Purpose : Build the druggable gene universe (Finan 2017) with GRCh37
#           coordinates + DGIdb 5.0 annotation.
# Output  : 04_druggable_mr/results/TaskD_01_druggable_gene_master_v1.csv
# Spec    : WEEK2_TASK_SPEC_v1.0.md, Section 5, Task D.1
#
# DATA SOURCES (verified):
#   Finan 2017 Sci Transl Med 9:eaag1166 (PMID 28356508) - 4,479 druggable genes
#   DGIdb 5.0 Cannon 2024 NAR 52(D1):D1227 (PMID 37953380) - drug-gene annotation
#   Ensembl GRCh37 archive via biomaRt
#
# PREREQUISITE (manual download):
#   1. Finan 2017 supplementary gene list (Table S1) -> data/finan2017_druggable.csv
#      Expected columns: gene_symbol, druggability_tier
#   2. DGIdb 5.0 interactions TSV (optional, for annotation) from dgidb.org
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
})

WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
FINAN_FILE <- file.path(WORK_DIR, "04_druggable_mr/data/finan2017_druggable.csv")
DGIDB_FILE <- file.path(WORK_DIR, "04_druggable_mr/data/dgidb_interactions.tsv")
OUTPUT_FILE <- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_01_druggable_gene_master_v1.csv")
LOG_FILE <- file.path(WORK_DIR, "00_meta/analytical_log.md")

# -----------------------------------------------------------------------------
# Step 1 — Read Finan 2017 druggable gene list
# -----------------------------------------------------------------------------
if (!file.exists(FINAN_FILE)) {
  stop("Finan 2017 druggable gene list not found at: ", FINAN_FILE,
       "\nDownload Table S1 from Finan et al. 2017 Sci Transl Med 9:eaag1166")
}
finan <- fread(FINAN_FILE)
cat("Finan 2017 genes:", nrow(finan), "\n")
# Expected ~4,479 genes. Confirm column names: gene_symbol, druggability_tier

# -----------------------------------------------------------------------------
# Step 2 — Get GRCh37 coordinates via biomaRt
# -----------------------------------------------------------------------------
cat("Connecting to Ensembl GRCh37 archive...\n")
grch37 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl",
                  host = "http://feb2014.archive.ensembl.org")

coords <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name",
                 "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = finan$gene_symbol,
  mart = grch37
)
setDT(coords)

# Keep only standard chromosomes (1-22, X, Y)
coords <- coords[chromosome_name %in% c(as.character(1:22), "X", "Y")]
cat("Genes mapped to GRCh37 coordinates:", nrow(coords), "\n")

# -----------------------------------------------------------------------------
# Step 3 — DGIdb 5.0 annotation (optional, annotation only - does NOT expand set)
# -----------------------------------------------------------------------------
dgidb_annot <- NULL
if (file.exists(DGIDB_FILE)) {
  dgidb <- fread(DGIDB_FILE)
  # DGIdb interactions TSV columns include: gene_name, drug_name, ...
  dgidb_annot <- dgidb[, .(
    dgidb_n_interactions = .N,
    dgidb_drugs = paste(unique(drug_name), collapse = ";")
  ), by = .(gene_symbol = gene_name)]
  cat("DGIdb annotation available for", nrow(dgidb_annot), "genes\n")
} else {
  cat("DGIdb file not found — proceeding without drug annotation (can add later)\n")
}

# -----------------------------------------------------------------------------
# Step 4 — Merge into master list
# -----------------------------------------------------------------------------
master <- merge(coords, finan, by.x = "hgnc_symbol", by.y = "gene_symbol", all.x = TRUE)
setnames(master,
         c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
         c("gene_symbol", "chr", "start_hg19", "end_hg19"))

if (!is.null(dgidb_annot)) {
  master <- merge(master, dgidb_annot, by = "gene_symbol", all.x = TRUE)
} else {
  master[, dgidb_n_interactions := NA_integer_]
  master[, dgidb_drugs := NA_character_]
}

master[, finan_tier := druggability_tier]
master[, in_mr_screen := TRUE]   # all genes proceed to instrument extraction
master[, notes := ""]

# Reorder to spec
out_cols <- c("gene_symbol", "ensembl_gene_id", "chr", "start_hg19", "end_hg19",
              "finan_tier", "dgidb_n_interactions", "dgidb_drugs",
              "in_mr_screen", "notes")
master <- master[, ..out_cols]

# Deduplicate (some symbols map to multiple Ensembl IDs)
master <- unique(master, by = "ensembl_gene_id")

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
fwrite(master, OUTPUT_FILE)
cat("Wrote:", OUTPUT_FILE, "(", nrow(master), "genes )\n")
cat("Genes with coordinates:", sum(!is.na(master$start_hg19)), "\n")

# Sanity: confirm IGF1R, IGF1, INSR, TSHR are present (they should be, as druggable)
for (g in c("TSHR", "IGF1R", "IGF1", "INSR")) {
  present <- g %in% master$gene_symbol
  cat(sprintf("  %s in druggable list: %s\n", g, present))
}

log_entry <- sprintf(
  "## %s — Task D.1\n- Finan 2017 input: %d genes\n- Mapped to GRCh37: %d\n- Final master list: %d genes\n- Output: %s\n\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  nrow(finan), nrow(coords), nrow(master), basename(OUTPUT_FILE)
)
cat(log_entry, file = LOG_FILE, append = TRUE)
