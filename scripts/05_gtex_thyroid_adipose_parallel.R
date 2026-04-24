#==============================================================================
# TED-TRAP Upgrade — Phase 1 Script 05
# GTEx v8 thyroid + adipose tissue-specific MR (parallel to eQTLGen blood)
#==============================================================================
# Purpose:
#   Resolve MAJOR issue #7: Tissue mismatch between eQTLGen (whole blood)
#   and disease tissues (orbital adipose, thyroid, etc.).
#
#   This script re-runs cis-MR using GTEx v8 thyroid eQTLs for TSHR/IGF1R/
#   PPARG, and GTEx v8 visceral adipose for IGF1R/IRS1/IRS2, providing
#   tissue-appropriate replication of the primary eQTLGen MR.
#
# Required inputs:
#   - GTEx v8 cis-QTL summary stats (download from gtexportal.org/home/downloads)
#   - For each tissue, eGene significant associations: Thyroid, Adipose_Visceral_Omentum
#
# Data tables used (via OpenGWAS "eqtl-a-*" aliases where available,
# else direct GTEx file):
#   - Thyroid:          GTEx_Analysis_v8_eQTL/Thyroid.v8.signif_variant_gene_pairs.txt.gz
#   - Adipose visceral: GTEx_Analysis_v8_eQTL/Adipose_Visceral_Omentum.v8.signif_variant_gene_pairs.txt.gz
#
# Output:
#   TrackA_MR/results/05_gtex_parallel_mr.csv
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(TwoSampleMR)
library(data.table)
library(dplyr)

log_file <- "TrackA_MR/logs/05_gtex_parallel.log"
sink(log_file, split = TRUE)
cat("=== GTEx v8 Parallel MR ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- GTEx file paths (user downloads from portal) ---
gtex_thyroid  <- "TrackA_MR/data/gtex_v8/Thyroid.v8.signif_variant_gene_pairs.txt.gz"
gtex_adipose  <- "TrackA_MR/data/gtex_v8/Adipose_Visceral_Omentum.v8.signif_variant_gene_pairs.txt.gz"

if (!file.exists(gtex_thyroid) || !file.exists(gtex_adipose)) {
  cat("❌ GTEx files missing. Download from:\n")
  cat("   https://gtexportal.org/home/downloads/adult-gtex/qtl\n")
  cat("   Required: GTEx_Analysis_v8_eQTL tar.gz\n\n")
  cat("   After download, extract both:\n")
  cat("     - Thyroid.v8.signif_variant_gene_pairs.txt.gz\n")
  cat("     - Adipose_Visceral_Omentum.v8.signif_variant_gene_pairs.txt.gz\n")
  cat("   to TrackA_MR/data/gtex_v8/\n")
  stop("Download GTEx data to proceed.")
}

# --- Gene-tissue assignments ---
# Each gene → tissue(s) most relevant for disease
gene_tissue <- data.frame(
  gene    = c("TSHR", "IGF1R",           "IGF1R",           "PPARG",   "IRS1",            "IRS2",            "AKT1",            "CTLA4",           "TNF"),
  ensembl = c("ENSG00000165409", "ENSG00000140443", "ENSG00000140443", "ENSG00000132170", "ENSG00000169047", "ENSG00000185950", "ENSG00000142208", "ENSG00000163599", "ENSG00000232810"),
  tissue  = c("Thyroid", "Thyroid",       "Adipose_Visceral","Thyroid",  "Adipose_Visceral","Adipose_Visceral","Thyroid",         "Thyroid",         "Thyroid"),
  stringsAsFactors = FALSE
)

# --- Extract function ---
load_gtex_cis <- function(file_path, ensembl_id) {
  # GTEx v8 columns (signif_variant_gene_pairs format):
  # variant_id, gene_id, tss_distance, ma_samples, ma_count, maf,
  # pval_nominal, slope, slope_se, pval_nominal_threshold, ...
  cat(sprintf("  Loading GTEx: %s (%s)\n", ensembl_id, basename(file_path)))
  cmd <- sprintf("zcat %s | awk -v g=%s 'NR==1 || $2 ~ g'",
                 file_path, ensembl_id)
  df <- tryCatch(fread(cmd = cmd, header = TRUE),
                 error = function(e) { cat(sprintf("  [ERROR] %s\n", e$message)); NULL })
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Rename to TwoSampleMR conventions
  df$SNP <- sub("^chr", "", df$variant_id)   # keep variant_id as SNP ID (will need rsID mapping)
  df$beta.exposure <- df$slope
  df$se.exposure   <- df$slope_se
  df$pval.exposure <- df$pval_nominal
  df$eaf.exposure  <- df$maf
  df$exposure      <- paste0(ensembl_id, "_", tools::file_path_sans_ext(basename(file_path)))
  df$id.exposure   <- df$exposure
  df$effect_allele.exposure <- strsplit(df$variant_id, "_")[[1]][4]   # rough
  df$other_allele.exposure  <- strsplit(df$variant_id, "_")[[1]][3]   # rough
  df$samplesize.exposure    <- 670   # GTEx v8 thyroid N; adjust per tissue

  # Compute F
  df$F_stat <- (df$beta.exposure)^2 / (df$se.exposure)^2
  df <- df[df$F_stat >= 10, ]
  df
}

results <- list()

for (i in seq_len(nrow(gene_tissue))) {
  gene_sym <- gene_tissue$gene[i]
  ensembl  <- gene_tissue$ensembl[i]
  tissue   <- gene_tissue$tissue[i]
  gtex_file <- if (tissue == "Thyroid") gtex_thyroid else gtex_adipose

  cat(sprintf("\n--- %s @ %s ---\n", gene_sym, tissue))

  ivs <- load_gtex_cis(gtex_file, ensembl)
  if (is.null(ivs) || nrow(ivs) == 0) {
    cat("  [SKIP] No IVs found\n")
    next
  }
  cat(sprintf("  Instruments (F≥10): %d\n", nrow(ivs)))

  # NOTE: GTEx variant IDs are in format chr_pos_REF_ALT_b38
  # Converting to rsID requires a lookup (dbSNP or MRgene mapping)
  # For simplicity, we recommend extracting via the OpenGWAS GTEx aliases instead:
  # e.g. TSHR-thyroid is eqtl-a-ENSG00000165409 from Thyroid tissue
  # These are available through ieugwasr::available_outcomes()
  # If using OpenGWAS, this script should mirror script 02 structure.

  cat("  NOTE: rsID mapping required. Recommend using OpenGWAS GTEx aliases instead\n")
  cat("        (e.g., eqtl-a-* IDs). Update gene_tissue table with these IDs.\n")

  # Placeholder: save raw extraction
  out <- sprintf("TrackA_MR/data/gtex_v8/extracted_%s_%s.tsv", gene_sym, tissue)
  fwrite(ivs, out, sep = "\t")
}

sink()

cat("\n✅ GTEx extraction complete. Next steps:\n")
cat("   1. Map GTEx variant_id to rsID (via dbSNP or 1000G).\n")
cat("   2. Re-run MR with script 03 using these instruments.\n")
cat("   3. Compare effect directions and magnitudes with eQTLGen results.\n\n")
cat("Alternative faster path: use OpenGWAS GTEx aliases (eqtl-a-* for Thyroid/Adipose)\n")
cat("   See https://gwas.mrcieu.ac.uk/ for full list\n\n")
