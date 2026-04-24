#==============================================================================
# TED-TRAP Upgrade — Phase 3 Script 17
# MAGMA gene-level analysis on FinnGen R12 Graves ophthalmopathy GWAS
#==============================================================================
# Purpose:
#   Use MAGMA (de Leeuw et al. 2015, PLoS Comput Biol 11:e1004219) to perform
#   gene-based association testing on FinnGen R12 Graves ophthalmopathy
#   summary statistics. Genes with MAGMA p < 0.05 / 20000 (Bonferroni) OR
#   p < 1e-4 (suggestive) are strong candidates for the TED disease gene set.
#
#   This provides a purely GWAS-driven, reproducible gene list that
#   complements the other source (Open Targets, DisGeNET, PubTator).
#
# Required inputs:
#   - FinnGen R12 GO summary statistics (SNP-level P values)
#   - MAGMA v1.10+ executable (https://ctg.cncr.nl/software/magma)
#   - Reference files: MAGMA-provided EUR 1000G LD reference
#   - Gene location file: MAGMA NCBI38.gene.loc
#
# This script writes a shell runner that calls MAGMA binaries.
# MAGMA itself is a command-line tool, not an R package; R is used for
# pre/post-processing summary statistics.
#
# Output:
#   TrackB_Network/results/17_magma_finngen_go.genes.out (MAGMA output)
#   TrackB_Network/results/17_magma_TED_genes.csv (filtered gene list)
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(data.table)

log_file <- "TrackA_MR/logs/17_magma.log"
sink(log_file, split = TRUE)
cat("=== MAGMA Gene-Level Analysis ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# --- Required file paths ---
finngen_ss     <- "TrackA_MR/data/finngen/finngen_R12_E4_GRAVES_OPHT.tsv.gz"
magma_bin      <- "TrackA_MR/tools/magma/magma"    # Linux binary; or .exe on Windows
ref_prefix     <- "TrackA_MR/tools/magma/g1000_eur"     # 1000G EUR reference
gene_loc_file  <- "TrackA_MR/tools/magma/NCBI38.gene.loc"

for (f in c(magma_bin, paste0(ref_prefix, ".bed"), gene_loc_file, finngen_ss)) {
  if (!file.exists(f)) {
    cat(sprintf("❌ MISSING: %s\n", f))
  }
}

cat("\nPre-requirements:\n")
cat("  1. Download MAGMA: https://cncr.nl/research/magma/\n")
cat("  2. Download reference: 1000G European (https://ctg.cncr.nl/software/magma)\n")
cat("  3. Download NCBI38.gene.loc file\n")
cat("  4. Place in TrackA_MR/tools/magma/\n\n")

# --- Step 1: Prepare SNP-level P values ---
cat("--- Step 1: Prepare SNP file ---\n")
if (file.exists(finngen_ss)) {
  fg <- fread(finngen_ss)
  cat(sprintf("FinnGen GO: %d rows\n", nrow(fg)))

  # MAGMA expects: SNP, CHR, BP, P, N
  snp_table <- data.frame(
    SNP = fg$rsids,
    CHR = fg$`#chrom`,
    BP  = fg$pos,
    P   = fg$pval,
    N   = 500348,   # FinnGen R12 GO total
    stringsAsFactors = FALSE
  )
  snp_table <- snp_table[!is.na(snp_table$SNP) & snp_table$SNP != "", ]

  snp_file <- "TrackB_Network/results/17_magma_snp_pvals.tsv"
  fwrite(snp_table, snp_file, sep = "\t")
  cat(sprintf("  Saved: %s (%d SNPs)\n", snp_file, nrow(snp_table)))
} else {
  cat("Summary stats not yet downloaded.\n")
  cat("Get from: https://r12.finngen.fi/pheno/E4_GRAVES_OPHT\n\n")
}

# --- Step 2: Write MAGMA shell command ---
cat("\n--- Step 2: MAGMA shell commands ---\n")
cmds <- sprintf('
# Step 2a: SNP-to-gene annotation
%s --annotate \\
  --snp-loc %s \\
  --gene-loc %s \\
  --out TrackB_Network/results/17_magma_annot

# Step 2b: Gene-based analysis
%s --bfile %s \\
  --gene-annot TrackB_Network/results/17_magma_annot.genes.annot \\
  --pval TrackB_Network/results/17_magma_snp_pvals.tsv ncol=N \\
  --out TrackB_Network/results/17_magma_finngen_go
',
  magma_bin, "TrackB_Network/results/17_magma_snp_pvals.tsv", gene_loc_file,
  magma_bin, ref_prefix
)
cat(cmds, "\n")

shell_file <- "TrackB_Network/scripts/17_run_magma.sh"
writeLines(cmds, shell_file)
Sys.chmod(shell_file, mode = "0755")
cat(sprintf("Shell script saved: %s\n", shell_file))
cat("Run manually: bash TrackB_Network/scripts/17_run_magma.sh\n\n")

# --- Step 3: Post-process MAGMA output (if exists) ---
genes_out_file <- "TrackB_Network/results/17_magma_finngen_go.genes.out"
if (file.exists(genes_out_file)) {
  cat("--- Step 3: Post-processing MAGMA gene results ---\n")
  genes <- fread(genes_out_file)
  cat(sprintf("MAGMA gene results: %d genes tested\n", nrow(genes)))
  print(head(genes))

  # Bonferroni threshold: 0.05 / 20000 = 2.5e-6
  bonf_thresh <- 0.05 / 20000
  # Suggestive: 1e-4
  sugg_thresh <- 1e-4

  bonf_genes <- genes[P < bonf_thresh, ]
  sugg_genes <- genes[P < sugg_thresh & P >= bonf_thresh, ]

  cat(sprintf("\nBonferroni-significant genes (P < %g): %d\n",
              bonf_thresh, nrow(bonf_genes)))
  cat(sprintf("Suggestive genes (P < %g): %d\n\n", sugg_thresh, nrow(sugg_genes)))

  # Combine and annotate
  final <- rbind(
    cbind(bonf_genes, tier = "Bonferroni"),
    cbind(sugg_genes, tier = "Suggestive")
  )

  out_csv <- "TrackB_Network/results/17_magma_TED_genes.csv"
  fwrite(final, out_csv)
  cat(sprintf("✅ Saved: %s\n", out_csv))
} else {
  cat("⚠️  MAGMA output not yet generated. Run the shell script first.\n")
}

sink()
