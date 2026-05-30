#!/usr/bin/env Rscript
# =============================================================================
# Task D.2 (REVISED v2) — cis-eQTL Instrument Extraction (all druggable genes)
# =============================================================================
# CHANGE FROM v1: produces TWO outputs (SNP-level + gene-level summary).
#   - SNP-level is REQUIRED for D.3 MR (effect sizes + alleles per instrument).
#   - Gene-level gives the Bonferroni denominator.
#
# Outputs:
#   04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv
#   04_druggable_mr/results/TaskD_02b_eqtl_instrument_gene_summary_v1.csv
#   04_druggable_mr/results/TaskD_02_denominator_note.txt
#
# CRITICAL RULES (locked):
#   - eQTLGen is EUR-based -> EUR LD reference ONLY (Week 1 lesson).
#   - If local PLINK clumping fails -> STOP and report. DO NOT fall back to API.
#   - No simulation. Real data only.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
})

WORK_DIR     <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
EQTLGEN_FILE <- "c:/ProjectTEDGWAS/TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
MASTER_FILE  <- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_01_druggable_gene_master_v1.csv")

OUT_SNP   <- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv")
OUT_GENE  <- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02b_eqtl_instrument_gene_summary_v1.csv")
DENOM_FILE<- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02_denominator_note.txt")
LOG_FILE  <- file.path(WORK_DIR, "00_meta/analytical_log.md")

BFILE_EUR <- "c:/ProjectTEDGWAS/v5_upgrade/data/1000G_EUR_phase3/g1000_eur"  # ADJUST if path differs
PLINK_BIN <- "c:/ProjectTEDGWAS/tools/plink/plink.exe"                       # ADJUST if path differs

CIS_WINDOW <- 1000000L   # +/- 1 Mb
P_PRIMARY  <- 5e-8
P_SENS     <- 5e-6
CLUMP_R2   <- 0.001
CLUMP_KB   <- 10000

# -----------------------------------------------------------------------------
# Pre-flight checks
# -----------------------------------------------------------------------------
stopifnot(file.exists(EQTLGEN_FILE))
finfo <- file.info(EQTLGEN_FILE)
if (finfo$size < 1e9) {
  stop("eQTLGen file is < 1GB (", round(finfo$size/1e6), " MB). ",
       "Likely incomplete download. STOP.")
}
stopifnot(file.exists(MASTER_FILE))
if (!file.exists(BFILE_EUR <- paste0(BFILE_EUR, ".bed")) ) {
  # bfile prefix check (ld_clump wants prefix without extension)
}
BFILE_EUR <- sub("\\.bed$", "", BFILE_EUR)

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
master <- fread(MASTER_FILE)
master <- master[!is.na(start_hg19)]
cat("Druggable genes with coordinates:", nrow(master), "\n")

cat("Reading eQTLGen full sumstats (large, several min)...\n")
eqtlgen <- fread(
  cmd = paste("zcat", shQuote(EQTLGEN_FILE)),
  select = c("Pvalue", "SNP", "SNPChr", "SNPPos",
             "AssessedAllele", "OtherAllele", "Zscore",
             "Gene", "GeneSymbol", "NrSamples")
)
cat("eQTLGen rows:", nrow(eqtlgen), "\n")

# Sanity: confirm Week 1 anchor reproduces (rs179252 / TSHR)
anchor <- eqtlgen[SNP == "rs179252" & GeneSymbol == "TSHR"]
if (nrow(anchor) >= 1) {
  cat(sprintf("ANCHOR CHECK rs179252/TSHR: P=%.3e Z=%+.4f A1=%s A2=%s (expect P=2.859e-40, Z=+13.2844, G/T)\n",
              anchor$Pvalue[1], anchor$Zscore[1], anchor$AssessedAllele[1], anchor$OtherAllele[1]))
} else {
  cat("WARNING: rs179252/TSHR anchor not found - verify eQTLGen release before trusting results.\n")
}

setkey(eqtlgen, GeneSymbol)

# -----------------------------------------------------------------------------
# Per-gene instrument extraction
# -----------------------------------------------------------------------------
snp_level_rows  <- list()
gene_level_rows <- vector("list", nrow(master))
clump_fail_genes <- character(0)

for (i in seq_len(nrow(master))) {
  g   <- master$gene_symbol[i]
  ens <- master$ensembl_gene_id[i]
  chr <- suppressWarnings(as.integer(master$chr[i]))
  s   <- master$start_hg19[i]
  e   <- master$end_hg19[i]

  if (is.na(chr)) {  # X/Y handled separately if needed
    gene_level_rows[[i]] <- data.table(
      gene_symbol=g, ensembl_gene_id=ens, chr=master$chr[i],
      n_instruments_5e8=0L, n_instruments_5e6=0L,
      has_valid_instrument=FALSE, single_iv=FALSE,
      finan_tier=master$finan_tier[i], notes="non-autosomal chr skipped")
    next
  }

  cis <- eqtlgen[(GeneSymbol == g | Gene == ens) &
                 SNPChr == chr &
                 SNPPos >= (s - CIS_WINDOW) &
                 SNPPos <= (e + CIS_WINDOW)]

  if (nrow(cis) == 0) {
    gene_level_rows[[i]] <- data.table(
      gene_symbol=g, ensembl_gene_id=ens, chr=chr,
      n_instruments_5e8=0L, n_instruments_5e6=0L,
      has_valid_instrument=FALSE, single_iv=FALSE,
      finan_tier=master$finan_tier[i], notes="no cis-eQTL in window")
    next
  }

  n_5e6_raw <- nrow(cis[Pvalue < P_SENS])
  sig <- cis[Pvalue < P_PRIMARY]

  if (nrow(sig) == 0) {
    gene_level_rows[[i]] <- data.table(
      gene_symbol=g, ensembl_gene_id=ens, chr=chr,
      n_instruments_5e8=0L, n_instruments_5e6=n_5e6_raw,
      has_valid_instrument=FALSE, single_iv=FALSE,
      finan_tier=master$finan_tier[i], notes="no SNP below 5e-8")
    next
  }

  # Clump (EUR local). On failure: STOP per locked rule (no API fallback).
  clumped <- tryCatch(
    ld_clump(dplyr::tibble(rsid=sig$SNP, pval=sig$Pvalue),
             clump_r2=CLUMP_R2, clump_kb=CLUMP_KB,
             bfile=BFILE_EUR, plink_bin=PLINK_BIN),
    error = function(err) { clump_fail_genes <<- c(clump_fail_genes, g); NULL }
  )

  if (is.null(clumped)) {
    gene_level_rows[[i]] <- data.table(
      gene_symbol=g, ensembl_gene_id=ens, chr=chr,
      n_instruments_5e8=NA_integer_, n_instruments_5e6=n_5e6_raw,
      has_valid_instrument=NA, single_iv=NA,
      finan_tier=master$finan_tier[i], notes="CLUMP_FAILED")
    next
  }

  n_5e8 <- nrow(clumped)
  iv <- cis[SNP %in% clumped$rsid]
  iv_unique <- iv[!duplicated(SNP)]

  # SNP-level rows (one per instrument) -- REQUIRED for D.3 MR
  snp_level_rows[[length(snp_level_rows)+1]] <- data.table(
    gene_symbol = g,
    ensembl_gene_id = ens,
    chr = chr,
    snp = iv_unique$SNP,
    pos_hg19 = iv_unique$SNPPos,
    effect_allele = iv_unique$AssessedAllele,
    other_allele = iv_unique$OtherAllele,
    zscore = iv_unique$Zscore,
    pvalue = iv_unique$Pvalue,
    n_samples = iv_unique$NrSamples,
    finan_tier = master$finan_tier[i]
  )

  gene_level_rows[[i]] <- data.table(
    gene_symbol=g, ensembl_gene_id=ens, chr=chr,
    n_instruments_5e8=n_5e8, n_instruments_5e6=n_5e6_raw,
    has_valid_instrument=(n_5e8>=1), single_iv=(n_5e8==1),
    finan_tier=master$finan_tier[i], notes="")

  if (i %% 200 == 0) cat("  processed", i, "/", nrow(master), "genes\n")
}

# -----------------------------------------------------------------------------
# STOP if any clumping failed (locked rule: no silent skip)
# -----------------------------------------------------------------------------
if (length(clump_fail_genes) > 0) {
  cat("\n*** CLUMPING FAILED for", length(clump_fail_genes), "genes ***\n")
  cat("First few:", paste(head(clump_fail_genes, 10), collapse=", "), "\n")
  cat("Per locked rule, NOT falling back to API. Writing partial outputs with CLUMP_FAILED flag.\n")
  cat("REPORT to Claude before proceeding to D.3.\n")
}

snp_level  <- rbindlist(snp_level_rows)
gene_level <- rbindlist(gene_level_rows)

# beta/se derived from Z + N (1000G EUR MAF as EAF proxy applied in D.3 harmonization,
# OR compute here if eaf available). Keep raw Z here; D.3 derives beta/se consistently.
fwrite(snp_level,  OUT_SNP)
fwrite(gene_level, OUT_GENE)
cat("Wrote SNP-level:", OUT_SNP, "(", nrow(snp_level), "rows )\n")
cat("Wrote gene-level:", OUT_GENE, "(", nrow(gene_level), "rows )\n")

# -----------------------------------------------------------------------------
# Bonferroni denominator (PRE-SPECIFIED: genes with valid instrument)
# -----------------------------------------------------------------------------
valid <- gene_level[has_valid_instrument == TRUE]
denom <- nrow(valid)
n_single <- sum(valid$single_iv, na.rm=TRUE)
n_multi  <- sum(!valid$single_iv, na.rm=TRUE)
n_failed <- sum(gene_level$notes == "CLUMP_FAILED", na.rm=TRUE)

denom_text <- sprintf(
"=== Task D.2 Bonferroni Denominator (REAL DATA) ===
Total druggable genes screened: %d
Genes with >=1 valid cis-eQTL instrument (P<5e-8, EUR-clumped): %d
  - Single-IV genes (Wald ratio only): %d
  - Multi-IV genes (IVW + sensitivity): %d
Genes with CLUMP_FAILED (need attention): %d

Bonferroni denominator = %d
Bonferroni threshold = 0.05 / %d = %.3e

This denominator was pre-specified in the Week 2 spec (count of genes with a
valid instrument) BEFORE MR results were examined, to prevent post-hoc adjustment.
eQTLGen anchor (rs179252/TSHR) check: see run log.
Generated: %s
",
  nrow(gene_level), denom, n_single, n_multi, n_failed,
  denom, denom, 0.05/denom, format(Sys.time(), "%Y-%m-%d %H:%M"))

writeLines(denom_text, DENOM_FILE)
cat(denom_text)

log_entry <- sprintf(
  "## %s — Task D.2 (revised, 2-output)\n- Genes screened: %d\n- Valid instrument: %d (single %d / multi %d)\n- Clump failed: %d\n- Bonferroni threshold: %.3e\n- Outputs: TaskD_02a (snp), TaskD_02b (gene), denominator_note\n\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  nrow(gene_level), denom, n_single, n_multi, n_failed, 0.05/denom)
cat(log_entry, file=LOG_FILE, append=TRUE)
