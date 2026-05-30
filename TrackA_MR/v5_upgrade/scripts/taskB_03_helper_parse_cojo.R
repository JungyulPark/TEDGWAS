#!/usr/bin/env Rscript
# =============================================================================
# Task B.3 Helper — Parse GCTA-COJO raw output into spec-formatted TSV
# =============================================================================
# Usage  : Rscript taskB_03_helper_parse_cojo.R <output_dir>
# Reads  : <dir>/TaskB_03_tshr_cojo_EUR_raw.jma.cojo
#          <dir>/TaskB_03_tshr_cojo_EAS_raw.jma.cojo
# Writes : <dir>/TaskB_03_tshr_cojo_independent_signals_v1.tsv
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript taskB_03_helper_parse_cojo.R <output_dir>")
}

OUT_DIR <- args[1]

suppressPackageStartupMessages({
  library(data.table)
  library(ieugwasr)
})

# -----------------------------------------------------------------------------
# Read both COJO outputs
# -----------------------------------------------------------------------------
# COJO .jma.cojo columns: Chr SNP bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r
# bJ = joint effect after conditioning on other selected SNPs

parse_cojo <- function(file, ld_id) {
  if (!file.exists(file)) {
    warning("COJO output not found: ", file)
    return(NULL)
  }
  d <- fread(file)
  if (nrow(d) == 0) return(NULL)

  out <- data.frame(
    ld_reference_id = ld_id,
    signal_rank = seq_len(nrow(d)),
    lead_snp = d$SNP,
    chr = d$Chr,
    pos_hg19 = d$bp,
    effect_allele = d$refA,
    other_allele = NA_character_,  # COJO doesn't output A2 directly; merge from .ma later
    eaf = d$freq,
    marginal_beta = d$b,
    marginal_p = d$p,
    conditional_beta = d$bJ,
    conditional_se = d$bJ_se,
    conditional_p = d$pJ,
    stringsAsFactors = FALSE
  )
  out
}

eur_results <- parse_cojo(
  file.path(OUT_DIR, "TaskB_03_tshr_cojo_EUR_raw.jma.cojo"),
  "1000G_phase3_EUR_hg19"
)
eas_results <- parse_cojo(
  file.path(OUT_DIR, "TaskB_03_tshr_cojo_EAS_raw.jma.cojo"),
  "1000G_phase3_EAS_hg19"
)

combined <- rbind(eur_results, eas_results)

# -----------------------------------------------------------------------------
# Compute LD with rs179252 in both ancestries
# -----------------------------------------------------------------------------
# Compute LD locally via PLINK rather than relying on online OpenGWAS APIs
get_r2_with_rs179252 <- function(target_snp, pop) {
  if (target_snp == "rs179252") return(1.0)
  
  bfile <- if (pop == "EUR") {
    "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/EUR_chr14"
  } else {
    "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3/EAS_chr14"
  }
  
  temp_extract <- tempfile(fileext = ".txt")
  writeLines(c("rs179252", target_snp), temp_extract)
  
  temp_out <- tempfile()
  cmd <- sprintf('c:/ProjectTEDGWAS/tools/plink.exe --bfile %s --r2 --ld-snp rs179252 --extract %s --ld-window-r2 0 --ld-window 99999 --ld-window-kb 99999 --out %s', bfile, temp_extract, temp_out)
  
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  ld_file <- paste0(temp_out, ".ld")
  r2_val <- NA_real_
  if (file.exists(ld_file)) {
    ld_data <- fread(ld_file)
    target_row <- ld_data[SNP_B == target_snp]
    if (nrow(target_row) > 0) {
      r2_val <- target_row$R2[1]
    }
    unlink(c(ld_file, paste0(temp_out, ".log")))
  }
  unlink(temp_extract)
  return(r2_val)
}

combined$ld_with_rs179252_eur <- vapply(combined$lead_snp,
                                         get_r2_with_rs179252,
                                         numeric(1), pop = "EUR")
combined$ld_with_rs179252_eas <- vapply(combined$lead_snp,
                                         get_r2_with_rs179252,
                                         numeric(1), pop = "EAS")


# -----------------------------------------------------------------------------
# Interpretation column
# -----------------------------------------------------------------------------
combined$interpretation <- with(combined, ifelse(
  signal_rank == 1, "primary",
  ifelse(pmax(ld_with_rs179252_eur, ld_with_rs179252_eas, na.rm = TRUE) > 0.6,
         "collinear with rs179252",
         "secondary independent")
))

# -----------------------------------------------------------------------------
# Reorder to spec and merge other_allele
# -----------------------------------------------------------------------------
ma_file <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/02_tshr_finemap/data/TSHR_eqtlgen_locus.ma"
if (file.exists(ma_file)) {
  ma_data <- fread(ma_file, select = c("SNP", "A2"))
  combined <- merge(combined, ma_data, by.x = "lead_snp", by.y = "SNP", all.x = TRUE, sort = FALSE)
  combined$other_allele <- combined$A2
  combined$A2 <- NULL
}

ordered_cols <- c("ld_reference_id", "signal_rank", "lead_snp", "chr", "pos_hg19",
                  "effect_allele", "other_allele", "eaf",
                  "marginal_beta", "marginal_p",
                  "conditional_beta", "conditional_se", "conditional_p",
                  "ld_with_rs179252_eur", "ld_with_rs179252_eas",
                  "interpretation")
combined <- combined[, ordered_cols]

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
out_file <- file.path(OUT_DIR, "TaskB_03_tshr_cojo_independent_signals_v1.tsv")
fwrite(combined, out_file, sep = "\t")

cat("Wrote:", out_file, "(", nrow(combined), "rows )\n")
cat("\nSummary:\n")
print(table(combined$ld_reference_id, combined$interpretation))
