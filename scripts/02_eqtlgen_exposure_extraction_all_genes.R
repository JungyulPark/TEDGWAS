#==============================================================================
# TED-TRAP Upgrade — Phase 1 Script 02
# Extract cis-eQTL instruments for all candidate genes from eQTLGen
#==============================================================================
# Purpose:
#   - Pull cis-eQTL IVs for 8 candidate genes from eQTLGen (via OpenGWAS)
#   - Clump, filter for F>10, save standardized instrument files
#   - Prepare for multi-method MR in script 03
#
# Candidate genes (pre-specified):
#   TSHR, IGF1R, TNF, PPARG, ARRB1, IRS1, AKT1, CTLA4
#
# Output:
#   TrackA_MR/data/instruments/<GENE>_instruments.tsv
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(data.table)

# --- Candidate genes with ensemble IDs ---
candidate_genes <- data.frame(
  gene      = c("TSHR",     "IGF1R",    "TNF",      "PPARG",
                "ARRB1",    "IRS1",     "AKT1",     "CTLA4"),
  ensembl   = c("ENSG00000165409", "ENSG00000140443", "ENSG00000232810", "ENSG00000132170",
                "ENSG00000137486", "ENSG00000169047", "ENSG00000142208", "ENSG00000163599"),
  stringsAsFactors = FALSE
)
candidate_genes$opengwas_id <- paste0("eqtl-a-", candidate_genes$ensembl)

# --- Parameters ---
# Note: For cis-MR, standard practice is p < 5e-8 (strict) OR p < 5e-6 (permissive)
# We will extract at BOTH thresholds and save separately for sensitivity
p_strict    <- 5e-8
p_permissive <- 5e-6
r2_clump    <- 0.001
kb_clump    <- 10000
pop_clump   <- "EUR"

# --- Output directory ---
inst_dir <- "TrackA_MR/data/instruments"
if (!dir.exists(inst_dir)) dir.create(inst_dir, recursive = TRUE)

# --- Log file ---
log_file <- "TrackA_MR/logs/02_extraction.log"
sink(log_file, split = TRUE)

cat("=== eQTLGen cis-eQTL Instrument Extraction ===\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Target genes (n=%d): %s\n", nrow(candidate_genes),
            paste(candidate_genes$gene, collapse = ", ")))
cat(sprintf("P-value thresholds: strict=%g, permissive=%g\n", p_strict, p_permissive))
cat(sprintf("Clumping: r2 < %g, kb = %d, pop = %s\n\n", r2_clump, kb_clump, pop_clump))

# --- Extraction loop ---
extraction_summary <- list()

for (i in seq_len(nrow(candidate_genes))) {
  gene_sym <- candidate_genes$gene[i]
  gwas_id  <- candidate_genes$opengwas_id[i]

  cat(sprintf("\n--- %s (%s) ---\n", gene_sym, gwas_id))

  # --- Strict threshold (preferred) ---
  strict_ivs <- tryCatch(
    extract_instruments(
      outcomes = gwas_id,
      p1       = p_strict,
      clump    = TRUE,
      r2       = r2_clump,
      kb       = kb_clump
    ),
    error = function(e) {
      cat(sprintf("  [ERROR @ strict]: %s\n", e$message))
      NULL
    }
  )

  n_strict <- if (is.null(strict_ivs) || nrow(strict_ivs) == 0) 0 else nrow(strict_ivs)
  cat(sprintf("  Strict (P<%g): %d IVs\n", p_strict, n_strict))

  # --- Permissive threshold (fallback) ---
  perm_ivs <- tryCatch(
    extract_instruments(
      outcomes = gwas_id,
      p1       = p_permissive,
      clump    = TRUE,
      r2       = r2_clump,
      kb       = kb_clump
    ),
    error = function(e) {
      cat(sprintf("  [ERROR @ permissive]: %s\n", e$message))
      NULL
    }
  )

  n_perm <- if (is.null(perm_ivs) || nrow(perm_ivs) == 0) 0 else nrow(perm_ivs)
  cat(sprintf("  Permissive (P<%g): %d IVs\n", p_permissive, n_perm))

  # --- Compute F-statistic per IV ---
  # F = (beta^2 * N) / (SE^2 * N) ≈ (beta/SE)^2 for large N
  # Or the more precise: F = R^2 * (N-2) / (1-R^2), where R^2 per SNP = 2*MAF*(1-MAF)*beta^2
  compute_F <- function(ivs) {
    if (is.null(ivs) || nrow(ivs) == 0) return(ivs)
    # Using chi-square approximation
    ivs$F_stat <- (ivs$beta.exposure)^2 / (ivs$se.exposure)^2
    ivs
  }
  strict_ivs <- compute_F(strict_ivs)
  perm_ivs   <- compute_F(perm_ivs)

  # Drop weak instruments (F<10)
  if (!is.null(strict_ivs) && nrow(strict_ivs) > 0) {
    n_before <- nrow(strict_ivs)
    strict_ivs <- strict_ivs[strict_ivs$F_stat >= 10, ]
    cat(sprintf("  Strict after F≥10 filter: %d IVs (removed %d)\n",
                nrow(strict_ivs), n_before - nrow(strict_ivs)))
  }
  if (!is.null(perm_ivs) && nrow(perm_ivs) > 0) {
    n_before <- nrow(perm_ivs)
    perm_ivs <- perm_ivs[perm_ivs$F_stat >= 10, ]
    cat(sprintf("  Permissive after F≥10 filter: %d IVs (removed %d)\n",
                nrow(perm_ivs), n_before - nrow(perm_ivs)))
  }

  # --- Save files ---
  if (!is.null(strict_ivs) && nrow(strict_ivs) > 0) {
    out_strict <- file.path(inst_dir, sprintf("%s_instruments_strict.tsv", gene_sym))
    fwrite(strict_ivs, out_strict, sep = "\t")
    cat(sprintf("  Saved strict → %s\n", out_strict))
  }
  if (!is.null(perm_ivs) && nrow(perm_ivs) > 0) {
    out_perm <- file.path(inst_dir, sprintf("%s_instruments_permissive.tsv", gene_sym))
    fwrite(perm_ivs, out_perm, sep = "\t")
    cat(sprintf("  Saved permissive → %s\n", out_perm))
  }

  # --- Decide primary IV set ---
  # Rule: if strict has ≥3, use strict. Else if permissive has ≥3, use permissive with disclosure.
  # For TSHR specifically (known to have few strong IVs), permissive is acceptable per cis-MR norm.
  primary_mode <- if (n_strict >= 3) "strict" else if (n_perm >= 3) "permissive" else "insufficient"
  cat(sprintf("  Primary IV mode: %s\n", primary_mode))

  extraction_summary[[gene_sym]] <- data.frame(
    gene            = gene_sym,
    ensembl         = candidate_genes$ensembl[i],
    opengwas_id     = gwas_id,
    n_ivs_strict    = n_strict,
    n_ivs_permissive = n_perm,
    primary_mode    = primary_mode,
    min_F_strict    = ifelse(n_strict > 0, round(min(strict_ivs$F_stat, na.rm=TRUE), 1), NA),
    min_F_permissive = ifelse(n_perm > 0, round(min(perm_ivs$F_stat, na.rm=TRUE), 1), NA),
    stringsAsFactors = FALSE
  )

  Sys.sleep(1)  # rate-limit courtesy
}

# --- Summary table ---
summary_df <- do.call(rbind, extraction_summary)
cat("\n\n=== Extraction Summary ===\n")
print(summary_df, row.names = FALSE)

summary_file <- "TrackA_MR/results/02_extraction_summary.csv"
fwrite(summary_df, summary_file)
cat(sprintf("\n✅ Summary saved: %s\n", summary_file))

# --- Readiness check for downstream MR ---
cat("\n=== Readiness for MR ===\n")
ready_genes     <- summary_df$gene[summary_df$primary_mode != "insufficient"]
insufficient    <- summary_df$gene[summary_df$primary_mode == "insufficient"]
cat(sprintf("Ready for MR (n=%d): %s\n", length(ready_genes),
            paste(ready_genes, collapse = ", ")))
if (length(insufficient) > 0) {
  cat(sprintf("⚠️  Insufficient IVs (n=%d): %s\n", length(insufficient),
              paste(insufficient, collapse = ", ")))
  cat("   These genes will be dropped from the main MR analysis.\n")
}

sink()
cat("\n✅ Extraction complete. Next: run 03_mr_all_genes_full_sensitivity.R\n\n")
