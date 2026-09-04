# =============================================================================
# Task 3 — TSHR SuSiE re-run with full diagnostics
# =============================================================================
# The original taskB_04_susie_finemap.R reports only per-credible-set summary
# rows. This script reproduces the identical fit and additionally reports what
# the Task 3 brief requires:
#
#   - every SNP's PIP (full table written to CSV)
#   - credible set membership, coverage, purity, convergence
#   - the actual PIP of rs179252
#   - susie_rss LD-vs-summary-statistics consistency diagnostics
#     (estimate_s_rss, kriging_rss) -- the original run emits
#     "IBSS algorithm did not converge" for the East Asian panel
#   - whether the five variants claimed in the manuscript's Table S3 PIP block
#     appear anywhere in the fit
#
# The original script is NOT modified and its v1 output is NOT overwritten.
# Outputs go to new dated filenames.
#
# Usage : Rscript internal/scripts/taskB_04b_susie_rerun_diagnostics.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

WORK_DIR   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
INPUT_FILE <- file.path(WORK_DIR, "02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv")
TMP_DIR    <- file.path(WORK_DIR, "02_tshr_finemap/data/susie_tmp")
RES_DIR    <- file.path(WORK_DIR, "02_tshr_finemap/results")

RUN_DATE       <- format(Sys.Date(), "%Y%m%d")
OUT_PIP        <- file.path(RES_DIR, paste0("TaskB_04_susie_rerun_", RUN_DATE, ".csv"))

# Identical parameters to taskB_04_susie_finemap.R
SUSIE_L        <- 10
SUSIE_COVERAGE <- 0.95
EQTLGEN_N      <- 31684L

ANCHOR <- "rs179252"

# The five variants the manuscript's Table S3 PIP block claims are the
# European credible set (commit e72f84c).
CLAIMED <- c("rs11603529", "rs179248", "rs2284720", "rs10137255", "rs2284722")

cat("=========================================================\n")
cat("TSHR SuSiE re-run with diagnostics\n")
cat("=========================================================\n")
cat("run date      :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version     :", R.version.string, "\n")
cat("susieR version:", as.character(packageVersion("susieR")), "\n")
cat("data.table    :", as.character(packageVersion("data.table")), "\n")
cat("L =", SUSIE_L, " coverage =", SUSIE_COVERAGE, " n =", EQTLGEN_N, "\n")
cat("seed          : not set (susie_rss is deterministic given z, R, n)\n\n")

# -----------------------------------------------------------------------------
locus <- fread(INPUT_FILE)
locus <- locus[!is.na(zscore) & !is.na(pos_hg19)]
cat("locus variants after filtering :", nrow(locus), "\n")
cat("locus span (hg19, chr14)       :",
    format(min(locus$pos_hg19), big.mark = ","), "-",
    format(max(locus$pos_hg19), big.mark = ","),
    sprintf("(%.2f Mb)\n\n", (max(locus$pos_hg19) - min(locus$pos_hg19)) / 1e6))

# -----------------------------------------------------------------------------
run_one <- function(ld_id, tmp_prefix) {
  cat("---------------------------------------------------------\n")
  cat("=== ", ld_id, " ===\n", sep = "")

  ld_file      <- paste0(tmp_prefix, ".ld")
  snplist_file <- paste0(tmp_prefix, ".snplist")
  if (!file.exists(ld_file) || !file.exists(snplist_file)) {
    cat("  LD matrix or snplist missing; run taskB_04_susie_finemap.R first.\n")
    return(NULL)
  }

  R_ld <- as.matrix(fread(ld_file))
  retained <- readLines(snplist_file)
  rownames(R_ld) <- colnames(R_ld) <- retained

  sub <- locus[snp %in% retained]
  sub <- sub[match(retained, snp)]
  cat("  variants after LD intersection :", nrow(sub), "\n")
  cat("  LD matrix                      :", nrow(R_ld), "x", ncol(R_ld), "\n")

  fit <- susie_rss(z = sub$zscore, R = R_ld, n = EQTLGEN_N,
                   L = SUSIE_L, coverage = SUSIE_COVERAGE)

  pips <- susie_get_pip(fit)
  cs   <- fit$sets$cs

  cat("  converged                      :", fit$converged, "\n")
  cat("  credible sets (coverage 0.95)  :", length(cs), "\n")

  cs_of <- rep(NA_integer_, nrow(sub))
  if (length(cs)) {
    for (i in seq_along(cs)) {
      idx <- cs[[i]]
      cs_of[idx] <- i
      pur <- if (!is.null(fit$sets$purity)) fit$sets$purity[i, "min.abs.corr"] else NA_real_
      top <- idx[which.max(pips[idx])]
      cat(sprintf("    CS_%d : n=%d  purity(min.abs.corr)=%.6f  top=%s  PIP=%.6f  pos=%s\n",
                  i, length(idx), pur, sub$snp[top], pips[top],
                  format(sub$pos_hg19[top], big.mark = ",")))
      cat("           members: ", paste(sub$snp[idx], collapse = ", "), "\n")
    }
  }

  # anchor PIP
  ai <- which(sub$snp == ANCHOR)
  if (length(ai)) {
    cat(sprintf("  %s PIP                   : %.10f  (rank %d of %d)\n",
                ANCHOR, pips[ai], as.integer(rank(-pips)[ai]), length(pips)))
  } else {
    cat("  ", ANCHOR, " not present in this panel\n", sep = "")
  }

  # top 20 PIPs
  ord <- order(-pips)[1:min(20, length(pips))]
  cat("\n  --- top 20 PIPs ---\n")
  print(data.frame(
    rank = seq_along(ord),
    snp  = sub$snp[ord],
    pos  = sub$pos_hg19[ord],
    z    = round(sub$zscore[ord], 4),
    pip  = signif(pips[ord], 6),
    cs   = cs_of[ord]
  ), row.names = FALSE)

  # claimed variants
  cat("\n  --- manuscript Table S3 claimed variants ---\n")
  for (v in CLAIMED) {
    j <- which(sub$snp == v)
    if (length(j)) {
      cat(sprintf("    %-12s PRESENT  pos=%-12s PIP=%.6f  CS=%s\n",
                  v, format(sub$pos_hg19[j], big.mark = ","), pips[j],
                  ifelse(is.na(cs_of[j]), "none", cs_of[j])))
    } else {
      cat(sprintf("    %-12s ABSENT from this panel's SuSiE input\n", v))
    }
  }

  # LD-vs-summary-statistics consistency diagnostics
  cat("\n  --- susie_rss diagnostics ---\n")
  s_hat <- tryCatch(estimate_s_rss(z = sub$zscore, R = R_ld, n = EQTLGEN_N),
                    error = function(e) { cat("    estimate_s_rss failed:",
                                              conditionMessage(e), "\n"); NA_real_ })
  cat(sprintf("    estimate_s_rss (s)           : %.6f\n", s_hat))
  cat("      s near 0 = z and LD consistent; s near 1 = strong inconsistency\n")

  cond <- tryCatch({
    k <- kriging_rss(z = sub$zscore, R = R_ld, n = EQTLGEN_N)
    k$conditional_dist
  }, error = function(e) { cat("    kriging_rss failed:",
                              conditionMessage(e), "\n"); NULL })
  if (!is.null(cond)) {
    cond$snp <- sub$snp
    cond$logLR[!is.finite(cond$logLR)] <- NA
    flagged <- cond[!is.na(cond$logLR) & cond$logLR > 2 & abs(cond$z) > 2, ]
    cat(sprintf("    kriging_rss outliers (logLR>2 & |z|>2) : %d\n", nrow(flagged)))
    if (nrow(flagged)) {
      f <- flagged[order(-flagged$logLR), ][1:min(10, nrow(flagged)), ]
      print(data.frame(snp = f$snp, z = round(f$z, 3),
                       cond_mean = round(f$condmean, 3),
                       logLR = round(f$logLR, 3)), row.names = FALSE)
    }
  }

  data.table(
    ld_reference_id = ld_id,
    snp             = sub$snp,
    pos_hg19        = sub$pos_hg19,
    zscore          = sub$zscore,
    pvalue          = sub$pvalue,
    pip             = pips,
    credible_set    = cs_of,
    converged       = fit$converged,
    estimate_s_rss  = s_hat
  )
}

res <- rbind(
  run_one("1000G_phase3_EUR_hg19", file.path(TMP_DIR, "eur")),
  run_one("1000G_phase3_EAS_hg19", file.path(TMP_DIR, "eas"))
)

if (!is.null(res)) {
  fwrite(res[order(ld_reference_id, -pip)], OUT_PIP)
  cat("\n---------------------------------------------------------\n")
  cat("Wrote full PIP table:", OUT_PIP, "\n")
  cat("rows:", nrow(res), "\n")
}
