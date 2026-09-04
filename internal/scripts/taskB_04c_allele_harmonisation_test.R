# =============================================================================
# Task 3 — Does allele harmonisation change the TSHR fine-mapping result?
# =============================================================================
# taskB_04_susie_finemap.R passes z-scores coded to the eQTLGen effect allele
# together with an LD matrix computed by PLINK --keep-allele-order, i.e. coded
# to the 1000G .bim A1. It never reconciles the two codings.
#
# A uniform flip across all SNPs would be harmless (r_ij picks up f_i*f_j, so
# two flips cancel). The damage comes from SNPs whose flip status disagrees.
# rs179252 (eQTLGen G/T vs bim T/G -> flipped) and rs179248 (C/T vs C/T -> not
# flipped) disagree, so the signed LD between the two European credible-set
# leads is passed to SuSiE with the wrong sign.
#
# This script re-runs SuSiE with z harmonised to the LD reference coding and
# compares credible sets, PIPs, and estimate_s_rss against the unharmonised run.
#
# Usage : Rscript internal/scripts/taskB_04c_allele_harmonisation_test.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

WORK_DIR   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
INPUT_FILE <- file.path(WORK_DIR, "02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv")
TMP_DIR    <- file.path(WORK_DIR, "02_tshr_finemap/data/susie_tmp")
RES_DIR    <- file.path(WORK_DIR, "02_tshr_finemap/results")
RUN_DATE   <- format(Sys.Date(), "%Y%m%d")

BIM <- list(
  EUR = "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/EUR_chr14.bim",
  EAS = "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3/EAS_chr14.bim"
)
PREFIX <- list(EUR = file.path(TMP_DIR, "eur"), EAS = file.path(TMP_DIR, "eas"))

SUSIE_L <- 10; SUSIE_COVERAGE <- 0.95; EQTLGEN_N <- 31684L
ANCHOR <- "rs179252"

cat("R:", R.version.string, " susieR:", as.character(packageVersion("susieR")), "\n")
cat("run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

locus <- fread(INPUT_FILE)
locus <- locus[!is.na(zscore) & !is.na(pos_hg19)]

flip_report <- list()

run_pair <- function(anc) {
  cat("=========================================================\n")
  cat("=== ", anc, " ===\n", sep = "")

  R_ld     <- as.matrix(fread(paste0(PREFIX[[anc]], ".ld")))
  retained <- readLines(paste0(PREFIX[[anc]], ".snplist"))
  rownames(R_ld) <- colnames(R_ld) <- retained

  sub <- locus[snp %in% retained][match(retained, snp)]

  bim <- fread(BIM[[anc]], header = FALSE,
               col.names = c("chr", "snp", "cm", "bp", "a1", "a2"))
  setkey(bim, snp)
  m <- bim[J(sub$snp)]

  # allele pairs must match as a set, otherwise the SNP is not comparable
  same_set <- (sub$effect_allele == m$a1 & sub$other_allele == m$a2) |
              (sub$effect_allele == m$a2 & sub$other_allele == m$a1)
  flip <- ifelse(sub$effect_allele == m$a1, 1, -1)
  flip[!same_set | is.na(same_set)] <- NA

  n_match <- sum(flip == 1, na.rm = TRUE)
  n_flip  <- sum(flip == -1, na.rm = TRUE)
  n_bad   <- sum(is.na(flip))
  cat(sprintf("  SNPs: %d | EA==A1: %d | EA==A2 (flip): %d | allele set mismatch: %d\n",
              nrow(sub), n_match, n_flip, n_bad))
  cat(sprintf("  anchor %s flip = %s\n", ANCHOR,
              ifelse(is.na(flip[sub$snp == ANCHOR]), "NA",
                     ifelse(flip[sub$snp == ANCHOR] == 1, "+1 (EA==A1)", "-1 (EA==A2)"))))

  flip_report[[anc]] <<- data.table(snp = sub$snp, flip = flip)

  z_raw  <- sub$zscore
  z_harm <- sub$zscore * ifelse(is.na(flip), 1, flip)   # NA -> leave as-is

  report <- function(tag, z) {
    fit  <- susie_rss(z = z, R = R_ld, n = EQTLGEN_N, L = SUSIE_L,
                      coverage = SUSIE_COVERAGE)
    pips <- susie_get_pip(fit)
    cs   <- fit$sets$cs
    s    <- tryCatch(estimate_s_rss(z = z, R = R_ld, n = EQTLGEN_N),
                     error = function(e) NA_real_)
    cat(sprintf("\n  [%s]  converged=%s  credible sets=%d  estimate_s_rss=%.6f\n",
                tag, fit$converged, length(cs), s))
    if (length(cs)) {
      for (i in seq_along(cs)) {
        idx <- cs[[i]]
        pur <- if (!is.null(fit$sets$purity)) fit$sets$purity[i, "min.abs.corr"] else NA_real_
        top <- idx[which.max(pips[idx])]
        cat(sprintf("    CS_%d n=%-3d purity=%.4f top=%-12s PIP=%.6f pos=%s\n",
                    i, length(idx), pur, sub$snp[top], pips[top],
                    format(sub$pos_hg19[top], big.mark = ",")))
      }
    }
    ai <- which(sub$snp == ANCHOR)
    if (length(ai)) cat(sprintf("    %s PIP = %.8f\n", ANCHOR, pips[ai]))
    invisible(data.table(ld_reference = anc, coding = tag, snp = sub$snp,
                         pos_hg19 = sub$pos_hg19, z = z, pip = pips,
                         credible_set = {
                           v <- rep(NA_integer_, nrow(sub))
                           if (length(cs)) for (i in seq_along(cs)) v[cs[[i]]] <- i
                           v
                         },
                         converged = fit$converged, estimate_s_rss = s))
  }

  rbind(report("unharmonised (as published)", z_raw),
        report("harmonised to LD coding",     z_harm))
}

out <- rbind(run_pair("EUR"), run_pair("EAS"))

f <- file.path(RES_DIR, paste0("TaskB_04c_allele_harmonisation_", RUN_DATE, ".csv"))
fwrite(out, f)
cat("\nWrote:", f, "\n")

# LD sign actually passed between the two European credible-set leads
R_eur <- as.matrix(fread(paste0(PREFIX$EUR, ".ld")))
ret   <- readLines(paste0(PREFIX$EUR, ".snplist"))
rownames(R_eur) <- colnames(R_eur) <- ret
if (all(c("rs179252", "rs179248") %in% ret)) {
  r <- R_eur["rs179252", "rs179248"]
  fl <- flip_report$EUR
  f1 <- fl[snp == "rs179252", flip]; f2 <- fl[snp == "rs179248", flip]
  cat(sprintf("\nLD passed to SuSiE: r(rs179252, rs179248) = %+.6f  (r2 = %.4f)\n", r, r^2))
  cat(sprintf("flip(rs179252)=%+d  flip(rs179248)=%+d  ->  product = %+d\n", f1, f2, f1 * f2))
  cat(sprintf("sign of r in eQTLGen coding should be %+.6f\n", r * f1 * f2))
}
