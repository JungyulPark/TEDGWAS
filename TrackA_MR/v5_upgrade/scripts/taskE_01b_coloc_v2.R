#!/usr/bin/env Rscript
# =============================================================================
# Task E.1b — Colocalization v2
# =============================================================================
# Revision of taskE_01_coloc.R. The original is NOT deleted or modified; it is
# kept as the provenance record of the published v1 results.
#
# What v1 did and why it is being redone:
#   1. Passed p-values, not beta/varbeta, so MAF, N and s entered the posterior
#      through coloc's variance approximation instead of being taken from the
#      studies' own standard errors.
#   2. Applied one 1000G EUR MAF vector to BOTH datasets and BOTH outcomes,
#      including the East Asian BBJ outcome.
#   3. Ran only BBJ and FinnGen. UKB hyperthyroidism — the outcome whose
#      ancestry best matches European eQTLGen, and where IGF1R is nominally
#      significant — had no colocalization at all.
#   4. Reported only PP.H0-H4. SNP.PP.H4 was computed but discarded.
#   5. Used FinnGen N = 520,387, an unreplaced placeholder (see Task 2).
#
# What v2 does:
#   - beta/varbeta for both datasets. eQTLGen beta/se are reconstructed from
#     Z and per-SNP NrSamples with European MAF (Zhu 2016), which is the correct
#     ancestry for eQTLGen; the outcome side uses each study's own beta and
#     standard error, so no MAF is needed there at all and the ancestry-mismatch
#     problem disappears rather than being patched.
#   - All three outcomes, with case/control counts confirmed in Task 4 against
#     the GWAS Catalog and, for FinnGen, against the file's own allele
#     frequencies (Task 2).
#   - UKB betas are rescaled from the BOLT-LMM linear scale to log-odds
#     (Lloyd-Jones 2018) before use, matching what the MR pipeline does.
#   - Explicit allele harmonisation between exposure and outcome, with
#     palindromic and allele-set-mismatch counts reported. (coloc.abf posteriors
#     are sign-invariant, so this is for interpretability and as a guard against
#     variant mismatches -- it is not the Task 3 SuSiE bug, which needed signed
#     LD. See internal/SUSIE_rerun_report.md section 5.5.)
#   - Reports PP.H0-H4, the top SNP.PP.H4 and its value, the H4-conditional 95%
#     credible set, the outcome's minimum P in the window, and overlap counts.
#   - Prior sensitivity over p12 = 1e-6, 5e-6, 1e-5 (Wallace 2020).
#
# Inputs are pre-filtered extracts (see internal/scripts/prepare_coloc_v2_inputs.sh).
#
# Usage : Rscript TrackA_MR/v5_upgrade/scripts/taskE_01b_coloc_v2.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

IN_DIR  <- Sys.getenv("COLOC_V2_IN",
                      "C:/Users/ophjy/AppData/Local/Temp/claude/C--ProjectTEDGWAS/37bf3923-e011-426d-bcdf-103398a6022f/scratchpad/coloc_v2")
RES_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/05_coloc_candidates"
RUN_DATE <- format(Sys.Date(), "%Y%m%d")

N_EQTL <- 31684L          # eQTLGen max; per-SNP NrSamples used for beta/se
WINDOW <- 1000000L
P12_GRID <- c(1e-6, 5e-6, 1e-5)
P12_MAIN <- 1e-5          # coloc default, used for the headline table

GENES <- list(
  TSHR     = list(chr=14L, start=81421333L,  end=81612646L),
  IGF1R    = list(chr=15L, start=99192200L,  end=99507759L),
  CTLA4    = list(chr=2L,  start=204732511L, end=204738683L),
  TNFSF14  = list(chr=19L, start=6669934L,   end=6669934L),
  IFNGR1   = list(chr=6L,  start=137239438L, end=137540335L),
  MAPKAPK5 = list(chr=12L, start=111994852L, end=112678697L),
  HSD3B7   = list(chr=16L, start=30992922L,  end=30997065L),
  VKORC1   = list(chr=16L, start=31102175L,  end=31106118L),
  PRSS36   = list(chr=16L, start=31142978L,  end=31149401L)
)

# Case/control counts confirmed in Task 2 (FinnGen) and Task 4 (BBJ, UKB).
# BBJ  : GWAS Catalog GCST90018627 "2,809 East Asian ancestry cases,
#        172,656 East Asian ancestry controls", 175,465 East Asian (Japan).
# UKB  : GWAS Catalog GCST90038636 "3,731 cases, 480,867 controls",
#        484,598 NR (U.K.). NOTE: this study is Donertas 2021 Nat Aging
#        (PMID 33959723), NOT Sakaue 2021 -- see internal/COLOC_v2_report.md.
# FinnGen: r12.finngen.fi/pheno/GRAVES_OPHT "858 cases / 499490 controls".
OUTCOMES <- list(
  BBJ_Graves = list(
    file = "BBJ_byrsid.tsv", cases = 2809L, controls = 172656L,
    N = 175465L, ancestry = "EAS", scale = "log_odds", label = "BBJ Graves disease"),
  UKB_hyperthyroid = list(
    file = "UKB_byrsid.tsv", cases = 3731L, controls = 480867L,
    N = 484598L, ancestry = "EUR", scale = "linear_lmm", label = "UKB hyperthyroidism"),
  FinnGen_GO = list(
    file = "FinnGen_byrsid.tsv", cases = 858L, controls = 499490L,
    N = 500348L, ancestry = "EUR", scale = "log_odds", label = "FinnGen R12 Graves ophthalmopathy")
)
for (nm in names(OUTCOMES)) OUTCOMES[[nm]]$s <- OUTCOMES[[nm]]$cases / OUTCOMES[[nm]]$N

cat("=========================================================\n")
cat("Task E.1b - colocalization v2\n")
cat("=========================================================\n")
cat("run           :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R             :", R.version.string, "\n")
cat("coloc         :", as.character(packageVersion("coloc")), "\n\n")
for (nm in names(OUTCOMES)) {
  o <- OUTCOMES[[nm]]
  cat(sprintf("  %-17s cases=%6d controls=%7d N=%7d s=%.8f  %s  %s\n",
              nm, o$cases, o$controls, o$N, o$s, o$ancestry, o$scale))
}
cat("\n")

# ---- inputs ----------------------------------------------------------------
eqtl <- fread(file.path(IN_DIR, "eqtlgen_9genes.tsv"))
frq  <- fread(file.path(IN_DIR, "eur_freq_subset.tsv"))
setkey(frq, SNP)
cat("eQTLGen rows:", nrow(eqtl), " | EUR freq rows:", nrow(frq), "\n")

oc_dt <- list()
for (nm in names(OUTCOMES)) {
  d <- fread(file.path(IN_DIR, OUTCOMES[[nm]]$file))
  d <- d[!is.na(beta) & !is.na(se) & se > 0 & rsid != "" & !is.na(rsid)]
  if (OUTCOMES[[nm]]$scale == "linear_lmm") {
    mu <- OUTCOMES[[nm]]$s              # case prevalence
    k  <- mu * (1 - mu)
    d[, `:=`(beta = beta / k, se = se / k)]   # Lloyd-Jones 2018; Z and P preserved
    cat(sprintf("  %s: rescaled LMM -> log-odds, mu=%.6f, divisor=%.8f\n", nm, mu, k))
  }
  oc_dt[[nm]] <- d
  cat(sprintf("  %-17s usable rows %d\n", nm, nrow(d)))
}
cat("\n")

# ---- exposure: eQTLGen beta/se from Z, per-SNP N, EUR MAF -------------------
prep_exposure <- function(dt) {
  m <- merge(dt, frq[, .(SNP, fA1 = A1, fA2 = A2, MAF)],
             by.x = "SNP", by.y = "SNP", all.x = TRUE)
  m[, eaf := fifelse(AssessedAllele == fA1, MAF,
              fifelse(AssessedAllele == fA2, 1 - MAF, NA_real_))]
  m <- m[!is.na(eaf) & eaf > 0 & eaf < 1]
  m[, se_e   := 1 / sqrt(2 * eaf * (1 - eaf) * (NrSamples + Zscore^2))]
  m[, beta_e := Zscore * se_e]
  m
}

flip_report <- list()

harmonise <- function(e, o) {
  m <- merge(e, o, by.x = "SNP", by.y = "rsid", suffixes = c("_e", "_o"))
  if (!nrow(m)) return(m)
  same <- (m$AssessedAllele == m$ea & m$OtherAllele == m$oa)
  swap <- (m$AssessedAllele == m$oa & m$OtherAllele == m$ea)
  m[, keep := same | swap]
  m[, beta_o := fifelse(swap, -beta, beta)]   # align outcome to exposure allele
  m[, palindromic := (AssessedAllele == "A" & OtherAllele == "T") |
                     (AssessedAllele == "T" & OtherAllele == "A") |
                     (AssessedAllele == "C" & OtherAllele == "G") |
                     (AssessedAllele == "G" & OtherAllele == "C")]
  m
}

run_one <- function(g, ginfo, oc_name, p12) {
  e <- eqtl[GeneSymbol == g]
  if (!nrow(e)) return(NULL)
  e <- prep_exposure(e)
  # Match on rsID only. The outcome files are on GRCh38 while eQTLGen and the
  # 1000G panels are on GRCh37, so a position window defined in hg19 would be
  # applied to hg38 coordinates -- an offset of ~466 kb on chr14 and ~865 kb on
  # chr2. v1 did exactly that. The analysis set is in any case the intersection
  # with eQTLGen's own cis window, so no positional filter is needed here.
  o <- oc_dt[[oc_name]][chr == ginfo$chr]
  if (!nrow(o)) return(NULL)

  m <- harmonise(e, o)
  if (!nrow(m)) return(NULL)
  n_pre <- nrow(m); n_pal <- sum(m$palindromic, na.rm = TRUE)
  n_bad <- sum(!m$keep, na.rm = TRUE)
  m <- m[keep == TRUE]
  m <- m[!duplicated(SNP)]
  if (nrow(m) < 50) return(NULL)

  oc <- OUTCOMES[[oc_name]]
  d_gwas <- list(snp = m$SNP, beta = m$beta_o, varbeta = m$se^2,
                 type = "cc", s = oc$s, N = oc$N)
  d_eqtl <- list(snp = m$SNP, beta = m$beta_e, varbeta = m$se_e^2,
                 type = "quant", N = N_EQTL, sdY = 1)

  res <- tryCatch(
    suppressWarnings(coloc.abf(d_gwas, d_eqtl, p12 = p12)),
    error = function(err) { cat("  coloc error:", conditionMessage(err), "\n"); NULL })
  if (is.null(res)) return(NULL)

  r <- as.data.table(res$results)
  setorder(r, -SNP.PP.H4)
  r[, cum := cumsum(SNP.PP.H4)]
  cs_n <- which(r$cum >= 0.95)[1]
  if (is.na(cs_n)) cs_n <- nrow(r)
  cs <- r[seq_len(cs_n)]

  data.table(
    gene = g, outcome = oc_name, p12 = p12,
    n_overlap = nrow(m), n_allele_mismatch = n_bad, n_palindromic = n_pal,
    PP.H0 = res$summary[["PP.H0.abf"]], PP.H1 = res$summary[["PP.H1.abf"]],
    PP.H2 = res$summary[["PP.H2.abf"]], PP.H3 = res$summary[["PP.H3.abf"]],
    PP.H4 = res$summary[["PP.H4.abf"]],
    top_snp = r$snp[1], top_SNP.PP.H4 = r$SNP.PP.H4[1],
    cs95_size = cs_n, cs95_snps = paste(head(cs$snp, 25), collapse = ";"),
    outcome_min_p = min(m$pval, na.rm = TRUE),
    eqtl_min_p = min(m$Pvalue, na.rm = TRUE)
  )
}

out <- list()
for (g in names(GENES)) {
  for (oc_name in names(OUTCOMES)) {
    for (p12 in P12_GRID) {
      r <- run_one(g, GENES[[g]], oc_name, p12)
      if (!is.null(r)) {
        out[[length(out) + 1]] <- r
        if (p12 == P12_MAIN) {
          cat(sprintf("  %-9s x %-17s n=%5d  H2=%.3f H3=%.3f H4=%.3f  top=%s (SNP.PP.H4=%.3f) cs95=%d  minP_out=%.2e\n",
                      g, oc_name, r$n_overlap, r$PP.H2, r$PP.H3, r$PP.H4,
                      r$top_snp, r$top_SNP.PP.H4, r$cs95_size, r$outcome_min_p))
        }
      }
    }
  }
}

res_all <- rbindlist(out)
f <- file.path(RES_DIR, paste0("TaskE_01b_coloc_v2_", RUN_DATE, ".csv"))
fwrite(res_all, f)
cat("\nWrote:", f, "  rows:", nrow(res_all), "\n")
