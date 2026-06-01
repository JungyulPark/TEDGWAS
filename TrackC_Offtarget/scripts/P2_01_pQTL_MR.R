#!/usr/bin/env Rscript
# =============================================================================
# Paper 2 (Track C) — STAGE 1: pQTL-MR + colocalization for the off-target axis
# Run ONLY after P2_00 audit confirms pQTL + outcomes + EUR LD are present.
# Hypothesis-honest: INSR is expected NULL (supports "off-target, not
# susceptibility"). Report whatever the data show; never invent estimates.
#
# Method (mirrors Paper 1, locked rules):
#   - cis-pQTL instruments: within +/-1 Mb of gene TSS, P < 5e-8
#   - LD clump r2 < 0.001 (10 Mb) using 1000G EUR (EUR-only locked rule)
#   - F = Z^2 (report min F; flag weak instruments F<10)
#   - MR: Wald ratio (1 IV) / IVW primary (+ weighted median, MR-Egger)
#   - coloc.abf: pQTL region vs each outcome (PP.H4 >= 0.8 = shared variant)
#   - Outcomes: BBJ Graves (discovery), UKB hyperthyroid, FinnGen GO (TED)
# Output: P2_01_pQTL_MR_results.csv  +  P2_01_coloc_results.csv
# =============================================================================
suppressWarnings(suppressMessages({
  library(TwoSampleMR); library(coloc); library(data.table); library(dplyr)
}))
options(stringsAsFactors = FALSE)
WORK <- "c:/ProjectTEDGWAS"                      # <-- edit if needed
LD   <- file.path(WORK,"TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/g1000_eur")
PLINK<- "plink"                                  # or full path to plink.exe
OUTD <- file.path(WORK,"TrackC_Offtarget/results"); dir.create(OUTD, showWarnings=FALSE)

# --- targets with plausible plasma pQTL (audit confirms availability) -------
targets <- c("INSR","IGF1R")     # extend with INS/IGF1 if pQTL files exist
gene_tss_hg19 <- list(           # hg19 TSS (approx; used only for +/-1Mb cis window)
  INSR  = c(chr=19, pos=7112266),
  IGF1R = c(chr=15, pos=99192200)
)

# --- outcome loader (harmonize column names to TwoSampleMR) -----------------
outcome_files <- list(
  BBJ_Graves       = file.path(WORK,"data/outcomes/GCST90018627_BBJ_Graves.tsv.gz"),
  UKB_hyperthyroid = file.path(WORK,"data/outcomes/GCST90038636_UKB_hyperthyroid.tsv.gz"),
  FinnGen_GO       = file.path(WORK,"data/outcomes/finngen_R12_GRAVES_OPHT.gz")
)

read_pqtl <- function(gene) {
  # EXPECTED: a per-protein cis-pQTL summary file from UKB-PPP or deCODE.
  # Must expose: SNP, chr, pos(hg19), effect_allele, other_allele, beta, se, pval, eaf
  cand <- Sys.glob(file.path(WORK,"data/pqtl","*",paste0("*",gene,"*")))
  if (!length(cand)) { cat("  [pQTL MISSING]", gene, "- skip\n"); return(NULL) }
  cat("  [pQTL] ", gene, "->", basename(cand[1]), "\n")
  fread(cand[1])   # NOTE: rename columns to the standard set above before MR
}

clump_eur <- function(snps) {
  # local PLINK LD clump on EUR panel; returns kept SNPs (r2<0.001, 10Mb)
  tf <- tempfile(); writeLines(c("SNP P", paste(snps$SNP, snps$pval)), tf)
  of <- tempfile()
  system2(PLINK, c("--bfile",LD,"--clump",tf,"--clump-p1","5e-8",
                   "--clump-r2","0.001","--clump-kb","10000","--out",of),
          stdout=FALSE, stderr=FALSE)
  cf <- paste0(of,".clumped")
  if (!file.exists(cf)) return(character(0))
  fread(cf)$SNP
}

all_mr <- list(); all_coloc <- list()
for (g in targets) {
  cat("\n=== ", g, " ===\n")
  ex <- read_pqtl(g); if (is.null(ex)) next
  # ... user/Gemini: standardize ex columns here (SNP/beta/se/pval/EA/OA/eaf/chr/pos)
  win <- gene_tss_hg19[[g]]
  ex <- ex[ex$chr==win["chr"] & abs(ex$pos-win["pos"])<=1e6 & ex$pval<5e-8, ]
  if (!nrow(ex)) { cat("  no genome-wide cis-pQTL in +/-1Mb — skip\n"); next }
  keep <- clump_eur(ex); ex <- ex[ex$SNP %in% keep, ]
  ex$F <- (ex$beta/ex$se)^2
  cat("  instruments:", nrow(ex), " min F =", round(min(ex$F),1),
      if (min(ex$F)<10) "(WEAK!)" else "", "\n")
  exp_dat <- format_data(as.data.frame(ex), type="exposure",
              snp_col="SNP", beta_col="beta", se_col="se", pval_col="pval",
              effect_allele_col="effect_allele", other_allele_col="other_allele",
              eaf_col="eaf"); exp_dat$exposure <- g

  for (oc in names(outcome_files)) {
    of <- outcome_files[[oc]]; if (!file.exists(of)) { cat("  [outcome MISSING]",oc,"\n"); next }
    out_raw <- fread(of)
    # ... standardize outcome columns to TwoSampleMR before read_outcome_data
    out_dat <- format_data(as.data.frame(out_raw), type="outcome",
                snps=exp_dat$SNP, snp_col="SNP", beta_col="beta", se_col="se",
                pval_col="pval", effect_allele_col="effect_allele",
                other_allele_col="other_allele", eaf_col="eaf")
    dat <- harmonise_data(exp_dat, out_dat, action=2)
    res <- mr(dat, method_list=if(sum(dat$mr_keep)>1)
                c("mr_ivw","mr_weighted_median","mr_egger_regression")
              else "mr_wald_ratio")
    res$gene <- g; res$outcome <- oc; all_mr[[paste(g,oc)]] <- res
    print(res[,c("gene","outcome","method","nsnp","b","se","pval")])

    # colocalization (needs region-wide pQTL + outcome; user supplies region frames)
    # all_coloc[[paste(g,oc)]] <- coloc.abf(dataset1=..., dataset2=...)
  }
}

if (length(all_mr)) {
  mr_tab <- bind_rows(all_mr)
  mr_tab$OR <- exp(mr_tab$b)
  write.csv(mr_tab, file.path(OUTD,"P2_01_pQTL_MR_results.csv"), row.names=FALSE)
  cat("\n== wrote P2_01_pQTL_MR_results.csv ==\n")
} else cat("\n== no MR computed (check audit / pQTL availability) ==\n")

cat("== Report full console output + the two CSVs back to Claude. ==\n")
cat("== Reminder: a NULL INSR result is expected and supportive (off-target). ==\n")
