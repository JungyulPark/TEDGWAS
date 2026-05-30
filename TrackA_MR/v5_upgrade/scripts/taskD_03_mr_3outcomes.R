#!/usr/bin/env Rscript
# =============================================================================
# Task D.3 — Druggable-Genome-Wide MR across 3 outcomes (LOCAL outcomes)
# =============================================================================
# Reads SNP-level instruments (D.2a) and runs MR against THREE LOCAL outcome
# files (Option A — no IEU API dependency, full reproducibility):
#   1. BBJ Graves disease       (discovery,   EAS, ~2,809 cases) hg19
#   2. UKB hyperthyroidism       (replication, EUR, ~3,731 cases) hg19
#   3. FinnGen R12 Graves ophth  (TED sensitivity, EUR, 858 cases) local .gz
#
# MODES:
#   DRY_RUN <- TRUE  -> only 5 anchor genes (TSHR/IGF1R/INSR/CTLA4/ARRB1),
#                       prints harmonization diagnostics, writes nothing final.
#   DRY_RUN <- FALSE -> full 2,545-gene run, writes all outputs.
#
# CRITICAL:
#   - Cross-ancestry: eQTLGen instruments EUR-derived; BBJ outcome EAS.
#     Harmonize on allele; flag palindromic; document as limitation.
#   - Bonferroni denominator READ from D.2 (never hardcoded blindly): 2545.
#   - EAF for exposure from 1000G EUR .frq (eQTLGen has no freq column).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)   # Hemani 2018 eLife 7:e34408 (PMID 29846171)
})

# ====================== CONFIG ======================
DRY_RUN <- FALSE   # FULL RUN approved 2026-05-26 after re-dry-run 6/6 pass

WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
SNP_FILE <- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv")
GENE_FILE<- file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02b_eqtl_instrument_gene_summary_v1.csv")
DENOM_FILE<-file.path(WORK_DIR, "04_druggable_mr/results/TaskD_02_denominator_note.txt")
FRQ_FILE <- file.path(WORK_DIR, "04_druggable_mr/data/g1000_eur_freq.frq")

OUT_DIR  <- file.path(WORK_DIR, "04_druggable_mr/results")
LOG_FILE <- file.path(WORK_DIR, "00_meta/analytical_log.md")

# Local outcome files (Option A) -- ADJUST paths to where Gemini downloaded them
OUTCOMES <- list(
  BBJ_Graves = list(
    path = file.path(WORK_DIR, "04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz"),
    ancestry = "EAS", role = "discovery", n_case = 2809),
  UKB_hyperthyroid = list(
    path = file.path(WORK_DIR, "04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz"),
    ancestry = "EUR", role = "replication", n_case = 3731),
  FinnGen_GO = list(
    path = "c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz",
    ancestry = "EUR", role = "sensitivity", n_case = 858)
)

DRY_RUN_GENES <- c("TSHR","IGF1R","INSR","CTLA4","ARRB1")

# Anchor expectations (v4.3) for dry-run validation
ANCHOR_EXPECT <- list(
  TSHR  = "BBJ Wald ~ beta=-1.408 P~9.1e-15 (single-IV)",
  IGF1R = "BBJ IVW P~0.086 null; UKB IVW P~0.052; 4 IVs pre-harmonization",
  ARRB1 = "v4.3 beta=-0.140 P=0.239 null"
)

# ====================== READ INPUTS ======================
snp_level  <- fread(SNP_FILE)
gene_level <- fread(GENE_FILE)

denom_lines <- readLines(DENOM_FILE)
denom_val <- as.integer(sub(".*= *", "", grep("Bonferroni denominator", denom_lines, value=TRUE)[1]))
if (is.na(denom_val) || denom_val < 1) stop("Cannot read Bonferroni denominator from D.2. STOP.")
BONF <- 0.05 / denom_val
cat(sprintf("Bonferroni denominator (from D.2): %d -> threshold %.3e\n", denom_val, BONF))

# EAF from 1000G EUR .frq
frq <- fread(FRQ_FILE)  # cols: CHR SNP A1 A2 MAF NCHROBS
setkey(frq, SNP)

# ====================== EXPOSURE PREP ======================
# eQTLGen gives Z + N; derive beta/se using EUR MAF as EAF proxy.
# Merge MAF; align so eaf corresponds to effect_allele.
prep_exposure <- function(dt) {
  m <- merge(dt, frq[, .(SNP, frq_A1=A1, frq_A2=A2, MAF)], by.x="snp", by.y="SNP", all.x=TRUE)
  # eaf for effect_allele: if effect_allele == frq_A1 -> MAF is freq(A1)=eaf; else 1-MAF
  m[, eaf := fifelse(effect_allele == frq_A1, MAF,
              fifelse(effect_allele == frq_A2, 1 - MAF, NA_real_))]
  # se from Z and N:  se = 1 / sqrt(2*eaf*(1-eaf)*(N + Z^2)); beta = Z*se
  m[, se   := 1 / sqrt(2*eaf*(1-eaf)*(n_samples + zscore^2))]
  m[, beta := zscore * se]
  m[, pval := pvalue]
  m
}

# ====================== OUTCOME LOADER ======================
# Each outcome may have different column names; map per file.
# GWAS Catalog harmonised TSV typically: hm_rsid/variant_id, effect_allele,
#   other_allele, beta, standard_error, p_value, effect_allele_frequency.
# FinnGen R12: #chrom, pos, ref, alt, rsids, pval, beta, sebeta, af_alt ...
# Gemini: confirm exact headers on first load; this maps common cases.
load_outcome <- function(oc) {
  d <- fread(oc$path)
  nm <- tolower(names(d))
  setnames(d, names(d), nm)
  # standardize
  get1 <- function(cands) { hit <- intersect(cands, nm); if(length(hit)) hit[1] else NA }
  snp_c  <- get1(c("hm_rsid","rsid","rsids","variant_id","snp","snpid"))
  ea_c   <- get1(c("effect_allele","alt","hm_effect_allele"))
  oa_c   <- get1(c("other_allele","ref","hm_other_allele"))
  beta_c <- get1(c("beta","hm_beta","effect"))
  se_c   <- get1(c("standard_error","sebeta","se"))
  p_c    <- get1(c("p_value","pval","p","pvalue"))
  eaf_c  <- get1(c("effect_allele_frequency","af_alt","eaf","maf"))
  if (any(is.na(c(snp_c,ea_c,oa_c,beta_c,se_c,p_c)))) {
    cat("OUTCOME HEADER MAP INCOMPLETE for", oc$path, "\n")
    cat("Columns found:", paste(nm, collapse=", "), "\n")
    stop("Map outcome columns manually before proceeding.")
  }
  out <- d[, .(SNP=get(snp_c), effect_allele.outcome=toupper(get(ea_c)),
               other_allele.outcome=toupper(get(oa_c)),
               beta.outcome=as.numeric(get(beta_c)),
               se.outcome=as.numeric(get(se_c)),
               pval.outcome=as.numeric(get(p_c)))]
  if (!is.na(eaf_c)) out[, eaf.outcome := as.numeric(d[[eaf_c]])]
  out <- out[!is.na(SNP) & SNP != ""]

  # --- UKB LMM -> log-odds rescaling (Method B) ---
  # Only for UKB hyperthyroidism (binary trait reported on LMM scale).
  # beta_logodds = beta_LMM / (mu * (1-mu));  same for se.
  # Z-score = beta/se is preserved (same constant), so P is unchanged.
  # BBJ Graves and FinnGen GO are already on log-odds scale — DO NOT rescale.
  if (!is.null(oc$role) && oc$role == "replication" && !is.null(oc$n_case) && oc$n_case == 3731) {
    ukb_mu    <- 3731 / 484598   # case prevalence
    ukb_scale <- ukb_mu * (1 - ukb_mu)  # ≈ 0.00764
    cat(sprintf("  [UKB RESCALE] LMM -> log-odds: mu=%.5f, scale=%.5f\n", ukb_mu, ukb_scale))
    out[, beta.outcome := beta.outcome / ukb_scale]
    out[, se.outcome   := se.outcome   / ukb_scale]
    # pval.outcome unchanged (z preserved)
  }
  out
}

# ====================== MR PER GENE x OUTCOME ======================
run_gene_outcome <- function(gene_sym, exp_g, outcome_dt, oc_name, oc_meta) {
  exp_fmt <- data.frame(
    SNP = exp_g$snp,
    beta.exposure = exp_g$beta,
    se.exposure = exp_g$se,
    effect_allele.exposure = exp_g$effect_allele,
    other_allele.exposure = exp_g$other_allele,
    eaf.exposure = exp_g$eaf,
    pval.exposure = exp_g$pval,
    exposure = gene_sym, id.exposure = gene_sym
  )
  out_sub <- outcome_dt[SNP %in% exp_g$snp]
  if (nrow(out_sub) == 0) return(NULL)
  out_fmt <- as.data.frame(out_sub)
  out_fmt$outcome <- oc_name; out_fmt$id.outcome <- oc_name

  h <- tryCatch(harmonise_data(exp_fmt, out_fmt, action=2), error=function(e) NULL)
  if (is.null(h) || nrow(h) == 0) return(NULL)
  h_keep <- h[h$mr_keep, ]
  n_iv <- nrow(h_keep)
  if (n_iv == 0) return(NULL)

  # method selection
  if (n_iv == 1) {
    res <- mr(h, method_list="mr_wald_ratio")
    egg_i <- NA; egg_p <- NA; q <- NA; q_p <- NA
  } else {
    res <- mr(h, method_list=c("mr_ivw","mr_egger_regression",
                               "mr_weighted_median","mr_weighted_mode"))
    plei <- tryCatch(mr_pleiotropy_test(h), error=function(e) NULL)
    egg_i <- if(!is.null(plei)) plei$egger_intercept else NA
    egg_p <- if(!is.null(plei)) plei$pval else NA
    het <- tryCatch(mr_heterogeneity(h, method_list="mr_ivw"), error=function(e) NULL)
    q   <- if(!is.null(het) && nrow(het)>0) het$Q[1] else NA
    q_p <- if(!is.null(het) && nrow(het)>0) het$Q_pval[1] else NA
  }
  steiger <- tryCatch(directionality_test(h), error=function(e) NULL)

  dt <- as.data.table(res)[, .(method, n_iv=nsnp, beta=b, se, pvalue=pval)]
  dt[, `:=`(gene_symbol=gene_sym, ensembl_gene_id=exp_g$ensembl_gene_id[1],
            outcome=oc_name,
            or=exp(beta), or_ci_lower=exp(beta-1.96*se), or_ci_upper=exp(beta+1.96*se),
            egger_intercept=egg_i, egger_intercept_p=egg_p,
            cochran_q=q, cochran_q_p=q_p,
            steiger_direction=if(!is.null(steiger)) steiger$correct_causal_direction[1] else NA,
            steiger_p=if(!is.null(steiger)) steiger$steiger_pval[1] else NA,
            n_iv_preharm=nrow(exp_g), n_iv_postharm=n_iv,
            n_palindromic=sum(h$palindromic, na.rm=TRUE),
            n_removed=nrow(exp_g)-n_iv,
            outcome_scale_original=fifelse(oc_name=="UKB_hyperthyroid","LMM_binary","log_odds"),
            outcome_scale_used="log_odds",
            ukb_rescaled=(oc_name=="UKB_hyperthyroid"))]
  dt
}

# ====================== MAIN ======================
genes_all <- gene_level[has_valid_instrument == TRUE, gene_symbol]
genes_run <- if (DRY_RUN) intersect(DRY_RUN_GENES, genes_all) else genes_all
cat(sprintf("MODE: %s | genes to run: %d\n", if(DRY_RUN) "DRY_RUN" else "FULL", length(genes_run)))

# Load outcomes once
cat("Loading outcomes...\n")
outcome_loaded <- lapply(OUTCOMES, load_outcome)
names(outcome_loaded) <- names(OUTCOMES)
for (nm in names(outcome_loaded)) cat(sprintf("  %s: %d rows\n", nm, nrow(outcome_loaded[[nm]])))

results <- list()
qc_rows <- list()
for (g in genes_run) {
  exp_g <- prep_exposure(snp_level[gene_symbol == g])
  if (any(is.na(exp_g$eaf))) cat(sprintf("  [warn] %s: %d SNPs missing EAF\n", g, sum(is.na(exp_g$eaf))))
  for (oc_name in names(OUTCOMES)) {
    r <- run_gene_outcome(g, exp_g, outcome_loaded[[oc_name]], oc_name, OUTCOMES[[oc_name]])
    if (!is.null(r)) {
      results[[length(results)+1]] <- r
      if (DRY_RUN) {
        cat(sprintf("  %s x %s: n_iv %d->%d, palindromic=%d | %s P=%.3e beta=%+.3f\n",
            g, oc_name, r$n_iv_preharm[1], r$n_iv_postharm[1], r$n_palindromic[1],
            r$method[which.min(ifelse(r$method=="Inverse variance weighted"|r$method=="Wald ratio",0,1))],
            r$pvalue[1], r$beta[1]))
      }
    }
  }
}

all_res <- rbindlist(results, fill=TRUE)

if (DRY_RUN) {
  cat("\n=== DRY-RUN ANCHOR CHECK ===\n")
  for (g in names(ANCHOR_EXPECT)) {
    cat(sprintf("%s expected: %s\n", g, ANCHOR_EXPECT[[g]]))
    print(all_res[gene_symbol==g & outcome=="BBJ_Graves",
                  .(method, n_iv, beta, pvalue, or)])
  }
  cat("\nDRY-RUN complete. Review anchors (esp. TSHR BBJ beta~-1.408 P~9.1e-15).\n")
  cat("If correct, set DRY_RUN <- FALSE and re-run for full 2,545-gene analysis.\n")
  fwrite(all_res, file.path(OUT_DIR, "TaskD_03_DRYRUN_5gene_v1.csv"))
} else {
  # Split + merged outputs
  fwrite(all_res[outcome=="BBJ_Graves"],      file.path(OUT_DIR,"TaskD_03a_MR_BBJ_discovery_v1.csv"))
  fwrite(all_res[outcome=="UKB_hyperthyroid"],file.path(OUT_DIR,"TaskD_03b_MR_UKB_replication_v1.csv"))
  fwrite(all_res[outcome=="FinnGen_GO"],      file.path(OUT_DIR,"TaskD_03c_MR_FinnGen_GO_sensitivity_v1.csv"))
  fwrite(all_res,                             file.path(OUT_DIR,"TaskD_03d_MR_all_outcomes_merged_v1.csv"))
  # QC
  qc <- all_res[, .(gene_symbol, outcome, method, n_iv_preharm, n_iv_postharm,
                    n_palindromic, n_removed)]
  fwrite(qc, file.path(OUT_DIR,"TaskD_03e_harmonization_qc_v1.csv"))
  cat("FULL D.3 complete. Wrote 03a-03e.\n")
  log_entry <- sprintf("## %s — Task D.3 FULL\n- Genes: %d | rows: %d | Bonf threshold: %.3e\n- Outputs: 03a-03e\n\n",
                       format(Sys.time(),"%Y-%m-%d %H:%M"), length(genes_run), nrow(all_res), BONF)
  cat(log_entry, file=LOG_FILE, append=TRUE)
}
