#!/usr/bin/env Rscript
# =============================================================================
# Task E.1 — Colocalization (coloc.abf) for backbone genes
# =============================================================================
# eQTL (eQTLGen, quant) vs GWAS (cc) per gene x outcome, gene +/-1Mb.
# Output: 05_coloc_candidates/TaskE_01_coloc_results_v1.csv
#
# coloc.abf single-causal-variant assumption (no LD matrix needed).
#   eQTL : type="quant", N=31684, MAF from 1000G EUR
#   GWAS : type="cc", s=case_proportion, N=total
# Refs: Giambartolomei 2014 PLoS Genet 10:e1004383; Wallace 2020 16:e1008720
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(coloc) })

WORK_DIR <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
EQTLGEN  <- "c:/ProjectTEDGWAS/TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
FRQ      <- file.path(WORK_DIR, "04_druggable_mr/data/g1000_eur_freq.frq")
BBJ      <- file.path(WORK_DIR, "04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz")
FINNGEN  <- "c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz"
OUT      <- file.path(WORK_DIR, "05_coloc_candidates/TaskE_02_candidate_coloc_v1.csv")
LOG      <- file.path(WORK_DIR, "00_meta/analytical_log.md")

N_EQTL <- 31684L
WINDOW <- 1000000L

# gene cis windows (hg19) — from TaskD_01 master
GENES <- list(
  TNFSF14 = list(chr=19L, start=6669934L, end=6669934L),
  IFNGR1  = list(chr=6L, start=137239438L, end=137540335L),
  MAPKAPK5= list(chr=12L, start=111994852L, end=112678697L),
  HSD3B7  = list(chr=16L, start=30992922L, end=30997065L),
  PRSS36  = list(chr=16L, start=31142978L, end=31149401L),
  VKORC1  = list(chr=16L, start=31102175L, end=31106118L)
)
PRIORITY_GENES <- c("TNFSF14","HSD3B7","PRSS36","IFNGR1","VKORC1","MAPKAPK5")

# outcome case proportions
OUTCOMES <- list(
  BBJ_Graves = list(file=BBJ,     s=2809/175465, N=175465, fmt="gwascat"),
  FinnGen_GO = list(file=FINNGEN, s=NA,           N=NA,     fmt="finngen")  # s/N set after reading
)

# ---- Load MAF reference once ----
frq <- fread(FRQ)  # CHR SNP A1 A2 MAF NCHROBS
setkey(frq, SNP)

# ---- Load eQTLGen once ----
cat("Loading eQTLGen (several min)...\n")
eqtl_all <- fread(EQTLGEN, select=c("Pvalue","SNP","SNPChr","SNPPos","GeneSymbol","NrSamples"))

# ---- Outcome loaders ----
load_bbj <- function() {
  d <- fread(BBJ)
  setnames(d, tolower(names(d)))
  # cols: chromosome base_pair_location effect_allele other_allele beta standard_error
  #       effect_allele_frequency p_value variant_id ... rsid
  d[, .(SNP=rsid, chr=as.integer(chromosome), pos=as.integer(base_pair_location),
        pval=as.numeric(p_value))]
}
load_finngen <- function() {
  d <- fread(FINNGEN)
  setnames(d, tolower(gsub("^#","",names(d))))
  # cols: chrom pos ref alt rsids ... pval ... af_alt ...
  list(dt=d[, .(SNP=rsids, chr=as.integer(chrom), pos=as.integer(pos),
                pval=as.numeric(pval))],
       n_case=NA, n_ctrl=NA)  # Gemini: fill s from FinnGen meta (858 cases / controls)
}

bbj_dt <- load_bbj()
fin_obj <- load_finngen()
fin_dt  <- fin_obj$dt

# FinnGen s/N: CONFIRM control count. R12 GO ~858 cases. Total from file rows or meta.
# Placeholder — Gemini sets these from FinnGen R12 GRAVES_OPHT documentation:
FINNGEN_S <- 858 / 520387
FINNGEN_N <- 520387

results <- list()

run_coloc <- function(g, ginfo, oc_name, oc_dt, s_val, n_val) {
  c0 <- ginfo$chr; lo <- ginfo$start - WINDOW; hi <- ginfo$end + WINDOW
  e <- eqtl_all[GeneSymbol==g & SNPChr==c0 & SNPPos>=lo & SNPPos<=hi]
  if (nrow(e) < 50) { cat(sprintf("  %s x %s: <50 eQTL SNPs (%d), skip\n", g, oc_name, nrow(e))); return(NULL) }
  o <- oc_dt[chr==c0 & pos>=lo & pos<=hi]
  m <- merge(e[, .(SNP, p_eqtl=Pvalue)], o[, .(SNP, p_gwas=pval)], by="SNP")
  m <- merge(m, frq[, .(SNP, MAF)], by="SNP")
  m <- m[!is.na(MAF) & MAF>0 & MAF<1 & p_eqtl>0 & p_gwas>0]
  m <- m[!duplicated(SNP)]
  if (nrow(m) < 50) { cat(sprintf("  %s x %s: <50 overlap (%d), skip\n", g, oc_name, nrow(m))); return(NULL) }

  res <- tryCatch(
    coloc.abf(
      dataset1=list(snp=m$SNP, pvalues=m$p_gwas, type="cc", s=s_val, N=n_val),
      dataset2=list(snp=m$SNP, pvalues=m$p_eqtl, type="quant", N=N_EQTL),
      MAF=m$MAF),
    error=function(err){cat("  coloc error:",conditionMessage(err),"\n"); NULL})
  if (is.null(res)) return(NULL)

  pp <- res$summary
  topsnp <- res$results[which.max(res$results$SNP.PP.H4), "snp"]
  h3 <- pp["PP.H3.abf"]; h4 <- pp["PP.H4.abf"]
  verdict <- if (h4>0.8) "strong" else if (h3>0.8) "distinct" else if (abs(h3-h4)<0.3 & h4>0.3) "ambiguous" else "weak"

  data.table(gene_symbol=g, outcome=oc_name, n_snps_overlap=nrow(m),
             pp_h0=pp["PP.H0.abf"], pp_h1=pp["PP.H1.abf"], pp_h2=pp["PP.H2.abf"],
             pp_h3=h3, pp_h4=h4, top_snp=topsnp, coloc_verdict=verdict, notes="")
}

for (g in names(GENES)) {
  if (!(g %in% PRIORITY_GENES)) next  # remove this line to run optional chr16 genes
  r1 <- run_coloc(g, GENES[[g]], "BBJ_Graves", bbj_dt, 2809/175465, 175465)
  if (!is.null(r1)) results[[length(results)+1]] <- r1
  r2 <- run_coloc(g, GENES[[g]], "FinnGen_GO", fin_dt, FINNGEN_S, FINNGEN_N)
  if (!is.null(r2)) results[[length(results)+1]] <- r2
}

out <- rbindlist(results)
dir.create(dirname(OUT), recursive=TRUE, showWarnings=FALSE)
fwrite(out, OUT)
cat("\n=== COLOC RESULTS ===\n"); print(out)

cat(sprintf("## %s — Task E.2 candidate coloc\n%s\n\n", format(Sys.time(),"%Y-%m-%d %H:%M"),
            paste(capture.output(print(out)), collapse="\n")), file=LOG, append=TRUE)
