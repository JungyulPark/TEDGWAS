#!/usr/bin/env Rscript
# =============================================================================
# Paper 2 (Track C) — STAGE 0: DATA AUDIT (run FIRST in Antigravity)
# Purpose: report what is runnable BEFORE any MR. Do NOT guess; if a file is
#          missing, print MISSING and continue. No numbers are invented here.
# Output : prints an audit manifest; writes P2_00_audit_manifest.csv
# =============================================================================
suppressWarnings(suppressMessages({}))
options(stringsAsFactors = FALSE)

WORK <- "c:/ProjectTEDGWAS"          # <-- Antigravity local root; edit if needed
cat("== Paper2 STAGE0 audit ==\n WORK =", WORK, "\n\n")

# ---- Target proteins/genes for the off-target axis -------------------------
# Membrane/secreted (plausible plasma pQTL) vs intracellular (often NO plasma pQTL)
targets <- data.frame(
  gene = c("INSR","IGF1R","IRS1","IRS2","PIK3R1","AKT2","PDPK1","INS","IGF1"),
  layer= c("receptor(primary)","receptor(comparator)","adaptor","adaptor",
           "intracellular","intracellular","intracellular","hormone","ligand"),
  plasma_pqtl_expected = c("maybe","maybe","unlikely","unlikely",
                           "unlikely","unlikely","unlikely","likely","likely")
)
cat("Targets:\n"); print(targets); cat("\n")

# ---- 1) pQTL sources (UKB-PPP / deCODE) ------------------------------------
pqtl_dirs <- c(
  ukb_ppp = file.path(WORK,"data/pqtl/UKB-PPP"),
  decode  = file.path(WORK,"data/pqtl/deCODE")
)
cat("[1] pQTL source directories:\n")
for (nm in names(pqtl_dirs)) {
  d <- pqtl_dirs[[nm]]
  cat(sprintf("  %-8s %s : %s\n", nm, d, if (dir.exists(d)) "FOUND" else "MISSING"))
  if (dir.exists(d)) {
    f <- list.files(d, pattern="(INSR|IGF1R|INS_|IGF1)", ignore.case=TRUE)
    cat("     candidate files:", if(length(f)) paste(head(f,8),collapse=", ") else "(none matched targets)","\n")
  }
}

# ---- 2) Outcome GWAS (same as Paper 1) -------------------------------------
outcomes <- c(
  BBJ_Graves        = file.path(WORK,"data/outcomes/GCST90018627_BBJ_Graves.tsv.gz"),
  UKB_hyperthyroid  = file.path(WORK,"data/outcomes/GCST90038636_UKB_hyperthyroid.tsv.gz"),
  FinnGen_GO        = file.path(WORK,"data/outcomes/finngen_R12_GRAVES_OPHT.gz")
)
cat("\n[2] Outcome GWAS files:\n")
for (nm in names(outcomes))
  cat(sprintf("  %-18s %s\n", nm, if (file.exists(outcomes[[nm]])) "FOUND" else "MISSING — edit path"))

# ---- 3) LD reference (EUR-only, locked rule) -------------------------------
ld <- file.path(WORK,"TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/g1000_eur")
cat("\n[3] EUR LD panel:", ld, ":", if (file.exists(paste0(ld,".bim"))) "FOUND" else "MISSING","\n")

# ---- 4) Public TED orbital snRNA-seq (Aim 1) -------------------------------
cat("\n[4] Public snRNA-seq for Aim1 — to be selected. Candidate: Kim 2024 JCI",
    "Insight TED orbital snRNA-seq (GEO). Report accession + n if available locally;",
    "otherwise we will fetch via GEOquery in P2_02.\n")

# ---- 5) FREEBIE: existing INSR eQTL-MR from Paper 1 ------------------------
cat("\n[5] Existing INSR blood eQTL-MR (Paper 1, reuse for Aim 2):\n")
mr <- file.path(WORK,"TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv")
if (file.exists(mr)) {
  d <- read.csv(mr)
  sub <- d[d$gene_symbol %in% c("INSR","IGF1R") &
           d$method %in% c("Inverse variance weighted","Wald ratio"),
           c("gene_symbol","outcome","n_iv","beta","se","pvalue")]
  print(sub, row.names=FALSE)
} else cat("  MISSING:", mr, "\n")

# ---- write manifest --------------------------------------------------------
man <- data.frame(
  item = c(names(pqtl_dirs), names(outcomes), "LD_EUR", "INSR_eQTL_MR_existing"),
  path = c(unlist(pqtl_dirs), unlist(outcomes), ld, mr),
  found= c(dir.exists(pqtl_dirs[1]),dir.exists(pqtl_dirs[2]),
           file.exists(outcomes[1]),file.exists(outcomes[2]),file.exists(outcomes[3]),
           file.exists(paste0(ld,".bim")), file.exists(mr))
)
out <- file.path(WORK,"TrackC_Offtarget/results/P2_00_audit_manifest.csv")
write.csv(man, out, row.names=FALSE)
cat("\n== manifest written:", out, "==\n")
cat("== NEXT: if pQTL + outcomes + LD all FOUND, run P2_01_pQTL_MR.R. ==\n")
cat("== Report this full console output back to Claude before proceeding. ==\n")
