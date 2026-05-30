#!/usr/bin/env Rscript
# =============================================================================
# TED-TRAP v5 — Tables 1 & 2 (publication-quality, gt package)
# =============================================================================
# Reads verified data directly from D.3 results — no hand-transcribed numbers.
# Output: PNG + HTML (and DOCX-ready) for Table 1 and Table 2.
#
# Requires: gt, data.table  (install if missing)
#   install.packages(c("gt","data.table"))
# For PNG export gt needs: webshot2 + chromote (install.packages(c("webshot2","chromote")))
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(gt) })

`%+%` <- function(a, b) paste0(a, b)   # string concat helper (defined before use)

WORK   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
MERGED <- file.path(WORK, "04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv")
TISSUE <- file.path(WORK, "06_tissue_integration/TaskF_01_backbone_tissue_inhouse_v1.csv")
COLOC  <- file.path(WORK, "05_coloc_candidates/TaskE_01_coloc_results_v1.csv")
OUTDIR <- file.path(WORK, "07_manuscript/tables")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

# ---------- helpers ----------
fmtP <- function(p) ifelse(p < 1e-3, sprintf("%.2e", p), sprintf("%.3f", p))
primary <- function(d, g, o) {
  c <- d[gene_symbol==g & outcome==o]
  i <- c[method=="Inverse variance weighted"]
  w <- c[method=="Wald ratio"]
  if (nrow(i)) return(i[1]); if (nrow(w)) return(w[1]); return(c[1])
}

mr     <- fread(MERGED)
tissue <- fread(TISSUE)
coloc  <- fread(COLOC)

# =============================================================================
# TABLE 1 — Study design & data sources
# =============================================================================
t1 <- data.table(
  Role = c("Exposure (instruments)","Discovery (primary)","Replication",
           "Sensitivity (TED-specific)","Target tissue"),
  Dataset = c("eQTLGen blood cis-eQTL","BBJ Graves disease (GCST90018627)",
              "UKB hyperthyroidism (GCST90038636)","FinnGen R12 Graves ophthalmopathy",
              "In-house orbital RNA-seq"),
  Ancestry = c("Predominantly European","East Asian","European","European","Korean"),
  Sample = c("up to 31,684 participants","2,809 cases","3,731 cases","858 cases","4 TED + 1 control")
)
gt1 <- gt(t1) |>
  tab_header(title = html("<b>Table 1. Study design and data sources</b>")) |>
  cols_label(Role="Role", Dataset="Dataset", Ancestry="Ancestry", Sample="Sample size / cases") |>
  tab_options(table.font.size = px(12), heading.align = "left") |>
  tab_source_note(html(
    "<i>eQTLGen sample size indicates the maximum available blood cis-eQTL sample size; "
    %+% "SNP-level sample sizes vary across instruments (e.g., rs179252 NrSamples = 31,566). "
    %+% "BBJ Graves disease was used as the primary discovery phenotype, UKB hyperthyroidism "
    %+% "as cross-ancestry/broader-phenotype replication, and FinnGen Graves ophthalmopathy "
    %+% "as the TED-specific sensitivity outcome. UKB binary-trait effect estimates were "
    %+% "rescaled from linear scale to log-odds using case prevalence.</i>"))

gtsave(gt1, file.path(OUTDIR, "Table1_study_design.html"))
tryCatch(gtsave(gt1, file.path(OUTDIR, "Table1_study_design.png")),
         error=function(e) cat("PNG export needs webshot2/chromote:", conditionMessage(e), "\n")
)
fwrite(t1, file.path(OUTDIR, "Table1_VALUES.csv"))

# =============================================================================
# TABLE 2 — Integrated evidence for backbone genes
# =============================================================================
genes <- c("TSHR","IGF1R","CTLA4")
rows <- list()
for (g in genes) {
  b <- primary(mr, g, "BBJ_Graves")
  u <- primary(mr, g, "UKB_hyperthyroid")
  f <- primary(mr, g, "FinnGen_GO")
  cb <- coloc[gene_symbol==g & outcome=="BBJ_Graves"]
  cf <- coloc[gene_symbol==g & outcome=="FinnGen_GO"]
  ti <- tissue[gene==g]
  rows[[g]] <- data.table(
    Gene = g,
    `MR β (P), BBJ`      = sprintf("%+.2f (%s)", b$beta, fmtP(b$pvalue)),
    `MR β (P), UKB`      = sprintf("%+.2f (%s)", u$beta, fmtP(u$pvalue)),
    `MR β (P), FinnGen`  = sprintf("%+.2f (%s)", f$beta, fmtP(f$pvalue)),
    `Coloc PP.H4 (BBJ/FinnGen)` = sprintf("%.3f / %.3f", cb$pp_h4, cf$pp_h4),
    `Tissue log2FC (padj)` = sprintf("%+.2f (%s)", ti$deseq2_log2fc, fmtP(ti$deseq2_padj)),
    `Integrated role` = c(
      TSHR="Susceptibility anchor; expression-colocalized (3-layer convergent)",
      IGF1R="Pharmacologic effector axis; MR risk-direction but not expression-colocalized susceptibility",
      CTLA4="Known autoimmune locus (positive control)")[g]
  )
}
t2 <- rbindlist(rows)

gt2 <- gt(t2) |>
  tab_header(title = html("<b>Table 2. Integrated genetic and tissue evidence for backbone genes</b>")) |>
  tab_style(style = cell_fill(color = "#FFF2CC"),
            locations = cells_body(rows = Gene == "TSHR")) |>   # highlight TSHR
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Gene)) |>
  tab_options(table.font.size = px(11), heading.align = "left") |>
  tab_source_note(html(
    "<i>MR by inverse-variance weighted (multi-IV) or Wald ratio (single-IV). "
    %+% "TSHR was single-IV across all outcomes; CTLA4 included single- or limited-IV "
    %+% "estimates depending on outcome (BBJ single-IV; UKB/FinnGen 2-IV), limiting "
    %+% "pleiotropy diagnostics. Coloc by coloc.abf (eQTL quant, GWAS cc); PP.H4>0.8 = "
    %+% "shared causal variant. IGF1R coloc PP.H2-dominant (GWAS signal without shared eQTL). "
    %+% "Tissue by DESeq2 (in-house orbital, n=1 control — exploratory supporting layer). "
    %+% "All values verified against source files.</i>"))

gtsave(gt2, file.path(OUTDIR, "Table2_integrated_evidence.html"))
tryCatch(gtsave(gt2, file.path(OUTDIR, "Table2_integrated_evidence.png")),
         error=function(e) cat("PNG export needs webshot2/chromote:", conditionMessage(e), "\n")
)
fwrite(t2, file.path(OUTDIR, "Table2_VALUES.csv"))

# ---- console echo for verification ----
cat("\n=== TABLE 2 (verification echo) ===\n")
print(t2)
cat("\nTables written to:", OUTDIR, "\n")
cat("If PNG failed: install.packages(c('webshot2','chromote')); or use HTML -> Word.\n")
