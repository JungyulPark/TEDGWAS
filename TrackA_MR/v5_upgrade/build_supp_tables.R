#!/usr/bin/env Rscript
# =============================================================================
# TED-TRAP v5 — Supplementary Tables (S1, S2, S5, S7)
# =============================================================================
# Output: CSV, HTML, and PNG formats for Tables S1, S2, S5, and S7.
# Saved under TrackA_MR/v5_upgrade/07_manuscript/supp_tables/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(gt)
})

`%+%` <- function(a, b) paste0(a, b)   # string concat helper

WORK   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
OUTDIR <- file.path(WORK, "07_manuscript/supp_tables")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

# Source data paths
S1_SRC <- file.path(WORK, "04_druggable_mr/results/TaskD_04_hit_classification_VERIFIED.csv")
S2_SRC <- file.path(WORK, "05_coloc_candidates/TaskE_02_candidate_coloc_v1.csv")
S7_SRC <- file.path(WORK, "04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv")
S4_SRC <- file.path(WORK, "04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv")

# Format P-value helper
fmtP <- function(p) {
  if (is.na(p) || p == "None" || p == "") return("-")
  p_num <- as.numeric(p)
  if (is.na(p_num)) return(p)
  if (p_num < 1e-3) sprintf("%.2e", p_num) else sprintf("%.3f", p_num)
}

fmtBeta <- function(b) {
  if (is.na(b) || b == "None" || b == "") return("-")
  b_num <- as.numeric(b)
  if (is.na(b_num)) return(b)
  sprintf("%+.2f", b_num)
}

# =============================================================================
# TABLE S1 — 13 Bonferroni hits classification
# =============================================================================
s1_raw <- fread(S1_SRC)

# Convert columns to numeric safely for direction/nominal support calculations
bbj_beta_num <- as.numeric(s1_raw$bbj_beta)
ukb_beta_num <- as.numeric(s1_raw$ukb_beta)
ukb_p_num <- as.numeric(s1_raw$ukb_p)
fin_beta_num <- as.numeric(s1_raw$fin_beta)
fin_p_num <- as.numeric(s1_raw$fin_p)

# Same direction = sign matches (and neither is NA)
ukb_same <- fifelse(is.na(ukb_beta_num) | is.na(bbj_beta_num), NA_character_,
                    fifelse(sign(ukb_beta_num) == sign(bbj_beta_num), "Yes", "No"))
ukb_nom  <- fifelse(is.na(ukb_same), NA_character_,
                    fifelse(ukb_same == "Yes" & !is.na(ukb_p_num) & ukb_p_num < 0.05, "Yes", "No"))

fin_same <- fifelse(is.na(fin_beta_num) | is.na(bbj_beta_num), NA_character_,
                    fifelse(sign(fin_beta_num) == sign(bbj_beta_num), "Yes", "No"))
fin_nom  <- fifelse(is.na(fin_same), NA_character_,
                    fifelse(fin_same == "Yes" & !is.na(fin_p_num) & fin_p_num < 0.05, "Yes", "No"))

s1_data <- s1_raw[, .(
  Gene = trimws(gene),
  Chr = chr,
  `Position (Mb)` = pos_mb,
  `Instruments (N)` = n_iv,
  Method = method,
  `BBJ Discovery Beta (P)` = paste0(sapply(bbj_beta, fmtBeta), " (", sapply(bbj_p, fmtP), ")"),
  `UKB Replication Beta (P)` = paste0(sapply(ukb_beta, fmtBeta), " (", sapply(ukb_p, fmtP), ")"),
  `UKB same dir` = fifelse(is.na(ukb_same), "NA", ukb_same),
  `UKB nominal` = fifelse(is.na(ukb_nom), "NA", ukb_nom),
  `FinnGen Sensitivity Beta (P)` = paste0(sapply(fin_beta, fmtBeta), " (", sapply(fin_p, fmtP), ")"),
  `FinnGen same dir` = fifelse(is.na(fin_same), "NA", fin_same),
  `FinnGen nominal` = fifelse(is.na(fin_nom), "NA", fin_nom),
  Class = class
)]

# Clean up any "NA ( - )" formats
for (col in c("BBJ Discovery Beta (P)", "UKB Replication Beta (P)", "FinnGen Sensitivity Beta (P)")) {
  s1_data[[col]] <- gsub("NA \\( - \\)", "-", s1_data[[col]])
}

# Pretty class names
class_map <- c(
  known_anchor = "Known anchor",
  MHC_spillover = "MHC spillover",
  chr16p11_LDcluster = "chr16p11.2 LD cluster",
  single_IV_unverifiable = "Single-IV (unverifiable)",
  discovery_only = "Discovery only"
)
s1_data[, Class := class_map[Class]]
s1_data[Gene == "CTLA4", Class := "Known anchor (positive control)"]

# Save S1 CSV
fwrite(s1_data, file.path(OUTDIR, "TableS1_mr_hits.csv"))

# Build S1 HTML/PNG
gt_s1 <- gt(s1_data) |>
  tab_header(title = html("<b>Table S1. Multi-outcome classification of 13 druggable-genome-wide discovery MR hits</b>")) |>
  cols_label(
    Gene = "Gene", Chr = "Chr", `Position (Mb)` = "Position (Mb)",
    `Instruments (N)` = "Instruments (N)", Method = "Method",
    `BBJ Discovery Beta (P)` = "BBJ Discovery β (P)",
    `UKB Replication Beta (P)` = "UKB Replication β (P)",
    `UKB same dir` = "UKB same dir",
    `UKB nominal` = "UKB nominal",
    `FinnGen Sensitivity Beta (P)` = "FinnGen Sensitivity β (P)",
    `FinnGen same dir` = "FinnGen same dir",
    `FinnGen nominal` = "FinnGen nominal",
    Class = "Classification Class"
  ) |>
  tab_options(table.font.size = px(11), heading.align = "left") |>
  tab_source_note(html(
    "<i>BBJ, Biobank Japan; UKB, UK Biobank. "
    %+% "Hits are classified based on location (MHC spillover, chr16p11.2 LD cluster), "
    %+% "replication support, and instrument size. "
    %+% "Same direction indicates whether the outcome effect sign matches that of the discovery analysis. "
    %+% "Nominal significance is defined as having the same direction and P < 0.05.</i>"))

gtsave(gt_s1, file.path(OUTDIR, "TableS1_mr_hits.html"))
tryCatch(gtsave(gt_s1, file.path(OUTDIR, "TableS1_mr_hits.png")),
         error=function(e) cat("Table S1 PNG export failed:", conditionMessage(e), "\n"))

# =============================================================================
# TABLE S5 — chr16p11.2 LD cluster genes
# =============================================================================
# Load S4 instruments to find lead eQTL SNP and position for each gene
inst_raw <- fread(S4_SRC)
inst_raw[, gene_symbol := trimws(gene_symbol)]
top_snps <- inst_raw[, .SD[which.min(pvalue)], by = gene_symbol]
top_snps_map <- top_snps[, .(
  Gene = gene_symbol,
  `Top eQTL SNP` = snp,
  `Top SNP position (Mb)` = round(pos_hg19 / 1e6, 6)
)]

s5_raw_data <- s1_raw[class == "chr16p11_LDcluster"]
s5_data <- s5_raw_data[, .(
  Gene = trimws(gene),
  Chr = chr,
  `Gene position (Mb)` = pos_mb,
  `Top eQTL SNP` = sapply(trimws(gene), function(g) top_snps_map[Gene == g, `Top eQTL SNP`]),
  `Top SNP position (Mb)` = sapply(trimws(gene), function(g) top_snps_map[Gene == g, `Top SNP position (Mb)`]),
  `BBJ Discovery Beta (P)` = paste0(sapply(bbj_beta, fmtBeta), " (", sapply(bbj_p, fmtP), ")"),
  `UKB Replication Beta (P)` = paste0(sapply(ukb_beta, fmtBeta), " (", sapply(ukb_p, fmtP), ")"),
  `FinnGen Sensitivity Beta (P)` = paste0(sapply(fin_beta, fmtBeta), " (", sapply(fin_p, fmtP), ")"),
  Class = "chr16p11.2 LD cluster"
)]

# Clean up any "NA ( - )" formats
for (col in c("BBJ Discovery Beta (P)", "UKB Replication Beta (P)", "FinnGen Sensitivity Beta (P)")) {
  s5_data[[col]] <- gsub("NA \\( - \\)", "-", s5_data[[col]])
}

# Save S5 CSV
fwrite(s5_data, file.path(OUTDIR, "TableS5_chr16_cluster.csv"))

# Build S5 HTML/PNG
gt_s5 <- gt(s5_data) |>
  tab_header(title = html("<b>Table S5. eQTL-mediated Mendelian Randomization findings in the chr16p11.2 LD cluster</b>")) |>
  cols_label(
    Gene = "Gene", Chr = "Chr", `Gene position (Mb)` = "Gene Position (Mb)",
    `Top eQTL SNP` = "Top eQTL SNP", `Top SNP position (Mb)` = "Top SNP Pos (Mb)",
    `BBJ Discovery Beta (P)` = "BBJ Discovery β (P)",
    `UKB Replication Beta (P)` = "UKB Replication β (P)",
    `FinnGen Sensitivity Beta (P)` = "FinnGen Sensitivity β (P)",
    Class = "Classification Class"
  ) |>
  tab_options(table.font.size = px(11), heading.align = "left") |>
  tab_source_note(html(
    "<i>All three genes (HSD3B7, VKORC1, PRSS36) reside within a tight linkage disequilibrium "
    %+% "block at the chr16p11.2 locus. The physical distance between the transcription start sites of "
    %+% "these genes is over 800 kb (spanning from 30.4 Mb to 31.2 Mb), but their lead cis-eQTL instruments "
    %+% "co-localize within a tight 143 kb genomic window (from 31.011 Mb to 31.154 Mb), suggesting limited "
    %+% "ability to distinguish independent gene-specific causal signals using blood cis-eQTL instruments alone.</i>"))

gtsave(gt_s5, file.path(OUTDIR, "TableS5_chr16_cluster.html"))
tryCatch(gtsave(gt_s5, file.path(OUTDIR, "TableS5_chr16_cluster.png")),
         error=function(e) cat("Table S5 PNG export failed:", conditionMessage(e), "\n"))


# =============================================================================
# TABLE S2 — Candidate colocalization results
# =============================================================================
s2_raw <- fread(S2_SRC)

outcome_map <- c(
  BBJ_Graves = "BBJ Graves' Disease",
  UKB_hyperthyroid = "UKB Hyperthyroidism",
  FinnGen_GO = "FinnGen Graves' Ophthalmopathy"
)
verdict_map <- c(
  strong = "Strong colocalization",
  weak = "Weak colocalization",
  distinct = "Distinct signals",
  ambiguous = "Ambiguous"
)

s2_data <- s2_raw[, .(
  Gene = gene_symbol,
  Outcome = outcome_map[outcome],
  `Overlapping SNPs` = n_snps_overlap,
  `PP.H0 (Null)` = sapply(pp_h0, fmtP),
  `PP.H1 (eQTL)` = sapply(pp_h1, fmtP),
  `PP.H2 (GWAS)` = sapply(pp_h2, fmtP),
  `PP.H3 (Distinct)` = sapply(pp_h3, fmtP),
  `PP.H4 (Shared)` = sapply(pp_h4, fmtP),
  `Top SNP` = top_snp,
  Verdict = verdict_map[coloc_verdict]
)]

# Save S2 CSV
fwrite(s2_data, file.path(OUTDIR, "TableS2_coloc.csv"))

# Build S2 HTML/PNG
gt_s2 <- gt(s2_data) |>
  tab_header(title = html("<b>Table S2. Colocalization probability profiles for candidate genes</b>")) |>
  cols_label(
    Gene = "Gene", Outcome = "Outcome", `Overlapping SNPs` = "Overlapping SNPs",
    `PP.H0 (Null)` = "PP.H0 (Null)", `PP.H1 (eQTL)` = "PP.H1 (eQTL)",
    `PP.H2 (GWAS)` = "PP.H2 (GWAS)", `PP.H3 (Distinct)` = "PP.H3 (Distinct)",
    `PP.H4 (Shared)` = "PP.H4 (Shared)", `Top SNP` = "Top SNP", Verdict = "Verdict"
  ) |>
  tab_options(table.font.size = px(11), heading.align = "left") |>
  tab_source_note(html(
    "<i>PP.H4 represents the posterior probability of a shared causal variant between "
    %+% "the eQTL dataset (eQTLGen) and the outcome GWAS. PP.H4 > 0.8 is considered "
    %+% "strong colocalization support. PP.H2-dominant profiles indicate GWAS association "
    %+% "without a shared cis-eQTL signal.</i>"))

gtsave(gt_s2, file.path(OUTDIR, "TableS2_coloc.html"))
tryCatch(gtsave(gt_s2, file.path(OUTDIR, "TableS2_coloc.png")),
         error=function(e) cat("Table S2 PNG export failed:", conditionMessage(e), "\n"))


# =============================================================================
# TABLE S6 — MR Sensitivity / Heterogeneity for Backbone Genes
# =============================================================================
s6_raw <- fread(S7_SRC)

# Filter for backbone genes
s6_sub <- s6_raw[gene_symbol %in% c("TSHR", "IGF1R", "CTLA4")]

# Pretty outcome names
s6_sub[, Outcome := outcome_map[outcome]]

# Select and rename columns
s6_data <- s6_sub[, .(
  Gene = gene_symbol,
  Outcome = Outcome,
  Method = method,
  `Instruments (N)` = n_iv_postharm,
  `Beta (SE)` = sprintf("%+.2f (%.2f)", beta, se),
  `P-value` = sapply(pvalue, fmtP),
  `Cochran Q (P)` = fifelse(is.na(cochran_q), "-", sprintf("%.2f (%s)", cochran_q, sapply(cochran_q_p, fmtP))),
  `Egger Intercept (P)` = fifelse(is.na(egger_intercept), "-", sprintf("%.3f (%s)", egger_intercept, sapply(egger_intercept_p, fmtP)))
)]

# Save S6 CSV
fwrite(s6_data, file.path(OUTDIR, "TableS6_mr_sensitivity.csv"))

# Build S6 HTML/PNG
gt_s6 <- gt(s6_data) |>
  tab_header(title = html("<b>Table S6. Mendelian Randomization sensitivity, heterogeneity, and pleiotropy tests</b>")) |>
  cols_label(
    Gene = "Gene", Outcome = "Outcome", Method = "Method",
    `Instruments (N)` = "Instruments (N)", `Beta (SE)` = "Beta (SE)",
    `P-value` = "P-value", `Cochran Q (P)` = "Cochran's Q (P)",
    `Egger Intercept (P)` = "Egger Intercept (P)"
  ) |>
  tab_options(table.font.size = px(11), heading.align = "left") |>
  tab_source_note(html(
    "<i>Cochran's Q-statistic tests for heterogeneity among instrument effects; "
    %+% "Egger intercept tests for directional pleiotropy. For single-instrument analyses "
    %+% "(Wald ratio), these diagnostic statistics are not calculable (-). "
    %+% "CTLA4 UKB hyperthyroidism Q-statistic shows significant heterogeneity (P = 0.001), "
    %+% "while IGF1R exhibits non-significant Egger intercept and Cochran's Q tests across all outcomes.</i>"))

gtsave(gt_s6, file.path(OUTDIR, "TableS6_mr_sensitivity.html"))
tryCatch(gtsave(gt_s6, file.path(OUTDIR, "TableS6_mr_sensitivity.png")),
         error=function(e) cat("Table S6 PNG export failed:", conditionMessage(e), "\n"))

cat("Supplementary tables generated successfully in:", OUTDIR, "\n")
