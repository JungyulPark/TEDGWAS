#!/usr/bin/env Rscript
# =============================================================================
# TED-TRAP v5 — Supplementary Table S3: TSHR SuSiE Credible Sets
# =============================================================================
# Output: CSV, HTML, and PNG formats for Table S3.
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

# Source data path
S3_SRC <- file.path(WORK, "02_tshr_finemap/results/TaskB_04_tshr_susie_credible_sets_v1.tsv")

# Load credible sets
s3_raw <- fread(S3_SRC)

# Format P-value helper
fmtP <- function(p) {
  if (is.na(p) || p == "") return("-")
  p_num <- as.numeric(p)
  if (is.na(p_num)) return(p)
  if (p_num < 1e-3) sprintf("%.2e", p_num) else sprintf("%.3f", p_num)
}

# Process data
s3_data <- s3_raw[, .(
  `LD reference panel` = fifelse(grepl("EUR", ld_reference_id), "European (1000G)", "East Asian (1000G)"),
  `Credible Set ID` = paste0("CS_", cs_id),
  `SNPs in CS` = n_snps_in_cs,
  Purity = round(cs_purity, 3),
  `Top SNP` = top_snp,
  `Top SNP PIP` = round(top_snp_pip, 3),
  `Position (hg19)` = top_snp_pos_hg19,
  `eQTL P-value` = sapply(top_snp_p, fmtP),
  `LD (r²) with rs179252 (EUR)` = round(ld_top_with_rs179252_eur, 3),
  `LD (r²) with rs179252 (EAS)` = round(ld_top_with_rs179252_eas, 3)
)]

# Save Table S3 CSV
fwrite(s3_data, file.path(OUTDIR, "TableS3_finemapping.csv"))

# Build S3 HTML/PNG
gt_s3 <- gt(s3_data) |>
  tab_header(title = html("<b>Table S3. TSHR cis-eQTL SuSiE fine-mapping credible sets</b>")) |>
  cols_label(
    `LD reference panel` = "LD reference panel",
    `Credible Set ID` = "CS ID",
    `SNPs in CS` = "N SNPs",
    Purity = "Purity",
    `Top SNP` = "Top SNP",
    `Top SNP PIP` = "Top SNP PIP",
    `Position (hg19)` = "Pos (hg19)",
    `eQTL P-value` = "eQTL P-value",
    `LD (r²) with rs179252 (EUR)` = "LD r² (EUR)",
    `LD (r²) with rs179252 (EAS)` = "LD r² (EAS)"
  ) |>
  tab_options(table.font.size = px(11), heading.align = "left") |>
  tab_source_note(html(
    "<i>TSHR cis-eQTL fine-mapping was performed using SuSiE with 1000 Genomes Phase 3 European "
    %+% "and East Asian LD reference panels. Credible sets resolved to single lead SNPs with PIP=1.0. "
    %+% "Secondary credible-set lead variants showed substantial LD with rs179252 in the corresponding "
    %+% "ancestry-matched LD reference panel, with minimum ancestry-matched r²=0.808 in EAS and r²=0.965 "
    %+% "in EUR, supporting retention of TSHR as a single-IV/locus-level anchor rather than a multi-IV exposure.</i>"
  ))

gtsave(gt_s3, file.path(OUTDIR, "TableS3_finemapping.html"))
tryCatch(gtsave(gt_s3, file.path(OUTDIR, "TableS3_finemapping.png")),
         error=function(e) cat("Table S3 PNG export failed:", conditionMessage(e), "\n"))

cat("Table S3 generated successfully in:", OUTDIR, "\n")
