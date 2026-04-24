# =============================================================================
# TED-TRAP — Table 2 Draft Generator
# Automatically builds Manuscript Table 2 (MR + Coloc + MVMR combined)
# Output: DOCX-ready CSV + formatted text
# =============================================================================

# Install dependencies if missing
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr", repos = "https://cloud.r-project.org")
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer", repos = "https://cloud.r-project.org")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable", repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# Check for optional flextable / officer (for DOCX)
have_officer <- requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)
if (!have_officer) {
    cat("Note: officer/flextable not installed. Will output CSV + markdown only.\n")
    cat("To enable DOCX export: install.packages(c('officer','flextable'))\n\n")
}

for (d in c("TrackA_MR/manuscript")) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# =============================================================================
# Load all analysis results
# =============================================================================

mr_rescaled <- read.csv("TrackA_MR/results/02b_mr_rescaled.csv", stringsAsFactors = FALSE)
sens <- read.csv("TrackA_MR/results/02_mr_sensitivity.csv", stringsAsFactors = FALSE)
mvmr <- read.csv("TrackA_MR/results/04v3_mvmr_finngen_summary.csv", stringsAsFactors = FALSE)
# Coloc is partly in 03_coloc (BBJ+UKB from Step 3) - reconstruct if needed

# =============================================================================
# PANEL A — Primary MR for Graves disease (BBJ)
# =============================================================================

cat("=== Building Panel A (Primary MR) ===\n")

pa_rows <- mr_rescaled[mr_rescaled$outcome_role == "Primary", ]
pa_focus <- pa_rows[
    pa_rows$gene %in% c("TSHR", "IGF1R"),
    c("gene", "method", "nsnp", "b", "se", "pval")
]

# Calculate OR
pa_focus$OR <- exp(pa_focus$b)
pa_focus$OR_CI_lo <- exp(pa_focus$b - 1.96 * pa_focus$se)
pa_focus$OR_CI_hi <- exp(pa_focus$b + 1.96 * pa_focus$se)

panelA <- data.frame(
    Gene = pa_focus$gene,
    Method = pa_focus$method,
    N_IV = pa_focus$nsnp,
    Beta = sprintf("%+.3f", pa_focus$b),
    SE = sprintf("%.3f", pa_focus$se),
    OR = sprintf("%.2f (%.2f-%.2f)", pa_focus$OR, pa_focus$OR_CI_lo, pa_focus$OR_CI_hi),
    P_value = ifelse(pa_focus$pval < 0.001,
        formatC(pa_focus$pval, format = "e", digits = 2),
        sprintf("%.3f", pa_focus$pval)
    ),
    Sig = ifelse(pa_focus$pval < 2.08e-3, "***",
        ifelse(pa_focus$pval < 0.05, "*", "")
    ),
    stringsAsFactors = FALSE
)

cat("\nPanel A preview:\n")
print(panelA, row.names = FALSE)

# =============================================================================
# PANEL B — Replication in UKB (rescaled to log-odds)
# =============================================================================

cat("\n=== Building Panel B (Replication) ===\n")

pb_rows <- mr_rescaled[mr_rescaled$outcome_role == "Replication" &
    mr_rescaled$gene %in% c("TSHR", "IGF1R") &
    mr_rescaled$method %in% c("Inverse variance weighted", "Wald ratio"), ]

panelB <- data.frame(
    Gene = pb_rows$gene,
    Method = pb_rows$method,
    N_IV = pb_rows$nsnp,
    Beta_rescaled = sprintf("%+.3f", pb_rows$b_rescaled),
    SE_rescaled = sprintf("%.3f", pb_rows$se_rescaled),
    P_value = ifelse(pb_rows$pval < 0.001,
        formatC(pb_rows$pval, format = "e", digits = 2),
        sprintf("%.3f", pb_rows$pval)
    ),
    Sig = ifelse(pb_rows$pval < 2.08e-3, "***",
        ifelse(pb_rows$pval < 0.05, "*", "")
    ),
    stringsAsFactors = FALSE
)

cat("\nPanel B preview:\n")
print(panelB, row.names = FALSE)

# =============================================================================
# PANEL C — MVMR (TED-specific, FinnGen outcome)
# =============================================================================

cat("\n=== Building Panel C (MVMR) ===\n")

panelC <- data.frame(
    Exposure = c(
        "TSHR (adj. GD liability)", "GD liability (adj. TSHR)",
        "TSHR — Univariable reference"
    ),
    Method = c("MVMR-IVW", "MVMR-IVW", "Univariable IVW"),
    N_IV = c(mvmr$N_SNPs[1], mvmr$N_SNPs[2], mvmr$N_SNPs[3]),
    Beta = sprintf("%+.3f", mvmr$Beta),
    SE = sprintf("%.3f", mvmr$SE),
    Conditional_F = ifelse(is.na(mvmr$Cond_F), "—", sprintf("%.1f", mvmr$Cond_F)),
    P_value = ifelse(mvmr$P < 0.001,
        formatC(mvmr$P, format = "e", digits = 2),
        sprintf("%.3f", mvmr$P)
    ),
    Sig = ifelse(mvmr$P < 0.001, "***",
        ifelse(mvmr$P < 0.05, "*", "")
    ),
    stringsAsFactors = FALSE
)

cat("\nPanel C preview:\n")
print(panelC, row.names = FALSE)

# =============================================================================
# PANEL D — Colocalization at TSHR locus
# =============================================================================

cat("\n=== Building Panel D (Coloc) ===\n")

# Hard-code from Step 3 + 3b results (reader can reconstruct from logs if needed)
panelD <- data.frame(
    Outcome = c(
        "Graves disease (BBJ)",
        "Hyperthyroidism (UKB)",
        "TED/GO (FinnGen R12)"
    ),
    Cohort_N = c(
        "212,453 (2,809 cases)",
        "484,598 (3,731 cases)",
        "500,348 (858 cases)"
    ),
    Common_SNPs = c("2,707", "3,376", "Technical limitation"),
    PP_H3 = c("0.049", "0.774", "NA"),
    PP_H4 = c("0.951 ⭐⭐⭐", "0.227", "NA"),
    Top_SNP = c("rs179252", "rs1023586", "—"),
    Note = c(
        "Shared causal variant",
        sprintf("Different lead tag; LD r²=0.85 in EAS with rs179252"),
        "Technical limitation (attempted)"
    ),
    stringsAsFactors = FALSE
)

cat("\nPanel D preview:\n")
print(panelD, row.names = FALSE)

# =============================================================================
# SAVE — CSV outputs
# =============================================================================

write.csv(panelA, "TrackA_MR/manuscript/Table2_PanelA_PrimaryMR.csv", row.names = FALSE)
write.csv(panelB, "TrackA_MR/manuscript/Table2_PanelB_Replication.csv", row.names = FALSE)
write.csv(panelC, "TrackA_MR/manuscript/Table2_PanelC_MVMR.csv", row.names = FALSE)
write.csv(panelD, "TrackA_MR/manuscript/Table2_PanelD_Coloc.csv", row.names = FALSE)

cat("\n💾 Saved CSVs to TrackA_MR/manuscript/\n")

# =============================================================================
# MARKDOWN version (for direct paste into Word/Google Docs)
# =============================================================================

md_text <- c(
    "",
    "## Table 2. Mendelian randomization and colocalization evidence for TSHR versus IGF1R causal effects on Graves disease and thyroid eye disease.",
    "",
    "### Panel A — Primary cis-MR (Graves disease, Biobank Japan, East Asian)",
    "",
    knitr::kable(panelA, format = "markdown"),
    "",
    "### Panel B — Replication (UK Biobank, European; β rescaled to log-odds using Lloyd-Jones 2018 method)",
    "",
    knitr::kable(panelB, format = "markdown"),
    "",
    "### Panel C — Multivariable MR for TED-specific signal (outcome: FinnGen R12 Graves ophthalmopathy)",
    "",
    knitr::kable(panelC, format = "markdown"),
    "",
    "### Panel D — Bayesian colocalization at TSHR locus (chr14:80.9-82.1 Mb)",
    "",
    knitr::kable(panelD, format = "markdown"),
    "",
    "### Footnotes",
    "",
    "*** P < 2.08×10⁻³ (Bonferroni-adjusted for 8 genes × 3 outcomes = 24 tests).",
    "* P < 0.05 (nominal).",
    "β: log-odds effect per 1-SD increase in gene expression.",
    "OR (95% CI): odds ratio computed as exp(β).",
    "N_IV: number of independent instruments (r² < 0.001, clumping window 10 Mb).",
    "Conditional F: Sanderson 2019 conditional F-statistic for instrument strength (>10 = strong).",
    "PP.H4: posterior probability of shared causal variant (Giambartolomei 2014; default priors).",
    "Panel B β values rescaled from UKB linear LMM to log-odds scale using factor 1/[p(1-p)] where p = case prevalence.",
    "",
    "Primary outcome: Graves disease (Sakaue et al., 2021; GCST90018627).",
    "Replication outcome: Hyperthyroidism (Dönertaş et al., 2021; GCST90038636).",
    "Sensitivity/TED outcome: FinnGen Release 12, Graves ophthalmopathy endpoint (858 cases)."
)

writeLines(md_text, "TrackA_MR/manuscript/Table2_markdown.md")
cat("💾 Saved: TrackA_MR/manuscript/Table2_markdown.md\n")

# =============================================================================
# DOCX (optional)
# =============================================================================

if (have_officer) {
    library(officer)
    library(flextable)

    doc <- read_docx()
    doc <- body_add_par(doc, "Table 2. Mendelian randomization and colocalization evidence",
        style = "heading 1"
    )

    doc <- body_add_par(doc, "Panel A — Primary cis-MR (Graves disease, BBJ)", style = "heading 2")
    ft_a <- flextable(panelA) %>%
        autofit() %>%
        theme_booktabs()
    doc <- body_add_flextable(doc, ft_a)
    doc <- body_add_par(doc, "", style = "Normal")

    doc <- body_add_par(doc, "Panel B — Replication (UKB, rescaled)", style = "heading 2")
    ft_b <- flextable(panelB) %>%
        autofit() %>%
        theme_booktabs()
    doc <- body_add_flextable(doc, ft_b)
    doc <- body_add_par(doc, "", style = "Normal")

    doc <- body_add_par(doc, "Panel C — MVMR for TED-specific signal", style = "heading 2")
    ft_c <- flextable(panelC) %>%
        autofit() %>%
        theme_booktabs()
    doc <- body_add_flextable(doc, ft_c)
    doc <- body_add_par(doc, "", style = "Normal")

    doc <- body_add_par(doc, "Panel D — Colocalization at TSHR locus", style = "heading 2")
    ft_d <- flextable(panelD) %>%
        autofit() %>%
        theme_booktabs()
    doc <- body_add_flextable(doc, ft_d)

    print(doc, target = "TrackA_MR/manuscript/Table2_draft.docx")
    cat("\n💾 Saved: TrackA_MR/manuscript/Table2_draft.docx\n")
} else {
    cat("\n[INFO] DOCX export skipped. To enable:\n")
    cat("  install.packages(c('officer','flextable','knitr'))\n")
}

cat("\n✅ Table 2 draft complete.\n")
cat("   4 Panels saved as CSV + Markdown (+ DOCX if available)\n\n")
