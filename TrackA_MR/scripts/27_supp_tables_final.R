# =============================================================================
# TED-TRAP — Supp Table S1 (Full MR) + S2 (Skin) final compilation
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(officer)
    library(flextable)
})

setwd("c:/ProjectTEDGWAS")
if (!dir.exists("TrackA_MR/manuscript")) dir.create("TrackA_MR/manuscript", recursive = TRUE)

# =============================================================================
# STEP 1: Explore all MR-related result files
# =============================================================================
cat("=== MR-related result files available ===\n")
res_files <- list.files("TrackA_MR/results",
    pattern = "\\.csv$",
    full.names = TRUE, recursive = TRUE
)
cat("All CSV files:\n")
print(basename(res_files))

# Identify full-MR candidates
mr_candidates <- res_files[grepl("mr_main|mr_full|mr_sens|02_|mr_primary|mr_rep|mr_repl",
    res_files,
    ignore.case = TRUE
)]
cat("\nLikely full MR files:\n")
print(basename(mr_candidates))

# =============================================================================
# STEP 2: Build Supp Table S1 — Full MR
# =============================================================================
# Load primary
if (file.exists("TrackA_MR/results/02_mr_main.csv")) {
    mr_main <- read.csv("TrackA_MR/results/02_mr_main.csv", stringsAsFactors = FALSE)
    cat("\n=== mr_main columns ===\n")
    print(names(mr_main))
    cat("=== Outcome roles in mr_main ===\n")
    print(unique(mr_main$outcome_role))

    # Format full MR table
    full_mr <- mr_main %>%
        mutate(
            Beta_fmt = sprintf("%+.3f", b),
            SE_fmt = sprintf("%.3f", se),
            OR_CI = sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)),
            P_fmt = ifelse(pval < 0.001, formatC(pval, format = "e", digits = 2),
                sprintf("%.3f", pval)
            ),
            Sig = case_when(
                pval < 2.08e-3 ~ "***", # Bonferroni for 8 genes × 3 outcomes
                pval < 0.05 ~ "*",
                TRUE ~ ""
            )
        ) %>%
        select(
            Gene = gene, Outcome_role = outcome_role, Method = method,
            N_IV = nsnp, Beta = Beta_fmt, SE = SE_fmt, OR_95CI = OR_CI,
            P_value = P_fmt, Significance = Sig
        )

    # Sort: gene → outcome_role → method (IVW/Wald first)
    method_order <- c(
        "Inverse variance weighted", "Wald ratio",
        "MR Egger", "Weighted median", "Weighted mode", "Simple mode"
    )
    full_mr$Method <- factor(full_mr$Method, levels = method_order)
    full_mr <- full_mr %>%
        arrange(Gene, Outcome_role, Method) %>%
        mutate(Method = as.character(Method))

    cat("\n=== Full MR table preview ===\n")
    print(head(full_mr, 15))
    cat("Total rows:", nrow(full_mr), "\n")

    write.csv(full_mr, "TrackA_MR/manuscript/TableS1_full_MR.csv", row.names = FALSE)

    # -----------------------------------------------------------------------------
    # Supp Table S1 DOCX
    # -----------------------------------------------------------------------------
    doc <- read_docx()
    doc <- body_add_par(doc,
        "Supplementary Table S1. Full Mendelian randomization results across all pre-specified candidate genes, outcomes, and methods.",
        style = "heading 1"
    )
    doc <- body_add_par(doc, "")

    ft <- flextable(full_mr) %>%
        autofit() %>%
        theme_booktabs() %>%
        fontsize(size = 8, part = "all") %>%
        bold(part = "header") %>%
        # Highlight Bonferroni-significant rows
        bg(i = which(full_mr$Significance == "***"), bg = "#E2EFDA", part = "body") %>%
        # Highlight nominal-significant rows
        bg(i = which(full_mr$Significance == "*"), bg = "#FFF2CC", part = "body")

    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "")

    # S1 footnote
    doc <- body_add_par(doc, paste0(
        "Full two-sample Mendelian randomization results for 8 genes (TSHR, IGF1R, TNF, PPARG, ARRB1, IRS1, AKT1, CTLA4) ",
        "for which at least one cis-eQTL instrument could be extracted from the eQTLGen Consortium (n = 31,684). ",
        "Outcome roles: Primary = Graves disease, Biobank Japan (GCST90018627, East Asian, 2,809 cases); ",
        "Replication = Hyperthyroidism, UK Biobank (GCST90038636, European, 3,731 cases); ",
        "Sensitivity = Graves ophthalmopathy, FinnGen Release 12 (858 cases). ",
        "Methods include inverse variance weighted (IVW) and Wald ratio as primary estimators, ",
        "with MR-Egger, weighted median, weighted mode, and simple mode as sensitivity analyses where instrument count permitted. ",
        "Effect sizes (β) are in log-odds per 1-SD increase in gene expression for Primary and Sensitivity outcomes, ",
        "and in linear mixed model units for the UK Biobank Replication outcome; ",
        "the Replication β has been rescaled to a comparable log-odds scale in Table 2 Panel B. ",
        "Significance thresholds: *** P < 2.08×10⁻³ (Bonferroni for 24 gene-outcome tests); * P < 0.05 (nominal). ",
        "Genes with no instruments available (IGF1, IL6, TGFB1, HAS2, PPARG, INSR, IRS2, FOXO1, PIK3R1, PDPK1, ",
        "HAS1, HAS3, ADIPOQ, FABP4, CEBPA) are not listed here; their MR non-testing is documented in Supplementary Figure S1."
    ), style = "Normal")

    print(doc, target = "TrackA_MR/manuscript/TableS1_full_MR.docx")
    cat("\n✅ Supp Table S1 saved\n")
} else {
    cat("\n❌ ERROR: 02_mr_main.csv not found!\n")
}

# =============================================================================
# STEP 3: Supp Table S2 — Skin color (relabel existing file)
# =============================================================================
# Existing file: TableS_skincolor_specificity.docx → rename to TableS2
existing_skin <- list.files("TrackA_MR/manuscript",
    pattern = "skin|Skin|specificity",
    full.names = TRUE, ignore.case = TRUE
)
cat("\n=== Existing Skin files ===\n")
print(basename(existing_skin))

# Load and rebuild as Supp Table S2
skin_csv <- existing_skin[grepl("\\.csv$", existing_skin)]
if (length(skin_csv) > 0) {
    skin_data <- read.csv(skin_csv[1], stringsAsFactors = FALSE)
    cat("Skin data:\n")
    print(skin_data)

    doc_s2 <- read_docx()
    doc_s2 <- body_add_par(doc_s2,
        "Supplementary Table S2. Specificity probe: IGF1R cis-eQTL Mendelian randomization against skin pigmentation.",
        style = "heading 1"
    )
    doc_s2 <- body_add_par(doc_s2, "")

    ft_s2 <- flextable(skin_data) %>%
        autofit() %>%
        theme_booktabs() %>%
        fontsize(size = 9, part = "all") %>%
        bold(part = "header")
    doc_s2 <- body_add_flextable(doc_s2, ft_s2)

    doc_s2 <- body_add_par(doc_s2, "")
    doc_s2 <- body_add_par(doc_s2, paste0(
        "As a specificity probe beyond the metabolic off-target phenotypes reported in Table 3, ",
        "IGF1R cis-eQTL instruments were tested against skin pigmentation (UK Biobank). ",
        "A nominal weighted median association was observed (β = +0.007, P_WM = 0.003), ",
        "which is biologically plausible given the documented expression of IGF-1R on melanocytes ",
        "and its role in pigmentation regulation. This finding is reported for full disclosure ",
        "and does not alter the interpretation of the primary metabolic off-target analyses, ",
        "as the skin pigmentation axis is mechanistically distinct from glycemic phenotypes. ",
        "IVW: inverse variance weighted; WM: weighted median."
    ), style = "Normal")

    # Save as standardized S2 filename
    print(doc_s2, target = "TrackA_MR/manuscript/TableS2_skincolor_specificity.docx")
    write.csv(skin_data, "TrackA_MR/manuscript/TableS2_skincolor_specificity.csv",
        row.names = FALSE
    )
    cat("\n✅ Supp Table S2 saved (Skin color specificity probe)\n")
} else {
    cat("\n❌ ERROR: Skin color CSV not found!\n")
}

# =============================================================================
# STEP 4: Update Main Table 3 footnote to reference Supp Table S2
# =============================================================================
# Reload Main Table 3 CSV
if (file.exists("TrackA_MR/manuscript/Table3_crossphenotype.csv")) {
    t3 <- read.csv("TrackA_MR/manuscript/Table3_crossphenotype.csv", stringsAsFactors = FALSE)

    doc_t3 <- read_docx()
    doc_t3 <- body_add_par(doc_t3,
        "Table 3. Cross-phenotype Mendelian randomization of IGF1R on teprotumumab-associated adverse event phenotypes and a biological positive control.",
        style = "heading 1"
    )
    doc_t3 <- body_add_par(doc_t3, "")

    ft_t3 <- flextable(t3) %>%
        autofit() %>%
        theme_booktabs() %>%
        fontsize(size = 9, part = "all") %>%
        bold(part = "header") %>%
        bg(i = which(t3$Verdict == "WM only"), bg = "#FFF2CC", part = "body")
    doc_t3 <- body_add_flextable(doc_t3, ft_t3)

    doc_t3 <- body_add_par(doc_t3, "")
    doc_t3 <- body_add_par(doc_t3, paste0(
        "Cross-phenotype MR was performed using IGF1R cis-eQTL instruments (eQTLGen Consortium, n = 31,684) ",
        "tested against outcome phenotypes grouped by hypothesis. Primary off-target (Metabolic) outcomes — ",
        "fasting glucose, HbA1c, and type 2 diabetes — address the hypothesis that teprotumumab-associated ",
        "hyperglycemia (~15% of treated patients) reflects on-target IGF1R perturbation; all three were null ",
        "under both IVW and weighted median estimators. Primary off-target (Hearing) outcomes ",
        "(self-reported hearing difficulty and hearing aid use, UKB) had no instrument overlap with the ",
        "IGF1R cis-eQTL set in available summary statistics and could not be formally tested. Height was ",
        "included as a biological positive-control anchor given the established role of the IGF-1/growth ",
        "hormone axis in stature regulation; weighted median P = 0.019 confirms that the IGF1R instrument ",
        "exerts detectable causal effects on a pathway-appropriate phenotype, arguing against generalized ",
        "instrument weakness as an explanation for the metabolic null findings. ",
        "An additional specificity probe against skin pigmentation is reported in Supplementary Table S2. ", # PATCH
        "IVW: inverse variance weighted. WM: weighted median."
    ), style = "Normal")

    print(doc_t3, target = "TrackA_MR/manuscript/Table3_crossphenotype.docx")
    cat("\n✅ Main Table 3 updated with Supp Table S2 reference\n")
} else {
    cat("\n❌ ERROR: Table3_crossphenotype.csv not found!\n")
}

cat("\n\n========== ALL SUPP TABLES FINALIZED ==========\n")
cat("  ✅ Supp Table S1: Full MR results (8 genes × 3 outcomes × all methods)\n")
cat("  ✅ Supp Table S2: Skin color specificity probe\n")
cat("  ✅ Main Table 3: Footnote updated with S2 reference\n")
