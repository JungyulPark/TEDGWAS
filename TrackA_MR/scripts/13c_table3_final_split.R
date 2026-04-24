# =============================================================================
# TED-TRAP — Table 3 FINAL (3-category Main + separate Supplementary)
# NO re-analysis. Splits existing CSV into Main + Supplementary docx.
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(officer)
    library(flextable)
})

# Load current CSV (whatever state it's in — 3-cat or 4-cat)
tbl3_all <- read.csv("c:/ProjectTEDGWAS/TrackA_MR/manuscript/Table3_crossphenotype.csv",
    stringsAsFactors = FALSE
)

cat("=== Current state ===\n")
print(tbl3_all[, c("Outcome", "Category", "Verdict")])

# -----------------------------------------------------------------------------
# Re-map categories (force 3-category scheme for Main + Skin→Supp)
# -----------------------------------------------------------------------------
tbl3_all$Category_new <- sapply(tbl3_all$Outcome, function(x) {
    if (grepl("hearing", x, ignore.case = TRUE)) {
        "Primary off-target (Hearing)"
    } else if (grepl("glucose|hba1c|diabetes", x, ignore.case = TRUE)) {
        "Primary off-target (Metabolic)"
    } else if (grepl("height", x, ignore.case = TRUE)) {
        "Positive control (Biological anchor)"
    } else if (grepl("skin", x, ignore.case = TRUE)) {
        "__TO_SUPP__"
    } else {
        "Other"
    }
})

# Split
main_rows <- tbl3_all[tbl3_all$Category_new != "__TO_SUPP__", ]
supp_rows <- tbl3_all[tbl3_all$Category_new == "__TO_SUPP__", ]

main_rows$Category <- main_rows$Category_new
main_rows$Category_new <- NULL
supp_rows$Category <- "Specificity probe (Pigmentation)"
supp_rows$Category_new <- NULL

# Main order
main_order <- c(
    "Primary off-target (Hearing)",
    "Primary off-target (Metabolic)",
    "Positive control (Biological anchor)"
)
main_rows$Category <- factor(main_rows$Category, levels = main_order)
main_rows <- main_rows[order(main_rows$Category), ]
main_rows$Category <- as.character(main_rows$Category)

cat("\n=== MAIN Table 3 ===\n")
print(main_rows[, c("Outcome", "Category", "N_IV", "IVW_P", "WM_P", "Verdict")])
cat("\n=== SUPPLEMENTARY Skin row ===\n")
print(supp_rows[, c("Outcome", "Category", "N_IV", "IVW_P", "WM_P", "Verdict")])

# -----------------------------------------------------------------------------
# Save CSVs
# -----------------------------------------------------------------------------
write.csv(main_rows, "c:/ProjectTEDGWAS/TrackA_MR/manuscript/Table3_crossphenotype.csv", row.names = FALSE)
write.csv(supp_rows, "c:/ProjectTEDGWAS/TrackA_MR/manuscript/TableS_skincolor_specificity.csv", row.names = FALSE)
cat("\n💾 CSVs saved (Main + Supplementary separately)\n")

# -----------------------------------------------------------------------------
# Regenerate MAIN docx
# -----------------------------------------------------------------------------
doc <- read_docx()
doc <- body_add_par(doc,
    "Table 3. Cross-phenotype Mendelian randomization of IGF1R on teprotumumab-associated adverse event phenotypes and a biological positive control.",
    style = "heading 1"
)
doc <- body_add_par(doc, "")

ft <- flextable(main_rows) %>%
    autofit() %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    bold(part = "header") %>%
    bg(i = which(main_rows$Verdict == "WM only"), bg = "#FFF2CC", part = "body")
doc <- body_add_flextable(doc, ft)

doc <- body_add_par(doc, "")
doc <- body_add_par(doc,
    paste0(
        "Cross-phenotype MR was performed using IGF1R cis-eQTL instruments (eQTLGen ",
        "Consortium, n = 31,684) tested against outcome phenotypes grouped by hypothesis. ",
        "Primary off-target (Metabolic) outcomes — fasting glucose, HbA1c, and type 2 ",
        "diabetes — address the hypothesis that teprotumumab-associated hyperglycemia ",
        "(~15% of treated patients) reflects on-target IGF1R perturbation; all three ",
        "were null under both IVW and weighted median estimators. Primary off-target ",
        "(Hearing) outcomes (self-reported hearing difficulty and hearing aid use, UKB) ",
        "had no instrument overlap with the IGF1R cis-eQTL set in available summary ",
        "statistics and could not be formally tested. Height was included as a ",
        "biological positive-control anchor given the established role of the IGF-1/",
        "growth hormone axis in stature regulation; weighted median P = 0.019 confirms ",
        "that the IGF1R instrument exerts detectable causal effects on a pathway-",
        "appropriate phenotype, arguing against generalized instrument weakness as an ",
        "explanation for the metabolic null findings. A specificity probe against skin ",
        "pigmentation is reported in Supplementary Table Sx. IVW: inverse variance ",
        "weighted. WM: weighted median."
    ),
    style = "Normal"
)
print(doc, target = "c:/ProjectTEDGWAS/TrackA_MR/manuscript/Table3_crossphenotype.docx")
cat("💾 Main Table 3 docx saved\n")

# -----------------------------------------------------------------------------
# Supplementary Skin color docx
# -----------------------------------------------------------------------------
doc_s <- read_docx()
doc_s <- body_add_par(doc_s,
    "Supplementary Table Sx. Specificity probe: IGF1R cross-phenotype MR against skin pigmentation.",
    style = "heading 1"
)
doc_s <- body_add_par(doc_s, "")

ft_s <- flextable(supp_rows) %>%
    autofit() %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    bold(part = "header")
doc_s <- body_add_flextable(doc_s, ft_s)

doc_s <- body_add_par(doc_s, "")
doc_s <- body_add_par(doc_s,
    paste0(
        "As a transparency measure, IGF1R cis-eQTL instruments were also tested against ",
        "skin colour as a specificity probe. A nominal weighted median signal (P = 0.003) ",
        "was observed, consistent with known IGF-1R expression in melanocytes and its ",
        "documented role in pigmentation biology. This finding is reported here for full ",
        "disclosure and does not alter the interpretation of the primary metabolic off-",
        "target analyses in Table 3, as the skin pigmentation axis is mechanistically ",
        "distinct from the glycemic outcomes of interest."
    ),
    style = "Normal"
)
print(doc_s, target = "c:/ProjectTEDGWAS/TrackA_MR/manuscript/TableS_skincolor_specificity.docx")
cat("💾 Supplementary Skin Color docx saved\n")

cat("\n✅ Table 3 FINAL split complete — Main (3-cat) + Supp (Skin).\n")
