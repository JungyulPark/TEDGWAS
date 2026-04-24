# =============================================================================
# TED-TRAP — Table 3 라벨 재매핑 — Skin color 처리 추가
# (Height = biological anchor, Skin color = specificity probe 둘 다 분리)
# =============================================================================

# CSV 원본 load
tbl3_final <- read.csv("c:/ProjectTEDGWAS/TrackA_MR/manuscript/Table3_crossphenotype.csv",
    stringsAsFactors = FALSE
)

# Gene-specific mapping — outcome명으로 직접 분류
tbl3_final$Category <- sapply(tbl3_final$Outcome, function(x) {
    if (grepl("hearing|Hearing", x, ignore.case = TRUE)) {
        "Primary off-target (Hearing)"
    } else if (grepl("Fasting_Glucose|HbA1c|Type2_Diabetes|Glucose|HbA1C|diabetes",
        x,
        ignore.case = TRUE
    )) {
        "Primary off-target (Metabolic)"
    } else if (grepl("Height|height", x, ignore.case = TRUE)) {
        "Positive control (Biological anchor)"
    } else if (grepl("Skin|skin", x, ignore.case = TRUE)) {
        "Specificity probe (Pigmentation)"
    } else {
        "Other"
    }
})

# 출력 순서 재정렬
cat_order_new <- c(
    "Primary off-target (Hearing)",
    "Primary off-target (Metabolic)",
    "Positive control (Biological anchor)",
    "Specificity probe (Pigmentation)"
)
tbl3_final$Category <- factor(tbl3_final$Category, levels = cat_order_new)
tbl3_final <- tbl3_final[order(tbl3_final$Category), ]

cat("=== Table 3 re-labeled ===\n")
print(tbl3_final, row.names = FALSE)

# Re-save CSV
write.csv(tbl3_final, "c:/ProjectTEDGWAS/TrackA_MR/manuscript/Table3_crossphenotype.csv",
    row.names = FALSE
)

# Re-generate DOCX with improved footnote
if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)) {
    library(officer)
    library(flextable)

    doc <- read_docx()
    doc <- body_add_par(doc,
        "Table 3. Cross-phenotype Mendelian randomization of IGF1R genetic effects on teprotumumab-associated adverse event phenotypes.",
        style = "heading 1"
    )
    doc <- body_add_par(doc, "")

    ft <- flextable(tbl3_final) %>%
        autofit() %>%
        theme_booktabs() %>%
        fontsize(size = 9, part = "all")
    doc <- body_add_flextable(doc, ft)

    doc <- body_add_par(doc, "")
    doc <- body_add_par(doc,
        paste0(
            "IVW: Inverse variance weighted. WM: Weighted median. Egger int. P: MR-Egger intercept P-value ",
            "for directional pleiotropy. Cochran Q P: heterogeneity test. ",
            "Hearing phenotypes (ukb-b-12002, ukb-b-14135) had no instrument overlap with IGF1R IVs in ",
            "OpenGWAS summary statistics and could not be formally tested. ",
            "Height is included as a biological positive-control anchor reflecting the established role of ",
            "the IGF-1 axis in growth regulation; a nominal weighted-median association is expected and ",
            "indicates functional instrument validity. Skin color is included as a specificity probe to ",
            "assess potential pleiotropic effects on pigmentation biology (melanocytes express IGF-1R); ",
            "the nominal signal is interpreted as exploratory and disclosed transparently."
        ),
        style = "Normal"
    )

    print(doc, target = "c:/ProjectTEDGWAS/TrackA_MR/manuscript/Table3_crossphenotype.docx")
    cat("\n💾 Saved revised: Table3_crossphenotype.docx\n")
}

cat("\n✅ Table 3 re-labeling complete.\n")
