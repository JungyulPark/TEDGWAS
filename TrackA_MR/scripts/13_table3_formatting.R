# =============================================================================
# TED-TRAP — Table 3: Cross-phenotype MR (SE2) — Publication Formatting
# Uses already-completed Step 5 results (05_cross_phenotype_mr.csv)
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
})

have_officer <- requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)

if (!dir.exists("TrackA_MR/manuscript")) dir.create("TrackA_MR/manuscript", recursive = TRUE)

# Load existing Step 5 results (DO NOT re-run analysis)
mr <- read.csv("TrackA_MR/results/05_cross_phenotype_mr.csv", stringsAsFactors = FALSE)
sens <- read.csv("TrackA_MR/results/05_cross_phenotype_sensitivity.csv", stringsAsFactors = FALSE)

# =============================================================================
# Build Table 3 — IVW primary + WM + sensitivity summary
# =============================================================================

ivw <- mr[mr$method == "Inverse variance weighted", ]
wm <- mr[mr$method == "Weighted median", c("outcome_name", "b", "pval")]
colnames(wm) <- c("outcome_name", "WM_b", "WM_pval")

eg <- mr[mr$method == "MR Egger", c("outcome_name", "b", "pval")]
colnames(eg) <- c("outcome_name", "Egger_b", "Egger_pval")

tbl3 <- ivw %>%
    select(outcome_name, outcome_role, nsnp, b, se, pval) %>%
    left_join(wm, by = "outcome_name") %>%
    left_join(eg, by = "outcome_name") %>%
    left_join(sens %>% select(outcome_name, Egger_int_p, Cochran_Q_p),
        by = "outcome_name"
    )

tbl3_fmt <- data.frame(
    Outcome = tbl3$outcome_name,
    Category = tbl3$outcome_role,
    N_IV = tbl3$nsnp,
    IVW_beta = sprintf("%+.3f", tbl3$b),
    IVW_SE = sprintf("%.3f", tbl3$se),
    IVW_P = ifelse(tbl3$pval < 0.001,
        formatC(tbl3$pval, format = "e", digits = 2),
        sprintf("%.3f", tbl3$pval)
    ),
    WM_beta = sprintf("%+.3f", tbl3$WM_b),
    WM_P = ifelse(tbl3$WM_pval < 0.001,
        formatC(tbl3$WM_pval, format = "e", digits = 2),
        sprintf("%.3f", tbl3$WM_pval)
    ),
    Egger_int_P = sprintf("%.3f", tbl3$Egger_int_p),
    Cochran_Q_P = ifelse(tbl3$Cochran_Q_p < 0.001,
        formatC(tbl3$Cochran_Q_p, format = "e", digits = 2),
        sprintf("%.3f", tbl3$Cochran_Q_p)
    ),
    Verdict = sapply(seq_len(nrow(tbl3)), function(i) {
        if (tbl3$pval[i] < 0.05) {
            return("Signal")
        }
        if (!is.na(tbl3$WM_pval[i]) && tbl3$WM_pval[i] < 0.05) {
            return("WM only")
        }
        return("Null")
    }),
    stringsAsFactors = FALSE
)

# Add Hearing rows (IV overlap = 0)
hearing_rows <- data.frame(
    Outcome = c(
        "Self-reported hearing difficulty (UKB)",
        "Hearing aid use (UKB)"
    ),
    Category = "primary_offtarget",
    N_IV = 0,
    IVW_beta = "—", IVW_SE = "—", IVW_P = "—",
    WM_beta = "—", WM_P = "—",
    Egger_int_P = "—", Cochran_Q_P = "—",
    Verdict = "N/A"
)

tbl3_final <- rbind(tbl3_fmt, hearing_rows)

# Order by category
cat_order <- c("primary_offtarget", "secondary_offtarget", "negative_control")
tbl3_final$Category <- factor(tbl3_final$Category, levels = cat_order)
tbl3_final <- tbl3_final[order(tbl3_final$Category), ]

cat("=== Table 3: Cross-phenotype MR ===\n")
print(tbl3_final, row.names = FALSE)

write.csv(tbl3_final, "TrackA_MR/manuscript/Table3_crossphenotype.csv", row.names = FALSE)

# DOCX export
if (have_officer) {
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
        "IVW: Inverse variance weighted. WM: Weighted median. Egger int. P: MR-Egger intercept P-value for directional pleiotropy. Cochran Q P: heterogeneity test. Hearing phenotypes (ukb-b-12002, ukb-b-14135) had no instrument overlap with IGF1R IVs in OpenGWAS summary statistics and could not be formally tested. Height serves as a biological positive control given the established role of IGF axis in growth.",
        style = "Normal"
    )

    print(doc, target = "TrackA_MR/manuscript/Table3_crossphenotype.docx")
    cat("\n💾 Saved: Table3_crossphenotype.docx\n")
} else {
    cat("\n⚠️  'officer' or 'flextable' not installed. Could not generate .docx, skipping docx export.\n")
}

cat("\n✅ Table 3 formatting complete (no re-analysis, uses Step 5 results).\n")
