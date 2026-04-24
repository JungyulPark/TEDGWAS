# =============================================================================
# TED-TRAP Step 2b — Lloyd-Jones β Rescale + Consolidated Summary
# Purpose:
#   1. Rescale Replication (UKB linear LMM) β to log-odds scale
#   2. Confirm cross-ancestry concordance
#   3. Build final Primary Endpoint table
# =============================================================================

library(dplyr)

mr_main <- read.csv("TrackA_MR/results/02_mr_main.csv")
sens <- read.csv("TrackA_MR/results/02_mr_sensitivity.csv")

cat("=== Step 2b: Scale harmonization ===\n\n")

# --- Lloyd-Jones 2018 rescaling for Replication (UKB linear LMM) ---
# β_logOR ≈ β_linear / [p(1-p)]
# p = case prevalence in UKB
p_rep <- 3731 / (3731 + 480867)
rescale_factor_rep <- 1 / (p_rep * (1 - p_rep))
cat(sprintf("Replication case prevalence: p = %.5f\n", p_rep))
cat(sprintf("Rescale factor: 1 / [p(1-p)] = %.1f\n\n", rescale_factor_rep))

# Apply rescale only to Replication rows
mr_rescaled <- mr_main
rep_mask <- mr_rescaled$outcome_role == "Replication"
mr_rescaled$b_rescaled <- ifelse(rep_mask,
    mr_rescaled$b * rescale_factor_rep,
    mr_rescaled$b
)
mr_rescaled$se_rescaled <- ifelse(rep_mask,
    mr_rescaled$se * rescale_factor_rep,
    mr_rescaled$se
)
mr_rescaled$scale_note <- ifelse(rep_mask,
    sprintf("rescaled (×%.1f)", rescale_factor_rep),
    "as reported (log-odds)"
)

# --- Verify rescaling worked ---
cat("=== TSHR Cross-Ancestry Concordance After Rescale ===\n")
tshr_check <- mr_rescaled[
    mr_rescaled$gene == "TSHR",
    c(
        "outcome_role", "method", "nsnp",
        "b", "b_rescaled", "se_rescaled", "pval", "scale_note"
    )
]
print(tshr_check, row.names = FALSE)

cat("\n=== IGF1R Cross-Ancestry Concordance After Rescale ===\n")
igf_check <- mr_rescaled[
    mr_rescaled$gene == "IGF1R" &
        mr_rescaled$method == "Inverse variance weighted",
    c(
        "outcome_role", "method", "nsnp",
        "b", "b_rescaled", "se_rescaled", "pval", "scale_note"
    )
]
print(igf_check, row.names = FALSE)

# --- Save rescaled results ---
write.csv(mr_rescaled, "TrackA_MR/results/02b_mr_rescaled.csv", row.names = FALSE)
cat("\n💾 Saved: TrackA_MR/results/02b_mr_rescaled.csv\n")

# =============================================================================
# PRIMARY ENDPOINT — Final consolidated table (Graves disease)
# =============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("FINAL PRIMARY ENDPOINT TABLE (for Table 2 Panel A)\n")
cat(rep("=", 70), "\n\n", sep = "")

primary <- mr_rescaled[
    mr_rescaled$outcome_role == "Primary",
    c("gene", "method", "nsnp", "b", "se", "pval")
]
primary$OR <- exp(primary$b)
primary$OR_CI_lo <- exp(primary$b - 1.96 * primary$se)
primary$OR_CI_hi <- exp(primary$b + 1.96 * primary$se)
primary$signif <- ifelse(primary$pval < 2.08e-3, "***",
    ifelse(primary$pval < 0.05, "*", "")
)

primary_display <- primary
primary_display$b <- round(primary_display$b, 3)
primary_display$se <- round(primary_display$se, 3)
primary_display$pval <- signif(primary_display$pval, 3)
primary_display$OR <- round(primary_display$OR, 2)
primary_display$OR_CI_lo <- round(primary_display$OR_CI_lo, 2)
primary_display$OR_CI_hi <- round(primary_display$OR_CI_hi, 2)

print(primary_display, row.names = FALSE)

# =============================================================================
# REPLICATION ENDPOINT — with rescaled β for direct comparison
# =============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("REPLICATION ENDPOINT (UKB rescaled to log-odds)\n")
cat(rep("=", 70), "\n\n", sep = "")

rep_tab <- mr_rescaled[
    mr_rescaled$outcome_role == "Replication",
    c(
        "gene", "method", "nsnp",
        "b_rescaled", "se_rescaled", "pval"
    )
]
rep_tab$b_rescaled <- round(rep_tab$b_rescaled, 3)
rep_tab$se_rescaled <- round(rep_tab$se_rescaled, 3)
rep_tab$pval <- signif(rep_tab$pval, 3)
print(rep_tab, row.names = FALSE)

# =============================================================================
# KEY COMPARISON — TSHR direction concordance across ancestries
# =============================================================================

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("KEY FINDING: Cross-Ancestry Consistency (Primary BBJ vs Replication UKB)\n")
cat(rep("=", 70), "\n\n", sep = "")

for (g in c("TSHR", "IGF1R", "TNF", "PPARG", "ARRB1", "IRS1", "AKT1", "CTLA4")) {
    pri <- mr_rescaled[mr_rescaled$gene == g &
        mr_rescaled$outcome_role == "Primary" &
        (mr_rescaled$method %in% c("Inverse variance weighted", "Wald ratio")), ]
    rep <- mr_rescaled[mr_rescaled$gene == g &
        mr_rescaled$outcome_role == "Replication" &
        (mr_rescaled$method %in% c("Inverse variance weighted", "Wald ratio")), ]

    if (nrow(pri) > 0 && nrow(rep) > 0) {
        dir_pri <- ifelse(pri$b[1] > 0, "+", "-")
        dir_rep <- ifelse(rep$b_rescaled[1] > 0, "+", "-")
        concord <- ifelse(dir_pri == dir_rep, "✅", "❌")

        cat(sprintf(
            "%-6s | Primary (BBJ)  β=%+.3f  (P=%.2e) %s\n",
            g, pri$b[1], pri$pval[1],
            ifelse(pri$pval[1] < 0.05, "SIG", "ns")
        ))
        cat(sprintf(
            "       | Repl. (UKB)*   β=%+.3f  (P=%.2e) %s  %s Direction concordance\n\n",
            rep$b_rescaled[1], rep$pval[1],
            ifelse(rep$pval[1] < 0.05, "SIG", "ns"), concord
        ))
    }
}
cat("* UKB β rescaled from linear to log-odds scale (Lloyd-Jones 2018 method)\n")

cat("\n✅ Step 2b complete. Ready for Step 3 (Colocalization).\n")
