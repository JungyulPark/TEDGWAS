# =============================================================================
# Step 1b — IV Availability WITH Clumping (실제 MR에 쓸 숫자)
# =============================================================================

library(ieugwasr)

genes <- c(
    TSHR  = "eqtl-a-ENSG00000165409",
    IGF1R = "eqtl-a-ENSG00000140443",
    TNF   = "eqtl-a-ENSG00000232810",
    PPARG = "eqtl-a-ENSG00000132170",
    ARRB1 = "eqtl-a-ENSG00000137486",
    IRS1  = "eqtl-a-ENSG00000169047",
    AKT1  = "eqtl-a-ENSG00000142208",
    CTLA4 = "eqtl-a-ENSG00000163599"
)

token <- Sys.getenv("OPENGWAS_JWT")

cat("=== IV counts BEFORE and AFTER clumping (r² < 0.001, kb = 10000) ===\n\n")
cat(sprintf(
    "%-6s %-15s %-20s %-20s\n",
    "Gene", "Before clump", "After strict (5e-8)", "After perm (5e-6)"
))
cat(rep("-", 70), sep = "", "\n")

results_v2 <- data.frame(
    gene = character(),
    before_clump = integer(),
    after_strict = integer(),
    after_permissive = integer(),
    stringsAsFactors = FALSE
)

for (g in names(genes)) {
    # Unclumped count at strict threshold
    before <- tryCatch(
        nrow(tophits(id = genes[[g]], pval = 5e-8, clump = 0, opengwas_jwt = token)),
        error = function(e) NA
    )

    # With clumping at strict
    strict_clumped <- tryCatch(
        tophits(
            id = genes[[g]], pval = 5e-8, clump = 1,
            r2 = 0.001, kb = 10000, opengwas_jwt = token
        ),
        error = function(e) NULL
    )
    n_strict_clumped <- if (!is.null(strict_clumped)) nrow(strict_clumped) else 0

    # With clumping at permissive
    perm_clumped <- tryCatch(
        tophits(
            id = genes[[g]], pval = 5e-6, clump = 1,
            r2 = 0.001, kb = 10000, opengwas_jwt = token
        ),
        error = function(e) NULL
    )
    n_perm_clumped <- if (!is.null(perm_clumped)) nrow(perm_clumped) else 0

    cat(sprintf(
        "%-6s %-15s %-20d %-20d\n",
        g, before, n_strict_clumped, n_perm_clumped
    ))

    results_v2 <- rbind(
        results_v2,
        data.frame(
            gene = g,
            before_clump = before,
            after_strict = n_strict_clumped,
            after_permissive = n_perm_clumped,
            stringsAsFactors = FALSE
        )
    )
    Sys.sleep(1)
}

# --- Interpretation ---
cat("\n=== Final IV verdict ===\n")
for (i in seq_len(nrow(results_v2))) {
    g <- results_v2$gene[i]
    ns <- results_v2$after_strict[i]
    np <- results_v2$after_permissive[i]

    verdict <- if (ns >= 3) {
        sprintf("✅ Strict OK: %d IVs → Full sensitivity (Egger/WM/WMode/PRESSO if ≥4)", ns)
    } else if (ns == 2) {
        sprintf("🟡 Strict borderline: 2 IVs → IVW + per-SNP Wald; consider permissive (%d)", np)
    } else if (ns == 1) {
        sprintf("🟠 Single IV at strict → Wald ratio only; try permissive (%d)", np)
    } else {
        sprintf("🔴 No strict IV → Must use permissive (%d IVs) and disclose", np)
    }
    cat(sprintf("  %-6s: %s\n", g, verdict))
}

# --- Save ---
if (!dir.exists("TrackA_MR/results")) dir.create("TrackA_MR/results", recursive = TRUE)
write.csv(results_v2, "TrackA_MR/results/01b_iv_clumped.csv", row.names = FALSE)
cat("\n💾 Saved: TrackA_MR/results/01b_iv_clumped.csv\n")
