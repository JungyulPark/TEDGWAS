# =============================================================================
# TED-TRAP Step 1 — IV Availability Pre-Check
# Purpose: Count cis-eQTL instruments available for each candidate gene
# =============================================================================

library(ieugwasr)

# --- Candidate genes ---
genes <- c(
    TSHR  = "eqtl-a-ENSG00000165409", # Primary focus (Primary EP)
    IGF1R = "eqtl-a-ENSG00000140443", # Primary focus (Primary EP)
    TNF   = "eqtl-a-ENSG00000232810", # Inflammatory
    PPARG = "eqtl-a-ENSG00000132170", # Adipogenesis
    ARRB1 = "eqtl-a-ENSG00000137486", # Scaffold
    IRS1  = "eqtl-a-ENSG00000169047", # IGF-1R adaptor
    AKT1  = "eqtl-a-ENSG00000142208", # Downstream
    CTLA4 = "eqtl-a-ENSG00000163599" # Immune regulation
)

# --- Verify token ---
token <- Sys.getenv("OPENGWAS_JWT")
cat(sprintf("Token length: %d characters\n", nchar(token)))
if (nchar(token) < 100) stop("Token not loaded. Check .Renviron.")

# --- Pre-check loop ---
cat("\n=== IV Availability at two thresholds ===\n\n")
cat(sprintf(
    "%-6s %-35s %-20s %-20s\n",
    "Gene", "OpenGWAS ID", "Strict P<5e-8", "Permissive P<5e-6"
))
cat(rep("-", 85), sep = "", "\n")

results <- data.frame(
    gene = character(), opengwas_id = character(),
    n_strict = integer(), n_permissive = integer(),
    stringsAsFactors = FALSE
)

for (g in names(genes)) {
    # Strict threshold (preferred for standard MR)
    strict <- tryCatch(
        tophits(
            id = genes[[g]], pval = 5e-8, clump = 0,
            opengwas_jwt = token
        ),
        error = function(e) NULL
    )
    n_strict <- if (!is.null(strict)) nrow(strict) else 0

    # Permissive threshold (used for cis-MR when strict is insufficient)
    perm <- tryCatch(
        tophits(
            id = genes[[g]], pval = 5e-6, clump = 0,
            opengwas_jwt = token
        ),
        error = function(e) NULL
    )
    n_perm <- if (!is.null(perm)) nrow(perm) else 0

    cat(sprintf("%-6s %-35s %-20d %-20d\n", g, genes[[g]], n_strict, n_perm))
    results <- rbind(
        results,
        data.frame(
            gene = g, opengwas_id = genes[[g]],
            n_strict = n_strict, n_permissive = n_perm,
            stringsAsFactors = FALSE
        )
    )
    Sys.sleep(0.3)
}

# --- Interpretation guide ---
cat("\n=== Interpretation ===\n")
for (i in seq_len(nrow(results))) {
    g <- results$gene[i]
    ns <- results$n_strict[i]
    np <- results$n_permissive[i]

    if (ns >= 3) {
        verdict <- "✅ Strict threshold sufficient (full sensitivity suite)"
    } else if (np >= 3) {
        verdict <- "🟡 Use permissive threshold (disclose as cis-MR)"
    } else if (np >= 2) {
        verdict <- "⚠️  cis-MR only (per-SNP Wald ratios as sensitivity)"
    } else {
        verdict <- "❌ Insufficient — consider excluding from analysis"
    }
    cat(sprintf("  %-6s: %s\n", g, verdict))
}

# --- Save summary ---
if (!dir.exists("TrackA_MR/results")) dir.create("TrackA_MR/results", recursive = TRUE)
write.csv(results, "TrackA_MR/results/01_iv_availability.csv", row.names = FALSE)
cat("\n💾 Saved: TrackA_MR/results/01_iv_availability.csv\n")
