# =============================================================================
# Step 3b: LD Proxy Check — rs179252 (BBJ top) vs rs1023586 (UKB top)
# =============================================================================

library(ieugwasr)
token <- Sys.getenv("OPENGWAS_JWT")

# LD check in both populations
cat("=== LD check: rs179252 (BBJ top) vs rs1023586 (UKB top) ===\n\n")

# European (1000G EUR)
cat("--- European (EUR) LD ---\n")
ld_eur <- tryCatch(
    ld_matrix(
        variants = c("rs179252", "rs1023586"),
        pop = "EUR", opengwas_jwt = token
    ),
    error = function(e) {
        cat(sprintf("EUR error: %s\n", e$message))
        NULL
    }
)
if (!is.null(ld_eur)) {
    print(ld_eur)
    cat(sprintf("\nEUR r² = %.3f\n", ld_eur[1, 2]^2))
}

# East Asian (1000G EAS)
cat("\n--- East Asian (EAS) LD ---\n")
ld_eas <- tryCatch(
    ld_matrix(
        variants = c("rs179252", "rs1023586"),
        pop = "EAS", opengwas_jwt = token
    ),
    error = function(e) {
        cat(sprintf("EAS error: %s\n", e$message))
        NULL
    }
)
if (!is.null(ld_eas)) {
    print(ld_eas)
    cat(sprintf("\nEAS r² = %.3f\n", ld_eas[1, 2]^2))
}

# Position check
cat("\n--- Position info (GRCh37) ---\n")
cat("rs179252:  chr14:81,606,536\n") # Looked up standard position
cat("rs1023586: chr14:81,598,349\n") # Looked up standard position
