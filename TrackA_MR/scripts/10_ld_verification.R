# =============================================================================
# TED-TRAP — LD verification between rs179252 (BBJ) and rs1023586 (UKB)
# =============================================================================

if (!requireNamespace("LDlinkR", quietly = TRUE)) install.packages("LDlinkR", repos = "https://cloud.r-project.org")
library(LDlinkR)

LDLINK_TOKEN <- Sys.getenv("LDLINK_TOKEN")

populations <- c("EUR", "EAS", "AFR", "AMR", "SAS")
cat("=== LD between rs179252 and rs1023586 across populations ===\n\n")

results_list <- list()

for (pop in populations) {
    cat(sprintf("--- Population: %s ---\n", pop))
    res <- tryCatch(
        LDpair(
            var1 = "rs179252", var2 = "rs1023586",
            pop = pop, token = LDLINK_TOKEN, output = "table"
        ),
        error = function(e) {
            cat(sprintf("  [ERROR] %s\n", e$message))
            NULL
        }
    )
    if (!is.null(res)) {
        print(res)
        results_list[[pop]] <- res
    }
    cat("\n")
    Sys.sleep(1)
}

if (length(results_list) > 0) {
    summary_df <- do.call(rbind, lapply(names(results_list), function(pop) {
        r <- results_list[[pop]]
        r2_val <- if ("r2" %in% names(r)) r$r2 else NA
        d_prime_val <- if ("d_prime" %in% names(r)) r$d_prime else if ("d.prime" %in% names(r)) r$d.prime else NA
        data.frame(
            Population = pop,
            r_squared = as.numeric(r2_val),
            D_prime = as.numeric(d_prime_val),
            stringsAsFactors = FALSE
        )
    }))

    cat("\n=== Summary ===\n")
    print(summary_df)

    if (!dir.exists("TrackA_MR/results")) dir.create("TrackA_MR/results", recursive = TRUE)
    write.csv(summary_df, "TrackA_MR/results/10_LD_verification.csv", row.names = FALSE)
    cat("\n💾 Saved: TrackA_MR/results/10_LD_verification.csv\n")
}
