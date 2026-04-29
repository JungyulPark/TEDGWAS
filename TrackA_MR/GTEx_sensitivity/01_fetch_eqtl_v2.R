setwd("c:/ProjectTEDGWAS/TrackA_MR/GTEx_sensitivity/")

options(repos = c(CRAN = "https://cloud.r-project.org"), pkgType = "binary")

tryCatch(
    {
        if (!require("gtexr")) install.packages("gtexr")
    },
    error = function(e) cat("gtexr install failed\n")
)

tryCatch(
    {
        if (!require("data.table")) install.packages("data.table")
    },
    error = function(e) cat("data.table install failed\n")
)

tryCatch(
    {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        if (!require("xQTLbiolinks")) BiocManager::install("xQTLbiolinks", ask = FALSE, update = FALSE)
    },
    error = function(e) cat("xQTLbiolinks install failed\n")
)

library(data.table)

# Try gtexr first
cat("Attempting gtexr...\n")
success_gtexr <- FALSE
if (require("gtexr")) {
    tryCatch(
        {
            # Directly use versioned ID first
            TSHR_ENS_confirmed <- "ENSG00000165409.13"
            cat("Fetching eQTLs for", TSHR_ENS_confirmed, "in Thyroid...\n")
            sig_eqtls <- get_significant_single_tissue_eqtls(
                gencodeIds = TSHR_ENS_confirmed,
                tissueSiteDetailIds = "Thyroid",
                datasetId = "gtex_v8"
            )
            if (!is.null(sig_eqtls) && nrow(sig_eqtls) > 0) {
                fwrite(sig_eqtls, "GTEx_v8_Thyroid_TSHR_significant_eQTLs.csv")
                success_gtexr <- TRUE
                cat("Success with gtexr.\n")
            }
        },
        error = function(e) {
            cat("Error in gtexr:", e$message, "\n")
        }
    )
}

# Fallback to xQTLbiolinks
if (!success_gtexr) {
    cat("Falling back to xQTLbiolinks...\n")
    if (require("xQTLbiolinks")) {
        tryCatch(
            {
                eqtl_all <- xQTLquery_eqtl(gene = "TSHR")
                if (!is.null(eqtl_all) && nrow(eqtl_all) > 0) {
                    gtex_thyroid_tshr <- eqtl_all[tissueSiteDetail == "Thyroid"]
                    if (nrow(gtex_thyroid_tshr) > 0) {
                        fwrite(gtex_thyroid_tshr, "GTEx_v8_Thyroid_TSHR_xQTLbiolinks.csv")
                        cat("Success with xQTLbiolinks. Got", nrow(gtex_thyroid_tshr), "rows.\n")
                        success_xQTL <- TRUE
                    } else {
                        cat("No Thyroid specific results in xQTLbiolinks.\n")
                    }
                } else {
                    cat("xQTLquery_eqtl returned empty.\n")
                }
            },
            error = function(e) {
                cat("Error in xQTLbiolinks:", e$message, "\n")
            }
        )
    } else {
        cat("xQTLbiolinks package not available.\n")
    }
}
