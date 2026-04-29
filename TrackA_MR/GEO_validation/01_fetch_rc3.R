setwd("c:/ProjectTEDGWAS/TrackA_MR/GEO_validation/")
options(repos = c(CRAN = "https://cloud.r-project.org"), pkgType = "binary", timeout = 300)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("recount3")) BiocManager::install("recount3", ask = FALSE, update = FALSE)

library(recount3)

cat("Fetching available projects...\n")
human_projects <- available_projects()
gse_proj <- human_projects[human_projects$project == "SRP125007", ]

if (nrow(gse_proj) > 0) {
    cat("Downloading GSE105149 (SRP125007)...\n")
    rse <- create_rse(gse_proj)

    col_info <- colData(rse)
    cat("Cols available:", paste(colnames(col_info), collapse = ", "), "\n")

    # Print identifiers to allow mapping
    print(col_info[, c("external_id", "sra.sample_attributes")])

    saveRDS(rse, "GSE105149_raw_rse.rds")
} else {
    cat("SRP125007 not found in recount3 projects.\n")
}
