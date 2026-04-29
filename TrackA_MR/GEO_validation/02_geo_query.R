setwd("c:/ProjectTEDGWAS/TrackA_MR/GEO_validation/")
options(repos = c(CRAN = "https://cloud.r-project.org"), pkgType = "binary")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GEOquery")) BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
if (!require("limma")) BiocManager::install("limma", ask = FALSE, update = FALSE)

library(GEOquery)
library(limma)
library(data.table)

cat("Fetching GSE105149 via GEOquery...\n")
gse <- getGEO("GSE105149", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gse) > 0) {
    eset <- gse[[1]]
    exprs_mat <- exprs(eset)
    metadata <- pData(eset)

    cat("Matrix dimensions:", dim(exprs_mat), "\n")
    cat("Sample IDs:", paste(rownames(metadata), collapse = ", "), "\n")

    # Determine case/control (usually in title or characteristics_ch1)
    cols_to_check <- c("title", "characteristics_ch1", "characteristics_ch1.1")
    available_cols <- intersect(cols_to_check, names(metadata))
    print(metadata[, available_cols, drop = FALSE])

    # Heuristic assignment
    metadata$group <- "Control"
    for (i in 1:nrow(metadata)) {
        row_text <- paste(metadata[i, ], collapse = " ")
        if (grepl("active ted|graves|patient", row_text, ignore.case = TRUE)) {
            metadata$group[i] <- "TED"
        } else if (grepl("inactive", row_text, ignore.case = TRUE)) {
            metadata$group[i] <- "InactiveTED"
        }
    }

    cat("Groups assigned:\n")
    print(table(metadata$group))

    # Save for further analysis
    saveRDS(eset, "GSE105149_eset.rds")
    saveRDS(metadata, "GSE105149_metadata.rds")
} else {
    cat("Failed to retrieve GSE105149.\n")
}
