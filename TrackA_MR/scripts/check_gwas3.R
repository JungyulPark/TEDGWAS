library(ieugwasr)
info <- tryCatch(gwasinfo(), error = function(e) NULL)
if (is.null(info)) {
    cat("Failed to get gwasinfo()\n")
} else {
    graves <- info[grepl("Graves", info$trait, ignore.case = TRUE), ]
    cat("--- Graves ---\n")
    if (nrow(graves) > 0) print(graves[, c("id", "trait", "sample_size", "ncase", "nsnp", "year")]) else cat("None\n")

    hyper <- info[grepl("hyperthyroidism", info$trait, ignore.case = TRUE), ]
    cat("--- Hyperthyroidism ---\n")
    if (nrow(hyper) > 0) print(hyper[, c("id", "trait", "sample_size", "ncase", "nsnp", "year")]) else cat("None\n")

    thyro <- info[grepl("thyrotoxicosis", info$trait, ignore.case = TRUE), ]
    cat("--- Thyrotoxicosis ---\n")
    if (nrow(thyro) > 0) print(thyro[, c("id", "trait", "sample_size", "ncase", "nsnp", "year")]) else cat("None\n")
}
