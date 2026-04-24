res <- tryCatch(readRDS("c:/ProjectTEDGWAS/TrackA_MR/results/coloc_finngen.rds"), error = function(e) NULL)
cat("\n=== Check 2: Coloc Original Data ===\n")
if (!is.null(res)) {
    print(res$summary)
} else {
    cat("RDS loading failed.\n")
}
