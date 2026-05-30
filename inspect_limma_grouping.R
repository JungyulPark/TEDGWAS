metadata <- readRDS("C:/ProjectTEDGWAS/TrackA_MR/GEO_validation/GSE105149_metadata.rds")
group <- rep(NA, nrow(metadata))
for (i in 1:nrow(metadata)) {
    row_text <- tolower(paste(metadata[i, ], collapse = " "))
    if (grepl("diagnosis: ted", row_text) || grepl("active ted", row_text)) {
        group[i] <- "TED"
    } else if (grepl("diagnosis: normal", row_text) || grepl("normal|control", row_text)) {
        group[i] <- "Control"
    }
}
sink("limma_grouping_inspect.txt")
cat("=== Group Counts ===\n")
print(table(group, useNA="always"))
cat("\n=== Control Samples ===\n")
print(metadata[which(group == "Control"), c("title", "diagnosis:ch1")])
cat("\n=== TED Samples ===\n")
print(metadata[which(group == "TED"), c("title", "diagnosis:ch1")])
sink()
cat("Done!\n")
