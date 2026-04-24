library(ieugwasr)
info <- tryCatch(gwasinfo(), error = function(e) NULL)
graves <- info[grepl("Graves", info$trait, ignore.case = TRUE), c("id", "trait", "sample_size", "ncase", "nsnp", "year")]
hyper <- info[grepl("hyperthyroidism", info$trait, ignore.case = TRUE), c("id", "trait", "sample_size", "ncase", "nsnp", "year")]
write.csv(graves, "c:/ProjectTEDGWAS/TrackA_MR/results/graves_gwas.csv", row.names = FALSE)
write.csv(hyper, "c:/ProjectTEDGWAS/TrackA_MR/results/hyper_gwas.csv", row.names = FALSE)
