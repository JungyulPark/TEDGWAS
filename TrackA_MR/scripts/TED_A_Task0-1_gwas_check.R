library(ieugwasr)

cat("1. OpenGWAS Graves disease search:\n")
graves_search <- gwassearch("Graves")
print(graves_search[, c("id", "trait", "sample_size", "year")])

cat("\n2. Candidate ID info:\n")
gwas_info <- gwasinfo(id = c("ieu-a-1098"))
print(gwas_info[, c("id", "trait", "sample_size", "ncase", "ncontrol", "population")])

cat("\n3. Hyperthyroidism search:\n")
hyper_search <- gwassearch("hyperthyroidism")
print(hyper_search[, c("id", "trait", "sample_size")])
