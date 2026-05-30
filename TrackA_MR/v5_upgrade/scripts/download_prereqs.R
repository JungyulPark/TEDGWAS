# R script to download and setup GCTA, PLINK, and 1000G Reference Panels

# Setup folders
dir.create("c:/ProjectTEDGWAS/tools", showWarnings=FALSE, recursive=TRUE)
dir.create("c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data", showWarnings=FALSE, recursive=TRUE)

# Download and unzip PLINK
plink_url <- "https://www.cog-genomics.org/static/bin/plink190207/plink_win64.zip"
plink_zip <- "c:/ProjectTEDGWAS/tools/plink_win64.zip"
cat("Downloading PLINK from:", plink_url, "\n")
tryCatch({
  download.file(plink_url, plink_zip, mode="wb")
  cat("Extracting PLINK...\n")
  unzip(plink_zip, exdir="c:/ProjectTEDGWAS/tools")
  file.remove(plink_zip)
  cat("PLINK installed.\n")
}, error = function(e) {
  cat("ERROR downloading PLINK:", conditionMessage(e), "\n")
})

# Download and unzip GCTA
gcta_url <- "https://yanglab.westlake.edu.cn/software/gcta/download/gcta-1.94.1-Win-x86_64.zip"
gcta_zip <- "c:/ProjectTEDGWAS/tools/gcta_win64.zip"
cat("\nDownloading GCTA from:", gcta_url, "\n")
tryCatch({
  download.file(gcta_url, gcta_zip, mode="wb")
  cat("Extracting GCTA...\n")
  unzip(gcta_zip, exdir="c:/ProjectTEDGWAS/tools")
  file.remove(gcta_zip)
  # Check if directory matches expected path and rename if needed
  extracted_dirs <- list.dirs("c:/ProjectTEDGWAS/tools", recursive=FALSE)
  gcta_extracted <- grep("gcta-1.94.1-Win", extracted_dirs, value=TRUE)
  if (length(gcta_extracted) > 0) {
    file.rename(gcta_extracted, "c:/ProjectTEDGWAS/tools/gcta-1.94.1-win")
  }
  cat("GCTA installed.\n")
}, error = function(e) {
  cat("ERROR downloading GCTA:", conditionMessage(e), "\n")
})

# Download and unzip 1000G EUR reference
eur_url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip"
eur_zip <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/g1000_eur.zip"
cat("\nDownloading 1000G EUR reference from:", eur_url, "\n")
tryCatch({
  download.file(eur_url, eur_zip, mode="wb")
  cat("Extracting EUR reference...\n")
  unzip(eur_zip, exdir="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3")
  file.remove(eur_zip)
  cat("1000G EUR reference installed.\n")
}, error = function(e) {
  cat("ERROR downloading 1000G EUR reference:", conditionMessage(e), "\n")
})

# Download and unzip 1000G EAS reference
eas_url <- "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eas.zip"
eas_zip <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/g1000_eas.zip"
cat("\nDownloading 1000G EAS reference from:", eas_url, "\n")
tryCatch({
  download.file(eas_url, eas_zip, mode="wb")
  cat("Extracting EAS reference...\n")
  unzip(eas_zip, exdir="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3")
  file.remove(eas_zip)
  cat("1000G EAS reference installed.\n")
}, error = function(e) {
  cat("ERROR downloading 1000G EAS reference:", conditionMessage(e), "\n")
})

cat("\nAll downloads and extractions complete.\n")
