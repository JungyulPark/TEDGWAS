# R script to download and setup GCTA, PLINK, and 1000G Reference Panels with user-agent

options(timeout = 1800) # Set timeout to 30 mins for large downloads

download_with_agent <- function(url, dest) {
  cat("Downloading:", url, "->", dest, "\n")
  # Use curl method which handles headers well in R
  # Set user agent header
  headers <- c("User-Agent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36")
  
  res <- tryCatch({
    download.file(url, dest, mode="wb", headers=headers, method="auto")
    sz <- file.size(dest)
    cat("Download finished. File size:", sz, "bytes\n")
    if (!is.na(sz) && sz > 1000000) {
      TRUE
    } else {
      cat("  WARNING: downloaded file is too small (might be HTML page).\n")
      FALSE
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    FALSE
  })
  res
}

# Create directories
dir.create("c:/ProjectTEDGWAS/tools", showWarnings=FALSE, recursive=TRUE)
dir.create("c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data", showWarnings=FALSE, recursive=TRUE)

# 1. Download PLINK
plink_url <- "https://s3.amazonaws.com/plink1-assets/plink_win64_20231211.zip"
plink_zip <- "c:/ProjectTEDGWAS/tools/plink_win64.zip"
if (download_with_agent(plink_url, plink_zip)) {
  cat("Extracting PLINK...\n")
  unzip(plink_zip, exdir="c:/ProjectTEDGWAS/tools")
  file.remove(plink_zip)
  cat("PLINK setup complete.\n")
}

# 2. Download GCTA
gcta_url <- "https://yanglab.westlake.edu.cn/software/gcta/download/gcta-1.94.1-Win-x86_64.zip"
gcta_zip <- "c:/ProjectTEDGWAS/tools/gcta_win64.zip"
if (download_with_agent(gcta_url, gcta_zip)) {
  cat("Extracting GCTA...\n")
  unzip(gcta_zip, exdir="c:/ProjectTEDGWAS/tools")
  file.remove(gcta_zip)
  extracted_dirs <- list.dirs("c:/ProjectTEDGWAS/tools", recursive=FALSE)
  gcta_extracted <- grep("gcta-1.94.1-Win", extracted_dirs, value=TRUE)
  if (length(gcta_extracted) > 0) {
    # If target dir already exists, remove it first to avoid rename errors
    target_dir <- "c:/ProjectTEDGWAS/tools/gcta-1.94.1-win"
    if (dir.exists(target_dir)) unlink(target_dir, recursive=TRUE)
    file.rename(gcta_extracted, target_dir)
    cat("GCTA renamed to expected path.\n")
  }
  cat("GCTA setup complete.\n")
}

# 3. Download EUR reference
eur_url <- "https://zenodo.org/records/14142263/files/g1000_eur.zip?download=1"
eur_zip <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/g1000_eur.zip"
if (download_with_agent(eur_url, eur_zip)) {
  cat("Extracting EUR reference...\n")
  unzip(eur_zip, exdir="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3")
  file.remove(eur_zip)
  cat("EUR reference setup complete.\n")
}

# 4. Download EAS reference
eas_url <- "https://zenodo.org/records/14108873/files/g1000_eas.zip?download=1"
eas_zip <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/g1000_eas.zip"
if (download_with_agent(eas_url, eas_zip)) {
  cat("Extracting EAS reference...\n")
  unzip(eas_zip, exdir="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3")
  file.remove(eas_zip)
  cat("EAS reference setup complete.\n")
}

cat("\nAll pre-requisites downloaded and extracted.\n")
