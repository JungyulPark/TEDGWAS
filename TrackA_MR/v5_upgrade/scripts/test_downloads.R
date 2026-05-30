# Test script for downloads

test_url <- function(url, dest) {
  cat("Testing:", url, "\n")
  res <- tryCatch({
    download.file(url, dest, mode="wb", quiet=TRUE)
    sz <- file.size(dest)
    if (!is.na(sz) && sz > 1000000) {
      cat("  SUCCESS! Size:", sz, "bytes\n")
      TRUE
    } else {
      cat("  FAILED: file too small or empty (HTML page). Size:", sz, "\n")
      if (file.exists(dest)) file.remove(dest)
      FALSE
    }
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    FALSE
  })
  res
}

# Test PLINK URLs
plink_urls <- c(
  "https://s3.amazonaws.com/plink1-assets/plink_win64_20231211.zip",
  "https://s3.amazonaws.com/plink1-assets/plink_win64_latest.zip",
  "https://www.cog-genomics.org/static/bin/plink190207/plink_win64.zip"
)
for (url in plink_urls) {
  if (test_url(url, "plink_test.zip")) break
}

# Test MAGMA URLs
magma_urls <- c(
  "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",
  "https://ctg.cncr.nl/software/magma/ref_data/g1000_eur.zip",
  "http://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",
  "http://ctg.cncr.nl/software/magma/ref_data/g1000_eur.zip"
)
for (url in magma_urls) {
  if (test_url(url, "magma_test.zip")) break
}
