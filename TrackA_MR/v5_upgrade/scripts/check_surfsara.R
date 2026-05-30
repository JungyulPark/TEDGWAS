# R script to find the filenames for the SurfSara MAGMA download links

urls <- c(
  "https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe",
  "https://vu.data.surfsara.nl/index.php/s/ePXET6IWVTwTes4",
  "https://vu.data.surfsara.nl/index.php/s/dz6PYdKOi3xVqHn",
  "https://vu.data.surfsara.nl/index.php/s/C6UkTV5nuFo8cJC",
  "https://vu.data.surfsara.nl/index.php/s/TXDEm70eEO7AgOb",
  "https://vu.data.surfsara.nl/index.php/s/alUE9KnBwCeml2S"
)

# Append /download to get the direct download link from SurfSara
for (url in urls) {
  direct_url <- paste0(url, "/download")
  cat("Checking URL:", direct_url, "\n")
  h <- curlGetHeaders(direct_url)
  disp <- grep("content-disposition", h, ignore.case=TRUE, value=TRUE)
  if (length(disp) > 0) {
    cat("  Disposition:", disp, "\n")
  } else {
    # Check location header for redirect
    loc <- grep("Location:", h, ignore.case=TRUE, value=TRUE)
    cat("  No disposition. Headers: ", paste(head(h, 5), collapse=" | "), "\n")
  }
}
