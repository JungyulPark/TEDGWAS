# R script to fix paths and tool configurations in the Week 1 scripts

scripts_dir <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/scripts"
files <- list.files(scripts_dir, pattern="\\.(R|sh)$", full.names=TRUE)

for (f in files) {
  content <- readLines(f, warn=FALSE)
  
  # Replace incorrect workspace paths
  new_content <- gsub("c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade", "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade", content, fixed=TRUE)
  new_content <- gsub("c:\\ProjectTEDGWAS\\v5_upgrade", "c:\\ProjectTEDGWAS\\TrackA_MR\\v5_upgrade", new_content, fixed=TRUE)
  
  # Update GCTA binary path in taskB_03_gcta_cojo.sh
  if (basename(f) == "taskB_03_gcta_cojo.sh") {
    new_content <- gsub(
      "GCTA_BIN=\"c:/ProjectTEDGWAS/tools/gcta-1.94.1-win/gcta64.exe\"",
      "GCTA_BIN=\"c:/ProjectTEDGWAS/tools/gcta-1.95.1-windows-x86_64/bin/gcta64.exe\"",
      new_content, fixed=TRUE
    )
  }
  
  if (!all(content == new_content)) {
    writeLines(new_content, f)
    cat("Updated paths in:", basename(f), "\n")
  }
}
cat("Path fixing completed.\n")
