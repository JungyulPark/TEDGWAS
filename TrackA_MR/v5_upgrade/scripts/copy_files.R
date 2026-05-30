src_dir <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/temp_extracted/"
dest_dir <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/"

# Create folders if not exist
dir.create(file.path(dest_dir, "00_meta"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(dest_dir, "02_tshr_finemap/results"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(dest_dir, "scripts"), showWarnings=FALSE, recursive=TRUE)

# Copy files
files_to_copy <- c(
  "data_sources_v1.2_final.csv" = "00_meta/data_sources.csv",
  "README.md" = "00_meta/README.md",
  "run_week1.sh" = "scripts/run_week1.sh",
  "taskA_01_pqtl_check.R" = "scripts/taskA_01_pqtl_check.R",
  "taskB_02_ld_clumping.R" = "scripts/taskB_02_ld_clumping.R",
  "taskB_03_gcta_cojo.sh" = "scripts/taskB_03_gcta_cojo.sh",
  "taskB_03_helper_make_ma.R" = "scripts/taskB_03_helper_make_ma.R",
  "taskB_03_helper_parse_cojo.R" = "scripts/taskB_03_helper_parse_cojo.R",
  "taskB_04_susie_finemap.R" = "scripts/taskB_04_susie_finemap.R",
  "taskB_05_integration.R" = "scripts/taskB_05_integration.R",
  "TaskB_00_ld_reference_metadata_template.csv" = "02_tshr_finemap/results/TaskB_00_ld_reference_metadata_v1.csv"
)

for (f in names(files_to_copy)) {
  src <- file.path(src_dir, f)
  dest <- file.path(dest_dir, files_to_copy[f])
  if (file.exists(src)) {
    file.copy(src, dest, overwrite=TRUE)
    cat("Copied:", f, "->", files_to_copy[f], "\n")
  } else {
    cat("ERROR: Source not found:", src, "\n")
  }
}
