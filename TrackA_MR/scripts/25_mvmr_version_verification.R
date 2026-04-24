setwd("c:/ProjectTEDGWAS")

cat("========== MVMR VERSION VERIFICATION ==========\n\n")

# -----------------------------------------------------------------------------
# 1. File timestamps
# -----------------------------------------------------------------------------
cat("### File timestamps ###\n")
files_to_check <- c(
    "TrackA_MR/results/04_mvmr_summary.csv",
    "TrackA_MR/results/04v3_mvmr_finngen_summary.csv",
    "TrackA_MR/results/04_mvmr.log"
)
for (f in files_to_check) {
    if (file.exists(f)) {
        info <- file.info(f)
        cat(sprintf("%-60s mtime: %s\n", basename(f), info$mtime))
    } else {
        cat(sprintf("%-60s NOT FOUND\n", basename(f)))
    }
}

# -----------------------------------------------------------------------------
# 2. v3 file contents (full)
# -----------------------------------------------------------------------------
cat("\n### 04v3_mvmr_finngen_summary.csv (full content) ###\n")
if (file.exists("TrackA_MR/results/04v3_mvmr_finngen_summary.csv")) {
    v3 <- read.csv("TrackA_MR/results/04v3_mvmr_finngen_summary.csv", stringsAsFactors = FALSE)
    print(v3)
}

# -----------------------------------------------------------------------------
# 3. Find v3 generator script
# -----------------------------------------------------------------------------
cat("\n### Find script that generates 04v3 ###\n")
all_scripts <- list.files("TrackA_MR/scripts",
    pattern = "\\.R$",
    full.names = TRUE, recursive = TRUE
)
for (f in all_scripts) {
    txt <- readLines(f, warn = FALSE)
    # Look for writes to 04v3 or specific variables
    hits <- grep("04v3|mvmr_finngen|GRAVES_OPHT|E4_GRAVES_OPHT|GO_outcome",
        txt,
        ignore.case = TRUE
    )
    if (length(hits) > 0) {
        cat("\n--- ", basename(f), " ---\n")
        for (h in hits) {
            lo <- max(1, h - 2)
            hi <- min(length(txt), h + 3)
            cat(paste0("[L", lo:hi, "] ", txt[lo:hi]), sep = "\n")
            cat("...\n")
        }
    }
}

# -----------------------------------------------------------------------------
# 4. Read v3's log if exists
# -----------------------------------------------------------------------------
cat("\n### Any v3-specific log? ###\n")
v3logs <- list.files("TrackA_MR",
    pattern = "04v3|mvmr_finngen",
    full.names = TRUE, recursive = TRUE, ignore.case = TRUE
)
v3logs <- v3logs[grepl("\\.log$", v3logs, ignore.case = TRUE)]
cat("Logs:\n")
print(basename(v3logs))
for (f in v3logs) {
    cat("\n--- ", basename(f), " ---\n")
    cat(head(readLines(f, warn = FALSE), 80), sep = "\n")
}

# -----------------------------------------------------------------------------
# 5. 04_mvmr.log — full dump
# -----------------------------------------------------------------------------
cat("\n### 04_mvmr.log (full) ###\n")
if (file.exists("TrackA_MR/results/04_mvmr.log")) {
    cat(readLines("TrackA_MR/results/04_mvmr.log"), sep = "\n")
}

# -----------------------------------------------------------------------------
# 6. DIAGNOSIS
# -----------------------------------------------------------------------------
cat("\n\n========== DIAGNOSIS ==========\n")
cat("Key questions:\n")
cat("1. Which file is NEWER (timestamps)?\n")
cat("2. What outcome did v3 actually use (from its generator script)?\n")
cat("3. Why was 04_mvmr re-run with E4_THYROID if v3 already had E4_GRAVES_OPHT?\n")
cat("4. Is 04v3 the legitimate final MVMR, or a discarded intermediate?\n")
