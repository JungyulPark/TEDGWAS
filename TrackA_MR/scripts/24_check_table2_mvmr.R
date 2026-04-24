setwd("c:/ProjectTEDGWAS")

# 의심 값들
suspicious <- c(
    "0.127", "0.962", # TSHR / GD_liability β_adj (Table 2 Panel C)
    "0.277", "0.074", # SEs
    "0.657", # TSHR P_adj
    "3.67e-07", "3.67e", # GD liability P_adj
    "25.4", "99.3", # Conditional F values
    "1.566", "-1.566", # Univariable TSHR β
    "2.82e-07", "2.82e" # Univariable TSHR P
)

# 모든 결과 파일 훑어서 이 숫자들 찾기
all_files <- c(
    list.files("TrackA_MR/results", full.names = TRUE, recursive = TRUE),
    list.files("TrackA_MR/manuscript", full.names = TRUE, recursive = TRUE),
    list.files("TrackA_MR/scripts", full.names = TRUE, recursive = TRUE)
)
all_files <- all_files[grepl("\\.(csv|txt|log|R|md|docx)$", all_files, ignore.case = TRUE)]

cat("=== Searching for suspicious numbers across all files ===\n\n")
for (val in suspicious) {
    cat("--- Searching for:", val, "---\n")
    found_any <- FALSE
    for (f in all_files) {
        if (grepl("\\.docx$", f, ignore.case = TRUE)) next # Skip docx parsing via readLines
        txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
        hits <- grep(val, txt, fixed = TRUE)
        if (length(hits) > 0) {
            cat("  FOUND in:", basename(f), "(lines:", paste(head(hits, 3), collapse = ","), ")\n")
            # Show context (first hit)
            ctx <- txt[max(1, hits[1] - 1):min(length(txt), hits[1] + 1)]
            cat("    Context:", paste(ctx, collapse = " | "), "\n")
            found_any <- TRUE
        }
    }
    if (!found_any) cat("  NOT FOUND anywhere\n")
    cat("\n")
}

# -----------------------------------------------------------------------------
# Step 2: Also grep MVMR scripts for how β was calculated
# -----------------------------------------------------------------------------
cat("\n=== MVMR-related scripts ===\n")
mvmr_scripts <- list.files("TrackA_MR/scripts",
    pattern = "mvmr|MVMR|04_",
    full.names = TRUE, ignore.case = TRUE
)
cat("Scripts found:\n")
print(basename(mvmr_scripts))

for (f in mvmr_scripts) {
    cat("\n--- ", basename(f), " ---\n")
    txt <- readLines(f, warn = FALSE)
    # Show last 50 lines (usually contains output)
    cat(tail(txt, 50), sep = "\n")
}

# -----------------------------------------------------------------------------
# Step 3: Any MVMR log files?
# -----------------------------------------------------------------------------
cat("\n=== MVMR log files ===\n")
logs <- list.files("TrackA_MR",
    pattern = "mvmr|04_",
    full.names = TRUE, recursive = TRUE, ignore.case = TRUE
)
logs <- logs[grepl("\\.(log|txt)$", logs)]
cat("Log files:\n")
print(basename(logs))
for (f in logs) {
    cat("\n--- ", basename(f), " ---\n")
    cat(readLines(f, warn = FALSE)[1:min(100, length(readLines(f, warn = FALSE)))], sep = "\n")
}
