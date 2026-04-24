# =============================================================================
# FinnGen R12 파일 검증 — 실행 전 반드시 확인
# =============================================================================

library(data.table)

FG_DIR <- "c:/ProjectTEDGWAS" # 선생님이 두신 경로

# 예상 파일 목록
expected_files <- c(
    "finngen_R12_GRAVES_OPHT.gz",
    "finngen_R12_E4_GRAVES_STRICT.gz",
    "finngen_R12_E4_GRAVES_OPHT_STRICT.gz"
)

cat("=== FinnGen R12 file verification ===\n\n")

for (f in expected_files) {
    full_path <- file.path(FG_DIR, f)

    if (!file.exists(full_path)) {
        cat(sprintf("❌ NOT FOUND: %s\n", f))
        next
    }

    size_mb <- file.info(full_path)$size / 1024^2
    cat(sprintf("✅ %s (%.1f MB)\n", f, size_mb))

    # Read first 3 lines
    cat("   Columns:\n")
    head_df <- tryCatch(
        fread(full_path, nrows = 3),
        error = function(e) {
            cat(sprintf("   [ERROR] %s\n", e$message))
            return(NULL)
        }
    )

    if (!is.null(head_df)) {
        cat(sprintf("     %s\n", paste(colnames(head_df), collapse = " | ")))
        cat("   First data row:\n")
        cat(sprintf("     %s\n\n", paste(as.character(head_df[1, ]), collapse = " | ")))
    }
}

cat("=== Verification complete ===\n")
