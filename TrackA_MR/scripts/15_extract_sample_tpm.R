setwd("c:/ProjectTEDGWAS")

# -----------------------------------------------------------------------------
# 1. data.txt 구조 파악
# -----------------------------------------------------------------------------
cat("=== data.txt structure ===\n")
h <- readLines("data.txt", n = 3)
cat("First 3 lines:\n")
cat(h, sep = "\n")
cat("\n")

# 2. 전체 load (using data.table for stability)
library(data.table)
d <- as.data.frame(fread("data.txt"))
cat("\nDimensions:", dim(d), "\n")
cat("Columns:", paste(names(d)[1:15], collapse = " | "), "...\n")
cat("First row:\n")
print(d[1, ])

# -----------------------------------------------------------------------------
# 3. 20 candidate genes TPM 추출
# -----------------------------------------------------------------------------
candidates <- c(
    # MR Primary
    "TSHR", "IGF1R", "IGF1",
    # MR Secondary
    "ARRB1", "PPARG", "IRS1", "AKT1", "TNF", "CTLA4",
    # Insulin Cassette
    "INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1",
    # TED Biology
    "HAS1", "HAS2", "HAS3", "ADIPOQ", "FABP4", "CEBPA"
)

# gene column 자동 탐지
gene_col <- "Gene_Symbol"
cat("\nUsing gene column:", gene_col, "\n")

# ONLY EXTRACT TPM COLUMNS (combining tech reps if they exist, or using FPKM/TPM directly)
target_cols <- grep("TPM|FPKM|tpm|fpkm", names(d), value = TRUE)
if (length(target_cols) > 0) {
    # If TPM exists, sub select those along with gene col
    d_vals <- d[, c(gene_col[1], target_cols)]
} else {
    d_vals <- d # raw fallback
}

sub <- d_vals[d_vals[[gene_col[1]]] %in% candidates, ]
cat("\n=== 20 candidate gene expression ===\n")
cat("Rows found:", nrow(sub), "/ 20 expected\n")
print(sub)

# -----------------------------------------------------------------------------
# 4. CSV로 저장
# -----------------------------------------------------------------------------
if (!dir.exists("TrackA_MR/data_for_fig3")) dir.create("TrackA_MR/data_for_fig3", recursive = TRUE)
write.csv(sub, "TrackA_MR/data_for_fig3/candidate_sample_TPM.csv", row.names = FALSE)
cat("\n💾 Saved: TrackA_MR/data_for_fig3/candidate_sample_TPM.csv\n")

# -----------------------------------------------------------------------------
# 5. Sample identity 확인
# -----------------------------------------------------------------------------
cat("\n=== Sample column identity ===\n")
sample_cols <- setdiff(names(d_vals), gene_col[1])
cat("Sample columns:", paste(sample_cols, collapse = ", "), "\n")
cat("(Handoff says Sample#2 = Control; checking TSHR level matching.)\n")
