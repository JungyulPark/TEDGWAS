setwd("c:/ProjectTEDGWAS")

cat("========== EVIDENCE MATRIX COVERAGE CHECK ==========\n\n")

# -----------------------------------------------------------------------------
# 1. MR Primary — BBJ Graves outcome
# -----------------------------------------------------------------------------
cat("### 1. MR PRIMARY ###\n")
mr_files <- list.files("TrackA_MR/results",
    pattern = "01_.*primary|primary.*mr|mr.*primary",
    full.names = TRUE, ignore.case = TRUE, recursive = TRUE
)
cat("MR primary candidate files:\n")
print(mr_files)

if (length(mr_files) > 0) {
    for (f in mr_files[grepl("\\.csv$", mr_files)]) {
        d <- read.csv(f)
        cat("\n--- File:", basename(f), "---\n")
        cat("Columns:", paste(names(d), collapse = ", "), "\n")
        gene_col <- intersect(c("exposure", "gene", "exposure_name", "Gene", "gene_symbol"), names(d))
        if (length(gene_col) > 0) {
            cat("Genes tested:\n")
            print(unique(d[[gene_col[1]]]))
        } else {
            print(head(d, 3))
        }
    }
}

# -----------------------------------------------------------------------------
# 2. Colocalization
# -----------------------------------------------------------------------------
cat("\n\n### 2. COLOCALIZATION ###\n")
coloc_files <- list.files("TrackA_MR/results",
    pattern = "coloc",
    full.names = TRUE, ignore.case = TRUE, recursive = TRUE
)
cat("Coloc files found:\n")
print(coloc_files)

for (f in coloc_files[grepl("\\.csv$", coloc_files)]) {
    d <- read.csv(f)
    cat("\n--- File:", basename(f), "---\n")
    cat("Columns:", paste(names(d), collapse = ", "), "\n")
    print(d) # coloc 결과는 row 수 적으므로 전체 출력
}

# -----------------------------------------------------------------------------
# 3. MVMR
# -----------------------------------------------------------------------------
cat("\n\n### 3. MVMR (Mediation) ###\n")
mvmr_files <- list.files("TrackA_MR/results",
    pattern = "mvmr|multivariable|mediation",
    full.names = TRUE, ignore.case = TRUE, recursive = TRUE
)
cat("MVMR files found:\n")
print(mvmr_files)

for (f in mvmr_files[grepl("\\.csv$", mvmr_files)]) {
    d <- read.csv(f)
    cat("\n--- File:", basename(f), "---\n")
    cat("Columns:", paste(names(d), collapse = ", "), "\n")
    print(d)
}

# -----------------------------------------------------------------------------
# 4. Tissue RNA-seq (candidate DESeq2 on 20 genes)
# -----------------------------------------------------------------------------
cat("\n\n### 4. TISSUE RNA-seq (Candidate DESeq2) ###\n")
tissue_files <- list.files("TrackA_MR/results",
    pattern = "rnaseq|deseq|candidate|tissue|06_",
    full.names = TRUE, ignore.case = TRUE, recursive = TRUE
)
cat("Tissue files found:\n")
print(tissue_files)

for (f in tissue_files[grepl("\\.csv$", tissue_files)]) {
    d <- tryCatch(read.csv(f), error = function(e) NULL)
    if (is.null(d)) next
    cat("\n--- File:", basename(f), "---\n")
    cat("Rows:", nrow(d), "| Columns:", paste(names(d), collapse = ", "), "\n")
    gene_col <- intersect(c("gene", "Gene", "gene_symbol", "symbol"), names(d))
    if (length(gene_col) > 0) {
        cat("Genes tested:\n")
        print(unique(d[[gene_col[1]]]))
        if (nrow(d) <= 25) {
            cat("\nFull content:\n")
            print(d)
        }
    } else {
        print(head(d, 5))
    }
}

cat("\n\n========== END COVERAGE CHECK ==========\n")
