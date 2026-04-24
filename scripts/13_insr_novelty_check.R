#==============================================================================
# TED-TRAP Upgrade — Phase 2 Script 13
# INSR upregulation novelty verification against prior TED RNA-seq reports
#==============================================================================
# Purpose:
#   Resolve CRITICAL issue #6: The v1 manuscript claimed INSR upregulation
#   as "first report of INSR in TED orbital tissue." This script systematically
#   verifies this claim against two key prior reports:
#
#     - Kim et al. 2021, IOVS 62(9):24 (doi:10.1167/iovs.62.9.24)
#       Bulk RNA-seq, n=6 TAO vs n=6 Ctrl orbital fat
#
#     - Kim et al. 2024, JCI Insight 9(24):e182352 (doi:10.1172/jci.insight.182352)
#       Single-nucleus RNA-seq, n=3 TED orbital fat
#
#   If INSR appears in either supplementary DEG table, the "first report" claim
#   must be softened to "direct transcriptomic confirmation complementing prior
#   pathway-level reports."
#
# Approach:
#   1. Download supplementary tables from both papers
#   2. Parse DEG lists / adjacent summary files
#   3. Search for "INSR" in gene columns
#   4. Compare effect direction
#
# Output:
#   TrackA_MR/results/13_insr_novelty_verification.md
#==============================================================================

setwd("c:/ProjectTEDGWAS")
library(data.table)
library(dplyr)

log_file <- "TrackA_MR/logs/13_insr_novelty.log"
sink(log_file, split = TRUE)
cat("=== INSR Novelty Verification ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

report_path <- "TrackA_MR/results/13_insr_novelty_verification.md"

# --- Search paths ---
kim2021_supp  <- "TrackA_MR/data/literature_checks/kim2021_iovs_supplementary.xlsx"
kim2024_supp  <- "TrackA_MR/data/literature_checks/kim2024_jci_insight_supplementary.xlsx"

cat("Looking for downloaded supplementary tables at:\n")
cat(sprintf("  Kim 2021 IOVS: %s\n", kim2021_supp))
cat(sprintf("  Kim 2024 JCI Insight: %s\n\n", kim2024_supp))

# --- Step 1: Kim 2021 IOVS check ---
result_2021 <- list(
  found = NA,
  effect = NA,
  source = "Kim et al. 2021, IOVS 62(9):24",
  notes = ""
)

if (file.exists(kim2021_supp)) {
  cat("Checking Kim 2021...\n")
  tryCatch({
    library(readxl)
    sheets <- excel_sheets(kim2021_supp)
    cat(sprintf("  Sheets found: %s\n", paste(sheets, collapse=", ")))

    for (sh in sheets) {
      df <- read_excel(kim2021_supp, sheet = sh)
      gene_cols <- grep("gene|symbol", names(df), ignore.case = TRUE, value = TRUE)
      if (length(gene_cols) == 0) next

      for (gc in gene_cols) {
        hit <- df[toupper(df[[gc]]) == "INSR", ]
        if (nrow(hit) > 0) {
          cat(sprintf("  ✅ INSR found in sheet '%s'\n", sh))
          cat("     Row(s):\n")
          print(hit)
          result_2021$found <- TRUE
          result_2021$notes <- paste0("Found in sheet: ", sh)
          # Extract effect direction if column exists
          fc_col <- grep("log2|logfc|fold", names(hit), ignore.case = TRUE, value = TRUE)
          if (length(fc_col) > 0) {
            result_2021$effect <- as.numeric(hit[[fc_col[1]]][1])
          }
        }
      }
    }
    if (is.na(result_2021$found)) {
      result_2021$found <- FALSE
      result_2021$notes <- "INSR not found in any sheet/column"
      cat("  ⚠️  INSR NOT found in Kim 2021 supplementary.\n")
    }
  }, error = function(e) {
    result_2021$notes <- sprintf("Error reading file: %s", e$message)
    cat(sprintf("  [ERROR] %s\n", e$message))
  })
} else {
  result_2021$notes <- "Supplementary file not downloaded yet"
  cat("⚠️  Kim 2021 supplementary not found locally.\n")
  cat("   Download from: https://iovs.arvojournals.org/article.aspx?articleid=2777099\n")
  cat("   (Supplementary materials tab)\n")
}

# --- Step 2: Kim 2024 JCI Insight check ---
result_2024 <- list(
  found = NA,
  effect = NA,
  source = "Kim et al. 2024, JCI Insight 9:e182352",
  notes = ""
)

if (file.exists(kim2024_supp)) {
  cat("\nChecking Kim 2024...\n")
  tryCatch({
    library(readxl)
    sheets <- excel_sheets(kim2024_supp)
    cat(sprintf("  Sheets found: %s\n", paste(sheets, collapse=", ")))

    for (sh in sheets) {
      df <- read_excel(kim2024_supp, sheet = sh)
      gene_cols <- grep("gene|symbol", names(df), ignore.case = TRUE, value = TRUE)
      if (length(gene_cols) == 0) next

      for (gc in gene_cols) {
        hit <- df[toupper(df[[gc]]) == "INSR", ]
        if (nrow(hit) > 0) {
          cat(sprintf("  ✅ INSR found in sheet '%s'\n", sh))
          print(hit)
          result_2024$found <- TRUE
          result_2024$notes <- paste0("Found in sheet: ", sh)
          fc_col <- grep("log2|logfc|fold", names(hit), ignore.case = TRUE, value = TRUE)
          if (length(fc_col) > 0) {
            result_2024$effect <- as.numeric(hit[[fc_col[1]]][1])
          }
        }
      }
    }
    if (is.na(result_2024$found)) {
      result_2024$found <- FALSE
      result_2024$notes <- "INSR not found in any sheet/column"
      cat("  ⚠️  INSR NOT found in Kim 2024 supplementary.\n")
    }
  }, error = function(e) {
    result_2024$notes <- sprintf("Error reading file: %s", e$message)
    cat(sprintf("  [ERROR] %s\n", e$message))
  })
} else {
  result_2024$notes <- "Supplementary file not downloaded yet"
  cat("\n⚠️  Kim 2024 supplementary not found locally.\n")
  cat("   Download from: https://insight.jci.org/articles/view/182352\n")
  cat("   (Supplementary materials section)\n")
}

# --- Step 3: Consolidated verdict ---
any_prior <- isTRUE(result_2021$found) || isTRUE(result_2024$found)

verdict <- if (any_prior) {
  "DOWNGRADE claim: INSR has been reported before."
} else if (is.na(result_2021$found) && is.na(result_2024$found)) {
  "INCONCLUSIVE: Supplementary files not yet verified. Use cautious framing."
} else {
  "UPHOLD claim: INSR appears to be a novel addition to TED RNA-seq literature."
}

cat("\n\n=== VERDICT ===\n")
cat(verdict, "\n\n")

# --- Write markdown report ---
md <- c(
  "# INSR Novelty Verification Report",
  sprintf("Date: %s\n", Sys.time()),
  "",
  "## Purpose",
  "",
  "The v1 manuscript claimed INSR upregulation as 'first report of INSR in TED orbital tissue'.",
  "This report verifies that claim against prior TED/TAO transcriptomic studies.",
  "",
  "## Sources checked",
  "",
  sprintf("- **%s**", result_2021$source),
  sprintf("  - Found: %s", ifelse(is.na(result_2021$found), "UNKNOWN", result_2021$found)),
  sprintf("  - Effect: %s", ifelse(is.na(result_2021$effect), "N/A", result_2021$effect)),
  sprintf("  - Notes: %s", result_2021$notes),
  "",
  sprintf("- **%s**", result_2024$source),
  sprintf("  - Found: %s", ifelse(is.na(result_2024$found), "UNKNOWN", result_2024$found)),
  sprintf("  - Effect: %s", ifelse(is.na(result_2024$effect), "N/A", result_2024$effect)),
  sprintf("  - Notes: %s", result_2024$notes),
  "",
  "## Verdict",
  "",
  paste0("**", verdict, "**"),
  "",
  "## Recommended manuscript language",
  ""
)

if (any_prior) {
  md <- c(md,
    "> *We provide direct transcriptomic evidence of INSR upregulation in TED orbital adipose tissue, complementing prior pathway-level reports of insulin-signaling enrichment in TED/TAO orbital transcriptomes (Kim 2021; Kim 2024).*",
    "",
    "**Do NOT use** 'first report' or 'novel' framing for INSR.")
} else if (is.na(result_2021$found) && is.na(result_2024$found)) {
  md <- c(md,
    "> *We report INSR upregulation in TED orbital adipose tissue, an observation that, to our knowledge, has not been previously highlighted in TED-specific transcriptomic reports, though prior studies have documented broader insulin-signaling pathway enrichment.*",
    "",
    "**Use cautious, non-absolute language**.")
} else {
  md <- c(md,
    "> *We report INSR upregulation in TED orbital adipose tissue; this gene-level finding has not been previously reported in available TED transcriptomic studies that we reviewed.*",
    "",
    "**'First report' language permissible but soften to 'to our knowledge' / 'in available studies'.**")
}

writeLines(md, report_path)
cat(sprintf("\n📝 Report saved: %s\n", report_path))

sink()
