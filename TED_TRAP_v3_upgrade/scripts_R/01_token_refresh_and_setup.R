#==============================================================================
# TED-TRAP Upgrade — Phase 1 Script 01
# OpenGWAS JWT Token Refresh + Package Installation
#==============================================================================
# Purpose:
#   - Install/load all R packages required for upgraded MR analysis
#   - Refresh OpenGWAS JWT token (resolves 401 error seen in terminal)
#   - Verify API access before running downstream scripts
#
# Author: TED-TRAP team
# Output: Token stored in .Renviron; logs to logs/01_setup.log
#==============================================================================

# --- Configuration ---
PROJECT_ROOT <- "c:/ProjectTEDGWAS"   # adjust if needed
setwd(PROJECT_ROOT)

# Create directory structure if not already present
dirs_needed <- c(
  "TrackA_MR/scripts", "TrackA_MR/data", "TrackA_MR/results",
  "TrackA_MR/logs", "TrackA_MR/figures",
  "TrackB_Network/scripts", "TrackB_Network/results",
  "TrackC_Offtarget/scripts", "TrackC_Offtarget/results",
  "Manuscript/archive_v1", "Manuscript/figures_v3", "Manuscript/tables_v3"
)
for (d in dirs_needed) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("✅ Directory structure verified.\n\n")

# --- Step 1: Package installation ---
# Required packages (core MR + sensitivity suite + coloc + DESeq2 prep)
required_CRAN <- c(
  "devtools", "remotes",
  "TwoSampleMR", "ieugwasr",
  "MendelianRandomization",
  "coloc",
  "susieR",
  "data.table", "dplyr", "tidyr", "readr", "stringr",
  "ggplot2", "cowplot", "ggrepel",
  "openxlsx"
)

# GitHub packages (not on CRAN)
required_GH <- list(
  "rondolab/MR-PRESSO" = "MRPRESSO",
  "n-mounier/MRlap"    = "MRlap",
  "WSpiller/MVMR"      = "MVMR",
  "jean997/cause"      = "cause"
)

# Bioconductor
required_Bioc <- c("DESeq2", "edgeR", "limma", "apeglm", "IHW")

cat("📦 Installing CRAN packages...\n")
for (pkg in required_CRAN) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("   Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("   [OK] %s already installed\n", pkg))
  }
}

cat("\n📦 Installing GitHub packages...\n")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
for (repo in names(required_GH)) {
  pkg_name <- required_GH[[repo]]
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat(sprintf("   Installing %s from %s...\n", pkg_name, repo))
    tryCatch(
      remotes::install_github(repo, upgrade = "never"),
      error = function(e) cat(sprintf("   [WARN] Failed to install %s: %s\n", pkg_name, e$message))
    )
  } else {
    cat(sprintf("   [OK] %s already installed\n", pkg_name))
  }
}

cat("\n📦 Installing Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_Bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("   Installing %s from Bioconductor...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("   [OK] %s already installed\n", pkg))
  }
}

# --- Step 2: OpenGWAS JWT Token Setup ---
cat("\n" , rep("=", 70), "\n", sep = "")
cat("🔑 OpenGWAS JWT Token Setup\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Current error in your log: 'Status code from OpenGWAS API: 401'\n")
cat("This means your JWT token has expired.\n\n")

cat("To refresh:\n")
cat("1. Open browser: https://api.opengwas.io/profile/\n")
cat("2. Sign in (ORCID or Google)\n")
cat("3. Click 'Generate new token'\n")
cat("4. Copy the token (very long string starting with 'eyJ...')\n")
cat("5. Paste it when prompted below.\n\n")

# Interactive token entry
if (interactive()) {
  new_token <- readline(prompt = "Paste your new OpenGWAS JWT token (or press ENTER to skip): ")

  if (nchar(new_token) > 20) {
    # Write to .Renviron in project root
    renviron_path <- file.path(PROJECT_ROOT, ".Renviron")
    existing_lines <- if (file.exists(renviron_path)) readLines(renviron_path) else character(0)
    existing_lines <- existing_lines[!grepl("^OPENGWAS_JWT=", existing_lines)]
    existing_lines <- c(existing_lines, sprintf('OPENGWAS_JWT=%s', new_token))
    writeLines(existing_lines, renviron_path)

    cat(sprintf("✅ Token saved to: %s\n", renviron_path))
    cat("⚠️  Restart R session (Ctrl+Shift+F10 in RStudio, or close+reopen R Interactive in VS Code)\n")
    cat("    to load the token into environment.\n")

    # Immediate in-session activation
    Sys.setenv(OPENGWAS_JWT = new_token)
    cat("✅ Token also activated for the current R session.\n")
  } else {
    cat("⚠️  No token entered. You can run this script again anytime.\n")
  }
}

# --- Step 3: Verify API access ---
cat("\n🔌 Verifying OpenGWAS API access...\n")
library(ieugwasr)
token <- Sys.getenv("OPENGWAS_JWT")

if (nchar(token) < 20) {
  cat("⚠️  No token in environment. Please set OPENGWAS_JWT in .Renviron and restart R.\n")
} else {
  api_test <- tryCatch(
    {
      # Simple API test: list user info
      user_status <- ieugwasr::user()
      cat("✅ OpenGWAS API access verified.\n")
      cat(sprintf("   User: %s\n", user_status$user))
      cat(sprintf("   Rate limit remaining: %s\n", user_status$allowance))
      TRUE
    },
    error = function(e) {
      cat(sprintf("❌ API test failed: %s\n", e$message))
      FALSE
    }
  )

  # Test a lightweight query
  if (api_test) {
    cat("\n🧪 Testing GWAS query (fetching TSHR info)...\n")
    test_query <- tryCatch(
      ieugwasr::gwasinfo(id = "eqtl-a-ENSG00000165409"),   # TSHR
      error = function(e) NULL
    )
    if (!is.null(test_query)) {
      cat(sprintf("✅ Test query successful. Trait: %s\n",
                  test_query$trait))
    } else {
      cat("⚠️  Test query failed. Check token validity.\n")
    }
  }
}

# --- Step 4: Session info log ---
log_path <- "TrackA_MR/logs/01_setup.log"
sink(log_path)
cat("=== TED-TRAP Setup Log ===\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))
cat("\n=== Installed key packages ===\n")
for (pkg in c("TwoSampleMR", "ieugwasr", "coloc", "susieR",
              "MRPRESSO", "MRlap", "MVMR", "cause", "DESeq2")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
  } else {
    cat(sprintf("  %s: NOT INSTALLED\n", pkg))
  }
}
sink()

cat(sprintf("\n📝 Setup log saved: %s\n", log_path))

cat("\n", rep("=", 70), "\n", sep = "")
cat("✅ Phase 1 Setup Complete\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nNext step: Run 02_eqtlgen_exposure_extraction_all_genes.R\n\n")
