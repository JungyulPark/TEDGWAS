# =============================================================================
# TED-TRAP Table 3 — IGF1R Cross-phenotype MR (off-target probe)
# Purpose: Test whether IGF1R germline perturbation causes metabolic phenotypes
# Positive control: Height (IGF-1/growth axis)
# =============================================================================

suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(dplyr)
    library(ieugwasr)
    library(writexl)
})

setwd("c:/ProjectTEDGWAS/TrackA_MR")

if (!dir.exists("results/Table3")) dir.create("results/Table3", recursive = TRUE)
if (!dir.exists("logs")) dir.create("logs", recursive = TRUE)

log_file <- file.path("logs", paste0("Table3_IGF1R_xphen_", format(Sys.Date(), "%Y%m%d"), ".log"))
sink(log_file, split = TRUE)

cat("=== IGF1R Cross-phenotype MR — started", as.character(Sys.time()), "===\n\n")

# -----------------------------------------------------------------------------
# STEP 1: Load IGF1R exposure (재사용 — Primary MR에서 쓴 것과 동일)
# -----------------------------------------------------------------------------
igf1r_iv_file <- "data/instruments/IGF1R_instruments.tsv" # From Script 02 Core

if (file.exists(igf1r_iv_file)) {
    exp_igf1r <- read.table(igf1r_iv_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    cat("Loaded existing IGF1R IVs from TSV:", nrow(exp_igf1r), "SNPs\n")
} else {
    stop("Instrument file missing")
}

cat("IGF1R instruments:\n")
print(exp_igf1r %>% select(SNP, beta.exposure, pval.exposure))

# -----------------------------------------------------------------------------
# STEP 2: Define 4 outcomes
# -----------------------------------------------------------------------------
outcomes_map <- tibble::tribble(
    ~label,              ~category,           ~opengwas_id,              ~trait_type,
    "Height",            "Positive control",  "ieu-a-89",                "continuous", # GIANT 2014
    "Fasting glucose",   "Glycemic",          "ebi-a-GCST90002232",      "continuous", # MAGIC 2021
    "HbA1c",             "Glycemic",          "ebi-a-GCST90002244",      "continuous", # MAGIC 2021
    "Type 2 diabetes",   "Glycemic",          "ebi-a-GCST90018926",      "binary" # DIAGRAM
)

cat("\nOutcomes:\n")
print(outcomes_map)

# -----------------------------------------------------------------------------
# STEP 3: MR loop
# -----------------------------------------------------------------------------
results_all <- list()
heterogeneity_all <- list()
pleiotropy_all <- list()
steiger_all <- list()

for (i in seq_len(nrow(outcomes_map))) {
    lbl <- outcomes_map$label[i]
    gid <- outcomes_map$opengwas_id[i]
    cat("\n==================================================\n")
    cat("Outcome:", lbl, "(", gid, ")\n")
    cat("==================================================\n")

    # Outcome extract
    out_dat <- tryCatch(
        extract_outcome_data(snps = exp_igf1r$SNP, outcomes = gid, proxies = TRUE),
        error = function(e) {
            cat("  ❌ extract_outcome_data failed:", e$message, "\n")
            NULL
        }
    )

    if (is.null(out_dat) || nrow(out_dat) == 0) {
        # Attempt alternative ID
        alt_gid <- switch(lbl,
            "Height" = "ukb-b-10787",
            "Fasting glucose" = "ieu-b-113",
            "HbA1c" = "ukb-b-4424",
            "Type 2 diabetes" = "ieu-a-26"
        )
        cat("  ⚠️  Failed. Attempting alt ID:", alt_gid, "\n")
        out_dat <- tryCatch(
            extract_outcome_data(snps = exp_igf1r$SNP, outcomes = alt_gid, proxies = TRUE),
            error = function(e) {
                NULL
            }
        )
    }

    if (is.null(out_dat) || nrow(out_dat) == 0) {
        cat("  ⚠️  No outcome SNPs — skipped\n")
        next
    }

    # Harmonise
    dat <- harmonise_data(exp_igf1r, out_dat, action = 2)
    dat <- dat[dat$mr_keep == TRUE, ]
    cat("  Harmonised IVs:", nrow(dat), "\n")
    if (nrow(dat) == 0) next

    # MR
    mr_res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_wald_ratio"))
    mr_res$outcome_label <- lbl
    mr_res$category <- outcomes_map$category[i]
    mr_res$trait_type <- outcomes_map$trait_type[i]
    results_all[[lbl]] <- mr_res

    cat("\nMR results:\n")
    print(mr_res %>% select(method, nsnp, b, se, pval))

    # Sensitivity
    if (nrow(dat) >= 3) {
        het <- mr_heterogeneity(dat)
        het$outcome_label <- lbl
        heterogeneity_all[[lbl]] <- het
        plt <- mr_pleiotropy_test(dat)
        plt$outcome_label <- lbl
        pleiotropy_all[[lbl]] <- plt
    }
    stg <- directionality_test(dat)
    stg$outcome_label <- lbl
    steiger_all[[lbl]] <- stg

    Sys.sleep(1) # OpenGWAS rate limit
}

# -----------------------------------------------------------------------------
# STEP 4: Consolidate
# -----------------------------------------------------------------------------
results_df <- bind_rows(results_all)
heterogeneity_df <- bind_rows(heterogeneity_all)
pleiotropy_df <- bind_rows(pleiotropy_all)
steiger_df <- bind_rows(steiger_all)

# Main IVW (or Wald if nsnp=1) row per outcome for Table 3
primary_rows <- results_df %>%
    group_by(outcome_label) %>%
    filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
        OR = ifelse(trait_type == "binary", sprintf("%.2f", exp(b)), NA),
        OR_CI = ifelse(trait_type == "binary",
            sprintf("%.2f (%.2f–%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)),
            NA
        ),
        beta_CI = sprintf("%.3f (%.3f to %.3f)", b, b - 1.96 * se, b + 1.96 * se)
    )

# WM + Egger per outcome
wm_rows <- results_df %>%
    filter(method == "Weighted median") %>%
    select(outcome_label, pval_WM = pval)
eg_rows <- results_df %>%
    filter(method == "MR Egger") %>%
    select(outcome_label, pval_Egger = pval)

table3 <- primary_rows %>%
    select(outcome_label, category, trait_type, nsnp, b, se, OR_CI, beta_CI, pval_IVW = pval) %>%
    left_join(wm_rows, by = "outcome_label") %>%
    left_join(eg_rows, by = "outcome_label") %>%
    left_join(steiger_df %>% select(outcome_label, correct_causal_direction, steiger_pval), by = "outcome_label")

cat("\n\n=== TABLE 3 (condensed) ===\n")
print(table3)

# -----------------------------------------------------------------------------
# STEP 5: Save
# -----------------------------------------------------------------------------
write.csv(results_df, "results/Table3/MR_all_methods_IGF1R_xphen.csv", row.names = FALSE)
write.csv(heterogeneity_df, "results/Table3/MR_heterogeneity_IGF1R_xphen.csv", row.names = FALSE)
write.csv(pleiotropy_df, "results/Table3/MR_pleiotropy_IGF1R_xphen.csv", row.names = FALSE)
write.csv(steiger_df, "results/Table3/MR_steiger_IGF1R_xphen.csv", row.names = FALSE)
write.csv(table3, "results/Table3/Table3_condensed.csv", row.names = FALSE)

write_xlsx(list(
    Table3_Main        = table3,
    Full_MR_methods    = results_df,
    Heterogeneity      = heterogeneity_df,
    Pleiotropy         = pleiotropy_df,
    Steiger            = steiger_df
), "results/Table3/Table3_IGF1R_xphen.xlsx")

cat("\n✅ Saved to results/Table3/ — XLSX + CSVs\n")
cat("\n=== DONE", as.character(Sys.time()), "===\n")
sink()
