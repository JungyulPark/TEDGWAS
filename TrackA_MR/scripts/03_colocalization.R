# =============================================================================
# TED-TRAP Step 3 — Colocalization at TSHR locus
# Purpose:
#   Resolve single-IV limitation of TSHR by showing shared causal variant
#   between eQTLGen blood TSHR expression and Graves disease GWAS.
# =============================================================================

required_pkgs <- c("coloc", "ieugwasr", "dplyr", "ggplot2")
for (p in required_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p, repos = "http://cran.us.r-project.org")
    }
}
library(coloc)
library(ieugwasr)
library(dplyr)
library(ggplot2)

token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(token) < 100) stop("Token not loaded.")

for (d in c("TrackA_MR/results", "TrackA_MR/figures", "TrackA_MR/logs")) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
log_file <- "TrackA_MR/logs/03_colocalization.log"
sink(log_file, split = TRUE)

cat("=== Step 3: Colocalization at TSHR locus ===\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

LOCUS_CHR <- 14
LOCUS_START <- 80921085
LOCUS_END <- 82081779

cat(sprintf(
    "Locus: chr%d : %s - %s (%d kb window)\n",
    LOCUS_CHR, format(LOCUS_START, big.mark = ","),
    format(LOCUS_END, big.mark = ","),
    round((LOCUS_END - LOCUS_START) / 1000)
))

gwas_ids <- list(
    eQTL_TSHR       = "eqtl-a-ENSG00000165409",
    GD_Primary_BBJ  = "ebi-a-GCST90018627",
    GD_Repl_UKB     = "ebi-a-GCST90038636",
    GD_Sens_FinnGen = "finn-b-E4_GRAVES_OPHT"
)

sample_info <- list(
    eQTL_TSHR       = list(N = 31684, type = "quant", s = NA),
    GD_Primary_BBJ  = list(N = 212453, type = "cc", s = 2809 / 212453),
    GD_Repl_UKB     = list(N = 484598, type = "cc", s = 3731 / 484598),
    GD_Sens_FinnGen = list(N = 218792, type = "cc", s = 1828 / 218792)
)

fetch_region <- function(id, chr, start, end, label) {
    cat(sprintf("\n  Fetching: %s [%s]\n", label, id))
    dat <- tryCatch(
        {
            associations(
                variants = sprintf("%d:%d-%d", chr, start, end),
                id = id,
                opengwas_jwt = token
            )
        },
        error = function(e) {
            NULL
        }
    )

    if (is.null(dat) || nrow(dat) == 0) {
        Sys.sleep(5)
        dat <- tryCatch(
            {
                associations(
                    variants = sprintf("%d:%d-%d", chr, start, end),
                    id = id,
                    opengwas_jwt = token
                )
            },
            error = function(e) {
                NULL
            }
        )

        if (is.null(dat) || nrow(dat) == 0) {
            cat(sprintf("    ⚠️ No data returned\n"))
            return(NULL)
        }
    }
    cat(sprintf("    ✅ %d SNPs retrieved\n", nrow(dat)))
    dat
}

cat("\n=== Fetching locus data via OpenGWAS API ===\n")
data_list <- list()
for (nm in names(gwas_ids)) {
    data_list[[nm]] <- fetch_region(gwas_ids[[nm]], LOCUS_CHR, LOCUS_START, LOCUS_END, nm)
    Sys.sleep(1)
}

prepare_coloc_input <- function(df, info, dataset_label) {
    df <- df[!is.na(df$beta) & !is.na(df$se) & df$se > 0, ]
    df <- df[!duplicated(df$rsid), ]
    if (!"eaf" %in% names(df)) df$eaf <- NA
    df$eaf[is.na(df$eaf)] <- 0.3

    input <- list(
        snp = df$rsid,
        beta = df$beta,
        varbeta = df$se^2,
        MAF = pmin(df$eaf, 1 - df$eaf),
        position = df$position,
        N = info$N,
        type = info$type
    )
    if (info$type == "cc") input[["s"]] <- info$s
    cat(sprintf("  %-20s: %d SNPs prepared\n", dataset_label, length(input$snp)))
    input
}

cat("\n=== Preparing coloc inputs ===\n")
if (!is.null(data_list$eQTL_TSHR)) {
    eqtl_input <- prepare_coloc_input(data_list$eQTL_TSHR, sample_info$eQTL_TSHR, "eQTL (TSHR)")
} else {
    eqtl_input <- NULL
}

gwas_inputs <- list()
for (nm in c("GD_Primary_BBJ", "GD_Repl_UKB", "GD_Sens_FinnGen")) {
    if (!is.null(data_list[[nm]])) {
        gwas_inputs[[nm]] <- prepare_coloc_input(data_list[[nm]], sample_info[[nm]], nm)
    }
}

if (is.null(eqtl_input)) stop("eQTL data failed to load. Aborting sumstats run.")

cat("\n=== Running coloc.abf ===\n\n")
coloc_results <- list()

for (nm in names(gwas_inputs)) {
    cat(sprintf("\n%s\nColoc: eQTL TSHR × %s\n%s\n", paste(rep("-", 50), collapse = ""), nm, paste(rep("-", 50), collapse = "")))
    common_snps <- intersect(eqtl_input$snp, gwas_inputs[[nm]]$snp)
    cat(sprintf("Common SNPs: %d\n", length(common_snps)))

    if (length(common_snps) < 50) {
        cat("  ⚠️ <50 common SNPs — coloc unreliable\n")
        next
    }

    subset_dataset <- function(d, snps) {
        idx <- match(snps, d$snp)
        res <- list(snp = snps, beta = d$beta[idx], varbeta = d$varbeta[idx], MAF = d$MAF[idx], position = d$position[idx], N = d$N, type = d$type)
        if ("s" %in% names(d)) res[["s"]] <- d[["s"]]
        res
    }

    e_sub <- subset_dataset(eqtl_input, common_snps)
    g_sub <- subset_dataset(gwas_inputs[[nm]], common_snps)

    result <- tryCatch(coloc.abf(dataset1 = e_sub, dataset2 = g_sub), error = function(e) {
        cat(sprintf("  [ERROR] %s\n", e$message))
        NULL
    })
    if (is.null(result)) next

    pp <- result$summary
    cat(sprintf("\n  PP.H0 (no causal):           %.4f\n  PP.H1 (eQTL only):           %.4f\n  PP.H2 (GWAS only):           %.4f\n  PP.H3 (different causals):   %.4f\n  PP.H4 (SHARED causal):       %.4f  %s\n", pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"], pp["PP.H3.abf"], pp["PP.H4.abf"], ifelse(pp["PP.H4.abf"] > 0.95, "⭐⭐⭐ VERY STRONG", ifelse(pp["PP.H4.abf"] > 0.80, "⭐⭐ STRONG", ifelse(pp["PP.H4.abf"] > 0.50, "⭐ SUGGESTIVE", "❌ NO EVIDENCE")))))

    top_snp_idx <- which.max(result$results$SNP.PP.H4)
    cat(sprintf("\n  Top colocalized SNP: %s (SNP-level PP.H4 = %.3f)\n", result$results$snp[top_snp_idx], result$results$SNP.PP.H4[top_snp_idx]))

    coloc_results[[nm]] <- list(summary = pp, top_snp = result$results$snp[top_snp_idx], n_common_snps = length(common_snps), full = result)
}

cat("\n\n", rep("=", 70), "\nFINAL COLOCALIZATION TABLE (Table 2 Panel B)\n", rep("=", 70), "\n\n", sep = "")
coloc_table <- data.frame()
for (nm in names(coloc_results)) {
    r <- coloc_results[[nm]]
    pp <- r$summary
    coloc_table <- rbind(coloc_table, data.frame(Outcome = nm, N_SNPs = r$n_common_snps, PP_H0 = round(pp["PP.H0.abf"], 4), PP_H1 = round(pp["PP.H1.abf"], 4), PP_H2 = round(pp["PP.H2.abf"], 4), PP_H3 = round(pp["PP.H3.abf"], 4), PP_H4 = round(pp["PP.H4.abf"], 4), Top_SNP = r$top_snp, Interpretation = ifelse(pp["PP.H4.abf"] > 0.95, "Very Strong", ifelse(pp["PP.H4.abf"] > 0.80, "Strong", ifelse(pp["PP.H4.abf"] > 0.50, "Suggestive", "No Evidence"))), stringsAsFactors = FALSE))
}
if (nrow(coloc_table) > 0) print(coloc_table, row.names = FALSE)
write.csv(coloc_table, "TrackA_MR/results/03_coloc_TSHR.csv", row.names = FALSE)
saveRDS(coloc_results, "TrackA_MR/results/03_coloc_full_results.rds")
cat("\n💾 Saved: TrackA_MR/results/03_coloc_TSHR.csv\n💾 Saved: TrackA_MR/results/03_coloc_full_results.rds\n")

cat("\n\n=== Prior Sensitivity Analysis (Wallace 2020) ===\nTesting robustness of PP.H4 > 0.8 to different prior choices.\n\n")
for (nm in names(coloc_results)) {
    cat(sprintf("\n-- %s --\n", nm))
    fig_path <- sprintf("TrackA_MR/figures/03_sensitivity_%s.pdf", nm)
    pdf(fig_path, width = 9, height = 6)
    tryCatch(
        {
            sensitivity(coloc_results[[nm]]$full, "H4 > 0.8")
            cat(sprintf("  ✅ Sensitivity plot saved: %s\n", fig_path))
        },
        error = function(e) {
            cat(sprintf("  [WARN] Sensitivity plot failed: %s\n", e$message))
        }
    )
    dev.off()
}

cat("\n\n", rep("=", 70), "\nSCIENTIFIC INTERPRETATION\n", rep("=", 70), "\n\n", sep = "")
for (nm in names(coloc_results)) {
    pph4 <- coloc_results[[nm]]$summary["PP.H4.abf"]
    pph3 <- coloc_results[[nm]]$summary["PP.H3.abf"]
    if (pph4 > 0.80) {
        cat(sprintf("✅ %s: eQTL and GWAS share a SINGLE causal variant (PP.H4 = %.3f).\n   This overcomes the single-IV limitation of Wald ratio MR.\n", nm, pph4))
    } else if (pph3 > 0.50) {
        cat(sprintf("⚠️  %s: Evidence for DIFFERENT causal variants (PP.H3 = %.3f).\n   MR signal may reflect locus-specific pleiotropy, not direct causality.\n", nm, pph3))
    } else {
        cat(sprintf("❓ %s: Inconclusive (PP.H4 = %.3f, PP.H3 = %.3f).\n", nm, pph4, pph3))
    }
    cat("\n")
}
sink()
cat("\n✅ Step 3 complete.\n   Log: TrackA_MR/logs/03_colocalization.log\n   Main result: TrackA_MR/results/03_coloc_TSHR.csv\n   Plots: TrackA_MR/figures/03_sensitivity_*.pdf\n\n")
