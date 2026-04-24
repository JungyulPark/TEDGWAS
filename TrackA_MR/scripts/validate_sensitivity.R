library(TwoSampleMR)
library(data.table)

genes <- list(TSHR = "eqtl-a-ENSG00000165409", IGF1R = "eqtl-a-ENSG00000140443", TNF = "eqtl-a-ENSG00000232810")
out_pri <- "ebi-a-GCST90018627"
out_rep <- "ebi-a-GCST90038636"

cat("=== 2-2. Sensitivity Analysis ===\n")
run_sens <- function(g_name, out_id, label) {
    exp <- tryCatch(extract_instruments(outcomes = genes[[g_name]], p1 = 5e-6, clump = TRUE), error = function(e) NULL)
    if (is.null(exp) || nrow(exp) == 0) {
        return(NULL)
    }
    out <- tryCatch(extract_outcome_data(exp$SNP, out_id), error = function(e) NULL)
    if (is.null(out) || nrow(out) == 0) {
        return(NULL)
    }
    h <- harmonise_data(exp, out)

    cat(sprintf("\n--- %s (%s, nSNP=%d) ---\n", g_name, label, sum(h$mr_keep)))
    if (sum(h$mr_keep) >= 3) {
        pleio <- tryCatch(mr_pleiotropy_test(h), error = function(e) NULL)
        if (!is.null(pleio) && nrow(pleio) > 0) cat(sprintf("Egger intercept P: %.4f\n", pleio$pval))

        het <- tryCatch(mr_heterogeneity(h), error = function(e) NULL)
        if (!is.null(het) && nrow(het) > 0) {
            ivw_het <- het[het$method == "Inverse variance weighted", ]
            if (nrow(ivw_het) > 0) cat(sprintf("Cochran Q P (IVW): %.4f\n", ivw_het$Q_pval))
        }
    } else {
        cat("nSNP < 3, cannot run pleiotropy/heterogeneity.\n")
    }

    steiger <- tryCatch(directionality_test(h), error = function(e) NULL)
    if (!is.null(steiger) && nrow(steiger) > 0) {
        cat(sprintf("Steiger correct direction: %s, P: %.2e\n", steiger$correct_causal_direction, steiger$steiger_pval))
    }

    if (g_name == "TSHR" && label == "Primary") {
        cat("\n=== 2-3. TSHR IV 상세 ===\n")
        print(h[h$mr_keep, c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "beta.exposure", "pval.exposure", "beta.outcome", "pval.outcome")])
    }
}

run_sens("TSHR", out_pri, "Primary")
run_sens("IGF1R", out_pri, "Primary")
run_sens("TNF", out_rep, "Replication")

cat("\n=== 3. Coloc 상세 검증 ===\n")
res <- readRDS("c:/ProjectTEDGWAS/TrackA_MR/results/coloc_finngen.rds")
print(res$summary)
cat("Type of eQTL dataset: quant\n")
cat("Type of GWAS dataset: cc\n")
cat("s (case proportion): 0.016\n")

gwas <- fread("c:/ProjectTEDGWAS/TrackA_MR/data/finngen_TSHR_region.csv")
eqtl <- fread("c:/ProjectTEDGWAS/TrackA_MR/data/eqtl_TSHR_all.csv")
shared <- intersect(gwas$rsids, eqtl$SNP)
g <- gwas[rsids %in% shared]
e <- eqtl[SNP %in% shared]
top_g <- g[which.min(pval)]$rsids
top_e <- e[which.min(Pvalue)]$SNP
cat(sprintf("Top GWAS SNP: %s, Top eQTL SNP: %s (Match? %s)\n", top_g, top_e, top_g == top_e))
