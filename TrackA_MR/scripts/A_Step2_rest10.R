library(TwoSampleMR)

genes <- list(
    IL6 = "eqtl-a-ENSG00000136244", TGFB1 = "eqtl-a-ENSG00000105329",
    HAS2 = "eqtl-a-ENSG00000153446", PPARG = "eqtl-a-ENSG00000132170",
    ARRB1 = "eqtl-a-ENSG00000137486", IGF1 = "eqtl-a-ENSG00000017427",
    IRS1 = "eqtl-a-ENSG00000169047", AKT1 = "eqtl-a-ENSG00000142208",
    TNF = "eqtl-a-ENSG00000232810", CTLA4 = "eqtl-a-ENSG00000163599"
)

out_pri <- "ebi-a-GCST90018627"
out_rep <- "ebi-a-GCST90038636"

results_df <- data.frame()

for (g in names(genes)) {
    cat("Running", g, "...\n")
    exp <- tryCatch(extract_instruments(outcomes = genes[[g]], p1 = 5e-6, clump = TRUE), error = function(e) NULL)
    if (is.null(exp) || nrow(exp) == 0) next
    cat(" Extracted IVs\n")
    out1 <- tryCatch(extract_outcome_data(exp$SNP, out_pri), error = function(e) NULL)
    if (!is.null(out1) && nrow(out1) > 0) {
        h1 <- harmonise_data(exp, out1)
        if (sum(h1$mr_keep) > 0) {
            mr1 <- tryCatch(mr(h1), error = function(e) NULL)
            if (!is.null(mr1)) {
                ivw1 <- mr1[mr1$method == "Inverse variance weighted", ]
                if (nrow(ivw1) > 0) results_df <- rbind(results_df, data.frame(Gene = g, GWAS = "Primary(Graves)", b = ivw1$b, p = ivw1$pval, nsnp = nrow(h1[h1$mr_keep, ])))
            }
        }
    }

    out2 <- tryCatch(extract_outcome_data(exp$SNP, out_rep), error = function(e) NULL)
    if (!is.null(out2) && nrow(out2) > 0) {
        h2 <- harmonise_data(exp, out2)
        if (sum(h2$mr_keep) > 0) {
            mr2 <- tryCatch(mr(h2), error = function(e) NULL)
            if (!is.null(mr2)) {
                ivw2 <- mr2[mr2$method == "Inverse variance weighted", ]
                if (nrow(ivw2) > 0) results_df <- rbind(results_df, data.frame(Gene = g, GWAS = "Replication(Hyper)", b = ivw2$b, p = ivw2$pval, nsnp = nrow(h2[h2$mr_keep, ])))
            }
        }
    }
}
write.csv(results_df, "c:/ProjectTEDGWAS/TrackA_MR/results/MR_rest10_summary.csv", row.names = F)
cat("\nDONE\n")
