library(data.table)

setwd("c:/ProjectTEDGWAS/TrackA_MR/GEO_validation/")
cand <- fread("GSE105149_candidate_genes_limma.csv")

modules <- list(
    MR_Primary       = c("TSHR", "IGF1R", "IGF1"),
    MR_Secondary     = c("ARRB1", "PPARG", "IRS1", "AKT1", "TNF", "CTLA4"),
    Insulin_Cassette = c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1"),
    TED_Biology      = c("HAS1", "HAS2", "HAS3", "ADIPOQ", "FABP4", "CEBPA")
)

# Using hgnc_symbol not gene
cand$module <- NA
for (mod_name in names(modules)) {
    cand$module[cand$hgnc_symbol %in% modules[[mod_name]]] <- mod_name
}

sink("TaskB_Deliverable.txt")
cat("[TASK B: GSE105149 additional detail]\n")
cat("\n=== HAS1/2/3 (active TED expected up) ===\n")
has_res <- cand[hgnc_symbol %in% c("HAS1", "HAS2", "HAS3")]
for (g in c("HAS1", "HAS2", "HAS3")) {
    r <- has_res[hgnc_symbol == g]
    if (nrow(r) > 0) {
        cat(g, ": logFC =", round(r$logFC, 3), ", P =", signif(r$P.Value, 3), "\n")
    } else {
        cat(g, ": NA\n")
    }
}
sanity_state <- if (any(has_res$logFC > 0, na.rm = T)) "Passed (upregulated in active TED)" else "Failed"
cat("Sanity check passed/failed:", sanity_state, "\n")

cat("\n=== IGF1R ===\n")
igf1r <- cand[hgnc_symbol == "IGF1R"]
if (nrow(igf1r) > 0) {
    cat("logFC =", round(igf1r$logFC, 3), ", P =", signif(igf1r$P.Value, 3), "\n")
    dir_vs_inhouse <- ifelse(igf1r$logFC > 0, "Concordant (assuming both upregulated)", "Discordant")
    cat("Direction vs in-house:", dir_vs_inhouse, "\n")
} else {
    cat("logFC = NA, P = NA\nDirection vs in-house: NA\n")
}

cat("\n=== Adipogenic markers ===\n")
for (g in c("PPARG", "CEBPA", "ADIPOQ", "FABP4")) {
    r <- cand[hgnc_symbol == g]
    if (nrow(r) > 0) cat(g, ": logFC =", round(r$logFC, 3), ", P =", signif(r$P.Value, 3), "\n") else cat(g, ": NA\n")
}

cat("\n=== Inflammation markers ===\n")
for (g in c("IL6", "TNF")) {
    r <- cand[hgnc_symbol == g]
    if (nrow(r) > 0) cat(g, ": logFC =", round(r$logFC, 3), ", P =", signif(r$P.Value, 3), "\n") else cat(g, ": NA\n")
}

all_candidates <- unlist(modules)
detected <- cand$hgnc_symbol
missing <- setdiff(all_candidates, detected)
cat("\n=== Genes NOT detected in GSE105149 ===\n")
cat("Missing:", paste(missing, collapse = ", "), "\n")
cat("Total mapped:", length(detected), "/", length(all_candidates), "\n")
sink()
