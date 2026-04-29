setwd("c:/ProjectTEDGWAS/TrackA_MR/GEO_validation/")
library(GEOquery)
library(limma)
library(data.table)

eset <- readRDS("GSE105149_eset.rds")
metadata <- readRDS("GSE105149_metadata.rds")

metadata$group <- NA
for (i in 1:nrow(metadata)) {
    row_text <- tolower(paste(metadata[i, ], collapse = " "))
    if (grepl("diagnosis: ted", row_text) || grepl("active ted", row_text)) {
        metadata$group[i] <- "TED"
    } else if (grepl("diagnosis: normal", row_text) || grepl("normal|control", row_text)) {
        metadata$group[i] <- "Control"
    }
}
keep_idx <- which(!is.na(metadata$group) & (metadata$group %in% c("TED", "Control")))
eset <- eset[, keep_idx]
metadata <- metadata[keep_idx, ]

exprs_mat <- exprs(eset)
if (max(exprs_mat, na.rm = TRUE) > 50) {
    exprs_mat <- log2(exprs_mat + 1)
}

group <- factor(metadata$group, levels = c("Control", "TED"))
design <- model.matrix(~group)

fit <- lmFit(exprs_mat, design)
fit <- eBayes(fit)

fdata <- fData(eset)
res <- topTable(fit, coef = "groupTED", number = Inf)
res$id_raw <- rownames(res)

if ("Gene symbol" %in% names(fdata)) {
    raw_sym <- as.character(fdata[rownames(res), "Gene symbol"])
    res$hgnc_symbol <- sapply(strsplit(raw_sym, "///"), function(x) trimws(x[1]))
}

fwrite(res, "GSE105149_limma_full_results.csv")

candidate_genes <- list(
    MR_Primary = c("TSHR", "IGF1R", "IGF1"),
    MR_Secondary = c("ARRB1", "PPARG", "IRS1", "AKT1", "TNF", "CTLA4"),
    Insulin_Cassette = c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1"),
    TED_Biology = c("HAS1", "HAS2", "HAS3", "ADIPOQ", "FABP4", "CEBPA")
)
all_candidates <- unlist(candidate_genes)

res_cand <- res[toupper(res$hgnc_symbol) %in% all_candidates | toupper(res$id_raw) %in% all_candidates, ]
res_cand$hgnc_symbol <- ifelse(toupper(res_cand$hgnc_symbol) %in% all_candidates, toupper(res_cand$hgnc_symbol), toupper(res_cand$id_raw))

# Aggregate by gene (lowest p-value)
res_cand <- res_cand[!is.na(res_cand$hgnc_symbol) & res_cand$hgnc_symbol != "", ]
res_cand <- res_cand[!is.na(res_cand$P.Value), ]
setDT(res_cand)
res_cand <- res_cand[order(P.Value)]
res_cand <- res_cand[, .SD[1], by = hgnc_symbol]

res_cand$padj_within20 <- p.adjust(res_cand$P.Value, method = "BH")
res_cand$module <- NA
for (mod in names(candidate_genes)) {
    res_cand$module[res_cand$hgnc_symbol %in% candidate_genes[[mod]]] <- mod
}

ins_res <- res_cand[module == "Insulin_Cassette"]
if (nrow(ins_res) > 0) {
    n_up <- sum(ins_res$logFC > 0)
    binom_p <- binom.test(n_up, nrow(ins_res), p = 0.5, alternative = "greater")$p.value
    cat("\nInsulin cassette (GSE105149):", n_up, "/", nrow(ins_res), "up, Binomial P:", binom_p, "\n")
} else {
    cat("\nInsulin cassette: 0 genes found.\n")
}

fwrite(res_cand, "GSE105149_candidate_genes_limma.csv")

inhouse <- data.table(
    hgnc_symbol = c(
        "TSHR", "INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1",
        "HAS1", "HAS2", "HAS3", "PPARG", "CEBPA", "ADIPOQ", "FABP4",
        "IL6", "TNF", "IGF1R", "IGF1", "IRS1", "AKT1", "ARRB1", "CTLA4"
    ),
    inhouse_logFC = c(
        2.33, 0.54, 0.65, 1.02, 0.44, 0.38,
        NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
    )
)

comp <- merge(inhouse, res_cand, by = "hgnc_symbol", all.x = TRUE)
comp$direction_match <- sign(comp$inhouse_logFC) == sign(comp$logFC)
fwrite(comp, "Cross_validation_inhouse_vs_GSE105149.csv")

cat("\n=== Final candidate summary ===\n")
print(comp[!is.na(logFC), .(hgnc_symbol, logFC, P.Value, padj_within20, inhouse_logFC, direction_match)])
