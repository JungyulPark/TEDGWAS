# =============================================================================
# TED-TRAP v5 — Table S8 external direction replication (BACKBONE GENES ONLY)
# STAGE 2 (CORRECTED): GSE58331 + GSE105149, STRICT TED-vs-Normal grouping
# -----------------------------------------------------------------------------
# SCOPE — backbone genes ONLY: TSHR, IGF1R, CTLA4.
#   DO NOT include insulin-cassette genes (INSR/IRS2/FOXO1/PIK3R1/PDPK1).
#   Those belong to Track C, NOT this v5 backbone manuscript.
#
# GROUPING — TED vs NORMAL/HEALTHY CONTROL ONLY.
#   Both datasets are multi-disease orbital/lacrimal collections, so a naive
#   "non-TED = control" grouping is WRONG. EXCLUDE sarcoidosis, GPA (Wegener),
#   nonspecific orbital inflammation (NSOI), lymphoma, and every non-TED disease.
#
# Confirmed citations (verified by Claude — use these, do not re-guess):
#   GSE58331  = Rosenbaum 2015 PLoS ONE 10(9):e0137654           PMID 26371757
#   GSE105149 = Rosenbaum 2017 JAMA Ophthalmol 135(11):1156-1162 PMID 28975236
#
# OUTPUT — for each dataset: explicit group counts + per-backbone-gene
#   log2FC, P, adjusted P (collapse multi-probe by highest mean expression).
# ACTION — paste FULL console output back to Claude; do not edit the manuscript.
# =============================================================================

suppressMessages({ library(GEOquery); library(limma) })

BACKBONE <- c("TSHR", "IGF1R", "CTLA4")

run_ds <- function(gse_id) {
  cat("\n############################## ", gse_id, " ##############################\n")
  g <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
  k <- which(grepl("GPL570", sapply(g, annotation))); if (!length(k)) k <- 1
  eset <- g[[k[1]]]
  pd <- pData(eset)

  # 1) locate diagnosis/disease column and PRINT raw values (auditable)
  dcol <- grep("diagnosis|disease|group|status|condition", colnames(pd),
               ignore.case = TRUE, value = TRUE)
  cat("[candidate group columns]:", paste(dcol, collapse = ", "), "\n")
  dx <- apply(pd[, dcol, drop=FALSE], 1, function(x) paste(x, collapse = " "))
  dx <- tolower(as.character(dx))
  cat("[raw combined diagnosis strings]:\n"); print(table(dx, useNA = "ifany"))

  # 2) STRICT grouping: TED vs Control only; everything else -> NA (excluded)
  grp <- rep(NA_character_, length(dx))
  grp[grepl("ted|thyroid|graves|orbitopathy|gravesdisease", dx)] <- "TED"
  grp[grepl("normal|healthy|control", dx)] <- "Control"
  # sarcoid / gpa / wegener / nsoi / pseudotumor / lymphoma -> remain NA (dropped)

  cat("\n[FINAL grouping — must be TED vs Control ONLY; all else excluded as NA]:\n")
  print(table(grp, useNA = "ifany"))

  keep <- !is.na(grp)
  eset2 <- eset[, keep]
  grp2  <- factor(grp[keep], levels = c("Control", "TED"))
  cat(sprintf("[n used] TED = %d ; Control = %d\n", sum(grp2 == "TED"), sum(grp2 == "Control")))

  # 3) limma DE: TED vs Control
  design <- model.matrix(~ grp2); colnames(design) <- c("Intercept", "TED_vs_Control")
  fit <- eBayes(lmFit(exprs(eset2), design))
  tt  <- topTable(fit, coef = "TED_vs_Control", number = Inf, adjust.method = "BH")
  tt$probe <- rownames(tt)

  fd <- fData(eset2)
  symcol <- grep("Gene.?symbol|^Symbol$|GENE_SYMBOL", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  tt$symbol <- fd[match(tt$probe, rownames(fd)), symcol]

  # 4) backbone genes — collapse multi-probe by highest mean expression
  ex_mean <- rowMeans(exprs(eset2))
  cat("\n[BACKBONE results — log2FC is TED vs Control; positive = UP in TED]:\n")
  for (gsym in BACKBONE) {
    rows <- tt[which(toupper(trimws(tt$symbol)) == gsym), , drop = FALSE]
    if (nrow(rows) == 0) { cat(sprintf("  %-6s : NO PROBE on platform\n", gsym)); next }
    rows$meanexpr <- ex_mean[rows$probe]
    best <- rows[which.max(rows$meanexpr), ]
    cat(sprintf("  %-6s | probe %-13s | log2FC=%+.3f | P=%.3g | adjP=%.3g  (n probes=%d)\n",
                gsym, best$probe, best$logFC, best$P.Value, best$adj.P.Val, nrow(rows)))
  }

  cat("\n[GEO record PubMedIds — should be 26371757 (GSE58331) / 28975236 (GSE105149)]:\n")
  print(experimentData(eset)@pubMedIds)
}

run_ds("GSE58331")
run_ds("GSE105149")
cat("\n\n==== DONE. Paste FULL output to Claude. Backbone only; no insulin cassette. ====\n")
