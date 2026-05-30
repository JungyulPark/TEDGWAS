# =============================================================================
# TED-TRAP v5  —  Table S8 external orbital-tissue direction replication
# STAGE 1 of 2 : METADATA + GROUPING + PMID AUDIT  (NO differential expression yet)
# -----------------------------------------------------------------------------
# Purpose : Before re-computing any log2FC / P value, audit the two external GEO
#           datasets so that (a) the TED-vs-control grouping is defined by the
#           data, not assumed, and (b) the PubMed ID for citation is read from
#           the GEO record itself.  Memory warns: DO NOT cite GSE105149 PMID
#           28910444 from memory — it must be confirmed here.
# Datasets: GSE58331  (expected: Rosenbaum 2015, PLoS One; GPL570)
#           GSE105149 (expected: Rosenbaum 2017, JAMA Ophthalmol; GPL570)
# Output  : prints sample annotation + platform + PubMedIds for both series.
# ACTION  : paste the FULL console output back to Claude. Claude will define the
#           exact grouping and send STAGE 2 (limma) for TSHR/IGF1R/CTLA4/HAS1.
# =============================================================================

# install.packages("BiocManager"); BiocManager::install(c("GEOquery","limma"))
suppressMessages({
  library(GEOquery)
  library(limma)
})

audit_geo <- function(gse_id){
  cat("\n#####################################################################\n")
  cat("##  ", gse_id, "\n")
  cat("#####################################################################\n")

  gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)

  # if multiple platforms, keep the GPL570 ExpressionSet(s); report all annotations
  anns <- sapply(gset, annotation)
  cat("\n[Platforms present in series]:\n"); print(anns)
  keep <- which(grepl("GPL570", anns))
  if (length(keep) == 0) keep <- seq_along(gset)   # fall back to all; flag in output

  for (i in keep){
    eset <- gset[[i]]
    cat("\n----- ExpressionSet:", names(gset)[i], " | platform:", annotation(eset),
        " | dim:", paste(dim(eset), collapse=" x "), "-----\n")

    pd <- pData(eset)
    info_cols <- grep("title|characteristic|source|disease|diagnos|group|status|tissue|type",
                      colnames(pd), ignore.case = TRUE)
    if (length(info_cols) == 0) info_cols <- seq_len(min(ncol(pd), 6))
    cat("\n[Sample annotation — define TED vs control from these]:\n")
    print(pd[, info_cols, drop = FALSE], right = FALSE)

    cat("\n[Expression value range — confirm log2-scale, not raw]:\n")
    ex <- exprs(eset)
    print(summary(as.vector(ex[, 1:min(3, ncol(ex))])))
    cat("  max overall:", round(max(ex, na.rm=TRUE), 2),
        " | negative values present:", any(ex < 0, na.rm=TRUE), "\n")

    cat("\n[Backbone-gene probes present on this platform]:\n")
    fd <- fData(eset)
    sym_col <- grep("Gene.?symbol|^Symbol$|GENE_SYMBOL", colnames(fd), ignore.case=TRUE)[1]
    if (!is.na(sym_col)){
      for (g in c("TSHR","IGF1R","CTLA4","HAS1")){
        hit <- which(toupper(trimws(fd[[sym_col]])) == g |
                     grepl(paste0("(^|[^A-Z])", g, "([^A-Z]|$)"), toupper(fd[[sym_col]])))
        cat(sprintf("  %-6s : %d probe(s)  [%s]\n", g, length(hit),
                    paste(head(rownames(fd)[hit], 5), collapse=", ")))
      }
    } else {
      cat("  (no gene-symbol column found in fData; STAGE 2 will map via GPL570 annotation)\n")
    }
  }

  cat("\n[Associated PubMedIds from GEO record — VERIFY against the GEO web page]:\n")
  pmid <- experimentData(gset[[keep[1]]])@pubMedIds
  print(pmid)
  cat("[Series title]:", experimentData(gset[[keep[1]]])@title, "\n")
}

audit_geo("GSE58331")
audit_geo("GSE105149")

cat("\n\n==== STAGE 1 COMPLETE — paste everything above back to Claude. ====\n")
cat("Do NOT run any limma / topTable yet; grouping must be confirmed first.\n")
