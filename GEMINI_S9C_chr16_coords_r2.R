# =============================================================================
# TED-TRAP v5 — Stage 9C : chr16p11.2 trio — bp coords + EUR pairwise r²
# -----------------------------------------------------------------------------
# Last round picked the WRONG reference (EAS chr14) for chr16 SNPs, so PLINK
# produced no .ld. Available references (from your own output) include:
#   c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/g1000_eur(.bim/.bed/.fam)
# That is the genome-wide EUR panel and DOES contain chr16. Use it.
#
# Trio: rs4889606 (HSD3B7), rs34649473 (VKORC1), rs78924645 (PRSS36)
# Goal: (a) hg19 bp position of each, (b) pairwise r² in 1000G EUR.
# DO NOT edit the manuscript. Paste full console output to Claude.
# =============================================================================

trio <- c("rs4889606", "rs34649473", "rs78924645")

# ---- (a) bp positions: read straight from the EUR .bim (cols: chr, rsid, cM, bp, a1, a2)
eur_bim <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/g1000_eur.bim"
cat("== bp positions from EUR .bim ==\n")
if (file.exists(eur_bim)) {
  bim <- read.table(eur_bim, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("chr","rsid","cM","bp","a1","a2"))
  sub <- bim[bim$rsid %in% trio, c("rsid","chr","bp","a1","a2")]
  print(sub)
  if (nrow(sub) > 0) {
    span <- max(sub$bp) - min(sub$bp)
    cat(sprintf("[span across trio]: %d bp (%.1f kb) on chr %s\n",
                span, span/1000, paste(unique(sub$chr), collapse=",")))
  }
} else {
  cat("!! EUR .bim not found at expected path; please list",
      "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/ and report.\n")
}

# ---- (b) pairwise r² in EUR (force the genome-wide EUR panel)
eur_ref <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/g1000_eur"
cat("\n== pairwise r² in 1000G EUR ==\n")
if (all(file.exists(paste0(eur_ref, c(".bim",".bed",".fam"))))) {
  tmp <- tempfile(fileext = ".txt"); writeLines(trio, tmp)
  outp <- tempfile()
  # locate a plink executable
  plink_candidates <- c("plink",
    "c:/ProjectTEDGWAS/tools/plink/plink.exe",
    "c:/ProjectTEDGWAS/tools/plink_win64/plink.exe",
    list.files("c:/ProjectTEDGWAS/tools", pattern="^plink\\.exe$",
                recursive=TRUE, full.names=TRUE))
  plink_candidates <- unique(plink_candidates)
  cat("[plink candidates]:\n"); print(plink_candidates)
  ran <- FALSE
  for (pk in plink_candidates) {
    cmd <- sprintf('"%s" --bfile "%s" --extract "%s" --r2 square --out "%s"',
                   pk, eur_ref, tmp, outp)
    cat("[try]:", cmd, "\n")
    st <- suppressWarnings(try(system(cmd), silent = TRUE))
    if (file.exists(paste0(outp, ".ld"))) { ran <- TRUE; break }
    # also try square0 / table forms if needed
    cmd2 <- sprintf('"%s" --bfile "%s" --extract "%s" --r2 inter-chr --ld-window-r2 0 --out "%s"',
                    pk, eur_ref, tmp, outp)
    st2 <- suppressWarnings(try(system(cmd2), silent = TRUE))
    if (file.exists(paste0(outp, ".ld"))) { ran <- TRUE; break }
  }
  ld_file <- paste0(outp, ".ld")
  if (file.exists(ld_file)) {
    cat("\n[pairwise r² (EUR)]:\n")
    print(read.table(ld_file, header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE))
  } else {
    cat("(still no .ld — please paste any PLINK stderr, and confirm a working",
        "plink.exe path under c:/ProjectTEDGWAS/tools/)\n")
  }
} else {
  cat("!! EUR genome-wide reference (g1000_eur.bed/.bim/.fam) not all present.\n",
      "   Files present:\n")
  print(file.exists(paste0(eur_ref, c(".bim",".bed",".fam"))))
}

cat("\n==== DONE — paste FULL output to Claude. ====\n")
