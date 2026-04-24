# === eQTLGen에서 TED target gene eQTL 추출 ===

library(data.table)

# eQTLGen 파일 경로
eqtl_file <- "c:/Projectbulid/CP3/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"


# Target genes
target_genes <- data.frame(
    gene = c(
        "IGF1R", "IGF1", "IRS1", "AKT1", "TSHR", "ARRB1",
        "IL6", "TNF", "HAS2", "PPARG", "TGFB1", "CTLA4"
    ),
    ensembl = c(
        "ENSG00000140443", "ENSG00000017427", "ENSG00000169047",
        "ENSG00000142208", "ENSG00000165409", "ENSG00000137486",
        "ENSG00000136244", "ENSG00000232810", "ENSG00000153446",
        "ENSG00000132170", "ENSG00000105329", "ENSG00000163599"
    ),
    module = c(
        rep("IGF1R_pathway", 4), rep("TSHR_pathway", 2),
        rep("Inflammatory", 2), "HA_synthesis", "Adipogenesis",
        "Fibrosis", "Immune"
    )
)

# 읽기 (대용량이므로 fread 사용)
cat("Reading eQTLGen data...\n")
eqtl <- fread(eqtl_file)
cat(sprintf("Total rows: %d\n", nrow(eqtl)))

# 각 유전자별 추출
for (i in 1:nrow(target_genes)) {
    gene_name <- target_genes$gene[i]
    ensembl_id <- target_genes$ensembl[i]

    gene_eqtl <- eqtl[Gene == ensembl_id | GeneSymbol == gene_name]

    cat(sprintf("\n%s (%s): %d cis-eQTLs total\n", gene_name, ensembl_id, nrow(gene_eqtl)))

    # Bonferroni significant만
    sig <- gene_eqtl[BonferroniP < 0.05]
    cat(sprintf("  Bonferroni significant: %d\n", nrow(sig)))

    # 저장
    outfile <- sprintf("c:/ProjectTEDGWAS/TrackA_MR/data/eqtl_%s_all.csv", gene_name)
    fwrite(gene_eqtl, outfile)

    outfile_sig <- sprintf("c:/ProjectTEDGWAS/TrackA_MR/data/eqtl_%s_sig.csv", gene_name)
    fwrite(sig, outfile_sig)
}

cat("\n=== Summary ===\n")
cat("Check: Do TSHR and IGF1R have sufficient IVs (>=3)?\n")
cat("If any gene has 0 significant eQTLs, exclude from MR\n")
