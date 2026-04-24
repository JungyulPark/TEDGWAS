library(data.table)
library(coloc)
FG_DIR <- "c:/ProjectTEDGWAS"
setwd(FG_DIR)

dt <- fread("finngen_R12_GRAVES_OPHT.gz", select = c("#chrom", "pos", "ref", "alt", "rsids", "beta", "sebeta", "pval", "af_alt"))
setnames(dt, "#chrom", "chrom")

fg_clean <- dt[chrom == 14 & pos >= 80921085 & pos <= 82081779]
fg_clean[, rsid := sub(",.*", "", rsids)]
fg_clean <- fg_clean[!is.na(beta) & !is.na(sebeta) & sebeta > 0 & !is.na(rsid) & rsid != ""]
fg_clean <- as.data.frame(fg_clean[!duplicated(rsid)])

cat("Fetching eQTL associations...\n")
e_loc <- ieugwasr::associations(variants = "14:80921085-82081779", id = "eqtl-a-ENSG00000165409", opengwas_jwt = Sys.getenv("OPENGWAS_JWT"))
e_clean <- as.data.frame(e_loc[!is.na(e_loc$beta) & !is.na(e_loc$se) & e_loc$se > 0, ])
e_clean <- e_clean[!duplicated(e_clean$rsid), ]

comm <- intersect(e_clean$rsid, fg_clean$rsid)
cat("Common SNPs:", length(comm), "\n")

e_sub <- e_clean[match(comm, e_clean$rsid), ]
f_sub <- fg_clean[match(comm, fg_clean$rsid), ]

me <- as.numeric(e_sub$eaf)
me[is.na(me) | me <= 0 | me >= 0.5] <- 0.3
mf <- as.numeric(f_sub$af_alt)
mf[is.na(mf) | mf <= 0 | mf >= 0.5] <- 0.3

d1 <- list(
    snp = comm,
    beta = as.numeric(e_sub$beta),
    varbeta = as.numeric(e_sub$se)^2,
    position = as.integer(e_sub$position),
    type = "quant",
    N = 31684L,
    MAF = me
)

d2 <- list(
    snp = comm,
    beta = as.numeric(f_sub$beta),
    varbeta = as.numeric(f_sub$sebeta)^2,
    position = as.integer(f_sub$pos),
    type = "cc",
    N = 500348L,
    s = 858 / 500348,
    MAF = mf
)

res <- coloc.abf(dataset1 = d1, dataset2 = d2)
print(res$summary)

top <- which.max(res$results$SNP.PP.H4)
cat(sprintf("\nTop shared SNP: %s (PP.H4 = %.4f)\n", res$results$snp[top], res$results$SNP.PP.H4[top]))
