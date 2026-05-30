#!/usr/bin/env Rscript
# =============================================================================
# TED-TRAP v5 — Figure 2: Druggable-genome-wide MR landscape (volcano)
# =============================================================================
# x = MR beta (BBJ discovery, primary method per gene), y = -log10(P)
# Bonferroni line at 1.965e-5. 2,545 genes. Highlight TSHR/IGF1R/CTLA4 + 13 hits.
# Reads directly from TaskD_03d — no hand numbers.
#
# Requires: ggplot2, data.table, ggrepel
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

WORK   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
MERGED <- file.path(WORK, "04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv")
SNP    <- file.path(WORK, "04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv")
OUTDIR <- file.path(WORK, "07_manuscript/figures")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

BONF <- 1.965e-5

mr <- fread(MERGED)
# primary method per gene for BBJ discovery
bbj <- mr[outcome=="BBJ_Graves" &
          ((n_iv==1 & method=="Wald ratio") |
           (n_iv>1  & method=="Inverse variance weighted"))]
bbj <- bbj[!is.na(beta) & !is.na(pvalue) & pvalue>0]
bbj[, neglog10P := -log10(pvalue)]

# gene positions for MHC/chr16 flagging
snp <- fread(SNP)
gpos <- snp[, .(chr=chr[1], pos=min(pos_hg19)), by=gene_symbol]
bbj <- merge(bbj, gpos, by="gene_symbol", all.x=TRUE)

# classify the 13 Bonferroni hits (verified)
known   <- c("TSHR","CTLA4")
mhc     <- function(c,p) (!is.na(c) & c==6 & p>=25e6 & p<=35e6)
c16     <- function(c,p) (!is.na(c) & c==16 & p>=30e6 & p<=32e6)
bbj[, class := "ns"]
bbj[pvalue < BONF, class := "other_hit"]
bbj[pvalue < BONF & mhc(chr,pos), class := "MHC"]
bbj[pvalue < BONF & c16(chr,pos), class := "chr16p11"]
bbj[pvalue < BONF & gene_symbol %in% known, class := "known_anchor"]
# label set: backbone always + all 13 hits
label_genes <- unique(c("TSHR","IGF1R","CTLA4",
                        bbj[pvalue < BONF, gene_symbol]))
bbj[, lab := ifelse(gene_symbol %in% label_genes, gene_symbol, NA_character_)]

# color palette
pal <- c(ns="grey80", other_hit="#888888", MHC="#9E9E9E",
         chr16p11="#C8A2C8", known_anchor="#C00000")
# IGF1R special (nominal, not Bonferroni) — mark distinctly
bbj[gene_symbol=="IGF1R", class := "IGF1R"]
pal["IGF1R"] <- "#2E75B6"

p <- ggplot(bbj, aes(beta, neglog10P)) +
  geom_point(aes(color=class), size=ifelse(bbj$class %in% c("ns","other_hit","MHC","chr16p11"),1.1,2.4),
             alpha=ifelse(bbj$class=="ns",0.35,0.9)) +
  geom_hline(yintercept=-log10(BONF), linetype="dashed", color="grey40") +
  annotate("text", x=max(bbj$beta,na.rm=TRUE)*0.7, y=-log10(BONF)+0.4,
           label="Bonferroni (1.97e-5; 2,545 genes)", size=3, color="grey40") +
  geom_vline(xintercept=0, linetype="dotted", color="grey60") +
  scale_color_manual(values=pal, name="Class",
                     breaks=c("known_anchor","IGF1R","MHC","chr16p11","other_hit"),
                     labels=c("Known anchor (TSHR/CTLA4)","IGF1R (effector)",
                              "MHC spillover","chr16p11.2 LD cluster","Other hit")) +
  geom_text_repel(aes(label=lab), size=3, max.overlaps=30,
                  segment.color="grey60", min.segment.length=0) +
  labs(x="MR effect on Graves disease (log-odds, BBJ)",
       y=expression(-log[10](italic(P))),
       title="Druggable-gene-wide BBJ discovery MR",
       subtitle="2,235 genes with BBJ estimates among 2,545 MR-testable druggable genes") +
  theme_bw(base_size=12) +
  theme(legend.position="right",
        plot.title=element_text(face="bold", size=13),
        panel.grid.minor=element_blank())

ggsave(file.path(OUTDIR,"Figure2_MR_volcano.png"), p, width=9, height=6.5, dpi=300)
ggsave(file.path(OUTDIR,"Figure2_MR_volcano.pdf"), p, width=9, height=6.5)

cat("Figure 2 written. N genes plotted:", nrow(bbj),
    "| Bonferroni hits:", nrow(bbj[pvalue<BONF]), "\n")
cat("Backbone positions:\n")
print(bbj[gene_symbol %in% c("TSHR","IGF1R","CTLA4"),
          .(gene_symbol, beta, pvalue, neglog10P, class)])
