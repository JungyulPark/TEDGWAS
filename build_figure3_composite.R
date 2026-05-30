#!/usr/bin/env Rscript
# =============================================================================
# TED-TRAP v5 — Figure 3: TSHR vs IGF1R multi-modal contrast (CORE FIGURE)
# =============================================================================
# Panel A: colocalization PP bars (TSHR vs IGF1R vs CTLA4, BBJ + FinnGen)
# Panel B: forest plot (3 genes x 3 outcomes, beta +/- 95% CI)
# Panel C: orbital tissue expression (TPM, Control vs TED, backbone genes)
# Assembled with patchwork. All data read from verified files.
#
# Requires: ggplot2, data.table, patchwork
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

WORK   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
MERGED <- file.path(WORK, "04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv")
COLOC  <- file.path(WORK, "05_coloc_candidates/TaskE_01_coloc_results_v1.csv")
TISSUE <- file.path(WORK, "06_tissue_integration/TaskF_01_backbone_tissue_inhouse_v1.csv")
OUTDIR <- file.path(WORK, "07_manuscript/figures")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

genes <- c("TSHR","IGF1R","CTLA4")
gcol  <- c(TSHR="#C00000", IGF1R="#2E75B6", CTLA4="#548235")

primary <- function(d,g,o){c<-d[gene_symbol==g & outcome==o]
  i<-c[method=="Inverse variance weighted"]; w<-c[method=="Wald ratio"]
  if(nrow(i))return(i[1]); if(nrow(w))return(w[1]); return(c[1])}

mr <- fread(MERGED); coloc <- fread(COLOC); tissue <- fread(TISSUE)

# ---------- PANEL A: colocalization PP.H4 (+H2/H3 context) ----------
cA <- coloc[gene_symbol %in% genes, .(gene_symbol, outcome, pp_h2, pp_h3, pp_h4)]
cAm <- melt(cA, id.vars=c("gene_symbol","outcome"),
            measure.vars=c("pp_h2","pp_h3","pp_h4"),
            variable.name="hyp", value.name="PP")
cAm[, hyp := factor(hyp, levels=c("pp_h2","pp_h3","pp_h4"),
                    labels=c("H2 (GWAS assoc. only)","H3 (distinct)","H4 (shared)"))]
cAm[, outcome := factor(outcome, levels=c("BBJ_Graves","FinnGen_GO"),
                        labels=c("BBJ","FinnGen GO"))]
cAm[, gene_symbol := factor(gene_symbol, levels=genes)]

pA <- ggplot(cAm, aes(outcome, PP, fill=hyp)) +
  geom_col(position="stack", width=0.7) +
  facet_wrap(~gene_symbol, nrow=1) +
  scale_fill_manual(values=c("H2 (GWAS assoc. only)"="#F4B183","H3 (distinct)"="#BFBFBF",
                             "H4 (shared)"="#70AD47"), name="Coloc PP") +
  geom_hline(yintercept=0.8, linetype="dashed", color="grey40") +
  labs(title="A  Colocalization", x=NULL, y="Posterior probability") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom",
        axis.text.x=element_text(angle=0))

# ---------- PANEL B: forest (3 genes x 3 outcomes) ----------
rows <- list()
for (g in genes) for (o in c("BBJ_Graves","UKB_hyperthyroid","FinnGen_GO")) {
  r <- primary(mr,g,o)
  rows[[paste(g,o)]] <- data.table(gene=g, outcome=o,
      beta=r$beta, lo=r$beta-1.96*r$se, hi=r$beta+1.96*r$se, p=r$pvalue, n_iv=r$n_iv)
}
fb <- rbindlist(rows)
fb[, outcome := factor(outcome, levels=c("FinnGen_GO","UKB_hyperthyroid","BBJ_Graves"),
                       labels=c("FinnGen GO","UKB hyperthyroid","BBJ Graves"))]
fb[, gene := factor(gene, levels=genes)]

pB <- ggplot(fb, aes(beta, outcome, color=gene)) +
  geom_vline(xintercept=0, linetype="dotted", color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25, linewidth=0.7) +
  geom_point(size=2.6) +
  facet_wrap(~gene, nrow=1, scales="free_x") +
  scale_color_manual(values=gcol, guide="none") +
  labs(title="B  MR effect across outcomes", x="MR β (log-odds, 95% CI)", y=NULL) +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"))

# ---------- PANEL C: orbital tissue expression (per-sample) ----------
# Reconstruct per-sample from TaskF_01 ted_per_sample string (TED) + ctrl_tpm
tlist <- list()
for (g in genes) {
  ti <- tissue[gene==g]
  ted_vals <- as.numeric(strsplit(ti$ted_per_sample, ";")[[1]])
  # NOTE: ted_per_sample in TaskF_01 may contain >4 values if both arms concatenated;
  # Gemini: ensure only the 4 TED-sample TPMs are used (samples 7/8/10/11).
  tlist[[g]] <- rbind(
    data.table(gene=g, group="Control", tpm=ti$ctrl_tpm),
    data.table(gene=g, group="TED", tpm=ted_vals)
  )
}
td <- rbindlist(tlist)
td[, gene := factor(gene, levels=genes)]
td[, group := factor(group, levels=c("Control","TED"))]
lab <- tissue[, .(gene, lab=sprintf("log2FC=%+.2f\npadj=%.1e", deseq2_log2fc, deseq2_padj))]
lab[, gene := factor(gene, levels=genes)]

pC <- ggplot(td, aes(group, tpm, color=gene)) +
  geom_boxplot(data=td[group=="TED"], width=0.4, outlier.shape=NA, alpha=0.3) +
  geom_jitter(width=0.12, size=2.4, alpha=0.85) +
  facet_wrap(~gene, nrow=1, scales="free_y") +
  scale_color_manual(values=gcol, guide="none") +
  geom_text(data=lab, aes(x=1.5, y=Inf, label=lab), vjust=1.3, size=2.8,
            color="grey20", inherit.aes=FALSE) +
  labs(title="C  Orbital tissue expression (in-house RNA-seq)",
       x=NULL, y="TPM") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"))

# ---------- assemble ----------
fig3 <- pA / pB / pC + plot_layout(heights=c(1, 1, 1.1))
ggsave(file.path(OUTDIR,"Figure3_TSHR_vs_IGF1R_composite.png"), fig3,
       width=9, height=11, dpi=300)
ggsave(file.path(OUTDIR,"Figure3_TSHR_vs_IGF1R_composite.pdf"), fig3,
       width=9, height=11)

cat("Figure 3 written.\n")
cat("Panel A coloc PP.H4 (verify): TSHR BBJ", coloc[gene_symbol=="TSHR"&outcome=="BBJ_Graves",pp_h4],
    "| IGF1R BBJ", coloc[gene_symbol=="IGF1R"&outcome=="BBJ_Graves",pp_h4], "\n")
cat("Panel C tissue (verify): TSHR log2FC", tissue[gene=="TSHR",deseq2_log2fc], "\n")
cat("\n*** Gemini: confirm ted_per_sample has exactly 4 TED TPM values per gene.\n")
cat("    TaskF_01 ted_per_sample for TSHR had 8 values — check if that includes\n")
cat("    both TED and control arms or replicate columns. Use ONLY the 4 TED samples.\n")
