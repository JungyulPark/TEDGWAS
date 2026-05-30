#!/usr/bin/env Rscript
# =============================================================================
# TED-TRAP v5 — Supplementary Figure S1: TSHR fine-mapping & LD rationale
# =============================================================================
# Purpose: show WHY TSHR is retained as a single-IV / locus-level anchor rather
# than a multi-IV exposure. Data from Table S3 (SuSiE credible sets).
#
# Panel A: credible-set lead SNP genomic positions (rs179252 = primary anchor)
# Panel B: ancestry-matched LD r2 of each secondary lead with rs179252
#
# Requires: ggplot2, data.table, patchwork
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

WORK   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
OUTDIR <- file.path(WORK, "07_manuscript/figures")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

# ---- data (from Table S3, verified) ----
cs <- data.table(
  panel   = c("European (1000G)","European (1000G)","East Asian (1000G)","East Asian (1000G)"),
  cs_id   = c("CS_1","CS_2","CS_1","CS_2"),
  snp     = c("rs179252","rs179248","rs3783950","rs3783949"),
  pos     = c(81435985, 81433038, 81448282, 81448382),
  r2_eur  = c(1.000, 0.965, 0.229, 0.207),
  r2_eas  = c(1.000, 0.992, 0.808, 0.808),
  is_primary = c(TRUE, FALSE, FALSE, FALSE)
)
# ancestry-matched r2 = r2 from the panel that detected the credible set
cs[, r2_matched := ifelse(panel=="European (1000G)", r2_eur, r2_eas)]
cs[, pos_kb := (pos - 81400000)/1000]   # position relative to 81.4 Mb, in kb
cols <- c("European (1000G)"="#2E75B6", "East Asian (1000G)"="#ED7D31")

# Numerical y-axis mapping to handle discrete levels with custom manual label nudges
cs[, y_num := ifelse(panel == "East Asian (1000G)", 1, 2)]
# Nudge rs179252 (up), rs179248 (down), rs3783950 (up), rs3783949 (down)
cs[, nudge_y := c(0.12, -0.12, 0.12, -0.12)]

# ---- Panel A: genomic positions ----
pA <- ggplot(cs, aes(pos_kb, y_num, color=panel)) +
  geom_segment(aes(x=min(pos_kb)-2, xend=max(pos_kb)+2, y=y_num, yend=y_num),
               color="grey85", linewidth=0.6) +
  geom_point(aes(size=is_primary, shape=is_primary)) +
  geom_text(aes(label=snp, y=y_num + nudge_y), size=3.2, show.legend=FALSE, fontface="bold") +
  scale_color_manual(values=cols, name="LD reference panel") +
  scale_size_manual(values=c(`TRUE`=4.5, `FALSE`=2.6), guide="none") +
  scale_shape_manual(values=c(`TRUE`=8, `FALSE`=16), guide="none") +  # star for primary
  scale_y_continuous(breaks=c(1, 2), labels=c("East Asian (1000G)", "European (1000G)"), limits=c(0.6, 2.4)) +
  labs(title="A  Credible-set lead SNP positions (TSHR locus, chr14)",
       x="Position relative to chr14:81.4 Mb (kb, hg19)", y=NULL) +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom",
        panel.grid.minor=element_blank())

# ---- Panel B: ancestry-matched LD with rs179252 ----
cs_sec <- cs[is_primary==FALSE]
cs_sec[, lab := sprintf("%s\n(%s)", snp, sub(" .*","",panel))]
pB <- ggplot(cs_sec, aes(reorder(snp, -r2_matched), r2_matched, fill=panel)) +
  geom_col(width=0.6) +
  geom_hline(yintercept=0.8, linetype="dashed", color="grey40") +
  annotate("text", x=0.7, y=0.83, label="r² = 0.8", size=3, color="grey40", hjust=0) +
  geom_text(aes(label=sprintf("%.3f", r2_matched)), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=cols, name="LD reference panel") +
  scale_y_continuous(limits=c(0,1.08), breaks=seq(0,1,0.25)) +
  labs(title="B  Ancestry-matched LD (r²) of secondary leads with rs179252",
       x="Secondary credible-set lead SNP", y="r² with rs179252\n(ancestry-matched panel)") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom",
        panel.grid.minor=element_blank())

sf1 <- pA / pB + plot_layout(heights=c(1, 1.15))
ggsave(file.path(OUTDIR,"FigureS1_TSHR_finemapping_LD.png"), sf1, width=8, height=8, dpi=300)
ggsave(file.path(OUTDIR,"FigureS1_TSHR_finemapping_LD.pdf"), sf1, width=8, height=8)

cat("Figure S1 written.\n")
cat("Secondary leads ancestry-matched r2 with rs179252 (all >= 0.808):\n")
print(cs_sec[, .(snp, panel, r2_matched)])
