#!/usr/bin/env Rscript
# =============================================================================
# TED-TRAP v5 — Supplementary Figure S2: Colocalization Profiles for Candidates
# =============================================================================
# Purpose: show why there are ZERO robust novel hits by demonstrating that
# all candidate genes fail to colocalize in FinnGen GO (replication/sensitivity).
# Data from Table S2 (candidate colocalization results).
#
# Requires: ggplot2, data.table, patchwork
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

WORK   <- "c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
S2_SRC <- file.path(WORK, "05_coloc_candidates/TaskE_02_candidate_coloc_v1.csv")
OUTDIR <- file.path(WORK, "07_manuscript/figures")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

# Load candidate coloc results
s2_raw <- fread(S2_SRC)

# Keep only columns of interest and melt for stacked bar chart
cA <- s2_raw[, .(gene_symbol, outcome, pp_h2, pp_h3, pp_h4)]
cAm <- melt(cA, id.vars=c("gene_symbol","outcome"),
            measure.vars=c("pp_h2","pp_h3","pp_h4"),
            variable.name="hyp", value.name="PP")

cAm[, hyp := factor(hyp, levels=c("pp_h2","pp_h3","pp_h4"),
                    labels=c("H2 (GWAS assoc. only)","H3 (distinct)","H4 (shared)"))]
cAm[, outcome := factor(outcome, levels=c("BBJ_Graves","FinnGen_GO"),
                        labels=c("BBJ","FinnGen GO"))]

# Define gene order for visual organization
cAm[, gene_symbol := factor(gene_symbol, levels=c("TNFSF14","IFNGR1","MAPKAPK5","HSD3B7","PRSS36","VKORC1"))]

# Match colors to Figure 3A
cols <- c("H2 (GWAS assoc. only)"="#F4B183", "H3 (distinct)"="#BFBFBF", "H4 (shared)"="#70AD47")

# Build the plot
sf2 <- ggplot(cAm, aes(outcome, PP, fill=hyp)) +
  geom_col(position="stack", width=0.65) +
  facet_wrap(~gene_symbol, nrow=2) +
  scale_fill_manual(values=cols, name="Coloc PP") +
  geom_hline(yintercept=0.8, linetype="dashed", color="grey40") +
  scale_y_continuous(limits=c(0, 1.05), breaks=seq(0, 1, 0.25)) +
  labs(title="Supplementary Figure S2. Colocalization probability profiles for candidate genes",
       x=NULL, y="Posterior probability") +
  theme_bw(base_size=11) +
  theme(plot.title=element_text(face="bold", size=12),
        legend.position="bottom",
        panel.grid.minor=element_blank(),
        strip.text=element_text(face="bold", size=10),
        axis.text.x=element_text(size=9.5))

# Save outputs
ggsave(file.path(OUTDIR, "FigureS2_candidate_coloc_failure.png"), sf2, width=9, height=7, dpi=300)
ggsave(file.path(OUTDIR, "FigureS2_candidate_coloc_failure.pdf"), sf2, width=9, height=7)

cat("Supplementary Figure S2 written successfully to:", OUTDIR, "\n")
