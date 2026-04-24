# =============================================================================
# TED-TRAP — Figure 2: MR Forest Plot (Publication-Quality)
# Panel A: BBJ Primary (East Asian) | Panel B: UKB Replication (European, rescaled)
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

# Install if needed
if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork", repos = "https://cloud.r-project.org")
}
library(patchwork)

if (!dir.exists("TrackA_MR/figures")) dir.create("TrackA_MR/figures", recursive = TRUE)

# =============================================================================
# Load rescaled MR results
# =============================================================================

mr <- read.csv("TrackA_MR/results/02b_mr_rescaled.csv", stringsAsFactors = FALSE)

# Focus on TSHR + IGF1R, and target methods
target_methods <- c(
    "Wald ratio", "Inverse variance weighted",
    "MR Egger", "Weighted median", "Weighted mode"
)

mr_focus <- mr[mr$gene %in% c("TSHR", "IGF1R") &
    mr$method %in% target_methods, ]

# =============================================================================
# Build BBJ Panel (Primary, log-odds native)
# =============================================================================

bbj <- mr_focus[mr_focus$outcome_role == "Primary", ]
bbj$OR <- exp(bbj$b)
bbj$OR_lo <- exp(bbj$b - 1.96 * bbj$se)
bbj$OR_hi <- exp(bbj$b + 1.96 * bbj$se)
bbj$label <- sprintf("%s — %s (N=%d)", bbj$gene, bbj$method, bbj$nsnp)

# Order: TSHR on top, then IGF1R methods
bbj$method <- factor(bbj$method, levels = target_methods)
bbj <- bbj[order(bbj$gene != "TSHR", bbj$method), ]
bbj$row_order <- seq_len(nrow(bbj))
bbj$label <- factor(bbj$label, levels = rev(bbj$label))

# Significance marker
bbj$sig_marker <- ifelse(bbj$pval < 2.08e-3, "***",
    ifelse(bbj$pval < 0.05, "*", "ns")
)

# P-value display
bbj$p_display <- ifelse(bbj$pval < 0.001,
    formatC(bbj$pval, format = "e", digits = 2),
    sprintf("P = %.3f", bbj$pval)
)

# =============================================================================
# Panel A — BBJ (Primary)
# =============================================================================

p_bbj <- ggplot(bbj, aes(x = OR, y = label, color = gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.25, linewidth = 0.7) +
    geom_point(aes(shape = sig_marker), size = 3.5, stroke = 1.2) +
    scale_shape_manual(
        values = c("***" = 19, "*" = 17, "ns" = 21),
        name = "Significance",
        labels = c(
            "*** P < 2.08×10⁻³ (Bonferroni)",
            "* P < 0.05 (nominal)",
            "Not significant"
        )
    ) +
    scale_color_manual(
        values = c("TSHR" = "#c0392b", "IGF1R" = "#2980b9"),
        name = "Gene"
    ) +
    # PATCH: x-axis expanded from [0.1, 4] → [0.1, 10]
    scale_x_log10(
        breaks = c(0.1, 0.25, 0.5, 1, 2, 4),
        limits = c(0.1, 10)
    ) +
    # PATCH: P-label position moved from x=3.2 → x=5.5 (outside widest CI)
    geom_text(aes(label = p_display, x = 5.5),
        hjust = 0, size = 3.2, color = "black",
        show.legend = FALSE
    ) +
    labs(
        title = "A. Primary cis-MR on Graves disease (Biobank Japan)",
        subtitle = "East Asian ancestry (n = 212,453; 2,809 cases)",
        x = "Odds ratio per 1-SD increase in gene expression (log scale)",
        y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "grey30", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        plot.margin = margin(10, 30, 10, 10) # PATCH: right margin reduced (70 → 30)
    ) +
    coord_cartesian(clip = "off")

# =============================================================================
# Panel B — UKB (Replication, rescaled)
# =============================================================================

ukb <- mr_focus[mr_focus$outcome_role == "Replication" &
    mr_focus$method %in% c("Wald ratio", "Inverse variance weighted"), ]
ukb$OR <- exp(ukb$b_rescaled)
ukb$OR_lo <- exp(ukb$b_rescaled - 1.96 * ukb$se_rescaled)
ukb$OR_hi <- exp(ukb$b_rescaled + 1.96 * ukb$se_rescaled)
ukb$label <- sprintf("%s — %s (N=%d)", ukb$gene, ukb$method, ukb$nsnp)
ukb <- ukb[order(ukb$gene != "TSHR"), ]
ukb$label <- factor(ukb$label, levels = rev(ukb$label))
ukb$sig_marker <- ifelse(ukb$pval < 2.08e-3, "***",
    ifelse(ukb$pval < 0.05, "*", "ns")
)
ukb$p_display <- ifelse(ukb$pval < 0.001,
    formatC(ukb$pval, format = "e", digits = 2),
    sprintf("P = %.3f", ukb$pval)
)

p_ukb <- ggplot(ukb, aes(x = OR, y = label, color = gene)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.25, linewidth = 0.7) +
    geom_point(aes(shape = sig_marker), size = 3.5, stroke = 1.2) +
    scale_shape_manual(values = c("***" = 19, "*" = 17, "ns" = 21), guide = "none") +
    scale_color_manual(values = c("TSHR" = "#c0392b", "IGF1R" = "#2980b9"), guide = "none") +
    # PATCH: same x-axis expansion
    scale_x_log10(
        breaks = c(0.1, 0.25, 0.5, 1, 2, 4),
        limits = c(0.1, 10)
    ) +
    # PATCH: same P-label repositioning
    geom_text(aes(label = p_display, x = 5.5),
        hjust = 0, size = 3.2, color = "black"
    ) +
    labs(
        title = "B. Replication on hyperthyroidism (UK Biobank, rescaled to log-odds)",
        subtitle = "European ancestry (n = 484,598; 3,731 cases)",
        x = "Odds ratio per 1-SD increase in gene expression (log scale)",
        y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "grey30", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(10, 30, 10, 10) # PATCH: consistent
    ) +
    coord_cartesian(clip = "off")

# =============================================================================
# Combine with patchwork
# =============================================================================

fig2 <- (p_bbj / p_ukb) +
    plot_layout(heights = c(6, 2)) +
    plot_annotation(
        title = "Figure 2. TSHR but not IGF1R is causally associated with Graves disease",
        subtitle = "Cis-MR across five methods, two independent cohorts, two ancestries",
        caption = "*** Bonferroni-significant (P < 2.08×10⁻³, 24 tests). Error bars: 95% CI.",
        theme = theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(color = "grey25", size = 11)
        )
    )

# Save at publication resolution
ggsave("TrackA_MR/figures/Figure2_MR_forest.pdf", fig2,
    width = 11, height = 7, device = grDevices::cairo_pdf
)
ggsave("TrackA_MR/figures/Figure2_MR_forest.png", fig2,
    width = 11, height = 7, dpi = 600
)
ggsave("TrackA_MR/figures/Figure2_MR_forest.tiff", fig2,
    width = 11, height = 7, dpi = 600, compression = "lzw"
)

cat("✅ Figure 2 PATCHED and re-saved\n")
