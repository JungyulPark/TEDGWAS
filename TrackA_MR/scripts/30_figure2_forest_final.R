# ============================================================
# Figure 2 v5 — TSHR vs IGF1R cis-MR
# v5 changes:
#  - Title: "TSHR but not IGF1R..." → "TSHR shows robust causal 
#    association, whereas IGF1R signals are phenotype-dependent"
#  - Subtitle: "two independent cohorts, two ancestries" 
#    → "two independent outcome datasets and ancestries"
#  - 3-tier marker scheme (preserved from v4)
# ============================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# ============================================================
# Verified ground truth (BBJ + UKB re-run)
# ============================================================

panel_a_data <- data.frame(
  label = c("TSHR — Wald ratio (N=1)",
            "IGF1R — Inverse-variance weighted (N=10)",
            "IGF1R — MR-Egger (N=10)",
            "IGF1R — Weighted median (N=10)",
            "IGF1R — Weighted mode (N=10)"),
  gene = c("TSHR", "IGF1R", "IGF1R", "IGF1R", "IGF1R"),
  beta = c(-1.408, +0.198, +0.311, +0.256, +0.253),
  se   = c( 0.182,  0.116,  0.233,  0.158,  0.537),
  pval = c(9.1e-15, 0.086, 0.219, 0.105, 0.649),
  stringsAsFactors = FALSE
)

panel_b_data <- data.frame(
  label = c("TSHR — Wald ratio (N=1)",
            "IGF1R — Inverse-variance weighted (N=10)",
            "IGF1R — MR-Egger (N=10)",
            "IGF1R — Weighted median (N=10)",
            "IGF1R — Weighted mode (N=10)"),
  gene = c("TSHR", "IGF1R", "IGF1R", "IGF1R", "IGF1R"),
  beta = c(-1.629, +0.165, +0.202, +0.273, +0.285),
  se   = c( 0.150,  0.085,  0.173,  0.097,  0.396),
  pval = c(1.07e-27, 0.052, 0.278, 0.005, 0.490),
  stringsAsFactors = FALSE
)

# ============================================================
# 3-TIER SIGNIFICANCE
# ============================================================

classify_sig <- function(pval) {
  if (pval < 2.08e-3) {
    return("Bonferroni-significant")
  } else if (pval < 0.05) {
    return("Nominal P<0.05")
  } else {
    return("Non-significant")
  }
}

panel_a_data$sig_class <- sapply(panel_a_data$pval, classify_sig)
panel_b_data$sig_class <- sapply(panel_b_data$pval, classify_sig)

prep_data <- function(d) {
  d$or       <- exp(d$beta)
  d$or_low   <- exp(d$beta - 1.96 * d$se)
  d$or_high  <- exp(d$beta + 1.96 * d$se)
  d$label    <- factor(d$label, levels = rev(d$label))
  d$pval_lab <- ifelse(d$pval < 1e-3,
                       sprintf("%.2e", d$pval),
                       sprintf("P = %.3f", d$pval))
  d$pval_lab <- ifelse(d$sig_class == "Bonferroni-significant",
                       paste0(d$pval_lab, " ***"),
                ifelse(d$sig_class == "Nominal P<0.05",
                       paste0(d$pval_lab, " *"),
                       d$pval_lab))
  d
}

panel_a_data <- prep_data(panel_a_data)
panel_b_data <- prep_data(panel_b_data)

# ============================================================
# Plot helper — 3-tier shape/fill scheme
# ============================================================

make_forest <- function(d, title, subtitle, x_max = 4.5) {
  d$fill_color <- with(d, ifelse(
    sig_class == "Non-significant", "white",
    ifelse(gene == "TSHR", "#C00000", "#4472C4")
  ))

  ggplot(d, aes(x = or, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = or_low, xmax = or_high, color = gene),
                   height = 0.25, linewidth = 0.7) +
    geom_point(aes(color = gene), shape = 21, size = 4.5, stroke = 1.5,
               fill = d$fill_color) +
    geom_text(aes(x = x_max * 1.05, label = pval_lab),
              color = "black", hjust = 0, size = 3.6) +
    scale_x_log10(limits = c(0.10, x_max * 1.6),
                  breaks = c(0.10, 0.25, 0.50, 1.00, 2.00, 4.00),
                  labels = sprintf("%.2f", c(0.10, 0.25, 0.50, 1.00, 2.00, 4.00))) +
    scale_color_manual(values = c("TSHR" = "#C00000", "IGF1R" = "#4472C4"),
                       name = "Gene") +
    labs(title    = title,
         subtitle = subtitle,
         x = "Odds ratio per 1-SD increase in gene expression (log scale)",
         y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0),
          plot.subtitle = element_text(size = 10, color = "gray40"),
          axis.text.y   = element_text(size = 10, color = "black"),
          axis.title.x  = element_text(size = 10),
          legend.position    = "right",
          legend.text        = element_text(size = 9),
          legend.title       = element_text(size = 10, face = "bold"),
          legend.key.size    = unit(0.5, "cm"),
          plot.margin   = margin(15, 12, 12, 10))
}

p_a <- make_forest(
  panel_a_data,
  title    = "A. Primary cis-MR on Graves disease (Biobank Japan)",
  subtitle = "East Asian ancestry (n = 212,453; 2,809 cases)"
)

p_b <- make_forest(
  panel_b_data,
  title    = "B. Replication on hyperthyroidism (UK Biobank, rescaled to log-odds)",
  subtitle = "European ancestry (n = 484,598; 3,731 cases)"
)

p_a <- p_a + theme(legend.position = "right",
                   plot.margin = margin(15, 100, 12, 10))
p_b <- p_b + theme(legend.position = "none",
                   plot.margin = margin(15, 100, 25, 10))

# ============================================================
# Combine — ★ NEW v5 TITLE + SUBTITLE
# ============================================================

final_plot <- p_a / p_b +
  plot_annotation(
    title    = "Figure 2. TSHR shows robust causal association, whereas IGF1R signals are phenotype-dependent",
    subtitle = "Cis-MR across two independent outcome datasets and ancestries",
    caption  = paste0(
      "Marker fill: filled = significant; open = non-significant. ",
      "*** Bonferroni-significant (P < 2.08\u00d710\u207b\u00b3, 24 tests). ",
      "* Nominal P<0.05 (not Bonferroni-significant). Error bars: 95% CI."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0),
      plot.caption  = element_text(size = 9, color = "gray40", hjust = 0,
                                   margin = margin(t = 10), lineheight = 1.2)
    )
  ) +
  plot_layout(heights = c(1.0, 1.0))

ggsave(
  filename = "Figure2.png",
  plot     = final_plot,
  width    = 14, height = 10, units = "in",
  dpi      = 600,
  bg       = "white"
)

cat("Figure 2 v5 saved.\n")
cat("Title: 'TSHR shows robust causal association, whereas IGF1R signals are phenotype-dependent'\n")
cat("Subtitle: 'Cis-MR across two independent outcome datasets and ancestries'\n")
