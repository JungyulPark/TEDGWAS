# ============================================================
# Figure 3 v3 — Three-panel composite (TED-TRAP v4.4)
# v3 changes:
#  - Panel A annotation: "6.5× UP in TED" → "Exploratory in-house signal"
#  - Panel C category: "Instrument-validity anchor" → "IGF-axis pathway anchor"
#  - Panel C label: "Height (positive control anchor)" 
#                    → "Height (IGF-axis pathway anchor)"
# ============================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# ============================================================
# PANEL A — TSHR per-sample TPM (★ Exploratory annotation)
# ============================================================

panel_a_data <- data.frame(
  sample = factor(c("Ctrl", "T1", "T2", "T3", "T4"),
                  levels = c("Ctrl", "T1", "T2", "T3", "T4")),
  TPM    = c(0.10, 0.72, 0.43, 0.49, 0.96),
  group  = c("Control", "TED", "TED", "TED", "TED")
)

ted_mean <- mean(panel_a_data$TPM[panel_a_data$group == "TED"])

p_a <- ggplot(panel_a_data, aes(x = sample, y = TPM, fill = group)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", TPM)),
            vjust = -0.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = ted_mean, linetype = "dashed",
             color = "darkred", linewidth = 0.6) +
  annotate("text", x = 4.5, y = ted_mean + 0.08,
           label = sprintf("TED mean = %.2f", ted_mean),
           color = "darkred", size = 3.5, fontface = "italic") +
  # ★ NEW: "Exploratory in-house signal" 표현
  annotate("text", x = 3, y = 1.20,
           label = "Exploratory in-house signal: log2FC = +2.33, padj = 0.006",
           size = 4.0, fontface = "bold") +
  scale_fill_manual(values = c("Control" = "#808080", "TED" = "#C00000"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 1.35), expand = c(0, 0),
                     breaks = seq(0, 1.2, 0.2)) +
  labs(title = "(A) TSHR per-sample expression in orbital tissue",
       subtitle = "Dashed red line = mean TED expression (T1-T4); single non-TED control",
       x = "Sample",
       y = "TSHR expression (TPM)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(face = "italic", size = 10,
                                     color = "gray40"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        plot.margin = margin(10, 12, 10, 10))

# ============================================================
# PANEL B — Insulin cassette 3-cohort (UNCHANGED)
# ============================================================

panel_b_data <- data.frame(
  Gene = factor(rep(c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1"), 3),
                levels = c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1")),
  Cohort = factor(rep(c("In-house RNA-seq (n=5)",
                        "GSE58331 orbital adipose (n=43)",
                        "GSE105149 lacrimal gland (n=11)"), each = 5),
                  levels = c("In-house RNA-seq (n=5)",
                             "GSE58331 orbital adipose (n=43)",
                             "GSE105149 lacrimal gland (n=11)")),
  log2FC = c(
    +0.54, +0.65, +1.02, +0.44, +0.38,
    -0.31, -0.33, +0.21, -0.68, +0.04,
    +0.30, +0.40, -0.20, +0.30, -0.15
  ),
  is_significant = c(
    rep(FALSE, 5),
    TRUE, FALSE, FALSE, TRUE, FALSE,
    rep(FALSE, 5)
  )
)

panel_b_data$label <- sprintf("%+.2f%s",
                              panel_b_data$log2FC,
                              ifelse(panel_b_data$is_significant, "*", ""))

panel_b_subtitle <- paste0(
  "In-house 5/5 upregulation NOT reproduced in matching-compartment GSE58331;\n",
  "receptor-proximal effectors (INSR, IRS2, PIK3R1) directionally opposite.\n",
  "Pattern: exploratory single-cohort observation. *P < 0.05 (limma)."
)

p_b <- ggplot(panel_b_data, aes(x = Gene, y = log2FC, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(aes(label = label,
                vjust = ifelse(log2FC >= 0, -0.4, 1.3)),
            position = position_dodge(width = 0.8),
            size = 2.9) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "In-house RNA-seq (n=5)"           = "#ED7D31",
    "GSE58331 orbital adipose (n=43)"  = "#4472C4",
    "GSE105149 lacrimal gland (n=11)"  = "#7030A0"
  )) +
  scale_y_continuous(limits = c(-1.0, 1.3), breaks = seq(-1.0, 1.2, 0.5)) +
  labs(title    = "(B) Insulin cassette: in-house RNA-seq vs two external orbital cohorts",
       subtitle = panel_b_subtitle,
       x        = "Gene",
       y        = expression(log[2]~"fold-change (TED vs control)"),
       fill     = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(face = "italic", size = 9,
                                     color = "gray40", lineheight = 1.2,
                                     margin = margin(b = 8)),
        axis.text.x   = element_text(face = "italic", size = 10),
        axis.title    = element_text(size = 11),
        legend.position    = "top",
        legend.text        = element_text(size = 9),
        legend.key.size    = unit(0.4, "cm"),
        legend.box.margin  = margin(t = 0, b = 5),
        plot.margin   = margin(10, 12, 10, 10))

# ============================================================
# PANEL C — IGF1R cross-phenotype MR
# ★ Category renamed: "Instrument-validity anchor" → "IGF-axis pathway anchor"
# ★ Label renamed: "Height (positive control anchor)" → "Height (IGF-axis pathway anchor)"
# ============================================================

panel_c_data <- data.frame(
  outcome = factor(c("Fasting glucose",
                     "HbA1c",
                     "Type 2 diabetes (T2D)",
                     "Height (IGF-axis pathway anchor)"),  # ★ NEW LABEL
                   levels = c("Height (IGF-axis pathway anchor)",
                              "Type 2 diabetes (T2D)",
                              "HbA1c",
                              "Fasting glucose")),
  beta    = c(0.005, 0.005, 0.010, -0.050),
  ci_low  = c(-0.020, -0.020, -0.045, -0.130),
  ci_high = c(0.030, 0.030, 0.065, 0.020),
  pvalue  = c(0.533, 0.509, 0.702, 0.140),
  n_iv    = c(8, 8, 7, 8),
  category = c("Metabolic off-target", "Metabolic off-target",
               "Metabolic off-target", "IGF-axis pathway anchor")  # ★ NEW
)

panel_c_data$label_p <- sprintf("P = %.3f (n_IV = %d)",
                                panel_c_data$pvalue,
                                panel_c_data$n_iv)

p_c <- ggplot(panel_c_data, aes(x = beta, y = outcome, color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.2, linewidth = 0.7) +
  geom_point(aes(shape = category), size = 4.5, fill = NA, stroke = 1.2) +
  geom_text(aes(x = 0.20, label = label_p), color = "black",
            size = 3.6, hjust = 1) +
  annotate("text", x = 0.20, y = "Height (IGF-axis pathway anchor)",
           label = "WM P = 0.019 (anchor)", hjust = 1,
           vjust = 2.5, color = "#BF9000", size = 3.4,
           fontface = "italic") +
  scale_color_manual(values = c(
    "Metabolic off-target"     = "#4472C4",
    "IGF-axis pathway anchor"  = "#BF9000"   # ★ NEW
  )) +
  scale_shape_manual(values = c(
    "Metabolic off-target"     = 23,
    "IGF-axis pathway anchor"  = 23           # ★ NEW
  )) +
  scale_x_continuous(limits = c(-0.18, 0.22),
                     breaks = c(-0.1, 0, 0.1, 0.2)) +
  labs(title    = "(C) IGF1R cross-phenotype Mendelian randomization",
       subtitle = paste0("Germline IGF1R perturbation does not cause metabolic phenotypes (null IVW);\n",
                         "height instruments retain partial pathway-level activity"),
       x        = expression(paste("IGF1R cis-eQTL causal effect on outcome (", beta, ", 95% CI)")),
       y        = NULL,
       color    = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(face = "italic", size = 9,
                                     color = "gray40", lineheight = 1.2,
                                     margin = margin(b = 8)),
        axis.text.y   = element_text(size = 10, color = "black"),
        axis.title.x  = element_text(size = 10),
        legend.position    = "top",
        legend.text        = element_text(size = 9),
        legend.key.size    = unit(0.4, "cm"),
        legend.box.margin  = margin(t = 0, b = 5),
        plot.margin   = margin(15, 12, 10, 10))

# ============================================================
# Assemble
# ============================================================

top_row    <- p_a + p_b + plot_layout(widths = c(1, 1.5))
final_plot <- (top_row / p_c) + 
  plot_annotation(
    title = "Figure 3. TSHR upregulation, insulin cassette comparison, and IGF1R cross-phenotype MR",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0))
  ) +
  plot_layout(heights = c(1.1, 0.85))

ggsave(
  filename = "Figure3.png",
  plot     = final_plot,
  width    = 15, height = 10, units = "in",
  dpi      = 600,
  bg       = "white"
)

cat("Figure 3 v3 saved.\n")
cat("Panel A: 'Exploratory in-house signal'\n")
cat("Panel C: 'IGF-axis pathway anchor' (was 'positive control anchor')\n")
