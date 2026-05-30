# ============================================================
# Figure 3 v7.1 — Three-panel vertical stack (TED-TRAP v4.5)
# v7.1 변경 (FROM v7):
#  - a/b/c/d footnote(legend)를 subtitle에서 PLOT 내부 annotation으로 이동
#  - 표준 figure 관례: panel 내부에 inset legend 배치
#  - Panel C subtitle은 핵심 메시지만 유지
# ============================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# ============================================================
# PANEL A — TSHR per-sample TPM
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
  geom_label(data = data.frame(x = 4.5, y = ted_mean,
                                label = sprintf("TED mean = %.2f", ted_mean)),
             aes(x = x, y = y, label = label),
             color = "darkred", size = 3.5, fontface = "italic",
             fill = "white", label.size = 0,
             label.padding = unit(0.15, "lines"),
             hjust = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Control" = "#808080", "TED" = "#C00000"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0),
                     breaks = seq(0, 1.0, 0.2)) +
  labs(title    = "(A) TSHR per-sample expression in orbital tissue",
       subtitle = "Exploratory in-house signal: log2FC = +2.33, padj = 0.006",
       x = "Sample",
       y = "TSHR expression (TPM)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0,
                                     margin = margin(b = 2)),
        plot.subtitle = element_text(face = "italic", size = 10,
                                     color = "gray40", hjust = 0,
                                     margin = margin(b = 4)),
        # ★ v7: plot.title.position 제거 (default 'panel' 사용)
        # → title이 차트와 일관된 위치 (안쪽으로 들어감)
        axis.title    = element_text(size = 11),
        axis.text     = element_text(size = 10),
        plot.margin   = margin(8, 15, 8, 5))

# ============================================================
# PANEL B — Insulin cassette 3-cohort
# ============================================================

panel_b_data <- data.frame(
  Gene = factor(rep(c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1"), 3),
                levels = c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1")),
  Cohort = factor(rep(c("In-house orbital adipose (n=5)\u2020",
                        "GSE58331 orbital adipose (n=43)",
                        "GSE105149 lacrimal gland (n=11)"), each = 5),
                  levels = c("In-house orbital adipose (n=5)\u2020",
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

panel_b_subtitle <- "\u2020 In-house = 4 TED orbital adipose vs 1 nasal fat control. *P < 0.05 (limma)."

p_b <- ggplot(panel_b_data, aes(x = Gene, y = log2FC, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(aes(label = label,
                vjust = ifelse(log2FC >= 0, -0.4, 1.3)),
            position = position_dodge(width = 0.8),
            size = 3.0) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "In-house orbital adipose (n=5)\u2020"  = "#ED7D31",
    "GSE58331 orbital adipose (n=43)"        = "#4472C4",
    "GSE105149 lacrimal gland (n=11)"        = "#7030A0"
  )) +
  scale_y_continuous(limits = c(-1.0, 1.3), breaks = seq(-1.0, 1.2, 0.5)) +
  labs(title    = "(B) Insulin cassette: in-house RNA-seq vs two external orbital cohorts",
       subtitle = panel_b_subtitle,
       x        = "Gene",
       y        = expression(log[2]~"fold-change (TED vs control)"),
       fill     = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0,
                                     margin = margin(b = 2)),
        plot.subtitle = element_text(face = "italic", size = 10,
                                     color = "gray40", hjust = 0,
                                     margin = margin(b = 4)),
        axis.text.x   = element_text(face = "italic", size = 10),
        axis.title    = element_text(size = 11),
        legend.position    = "top",
        legend.text        = element_text(size = 9),
        legend.title       = element_blank(),
        legend.key.size    = unit(0.4, "cm"),
        legend.box.margin  = margin(t = 0, b = 3),
        plot.margin        = margin(8, 15, 8, 5))

# ============================================================
# PANEL C — IGF1R cross-phenotype MR
# ★ v7: y-axis labels을 a/b/c/d로 (1글자!) → axis padding 최소
# ============================================================

panel_c_data <- data.frame(
  outcome = factor(c("a", "b", "c", "d"),
                   levels = c("d", "c", "b", "a")),  # bottom-up
  outcome_full = c("Fasting glucose", "HbA1c", "T2D",
                   "Height (IGF-axis pathway anchor)"),
  beta    = c(0.005, 0.005, 0.010, -0.050),
  ci_low  = c(-0.020, -0.020, -0.045, -0.130),
  ci_high = c(0.030, 0.030, 0.065, 0.020),
  pvalue  = c(0.533, 0.509, 0.702, 0.140),
  n_iv    = c(8, 8, 7, 8),
  category = c("Metabolic off-target", "Metabolic off-target",
               "Metabolic off-target", "IGF-axis pathway anchor")
)

panel_c_data$label_p <- sprintf("P = %.3f (n_IV = %d)",
                                panel_c_data$pvalue,
                                panel_c_data$n_iv)

# Subtitle: 핵심 메시지만 (a/b/c/d 정의는 plot 내부 annotation)
panel_c_subtitle <- "Germline IGF1R does not cause metabolic phenotypes (null IVW); height retains pathway-level activity"


p_c <- ggplot(panel_c_data, aes(x = beta, y = outcome, color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.2, linewidth = 0.7) +
  geom_point(aes(shape = category), size = 4.5, fill = NA, stroke = 1.2) +
  geom_text(aes(x = 0.20, label = label_p), color = "black",
            size = 3.4, hjust = 1) +
  annotate("text", x = 0.20, y = "d",
           label = "WM P = 0.019 (anchor)", hjust = 1,
           vjust = 2.2, color = "#BF9000", size = 3.2,
           fontface = "italic") +
  # ★ v7.1: a/b/c/d 정의를 plot 내부 좌측 상단에 inset annotation으로
  annotate("label", x = -0.18, y = "a",
           label = "a = Fasting glucose\nb = HbA1c\nc = T2D\nd = Height (IGF-axis pathway anchor)",
           hjust = 0, vjust = 1.0,
           color = "gray30", size = 3.0,
           fill = "white", label.size = 0.3,
           label.padding = unit(0.3, "lines")) +
  scale_color_manual(values = c(
    "Metabolic off-target"     = "#4472C4",
    "IGF-axis pathway anchor"  = "#BF9000"
  )) +
  scale_shape_manual(values = c(
    "Metabolic off-target"     = 23,
    "IGF-axis pathway anchor"  = 23
  )) +
  scale_x_continuous(limits = c(-0.20, 0.22),
                     breaks = c(-0.1, 0, 0.1, 0.2)) +
  labs(title    = "(C) IGF1R cross-phenotype Mendelian randomization",
       subtitle = panel_c_subtitle,
       x        = expression(paste("IGF1R cis-eQTL causal effect on outcome (", beta, ", 95% CI)")),
       y        = NULL,
       color    = NULL, shape = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0,
                                     margin = margin(b = 2)),
        plot.subtitle = element_text(face = "italic", size = 10,
                                     color = "gray40", hjust = 0,
                                     margin = margin(b = 4)),
        # ★ v7.1: y-axis labels이 a/b/c/d (1글자) → bold + 큰 size
        axis.text.y   = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x  = element_text(size = 10),
        legend.position    = "top",
        legend.text        = element_text(size = 9),
        legend.key.size    = unit(0.4, "cm"),
        legend.box.margin  = margin(t = 0, b = 3),
        plot.margin        = margin(8, 15, 8, 5))

# ============================================================
# Vertical stack — patchwork align 영향 최소
# ============================================================

final_plot <- (p_a / p_b / p_c +
  plot_layout(heights = c(1.0, 1.3, 1.0))) +
  plot_annotation(
    title    = "Figure 3. Tissue-level and cross-phenotype evidence for TSHR vs IGF-1R in TED",
    subtitle = "Three-panel composite: (A) TSHR per-sample, (B) insulin cassette across cohorts, (C) IGF1R cross-phenotype MR",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0,
                                   margin = margin(t = 5, b = 3)),
      plot.subtitle = element_text(face = "italic", size = 11,
                                   color = "gray30", hjust = 0,
                                   margin = margin(b = 8),
                                   lineheight = 1.2)
    )
  )

ggsave(
  filename = "Figure3.png",
  plot     = final_plot,
  width    = 11, height = 14, units = "in",
  dpi      = 470,
  bg       = "white"
)

cat("Figure 3 v7.1 saved.\n")
cat("핵심 변경 (FROM v7):\n")
cat("  ✅ a/b/c/d 정의를 plot 내부 좌측 상단 annotation으로 이동\n")
cat("  ✅ 표준 figure 관례: panel 내부 inset legend\n")
cat("Dimensions: 11 x 14 in @ 470 DPI = 5170 x 6580 px = 34 MP\n")
