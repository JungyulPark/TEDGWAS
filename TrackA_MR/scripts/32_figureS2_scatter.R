# ============================================================
# Figure S2 — Three-cohort cross-validation scatter
# 라벨 겹침 해결 + ggrepel 사용
# Sage Thyroid publication-quality
# ============================================================

# Required packages: ggplot2, ggrepel, patchwork, scales
# install.packages(c("ggplot2", "ggrepel", "patchwork", "scales"))

library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)

# ============================================================
# DATA — locked v4.3 ground truth values
# ============================================================

# Panel A: in-house vs GSE58331 (matching tissue compartment)
panel_a_data <- data.frame(
  Gene  = c("TSHR",  "FOXO1", "PDPK1", "INSR",  "IRS2",  "PIK3R1"),
  in_house = c(2.33,  1.02,    0.38,    0.54,    0.65,    0.44),
  external = c(0.249, 0.21,    0.04,   -0.31,   -0.33,   -0.68),  # GSE58331 log2FC
  category = c("primary", "insulin", "insulin", "insulin", "insulin", "insulin"),
  concordant = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
)

# Panel B: in-house vs GSE105149 (cross-compartment, lacrimal gland)
panel_b_data <- data.frame(
  Gene  = c("TSHR", "IGF1R", "INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1"),
  in_house = c(2.33, 0.55,    0.54,   0.65,   1.02,    0.44,     0.38),
  external = c(0.24, 0.64,    0.30,   0.40,  -0.20,    0.30,    -0.15),  # GSE105149 log2FC
  category = c("primary", "primary", "insulin", "insulin", "insulin", "insulin", "insulin"),
  concordant = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)
)

# ============================================================
# Helper function — make a clean scatter panel
# ============================================================

make_panel <- function(data, ylab, title, xrange, yrange) {
  ggplot(data, aes(x = in_house, y = external)) +
    # Concordance background shading
    annotate("rect", xmin = 0, xmax = xrange[2], ymin = 0, ymax = yrange[2], 
             fill = "#90EE90", alpha = 0.15) +  # upper right concordant
    annotate("rect", xmin = xrange[1], xmax = 0, ymin = yrange[1], ymax = 0, 
             fill = "#90EE90", alpha = 0.15) +  # lower left concordant
    annotate("rect", xmin = 0, xmax = xrange[2], ymin = yrange[1], ymax = 0, 
             fill = "#FFB6B6", alpha = 0.15) +  # lower right discordant
    annotate("rect", xmin = xrange[1], xmax = 0, ymin = 0, ymax = yrange[2], 
             fill = "#FFB6B6", alpha = 0.15) +  # upper left discordant
    
    # Reference lines
    geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "gray50", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "gray70", 
                linetype = "dashed", linewidth = 0.4) +
    
    # Points — TSHR/IGF1R as diamond, others as circle
    geom_point(aes(shape = category, fill = concordant), 
               size = 4.5, color = "black", stroke = 0.5) +
    scale_shape_manual(values = c("primary" = 23, "insulin" = 21),
                       guide = "none") +
    scale_fill_manual(values = c("TRUE" = "#90EE90", "FALSE" = "#FFB6B6"),
                      guide = "none") +
    
    # Labels using ggrepel — automatic non-overlapping placement
    geom_text_repel(aes(label = Gene, fontface = ifelse(category == "primary", "bold", "plain")),
                    size = 4.2,
                    box.padding = 0.5,
                    point.padding = 0.4,
                    segment.color = "gray60",
                    segment.size = 0.3,
                    min.segment.length = 0.2,
                    max.overlaps = Inf,
                    seed = 42) +
    
    # Quadrant labels
    annotate("text", x = xrange[2]*0.85, y = yrange[2]*0.92, 
             label = "Concordant\n(UP)", size = 3.3, 
             color = "darkgreen", fontface = "bold", lineheight = 0.9) +
    annotate("text", x = xrange[2]*0.85, y = yrange[1]*0.92, 
             label = "Discordant", size = 3.3, 
             color = "darkred", fontface = "bold") +
    
    coord_cartesian(xlim = xrange, ylim = yrange) +
    labs(x = expression(paste("In-house orbital adipose ", log[2], "FC (n=5)")),
         y = ylab,
         title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
          axis.title = element_text(size = 11),
          panel.grid.major = element_line(color = "gray95", linewidth = 0.2),
          plot.margin = margin(10, 10, 10, 10))
}

# ============================================================
# Build panels
# ============================================================

p_a <- make_panel(
  panel_a_data,
  ylab = expression(paste("GSE58331 orbital adipose ", log[2], "FC (23 TED vs 20 normal)")),
  title = "(A) Matching tissue compartment",
  xrange = c(-0.3, 2.7),
  yrange = c(-0.85, 0.55)
)

p_b <- make_panel(
  panel_b_data,
  ylab = expression(paste("GSE105149 lacrimal gland ", log[2], "FC (4 TED vs 7 normal)")),
  title = "(B) Cross-compartment direction test",
  xrange = c(-0.3, 2.7),
  yrange = c(-0.4, 0.85)
)

# ============================================================
# Combine + save
# ============================================================

combined <- (p_a | p_b) +
  plot_annotation(
    title = "Supplementary Figure S2. Three-cohort cross-validation: in-house orbital adipose RNA-seq vs two external microarray cohorts",
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0))
  ) +
  plot_layout(widths = c(1, 1))

# Save at 600 DPI for Sage Thyroid
ggsave(
  "FigureS2.png",
  combined,
  width = 14, height = 7, units = "in",
  dpi = 600,
  bg = "white"
)

cat("Figure S2 saved at 600 DPI: FigureS2.png\n")
cat("Dimensions: 14 × 7 inches at 600 DPI = 8400 × 4200 px\n")
