#!/usr/bin/env Rscript
# =============================================================================
# build_figure1_studydesign.R
# =============================================================================
# Drawing the Figure 1 study design schematic using ggplot2 and grid.
# Grayscale, high-contrast publication design based on Figure1_editable.pptx
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(grid)
})

# Create directory if not exists
outdir <- "C:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/07_manuscript/figures"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Base Plot Setup
p <- ggplot() +
    theme_void() +
    xlim(-1.0, 11.0) +
    ylim(-0.1, 11.6)

# =============================================================================
# BAND 1: DATA SOURCES
# =============================================================================

# Helper to draw standard Band 1 Header Box
draw_band1_box <- function(plot_obj, xmin, xmax, ymin, ymax, header_text, body_lines) {
    # Main Box
    plot_obj <- plot_obj +
        annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                 fill = "#FFFFFF", color = "black", linewidth = 0.6) +
        # Header Box
        annotate("rect", xmin = xmin, xmax = xmax, ymin = ymax - 0.35, ymax = ymax,
                 fill = "#1A1A1A", color = "black", linewidth = 0.6) +
        # Header Text
        annotate("text", x = (xmin + xmax)/2, y = ymax - 0.175,
                 label = header_text, color = "#FFFFFF", fontface = "bold", size = 3.3, hjust = 0.5)
    
    # Body Text (centered)
    y_start <- ymax - 0.62
    for (i in seq_along(body_lines)) {
        plot_obj <- plot_obj +
            annotate("text", x = (xmin + xmax)/2, y = y_start - (i-1)*0.24,
                     label = body_lines[i], color = "#000000",
                     fontface = if (i == 1) "bold" else "plain",
                     size = 3.0, hjust = 0.5)
    }
    return(plot_obj)
}

p <- draw_band1_box(p, 0.4, 3.4, 8.1, 9.7, "Genetic exposure",
                    c("eQTLGen", "blood cis-eQTL", "predominantly European", "ancestry"))

p <- draw_band1_box(p, 3.7, 6.7, 8.1, 9.7, "Disease outcomes (GWAS)",
                    c("BBJ Graves disease", "discovery", "UKB hyperthyroidism", "replication", "FinnGen Graves", "ophthalmopathy", "TED-specific sensitivity"))

p <- draw_band1_box(p, 7.0, 10.0, 8.1, 9.7, "Orbital tissue",
                    c("In-house RNA-seq", "TED vs control", "orbital tissue", "Korean cohort"))


# =============================================================================
# BAND 2: ANALYSIS
# =============================================================================

# Helper to draw Band 2 step box
draw_band2_box <- function(plot_obj, xmin, xmax, ymin, ymax, step_num, title, subtitle) {
    plot_obj <- plot_obj +
        # Main gray card
        annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                 fill = "#F7F7F7", color = "grey40", linewidth = 0.6) +
        # Number Circle background (as a point)
        annotate("point", x = xmin + 0.35, y = (ymin + ymax)/2, size = 6.5, color = "#1A1A1A") +
        # Number Circle text
        annotate("text", x = xmin + 0.35, y = (ymin + ymax)/2, label = step_num,
                 color = "#FFFFFF", fontface = "bold", size = 2.8, hjust = 0.5, vjust = 0.5) +
        # Title text
        annotate("text", x = xmin + 0.7, y = ymax - 0.28, label = title,
                 color = "#000000", fontface = "bold", size = 3.2, hjust = 0) +
        # Subtitle text
        annotate("text", x = xmin + 0.7, y = ymin + 0.24, label = subtitle,
                 color = "grey30", size = 2.8, hjust = 0)
    return(plot_obj)
}

p <- draw_band2_box(p, 0.4, 5.1, 6.4, 7.3, "1",
                    "Druggable-gene-wide Mendelian randomization",
                    "2,545 MR-testable genes - 13 BBJ Bonferroni hits")

p <- draw_band2_box(p, 5.3, 10.0, 6.4, 7.3, "2",
                    "Colocalization (shared causal variant)",
                    "shared vs distinct eQTL - GWAS signals")

p <- draw_band2_box(p, 0.4, 5.1, 5.2, 6.1, "3",
                    "Statistical fine-mapping (instrument structure)",
                    "credible sets & ancestry-matched LD")

p <- draw_band2_box(p, 5.3, 10.0, 5.2, 6.1, "4",
                    "Orbital expression (differential analysis)",
                    "exploratory tissue-level support")


# =============================================================================
# BAND 3: INTEGRATION
# =============================================================================

# Row 1: TSHR
p <- p +
    annotate("rect", xmin = 0.4, xmax = 10.0, ymin = 3.9, ymax = 4.8,
             fill = "#FFFFFF", color = "black", linewidth = 0.6) +
    # Accent bar on left
    annotate("rect", xmin = 0.4, xmax = 0.48, ymin = 3.9, ymax = 4.8,
             fill = "#111111", color = NA) +
    # TSHR susceptibility labels
    annotate("text", x = 0.65, y = 4.5, label = "TSHR", color = "#000000", fontface = "bold", size = 3.8, hjust = 0) +
    annotate("text", x = 0.65, y = 4.15, label = "susceptibility", color = "#555555", size = 3.0, hjust = 0) +
    # Vertical line separator inside card
    annotate("segment", x = 2.0, xend = 2.0, y = 3.95, yend = 4.75, color = "grey70", linewidth = 0.5) +
    # Right-side text
    annotate("text", x = 2.2, y = 4.5, label = "Genetically anchored susceptibility",
             color = "#000000", fontface = "bold", size = 3.3, hjust = 0) +
    annotate("text", x = 2.2, y = 4.15, label = "convergent MR, colocalization, and orbital expression evidence",
             color = "grey30", size = 3.0, hjust = 0)

# Row 2: IGF1R
p <- p +
    annotate("rect", xmin = 0.4, xmax = 10.0, ymin = 2.9, ymax = 3.8,
             fill = "#FFFFFF", color = "black", linewidth = 0.6) +
    # Accent bar on left
    annotate("rect", xmin = 0.4, xmax = 0.48, ymin = 2.9, ymax = 3.8,
             fill = "#666666", color = NA) +
    # IGF1R effector labels
    annotate("text", x = 0.65, y = 3.5, label = "IGF1R", color = "#000000", fontface = "bold", size = 3.8, hjust = 0) +
    annotate("text", x = 0.65, y = 3.15, label = "effector", color = "#555555", size = 3.0, hjust = 0) +
    # Vertical line separator inside card
    annotate("segment", x = 2.0, xend = 2.0, y = 2.95, yend = 3.75, color = "grey70", linewidth = 0.5) +
    # Right-side text
    annotate("text", x = 2.2, y = 3.5, label = "Pharmacologic effector axis",
             color = "#000000", fontface = "bold", size = 3.3, hjust = 0) +
    annotate("text", x = 2.2, y = 3.15, label = "MR risk-direction signal, but not expression-colocalized susceptibility",
             color = "grey30", size = 3.0, hjust = 0)

# Row 3: CTLA4
p <- p +
    annotate("rect", xmin = 0.4, xmax = 10.0, ymin = 2.2, ymax = 2.8,
             fill = "#FAFAFA", color = "black", linewidth = 0.6) +
    # Left text
    annotate("text", x = 0.65, y = 2.5, label = "CTLA4", color = "#000000", fontface = "bold", size = 3.8, hjust = 0) +
    # Vertical line separator inside card
    annotate("segment", x = 2.0, xend = 2.0, y = 2.25, yend = 2.75, color = "grey70", linewidth = 0.5) +
    # Right-side text
    annotate("text", x = 2.2, y = 2.5, label = "known autoimmune locus - included as a positive control",
             color = "grey30", size = 3.1, hjust = 0)

# Row 4: Filtering
p <- p +
    annotate("rect", xmin = 0.4, xmax = 10.0, ymin = 1.5, ymax = 2.1,
             fill = "#FAFAFA", color = "black", linewidth = 0.6) +
    # Left text
    annotate("text", x = 0.65, y = 1.8, label = "Filtering", color = "#000000", fontface = "bold", size = 3.8, hjust = 0) +
    # Vertical line separator inside card
    annotate("segment", x = 2.0, xend = 2.0, y = 1.55, yend = 2.05, color = "grey70", linewidth = 0.5) +
    # Right-side text
    annotate("text", x = 2.2, y = 1.8, label = "no robust novel target after candidate colocalization",
             color = "grey30", size = 3.1, hjust = 0)


# =============================================================================
# ARROWS (CONNECTIONS)
# =============================================================================

# ---- Band 1 to Band 2 arrows (vertical down) ----
# Arrow 1: x = 1.9, y = 8.1 to y = 7.35
p <- p + annotate("segment", x = 1.9, xend = 1.9, y = 8.1, yend = 7.35,
                  arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                  color = "grey20", linewidth = 0.5)

# Arrow 2: x = 5.2, y = 8.1 to y = 7.35
p <- p + annotate("segment", x = 5.2, xend = 5.2, y = 8.1, yend = 7.35,
                  arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                  color = "grey20", linewidth = 0.5)

# Arrow 3: x = 8.5, y = 8.1 to y = 7.35
p <- p + annotate("segment", x = 8.5, xend = 8.5, y = 8.1, yend = 7.35,
                  arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                  color = "grey20", linewidth = 0.5)

# ---- Band 2 to Band 3 converging arrows ----
# Step 3 and 4 bottom: y = 5.2. Band 3 top: y = 4.8.
p <- p +
    # Left vertical down
    annotate("segment", x = 2.75, xend = 2.75, y = 5.2, yend = 4.95, color = "grey20", linewidth = 0.5) +
    # Right vertical down
    annotate("segment", x = 7.65, xend = 7.65, y = 5.2, yend = 4.95, color = "grey20", linewidth = 0.5) +
    # Horizontal line
    annotate("segment", x = 2.75, xend = 7.65, y = 4.95, yend = 4.95, color = "grey20", linewidth = 0.5) +
    # Center vertical arrow
    annotate("segment", x = 5.2, xend = 5.2, y = 4.95, yend = 4.82,
              arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
              color = "grey20", linewidth = 0.5)


# =============================================================================
# BAND LABELS (LEFT SIDE)
# =============================================================================
p <- p +
    annotate("text", x = -0.5, y = 8.9, label = "DATA SOURCES",
             fontface = "bold", size = 3.6, color = "grey25", angle = 90) +
    annotate("text", x = -0.5, y = 6.25, label = "ANALYSIS",
             fontface = "bold", size = 3.6, color = "grey25", angle = 90) +
    annotate("text", x = -0.5, y = 3.15, label = "INTEGRATION",
             fontface = "bold", size = 3.6, color = "grey25", angle = 90)

# =============================================================================
# TITLE & SUBTITLE
# =============================================================================
p <- p +
    annotate("text", x = 5.2, y = 11.2,
             label = "Figure 1. Study design and analytical workflow",
             fontface = "bold", size = 5.0, color = "black") +
    annotate("text", x = 5.2, y = 10.8,
             label = "Multi-outcome mapping, locus colocalization, fine-mapping, and target-tissue transcriptomics",
             size = 3.2, color = "grey35", fontface = "italic")

# Bottom Subtitle
p <- p +
    annotate("text", x = 5.2, y = 0.8,
             label = "A druggable-gene-wide framework distinguishing genetically anchored susceptibility",
             size = 3.3, color = "black", fontface = "bold") +
    annotate("text", x = 5.2, y = 0.5,
             label = "from pharmacologic effector biology in Graves disease and thyroid eye disease",
             size = 3.3, color = "black", fontface = "bold")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
ggsave(file.path(outdir, "Figure1_study_design.pdf"), p,
    width = 11, height = 10.5, device = grDevices::cairo_pdf, bg = "white"
)
ggsave(file.path(outdir, "Figure1_study_design.png"), p,
    width = 11, height = 10.5, dpi = 600, bg = "white"
)
ggsave(file.path(outdir, "Figure1_study_design.tiff"), p,
    width = 11, height = 10.5, dpi = 600, compression = "lzw", bg = "white"
)

cat("Figure 1 study design schematic updated successfully to match Figure1_editable design.\n")
