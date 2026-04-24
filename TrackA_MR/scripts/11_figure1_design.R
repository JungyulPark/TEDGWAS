# =============================================================================
# TED-TRAP — Figure 1 PATCH v4
# Fix: Remove "SE" jargon, use descriptive analysis-based labels
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(grid)
})

if (!dir.exists("TrackA_MR/figures")) dir.create("TrackA_MR/figures", recursive = TRUE)

col_genetic <- "#2E75B6"
col_pathway <- "#ED7D31"
col_tissue <- "#7030A0"
col_method <- "#70AD47"
col_result <- "#D9E1F2"
col_integ <- "#A87900"

# BAND 1 — INPUT (unchanged)
band1 <- data.frame(
    x = c(1.5, 5, 8.5), y = 9, w = 2.8, h = 1.4,
    fill = c(col_genetic, col_pathway, col_tissue),
    title = c("GENETIC DATA", "CANDIDATE GENES", "TRANSCRIPTOMIC DATA"),
    content = c(
        "eQTLGen (n = 31,684)\nBBJ Graves (2,809 cases)\nUKB Hyperthyroid (3,731)\nFinnGen R12 GO (858)",
        "Pre-specified 20 genes:\nMR Primary (3)\nMR Secondary (6)\nInsulin Cassette (5)\nTED Biology (6)",
        "Bulk RNA-seq\n4 inactive TED\n+ 1 Control\nKorean cohort\nIRB-approved"
    ),
    stringsAsFactors = FALSE
)

# BAND 2 — METHODS (label naming polished)
band2 <- data.frame(
    x = c(1, 3.67, 6.33, 9), y = 6.5, w = 2.1, h = 1, fill = col_method,
    label = c(
        "Cis-MR", # PATCH: simplified
        "Bayesian\nColocalization",
        "Multivariable\nMR (MVMR)",
        "Candidate\nDESeq2"
    ),
    stringsAsFactors = FALSE
)

# BAND 3 — RESULTS (PATCH: descriptive labels, matching method order)
band3 <- data.frame(
    x = c(1, 3.67, 6.33, 9), y = 4, w = 2.1, h = 1.4, fill = col_result,
    title = c(
        "CAUSAL MR", # was: PRIMARY EP
        "COLOCALIZATION", # was: SE4: Coloc
        "MEDIATION (MVMR)", # was: SE1: MVMR
        "TISSUE EXPRESSION"
    ), # was: SE3: Tissue
    content = c(
        "TSHR P=9.1e-15\nIGF1R P=0.089\nCross-ancestry\nreplication", # PATCH: "concordance" → "replication"
        "BBJ PP.H4=0.951\nrs179252 shared\nLD proxy (EAS)\nr²=0.85",
        "TSHR β_adj=+0.13\nGD liability as\nupstream mediator", # PATCH: simplified
        "TSHR log2FC=+2.33\npadj=0.006\nInsulin cassette\n5/5 up (P=0.031)"
    ),
    stringsAsFactors = FALSE
)

# BAND 4 — INTEGRATION (unchanged content, wording slightly tightened)
band4 <- data.frame(
    x = 5, y = 1.3, w = 9, h = 1.3, fill = col_integ,
    title = "MULTI-MODAL INTEGRATION", # PATCH: "TRIANGULATION" → "INTEGRATION" for consistency
    content = "Convergent evidence for TSHR in TED — Population causality + Locus coloc + Orbital tissue",
    stringsAsFactors = FALSE
)

# =============================================================================
# Build plot
# =============================================================================

p <- ggplot() +
    theme_void() +
    xlim(-0.8, 10.8) +
    ylim(-0.2, 11.5)

# ---- BAND 1 ----
for (i in seq_len(nrow(band1))) {
    p <- p +
        annotate("rect",
            xmin = band1$x[i] - band1$w[i] / 2, xmax = band1$x[i] + band1$w[i] / 2,
            ymin = band1$y[i] - band1$h[i] / 2, ymax = band1$y[i] + band1$h[i] / 2,
            fill = band1$fill[i], alpha = 1.0, color = "black", linewidth = 0.6
        ) +
        annotate("text",
            x = band1$x[i], y = band1$y[i] + 0.5,
            label = band1$title[i], fontface = "bold", size = 3.5, color = "white"
        ) +
        annotate("text",
            x = band1$x[i], y = band1$y[i] - 0.15,
            label = band1$content[i], size = 3.0, color = "white",
            lineheight = 1.0, fontface = "bold"
        )
}

# ---- Band 1 → Band 2 arrows ----
input_method_edges <- data.frame(
    from_x = c(1.5, 1.5, 5, 5, 5, 8.5),
    from_y = rep(8.3, 6),
    to_x   = c(1, 3.67, 1, 3.67, 6.33, 9),
    to_y   = rep(7.0, 6)
)
for (i in seq_len(nrow(input_method_edges))) {
    p <- p + annotate("segment",
        x = input_method_edges$from_x[i],
        xend = input_method_edges$to_x[i],
        y = input_method_edges$from_y[i],
        yend = input_method_edges$to_y[i],
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        color = "grey20", linewidth = 0.5, alpha = 1.0
    )
}

# ---- BAND 2 ----
for (i in seq_len(nrow(band2))) {
    p <- p +
        annotate("rect",
            xmin = band2$x[i] - band2$w[i] / 2, xmax = band2$x[i] + band2$w[i] / 2,
            ymin = band2$y[i] - band2$h[i] / 2, ymax = band2$y[i] + band2$h[i] / 2,
            fill = band2$fill[i], alpha = 1.0, color = "black", linewidth = 0.6
        ) +
        annotate("text",
            x = band2$x[i], y = band2$y[i],
            label = band2$label[i], fontface = "bold", size = 3.3, color = "white",
            lineheight = 1.0
        )
}

# ---- Band 2 → Band 3 arrows ----
for (i in seq_len(nrow(band2))) {
    p <- p + annotate("segment",
        x = band2$x[i], xend = band3$x[i],
        y = 6, yend = 4.75,
        arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
        color = "grey20", linewidth = 0.55
    )
}

# ---- BAND 3 ----
for (i in seq_len(nrow(band3))) {
    p <- p +
        annotate("rect",
            xmin = band3$x[i] - band3$w[i] / 2, xmax = band3$x[i] + band3$w[i] / 2,
            ymin = band3$y[i] - band3$h[i] / 2, ymax = band3$y[i] + band3$h[i] / 2,
            fill = band3$fill[i], color = "black", linewidth = 0.6
        ) +
        annotate("text",
            x = band3$x[i], y = band3$y[i] + 0.52,
            label = band3$title[i], fontface = "bold", size = 3.2, color = "#1F4E79"
        ) + # PATCH: slightly smaller to fit longer labels
        annotate("text",
            x = band3$x[i], y = band3$y[i] - 0.15,
            label = band3$content[i], size = 2.8, color = "black",
            lineheight = 1.0
        )
}

# ---- Band 3 → Band 4 arrows ----
for (i in seq_len(nrow(band3))) {
    p <- p + annotate("segment",
        x = band3$x[i], xend = band4$x,
        y = 3.3, yend = 2,
        arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
        color = "grey20", linewidth = 0.6
    )
}

# ---- BAND 4 (integration hero) ----
p <- p +
    annotate("rect",
        xmin = band4$x - band4$w / 2, xmax = band4$x + band4$w / 2,
        ymin = band4$y - band4$h / 2, ymax = band4$y + band4$h / 2,
        fill = band4$fill, alpha = 1.0, color = "black", linewidth = 0.9
    ) +
    annotate("text",
        x = band4$x, y = band4$y + 0.35,
        label = band4$title, fontface = "bold", size = 4.8, color = "white"
    ) +
    annotate("text",
        x = band4$x, y = band4$y - 0.2,
        label = band4$content, size = 3.4, color = "white",
        fontface = "bold", lineheight = 1.0
    )

# ---- Band labels (left) ----
p <- p +
    annotate("text",
        x = -0.65, y = 9, label = "INPUT",
        fontface = "bold", size = 3.8, color = "grey15", angle = 90
    ) +
    annotate("text",
        x = -0.65, y = 6.5, label = "ANALYSIS",
        fontface = "bold", size = 3.8, color = "grey15", angle = 90
    ) +
    annotate("text",
        x = -0.65, y = 4, label = "RESULTS",
        fontface = "bold", size = 3.8, color = "grey15", angle = 90
    ) +
    annotate("text",
        x = -0.65, y = 1.3, label = "INTEGRATION",
        fontface = "bold", size = 3.8, color = "grey15", angle = 90
    )

# ---- Title ----
p <- p +
    annotate("text",
        x = 5, y = 11.2,
        label = "Figure 1. Study design and analytical workflow",
        fontface = "bold", size = 5.5, color = "black"
    ) +
    annotate("text",
        x = 5, y = 10.7,
        label = "Multi-modal causal inference for TSHR vs IGF1R in thyroid eye disease",
        size = 3.4, color = "grey30", fontface = "italic"
    )

# =============================================================================
# Save
# =============================================================================

ggsave("TrackA_MR/figures/Figure1_study_design.pdf", p,
    width = 12, height = 10, device = grDevices::cairo_pdf, bg = "white"
)
ggsave("TrackA_MR/figures/Figure1_study_design.png", p,
    width = 12, height = 10, dpi = 600, bg = "white"
)
ggsave("TrackA_MR/figures/Figure1_study_design.tiff", p,
    width = 12, height = 10, dpi = 600, compression = "lzw", bg = "white"
)

cat("✅ Figure 1 PATCH v4 saved (SE → descriptive labels, flow order)\n")
