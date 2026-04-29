# =============================================================================
# TED-TRAP — Figure 3 Main v2 (polish patch)
# Changes: remove "TED mean" text overlap, remove ★ stars, tighter layout
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(scales)
    library(grid)
})

setwd("c:/ProjectTEDGWAS")
base_family <- "sans"

# ---- Load data ----
tpm <- read.csv("TrackA_MR/data_for_fig3/candidate_sample_TPM.csv",
    stringsAsFactors = FALSE, check.names = FALSE
)
raw <- read.csv("TrackA_MR/results/05_cross_phenotype_mr.csv", stringsAsFactors = FALSE)

# =============================================================================
# Panel A — TSHR per-sample (PATCH: remove "TED mean" label, add to subtitle)
# =============================================================================
tpm_cols_TPM <- grep("_TPM$", names(tpm), value = TRUE)
tshr_tpm <- tpm %>%
    filter(Gene_Symbol == "TSHR") %>%
    select(all_of(tpm_cols_TPM))
samples <- c("2", "7", "8", "10", "11")
tshr_sample <- tibble(
    tpm = sapply(samples, function(s) {
        mean(as.numeric(tshr_tpm[, c(paste0(s, "__1_TPM"), paste0(s, "__2_TPM"))]))
    }),
    label = c("Ctrl", "T1", "T2", "T3", "T4"),
    group = c("Control", rep("TED", 4))
)
tshr_sample$label <- factor(tshr_sample$label, levels = c("Ctrl", "T1", "T2", "T3", "T4"))
ted_mean <- mean(tshr_sample$tpm[tshr_sample$group == "TED"])

panelA <- ggplot(tshr_sample, aes(x = label, y = tpm, fill = group)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.35) +
    # PATCH: remove annotate text — keep dashed line only
    geom_hline(
        yintercept = ted_mean, linetype = "22",
        color = "#8B0000", linewidth = 0.55
    ) +
    geom_text(aes(label = sprintf("%.2f", tpm)), vjust = -0.55, size = 2.9) +
    scale_fill_manual(
        values = c("Control" = "grey60", "TED" = "#B22222"),
        guide = "none"
    ) +
    annotate("text",
        x = 3, y = max(tshr_sample$tpm) * 1.3,
        label = "6.5× up in TED,  log2FC = +2.33,  padj = 0.006",
        size = 2.9, fontface = "bold"
    ) +
    coord_cartesian(ylim = c(0, max(tshr_sample$tpm) * 1.45)) +
    labs(
        title = "A. TSHR per-sample expression in orbital tissue",
        subtitle = "Dashed red line indicates mean TED expression (excluding control).",
        x = NULL, y = "TSHR expression (TPM)"
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(
            size = 8.5, face = "italic", color = "grey30",
            margin = margin(b = 6)
        ),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )

# =============================================================================
# Panel B — Insulin cassette (unchanged)
# =============================================================================
ic_genes <- c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1")
ic_fc <- lapply(ic_genes, function(g) {
    row <- tpm %>%
        filter(Gene_Symbol == g) %>%
        select(all_of(tpm_cols_TPM))
    vals <- sapply(samples, function(s) mean(as.numeric(row[, c(paste0(s, "__1_TPM"), paste0(s, "__2_TPM"))])))
    tibble(gene = g, log2FC = log2(mean(vals[c("7", "8", "10", "11")]) / vals["2"]))
}) %>% bind_rows()
ic_fc$gene <- factor(ic_fc$gene, levels = ic_fc$gene[order(-ic_fc$log2FC)])

panelB <- ggplot(ic_fc, aes(x = gene, y = log2FC, fill = log2FC)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.35) +
    geom_hline(yintercept = 0, color = "grey25", linewidth = 0.5) +
    geom_text(aes(label = sprintf("+%.2f", log2FC)),
        vjust = -0.5, size = 2.9, fontface = "bold"
    ) +
    scale_fill_gradient(
        low = "#FDD9B5", high = "#D95F0E",
        limits = c(0, 1.2), guide = "none"
    ) +
    coord_cartesian(ylim = c(0, 1.35)) +
    annotate("text",
        x = 3, y = 1.28,
        label = "5/5 directional up, Binomial P = 0.031",
        size = 2.9, fontface = "italic"
    ) +
    labs(
        title = "B. Insulin signaling cassette in orbital tissue",
        subtitle = "Fold-change relative to single non-TED control; all five effectors upregulated.",
        x = NULL, y = "log2 fold-change (TED / Ctrl)"
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(
            size = 8.5, face = "italic", color = "grey30",
            margin = margin(b = 6)
        ),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold")
    )

# =============================================================================
# Panel C — Forest plot (PATCH: remove ★, keep gold bold formatting only)
# =============================================================================
ivw <- raw %>%
    filter(method == "Inverse variance weighted") %>%
    select(outcome = outcome_name, role = outcome_role, nsnp, beta = b, se, pval) %>%
    mutate(lci = beta - 1.96 * se, uci = beta + 1.96 * se)

ivw_plot <- ivw %>%
    filter(!grepl("skin|hearing", outcome, ignore.case = TRUE)) %>%
    mutate(
        Outcome = case_when(
            grepl("Fasting.?Glucose|glucose", outcome, ignore.case = TRUE) ~ "Fasting glucose",
            grepl("HbA1c", outcome, ignore.case = TRUE) ~ "HbA1c",
            grepl("Type.?2.?Diabetes|T2D|diabetes", outcome, ignore.case = TRUE) ~ "Type 2 diabetes",
            grepl("Height|height", outcome, ignore.case = TRUE) ~ "Height (IGF-axis anchor)",
            TRUE ~ outcome
        ),
        Category = case_when(
            grepl("Height", Outcome) ~ "Instrument-validity anchor",
            TRUE ~ "Metabolic off-target"
        ),
        label_pval = ifelse(pval < 0.001, formatC(pval, format = "e", digits = 1),
            sprintf("%.3f", pval)
        )
    )

# WM P for Height
wm_height_val <- raw %>%
    filter(
        method == "Weighted median",
        grepl("Height", outcome_name, ignore.case = TRUE)
    ) %>%
    pull(pval)
wm_height <- if (length(wm_height_val) > 0) formatC(wm_height_val, format = "e", digits = 2) else "NA"

ivw_plot$Outcome <- factor(ivw_plot$Outcome,
    levels = rev(c(
        "Fasting glucose", "HbA1c", "Type 2 diabetes",
        "Height (IGF-axis anchor)"
    ))
)

panelC <- ggplot(ivw_plot, aes(x = beta, y = Outcome, color = Category)) +
    geom_vline(xintercept = 0, linetype = "22", color = "grey40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.15, linewidth = 0.7) +
    geom_point(size = 3.2, shape = 18) +
    geom_text(
        aes(
            x = max(ivw_plot$uci) * 1.15,
            label = paste0("P = ", label_pval, " (n_IV = ", nsnp, ")")
        ),
        color = "black", hjust = 0, size = 2.9, family = base_family
    ) +
    # PATCH: remove ★, keep gold bold
    annotate("text",
        x = max(ivw_plot$uci) * 1.15, y = 1,
        label = paste0("WM P = ", wm_height, " (anchor)"),
        color = "#BF9000", hjust = 0, size = 2.9, fontface = "bold",
        vjust = 2.3
    ) +
    scale_color_manual(
        values = c(
            "Metabolic off-target" = "#4472C4",
            "Instrument-validity anchor" = "#BF9000"
        ),
        name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.55))) +
    labs(
        title = "C. IGF1R cross-phenotype Mendelian randomization",
        subtitle = "Germline IGF1R perturbation does not cause metabolic phenotypes (null IVW); Height instruments retain detectable pathway activity.",
        x = "IGF1R cis-eQTL causal effect on outcome (beta, 95% CI)",
        y = NULL
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(
            size = 8.5, face = "italic", color = "grey25",
            margin = margin(b = 6)
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.4, "cm"),
        axis.text.y = element_text(face = "bold")
    )

# =============================================================================
# Assemble — PATCH: tighter Panel C height
# =============================================================================
combined <- (panelA | panelB) / panelC +
    plot_layout(heights = c(1, 0.85)) # PATCH: was c(1, 1.1)

ggsave("TrackA_MR/figures/Figure3_main_3panel_v3.pdf", combined,
    width = 11, height = 8, device = cairo_pdf, bg = "white"
) # PATCH: height 9→8
ggsave("TrackA_MR/figures/Figure3_main_3panel_v3.png", combined,
    width = 11, height = 8, dpi = 600, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_main_3panel_v3.tiff", combined,
    width = 11, height = 8, dpi = 600, compression = "lzw", bg = "white"
)

cat("\n✅ Figure 3 Main v3 saved (Instrument-validity anchor applied)\n")
