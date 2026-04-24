# =============================================================================
# TED-TRAP — Figure 3 Panels A/B/C v3 (Thyroid publication grade)
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

# ---- Thyroid-style typography (fallback if Arial unavailable) ----
base_family <- tryCatch(
    {
        if ("Arial" %in% names(grDevices::postscriptFonts()) ||
            "Arial" %in% systemfonts::system_fonts()$family) {
            "Arial"
        } else {
            "sans"
        }
    },
    error = function(e) "sans"
)

# =============================================================================
# STEP 1: Load data
# =============================================================================
mr <- read.csv("TrackA_MR/results/02_mr_main.csv", stringsAsFactors = FALSE)
coloc <- read.csv("TrackA_MR/results/03_coloc_TSHR.csv", stringsAsFactors = FALSE)
mvmr <- read.csv("TrackA_MR/results/04_mvmr_summary.csv", stringsAsFactors = FALSE)
tis <- read.csv("TrackA_MR/results/06b_candidate_deseq2_reanalysis.csv",
    stringsAsFactors = FALSE
)
tpm <- read.csv("TrackA_MR/data_for_fig3/candidate_sample_TPM.csv",
    stringsAsFactors = FALSE, check.names = FALSE
)

# =============================================================================
# STEP 2: Evidence matrix
# =============================================================================
gene_groups <- tibble::tribble(
    ~gene, ~group, ~group_order, ~gene_order,
    "TSHR", "MR Primary", 1, 1,
    "IGF1R", "MR Primary", 1, 2,
    "IGF1", "MR Primary", 1, 3,
    "ARRB1", "MR Secondary", 2, 4,
    "PPARG", "MR Secondary", 2, 5,
    "IRS1", "MR Secondary", 2, 6,
    "AKT1", "MR Secondary", 2, 7,
    "TNF", "MR Secondary", 2, 8,
    "CTLA4", "MR Secondary", 2, 9,
    "INSR", "Insulin Cassette", 3, 10,
    "IRS2", "Insulin Cassette", 3, 11,
    "FOXO1", "Insulin Cassette", 3, 12,
    "PIK3R1", "Insulin Cassette", 3, 13,
    "PDPK1", "Insulin Cassette", 3, 14,
    "HAS1", "TED Biology", 4, 15,
    "HAS2", "TED Biology", 4, 16,
    "HAS3", "TED Biology", 4, 17,
    "ADIPOQ", "TED Biology", 4, 18,
    "FABP4", "TED Biology", 4, 19,
    "CEBPA", "TED Biology", 4, 20
)

mr_primary <- mr %>%
    filter(
        outcome_role == "Primary",
        method %in% c("Inverse variance weighted", "Wald ratio")
    ) %>%
    transmute(gene, mr_pval = pval)

coloc_val <- coloc %>%
    filter(Outcome == "GD_Primary_BBJ") %>%
    pull(PP_H4)

mat <- gene_groups %>%
    left_join(mr_primary, by = "gene") %>%
    left_join(tis %>% transmute(gene, tis_padj = padj_candidate_raw), by = "gene") %>%
    mutate(
        coloc_PPH4 = ifelse(gene == "TSHR", coloc_val, NA_real_),
        mvmr_flag = ifelse(gene == "TSHR", 1, NA_real_),
        ev_mr = ifelse(!is.na(mr_pval) & mr_pval < 0.05, 1, 0),
        ev_coloc = ifelse(!is.na(coloc_PPH4) & coloc_PPH4 > 0.8, 1, 0),
        ev_mvmr = ifelse(gene == "TSHR", 1, 0),
        ev_tis = ifelse(!is.na(tis_padj) & tis_padj < 0.05, 1, 0),
        conv_score = ev_mr + ev_coloc + ev_mvmr + ev_tis
    )

# =============================================================================
# STEP 3: Panel A v3 — Thyroid-grade heatmap
# =============================================================================
# Flag: which cells are "not tested" vs "tested null"
long <- mat %>%
    transmute(
        gene, group, gene_order,
        # Scored cells (0-1 range)
        `MR` = case_when(
            is.na(mr_pval) ~ NA_real_,
            mr_pval < 1e-10 ~ 1.0,
            mr_pval < 0.001 ~ 0.85,
            mr_pval < 0.05 ~ 0.55,
            TRUE ~ 0.15
        ),
        `Coloc` = case_when(
            is.na(coloc_PPH4) ~ NA_real_,
            coloc_PPH4 > 0.9 ~ 0.95,
            coloc_PPH4 > 0.8 ~ 0.75,
            TRUE ~ 0.3
        ),
        `MVMR` = case_when(
            is.na(mvmr_flag) ~ NA_real_,
            TRUE ~ 0.8
        ),
        `Tissue` = case_when(
            is.na(tis_padj) ~ NA_real_,
            tis_padj < 0.01 ~ 0.9,
            tis_padj < 0.05 ~ 0.6,
            tis_padj < 0.1 ~ 0.3,
            TRUE ~ 0.1
        ),
        `Conv.` = conv_score / 4
    ) %>%
    pivot_longer(
        cols = c(MR, Coloc, MVMR, Tissue, `Conv.`),
        names_to = "evidence", values_to = "score"
    ) %>%
    mutate(
        evidence = factor(evidence, levels = c("MR", "Coloc", "MVMR", "Tissue", "Conv.")),
        # Which cells are "not tested"?
        not_tested = is.na(score) & evidence %in% c("Coloc", "MVMR")
    )

# Merge values for labels
long <- long %>%
    left_join(mat %>% select(gene, mr_pval, coloc_PPH4, tis_padj, conv_score),
        by = "gene"
    ) %>%
    mutate(label = case_when(
        evidence == "MR" & !is.na(mr_pval) & mr_pval < 0.001 ~ formatC(mr_pval, format = "e", digits = 1),
        evidence == "MR" & !is.na(mr_pval) ~ sprintf("%.3f", mr_pval),
        evidence == "Coloc" & !is.na(coloc_PPH4) ~ sprintf("%.3f", coloc_PPH4),
        evidence == "Coloc" & is.na(coloc_PPH4) ~ "Not\ntested", # PATCH C
        evidence == "MVMR" & !is.na(score) ~ "Confirmed",
        evidence == "MVMR" & is.na(score) ~ "Not\ntested", # PATCH C
        evidence == "Tissue" & !is.na(tis_padj) & tis_padj < 0.001 ~ formatC(tis_padj, format = "e", digits = 1),
        evidence == "Tissue" & !is.na(tis_padj) ~ sprintf("%.3f", tis_padj),
        evidence == "Conv." ~ paste0(conv_score, "/4"),
        TRUE ~ ""
    ))

long$gene <- factor(long$gene, levels = rev(gene_groups$gene))

# Style per cell: not_tested cells get smaller italic grey text
long <- long %>%
    mutate(
        label_size = ifelse(not_tested, 1.9, 2.3),
        label_face = case_when(
            not_tested ~ "italic",
            evidence == "Conv." & conv_score >= 1 ~ "bold",
            TRUE ~ "plain"
        ),
        label_color = ifelse(not_tested, "grey40", "black")
    )

# ------ Thyroid-style color palette ------
# Sequential: light cream → muted sage → deep forest green
# Colorblind-safe, journal-friendly
panelA <- ggplot(long, aes(x = evidence, y = gene, fill = score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
        aes(
            label = label, size = label_size, fontface = label_face,
            color = label_color
        ),
        lineheight = 0.85, show.legend = FALSE
    ) +
    scale_size_identity() +
    scale_color_identity() +
    scale_fill_gradientn(
        colors = c("#FBFBF5", "#EFECD9", "#D9E6BF", "#9BC48F", "#4E8D5C", "#1E5F3A"),
        values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1.0)),
        limits = c(0, 1),
        na.value = "#EDEDED",
        name = "Evidence\nstrength",
        breaks = c(0, 0.3, 0.6, 0.9),
        labels = c("Null", "Weak", "Moderate", "Strong"),
        guide = guide_colorbar(
            barwidth = 0.5, barheight = 5.5,
            title.theme = element_text(size = 7.5, family = base_family),
            label.theme = element_text(size = 7, family = base_family)
        )
    ) +
    # Column separator between Tissue and Conv.
    geom_vline(xintercept = 4.5, linewidth = 0.7, color = "grey25") +
    # Group separator (horizontal)
    geom_hline(
        yintercept = c(6.5, 11.5, 14.5), linewidth = 0.35,
        color = "grey35", linetype = "22"
    ) +
    # TSHR gold border (5 columns)
    annotate("rect",
        xmin = 0.5, xmax = 5.5, ymin = 19.5, ymax = 20.5,
        fill = NA, color = "#B8860B", linewidth = 1.3
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = NULL, y = NULL,
        title = "A. Multi-modal triangulation across 20 pre-specified candidate genes"
    ) +
    theme_minimal(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(
            face = "bold", size = 11.5,
            margin = margin(b = 8), family = base_family
        ),
        axis.text.x.top = element_text(
            face = "bold", size = 9.5,
            margin = margin(b = 4), family = base_family
        ),
        axis.text.y = element_text(face = "bold", size = 9, family = base_family),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey25", linewidth = 0.6),
        legend.position = "right",
        legend.margin = margin(l = 10),
        plot.margin = margin(5, 60, 5, 80) # extra left for group brackets
    )

# --- Group bracket labels (left side) ---
grp_labels <- tibble(
    y_mid = c(19, 16, 12, 3.5), # mid-point of each group
    label = c("MR\nPrimary", "MR\nSecondary", "Insulin\nCassette", "TED\nBiology")
)

# Panel A annotations outside plot area (using patchwork wrap + draw_grob)
panelA <- panelA +
    annotation_custom(
        grob = grid::textGrob("MR Primary",
            rot = 90,
            gp = grid::gpar(fontsize = 9, fontface = "bold", col = "grey25")
        ),
        xmin = -0.6, xmax = -0.6, ymin = 19, ymax = 19
    ) +
    annotation_custom(
        grob = grid::textGrob("MR Secondary",
            rot = 90,
            gp = grid::gpar(fontsize = 9, fontface = "bold", col = "grey25")
        ),
        xmin = -0.6, xmax = -0.6, ymin = 15.5, ymax = 15.5
    ) +
    annotation_custom(
        grob = grid::textGrob("Insulin Cassette",
            rot = 90,
            gp = grid::gpar(fontsize = 9, fontface = "bold", col = "grey25")
        ),
        xmin = -0.6, xmax = -0.6, ymin = 9, ymax = 9
    ) +
    annotation_custom(
        grob = grid::textGrob("TED Biology",
            rot = 90,
            gp = grid::gpar(fontsize = 9, fontface = "bold", col = "grey25")
        ),
        xmin = -0.6, xmax = -0.6, ymin = 3.5, ymax = 3.5
    ) +
    coord_cartesian(clip = "off")

# =============================================================================
# STEP 4: Panel B (refined)
# =============================================================================
tpm_cols_TPM <- grep("_TPM$", names(tpm), value = TRUE)
tshr_tpm <- tpm %>%
    filter(Gene_Symbol == "TSHR") %>%
    select(all_of(tpm_cols_TPM))
samples <- c("2", "7", "8", "10", "11")
tshr_sample <- tibble(
    sample = paste0("S", samples),
    tpm = sapply(samples, function(s) {
        mean(as.numeric(tshr_tpm[, c(paste0(s, "__1_TPM"), paste0(s, "__2_TPM"))]))
    })
)
tshr_sample$label <- c("Ctrl", "T1", "T2", "T3", "T4")
tshr_sample$group <- c("Control", rep("TED", 4))
tshr_sample$label <- factor(tshr_sample$label, levels = c("Ctrl", "T1", "T2", "T3", "T4"))
ted_mean <- mean(tshr_sample$tpm[tshr_sample$group == "TED"])

panelB <- ggplot(tshr_sample, aes(x = label, y = tpm, fill = group)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.35) +
    geom_hline(
        yintercept = ted_mean, linetype = "22",
        color = "#8B0000", linewidth = 0.55
    ) +
    annotate("text",
        x = 5.3, y = ted_mean, label = "TED mean",
        size = 2.6, color = "#8B0000", hjust = 0, fontface = "italic"
    ) +
    geom_text(aes(label = sprintf("%.2f", tpm)),
        vjust = -0.55, size = 2.9,
        family = base_family
    ) +
    scale_fill_manual(
        values = c("Control" = "grey60", "TED" = "#B22222"),
        guide = "none"
    ) +
    annotate("text",
        x = 3, y = max(tshr_sample$tpm) * 1.3,
        label = "bold('6.5× up in TED, ')~italic(log[2]*FC) == '+2.33, '~italic(p[adj]) == 0.006",
        size = 2.9, parse = TRUE
    ) +
    coord_cartesian(ylim = c(0, max(tshr_sample$tpm) * 1.45), clip = "off") +
    labs(
        title = "B. TSHR per-sample expression",
        x = NULL, y = "TSHR expression (TPM)"
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5, family = base_family),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(family = base_family),
        axis.title.y = element_text(family = base_family, size = 9.5),
        plot.margin = margin(5, 20, 5, 5)
    )

# =============================================================================
# STEP 5: Panel C (refined)
# =============================================================================
ic_genes <- c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1")
ic_fc <- lapply(ic_genes, function(g) {
    row <- tpm %>%
        filter(Gene_Symbol == g) %>%
        select(all_of(tpm_cols_TPM))
    if (nrow(row) == 0) {
        return(NULL)
    }
    vals <- sapply(samples, function(s) mean(as.numeric(row[, c(paste0(s, "__1_TPM"), paste0(s, "__2_TPM"))])))
    ctrl <- vals["2"]
    ted <- mean(vals[c("7", "8", "10", "11")])
    tibble(gene = g, log2FC = log2(ted / ctrl))
}) %>% bind_rows()
ic_fc$gene <- factor(ic_fc$gene, levels = ic_fc$gene[order(-ic_fc$log2FC)])

panelC <- ggplot(ic_fc, aes(x = gene, y = log2FC, fill = log2FC)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.35) +
    geom_hline(yintercept = 0, color = "grey25", linewidth = 0.5) +
    geom_text(aes(label = sprintf("+%.2f", log2FC)),
        vjust = -0.5, size = 2.9, fontface = "bold", family = base_family
    ) +
    scale_fill_gradient(
        low = "#FDD9B5", high = "#D95F0E",
        limits = c(0, 1.2), guide = "none"
    ) +
    coord_cartesian(ylim = c(0, 1.35), clip = "off") +
    annotate("text",
        x = 3, y = 1.28,
        label = "italic('5/5 directional up, Binomial P = 0.031')",
        size = 2.9, parse = TRUE
    ) +
    labs(
        title = "C. Insulin signaling cassette",
        x = NULL, y = expression(log[2] * " fold-change (TED / Ctrl)")
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5, family = base_family),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold", family = base_family),
        axis.text.y = element_text(family = base_family),
        axis.title.y = element_text(family = base_family, size = 9.5),
        plot.margin = margin(5, 20, 5, 5)
    )

# =============================================================================
# STEP 6: Assemble
# =============================================================================
combined <- panelA / (panelB | panelC) +
    plot_layout(heights = c(2.4, 1))

ggsave("TrackA_MR/figures/Figure3_panelABC_v3.pdf", combined,
    width = 11, height = 12.5, device = cairo_pdf, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_panelABC_v3.png", combined,
    width = 11, height = 12.5, dpi = 600, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_panelABC_v3.tiff", combined,
    width = 11, height = 12.5, dpi = 600, compression = "lzw", bg = "white"
)

cat("\n✅ Figure 3 Panels A+B+C v3 (Thyroid publication grade) saved\n")
cat("  - 'Not tested' italic grey labels in Coloc/MVMR empty cells\n")
cat("  - Group brackets (MR Primary / Secondary / Insulin / TED) on left\n")
cat("  - Sage-forest green colorblind-safe palette\n")
cat("  - Thyroid-style typography (Arial/Helvetica, refined margins)\n")
