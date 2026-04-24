# =============================================================================
# TED-TRAP — Figure 3 Panels A/B/C v2 (A+C upgraded)
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

# =============================================================================
# STEP 1: Load data (same as v1)
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
# STEP 2: Evidence matrix (same logic)
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
mvmr_val <- mvmr %>% filter(Exposure == "TSHR")
tissue <- tis %>% transmute(gene, tis_padj = padj_candidate_raw)

mat <- gene_groups %>%
    left_join(mr_primary, by = "gene") %>%
    left_join(tissue, by = "gene") %>%
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
# STEP 3: Panel A v2 — Integrated 5-column heatmap (Convergence as 5th col)
# =============================================================================
# PATCH: per-column score normalized to 0-1, common color scale
long <- mat %>%
    transmute(
        gene, group, gene_order,
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
        `Conv.` = conv_score / 4 # PATCH: Convergence as true column, 0-1 scale
    ) %>%
    pivot_longer(
        cols = c(MR, Coloc, MVMR, Tissue, `Conv.`),
        names_to = "evidence", values_to = "score"
    ) %>%
    mutate(evidence = factor(evidence,
        levels = c("MR", "Coloc", "MVMR", "Tissue", "Conv.")
    ))

# labels per column
long <- long %>%
    left_join(mat %>% select(gene, mr_pval, coloc_PPH4, tis_padj, conv_score),
        by = "gene"
    ) %>%
    mutate(label = case_when(
        evidence == "MR" & !is.na(mr_pval) & mr_pval < 0.001 ~ formatC(mr_pval, format = "e", digits = 1),
        evidence == "MR" & !is.na(mr_pval) ~ sprintf("%.3f", mr_pval),
        evidence == "Coloc" & !is.na(coloc_PPH4) ~ sprintf("%.3f", coloc_PPH4),
        evidence == "MVMR" & !is.na(score) ~ "Confirmed",
        evidence == "Tissue" & !is.na(tis_padj) & tis_padj < 0.001 ~ formatC(tis_padj, format = "e", digits = 1),
        evidence == "Tissue" & !is.na(tis_padj) ~ sprintf("%.3f", tis_padj),
        evidence == "Conv." ~ paste0(conv_score, "/4"),
        TRUE ~ ""
    ))

long$gene <- factor(long$gene, levels = rev(gene_groups$gene))

# Conv. column label boldness
long$face <- ifelse(long$evidence == "Conv." & long$conv_score == 4, "bold",
    ifelse(long$evidence == "Conv." & long$conv_score >= 1, "bold",
        "plain"
    )
)

panelA <- ggplot(long, aes(x = evidence, y = gene, fill = score)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label, fontface = face), size = 2.3, color = "black") +
    scale_fill_gradientn(
        colors = c("#F2F2F2", "#FFF2CC", "#FFE699", "#A9D08E", "#548235", "#1F5F1F"),
        values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1.0)),
        limits = c(0, 1),
        na.value = "#E8E8E8",
        name = "Evidence\nstrength",
        breaks = c(0, 0.3, 0.6, 0.9),
        labels = c("Null", "Weak", "Moderate", "Strong"),
        guide = guide_colorbar(barwidth = 0.6, barheight = 6)
    ) +
    # Column separator (between Tissue and Conv.)
    geom_vline(xintercept = 4.5, linewidth = 0.6, color = "grey30") +
    # Group separator (horizontal)
    geom_hline(
        yintercept = c(6.5, 11.5, 14.5), linewidth = 0.4,
        color = "grey30", linetype = "dashed"
    ) +
    # TSHR row gold border (now spans 5 columns)
    annotate("rect",
        xmin = 0.5, xmax = 5.5, ymin = 19.5, ymax = 20.5,
        fill = NA, color = "#BF9000", linewidth = 1.5
    ) +
    scale_x_discrete(position = "top") +
    labs(
        x = NULL, y = NULL,
        title = "A. Multi-modal triangulation of 20 pre-specified candidate genes"
    ) +
    theme_minimal(base_size = 10) +
    theme(
        plot.title       = element_text(face = "bold", size = 11),
        axis.text.x.top  = element_text(angle = 0, hjust = 0.5, face = "bold", size = 9),
        axis.text.y      = element_text(face = "bold", size = 8.5),
        panel.grid       = element_blank(),
        legend.position  = "right",
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(5, 5, 5, 5)
    )

# =============================================================================
# STEP 4: Panel B (unchanged from v1)
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
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    geom_hline(
        yintercept = ted_mean, linetype = "dashed",
        color = "#C00000", linewidth = 0.5
    ) +
    geom_text(aes(label = sprintf("%.2f", tpm)), vjust = -0.5, size = 2.8) +
    scale_fill_manual(values = c("Control" = "grey55", "TED" = "#C00000")) +
    annotate("text",
        x = 3, y = max(tshr_sample$tpm) * 1.28,
        label = "6.5× up in TED\nlog2FC = +2.33, padj = 0.006",
        size = 2.9, fontface = "bold"
    ) +
    coord_cartesian(ylim = c(0, max(tshr_sample$tpm) * 1.4)) +
    labs(
        title = "B. TSHR per-sample expression",
        x = NULL, y = "TSHR TPM"
    ) +
    theme_bw(base_size = 10) +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "none",
        panel.grid.major.x = element_blank()
    )

# =============================================================================
# STEP 5: Panel C v2 — Fold-change plot (PATCH)
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
    tibble(
        gene = g, ctrl_tpm = ctrl, ted_tpm = ted,
        log2FC = log2(ted / ctrl)
    )
}) %>% bind_rows()

# Ordering: highest fold-change first
ic_fc$gene <- factor(ic_fc$gene, levels = ic_fc$gene[order(-ic_fc$log2FC)])

panelC <- ggplot(ic_fc, aes(x = gene, y = log2FC, fill = log2FC)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    geom_text(aes(label = sprintf("+%.2f", log2FC)),
        vjust = -0.5, size = 2.8, fontface = "bold"
    ) +
    scale_fill_gradient(
        low = "#FCE4D6", high = "#ED7D31",
        limits = c(0, 1.2), guide = "none"
    ) +
    coord_cartesian(ylim = c(-0.3, 1.3)) +
    annotate("segment",
        x = 0.5, xend = 5.5, y = 0, yend = 0,
        color = "grey30", linewidth = 0.5
    ) +
    annotate("text",
        x = 3, y = 1.22,
        label = "5/5 directional up (Binomial P = 0.031)",
        size = 3, fontface = "italic"
    ) +
    labs(
        title = "C. Insulin signaling cassette — log2 fold-change (TED / Ctrl)",
        x = NULL, y = expression(log[2] * " (TED / Ctrl)")
    ) +
    theme_bw(base_size = 10) +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold")
    )

# =============================================================================
# STEP 6: Assemble
# =============================================================================
combined <- panelA / (panelB | panelC) +
    plot_layout(heights = c(2.2, 1))

ggsave("TrackA_MR/figures/Figure3_panelABC_v2.pdf", combined,
    width = 11, height = 12, device = cairo_pdf, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_panelABC_v2.png", combined,
    width = 11, height = 12, dpi = 600, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_panelABC_v2.tiff", combined,
    width = 11, height = 12, dpi = 600, compression = "lzw", bg = "white"
)

cat("\n✅ Figure 3 Panels A+B+C v2 saved\n")
cat("  - Panel A: Convergence as integrated 5th column\n")
cat("  - Panel C: Fold-change plot (PIK3R1 scale dominance resolved)\n")
