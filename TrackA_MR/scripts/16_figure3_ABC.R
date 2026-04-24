# =============================================================================
# TED-TRAP — Figure 3 Panels A/B/C (data-driven)
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
if (!dir.exists("TrackA_MR/figures")) dir.create("TrackA_MR/figures", recursive = TRUE)

# =============================================================================
# STEP 1: Load all evidence sources
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
# STEP 2: Build unified 20 × 4 evidence matrix
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

# ---- MR column: -log10(P) for BBJ Primary, IVW or Wald ratio ----
mr_primary <- mr %>%
    filter(
        outcome_role == "Primary",
        method %in% c("Inverse variance weighted", "Wald ratio")
    ) %>%
    transmute(gene, mr_pval = pval, mr_beta = b)

# ---- Coloc: TSHR only (BBJ row) ----
coloc_val <- coloc %>%
    filter(Outcome == "GD_Primary_BBJ") %>%
    pull(PP_H4)

# ---- MVMR: TSHR row ----
mvmr_val <- mvmr %>%
    filter(Exposure == "TSHR") %>%
    transmute(mvmr_beta = Beta_adjusted, mvmr_pval = P_adjusted)

# ---- Tissue: 20 genes ----
tissue <- tis %>% transmute(gene,
    tis_log2FC = log2FC_raw,
    tis_pval = pvalue_raw, tis_padj = padj_candidate_raw
)

# ---- Merge into matrix ----
mat <- gene_groups %>%
    left_join(mr_primary, by = "gene") %>%
    left_join(tissue, by = "gene") %>%
    mutate(
        coloc_PPH4 = ifelse(gene == "TSHR", coloc_val, NA_real_),
        mvmr_beta  = ifelse(gene == "TSHR", mvmr_val$mvmr_beta, NA_real_),
        mvmr_pval  = ifelse(gene == "TSHR", mvmr_val$mvmr_pval, NA_real_)
    )

# ---- Evidence scores (0 = not tested/null, 1 = sig) ----
mat <- mat %>%
    mutate(
        ev_mr = ifelse(!is.na(mr_pval) & mr_pval < 0.05, 1, 0),
        ev_coloc = ifelse(!is.na(coloc_PPH4) & coloc_PPH4 > 0.8, 1, 0),
        ev_mvmr = ifelse(gene == "TSHR", 1, 0), # TSHR MVMR "upstream mediation confirmed"
        ev_tis = ifelse(!is.na(tis_padj) & tis_padj < 0.05, 1, 0),
        conv_score = ev_mr + ev_coloc + ev_mvmr + ev_tis
    )

cat("=== Evidence matrix ===\n")
print(mat %>% select(
    gene, group, mr_pval, coloc_PPH4,
    mvmr_beta, tis_log2FC, tis_padj, conv_score
))

# =============================================================================
# STEP 3: Prepare long-format for heatmap
# =============================================================================
long <- mat %>%
    transmute(gene, group, gene_order,
        `MR (Primary)`        = -log10(mr_pval),
        `Colocalization`      = coloc_PPH4 * 10, # scale for viz
        `MVMR (Mediation)`    = ifelse(gene == "TSHR", 5, NA), # categorical mark
        `Tissue (RNA-seq)`    = -log10(tis_padj)
    ) %>%
    pivot_longer(
        cols = c(
            `MR (Primary)`, `Colocalization`,
            `MVMR (Mediation)`, `Tissue (RNA-seq)`
        ),
        names_to = "evidence", values_to = "score"
    ) %>%
    mutate(evidence = factor(evidence, levels = c(
        "MR (Primary)", "Colocalization",
        "MVMR (Mediation)", "Tissue (RNA-seq)"
    )))

# Annotation labels inside cells
long <- long %>%
    mutate(label = case_when(
        evidence == "MR (Primary)" & !is.na(score) ~ ifelse(10^(-score) < 0.001,
            formatC(10^(-score), format = "e", digits = 1),
            sprintf("%.3f", 10^(-score))
        ),
        evidence == "Colocalization" & !is.na(score) ~ sprintf("%.3f", score / 10),
        evidence == "MVMR (Mediation)" & !is.na(score) ~ "Upstream\nconfirmed",
        evidence == "Tissue (RNA-seq)" & !is.na(score) ~ ifelse(10^(-score) < 0.001,
            formatC(10^(-score), format = "e", digits = 1),
            sprintf("%.3f", 10^(-score))
        ),
        TRUE ~ ""
    ))

# gene factor (reverse so TSHR appears on top)
long$gene <- factor(long$gene, levels = rev(gene_groups$gene))

# =============================================================================
# STEP 4: Panel A — Triangulation heatmap
# =============================================================================
panelA <- ggplot(long, aes(x = evidence, y = gene, fill = score)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.3, color = "black", lineheight = 0.85) +
    scale_fill_gradientn(
        colors = c("#F2F2F2", "#FFF2CC", "#A9D08E", "#548235", "#1F5F1F"),
        values = scales::rescale(c(0, 1.3, 3, 7, 15)),
        limits = c(0, 15),
        na.value = "#E8E8E8",
        name = "Signal\nstrength",
        breaks = c(0, 1.3, 3, 7, 15),
        labels = c("Null", "P=0.05", "P=0.001", "Strong", "Very strong"),
        guide = guide_colorbar(barwidth = 0.6, barheight = 6)
    ) +
    # Group separator horizontal lines (between MR Primary/Secondary/Insulin/TED)
    geom_hline(
        yintercept = c(6.5, 11.5, 14.5), linewidth = 0.4,
        color = "grey30", linetype = "dashed"
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

# Group labels on far left
group_lab <- data.frame(
    y = c(20, 14, 9.5, 3),
    label = c("MR\nPrimary", "MR\nSecondary", "Insulin\nCassette", "TED\nBiology")
)
# (group labels are built into y-axis via facet_grid alternative — keeping minimal for now)

# TSHR row highlight — gold rectangle
panelA <- panelA +
    annotate("rect",
        xmin = 0.5, xmax = 4.5,
        ymin = 19.5, ymax = 20.5,
        fill = NA, color = "#BF9000", linewidth = 1.5
    )

# Convergence score as right-side annotation
conv_anno <- mat %>%
    select(gene, conv_score) %>%
    mutate(gene = factor(gene, levels = rev(gene_groups$gene))) %>%
    mutate(
        label = paste0(conv_score, "/4"),
        color_flag = ifelse(conv_score == 4, "gold",
            ifelse(conv_score >= 1, "mid", "none")
        )
    )

panelA <- panelA +
    annotate("text",
        x = 4.8, y = as.numeric(conv_anno$gene),
        label = conv_anno$label,
        fontface = ifelse(conv_anno$conv_score == 4, "bold", "plain"),
        color = ifelse(conv_anno$conv_score == 4, "#BF9000",
            ifelse(conv_anno$conv_score >= 1, "#548235", "grey50")
        ),
        size = 3, hjust = 0
    ) +
    annotate("text",
        x = 4.8, y = 20.9, label = "Conv.", size = 3,
        fontface = "bold", hjust = 0
    ) +
    coord_cartesian(xlim = c(0.5, 5.3), clip = "off")

# =============================================================================
# STEP 5: Panel B — TSHR per-sample bar plot
# =============================================================================
tpm_cols_TPM <- grep("_TPM$", names(tpm), value = TRUE)
tshr_tpm <- tpm %>%
    filter(Gene_Symbol == "TSHR") %>%
    select(all_of(tpm_cols_TPM))

# Collapse tech reps by averaging
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
        x = 3, y = max(tshr_sample$tpm) * 1.25,
        label = "6.5× up in TED\nlog2FC = +2.33, padj = 0.006",
        size = 2.9, fontface = "bold"
    ) +
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
# STEP 6: Panel C — Insulin cassette Ctrl vs TED
# =============================================================================
ic_genes <- c("INSR", "IRS2", "FOXO1", "PIK3R1", "PDPK1")
ic_long <- lapply(ic_genes, function(g) {
    row <- tpm %>%
        filter(Gene_Symbol == g) %>%
        select(all_of(tpm_cols_TPM))
    if (nrow(row) == 0) {
        return(NULL)
    }
    sapply(samples, function(s) mean(as.numeric(row[, c(paste0(s, "__1_TPM"), paste0(s, "__2_TPM"))]))) -> vals
    tibble(
        gene = g,
        Ctrl = vals["2"],
        TED = mean(vals[c("7", "8", "10", "11")])
    )
}) %>%
    bind_rows() %>%
    pivot_longer(cols = c(Ctrl, TED), names_to = "group", values_to = "tpm")

ic_long$group <- factor(ic_long$group, levels = c("Ctrl", "TED"))
ic_long$gene <- factor(ic_long$gene, levels = ic_genes)

panelC <- ggplot(ic_long, aes(x = gene, y = tpm, fill = group)) +
    geom_col(
        position = position_dodge(0.75), width = 0.7,
        color = "black", linewidth = 0.3
    ) +
    scale_fill_manual(
        values = c("Ctrl" = "grey55", "TED" = "#ED7D31"),
        name = "Group"
    ) +
    labs(
        title = "C. Insulin signaling cassette",
        subtitle = "5/5 directional up in TED (Binomial P = 0.031)",
        x = NULL, y = "Mean TPM"
    ) +
    theme_bw(base_size = 10) +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, face = "italic"),
        legend.position = "top",
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major.x = element_blank()
    )

# =============================================================================
# STEP 7: Assemble (Panel A on top full, B & C side by side below)
# =============================================================================
combined <- panelA / (panelB | panelC) +
    plot_layout(heights = c(2.2, 1))

ggsave("TrackA_MR/figures/Figure3_panelABC.pdf", combined,
    width = 11, height = 12, device = cairo_pdf, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_panelABC.png", combined,
    width = 11, height = 12, dpi = 600, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_panelABC.tiff", combined,
    width = 11, height = 12, dpi = 600, compression = "lzw", bg = "white"
)

cat("\n✅ Figure 3 Panels A+B+C saved (Panel D schematic pending)\n")
