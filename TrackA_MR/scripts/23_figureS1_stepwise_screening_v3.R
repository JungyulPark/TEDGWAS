# =============================================================================
# TED-TRAP — Supplementary Figure S1 FINAL (5-fix integrated patch)
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
    library(grid)
})

setwd("c:/ProjectTEDGWAS")
base_family <- "sans"

# ---- Load ----
mr <- read.csv("TrackA_MR/results/02_mr_main.csv", stringsAsFactors = FALSE)
coloc <- read.csv("TrackA_MR/results/03_coloc_TSHR.csv", stringsAsFactors = FALSE)
mvmr <- read.csv("TrackA_MR/results/04_mvmr_summary.csv", stringsAsFactors = FALSE)
tis <- read.csv("TrackA_MR/results/06b_candidate_deseq2_reanalysis.csv",
    stringsAsFactors = FALSE
)

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

# PATCH: Extract TSHR MVMR P-value (actual value from results)
mvmr_tshr_p <- mvmr %>%
    filter(Exposure == "TSHR") %>%
    pull(P_adjusted)
cat("TSHR MVMR P_adjusted =", mvmr_tshr_p, "\n")

# PATCH: Triangulation score re-computed — MVMR null does NOT count
mat <- gene_groups %>%
    left_join(mr_primary, by = "gene") %>%
    left_join(tis %>% transmute(gene, tis_padj = padj_candidate_raw), by = "gene") %>%
    mutate(
        coloc_PPH4 = ifelse(gene == "TSHR", coloc_val, NA_real_),
        mvmr_pval = ifelse(gene == "TSHR", mvmr_tshr_p, NA_real_),
        ev_mr = ifelse(!is.na(mr_pval) & mr_pval < 0.05, 1, 0),
        ev_coloc = ifelse(!is.na(coloc_PPH4) & coloc_PPH4 > 0.8, 1, 0),
        # PATCH: MVMR significance criterion — P < 0.05 (null TSHR fails this)
        ev_mvmr = ifelse(!is.na(mvmr_pval) & mvmr_pval < 0.05, 1, 0),
        ev_tis = ifelse(!is.na(tis_padj) & tis_padj < 0.05, 1, 0),
        triang_score = ev_mr + ev_coloc + ev_mvmr + ev_tis
    )

cat("\n=== Triangulation scores ===\n")
print(mat %>% select(gene, mr_pval, coloc_PPH4, mvmr_pval, tis_padj, triang_score))
cat("TSHR score:", mat$triang_score[mat$gene == "TSHR"], "/ 4\n")

# ---- Long format ----
long <- mat %>%
    transmute(
        gene, group, gene_order,
        # PATCH: Column keys now include metric label
        `MR (P)` = case_when(
            is.na(mr_pval) ~ NA_real_,
            mr_pval < 1e-10 ~ 1.0,
            mr_pval < 0.001 ~ 0.85,
            mr_pval < 0.05 ~ 0.55,
            TRUE ~ 0.15
        ),
        `Coloc (PP.H4)` = case_when(
            is.na(coloc_PPH4) ~ NA_real_,
            coloc_PPH4 > 0.9 ~ 0.95,
            coloc_PPH4 > 0.8 ~ 0.75,
            TRUE ~ 0.3
        ),
        # PATCH: MVMR now scored on P-value, not flag
        `MVMR (P)` = case_when(
            is.na(mvmr_pval) ~ NA_real_,
            mvmr_pval < 0.001 ~ 0.85,
            mvmr_pval < 0.05 ~ 0.55,
            TRUE ~ 0.1 # PATCH: null → pale yellow
        ),
        `Tissue (p_adj)` = case_when(
            is.na(tis_padj) ~ NA_real_,
            tis_padj < 0.01 ~ 0.9,
            tis_padj < 0.05 ~ 0.6,
            tis_padj < 0.1 ~ 0.3,
            TRUE ~ 0.1
        ),
        # PATCH: column renamed from Conv. to Triangulation
        `Triangulation` = triang_score / 4
    ) %>%
    pivot_longer(
        cols = c(
            `MR (P)`, `Coloc (PP.H4)`, `MVMR (P)`,
            `Tissue (p_adj)`, `Triangulation`
        ),
        names_to = "evidence", values_to = "score"
    ) %>%
    mutate(evidence = factor(evidence,
        levels = c(
            "MR (P)", "Coloc (PP.H4)",
            "MVMR (P)", "Tissue (p_adj)",
            "Triangulation"
        )
    ))

# ---- Labels ----
long <- long %>%
    left_join(mat %>% select(gene, mr_pval, coloc_PPH4, mvmr_pval, tis_padj, triang_score),
        by = "gene"
    ) %>%
    mutate(label = case_when(
        evidence == "MR (P)" & !is.na(mr_pval) & mr_pval < 0.001 ~ formatC(mr_pval, format = "e", digits = 1),
        evidence == "MR (P)" & !is.na(mr_pval) ~ sprintf("%.3f", mr_pval),
        evidence == "Coloc (PP.H4)" & !is.na(coloc_PPH4) ~ sprintf("%.3f", coloc_PPH4),
        # PATCH: MVMR now shows actual P-value
        evidence == "MVMR (P)" & !is.na(mvmr_pval) & mvmr_pval < 0.001 ~ formatC(mvmr_pval, format = "e", digits = 1),
        evidence == "MVMR (P)" & !is.na(mvmr_pval) ~ sprintf("%.3f", mvmr_pval),
        evidence == "Tissue (p_adj)" & !is.na(tis_padj) & tis_padj < 0.001 ~ formatC(tis_padj, format = "e", digits = 1),
        evidence == "Tissue (p_adj)" & !is.na(tis_padj) ~ sprintf("%.3f", tis_padj),
        evidence == "Triangulation" ~ paste0(triang_score, "/4"),
        TRUE ~ ""
    ))

long$gene <- factor(long$gene, levels = rev(gene_groups$gene))
long$face <- case_when(
    long$evidence == "Triangulation" & long$triang_score >= 1 ~ "bold",
    TRUE ~ "plain"
)

# ---- Subtitle ----
subtitle_text <- paste0(
    "Colocalization and MVMR were performed only on genes meeting the primary MR significance threshold (P < 0.05).\n",
    "Grey cells denote conditional non-testing rather than null results. Triangulation score counts layers with nominal significance."
)

# =============================================================================
# Build plot
# =============================================================================
suppS1 <- ggplot(long, aes(x = evidence, y = gene, fill = score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label, fontface = face), size = 2.4) +
    geom_vline(xintercept = 4.5, linewidth = 0.7, color = "grey25") +
    geom_hline(
        yintercept = c(6.5, 11.5, 14.5), linewidth = 0.35,
        color = "grey35", linetype = "22"
    ) +
    # TSHR gold border (on top)
    annotate("rect",
        xmin = 0.5, xmax = 5.5, ymin = 19.5, ymax = 20.5,
        fill = NA, color = "#B8860B", linewidth = 1.5
    ) +
    scale_fill_gradientn(
        colors = c("#FBFBF5", "#EFECD9", "#D9E6BF", "#9BC48F", "#4E8D5C", "#1E5F3A"),
        values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1.0)),
        limits = c(0, 1),
        na.value = "#EDEDED",
        name = "Evidence\nstrength",
        breaks = c(0, 0.3, 0.6, 0.9),
        labels = c("Null", "Weak", "Moderate", "Strong"),
        guide = guide_colorbar(
            barwidth = 0.55, barheight = 6,
            title.theme = element_text(size = 8),
            label.theme = element_text(size = 7.5)
        )
    ) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = NULL, y = NULL,
        title = "Supplementary Figure S1. Stepwise multi-modal screening across 20 pre-specified candidate genes",
        subtitle = subtitle_text
    ) +
    theme_minimal(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5, margin = margin(b = 3)),
        plot.subtitle = element_text(
            size = 9, face = "italic", color = "grey30",
            margin = margin(b = 10), lineheight = 1.15
        ),
        # PATCH: column header — slightly smaller to fit "Coloc (PP.H4)" etc.
        axis.text.x.top = element_text(face = "bold", size = 8.8, margin = margin(b = 4)),
        axis.text.y = element_text(face = "bold", size = 9),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey25", linewidth = 0.6),
        legend.position = "right",
        legend.margin = margin(l = 10),
        # PATCH: larger left margin to fit group brackets fully
        plot.margin = margin(5, 60, 5, 140)
    )

# PATCH: Group brackets — pulled closer to 0.1 to avoid truncation
suppS1 <- suppS1 +
    annotation_custom(
        grob = grid::textGrob("MR Primary",
            rot = 90,
            gp = grid::gpar(fontsize = 10, fontface = "bold", col = "grey20")
        ),
        xmin = -0.3, xmax = -0.3, ymin = 19, ymax = 19
    ) +
    annotation_custom(
        grob = grid::textGrob("MR Secondary",
            rot = 90,
            gp = grid::gpar(fontsize = 10, fontface = "bold", col = "grey20")
        ),
        xmin = -0.3, xmax = -0.3, ymin = 15.5, ymax = 15.5
    ) +
    annotation_custom(
        grob = grid::textGrob("Insulin Cassette",
            rot = 90,
            gp = grid::gpar(fontsize = 10, fontface = "bold", col = "grey20")
        ),
        xmin = -0.3, xmax = -0.3, ymin = 9, ymax = 9
    ) +
    annotation_custom(
        grob = grid::textGrob("TED Biology",
            rot = 90,
            gp = grid::gpar(fontsize = 10, fontface = "bold", col = "grey20")
        ),
        xmin = -0.3, xmax = -0.3, ymin = 3.5, ymax = 3.5
    ) +
    coord_cartesian(clip = "off")

# =============================================================================
# Save
# =============================================================================
ggsave("TrackA_MR/figures/FigureS1_stepwise_screening.pdf", suppS1,
    width = 12, height = 11, device = cairo_pdf, bg = "white"
)
ggsave("TrackA_MR/figures/FigureS1_stepwise_screening.png", suppS1,
    width = 12, height = 11, dpi = 600, bg = "white"
)
ggsave("TrackA_MR/figures/FigureS1_stepwise_screening.tiff", suppS1,
    width = 12, height = 11, dpi = 600, compression = "lzw", bg = "white"
)

cat("\n✅ Supp S1 FINAL integrated patch saved\n")
cat("   Changes applied:\n")
cat("   1. Conv. → Triangulation\n")
cat("   2. Column headers now include metric name (P / PP.H4 / P / p_adj)\n")
cat("   3. TSHR MVMR: 'Confirmed' → actual P = ", sprintf("%.3f", mvmr_tshr_p), "\n")
cat("   4. TSHR Triangulation: 4/4 → ", mat$triang_score[mat$gene == "TSHR"], "/4 (honest)\n")
cat("   5. Group bracket position fixed (margin expanded)\n")
