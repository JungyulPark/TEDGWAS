# =============================================================================
# TED-TRAP — Figure 3 FINAL (3-panel Main + Supp S1)
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

base_family <- "sans" # Arial if available

# =============================================================================
# STEP 1: Load data
# =============================================================================
tpm <- read.csv("TrackA_MR/data_for_fig3/candidate_sample_TPM.csv",
    stringsAsFactors = FALSE, check.names = FALSE
)

# Table 3 data (IGF1R cross-phenotype)
# ⚠️ Gemini: 이 파일이 Main Table 3 CSV임. 경로 확인 필요
tbl3 <- read.csv("TrackA_MR/manuscript/Table3_crossphenotype.csv",
    stringsAsFactors = FALSE
)

cat("=== Table 3 loaded ===\n")
print(tbl3)

# Also need IVW_SE for error bars — pull from original MR results
# ⚠️ 만약 Table 3 CSV에 SE가 없다면 원본 05_cross_phenotype_mr.csv에서 가져옴
if (!"IVW_SE" %in% names(tbl3) && !"IVW_beta" %in% names(tbl3)) {
    raw <- read.csv("TrackA_MR/results/05_cross_phenotype_mr.csv", stringsAsFactors = FALSE)
    cat("Raw MR columns:", paste(names(raw), collapse = ", "), "\n")
    ivw <- raw[
        raw$method == "Inverse variance weighted",
        c("outcome_name", "b", "se", "pval", "nsnp")
    ]
    names(ivw) <- c("Outcome", "beta", "se", "pval", "nsnp")
    # Merge / rebuild tbl3 if needed
    tbl3_plot <- ivw
} else {
    # Parse beta/SE from Table 3 formatted strings if needed
    tbl3_plot <- tbl3
}

# =============================================================================
# STEP 2: Panel A — TSHR per-sample expression (same as before)
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
    geom_hline(
        yintercept = ted_mean, linetype = "22",
        color = "#8B0000", linewidth = 0.55
    ) +
    annotate("text",
        x = 5.3, y = ted_mean, label = "TED mean",
        size = 2.6, color = "#8B0000", hjust = 0, fontface = "italic"
    ) +
    geom_text(aes(label = sprintf("%.2f", tpm)), vjust = -0.55, size = 2.9) +
    scale_fill_manual(
        values = c("Control" = "grey60", "TED" = "#B22222"),
        guide = "none"
    ) +
    annotate("text",
        x = 3, y = max(tshr_sample$tpm) * 1.3,
        label = "bold('6.5× up in TED,')~italic(log[2]*FC)*'='*'+2.33, '*italic(p[adj])*'='*'0.006'",
        size = 2.9, parse = TRUE
    ) +
    coord_cartesian(ylim = c(0, max(tshr_sample$tpm) * 1.45), clip = "off") +
    labs(
        title = "A. TSHR per-sample expression in orbital tissue",
        x = NULL, y = "TSHR expression (TPM)"
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )

# =============================================================================
# STEP 3: Panel B — Insulin cassette log2FC (same as before)
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
    coord_cartesian(ylim = c(0, 1.35), clip = "off") +
    annotate("text",
        x = 3, y = 1.28,
        label = "italic('5/5 directional up, Binomial P = 0.031')",
        size = 2.9, parse = TRUE
    ) +
    labs(
        title = "B. Insulin signaling cassette in orbital tissue",
        x = NULL, y = expression(log[2] * " fold-change (TED / Ctrl)")
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold")
    )

# =============================================================================
# STEP 4: Panel C — NEW — IGF1R cross-phenotype MR forest plot
# =============================================================================
# Use raw MR results for beta/SE/CI
raw <- read.csv("TrackA_MR/results/05_cross_phenotype_mr.csv", stringsAsFactors = FALSE)
ivw <- raw %>%
    filter(method == "Inverse variance weighted") %>%
    select(outcome = outcome_name, role = outcome_role, nsnp, beta = b, se, pval) %>%
    mutate(
        lci = beta - 1.96 * se,
        uci = beta + 1.96 * se
    )

cat("\n=== IVW for forest plot ===\n")
print(ivw)

# Re-label per Main Table 3 categories (exclude Skin color; it's in Supp)
ivw_plot <- ivw %>%
    filter(!grepl("skin|hearing", outcome, ignore.case = TRUE)) %>% # Skin→Supp, Hearing N_IV=0
    mutate(
        Outcome = case_when(
            grepl("Fasting.?Glucose|glucose", outcome, ignore.case = TRUE) ~ "Fasting glucose",
            grepl("HbA1c", outcome, ignore.case = TRUE) ~ "HbA1c",
            grepl("Type.?2.?Diabetes|T2D|diabetes", outcome, ignore.case = TRUE) ~ "Type 2 diabetes",
            grepl("Height|height", outcome, ignore.case = TRUE) ~ "Height (positive control)",
            TRUE ~ outcome
        ),
        Category = case_when(
            grepl("Height", Outcome) ~ "Positive control",
            TRUE ~ "Metabolic off-target"
        ),
        is_sig = pval < 0.05,
        label_pval = ifelse(pval < 0.001, formatC(pval, format = "e", digits = 1),
            sprintf("%.3f", pval)
        )
    )

# Add Weighted Median P for Height (anchor signal)
wm <- raw %>%
    filter(
        method == "Weighted median",
        grepl("Height", outcome_name, ignore.case = TRUE)
    ) %>%
    pull(pval)
wm_height <- if (length(wm) > 0) formatC(wm, format = "e", digits = 2) else "NA"

# Order: metabolic first, then positive control (bottom)
ivw_plot$Outcome <- factor(ivw_plot$Outcome,
    levels = rev(c(
        "Fasting glucose", "HbA1c", "Type 2 diabetes",
        "Height (positive control)"
    ))
)

panelC <- ggplot(ivw_plot, aes(x = beta, y = Outcome, color = Category)) +
    geom_vline(xintercept = 0, linetype = "22", color = "grey40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.15, linewidth = 0.7) +
    geom_point(size = 3.2, shape = 18) +
    # Annotate P-values to right of plot
    geom_text(
        aes(
            x = max(ivw_plot$uci) * 1.15,
            label = paste0("P = ", label_pval, " (n_IV=", nsnp, ")")
        ),
        color = "black", hjust = 0, size = 2.9, family = base_family
    ) +
    # Annotate WM P for Height (anchor signal)
    annotate("text",
        x = max(ivw_plot$uci) * 1.15, y = 1,
        label = paste0("WM P = ", wm_height, " \u2605"),
        color = "#BF9000", hjust = 0, size = 2.9, fontface = "bold",
        vjust = 2.3
    ) +
    scale_color_manual(
        values = c(
            "Metabolic off-target" = "#4472C4",
            "Positive control" = "#BF9000"
        ),
        name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.55))) +
    labs(
        title = "C. IGF1R cross-phenotype Mendelian randomization",
        subtitle = "Germline IGF1R perturbation does not cause metabolic phenotypes (null IVW); Height serves as biological anchor (WM significant).",
        x = expression("IGF1R cis-eQTL causal effect on outcome (" * beta * ", 95% CI)"),
        y = NULL
    ) +
    theme_bw(base_size = 10, base_family = base_family) +
    theme(
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(
            size = 9, face = "italic", color = "grey25",
            margin = margin(b = 8)
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.4, "cm"),
        axis.text.y = element_text(face = "bold")
    )

# =============================================================================
# STEP 5: Assemble Main Figure 3
# =============================================================================
combined <- (panelA | panelB) / panelC +
    plot_layout(heights = c(1, 1.1))

ggsave("TrackA_MR/figures/Figure3_main_3panel.pdf", combined,
    width = 11, height = 9, device = cairo_pdf, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_main_3panel.png", combined,
    width = 11, height = 9, dpi = 600, bg = "white"
)
ggsave("TrackA_MR/figures/Figure3_main_3panel.tiff", combined,
    width = 11, height = 9, dpi = 600, compression = "lzw", bg = "white"
)

cat("\n✅ Main Figure 3 (3-panel) saved\n")

# =============================================================================
# STEP 6: Rename existing Panel A heatmap → Supp Figure S1
# (Remove 'Not tested' labels per user feedback)
# =============================================================================
# 기존 Figure3_panelABC_v3 파일은 그대로 두고, heatmap only version을 Supp S1으로 새로 만들지
# 아니면 기존 파일을 그냥 rename만 할지 결정 필요.
# 간단 처리: 기존 v3 heatmap 부분만 추출해서 Supp S1로 별도 저장

# → 다음 patch에서 Supp S1 clean version 별도로 처리 (not-tested 텍스트 제거한 버전)
cat("\n⚠️  Note: Supp Figure S1 (20-gene heatmap without 'Not tested' text)\n")
cat("   → 별도 patch로 처리 예정\n")
