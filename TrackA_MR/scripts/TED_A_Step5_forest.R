library(ggplot2)
library(dplyr)

cat("Loading MR results...\n")
df_rest <- read.csv("c:/ProjectTEDGWAS/TrackA_MR/results/MR_rest10_summary.csv")
# Back-calculate SE from p-value and b
df_rest$se <- abs(df_rest$b) / qnorm(1 - df_rest$p / 2)

# Add primary genes
df_main <- data.frame(
    Gene = c("IGF1R", "TSHR", "IGF1R", "TSHR"),
    GWAS = c("Primary(Graves)", "Primary(Graves)", "Replication(Hyper)", "Replication(Hyper)"),
    b = c(0.217, -1.394, 0.001, -0.012),
    se = c(0.112, 0.167, 0.001, 0.001),
    p = c(5.27e-02, 6.79e-17, 4.12e-02, 3.88e-29),
    nsnp = c(11, 3, 11, 5)
)

df_all <- bind_rows(df_main, df_rest)

df_all$GWAS <- factor(df_all$GWAS, levels = c("Primary(Graves)", "Replication(Hyper)"))
df_all$LCI <- df_all$b - 1.96 * df_all$se
df_all$UCI <- df_all$b + 1.96 * df_all$se

ordered_genes <- c("TSHR", "IGF1R", "TNF", "CTLA4", "ARRB1", "PPARG", "IRS1", "AKT1")
df_all <- df_all[df_all$Gene %in% ordered_genes, ]
df_all$Gene <- factor(df_all$Gene, levels = rev(ordered_genes))

cat("Plotting...\n")
p <- ggplot(df_all, aes(x = b, y = Gene, color = GWAS)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(xmin = LCI, xmax = UCI), position = position_dodge(width = 0.5), width = 0.2) +
    facet_wrap(~GWAS, scales = "free_x") +
    theme_bw() +
    scale_color_manual(values = c("Primary(Graves)" = "#0072B2", "Replication(Hyper)" = "#D55E00")) +
    labs(
        title = "Mendelian Randomization: Target Genes vs Graves/Hyperthyroidism",
        x = "IVW Causal Estimate (Beta)", y = "Target Gene"
    ) +
    theme(text = element_text(size = 13), legend.position = "none")

ggsave("c:/ProjectTEDGWAS/TrackA_MR/results/MR_Summary_Forest.png", plot = p, width = 10, height = 6)
cat("Saved to MR_Summary_Forest.png\n")
