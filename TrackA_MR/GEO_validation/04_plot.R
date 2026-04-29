setwd("c:/ProjectTEDGWAS/TrackA_MR/GEO_validation/")

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(data.table)
library(ggplot2)

comp <- fread("Cross_validation_inhouse_vs_GSE105149.csv")
comp <- comp[!is.na(inhouse_logFC) & !is.na(logFC), ]

png("Figure_S2_cross_validation_scatter.png", width = 6, height = 6, res = 300, units = "in")
p <- ggplot(comp, aes(x = inhouse_logFC, y = logFC)) +
    geom_point(aes(color = direction_match), size = 3) +
    geom_text(aes(label = hgnc_symbol), vjust = -0.8, size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dotted") +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "darkgreen")) +
    labs(
        x = "In-house TED cohort log2FC",
        y = "GSE105149 (Rosenbaum 2018) log2FC",
        title = "Candidate gene cross-validation",
        color = "Direction\nconcordant"
    ) +
    theme_bw()
print(p)
invisible(dev.off())
