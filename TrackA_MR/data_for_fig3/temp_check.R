df <- read.csv('C:/ProjectTEDGWAS/TrackA_MR/data_for_fig3/candidate_sample_TPM.csv')
genes <- c('ARRB1', 'TNF', 'AKT1', 'PIK3R1', 'HAS1', 'HAS2', 'HAS3', 'PPARG', 'ADIPOQ', 'IGF1R', 'TSHR')
df_sub <- subset(df, Gene_Symbol %in% genes)

results <- data.frame(Gene=character(), Ctrl_TPM=numeric(), TED_Mean_TPM=numeric(), stringsAsFactors=FALSE)

for(i in 1:nrow(df_sub)) {
  g <- df_sub$Gene_Symbol[i]
  ctrl_cols <- c('X2__1_TPM', 'X2__2_TPM')
  ted_cols <- c('X10__1_TPM','X10__2_TPM','X11__1_TPM','X11__2_TPM','X7__1_TPM','X7__2_TPM','X8__1_TPM','X8__2_TPM')
  
  # Note: R reads columns starting with number as X...
  ctrl_mean <- mean(as.numeric(df_sub[i, ctrl_cols]))
  ted_mean <- mean(as.numeric(df_sub[i, ted_cols]))
  
  results <- rbind(results, data.frame(Gene=g, Ctrl_TPM=ctrl_mean, TED_Mean_TPM=ted_mean))
}
print(results)
