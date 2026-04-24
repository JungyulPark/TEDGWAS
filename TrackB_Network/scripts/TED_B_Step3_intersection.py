import pandas as pd
import os
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Load datasets
igf1r_df = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/data/drug_targets/IGF1R_pathway_genes.csv')
tshr_df = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/data/drug_targets/TSHR_pathway_genes.csv')
ted_df = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/data/disease_genes/TED_disease_genes.csv')

igf1r_genes = set(igf1r_df['Gene'].dropna())
tshr_genes = set(tshr_df['Gene'].dropna())
ted_genes = set(ted_df['Gene'].dropna())

# Venn diagram
plt.figure(figsize=(10, 8))
venn = venn3([igf1r_genes, tshr_genes, ted_genes], ('IGF1R Pathway', 'TSHR Pathway', 'TED Disease Genes'))
plt.title("TED-TRAP Gene Intersection")
plt.savefig('c:/ProjectTEDGWAS/TrackB_Network/figures/Venn_Intersection.png', dpi=300)
plt.close()

# Compute zones
igf1r_only = (igf1r_genes & ted_genes) - tshr_genes
tshr_only = (tshr_genes & ted_genes) - igf1r_genes
shared = igf1r_genes & tshr_genes & ted_genes

# Save to CSV
os.makedirs('c:/ProjectTEDGWAS/TrackB_Network/results', exist_ok=True)
pd.Series(list(igf1r_only), name="Gene").to_csv('c:/ProjectTEDGWAS/TrackB_Network/results/IGF1R_only_genes.csv', index=False)
pd.Series(list(tshr_only), name="Gene").to_csv('c:/ProjectTEDGWAS/TrackB_Network/results/TSHR_only_genes.csv', index=False)
pd.Series(list(shared), name="Gene").to_csv('c:/ProjectTEDGWAS/TrackB_Network/results/Shared_genes.csv', index=False)

print(f"IGF1R-only TED overlap: {len(igf1r_only)} genes")
print(f"TSHR-only TED overlap: {len(tshr_only)} genes")
print(f"Shared TED overlap: {len(shared)} genes")
