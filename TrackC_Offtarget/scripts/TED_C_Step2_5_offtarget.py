import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os

# Create curated sets
hearing = ['NTRN','GAP43','SOX2','ATOH1']
insulin = ['INSR','IRS1','IRS2','IGFBP1','IGFBP2','IGFBP3']
cns = ['BCL2','BAX','CASP3','CASP9','BDNF','ARC','CREB1','MBP','PLP1','OLIG2']

os.makedirs('c:/ProjectTEDGWAS/TrackC_Offtarget/data', exist_ok=True)
pd.DataFrame({'Gene': hearing, 'Tissue': 'Cochlea'}).to_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/hearing_loss_genes.csv', index=False)
pd.DataFrame({'Gene': insulin, 'Tissue': 'Pancreas'}).to_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/insulin_pathway_genes.csv', index=False)
pd.DataFrame({'Gene': cns, 'Tissue': 'CNS'}).to_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/cns_pathway_genes.csv', index=False)

# Network Figure (Fig 4)
G = nx.Graph()
G.add_node('IGF-1R')
G.add_edge('IGF-1R', 'Orbit (Target)')
G.add_edge('IGF-1R', 'Cochlea (Off-target)')
G.add_edge('IGF-1R', 'Pancreas (Off-target)')
G.add_edge('IGF-1R', 'CNS (Off-target)')

G.add_node('TSHR')
G.add_edge('TSHR', 'Orbit (Target)')

for g in hearing: G.add_edge('Cochlea (Off-target)', g)
for g in insulin: G.add_edge('Pancreas (Off-target)', g)
for g in cns: G.add_edge('CNS (Off-target)', g)

plt.figure(figsize=(12, 8))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color='lightgreen', font_size=10, font_weight='bold')
plt.title("TED-TRAP Integrated Off-target Network (Fig 4)")
os.makedirs('c:/ProjectTEDGWAS/TrackC_Offtarget/figures', exist_ok=True)
plt.savefig('c:/ProjectTEDGWAS/TrackC_Offtarget/figures/Offtarget_Network.png', dpi=300)
print("Track C complete.")
