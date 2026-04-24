import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

igf_df = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/IGF1R_only_genes.csv')
tshr_df = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/TSHR_only_genes.csv')
shared_df = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/Shared_genes.csv')

G = nx.Graph()

for g in igf_df['Gene'].dropna(): G.add_node(g, subset='IGF1R Only', color='blue')
for g in tshr_df['Gene'].dropna(): G.add_node(g, subset='TSHR Only', color='red')
for g in shared_df['Gene'].dropna(): G.add_node(g, subset='Shared', color='purple')

for center, df, c in [('Teprotumumab', igf_df, 'blue'), ('Teprotumumab', shared_df, 'blue'),
                      ('TSHR-Target', tshr_df, 'red'), ('TSHR-Target', shared_df, 'red')]:
    G.add_node(center, subset='Drug', color='yellow')
    for g in df['Gene'].dropna(): G.add_edge(center, g, color=c)

plt.figure(figsize=(14, 10))
pos = nx.spring_layout(G, seed=42)
colors = [n[1]['color'] for n in G.nodes(data=True)]
nx.draw(G, pos, with_labels=True, node_color=colors, edge_color='gray', node_size=800, font_size=8, font_weight='bold')
plt.title("TED-TRAP Differential Pathway Architecture")

import os
os.makedirs('c:/ProjectTEDGWAS/TrackB_Network/figures', exist_ok=True)
plt.savefig('c:/ProjectTEDGWAS/TrackB_Network/figures/Differential_Pathway_Network.png', dpi=300)
print("Network diagram generated.")
