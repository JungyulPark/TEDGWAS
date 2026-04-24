import networkx as nx
import matplotlib.pyplot as plt

G = nx.DiGraph()

canonical_nodes = ["TSH\n(Ligand)", "TSHR\n(GPCR)", "Gs\n(G-protein)", "AC\n(Adenylyl Cyclase)", "cAMP\n(Messenger)", "PKA\n(Kinase)", "CREB\n(Transcription Factor)"]
for i in range(len(canonical_nodes)-1):
    G.add_edge(canonical_nodes[i], canonical_nodes[i+1], color="blue", style="solid", label="canonical")

igf_nodes = ["IGF-1\n(Ligand)", "IGF-1R\n(RTK)", "PI3K\n(Kinase)", "AKT\n(Kinase)"]
for i in range(len(igf_nodes)-1):
    G.add_edge(igf_nodes[i], igf_nodes[i+1], color="green", style="solid", label="canonical")

G.add_edge("AKT\n(Kinase)", "CREB\n(Transcription Factor)", color="red", style="dashed", label="Not in hsa04024")

pos = {
    "TSH\n(Ligand)": (0, 4), "TSHR\n(GPCR)": (1.5, 4), "Gs\n(G-protein)": (3, 4), "AC\n(Adenylyl Cyclase)": (4.5, 4), 
    "cAMP\n(Messenger)": (6, 4), "PKA\n(Kinase)": (7.5, 4), "CREB\n(Transcription Factor)": (9, 2.5),
    "IGF-1\n(Ligand)": (0, 1), "IGF-1R\n(RTK)": (1.5, 1), "PI3K\n(Kinase)": (3.5, 1), "AKT\n(Kinase)": (5.5, 1)
}

plt.figure(figsize=(14, 6))

edges = G.edges()
colors = [G[u][v]['color'] for u,v in edges]
styles = [G[u][v]['style'] for u,v in edges]

nx.draw(G, pos, with_labels=True, node_color='#F0F0F0', edgecolors='black', node_size=3500, font_size=9, font_weight="bold", arrows=True, arrowsize=20, edge_color=colors, style=styles, width=2.5)

plt.text(1.5, 4.4, "KEGG: hsa04024 (cAMP signaling pathway)", fontsize=11, color="blue", weight="bold")
plt.text(1.5, 0.4, "KEGG: hsa04151 (PI3K-Akt signaling pathway)", fontsize=11, color="green", weight="bold")
plt.text(7, 1.5, "Indirect\n(No directed edge in hsa04024)", fontsize=10, color="red", weight="bold", ha="center")

plt.title("Track B Extension Layer: Biological Directionality (KEGG Overlay)", fontsize=14, pad=20)
plt.savefig("c:/ProjectTEDGWAS/TrackB_Network/results/KEGG_Directed_Extension.png", bbox_inches='tight', dpi=300)
plt.close()
print("KEGG Overlay Plot generated.")
