import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import os

print("=== Rebuilding Fig 4 ===")
hearing = pd.read_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/hearing_loss_genes.csv')['Gene'].dropna().tolist()
insulin = pd.read_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/insulin_pathway_genes.csv')['Gene'].dropna().tolist()
igf1r_exp = pd.read_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/igf1r_expanded_genes.csv')['Gene'].dropna().tolist()

coch_int = list(set(igf1r_exp) & set(hearing))
ins_int = list(set(igf1r_exp) & set(insulin))

G = nx.Graph()
G.add_node('IGF-1R\n(Target)', color='red', size=3000)
G.add_node('TSHR\n(Target)', color='blue', size=3000)

h_label = f'Hearing Loss\n({len(coch_int)} Shared Genes)'
i_label = f'Hyperglycemia\n({len(ins_int)} Shared Genes)'

G.add_node(h_label, color='orange', size=2000)
G.add_node(i_label, color='orange', size=2000)

G.add_edge('IGF-1R\n(Target)', h_label)
G.add_edge('IGF-1R\n(Target)', i_label)

plt.figure(figsize=(10,6))
pos = nx.spring_layout(G, seed=42)
colors = [nx.get_node_attributes(G, 'color')[node] for node in G.nodes()]
sizes = [nx.get_node_attributes(G, 'size')[node] for node in G.nodes()]
nx.draw(G, pos, with_labels=True, node_color=colors, node_size=sizes, font_size=10, font_weight='bold', width=2)
plt.title(f"Data-Driven Off-Target Network (Shared PI3K/AKT effectors)")
os.makedirs('c:/ProjectTEDGWAS/TrackC_Offtarget/figures', exist_ok=True)
plt.savefig('c:/ProjectTEDGWAS/TrackC_Offtarget/figures/Offtarget_Network_Rebuilt.png', dpi=300)
plt.close()

print("\n=== Enrichment ===")
def run_enrichr(genes, desc):
    if len(genes) < 3: 
        print(f"\n[{desc}] Too few genes")
        return
    url = 'https://maayanlab.cloud/Enrichr/addList'
    r = requests.post(url, files={'list': (None, '\n'.join(genes)), 'description': (None, desc)})
    if not r.ok: return
    user_list_id = r.json()['userListId']
    
    url_kegg = f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType=KEGG_2021_Human'
    kegg = requests.get(url_kegg).json()
    print(f"\n[{desc}] Top 3 KEGG:")
    if 'KEGG_2021_Human' in kegg:
        for t in kegg['KEGG_2021_Human'][:3]:
            print(f"  {t[1]}: P={t[2]:.2e}")
            
    url_lincs = f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType=LINCS_L1000_Chem_Pert_down'
    lincs = requests.get(url_lincs).json()
    print(f"[{desc}] Top 3 LINCS L1000 (Down):")
    if 'LINCS_L1000_Chem_Pert_down' in lincs:
        for t in lincs['LINCS_L1000_Chem_Pert_down'][:3]:
            print(f"  {t[1]}: P={t[2]:.2e}")

igf = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/IGF1R_only_genes.csv')['Gene'].dropna().tolist()
tshr = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/TSHR_only_genes.csv')['Gene'].dropna().tolist()
shared = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/Shared_genes.csv')['Gene'].dropna().tolist()

run_enrichr(igf, 'IGF1R-only')
run_enrichr(tshr, 'TSHR-only')
run_enrichr(shared, 'Shared')
print("Done.")
