import networkx as nx, pandas as pd, requests

def fetch_network(genes):
    url = "https://version-12-0.string-db.org/api/json/network"
    params = {"identifiers": "%0d".join(genes), "species": 9606, "caller_identity": "TED", "required_score": 400}
    r = requests.post(url, data=params)
    if not r.ok: return []
    return r.json()

igf1r_net = ['IGF1R','IGF1','IGF2','IRS1','IRS2','SHC1','GRB2','SOS1','PIK3CA','PIK3CB','PIK3R1','AKT1','AKT2','MTOR','RPS6KB1','KRAS','BRAF','MAP2K1','MAP2K2','MAPK1','MAPK3','JAK1','JAK2','STAT3','STAT5A']
tshr_net = ['TSHR','FOXO1','FOXO3'] # Pure TSHR source without the extension modules

ext_camp = ['GNAS', 'ADCY3', 'ADCY5', 'ADCY6', 'ADCY7', 'ADCY9', 'PRKACA', 'PRKACB', 'CREB1', 'CREB3', 'ATF1']
ext_arrb = ['ARRB1', 'ARRB2', 'GRK2', 'GRK5']
ext_gq = ['PLCB1', 'PLCB3', 'GNAQ', 'GNA11', 'PRKCB']

all_genes = list(set(igf1r_net + tshr_net + ext_camp + ext_arrb + ext_gq))
edges = fetch_network(all_genes)

G = nx.Graph()
for g in all_genes: G.add_node(g)
for e in edges: G.add_edge(e['preferredName_A'], e['preferredName_B'])

res = []
for name, mod in [('cAMP/PKA', ext_camp), ('beta-arrestin', ext_arrb), ('Gq/PLC', ext_gq)]:
    igf_edges = sum(1 for s in igf1r_net for t in mod if G.has_edge(s, t))
    tshr_edges = sum(1 for s in tshr_net for t in mod if G.has_edge(s, t))
    
    def get_avg_sp(src_net, tgt_net):
        sp = []
        for s in src_net:
            for t in tgt_net:
                if G.has_node(s) and G.has_node(t) and nx.has_path(G, s, t):
                    sp.append(nx.shortest_path_length(G, s, t))
        return round(sum(sp)/len(sp), 2) if sp else "No Path"
        
    res.append({
        'Extension': name, 
        'Direct Edges from IGF1R Net': igf_edges,
        'Direct Edges from TSHR Net': tshr_edges,
        'Avg SP from IGF1R': get_avg_sp(igf1r_net, mod),
        'Avg SP from TSHR': get_avg_sp(tshr_net, mod)
    })

df = pd.DataFrame(res)
print(df.to_string(index=False))
df.to_csv('c:/ProjectTEDGWAS/TrackB_Network/results/Extension_Layer_Connectivity_v2.csv', index=False)
