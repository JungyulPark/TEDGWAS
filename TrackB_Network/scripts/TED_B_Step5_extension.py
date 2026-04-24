import pandas as pd
import requests
import networkx as nx

def fetch_string_network(genes):
    STRING_API_URL = "https://version-12-0.string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(genes),
        "species": 9606,
        "caller_identity": "TED_TRAP_Project",
        "required_score": 700
    }
    response = requests.post(STRING_API_URL, data=params)
    if not response.ok: return []
    return response.json()

# Load targets
igf1r_genes = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/IGF1R_only_genes.csv')['Gene'].tolist()
tshr_genes = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/TSHR_only_genes.csv')['Gene'].tolist()
shared_genes = pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/Shared_genes.csv')['Gene'].tolist()

ext_arrb = ['ARRB1', 'ARRB2', 'GRK2', 'GRK5']
ext_cd34 = ['CD34', 'THY1', 'CXCR4', 'CXCL12']
all_ext = ext_arrb + ext_cd34

def analyze_connectivity(source_genes, name):
    print(f"\nAnalyzing {name} connectivity to Extension Layer...")
    query_genes = source_genes + all_ext
    edges = fetch_string_network(query_genes)
    
    G = nx.Graph()
    for e in edges:
        G.add_edge(e['preferredName_A'], e['preferredName_B'], weight=e['score'])
        
    for ext_set, set_name in zip([ext_arrb, ext_cd34], ["ARRB1/2 Scaffold", "CD34+ Fibrocyte"]):
        print(f"  [{set_name}]")
        total_edges = 0
        for s in source_genes:
            for t in ext_set:
                if G.has_edge(s, t):
                    total_edges += 1
                    print(f"    Direct edge: {s} - {t} (score: {G[s][t]['weight']})")
        print(f"    Total direct edges: {total_edges}")

analyze_connectivity(igf1r_genes + shared_genes, "IGF1R Network (Teprotumumab)")
analyze_connectivity(tshr_genes + shared_genes, "TSHR Network (Targeting)")
