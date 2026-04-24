import requests
import pandas as pd
import os

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
QUERY_URL = 'https://maayanlab.cloud/Enrichr/enrich'

gene_sets = ['LINCS_L1000_Chem_Pert_up', 'LINCS_L1000_Chem_Pert_down']

os.makedirs('c:/ProjectTEDGWAS/TrackB_Network/results/cmap', exist_ok=True)

print("Starting CMap Drug Repurposing via LINCS L1000...")
for name in ['IGF1R_only_genes', 'TSHR_only_genes', 'Shared_genes']:
    file_path = f'c:/ProjectTEDGWAS/TrackB_Network/results/{name}.csv'
    if not os.path.exists(file_path): continue

    df = pd.read_csv(file_path)
    if 'Gene' not in df.columns: continue
    
    genes_list = df['Gene'].dropna().tolist()
    if not genes_list: continue

    genes_str = '\n'.join(genes_list)
    print(f"Adding list for {name} ({len(genes_list)} genes)...")
    
    response = requests.post(ENRICHR_URL, files={'list': (None, genes_str), 'description': (None, name)})
    if not response.ok: continue
        
    user_list_id = response.json().get('userListId')
    
    for gs in gene_sets:
        res = requests.get(QUERY_URL, params={'userListId': user_list_id, 'backgroundType': gs})
        if not res.ok: continue
            
        data = res.json()
        if gs not in data: continue
            
        out_df = pd.DataFrame(data[gs], columns=['Rank', 'Term', 'P-value', 'Z-score', 'Combined Score', 'Overlapping Genes', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value'])
        out_file = f'c:/ProjectTEDGWAS/TrackB_Network/results/cmap/{name}_{gs}.csv'
        out_df.to_csv(out_file, index=False)
        print(f"Saved {out_file} ({len(out_df)} records)")

print("CMap Drug Repurposing complete.")
