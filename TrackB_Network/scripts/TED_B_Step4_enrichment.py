import requests
import pandas as pd
import os

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
QUERY_URL = 'https://maayanlab.cloud/Enrichr/enrich'

gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021']

# Make output directory
os.makedirs('c:/ProjectTEDGWAS/TrackB_Network/results/enrichment', exist_ok=True)

print("Starting Enrichr analysis...")
for name in ['IGF1R_only_genes', 'TSHR_only_genes', 'Shared_genes']:
    file_path = f'c:/ProjectTEDGWAS/TrackB_Network/results/{name}.csv'
    if not os.path.exists(file_path):
        print(f"Skipping {name}: file not found")
        continue

    df = pd.read_csv(file_path)
    if 'Gene' not in df.columns: continue
    
    genes_list = df['Gene'].dropna().tolist()
    if len(genes_list) == 0:
        print(f"Skipping {name}: no genes")
        continue

    genes_str = '\n'.join(genes_list)
    print(f"Adding list for {name} ({len(genes_list)} genes)...")
    
    response = requests.post(ENRICHR_URL, files={'list': (None, genes_str), 'description': (None, name)})
    if not response.ok:
        print(f"Failed to add list to Enrichr: {response.text}")
        continue
        
    user_list_id = response.json().get('userListId')
    
    for gs in gene_sets:
        res = requests.get(QUERY_URL, params={'userListId': user_list_id, 'backgroundType': gs})
        if not res.ok:
            print(f"Failed to fetch enrichment for {gs}")
            continue
            
        data = res.json()
        if gs not in data:
            continue
            
        out_df = pd.DataFrame(data[gs], columns=['Rank', 'Term', 'P-value', 'Z-score', 'Combined Score', 'Overlapping Genes', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value'])
        out_file = f'c:/ProjectTEDGWAS/TrackB_Network/results/enrichment/{name}_{gs}.csv'
        out_df.to_csv(out_file, index=False)
        print(f"Saved {out_file} ({len(out_df)} records)")

print("Enrichr analysis complete.")
