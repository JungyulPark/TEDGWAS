import requests
import json
import pandas as pd
import os

def fetch_enrichr(lib_name):
    url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib_name}"
    r = requests.get(url)
    data = {}
    if r.ok:
        for line in r.text.strip().split('\n'):
            parts = line.split('\t')
            term = parts[0]
            genes = [g for g in parts[2:] if g]
            data[term] = genes
    return data

print("Downloading KEGG 2021 Human...")
kegg = fetch_enrichr('KEGG_2021_Human')
insulin_genes = []
igf1r_kegg = []
for k, v in kegg.items():
    if 'Insulin signaling pathway' in k:
        insulin_genes.extend(v)
    if 'PI3K-Akt signaling pathway' in k:
        igf1r_kegg.extend(v)

print("Downloading DisGeNET...")
disgenet = fetch_enrichr('DisGeNET')
hearing_genes = set()
for k, v in disgenet.items():
    if 'hearing loss' in k.lower() or 'deafness' in k.lower() or 'sensorineural' in k.lower():
        for g in v: hearing_genes.add(g)
hearing_genes = list(hearing_genes)

igf1r_original = ['IGF1R','IGF1','IGF2','IRS1','IRS2','SHC1','GRB2','SOS1','PIK3CA','PIK3CB','PIK3R1','AKT1','AKT2','MTOR','RPS6KB1','KRAS','BRAF','MAP2K1','MAP2K2','MAPK1','MAPK3','JAK1','JAK2','STAT3','STAT5A']
igf1r_expanded = list(set(igf1r_original + igf1r_kegg))

print("\n=== Rebuilt Sets ===")
print("IGF1R Expanded downstream genes:", len(igf1r_expanded))
print("Hearing loss genes:", len(hearing_genes))
print("Insulin pathway genes:", len(insulin_genes))

coch_int = set(igf1r_expanded) & set(hearing_genes)
ins_int = set(igf1r_expanded) & set(insulin_genes)

print("\n=== Intersections ===")
print("Cochlear intersection:", len(coch_int))
print("Cochlear intersect genes:", list(coch_int)[:30])

print("Insulin intersection:", len(ins_int))
print("Insulin intersect genes:", list(ins_int)[:30])

os.makedirs('c:/ProjectTEDGWAS/TrackC_Offtarget/data', exist_ok=True)
pd.DataFrame({'Gene': list(hearing_genes), 'Tissue': 'Cochlea'}).to_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/hearing_loss_genes.csv', index=False)
pd.DataFrame({'Gene': list(insulin_genes), 'Tissue': 'Pancreas'}).to_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/insulin_pathway_genes.csv', index=False)
pd.DataFrame({'Gene': list(igf1r_expanded), 'Pathway': 'IGF1R_Expanded'}).to_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/igf1r_expanded_genes.csv', index=False)
