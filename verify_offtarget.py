import pandas as pd
import os

BASE = r'c:\ProjectTEDGWAS'

# 1. Intersection gene names check
igf1r = set(pd.read_csv(os.path.join(BASE, 'TrackB_Network/data/drug_targets/IGF1R_pathway_genes.csv'))['Gene'].dropna())
tshr = set(pd.read_csv(os.path.join(BASE, 'TrackB_Network/data/drug_targets/TSHR_pathway_genes.csv'))['Gene'].dropna())
ted = set(pd.read_csv(os.path.join(BASE, 'TrackB_Network/data/disease_genes/TED_disease_genes.csv'))['Gene'].dropna())

igf1r_only = sorted(igf1r & ted - tshr)
tshr_only = sorted(tshr & ted - igf1r)
shared = sorted(igf1r & tshr & ted)

print('=== 1. INTERSECTION GENE NAMES ===')
print(f'IGF1R-only (4 expected): {", ".join(igf1r_only)}')
print(f'TSHR-only (3 expected): {", ".join(tshr_only)}')
print(f'Shared (8 expected): {", ".join(shared)}')
print(f'Total: {len(igf1r_only) + len(tshr_only) + len(shared)}')

# 2. Off-target files check
print('\n=== 2. OFF-TARGET FILE COUNTS ===')
def check_len(path):
    if not os.path.exists(path): return f'FILE NOT FOUND'
    return len(pd.read_csv(path))

print(f'Cochlear intersection (38 expected): {check_len(os.path.join(BASE, "TrackC_Offtarget/results/cochlear_intersection_full.csv"))}')
print(f'Insulin intersection (51 expected): {check_len(os.path.join(BASE, "TrackC_Offtarget/results/insulin_intersection_full.csv"))}')
print(f'IGF1R expanded downstream (359 expected): {check_len(os.path.join(BASE, "TrackC_Offtarget/data/IGF1R_expanded_downstream_genes.csv"))}')
print(f'Hearing loss (855 expected): {check_len(os.path.join(BASE, "TrackC_Offtarget/data/DisGeNET_Hearing_Loss_genes.csv"))}')
print(f'Insulin signaling (137 expected): {check_len(os.path.join(BASE, "TrackC_Offtarget/data/KEGG_Insulin_Signaling_genes.csv"))}')
