import pandas as pd
import os

BASE = r'c:\ProjectTEDGWAS'

def check_len(path):
    if not os.path.exists(path): return 'Not Found'
    return len(pd.read_csv(path))

print(f'Hearing loss: {check_len(os.path.join(BASE, "TrackC_Offtarget/data/hearing_loss_genes.csv"))}')
print(f'IGF1R expanded downstream: {check_len(os.path.join(BASE, "TrackC_Offtarget/data/igf1r_expanded_genes.csv"))}')
print(f'Insulin signaling: {check_len(os.path.join(BASE, "TrackC_Offtarget/data/insulin_pathway_genes.csv"))}')
