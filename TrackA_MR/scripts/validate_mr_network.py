import pandas as pd
import codecs

print("\n=== Check 1: MR Results Original Data ===")
text = codecs.open('c:/ProjectTEDGWAS/TrackA_MR/results/clean_mr.txt', encoding='utf-8', errors='ignore').read()
for line in text.split('\n'):
    if ('TSHR' in line or 'IGF1R' in line) and ('IVW:' in line or '===' in line):
        print(line.strip())

df10 = pd.read_csv('c:/ProjectTEDGWAS/TrackA_MR/results/MR_rest10_summary.csv')
tnf = df10[(df10['Gene']=='TNF') & (df10['GWAS']=='Replication(Hyper)')]
print("\n[TNF Replication]")
print(tnf.to_string(index=False))

print("\n=== Check 3: Network Gene Lists ===")
def get_genes(path):
    import os
    if not os.path.exists(path): return []
    df = pd.read_csv(path)
    if 'Gene' in df.columns: return df['Gene'].dropna().tolist()
    return []

igf = get_genes('c:/ProjectTEDGWAS/TrackB_Network/results/IGF1R_only_genes.csv')
tshr = get_genes('c:/ProjectTEDGWAS/TrackB_Network/results/TSHR_only_genes.csv')
share = get_genes('c:/ProjectTEDGWAS/TrackB_Network/results/Shared_genes.csv')
print(f"IGF1R-only (N={len(igf)}): {igf}")
print(f"TSHR-only (N={len(tshr)}): {tshr}")
print(f"Shared (N={len(share)}): {share}")

print("\n=== Check 4: Off-target intersection ===")
cochlea = get_genes('c:/ProjectTEDGWAS/TrackC_Offtarget/results/Cochlear_Offtarget_Genes.csv')
insulin = get_genes('c:/ProjectTEDGWAS/TrackC_Offtarget/results/Insulin_Offtarget_Genes.csv')
print(f"Cochlear off-target (N={len(cochlea)}): {cochlea}")
print(f"Insulin off-target (N={len(insulin)}): {insulin}")
