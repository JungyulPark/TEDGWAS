import pandas as pd

genes_to_check = ['GNAS', 'ADCY3', 'ADCY5', 'ADCY6', 'ADCY7', 'ADCY9', 'PRKACA', 'PRKACB', 'CREB1', 'CREB3', 'ATF1']
zones = {
    'IGF1R_only': pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/IGF1R_only_genes.csv')['Gene'].tolist(),
    'TSHR_only': pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/TSHR_only_genes.csv')['Gene'].tolist(),
    'Shared': pd.read_csv('c:/ProjectTEDGWAS/TrackB_Network/results/Shared_genes.csv')['Gene'].tolist()
}

print("=== Gene Zone Check ===")
for g in genes_to_check:
    found = "Not in TED intersection (Means not in TED_disease_genes.csv)"
    for zname, zgenes in zones.items():
        if g in zgenes:
            found = zname
    print(f"{g}: {found}")
