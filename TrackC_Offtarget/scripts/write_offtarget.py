import pandas as pd

hearing = set(pd.read_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/hearing_loss_genes.csv')['Gene'].dropna())
insulin = set(pd.read_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/insulin_pathway_genes.csv')['Gene'].dropna())
igf1r = set(pd.read_csv('c:/ProjectTEDGWAS/TrackC_Offtarget/data/igf1r_expanded_genes.csv')['Gene'].dropna())

coch_int = list(igf1r & hearing)
ins_int = list(igf1r & insulin)

lines = []
lines.append("=== Rebuilt Sets ===")
lines.append(f"IGF1R downstream (expanded): {len(igf1r)}")
lines.append(f"Hearing loss genes (DisGeNET): {len(hearing)}")
lines.append(f"Insulin pathway genes (KEGG): {len(insulin)}")

lines.append("\n=== High-Confidence Intersections ===")
lines.append(f"Cochlear intersection: {len(coch_int)}")
lines.append(f"Genes: {coch_int[:20]} ...")

lines.append(f"\nInsulin intersection: {len(ins_int)}")
lines.append(f"Genes: {ins_int[:20]} ...")

with open('c:/ProjectTEDGWAS/TrackC_Offtarget/results/offtarget_report.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
