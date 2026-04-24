"""
Phase A-10: Off-target Intersection Extraction
Extract full gene lists for cochlear (38) and insulin (51) intersections.
"""
import csv
import os

DATA_DIR = r"c:\ProjectTEDGWAS\TrackC_Offtarget\data"
OUT_DIR = r"c:\ProjectTEDGWAS\TrackC_Offtarget\results"

def read_gene_set(filepath):
    genes = set()
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if row and row[0].strip():
                genes.add(row[0].strip())
    return genes

# Read all gene sets
igf1r = read_gene_set(os.path.join(DATA_DIR, 'igf1r_expanded_genes.csv'))
hearing = read_gene_set(os.path.join(DATA_DIR, 'hearing_loss_genes.csv'))
insulin = read_gene_set(os.path.join(DATA_DIR, 'insulin_pathway_genes.csv'))

print(f"IGF1R expanded: {len(igf1r)}")
print(f"Hearing loss: {len(hearing)}")
print(f"Insulin pathway: {len(insulin)}")

# Compute intersections
cochlear = sorted(igf1r & hearing)
insulin_inter = sorted(igf1r & insulin)

print(f"\nCochlear intersection: {len(cochlear)}")
print("Genes:", cochlear)

print(f"\nInsulin intersection: {len(insulin_inter)}")
print("Genes:", insulin_inter)

insr_check = "INSR" in insulin_inter
print(f"\nINSR in insulin intersection: {insr_check}")

# Save full lists
out1 = os.path.join(OUT_DIR, 'cochlear_intersection_full.csv')
with open(out1, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Gene'])
    for g in cochlear:
        writer.writerow([g])
print(f"Saved: {out1}")

out2 = os.path.join(OUT_DIR, 'insulin_intersection_full.csv')
with open(out2, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Gene'])
    for g in insulin_inter:
        writer.writerow([g])
print(f"Saved: {out2}")
