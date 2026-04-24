import csv

# Read CTD genes
ctd_genes = {}
with open(r"c:\ProjectTEDGWAS\CTD_D049970_genes_20260421243921.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene = row["Gene Symbol"]
        direct = row.get("Direct Evidence", "").strip()
        score = float(row.get("Inference Score", 0) or 0)
        ctd_genes[gene] = {"direct": direct, "score": score}

# Read our 59 disease genes
our_genes = []
with open(r"c:\ProjectTEDGWAS\TrackB_Network\data\disease_genes\TED_disease_genes.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        our_genes.append(row["Gene"])

print(f"CTD total: {len(ctd_genes)}")
print(f"Our total: {len(our_genes)}")

direct_count = sum(1 for g in ctd_genes.values() if g["direct"])
print(f"Direct evidence: {direct_count}")
print()

# Cross-reference
in_ctd = []
not_in_ctd = []
for g in our_genes:
    if g in ctd_genes:
        info = ctd_genes[g]
        in_ctd.append((g, info["direct"], info["score"]))
    else:
        not_in_ctd.append(g)

print(f"Our genes IN CTD: {len(in_ctd)}/{len(our_genes)}")
for g, direct, score in sorted(in_ctd, key=lambda x: -x[2]):
    d_flag = " ***DIRECT***" if direct else ""
    print(f"  {g}: score={score:.2f}{d_flag}")

print(f"\nOur genes NOT in CTD: {len(not_in_ctd)}")
for g in not_in_ctd:
    print(f"  {g}")

# Which genes were in the original "44" analysis (intersection analysis)
# Original intersection was: IGF1R pathway (25) + TSHR pathway (24) intersected with disease genes
# The 15 intersection genes are known. The remaining ~29 are non-intersecting disease genes
intersection_15 = ["IGF1R", "IGF1", "IL1B", "CXCL8", "TSHR", "ADIPOQ", "FABP4",
                    "HAS2", "HAS1", "TNF", "ARRB1", "IL6", "HAS3", "PPARG", "CEBPA"]

print(f"\n=== INTERSECTION ANALYSIS ===")
print(f"15 intersection genes: {intersection_15}")
print(f"\n59 disease genes total")
print(f"15 intersection genes (in disease + pathway)")
print(f"44 remaining disease genes (not in pathway intersection)")
remaining_44 = [g for g in our_genes if g not in intersection_15]
print(f"Remaining: {len(remaining_44)}")
for g in remaining_44:
    ctd_flag = f" (CTD score={ctd_genes[g]['score']:.1f})" if g in ctd_genes else " (NOT in CTD)"
    print(f"  {g}{ctd_flag}")
