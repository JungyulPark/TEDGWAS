our_genes = ['IGF1R','IGF1','IL1B','CXCL8',
             'TSHR','ADIPOQ','FABP4',
             'HAS2','HAS1','HAS3','TNF','ARRB1',
             'IL6','PPARG','CEBPA']

with open('c:/ProjectTEDGWAS/RNAseq_fulltext.txt', encoding='utf-8') as f:
    text = f.read()

lines = text.split('\n')

print("=== Full manuscript length:", len(lines), "lines ===\n")
print("=== Intersection gene mentions ===")
for gene in our_genes:
    matches = [l.strip() for l in lines if gene in l and l.strip()]
    print(f"\n[{gene}] ({len(matches)} mentions):")
    for m in matches[:4]:
        print(f"  > {m[:200]}")

# Also look for log2FC, DEG, differentially expressed keyword context
import re
print("\n=== DEG/FC numerical entries ===")
fc_lines = [l for l in lines if re.search(r'log2|fold.?change|FC|upregulat|downregulat', l, re.I)]
for l in fc_lines[:30]:
    print(l.strip()[:200])
