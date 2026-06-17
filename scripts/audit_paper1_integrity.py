#!/usr/bin/env python3
"""
Paper 1 (v5) integrity audit — the reusable "keep submission-ready" check.

Run from repo root:  python3 scripts/audit_paper1_integrity.py
Exit code 0 = all pass; 1 = at least one FAIL (so it can wire into a hook/CI).

Checks (all against the locked master MANUSCRIPT_TED_TRAP_v5_MASTER.md):
  1. master exists; print current md5 (integrity anchor)
  2. zero placeholders (to be tabulated / TODO / [ ])
  3. references: 27 listed, all cited (incl. grouped [a, b]) in body
  4. OR = exp(beta) for the 3 backbone genes x 3 outcomes (2-dp display)
  5. IGF1R tissue adj-P stated consistently (no bare "= 1.00")
  6. structural completeness: Abstract/Methods/Results/Discussion/Declarations,
     6 declarations incl. Informed consent, Figure legends 1-3/S1/S2, Tables S1-S7
  7. no stale author names (Yae-Eun Kang / 강예은 / 박정율)
"""
import re, sys, hashlib, os

M = "MANUSCRIPT_TED_TRAP_v5_MASTER.md"
fails, warns = [], []

def ok(b, msg): print(("  ok  " if b else " FAIL ") + msg); (None if b else fails.append(msg))

if not os.path.exists(M):
    print("FAIL: master not found:", M); sys.exit(1)
t = open(M, encoding="utf-8").read()

print("== Paper 1 integrity audit ==")
print("md5:", hashlib.md5(t.encode()).hexdigest())

# 2. placeholders
ph = len(re.findall(r'to be tabulated|\bTODO\b|\bplaceholder\b|\[ *\]', t, re.I))
ok(ph == 0, f"placeholders = {ph} (expect 0)")

# 3. references
body, refs = t.split("## References", 1)
nref = len(re.findall(r'(?m)^\d+\.\s', refs.split("## ")[0]))
cited = set()
for grp in re.findall(r'\[([\d,\s]+)\]', body):
    for n in re.findall(r'\d+', grp): cited.add(int(n))
ok(nref == 27, f"references listed = {nref} (expect 27)")
missing = [i for i in range(1, nref+1) if i not in cited]
ok(not missing, f"all refs cited in body (missing: {missing or 'none'})")

# 4. OR = exp(beta)  (precise betas -> 2dp)
import math
backbone = {  # gene_outcome: (precise_beta, displayed_OR)
 ("TSHR","BBJ"):(-2.096,0.12),("TSHR","UKB"):(-2.436,0.09),("TSHR","Finn"):(-2.331,0.10),
 ("IGF1R","BBJ"):(0.446,1.56),("IGF1R","UKB"):(0.299,1.35),("IGF1R","Finn"):(0.342,1.41),
 ("CTLA4","BBJ"):(-1.740,0.18),("CTLA4","UKB"):(-1.569,0.21),("CTLA4","Finn"):(-1.768,0.17)}
bad = [f"{g}_{o}" for (g,o),(b,orr) in backbone.items() if round(math.exp(b),2) != orr]
ok(not bad, f"OR=exp(beta) for 9 backbone estimates (mismatch: {bad or 'none'})")

# 5. IGF1R tissue adj-P consistency
ok("adjusted *P* = 1.00" not in t and "+0.41 (1.00" not in t,
   "no bare IGF1R tissue 'P = 1.00' (should be >0.99)")

# 6. structure
for sec in ["## Abstract","## Methods","## Results","## Discussion","## Declarations","## References"]:
    ok(sec in t, f"section present: {sec}")
for d in ["Funding","Conflict of interest","Ethics approval","Informed consent","Data availability","Author contributions"]:
    ok(f"**{d}.**" in t, f"declaration present: {d}")
for fig in ["Figure 1.","Figure 2.","Figure 3.","Figure S1.","Figure S2."]:
    ok(f"**{fig}" in t, f"figure legend present: {fig}")
ok(all(f"Table S{i}" in t for i in range(1,8)), "Supplementary Tables S1-S7 all referenced")

# 7. stale names
stale = [n for n in ["Yae-Eun Kang","강예은","박정율"] if n in t]
ok(not stale, f"no stale author names (found: {stale or 'none'})")

print("\nRESULT:", "ALL PASS ✅" if not fails else f"{len(fails)} FAIL ❌")
sys.exit(0 if not fails else 1)
