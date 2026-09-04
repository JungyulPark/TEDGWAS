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
  8. single master — no stale duplicate copy elsewhere in the repo
  9. PP.H2 never described as "no detectable disease association" (IGF1R has a
     nominal MR association; PP.H2 is about the shared-variant hypothesis only)
 10. power language hedged ("constrains"), never "excludes"
 11. no "identical instruments" claim (TSHR = 1 IV, IGF1R = 4 IVs)
 12. FinnGen Graves ophthalmopathy never called a "case series"
 13. the same banned phrases absent from README and the submission checklist
 14. outcome hierarchy described as functional (case counts 2,809 / 3,731 / 858 are
     not monotonically decreasing, and TED specificity dips at the UKB step)
 15. the fine-mapping layer stays withdrawn (the SuSiE run was invalid)
 16. no P value is reported for the single-control orbital tissue data
 17. the null is phrased as a criterion failure, not "no robust novel target"
 18. corrected screening figures present, superseded ones absent
 19. TSHR's UKB non-colocalization is disclosed wherever PP.H4 = 0.951 appears
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
ok(all(f"Table S{i}" in t for i in range(1,9)), "Supplementary Tables S1-S8 all referenced")

# 7. stale names
stale = [n for n in ["Yae-Eun Kang","강예은","박정율"] if n in t]
ok(not stale, f"no stale author names (found: {stale or 'none'})")

# 8. one master only -- a second copy silently goes stale (it did once, in submission/)
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dupes = []
for dirpath, dirnames, filenames in os.walk(root):
    dirnames[:] = [d for d in dirnames if d not in (".git", "archives", "node_modules")]
    for fn in filenames:
        full = os.path.join(dirpath, fn)
        if os.path.abspath(full) == os.path.abspath(M) or not fn.endswith(".md"):
            continue
        if fn.startswith("MANUSCRIPT_TED_TRAP") or "_MASTER.md" in fn:
            dupes.append(os.path.relpath(full, root))
ok(not dupes, f"single master, no stale duplicate copies (found: {dupes or 'none'})")

# 9. IGF1R: PP.H2 must never be described as absence of a disease association
overread = [q for q in ["without a detectable disease association",
                        "without a detectable outcome association",
                        "no detectable disease association"] if q in t]
ok(not overread, f"PP.H2 not over-read as 'no disease association' (found: {overread or 'none'})")

# 10. the null constrains, it does not exclude (only 35.6% of genes powered for OR>=2.0)
overclaim = [q for q in ["null excludes", "evidence against additional large"] if q in t]
ok(not overclaim, f"power language hedged, not 'excludes' (found: {overclaim or 'none'})")

# 11. TSHR and IGF1R do not share instruments
ok("identical instruments" not in t, "no 'identical instruments' claim")
# 12. FinnGen GO is a case-control GWAS, not a case series
ok("case series" not in t, "FinnGen outcome not called a 'case series'")

# 13. the same banned phrases in the reader-facing package docs, not just the master.
#     ("case series" once survived in SUBMISSION_CHECKLIST.md while this audit said PASS,
#      because checks 1-12 only ever read the master.)
# FIGURE_VERIFICATION.md is deliberately excluded: it is the defect log, and its job is to
# quote the wording that was removed. README and the checklist are read as current statements.
PACKAGE = ["README.md", "submission/SUBMISSION_CHECKLIST.md"]
BANNED = ["identical instruments", "TED-specific sensitivity",
          "no detectable disease association", "without a detectable outcome association"]
leaks = []
for rel in PACKAGE:
    fp = os.path.join(root, rel)
    if not os.path.exists(fp):
        continue
    body = open(fp, encoding="utf-8").read()
    for phrase in BANNED + ["case series"]:
        for ln, line in enumerate(body.splitlines(), 1):
            if phrase not in line:
                continue
            # these files legitimately quote the banned wording when recording its removal
            if any(m in line for m in ("Never call", "never call", "Stop calling",
                                       "it still carried", "previously said", "was being over-read",
                                       "| Was |", "read \"", "The footnote read")):
                continue
            leaks.append(f"{rel}:{ln} {phrase!r}")
ok(not leaks, f"package docs free of banned phrasing (found: {leaks or 'none'})")

# 14. the outcome hierarchy is FUNCTIONAL, not monotone. Case counts run
#     2,809 (BBJ) -> 3,731 (UKB) -> 858 (FinnGen), so sample size does not decrease
#     monotonically, and TED specificity dips at the UKB hyperthyroidism step.
counts = [("2,809", "BBJ Graves"), ("3,731", "UKB hyperthyroidism"), ("858", "FinnGen GO")]
missing_counts = [f"{n} ({who})" for n, who in counts if n not in t]
ok(not missing_counts, f"outcome case counts present (missing: {missing_counts or 'none'})")
monotone = [q for q in ["decreasing sample size", "increasing TED specificity"] if q in t]
ok(not monotone, f"hierarchy described functionally, not as monotone (found: {monotone or 'none'})")

# 15. the fine-mapping layer was withdrawn (invalid SuSiE run — no allele
#     harmonisation between the z-scores and the LD reference). Nothing in the
#     manuscript may reintroduce it.
finemap = [q for q in ["fine-map", "fine map", "SuSiE", "susieR", "credible set",
                       "credible-set", "posterior inclusion probability"] if q in t]
ok(not finemap, f"fine-mapping layer stays withdrawn (found: {finemap or 'none'})")

# 16. the orbital tissue data are descriptive only — no P value may reappear
#     alongside them (one control cannot support differential-expression inference)
tissue_p = [q for q in ["adjusted *P* = 0.032", "adjusted P=0.032", "log2FC=+2.33, adjusted",
                        "padj = 0.032"] if q in t]
ok(not tissue_p, f"no tissue P value reported (found: {tissue_p or 'none'})")

# 17. the null is a criterion failure, not an absence of targets
ok("robust novel target" not in t, "null stated as a criterion, not 'robust novel target'")

# 18. corrected screening figures, and the superseded ones absent
for n in ("2,544", "6,135", "2,234", "29.7"):
    ok(n in t, f"corrected screening figure present: {n}")
stale_counts = [n for n in ("2,545", "6,136", "2,235") if n in t]
ok(not stale_counts, f"superseded screening figures absent (found: {stale_counts or 'none'})")

# 19. TSHR colocalization must not be stated without the UKB counterexample
if "0.951" in t:
    ok("0.774" in t or "rs1023586" in t,
       "TSHR UKB non-colocalization disclosed alongside PP.H4 = 0.951")

print("\nRESULT:", "ALL PASS ✅" if not fails else f"{len(fails)} FAIL ❌")
sys.exit(0 if not fails else 1)
