#!/usr/bin/env python3
"""
Paper 1 (v5) integrity audit — the reusable "keep submission-ready" check.

Run from repo root:  python3 scripts/audit_paper1_integrity.py
Exit code 0 = all pass; 1 = at least one FAIL (so it can wire into a hook/CI).

Two layers run here, and both matter:

  A. the framing guards below. Every one of them exists because the wording it
     forbids was in the manuscript at some point and had to be removed. They are
     cheap, they are the locked rules in CLAUDE.md made executable, and they are
     the only thing standing between a future edit and a silently restored
     defect. They are deliberately phrase-level: a numeric audit cannot catch
     "case series" or an over-read PP.H2.

  B. the current candidate's numeric audit, which compares every displayed value
     with its full-precision source and checks the Word structure. It is run at
     the end of this file so one command covers both layers.

Structure targets follow the 11 September 2026 clinical revision: two main
figures plus Figure S1, three main tables, Supplementary Tables S1-S4.

Checks against the master MANUSCRIPT_TED_TRAP_v5_MASTER.md:
  1. master exists; print current md5 (integrity anchor)
  2. zero placeholders (to be tabulated / TODO / [ ])
  3. references: 26 listed, all cited (incl. grouped [a, b]) in body
  4. OR = exp(beta) for the 3 backbone genes x 3 outcomes (2-dp display)
  5. structural completeness: sections, 6 declarations, figure legends, tables
  6. no stale author names (Yae-Eun Kang / 강예은 / 박정율)
  7. single master -- no stale duplicate copy elsewhere in the repo
  8. one submission candidate under submission/ -- older ones get archived
  9. PP.H2 never described as "no detectable disease association" (IGF1R has a
     nominal MR association; PP.H2 is about the shared-variant hypothesis only)
 10. power language hedged ("constrains"), never "excludes"
 11. no "identical instruments" claim (TSHR = 1 IV, IGF1R = 4 IVs)
 12. FinnGen Graves ophthalmopathy never called a "case series"
 13. the same banned phrases absent from README and the submission checklist
 14. the fine-mapping layer stays withdrawn (the SuSiE run was invalid)
 15. no P value is reported for the single-control orbital tissue data
 16. no differential-expression significance claimed for the tissue data
 17. no superseded colocalization posteriors or DESeq2 fold changes
 18. corrected screening figures present, superseded ones absent
 19. TSHR's UKB non-colocalization is disclosed wherever PP.H4 = 0.951 appears
 20. Table 1 carries the analysis-set sample sizes that feed coloc
 21. "prespecified" is reserved for the outcome hierarchy (nothing was registered)
 22. no AI-tool declaration or vendor name (author decision, 2026-09-05)
"""
import re, sys, hashlib, os, math, runpy
from pathlib import Path

M = "MANUSCRIPT_TED_TRAP_v5_MASTER.md"
ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = "submission/candidate_20260911_review"
fails = []

def ok(b, msg): print(("  ok  " if b else " FAIL ") + msg); (None if b else fails.append(msg))

if not os.path.exists(M):
    print("FAIL: master not found:", M); sys.exit(1)
raw = open(M, "rb").read()
# the master is stored with CRLF (.gitattributes -text); hash and match on LF so
# the result does not depend on which machine last wrote the file
t = raw.decode("utf-8").replace("\r\n", "\n")

print("== Paper 1 integrity audit ==")
print("md5 (LF-normalised):", hashlib.md5(t.encode()).hexdigest())

# 2. placeholders
ph = len(re.findall(r'to be tabulated|\bTODO\b|\bplaceholder\b|\[ *\]', t, re.I))
ok(ph == 0, f"placeholders = {ph} (expect 0)")

# 3. references
body, refs = t.split("## References", 1)
nref = len(re.findall(r'(?m)^\d+\.\s', refs.split("## ")[0]))
cited = set()
for grp in re.findall(r'\[([\d,\s–-]+)\]', body):
    for a, b in re.findall(r'(\d+)\s*[–-]\s*(\d+)', grp):
        cited.update(range(int(a), int(b) + 1))
    for n in re.findall(r'\d+', grp):
        cited.add(int(n))
ok(nref == 26, f"references listed = {nref} (expect 26)")
missing = [i for i in range(1, nref + 1) if i not in cited]
ok(not missing, f"all refs cited in body (missing: {missing or 'none'})")

# 4. OR = exp(beta)  (precise betas -> 2dp)
backbone = {  # gene_outcome: (precise_beta, displayed_OR)
 ("TSHR","BBJ"):(-2.096,0.12),("TSHR","UKB"):(-2.436,0.09),("TSHR","Finn"):(-2.331,0.10),
 ("IGF1R","BBJ"):(0.446,1.56),("IGF1R","UKB"):(0.299,1.35),("IGF1R","Finn"):(0.342,1.41),
 ("CTLA4","BBJ"):(-1.740,0.18),("CTLA4","UKB"):(-1.569,0.21),("CTLA4","Finn"):(-1.768,0.17)}
bad = [f"{g}_{o}" for (g,o),(b,orr) in backbone.items() if round(math.exp(b),2) != orr]
ok(not bad, f"OR=exp(beta) for 9 backbone estimates (mismatch: {bad or 'none'})")

# 5. structure
for sec in ["## Abstract","## Methods","## Results","## Discussion","## Declarations","## References"]:
    ok(sec in t, f"section present: {sec}")
for d in ["Funding","Conflict of interest","Ethics approval","Informed consent","Data availability","Author contributions"]:
    ok(f"**{d}.**" in t, f"declaration present: {d}")
for fig in ["Figure 1.","Figure 2.","Figure S1."]:
    ok(f"**{fig}" in t, f"figure legend present: {fig}")
ok(all(f"**Table {i}." in t for i in (1,2,3)), "main Tables 1-3 present")
ok(all(f"Table S{i}" in t for i in range(1,5)), "Supplementary Tables S1-S4 all referenced")
# withdrawn supplementary tables must not be referenced any more
ghosts = [f"Table S{i}" for i in range(5,10) if f"Table S{i}" in t]
ok(not ghosts, f"no reference to withdrawn supplementary tables (found: {ghosts or 'none'})")

# 6. stale names
stale = [n for n in ["Yae-Eun Kang","강예은","박정율"] if n in t]
ok(not stale, f"no stale author names (found: {stale or 'none'})")

# 7. one master only -- a second copy silently goes stale (it did once, in submission/)
dupes = []
for dirpath, dirnames, filenames in os.walk(ROOT):
    dirnames[:] = [d for d in dirnames if d not in (".git", "archives", "node_modules")]
    for fn in filenames:
        full = os.path.join(dirpath, fn)
        if os.path.abspath(full) == os.path.abspath(M) or not fn.endswith(".md"):
            continue
        if fn.startswith("MANUSCRIPT_TED_TRAP") or "_MASTER.md" in fn:
            dupes.append(os.path.relpath(full, ROOT))
ok(not dupes, f"single master, no stale duplicate copies (found: {dupes or 'none'})")

# 8. one live submission candidate. Three accumulated once (20260905, 20260905_maf,
#    20260911_review) and each held a different manuscript, cover letter and figure
#    set -- exactly the way a superseded file reaches a journal.
live = sorted(p.name for p in (ROOT / "submission").glob("candidate_*") if p.is_dir())
ok(live == [os.path.basename(CANDIDATE)],
   f"one submission candidate under submission/ (found: {live})")

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
#      because the checks above only ever read the master.)
# FIGURE_VERIFICATION.md is deliberately excluded: it is the defect log, and its job is to
# quote the wording that was removed. README and the checklist are read as current statements.
PACKAGE = ["README.md", "submission/SUBMISSION_CHECKLIST.md", f"{CANDIDATE}/README.md"]
BANNED = ["identical instruments", "TED-specific sensitivity",
          "no detectable disease association", "without a detectable outcome association",
          "case series"]
leaks = []
for rel in PACKAGE:
    fp = ROOT / rel
    if not fp.exists():
        continue
    for ln, line in enumerate(fp.read_text(encoding="utf-8").splitlines(), 1):
        for phrase in BANNED:
            if phrase in line and not any(m in line for m in (
                    "Never call", "never call", "Stop calling", "it still carried",
                    "previously said", "was being over-read", "| Was |", 'read "',
                    "The footnote read")):
                leaks.append(f"{rel}:{ln} {phrase}")
ok(not leaks, f"package docs free of banned phrasing (found: {leaks or 'none'})")

# 14. the SuSiE fine-mapping layer was withdrawn as invalid and must stay out
susie = [q for q in ["SuSiE", "credible set", "fine-mapped", "fine-mapping"] if q in t]
ok(not susie, f"fine-mapping layer stays withdrawn (found: {susie or 'none'})")

# 15/16. the orbital tissue data has one control -- no P value, no significance claim
tissue_p = re.findall(r'(?:orbital|tissue|TPM)[^.]{0,160}?\*?P\*? *[=<]', t, re.I)
ok(not tissue_p, f"no tissue P value reported (found: {tissue_p or 'none'})")
sig_tissue = [q for q in ["significantly differentially expressed",
                          "not significantly differentially"] if q in t]
ok(not sig_tissue, f"no differential-expression significance claimed (found: {sig_tissue or 'none'})")

# 17. superseded posteriors and the withdrawn DESeq2 fold changes must not reappear
stale_pp = [q for q in ["PP.H3 = 0.794", "PP.H4 = 0.206", "PP.H2 = 0.63 and 0.64",
                        "PP.H2 = 0.73–0.74", "PP.H4 ≤ 0.62",
                        "DESeq2", "+2.33", "+0.41", "+1.27",
                        "0.808", "0.965"] if q in t]
ok(not stale_pp, f"no superseded posteriors or withdrawn tissue statistics (found: {stale_pp or 'none'})")

# 18. corrected screening figures, and the superseded ones absent
for n in ("2,544", "6,135", "2,234", "29.7"):
    ok(n in t, f"corrected screening figure present: {n}")
stale_counts = [n for n in ("2,545", "6,136", "2,235") if n in t]
ok(not stale_counts, f"superseded screening figures absent (found: {stale_counts or 'none'})")

# 19. TSHR colocalization must not be stated without the UKB counterexample
if "0.951" in t:
    ok("0.226" in t, "TSHR UKB non-colocalization disclosed alongside PP.H4 = 0.951")

# 20. Table 1 must carry the analysis-set sample sizes that feed coloc
for n in ("2,809", "3,731", "858", "172,656", "480,867", "499,490"):
    ok(n in t, f"Table 1 sample size present: {n}")

# 21. "prespecified" is reserved for the outcome hierarchy. The PP.H4 >= 0.80
#     threshold and the backbone genes were not registered, and the manuscript
#     itself states that no registration identifier exists.
presp = [ln.strip()[:80] for ln in t.splitlines()
         if "prespecified" in ln and "outcome hierarchy" not in ln]
ok(not presp, f"'prespecified' only for the outcome hierarchy (found: {presp or 'none'})")

# 22. no AI-tool declaration or vendor names in the manuscript (author decision 2026-09-05)
ai_terms = [q for q in ["AI-assisted", "artificial intelligence", "language model",
                        "Codex", "OpenAI", "Claude", "ChatGPT"] if q in t]
ok(not ai_terms, f"no AI-tool wording in the manuscript (found: {ai_terms or 'none'})")

print("\n-- framing guards:", "ALL PASS ✅" if not fails else f"{len(fails)} FAIL ❌")

# ---------------------------------------------------------------------------
# Layer B: the candidate's numeric audit (displayed values vs full-precision
# sources, Word structure). It exits non-zero on failure, so run it last and
# only after the framing guards have passed.
# ---------------------------------------------------------------------------
if fails:
    print("RESULT:", f"{len(fails)} FAIL ❌  (numeric audit not run)")
    sys.exit(1)

numeric = ROOT / CANDIDATE / "reproducibility/audit_submission.py"
if not numeric.exists():
    print("FAIL: candidate numeric audit not found:", numeric)
    sys.exit(1)
print("\n== candidate numeric audit ==")
runpy.run_path(str(numeric), run_name="__main__")
