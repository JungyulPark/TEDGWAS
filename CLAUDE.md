> **2026-09-05 reviewed revision:** The current master and `submission/candidate_20260905/` supersede the historical numerical claims, fixed scientific assertions, journal status and build instructions below. See the candidate revision report and source verification before editing. In particular, outdated SuSiE/tissue P-value claims must not be restored, the current IGF1R UKB posterior is split H2/H4, and eQTLGen frequency sensitivity remains unperformed **on the GitHub copy** (a run completed on the local machine on 2026-09-05 has not been pushed; see `internal/HANDOFF_local_reanalysis.md`, Task 6). The manuscript carries **no AI-tool declaration** by author decision (audit check 27), and "prespecified" is reserved for the outcome hierarchy (check 26). The only submission package is `submission/candidate_20260905/`; the older root package is archived under `archives/submission_pre20260905/`. The non-negotiable data-management rules below remain in force.

# CLAUDE.md — TED-TRAP project rules

This file orients any AI assistant working in this repository. Read it before
making changes. The guiding value of this project is **정확하고 진실한게 생명**
— accuracy and truth above all. No hallucination; verification first.

## What this project is
Druggable-gene-wide Mendelian randomization (MR) + colocalization +
orbital transcriptomics distinguishing **TSHR-anchored genetic susceptibility**
from **IGF1R pharmacologic effector biology** in Graves disease (GD) and
thyroid eye disease (TED). Current state: **v5**, prepared for submission to the
Journal of Endocrinological Investigation (JEI). Final package in `submission/`.

## Non-negotiable data rules
1. **NEVER commit IRB raw data.** The in-house orbital RNA-seq (`data.txt`,
   FASTQ/BAM/counts) is IRB-protected (Pusan National University Hospital
   2104-018-102) and is "available on request" per the manuscript. It must never
   enter git. `.gitignore` enforces this; do not weaken it. (Note: `data.txt`
   was committed once early and has since been **purged from all history** — see
   `internal/GITHUB_SECURITY_CLEANUP.md`.)
2. **NEVER commit TSHR-ATrap patent material** (sequences, structures,
   `.pdb`/`.fasta`, AlphaFold/HDOCK outputs).
3. **Do not redistribute license-bound data** (eQTLGen full sumstats, FinnGen
   `.gz`, 1000G panels). Commit the *download script*, not the data.
4. Commit only **scripts** and **aggregate results** (summary tables, figure
   code, final figures).

## Locked scientific ground truths (do not silently change)
- **EUR-only LD reference is a locked rule** for clumping/coloc (eQTLGen is
  EUR-based; EAS LD mismatch produced artifactual COJO signals).
- **TSHR** — single-instrument, locus-level anchor (rs179252, chr14:81,435,985
  hg19). SuSiE = 1 credible set; no independent secondary IV. MR β: BBJ −2.10,
  UKB −2.44, FinnGen −2.33. Coloc PP.H4: BBJ 0.951 / FinnGen 0.986.
- **IGF1R** — effector, NOT expression-colocalized. Multi-IV (4). Nominal risk
  direction (BBJ P=0.021, UKB P=0.012, FinnGen NS). Coloc H2-dominant
  (0.043/0.036). Tissue log2FC +0.41 (NS).
- **CTLA4** — positive control. BBJ H3-dominant (0.79), FinnGen H4=0.978.
- **robust_novel = 0** after MHC + chr16p11.2 LD-spillover + cross-outcome coloc
  filtering. This is an *informative* result, not a negative one — do not reframe
  it as a discovery.
- **Tissue n=4 TED + 1 control**, biological-sample level (technical replicates
  collapsed). Correct TSHR padj = 0.032; earlier pseudoreplicated padj values
  (6.6e-5 / 0.006) are WRONG.
- **No "first systematic" claim** — druggable GD MR is already published.
- **FinnGen Graves ophthalmopathy is NOT a TED-specific contrast.** Cases are ascertained among Graves
  disease patients and compared with population controls, so its associations substantially re-measure GD
  susceptibility. TSHR is an anchor **for GD**, "whose signal is also recovered in a TED-enriched
  Graves ophthalmopathy outcome" — never claim a TED-specific effect separable from GD liability. See
  `internal/MVMR_FORENSIC_VERDICT.md`.
- **MVMR results in the repo are UNREPORTABLE** (100% sample overlap with the outcome nested in the
  exposure, correlated/trans/MHC instruments, sign instability). Never import
  `04v3_mvmr_finngen_summary.csv`, `04_mvmr_summary.csv`, or `Table2_PanelC_MVMR.csv` into any manuscript.

## Framing rules (reviewer-proof, locked)
- "Replicated" → "directionally reproduced".
- IGF1R: present as *effector axis*, never as susceptibility anchor. Absence of
  coloc for IGF1R is *not* evidence against its therapeutic role.
- **IGF1R significance is outcome-specific**: nominal in BBJ (P=0.021) and UKB (P=0.012), **not** in
  FinnGen (P=0.182). Never write "nominal across outcomes" — say "directionally consistent, nominally
  significant in BBJ and UKB but not FinnGen".
- **Never call the FinnGen outcome a "case series"** — it is a case–control GWAS with population
  controls. Use "TED-enriched Graves ophthalmopathy outcome".
- **"Prespecified" applies ONLY to the outcome hierarchy** (documented in
  `TrackA_MR/v5_upgrade/03_decision/TaskC_pre_analysis_plan_v1.md`, which covers the BBJ/UKB/FinnGen
  hierarchy and TSHR — *not* IGF1R or CTLA4). The backbone genes were designated **a priori on
  biological and therapeutic grounds**; write it that way, never "prespecified backbone genes".
- **TSHR and IGF1R do NOT share instruments** (TSHR = 1 IV, IGF1R = 4 IVs). What is shared is the
  outcome hierarchy and the analytic framework — never write "identical/same instruments".
- **PP.H2 is NOT "no disease association".** IGF1R has a nominal MR association (BBJ P=0.021, UKB
  P=0.012). PP.H2 = 0.79/0.68 says the *cis*-eQTL signal does not resolve to a variant shared with
  the outcome. Never write "no detectable disease/outcome association" for IGF1R — it contradicts
  the paper's own Table 2 row. Say "does not resolve to a variant shared with the outcome".
- **The null CONSTRAINS, it does not EXCLUDE.** Only 35.6% of genes were powered for OR≥2.0 and
  14.6% for OR≥1.5. Write "constrains additional large expression-mediated effects, particularly
  among well-powered genes, but does not exclude moderate effects" — never "excludes" / "evidence
  against" / "rules out".
- **Fine-mapping was run only at TSHR.** Never write that all three backbone genes went through
  "every evidence layer" — the shared layers are MR, colocalization and orbital tissue.
- **One master, one copy.** `submission/` must not hold a second copy of the manuscript markdown; a
  mirror there went stale once and re-introduced fixed errors. `scripts/audit_paper1_integrity.py`
  now fails if a duplicate reappears.
- Tissue evidence is *exploratory* (single control), never confirmatory.
- **Call the FinnGen outcome "TED-enriched", never "TED-specific"** (it is GO cases vs population
  controls). "TED-specific" survives ONLY where it denotes the *concept* of a TED-specific effect or
  the field's future TED GWAS — never as a label for our outcome.
- **IGF1R wording is comparative, never categorical**: "more compatible with an effector-target
  interpretation than with the TSHR-like expression-colocalized susceptibility architecture" —
  do not write that IGF1R "is not a susceptibility locus".
- External GEO is **not** included (Option A): three external cohorts did not
  reproduce the TSHR tissue direction (lacrimal-enriched / inactive TED). The
  "not externally replicated" limitation is honest and stays. Record:
  `internal/INTERNAL_external_GEO_research_note.md`.

## Working style
- **One task at a time (하나씩).** Verify before moving on.
- When in doubt about a number, check it against the source analysis file under
  `TrackA_MR/v5_upgrade/` — do not fill from memory.
- Figures are built by the R/Python scripts in the repo; **verify figure inputs
  against the locked master** before trusting a render (`FIGURE_VERIFICATION.md`).
  Figure 3B β-axis scale MUST match Table 2.
- Master integrity is tracked by MD5. The master is stored with CRLF line endings
  (`.gitattributes -text`); `scripts/audit_paper1_integrity.py` prints the LF-normalised
  hash. Current master `MANUSCRIPT_TED_TRAP_v5_MASTER.md`: raw `d6b944ecedbe091b3453492ae536304d`,
  normalised `eb72ed12849aaa48236457163244aa05` (placeholders = 0).
- **Length budget (Endocrine Connections):** main text ≤ 5,000 words — currently 4,423 (Introduction–Discussion,
  headings excluded; 4,501 with headings). Abstract is a **single unstructured paragraph**, 245 words (≤250).
  Recount with `python3 scripts/26_wordcount_main_text.py` after any edit.

## Repo layout (actual)
```
submission/        # FINAL JEI package: docx, cover letter, 5 figures, checklist
TrackA_MR/         # v5 core: MR, coloc, fine-mapping, tissue
  v5_upgrade/      #   final analysis results + 07_manuscript/figures (canonical PNGs)
TrackB_Network/    # network / pathway analysis
TrackC_Offtarget/  # off-target + insulin cassette (SEPARATE paper — not v5)
scripts/           # numbered R analysis + download scripts (no data)
finn/  Literature/ Manuscript/  archives/
internal/          # notes NOT for submission (external-GEO, venue, security)
MANUSCRIPT_TED_TRAP_v5_MASTER.md   # markdown source of truth
*.docx             # final submission deliverables (only these are tracked)
```
Local-only (gitignored) raw/reference data lives under each track's `data/`.
