> **15 September 2026 compact supplementary review:** Leave-one-out is now reported in Methods/Results/Discussion, Figure S2, Data 6 and STROBE item 13. Fifteen omissions comprise 11 IGF1R and 4 CTLA4 results; all nine full-set estimates and every omission were independently checked in base R. The European CTLA4 dependence on rs13030124 is disclosed in Table 2. No new multi-signal colocalization was undertaken: its absence and the single-causal-variant assumption remain explicit limitations. Formal directionality and participant overlap remain unverified. Read the current REVIEW_REPORT_KO.md for review status and audit counts.

> **2026-09-11 reviewed revision:** The current master and `submission/candidate_20260911_review/` supersede the historical numerical claims, fixed scientific assertions, journal status and build instructions below. See the current revision report before editing. Outdated SuSiE/tissue P-value claims must not be restored; IGF1R UKB remains split H2/H4. The requested eQTLGen frequency sensitivity is completed (Tables S2–S3; Supplementary Data 2–4). This is a final author-review candidate, with author details/declarations and live journal checks pending. The non-negotiable data-management rules below remain in force. The 2026-09-05 author decision to omit the manuscript AI-tool declaration is retained; journal disclosure requirements remain an explicit final-review item. "Prespecified" must not describe unregistered thresholds. The old root package remains archived in `archives/submission_pre20260905/`.

# CLAUDE.md — TED-TRAP project rules

This file orients any AI assistant working in this repository. Read it before
making changes. The guiding value of this project is **정확하고 진실한게 생명**
— accuracy and truth above all. No hallucination; verification first.

## What this project is
Druggable-gene-wide Mendelian randomization (MR) + colocalization +
orbital transcriptomics distinguishing **TSHR-anchored genetic susceptibility**
from **IGF1R pharmacologic effector biology** in Graves disease (GD) and
thyroid eye disease (TED). Current state: **v5 leave-one-out integrated revision (14 Sep 2026)**,
prepared for **Endocrine Connections**. The one live package is
`submission/candidate_20260911_review/`; everything else is archived.

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
- **TSHR** — single-instrument locus (rs179252, chr14:81,435,985 hg19) under the
  European selection reference; the East Asian panel gives two. OR 0.12 (BBJ,
  *P* = 1.09×10⁻¹⁴) / 0.09 (UKB, 8.77×10⁻²⁸) / 0.10 (FinnGen, 2.82×10⁻⁷).
  PP.H4 0.951 BBJ / **0.226 UKB** / 0.986 FinnGen — the UKB non-colocalization is
  a headline finding and must never be dropped. The SuSiE fine-mapping layer was
  **withdrawn as invalid**; do not restore credible sets, PIPs or SuSiE r².
- **IGF1R** — 4 instruments (3 in FinnGen). OR 1.56 (BBJ, *P* = 0.0212) / 1.35
  (UKB, 0.0117) / 1.41 (FinnGen, **0.182, not significant**). PP.H4 0.073 / 0.404
  / 0.032 — UKB splits H2 0.400 vs H4 0.404, i.e. **unresolved**, not "eQTL only".
  Write it comparatively (weaker genetic support than *TSHR*), never as "IGF1R is
  not a susceptibility locus" and never as proof of an exclusive effector role.
- **CTLA4** — biologically selected comparator, not a validating positive control:
  PP.H4 0.201 BBJ (H3 0.799) / 0.953 UKB / 0.978 FinnGen, so it **fails** the
  combined BBJ-plus-FinnGen criterion. OR 0.18 / 0.21 / 0.17.
- **eQTLGen per-variant sample size**: max **31,684**, median **28,092**, min **232** across the 6,135
  instruments; 27.4% have n < 20,000. Methods reports the median and range, not just "up to 31,684".
  Instrument strength is unaffected (P < 5e-8 forces F >= 29.7 at any n), and every instrument behind
  the reported backbone genes has n 20,515-31,567.
- **eQTLGen allele-frequency sensitivity is DONE** (2026-09-11). All 13 discovery
  hits retained, no PP.H4 crossed 0.80 (max |Δ| 0.002142), *IGF1R* stayed nominal
  in BBJ/UKB and non-significant in FinnGen. Detection fell: BBJ OR 1.5 power
  14.6% → 12.5%. Never describe this analysis as unperformed again.
- **Leave-one-out is reported** (2026-09-14): 15 omissions = 11 IGF1R + 4 CTLA4; Figure S2, Data 6 and STROBE item 13. European CTLA4 support is concentrated in rs13030124; imprecision of the other SNP is not proof of no effect. New multi-signal coloc is outside this revision; retain model limitations. Do not restore the superseded "not reported" statement.
- **robust_novel = 0** after MHC + chr16p11.2 LD-spillover + cross-outcome coloc
  filtering. This is an *informative* result, not a negative one — do not reframe
  it as a discovery.
- **Tissue n=4 TED + 1 control**, biological-sample level (technical replicates
  collapsed). **Descriptive only** — with one control, control-side variance is not
  estimable, so no *P* value, no DESeq2 result and no significance language may
  appear anywhere. Earlier padj values (0.032, 6.6e-5, 0.006) are all withdrawn.
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
- IGF1R: the genetic evidence is **weaker and unresolved**, not absent. Absence of
  colocalization is *not* evidence against its therapeutic role, and the manuscript
  must never imply that the teprotumumab evidence is challenged by these data.
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
- **PP.H2 is NOT "no disease association".** IGF1R has a nominal MR association (BBJ P=0.0212, UKB
  P=0.0117). A high H2 says the *cis*-eQTL signal does not resolve to a variant shared with the
  outcome. Never write "no detectable disease/outcome association" for IGF1R — it contradicts the
  paper's own Table 2 row. Say "does not resolve to a variant shared with the outcome", and for UKB
  say the posterior is **split** between H2 (0.400) and H4 (0.404).
- **The null CONSTRAINS, it does not EXCLUDE.** Only 35.6% of genes were powered for OR≥2.0 and
  14.6% for OR≥1.5. Write "constrains additional large expression-mediated effects, particularly
  among well-powered genes, but does not exclude moderate effects" — never "excludes" / "evidence
  against" / "rules out".
- **Fine-mapping is withdrawn entirely.** The SuSiE run had an allele-harmonisation defect; no
  credible set, PIP, purity or SuSiE-derived r² may appear. The evidence layers are MR and
  colocalization; the orbital tissue is descriptive context only.
- **One master, one package.** `submission/` must not hold a second copy of the manuscript markdown
  (a mirror went stale once and re-introduced fixed errors), and must hold exactly **one**
  `candidate_*` directory — three accumulated once, each with a different manuscript, cover letter
  and figure set. Superseded packages go to `archives/submission_candidates/`.
  `scripts/audit_paper1_integrity.py` fails on either.
- Tissue evidence is *exploratory* (single control), never confirmatory.
- **Call the FinnGen outcome "TED-enriched", never "TED-specific"** (it is GO cases vs population
  controls). "TED-specific" survives ONLY where it denotes the *concept* of a TED-specific effect or
  the field's future TED GWAS — never as a label for our outcome.
- **IGF1R wording is comparative, never categorical**: "weaker genetic support for *IGF1R* than for
  *TSHR* in the two Graves disease outcomes; they do not show that *IGF1R* has no inherited
  contribution" — do not write that IGF1R "is not a susceptibility locus".
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
  The current set is Figures 1-3 plus Figures S1-S2; Figure 2's OR axis MUST match Table 2. All figure text must remain black; asterisks denote nominal P<0.05, with an additional dagger only for BBJ discovery P<0.05/2,544. Figure3 borders denote PP.H4>=0.80, not P-value significance.
- **Run `python3 scripts/audit_paper1_integrity.py` after every manuscript edit.** It runs the
  framing guards (each one encodes a defect that was actually removed) and then the candidate's
  numeric audit, which compares every displayed value with its full-precision source.
- Master integrity is tracked by MD5. The master is stored with CRLF line endings
  (`.gitattributes -text`); `scripts/audit_paper1_integrity.py` prints the LF-normalised
  hash. Current master `MANUSCRIPT_TED_TRAP_v5_MASTER.md`: raw `58f81a603dff5d1605c562fb7389e35d`,
  normalised `7361d2d9bdaa8dd02fc8626ec1b19a4c` (placeholders = 0). **33 references** since the 1000 Genomes panel got its own.
- ***P* is always italic** — in "*P* value" and "*P* <" too, not only before a number.
- **Length budget (Endocrine Connections):** main text ≤ 5,000 words — currently **3,796**
  (Introduction–Discussion; 3,828 with sub-headings, which is what Word reports). Abstract is a
  **single paragraph** with inline `Objective:/Methods:/Results:/Conclusions:` labels, **250** words
  (≤250). `scripts/26_wordcount_main_text.py` is the number of record.
  Recount with `python3 scripts/26_wordcount_main_text.py` after any edit.
- **Introduction is written for a first-time MR/colocalization reader** (medical-student/resident level,
  2026-09-15): MR and colocalization are each glossed with their core intuition on first use,
  "genetically proxied expression" and "druggable genes" are defined at first use, and one sentence
  previews that MR results are OR/95% CI/*P* while colocalization results are a posterior probability,
  not a *P* value. **Paragraph 3 stays impersonal** — it defines the two methods generically, with no
  "we used X" statement and only a generic "exposure"; every first-person statement of what this study
  did belongs in the final aims paragraph, which opens "Here, we...". That aims paragraph also
  discloses the in-house orbital RNA-seq as descriptive context only — the Abstract announces the
  dataset, so the Introduction must not stay silent about it, but it is never framed as an evidence
  layer, a validation or a result.
  Keep this register if the Introduction is edited again; do not revert to the terser
  original phrasing.

## Repo layout (actual)
```
submission/candidate_20260911_review/   # THE package: docx, cover letter,
                   #   STROBE-MR checklist, Figures 1-3 + S1-S2, Supplementary Data 1-6,
                   #   provenance/ (incl. maf/) and reproducibility/
archives/submission_candidates/         # superseded packages -- never submit from here
TrackA_MR/         # v5 core: MR, coloc, tissue (fine-mapping withdrawn)
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
