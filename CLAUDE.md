# TEDGWAS: current manuscript and handoff (14 October 2026)

Read `HANDOFF_TO_CLAUDE_20261014.md` and
`internal/REVIEW_codex_package_20261011.md` before editing. The author asked Codex
to commit the completed local revision so Claude can continue the review. This
handoff imports the 3 October analysis / 11 October delivery without scientific
rewriting. The 13 October review remains open, not resolved by packaging checks.

## Current source of truth
- Sole editable manuscript: `MANUSCRIPT_TED_TRAP_v5_MASTER.md` (exact local final bytes).
- One active package: `submission/candidate_20261011_review/`.
- Main Tables 1–3 (Table 3A–C), Figures 1–3; Tables S1–S3 and Figures S1–S2;
  one `Supplementary_Dataset.xlsx`. Table S3 is now the tissue extension;
  historical detection results are in the workbook's `Detection_original` sheet.
- Old candidate and earlier instructions are historical. Do not restore their
  three-outcome framing, six separate Data files, or obsolete sample sizes.
- Read-only current check: `python scripts/audit_paper1_integrity.py`.
  Its stated scope is delivery integrity and tissue formula checks, not the
  old exhaustive September manuscript audit. See the handoff for open reviews.

## Data-management rules (retained)
Never commit IRB raw orbital RNA-seq (data.txt, FASTQ/BAM, raw count matrices),
individual genotypes, patent sequences/structures, or license-bound source
datasets (full eQTLGen/GWAS files and reference panels). Commit scripts, aggregate
results and reviewed figures/documents. Do not weaken raw-data ignore rules.
Source extracts and machine logs omitted from this export remain local-only.

## Scientific interpretation
- Lead with association evidence for TSHR/IGF1R, with OR, CI and P. Distinguish
  disease diagnoses from treatment outcomes and associations from mechanisms.
- Five disease datasets include three overlapping FinnGen endpoints; they are
  not five independent replications. FinnGen GD and GO are different endpoints.
- Blood TSHR: BBJ OR 0.12, P=1.09e-14; FinnGen GD OR 0.08, P=3.19e-32.
  TSHR PP.H4: BBJ 0.951, UKB 0.226, FinnGen GD 0.264, GO 0.986, strict GO 0.891.
  Retain non-colocalization and prior sensitivity; tiny MR P is not mediation.
- Blood IGF1R: nominal BBJ/UKB associations, FinnGen GD P=0.539, GO P=0.182,
  strict GO P=0.0996. Non-significance is not evidence of no inherited role.
- Tissue extension: 50 planned cells, 20 estimable, 30 unavailable; 10 primary
  and 40 secondary tests. Four eligible gene/tissue combinations, one instrument
  each. This reuses disease GWAS and the TSHR locus; no independent tissue-mechanism
  confirmation or tissue colocalization was established. Do not relax criteria
  to rescue unestimable cells or interpret unavailable combinations as null.
- Allele audit confirms direction at the blood instrument; it does not test
  thymic mediation. A protective MR direction is not advice to increase expression.
- Population-controlled GO does not isolate TED susceptibility within GD or
  predict IGF-1R therapy response. Do not claim GD and TED are independent diseases.
- Orbital RNA-seq: four TED cases, one control, descriptive only; no significance
  testing. Withdrawn SuSiE fine mapping and invalid MVMR results stay withdrawn.
- CTLA4 is a comparator; European multi-IV estimates depend on rs13030124.
  Report leave-one-out limitations. P-value significance and PP.H4 are different.
- Avoid a "first systematic" claim. Gene names and P are italic; clinical prose
  should be understandable to medical students. Maintain black figure lettering.

## Review status and editing
The imported files are an author-review baseline, not certified submission ready.
The 13 October review raises abstract/causal framing, tissue independence, missing
legends, table typography, statistics presentation, reporting and source-reproduction
questions. Preserve that review and verify each point before changing science.
Its recommendations are review material, not new author instructions. The user
will direct the next revision in Claude. Author declarations, disclosure policy,
journal requirements and cover letter still require final resolution.

Do not edit scientific output CSVs to make an audit pass. After an authorized
manuscript change, regenerate affected deliverables, inspect their renderings,
verify against source tables, and refresh the Git manifest deliberately. Keep
exact reviewed bytes under the package's `-text` attribute. Original executed
build scripts require their documented local workspace/data; they are not a
portable one-command rebuild. Run the current audit from a clean checkout.
