# TEDGWAS: current manuscript and handoff (14 October 2026)

Read `HANDOFF_TO_CLAUDE_20261014.md` and
`internal/REVIEW_codex_package_20261011.md` before editing. The author asked Codex
to commit the completed local revision so Claude can continue the review. This
handoff imports the 3 October analysis / 11 October delivery without scientific
rewriting. The 13 October review items were then addressed on 14 October (see the
review-status section below).

## Current source of truth
- Sole editable manuscript: `MANUSCRIPT_TED_TRAP_v5_MASTER.md` (exact local final bytes).
- One active package: `submission/candidate_20261011_review/`.
- Main Tables 1–3 (Table 3A–C), Figures 1–3; Tables S1–S3 and Figures S1–S2;
  one `Supplementary_Dataset.xlsx`. Table S3 is now the tissue extension;
  historical detection results are in the workbook's `Detection_original` sheet.
- Old candidate and earlier instructions are historical. Do not restore their
  three-outcome framing, six separate Data files, or obsolete sample sizes.
- Read-only current check: `python scripts/audit_paper1_integrity.py`.
  It checks delivery integrity, tissue formulas and 52 framing guards (added 14 October).

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

## Review status and editing (updated 14 October 2026)
The 13 October review items were addressed on top of the imported package; see the
"Resolution" section of `internal/REVIEW_codex_package_20261011.md`. Remaining author items:
journal requirements and all-author approval.

Locked from this revision:
- The Abstract reports TSHR PP.H4 in UKB (0.226) and FinnGen GD (0.264) and that BBJ and strict
  GO support is prior-sensitive (0.661 and 0.451 at p12 = 1e-6). Never again reduce this to
  "support varied". Orbital data are "Mean ... exceeded the control" — IGF1R is higher in only 3 of 4
  TED samples, so never drop "Mean".
- FinnGen GD/strict GO colocalization and MR for TSHR/IGF1R/CTLA4 were independently reproduced
  (`Verification/Review_20261014/finngen_coloc_reproduction/`): PP.H4 within 1.2e-3, MR exact on
  the original exposure scale; IGF1R PP.H2/H3 only within 0.003 (frequency source). The FinnGen GD
  screen beyond these three genes is Codex-computed, not independently reproduced.
- Adipose TSHR instruments rs12101255/rs74064796 sit 50 bp apart in TSHR intron 1, ~15 kb from
  rs179252 (EUR r² 0.46), and are blood TSHR eQTLs (P ~1e-34). Single-SNP P = the variant's GWAS P.
  Report them as re-expressing the regional association, never as independent support.
- IGF1R: all ORs > 1, nominal only in BBJ/UKB; FinnGen GD OR 1.07 (0.86–1.33) is the most precise
  and near-null (detectable OR ~1.36 at P < 0.05). Never write "consistently indicated higher odds".
  Weighted median BBJ 1.56 (P 0.0288), UKB 1.42 (P 0.00244); Egger intercepts P 0.66/0.56/0.43.
- Per-allele TSHR ORs 0.80/0.77/0.78/0.77/0.82 (BBJ/UKB/GO/GD/strict GO) explain the per-SD ORs;
  call them moderate, not small. The main text writes "FinnGen ophthalmopathy"/"strict ophthalmopathy";
  "GO" is defined only in the Abstract and table footnotes.
- Thyroid IGF1R: only the primary FinnGen GD estimate is below 1 (0.76); the four secondary ones are
  above 1. Write "directions were inconsistent", never "opposite to blood" without that scope.
- Ji et al. [9]: PP.H4 > 0.5 is confirmed; their p12 exponent is unreadable in PMC — do not state it.
- Colocalization is not called "secondary"; the post hoc change of emphasis is disclosed in Methods.
- No AI-assisted wording in any submission-facing file (author decision 2026-09-05, reaffirmed 2026-10-14).
- Abstract 248/250; main text 4,951 (4,998 with sub-headings) of 5,000 — effectively no headroom.
- Author decisions (14 October): no grant funded this work, so "no specific grant" stands; the
  companion orbital whole-transcriptome manuscript is not under review; orbital tissue was
  collected during decompression from patients with inactive TED. Keep the orbital data brief and
  secondary (one short supplementary paragraph; descriptive only).

Editing workflow:
1. Edit `MANUSCRIPT_TED_TRAP_v5_MASTER.md` and/or the package Markdown sources
   (`Tables.md`, `Figure_Legends.md`, `Supplementary_Material.md`, `SupplementaryTables_S1_S2.md`,
   `TableS3_Tissue_extension.md`, `Cover_Letter_Draft.md`).
2. Rebuild: `python3 submission/candidate_20261011_review/reproducibility/build_documents_20261014.py`
   (pandoc + python-docx + LibreOffice + poppler) and inspect rendered pages.
3. If new numbers were added, extend `reproducibility/derive_review_values_20261014.py` and rerun it.
4. `python3 scripts/26_wordcount_main_text.py`; then `python3 scripts/29_refresh_manifest.py YYYY-MM-DD`;
   then `python3 scripts/audit_paper1_integrity.py` (manifest, tissue formulas, 52 framing guards).
Do not edit scientific output CSVs or workbook data sheets to make an audit pass; the workbook's
Readme sheet was changed by editing only `xl/worksheets/sheet1.xml`. Keep exact reviewed bytes under
the package's `-text` attribute.
