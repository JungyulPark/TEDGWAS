# Review of the Codex package TED_Final_Review_20261011 (internal)

Reviewed 2026-10-13 from the ZIP the author uploaded (`TED_Final_Review_20261011.zip`, 220 entries,
MANIFEST_SHA256 verified: 219 files, 0 missing, 0 mismatched). Seven independent reviewers, an
adversarial check of every critical/major finding, and a completeness critic. The package is **not on
GitHub**: as of this review the newest remote commit is still `c46a241` (2026-09-16).

Call the GitHub master **A** (`MANUSCRIPT_TED_TRAP_v5_MASTER.md` at `c46a241`) and the package **B**.

## Lineage
B descends from the unmerged branch `origin/codex/abstract-clinical-refinement-20260915` (`3468bf3`),
then adds FinnGen GD (`E4_GRAVES_STRICT`), strict Graves ophthalmopathy (`E4_GRAVES_OPHT_STRICT`), a
24-gene colocalization (1,440 settings), a TSHR allele audit and a fixed five-tissue GTEx v8 extension.
B lacks most of A's 2026-09-15/16 wording fixes but already contains the eQTLGen per-variant n sentence
and an impersonal MR/colocalization paragraph with the OR-vs-posterior preview.

## Verified (recomputed from primary inputs)
- All 20 estimable GTEx tissue Wald ratios, from raw GTEx significant-pair/eGene extracts and raw disease
  extracts with independent allele harmonisation: max relative difference 4e-14. F 32.1–61.1, the 50-cell
  grid, the 30 unestimable verdicts, the Bonferroni families (0.005 / 0.00125) all correct.
  Protocol freeze 11:08:50Z precedes GTEx retrieval 11:47Z (self-attested).
- TSHR rs179252 in all five outcomes from source rows: MR β −2.09595 / −2.43623 / −2.33149 / −2.52690 /
  −1.90408 (BBJ / UKB / FinnGen GO / FinnGen GD / strict GO). FinnGen GD OR 0.0799 (0.0525–0.1215),
  P 3.19×10⁻³²; strict GO OR 0.149 (0.057–0.386), P 9.0×10⁻⁵.
- Every number carried over from A is identical to A's provenance (7,219 MR rows, 13,039 estimator rows,
  30,115 frequency rows, 324 coloc rows, power, LOO, 6,135 instruments).
- 307 values in Tables 2, 3A–C, S1, S2 and Figures 2, 3, S2 match the workbook.
- Allele audit: signed LD rs179252-G / rs179247-G r = 0.978 (EUR) / 0.996 (EAS), G-G coupling.
  rs179252-G raises blood TSHR (Z +13.28) and lowers odds in all five outcomes (per-allele OR
  0.766–0.818). Marín-Sánchez 2019 (full text) names rs179247-G as GD-protective and preferentially
  transcribed in thymus, so A's conditional central-tolerance argument is directionally supported.
  (The Colobran 2011 allele label could not be read; egress blocked.)
- 33 references, all cited, numbered in first-citation order. B resolves A's OPEN [9, 10] TNFSF14/IFNGR1
  citation issue (the sentence is gone; [10] is cited only for inflammatory proteins).
- Words: Abstract 247, main text 4,555.

## Not verifiable here
- FinnGen GD / strict GO colocalization (1,116 of 1,440 settings) and the FinnGen GD full screen beyond
  TSHR: no E4 regional summary statistics in B or the repo. Reproduce before quoting in the Abstract.
- Five reference replacements and their attributed claims (Ji et al. p12 5×10⁻⁵ / PP.H4 > 0.50; new refs
  11, 29–32 incl. two 2026 papers). Endocrine Connections limits and AI-disclosure policy.

## Main problems in B (adversarially confirmed)
1. **Abstract hides the TSHR non-colocalization.** UKB PP.H4 0.226 (H3 0.774) and FinnGen GD 0.264
   (H3 0.736) — the two best-powered European GD/hyperthyroid datasets — became "Secondary colocalization
   support varied across datasets", while FinnGen GD MR P 3.19×10⁻³² and the adipose result were promoted.
2. **Adipose extension overstated in Abstract and Conclusions.** One SNP per tissue: rs12101255
   (subcutaneous) and rs74064796 (visceral, 50 bp away) in the TSHR intron-1 disease-lead cluster;
   rs12101255 has EUR r² 0.459 with rs179252 and is itself an eQTLGen blood TSHR eQTL (P 4.1×10⁻³⁴).
   Single-SNP MR P equals that SNP's GWAS P (1.28×10⁻³⁴ both), so "passed correction" was nearly
   guaranteed; no tissue colocalization; same disease GWAS reused. Not independent support.
3. **IGF1R "consistently indicated higher odds"** glosses over FinnGen GD OR 1.07 (0.86–1.33), P 0.539,
   the most precise IGF1R estimate (detectable OR at α 0.05 ≈ 1.36 vs 1.72 BBJ, 1.39 UKB).
4. HLA-DQA2 / TUBB given Abstract space although PP.H3 ≈ 1.000 in all five datasets (MHC, low-n
   instruments). "Five datasets" concordance counts three overlapping FinnGen endpoints
   (GO vs strict GO gene-level z correlation 0.947).
5. Colocalization relabelled "secondary" after TSHR failed in FinnGen GD. B discloses it as post hoc
   (verifier: minor), but a reviewer may still question it.

## Lost relative to A
Per-allele TSHR ORs and the 0.105 SD shift (OR 0.08 left unexplained); the teprotumumab-trial paragraph
and refs (phase 2/3, OPTIC-X, chronic TED) plus Hemani 2017 (Steiger); PP.H2 column in Table 2 (IGF1R UKB
H2/H4 split invisible); alternative MR estimators (named in Methods, reported nowhere); detection-limit
table (only in an XLSX sheet; 35.6% figure gone; nothing computed for FinnGen GD — recomputed 24.8% /
52.2% at OR 1.5 / 2.0); Introduction lay glosses and RNA-seq disclosure; Conclusions IGF1R guard;
"tissue-specific step remains untested" sentence; CTLA4 takeaway again before its caveat. Results cite
Table 3C/S2 for TSHR/IGF1R frequency stability but neither table contains those genes.

## Package defects
AI-assisted declaration added twice (lines 74 and 162, inconsistent wording) against the 2026-09-05 author
decision; Figure 1–3 legends missing from every Word/PDF deliverable; internal to-do/version-history notes
in the supplement, Figure S1 legend, STROBE items 18/20 and the workbook Readme; cover letter without
journal, date or disclosures; no italics in Tables.docx / Supplementary_Material.docx; lost SE subscript
("SEX"); Figure 2 dagger dropped; Figure 1 label collision; STROBE collapsed to 20 rows.
Data protection: no IRB raw data, patent material or full licensed summary statistics. Contains
open-access 1000 Genomes genotype dosages at two SNPs and local username/host paths in 29 log files —
do not commit those as-is.

## Placement judgement
Agree with Codex: no main table/figure for the tissue extension; allele audit in text plus a workbook
sheet; FinnGen GD in Table 2/Figure 2 is a real improvement (restore PP.H2). Disagree: tissue result in
Abstract/Conclusions; dropping detection limits from the supplement; Table 3 at three pages with
duplication (3C = S2; 3A/3B PP.H4 repeat Figure 3A).

## Recommended merge
Base on B's expanded, numerically verified analysis; port A's wording fixes; rebalance the Abstract and
Conclusions; restore the lost items above; reproduce FinnGen GD/strict GO colocalization for TSHR, IGF1R
and CTLA4 from FinnGen R12 regional data before quoting it; get Codex's working repository pushed so
scripts and provenance come from source rather than the ZIP.
