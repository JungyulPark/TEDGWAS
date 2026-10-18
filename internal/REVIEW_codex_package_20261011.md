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

## Resolution (2026-10-14, on top of the imported package 698e4b6)
The author asked to merge on B, reproduce the FinnGen colocalization, restore what was
seriously missing, correct overstatements and remove the AI-assisted wording.

- **FinnGen GD / strict GO colocalization and MR reproduced** for *TSHR*, *IGF1R*, *CTLA4*
  from re-downloaded FinnGen R12 regions with two independent implementations
  (`submission/candidate_20261011_review/Verification/Review_20261014/finngen_coloc_reproduction/`).
  All PP.H4 within 1.2e-3 (identity control 4.1e-4; gnomAD substituted for 1000G EUR),
  no setting crosses 0.80; *IGF1R* PP.H2/H3 within 0.003 (frequency source); MR exact on the
  original exposure scale. The FinnGen GD screen beyond these three genes remains Codex-only.
- **Abstract/Conclusions**: TSHR non-colocalization in UKB (0.226) and FinnGen GD (0.264) and BBJ
  prior sensitivity now in the Abstract; tissue result reduced to "single-variant tissue estimates
  added no independent evidence"; HLA-DQA2/TUBB removed from the Abstract; IGF1R guard restored.
- **Adipose overstatement corrected**: Results/Discussion/Table S3 name rs12101255/rs74064796,
  their intron-1 location, r² 0.46 with rs179252, blood-eQTL status and single-SNP P equivalence.
- **IGF1R**: "consistently indicated higher odds" replaced; FinnGen GD near-null and detectable OR
  ~1.36 stated; alternative estimators reported (weighted median BBJ/UKB, Egger/mode, intercepts);
  two-of-four LOO detail restored; thyroid estimate's opposite direction stated.
- **Restored**: per-allele TSHR ORs and the 0.11 SD scale; PP.H2 column in Table 2; Introduction
  gloss for genetically proxied expression and the RNA-seq disclosure; CTLA4 caveat-before-takeaway;
  "thymic step not tested"; detection benchmarks extended to FinnGen GD (24.8%/52.2%); MHC
  composition of the additional FinnGen GD genes; HLA-A heterogeneity; corrected cross-reference to
  MR_new_freq.
- **Ji et al.**: unverifiable p12 = 5×10⁻⁵ removed (PMC text drops the exponent); PP.H4 > 0.5 kept.
  New references 11 and 29–32 confirmed in PubMed with matching claims.
- **"Secondary" colocalization** wording removed; the post hoc change of emphasis is disclosed.
- **Package**: AI-assisted wording removed (manuscript ×2); Figure 1–3 legends now in
  Manuscript.docx/pdf; internal to-do/version notes removed (supplement, Figure S1 legend, STROBE
  18–20, workbook Readme via XML-only edit); italics and P notation in tables/legends; cover letter
  rewritten for Endocrine Connections; portable build script
  (`reproducibility/build_documents_20261014.py`); framing guards added to the audit.
- **Not resolved / author items**: teprotumumab phase 2/OPTIC-X/chronic-TED references not
  re-added (ref 3 retained with an explicit guard); orbital tissue collection/sequencing methods are
  absent from both A and B and must come from the tissue investigators; STROBE crosswalk remains
  condensed; Figure 1 label collision not re-plotted; companion-manuscript status, journal
  requirements and author approval remain for the authors.

### Verification of the revision (2026-10-14, second pass)
An independent four-lens check (numbers, claims, coverage, rendering; 17 agents) confirmed the new
numbers and found defects introduced or left by the first pass, now fixed: the Abstract had dropped
"Mean" from the orbital sentence (false for one IGF1R sample); prior sensitivity was flagged for BBJ
only (strict GO falls further, 0.891 → 0.451); "two largest datasets" is true only by cases; the
thyroid IGF1R direction was generalised; MR_new_freq was cited for BBJ evidence; antecedent and
naming problems ("its", "Their", undefined "GO"); reproduction-paragraph counts unexplained and the
TSHR strict-GO 0.449 vs 0.451 difference undisclosed; roman gene names, P and flat subscripts in the
supplement; cover letter overstated the STROBE crosswalk and restated unconfirmed declarations;
STROBE heading font and page size. The orbital tissue methods were added to the supplement from the
authors' companion manuscript, pending their confirmation; the companion manuscript's PNUH grant
versus this manuscript's "no specific grant" is flagged for the authors. Figure re-plots (Figure 2
dagger, Figure 1 label collision, in-image titles) remain open.

### Display review of tables and figures (2026-10-14, third pass)
Three independent reviewers checked Tables 1–3, S1–S3 and Figures 1–3, S2 against the workbook.
All displayed numbers matched. One plotting error, inherited from the imported design, was found:
Figure 2B stacked PP.H4 last, so the dotted line at 0.80 marked PP.H4 = 0.20; H4 is now stacked
from zero. Author decisions: Table 2 keeps the numbers (Q-test *P* column added; *CTLA4* IVW kept
with its rs13030124 dependence in the footnote) and Figure 2 shows none; Table 3A is replaced by
Table 3 plus Figure 3, Table 3B became Table S1, Table 3C (a duplicate of S2) was deleted, and the
*TSHR* regional table (S3) was expanded with matched variants, disease lead, rs179252 *P* and
per-allele OR. Supplementary tables are renumbered S1–S4 in first-citation order (S4 tissue).
Tables are now generated from the workbook (`reproducibility/build_tables_20261014.py`); figures
were redrawn at final width with text ≥7 pt. Open: a *TSHR* regional association plot with LD r²
needs the local reference panel.

An independent four-dimension check of this revision (numbers, figures, claims, rendering; each
finding adversarially re-verified) confirmed every displayed number and found 21 defects, all fixed:
the cover letter wrongly said BBJ/strict support held "only under the primary prior" (it holds at
p12 = 5×10⁻⁶ and fails at 10⁻⁶); Table 2 called all *CTLA4* multi-IV estimates heterogeneous (strict
ophthalmopathy Q *P* = 0.0755); Table S4 omitted that UKB single-SNP *P* values are z-test values;
"European" was used for the UKB *CTLA4* estimates; Figure 1 box borders were clipped; Figure 3 and
S2 had non-black lettering; the S1 legend misnamed the log~2~ fold change; missing abbreviations;
three roman *TSHR*/*P* in the supplement; Unicode exponents fell back to a non-Times font (now pandoc
superscripts); a near-empty supplement page; and page numbers restarted inside the merged review
PDFs (now continuous). Gene names in the title and headings stay roman by the earlier convention.

### Results review (2026-10-16)
Four independent reviewers (statistics, clinical reader, editorial, cross-section consistency), each
followed by an adversarial verifier, then one synthesis. Every Results number matched the workbook and
displays. 92 findings survived verification and were merged into 21 corrections and 19 decisions; the
author accepted every recommendation. Corrections included selective reporting (weighted median for
two of three datasets; CTLA4 omission omitting strict ophthalmopathy), the thyroid IGF1R "opposite to
the blood estimates" wording (now "directions were inconsistent"), bare "ophthalmopathy", "corrected
support" for the adipose estimates, and Discussion-style sentences. New in the Results: TUBB/C4A
single-variant dependence, IGF1R H2-dominant colocalization, per-allele ORs, detectable OR 1.36, inactive
TED and the control depot; the C4B/HLA-G reproduction sentence moved to the Supplementary Methods.
Results 1,487 → 1,384 words; main text 4,992 → 4,889 with sub-headings.

### Discussion review (2026-10-17)
Same design as the Results review (evidence and consistency, clinical reader, editorial, citations checked
in PubMed; adversarial verification; one synthesis). All Discussion numbers matched. 83 findings survived
and were merged into 29 corrections, 2 companion edits and 12 decisions; the author accepted every
recommendation. The one factual error: Ji et al. validated TSHR in FinnGen R10 E4_GRAVES_OPHT_STRICT, not
GRAVES_OPHT (Introduction and Discussion corrected). Other changes: removal of overstatement
("support for expression-mediated causality", "further limits confidence"), Results repetition and a
rebuttal of an unstated hypothesis; per-SD ORs explained as lifelong inherited differences; ophthalmopathy
estimates may largely reflect GD susceptibility; paragraph order consolidated around TSHR; [26] moved to
the Discussion and [33] to the Limitations (order 1–33 kept). Discussion 1,739 → 1,427 words; main text
4,889 → 4,573 with sub-headings.

### Results second pass and Discussion rewrite (2026-10-18)
At the author's request the Results were reduced to findings only (1,379 → 1,019 words): interpretation
moved to the Discussion, Limitations, footnotes and legends; numbers shown in displays are cited rather
than repeated, except the TSHR/IGF1R headline estimates; workbook-only values stay. The author then
approved a Discussion outline and the Discussion was rewritten to it (1,427 → 1,162 words), with
findings-only Conclusions without numbers and two new references read in PubMed: Zuber 2022 (abstract
only; supports "a non-zero MR estimate can arise from distinct but correlated variants") and Douglas
2024 (full text; randomized trial in long-duration, low-activity TED, primary outcome proptosis). The
HLA-DQA2 UKB sentence is supported by MR_estimators_original (two instruments, Q 220.1). References now
35, first-cited in order. Main text 3,948 words with sub-headings.
