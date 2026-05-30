# [Gemini → User] v5 Manuscript Update & Verification Confirmation (2026-05-27)

We have verified, implemented, and synchronized all required updates across the repository and manuscript folder. Below is the final verification checklist.

## 1. Table & Text Corrections (Reconciliation)
- **Table 3 *PSMB8* CI Upper Bound:** Updated to `1.25 (1.13–1.38)` (raw `1.3753` rounded to `1.38`) on line 253 of `MANUSCRIPT_TED_TRAP_SUBMISSION.md`.
- **Table S1 *IFNGR1* FinnGen Beta:** Updated to `+0.14 (0.68)` (raw `+0.145` rounded to `+0.14` under the R round-to-even behavior) on line 280 of `MANUSCRIPT_TED_TRAP_SUBMISSION.md`.
- **Table 2 Footnote Threshold:** Normalized to `PP.H4 ≥ 0.80` on line 236 of `MANUSCRIPT_TED_TRAP_SUBMISSION.md` to establish absolute consistency with the Methods section.

## 2. Table S7. STROBE-MR Checklist (TableS7_STROBE_MR.md)
- Replaced the older version in the manuscript folder (`TrackA_MR/v5_upgrade/07_manuscript/TableS7_STROBE_MR.md`) with the user's updated checklist (`SuppTable_S7_STROBE_MR.md`).
- **Checklist modifications verified:**
  - Title corrected to **Table S7**.
  - Sensitivity table references corrected to **S6** (4 locations).
  - Item 22 (Data availability) GitHub link removed and replaced with public repositories (*eQTLGen, GWAS Catalog, FinnGen*).
  - Software versions specified (*R 4.3.3; TwoSampleMR v0.7.4; ieugwasr; coloc v5.2.3; susieR v0.14.2; DESeq2 v1.42.1; PLINK v1.9*).
  - Colocalization findings updated to exact probabilities (*PP.H4 = 0.951 (BBJ) / 0.986 (FinnGen)*).
- Obsolete S8 files removed.

## 3. Acknowledgements & Declarations
- The complete **Declarations** section (including Acknowledgements with the FinnGen standard disclaimer, GWAS Catalog, and Sakaue 2021 cohort citations; Funding; Author contributions; and Ethics approval) is fully integrated into the manuscript (`MANUSCRIPT_TED_TRAP_SUBMISSION.md` starting at line 150).

## 4. Cover Letter (Thyroid/ETJ Restructuring Statement)
- Updated `COVER_LETTER.md` to explicitly state that this manuscript represents a completely restructured and refocused study (previously scoped for *Thyroid* or the *European Thyroid Journal*) refined to concentrate strictly on the genetic triangulation of druggable targets using high-rigor MR, colocalization, and target-relevant orbital tissue expression.

## 5. Main Text Word Count
- The word count of the main body (Introduction, Methods, Results, Discussion) is **4,380 words**, which is safely within the *Journal of Endocrinological Investigation* limit of **4,500 words** (leaving 120 words of margin).

## 6. Directory Synchronization & DOCX Compilation
- **Markdown files synced:** The updated master markdown files are synced in:
  - `c:\ProjectTEDGWAS\MANUSCRIPT_TED_TRAP_SUBMISSION.md`
  - `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\07_manuscript\MANUSCRIPT_TED_TRAP_SUBMISSION.md`
  - `c:\ProjectTEDGWAS\COVER_LETTER.md`
  - `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\07_manuscript\COVER_LETTER.md`
- **DOCX compiled:** Recompiled `.docx` files generated and synchronized in:
  - `c:\ProjectTEDGWAS\MANUSCRIPT_TED_TRAP_SUBMISSION.docx`
  - `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\07_manuscript\MANUSCRIPT_TED_TRAP_SUBMISSION.docx`
  - `c:\ProjectTEDGWAS\COVER_LETTER.docx`
  - `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\07_manuscript\COVER_LETTER.docx`
