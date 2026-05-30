# TED-TRAP v5 — Manuscript Assets Lock Registry

**Purpose:** Track which manuscript tables/figures are verified and locked, with the
exact values, source files, and lock status. Nothing enters the manuscript unless it
appears here as LOCKED.

**Lock procedure:** Claude verifies values against source data → user approves →
status set to LOCKED → Gemini records file path + MD5 in this registry → asset frozen.

---

## STATUS LEGEND
- 🔄 IN PROGRESS — being verified/edited
- ✅ LOCKED — verified, approved, frozen (do not change without re-opening)

---

## MAIN TABLES

### Table 1 — Study design & data sources
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** hand-built (study metadata), no computed values
- **File:** `07_manuscript/tables/Table1_study_design.{html,png}` and `Table1_VALUES.csv`
- **Verified values (all ✅ correct):**
  - eQTLGen: EUR, 31,684 (→ "up to 31,684 participants"; rs179252 NrSamples=31,566)
  - BBJ Graves GCST90018627: EAS, 2,809 cases (discovery)
  - UKB hyperthyroid GCST90038636: EUR, 3,731 cases (replication)
  - FinnGen R12 GO: EUR, 858 cases (TED-specific sensitivity)
  - In-house orbital RNA-seq: Korean, 4 TED + 1 control
- **MD5:**
  - HTML: f595d88bc2485c08dd935ba383bc978a (size: 12274 bytes)
  - PNG: a07958406e016140ed816a386067bfcb (size: 82986 bytes)
  - CSV: 879a58cfdf85c39d3ece905a1918185a (size: 425 bytes)

### Table 2 — Integrated evidence for backbone genes
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskD_03d` (MR), `TaskE_01` (coloc), `TaskF_01` (tissue) — all ✅ Claude-verified
- **File:** `07_manuscript/tables/Table2_integrated_evidence.{html,png}` and `Table2_VALUES.csv`
- **Verified values (all ✅ correct, re-confirmed):**
  | Gene | MR BBJ | MR UKB | MR FinnGen | Coloc H4 BBJ/FinnGen | Tissue log2FC (padj) |
  |------|--------|--------|------------|----------------------|----------------------|
  | TSHR | −2.10 (1.1e-14) | −2.44 (8.8e-28) | −2.33 (2.8e-7) | 0.951/0.986 | +2.33 (0.032) |
  | IGF1R | +0.45 (0.021) | +0.30 (0.012) | +0.34 (0.182) | 0.043/0.036 | +0.41 (1.00 NS) |
  | CTLA4 | −1.74 (5.5e-15) | −1.57 (0.002) | −1.77 (0.012) | 0.206/0.978 | +1.27 (0.815 NS) |
- **MD5:**
  - HTML: 7b45adcbfe6b49fa47c874bd6a5891dd (size: 13430 bytes)
  - PNG: e948f45ae01fa9d0d6e7080be8574a98 (size: 92296 bytes)
  - CSV: 923db5eeb5f225e319563e22eac9bb22 (size: 568 bytes)

---

## MAIN FIGURES

### Figure 1 — Study design and analytical workflow
- **Status:** ✅ LOCKED (2026-05-27)
- **Source:** `build_figure1_studydesign.R` — ✅ verified (3-band structure: input, analytical streams, integration)
- **File:** `07_manuscript/figures/Figure1_study_design.{png,pdf,tiff}`
- **Verified:** standard ASCII text labels (no hollow boxes), monochrome theme matching PPTX/DOCX template, accurate representation of inputs and analytical streams for TSHR, IGF1R, and CTLA4.
- **MD5:**
  - PNG: 03f318da6835b35e0c0f802c94f50733 (size: 678974 bytes)
  - PDF: 911bc8dbcc8dab85d8d4bfacee45dedb (size: 58946 bytes)
  - TIFF: a06d73628ad2b402ee26e7a0e2813b25 (size: 917636 bytes)

### Figure 2 — MR volcano
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskD_03d` — ✅ verified (2,235 genes plotted, 13 Bonferroni hits)
- **File:** `07_manuscript/figures/Figure2_MR_volcano.{png,pdf}`
- **Verified:** TSHR/CTLA4 upper-left (protective, strong); IGF1R mid-low (nominal, below Bonferroni); 13 hits labeled; class colors correct; Bonferroni line at 1.965e-5. Plotted title displays "Druggable-gene-wide BBJ discovery MR".
- **MD5:**
  - PNG: ecfa7368a8a20826fe1ae05634bbe919 (size: 220210 bytes)
  - PDF: 0e8cd58e099f17b0971041b90794f67c (size: 129913 bytes)

### Figure 3 — TSHR vs IGF1R multi-modal composite
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskE_01` (coloc), `TaskD_03d` (MR), `TaskF_01` (tissue n=5) — all ✅
- **File:** `07_manuscript/figures/Figure3_TSHR_vs_IGF1R_composite.{png,pdf}`
- **Verified:** Panel A coloc (TSHR H4 green / IGF1R H2 orange [labeled: H2 (GWAS assoc. only)] / CTLA4 mixed); Panel B forest (TSHR/CTLA4 negative, IGF1R positive); Panel C tissue (TED n=4 + Control n=1, TSHR padj=3.2e-2). All values match source.
- **MD5:**
  - PNG: b388760c8f7468c1274121c6d4b10d13 (size: 233142 bytes)
  - PDF: 2892925f9785012e921228e84015df59 (size: 9713 bytes)

---

## CAPTIONS & TEXT

### Main Captions & Methods Sentences
- **Status:** ✅ LOCKED (2026-05-26)
- **File:** `07_manuscript/MAIN_CAPTIONS.md`
- **MD5:** 3906c4294b42ca6d3a7c5709cf546a45 (size: 6944 bytes)
## SUPPLEMENTARY TABLES

### Table S1 — 13 discovery hits classification
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskD_04` — verified (13 Bonferroni hits classified)
- **File:** `07_manuscript/supp_tables/TableS1_mr_hits.{csv,html,png}`
- **MD5:**
  - CSV: 5c13b61649273805ad11a60dbeb049e2 (size: 1497 bytes)
  - HTML: 72753ccc7a753670dd8c59c11cccd9bc (size: 23100 bytes)
  - PNG: 61728b4c10de7a620e2ebdcd0bcc7dd2 (size: 183431 bytes)

### Table S2 — Candidate colocalization results
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskE_02` — verified ( colocalization PP for candidate genes )
- **File:** `07_manuscript/supp_tables/TableS2_coloc.{csv,html,png}`
- **MD5:**
  - CSV: fc997091f0bd595459ff6cfc51abe36a (size: 1393 bytes)
  - HTML: 2f34e90a138775fcb9c7126219f8b2c4 (size: 19375 bytes)
  - PNG: ffe12d9a6066e42560ab12c8662fb871 (size: 152971 bytes)

### Table S3 — TSHR fine-mapping credible set
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskB_04` — verified (SuSiE credible sets from EUR & EAS reference panels)
- **File:** `07_manuscript/supp_tables/TableS3_finemapping.{csv,html,png}`
- **MD5:**
  - CSV: 1465b422e423bd0283adf298298aa10a (size: 431 bytes)
  - HTML: 0724e385d3c35e343d0fbfef18a94aae (size: 14346 bytes)
  - PNG: a9d930ff7a541d56e25e5a62ff60b0d4 (size: 73816 bytes)

### Table S4 — eQTL instruments & F-statistics
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskD_02a` — verified (6,136 instruments, min F-stat = 14.20084, 0 weak instruments)
- **File:** `07_manuscript/supp_tables/TableS4_instruments.{csv,html,png}`
- **MD5:**
  - CSV: 9a8b74c09fc9d31aff744b1dd41fe4fa (size: 530433 bytes)
  - HTML: d758ee89178871d45ca6e7abc4319643 (size: 39201 bytes)
  - PNG: bff8ca5ecc1d7e5e00d8881fe94c7d61 (size: 285586 bytes)

### Table S5 — chr16p11.2 LD cluster
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskD_04` — subset of S1 with lead SNP co-localization (HSD3B7, VKORC1, PRSS36)
- **File:** `07_manuscript/supp_tables/TableS5_chr16_cluster.{csv,html,png}`
- **MD5:**
  - CSV: 3487116233b09979182daf5ba58a840f (size: 450 bytes)
  - HTML: 77eecf4ba19a4f168b7fbaa6a064abce (size: 13600 bytes)
  - PNG: 27fede03b678102cc22e29560121020e (size: 73904 bytes)

### Table S6 — MR Sensitivity / Heterogeneity for Backbone Genes
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** `TaskD_03d` — verified (heterogeneity & Egger intercept tests for TSHR, IGF1R, CTLA4)
- **File:** `07_manuscript/supp_tables/TableS6_mr_sensitivity.{csv,html,png}`
- **MD5:**
  - CSV: 8f876306fc3ee8d5b17ba85181677f2c (size: 1694 bytes)
  - HTML: 74e028f7f66d48d814d8139f6c62587b (size: 20882 bytes)
  - PNG: 5f4e09ef905c6f8dadf5fbc6d73bc36f (size: 190319 bytes)

### Table S7 — STROBE-MR checklist
- **Status:** ✅ LOCKED (2026-05-26)
- **Source:** hand-built (reporting standards)
- **File:** `07_manuscript/TableS7_STROBE_MR.md`
- **MD5:** 3edb86275c81b8ae0335a5c0788236b9 (size: 5962 bytes)

---

## SUPPLEMENTARY FIGURES

### Figure S1 — TSHR fine-mapping & LD rationale
- **Status:** ✅ LOCKED (2026-05-27)
- **Source:** `build_figureS1_finemapping.R` (based on Table S3 SuSiE credible sets)
- **File:** `07_manuscript/figures/FigureS1_TSHR_finemapping_LD.{png,pdf}`
- **MD5:**
  - PNG: 98f379bbd1796c2cc368e66baa65aa78 (size: 188335 bytes)
  - PDF: 5b9a7736cc66108b7ecede1051ca031f (size: 5803 bytes)

### Figure S2 — Candidate colocalization profiles
- **Status:** ✅ LOCKED (2026-05-27)
- **Source:** `build_figureS2_coloc.R` (based on Table S2 candidate colocalization results)
- **File:** `07_manuscript/figures/FigureS2_candidate_coloc_failure.{png,pdf}`
- **MD5:**
  - PNG: 5885f8129bce8d7ae7405a4cea23ae69 (size: 112466 bytes)
  - PDF: 3c66a8ed0b7660fa97af4db6482f37e2 (size: 5887 bytes)

---

## LOCK LOG
| Date | Asset | Action | By |
|------|-------|--------|-----|
| 2026-05-26 | Table 1 | LOCKED — wording updated, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Table 2 | LOCKED — wording updated, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Fig 2   | LOCKED — title updated, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Fig 3   | LOCKED — legend label updated, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Captions| LOCKED — MAIN_CAPTIONS.md written and MD5 registered | Claude+user |
| 2026-05-26 | Table S1| LOCKED — split to 4 same-dir/nominal columns, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Table S2| LOCKED — coloc table generated and MD5 registered | Claude+user |
| 2026-05-26 | Table S3| LOCKED — column names and ancestry-matched footnote updated, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Table S4| LOCKED — title and footnotes updated, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Table S5| LOCKED — lead SNP, lead SNP position columns, and softened footnote added, rebuilt and MD5 registered | Claude+user |
| 2026-05-26 | Table S6| LOCKED — sensitivity table generated, renumbered from S7, and MD5 registered | Claude+user |
| 2026-05-26 | Table S7| LOCKED — STROBE-MR checklist renumbered from S8, toned-down language updated and MD5 registered | Claude+user |
| 2026-05-27 | Fig 1   | LOCKED — R study design schematic generated, standard ASCII, MD5 registered | Claude+user |
| 2026-05-27 | Fig S1  | LOCKED — label overlap fixed, x-axis label improved, caption updated, rebuilt and MD5 registered | Claude+user |
| 2026-05-27 | Fig S2  | LOCKED — candidate colocalization stacked bar plot generated, caption refined and MD5 registered | Claude+user |
| 2026-05-27 | Title/Abs| LOCKED — Title + Abstract v2 locked (TED aligned, structured 227 words, JEI format) | Claude+user |
| 2026-05-27 | Results | LOCKED — Results 6 sections (R1–R6) FINAL prose locked and integrated | Claude+user |
| 2026-05-27 | Methods P1| LOCKED — Methods Part 1 (Data sources & instrument definition) drafted and locked | Claude+user |
| 2026-05-27 | Methods P2| LOCKED — Methods Part 2 (Outcome cohorts & data harmonization) drafted and locked | Claude+user |
| 2026-05-27 | Methods P3| LOCKED — Methods Part 3 (MR estimators & sensitivity analyses) drafted and locked | Claude+user |
| 2026-05-27 | Methods P4| LOCKED — Methods Part 4 (Fine-mapping & colocalization) drafted and locked | Claude+user |
| 2026-05-27 | Methods P5| LOCKED — Methods Part 5 (In-house RNA-seq) drafted and locked | Claude+user |
| 2026-05-27 | Methods   | LOCKED — Methods Section (M1–M8 / Parts 1–6) FINAL prose locked and integrated | Claude+user |
| 2026-05-27 | Declarations| LOCKED — Declarations (Acknowledgments, Funding, Contributions, Ethics) drafted and integrated | Claude+user |
| 2026-05-27 | Cover Letter| LOCKED — Cover Letter for JEI submission drafted and compiled | Claude+user |
| 2026-05-27 | docx Files| LOCKED — MANUSCRIPT_TED_TRAP_SUBMISSION.docx and COVER_LETTER.docx compiled and locked | Claude+user |

---

*Registry maintained as assets are verified. Manuscript draws only from LOCKED assets.*
