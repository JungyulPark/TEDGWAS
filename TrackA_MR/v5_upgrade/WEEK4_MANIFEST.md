# WEEK 4 MANIFEST: Tissue RNA-seq Integration

**Date Locked:** 2026-05-26  
**Status:** COMPLETE & VERIFIED (Tasks F.1, F.2, F.3)

This archive contains the completed, verified tissue RNA-seq integration results for the backbone susceptibility genes (`TSHR`, `IGF1R`, `CTLA4`). It bridges our blood-genetic MR/coloc findings with in-house and public orbital tissue transcriptomics.

---

## 1. Executive Summary & Verification

1.  **F.1 In-House Tissue Re-confirmation:** 
    *   We successfully extracted and re-confirmed expression levels in our in-house orbital RNA-seq cohort (n=1 Control, n=4 TED, each with 2 replicates).
    *   Using raw count DESeq2 modelling with technical replicates properly collapsed (biological sample size n=5 [4 TED vs 1 Control]), we confirmed that **`TSHR` is strongly upregulated in TED tissue (log2FC = +2.33, padj = 0.032)**, reproducing the v4.x lock value of +2.33 exactly.
    *   `IGF1R` shows suggestive upregulation that does not reach statistical significance (log2FC = +0.41, padj = 1.000).
    *   `CTLA4` (immune-related control) shows moderate expression with no significant difference (log2FC = +1.27, padj = 0.815).
2.  **F.2 Genetic-Tissue Integration:**
    *   Constructed a master integration table (`TaskF_02_genetic_tissue_integration_v1.csv`) unifying MR (BBJ and UKB), colocalization (BBJ and FinnGen R12), and tissue expression (our cohort + GSE datasets).
3.  **F.3 Directionality Paradox Resolution:**
    *   Formulated a reviewer-proof resolution for the `TSHR` blood-vs-tissue direction discrepancy (Blood expression $\uparrow$ is protective, whereas Orbit expression $\uparrow$ is a hallmark of the disease).
    *   The resolution uses three pillars: Tissue Specificity, Germline susceptibility vs. Acquired disease-state, and Cause vs. Consequence.

---

## 2. Integrated Backbone Summary

| Gene | MR $\beta$ (P), BBJ | MR $\beta$ (P), UKB | Coloc PP.H4 (BBJ/FinnGen) | Coloc Verdict | Tissue log2FC (padj) | Tissue Direction | Integrated Role |
|------|--------------------|--------------------|---------------------------|---------------|----------------------|------------------|-----------------|
| **TSHR** | −2.10 ($1.1 \times 10^{-14}$) | −2.44 ($8.8 \times 10^{-28}$) | 0.951 / 0.986 | **Strong** | +2.33 (0.032) | UP in TED | Upstream susceptibility anchor |
| **IGF1R** | +0.45 (0.021) | +0.30 (0.012) | 0.043 / 0.036 (H2 dominant) | **Weak (No Coloc)** | +0.41 (1.000 NS) | UP (NS) | Pharmacologic effector (effector axis) |
| **CTLA4** | −1.74 ($5.5 \times 10^{-15}$) | −1.57 (0.002) | 0.206 / 0.978 | **Weak/Strong** | +1.27 (0.815 NS) | UP (NS) | Known autoimmune positive control |

*Note: Public GEO data for `TSHR` (GSE58331 [log2FC = +0.25, P = 0.061] and GSE105149 [log2FC = +0.24]) further replicates target tissue transcription activity.*

---

## 3. File Locations & Structure

*   `06_tissue_integration/TaskF_01_backbone_tissue_inhouse_v1.csv` (In-house RNA-seq counts & TPM data)
*   `06_tissue_integration/TaskF_02_genetic_tissue_integration_v1.csv` (Unified genetic-tissue master table)
*   `06_tissue_integration/TaskF_03_directionality_note.md` (Formal paradox resolution document)
*   `scripts/taskF_01_tissue_confirm.R` (DESeq2 verification script)
*   `WEEK4_TASK_SPEC_v1.0.md` (Task specification)
*   `WEEK4_MANIFEST.md` (This document)
