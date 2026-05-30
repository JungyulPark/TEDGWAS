# TED-TRAP v5: Pre-Analysis Plan

## 1. Background and Rationale
The overarching goal of the TED-TRAP v5 upgrade is to execute a druggable-gene-wide genetic prioritization of Graves' Disease (GD) and Thyroid Eye Disease (TED) susceptibility architecture. Based on Week 1 fine-mapping and LD conditional analyses, the TSHR locus is determined to contain a single causal signal represented by rs179252.

## 2. Exposures
* **Primary Exposure:** TSHR expression (eQTLGen, cis-eQTL).
* **Instrument Variable (IV):** Single-IV model utilizing `rs179252` as the instrument. All other secondary signals detected via COJO or SuSiE were found to be either false positives due to LD reference mismatch or highly collinear proxies of rs179252 (r2 > 0.8). 
* **pQTL Exposure:** Currently marked as PENDING. pQTL availability is explicitly scheduled as the first task of Week 2. At the start of Week 2, the complete protein panel lists from UKB-PPP (Sun 2023) and deCODE (Ferkingstad 2021) will be downloaded to definitively resolve all 8 target cells. In particular, because IGF1 is a secreted protein, it has a high probability of being present on the plasma panels; verifying its pQTL MR feasibility is a mandatory next step to establish orthogonal evidence for the IGF axis.

## 3. Outcomes
The MR analyses will proceed according to the established outcome hierarchy:
1. **Primary Discovery:** BBJ Graves disease (EAS, 2,809 cases)
2. **Replication:** UK Biobank hyperthyroidism (EUR, 3,731 cases)
3. **TED-specific Sensitivity:** FinnGen R12 Graves ophthalmopathy (858 cases)

## 4. Analysis Strategy
* **Methodology:** Since TSHR yields a single independent cis-eQTL instrument, all causal estimates for TSHR will be derived using the **Wald Ratio** method. Multi-IV MR techniques (e.g., Inverse-Variance Weighted, MR-Egger) are structurally not applicable.
* **Sensitivity Analyses:** Steiger filtering to assess reverse causality. Colocalization will be considered if sufficient summary statistics are available.

## 5. Downstream Expansion (Week 2-3)
Given the single-IV constraint at the TSHR locus, the study will be explicitly framed as a **locus-level prioritization** for TSHR. To achieve the required methodological strength for high-impact journals (e.g., JCI Insight, EBioMedicine), the study must proceed to **Week 2-3 Druggable-Genome-Wide MR**. This expansion aims to discover novel, replicated druggable hits that extend beyond the TSHR locus to establish robust multi-IV causal pathways in TED pathogenesis.
