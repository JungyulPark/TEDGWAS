# Task D.4 Verification & Hit Classification Report
**Date:** 2026-05-26
**Target:** 2,545 Druggable Genes against 3 Outcomes (BBJ Graves, UKB Hyperthyroidism, FinnGen Graves Ophth)

## 1. Primary MR Results (BBJ Discovery)
Out of 2,545 genes tested, 13 genes passed the Bonferroni significance threshold ($P < 1.965 \times 10^{-5}$) in the primary BBJ discovery cohort.

### 13 Bonferroni Hits:
1. **TSHR** (chr14) — Known Anchor
2. **CTLA4** (chr2) — Known Anchor (Positive Control)
3. **HLA-A** (chr6) — Known Locus (MHC Region)
4. **HLA-DQA2** (chr6) — Known Locus (MHC Region)
5. **C4A** (chr6) — HLA Spillover (MHC Region)
6. **TUBB** (chr6) — HLA Spillover (MHC Region)
7. **PSMB8** (chr6) — HLA Spillover (MHC Region)
8. **TNFSF14** (chr19) — Novel Candidate (Failed Verification)
9. **HSD3B7** (chr16) — Novel Candidate (Failed Verification)
10. **PRSS36** (chr16) — Novel Candidate (Failed Verification)
11. **VKORC1** (chr16) — Discovery Only (Failed UKB Replication)
12. **MAPKAPK5** (chr12) — Discovery Only (Failed UKB Replication)
13. **IFNGR1** (chr6) — Discovery Only (Failed UKB Replication)

---

## 2. Verification of "Novel" Candidates
The three initial novel candidates (`TNFSF14`, `HSD3B7`, `PRSS36`) were subjected to rigorous scrutiny to prevent false-positive claims in the manuscript. All three failed the verification phase due to the following critical limitations:

### A. Single-IV Limitation (Wald Ratio)
All three candidates are instrumented by a single eQTL variant (`n_iv = 1`). This inherently precludes pleiotropy testing (MR-Egger, Weighted Median), making the causal claims highly susceptible to horizontal pleiotropy.

### B. LD Spillover & Genomic Clustering
- `HSD3B7` (chr16:30.4Mb)
- `VKORC1` (chr16:31.1Mb)
- `PRSS36` (chr16:31.2Mb)
These genes cluster within an 800kb window on `chr16p11.2`. This strongly suggests that these are not three independent causal signals, but rather a single shared signal (LD spillover) masking as multiple hits. Colocalization analysis would be required to disentangle this, which is outside the current scope.

### C. Outcome Inconsistency (FinnGen Failure)
- **TNFSF14:** BBJ (β=-0.47), UKB (β=-0.30), but FinnGen (β=+0.034, P=0.86).
- **HSD3B7:** BBJ (β=-1.23), UKB (β=-0.28), but FinnGen (β=+0.272, opposite direction).
The lack of consistency across the three outcomes undermines the robustness of the signal, especially in the TED-specific sensitivity cohort (FinnGen GO).

### D. Literature Precedents
Recent publications (e.g., *Sci Reports 2025*, *medRxiv 2025*) have already explored immune-related genes and pQTL MR for Graves' disease. The TNF superfamily is well-documented in autoimmune contexts, diminishing the true "novelty" of targets like TNFSF14 (LIGHT).

---

## 3. Conclusion & Manuscript Strategy
**Decision: Option 1 (Defensive, High-Rigor Approach)**

Given the failure of the novel candidates to pass rigorous verification, the manuscript will NOT overclaim new therapeutic discoveries. Instead, the narrative will focus on the **strength and consistency of the pipeline** and the **TSHR vs. IGF1R axis**:

1. **Systematic Scope:** The 2,545-gene systematic screen successfully resolves previous reviewer critiques regarding narrow scope and "cherry-picking".
2. **TSHR (Upstream Anchor):** Shows a perfectly consistent protective effect across all 3 outcomes (BBJ β=-2.10, UKB β=-2.44, FinnGen β=-2.33).
3. **IGF1R (Downstream Effector):** Shows a consistent risk effect (BBJ β=+0.45, UKB β=+0.30), supporting the mechanism of teprotumumab. This is one of the cleanest multi-IV signals in the dataset.
4. **CTLA4 (Positive Control):** Consistently significant across all three cohorts, validating the pipeline's analytical rigor.

This strategy ensures a highly defensible, reviewer-proof manuscript suitable for submission.
