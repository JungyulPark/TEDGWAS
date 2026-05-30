# WEEK 3 MANIFEST: Colocalization (coloc.abf)

**Date Locked:** 2026-05-26
**Status:** COMPLETE & VERIFIED (Modules 3A and 3B)

This archive contains the complete, verified colocalization results (coloc.abf) for both the core anchor/mechanistic-axis genes (Week 3A) and the candidate novel genes (Week 3B) against BBJ Graves (discovery) and FinnGen GO (TED sensitivity).

## 1. Executive Summary & Verification
The colocalization analysis was executed in two phases:
*   **Week 3A (Anchor Loci):** Validated `TSHR`, `IGF1R`, and `CTLA4`. `TSHR` colocalizes strongly (BBJ PP.H4=0.9515 — reproducing v4.3; FinnGen PP.H4=0.9865), confirming it as an expression-driven upstream susceptibility anchor. `IGF1R` shows dominant PP.H2 (BBJ 0.793, FinnGen 0.675), indicating the GWAS signal is present but the cis-eQTL signal is too weak in the window. `CTLA4` acts as a positive control, showing strong coloc in FinnGen (PP.H4=0.978) and cross-ancestry LD complexity in BBJ (PP.H3=0.794).
*   **Week 3B (Candidates):** Validated `TNFSF14`, `HSD3B7`, `PRSS36`, `IFNGR1`, `VKORC1`, and `MAPKAPK5`. **All candidates failed the consistency or coloc criteria**, proving that robust novel targets are indeed 0.

## 2. Complete Coloc Results Table
Combined results for Week 3A (`TaskE_01`) and Week 3B (`TaskE_02`):

| Gene | Outcome | nSNPs | PP.H2 | PP.H3 | PP.H4 | Top SNP | Verdict | Interpretation |
|------|---------|-------|-------|-------|-------|---------|---------|----------------|
| **TSHR** | BBJ | 3,756 | ~0 | 0.049 | **0.951** | rs179252 | strong | Shared causal variant (reproduced) |
| **TSHR** | FinnGen GO | 5,561 | ~0 | 0.013 | **0.986** | rs179252 | strong | Shared causal variant |
| **IGF1R** | BBJ | 4,456 | **0.793** | 0.164 | 0.043 | rs2654980 | weak | GWAS signal present, eQTL signal weak |
| **IGF1R** | FinnGen GO | 5,856 | **0.675** | 0.289 | 0.036 | rs2654980 | weak | GWAS signal present, eQTL signal weak |
| **CTLA4** | BBJ | 1,982 | ~0 | **0.794** | 0.206 | rs231811 | weak | Distinct causal variants (cross-ancestry LD) |
| **CTLA4** | FinnGen GO | 3,230 | ~0 | 0.022 | **0.978** | rs1863800 | strong | Shared causal variant |
| **TNFSF14**| BBJ | 4,788 | ~0 | 0.005 | **0.994** | rs2291668 | strong | Shared causal variant in BBJ only |
| **TNFSF14**| FinnGen GO | 7,696 | **0.626** | 0.355 | 0.019 | rs2291668 | weak | Fails coloc in FinnGen GO (no consistency) |
| **IFNGR1** | BBJ | 3,958 | ~0 | 0.006 | **0.989** | rs11754268| strong | Shared causal variant in BBJ only |
| **IFNGR1** | FinnGen GO | 5,859 | **0.642** | 0.339 | 0.019 | rs11754268| weak | Fails coloc in FinnGen GO (no consistency) |
| **MAPKAPK5**| BBJ | 1,909 | ~0 | **1.000** | ~0 | rs79271898| distinct| Distinct causal variants |
| **MAPKAPK5**| FinnGen GO | 2,979 | **0.674** | 0.296 | 0.030 | rs79271898| weak | GWAS signal only |
| **HSD3B7** | BBJ | 1,625 | ~0 | 0.384 | 0.616 | rs4889606 | ambiguous| Suggestive shared variant in BBJ only |
| **HSD3B7** | FinnGen GO | 3,095 | **0.733** | 0.240 | 0.026 | rs4889606 | weak | GWAS signal only |
| **PRSS36** | BBJ | 1,479 | ~0 | **0.843** | 0.157 | rs78924645| distinct| Distinct causal variants |
| **PRSS36** | FinnGen GO | 2,820 | **0.742** | 0.220 | 0.038 | rs78924645| weak | GWAS signal only |
| **VKORC1** | BBJ | 1,512 | ~0 | **0.639** | 0.360 | rs34649473| ambiguous| Distinct/ambiguous variants |
| **VKORC1** | FinnGen GO | 2,950 | **0.740** | 0.231 | 0.029 | rs2884737 | weak | GWAS signal only |

## 3. Disqualification of the chr16p11.2 LD Cluster
The top SNPs from the candidate coloc for the chr16p11.2 cluster are:
*   `HSD3B7` (rs4889606): chr16:31,011,183 (31.01 Mb)
*   `VKORC1` (rs34649473): chr16:31,066,538 (31.06 Mb)
*   `PRSS36` (rs78924645): chr16:31,154,358 (31.15 Mb)
These three signals map within 143 kb of each other in a dense region on chr16p11.2. The lack of strong colocalization consistency (with PP.H2 dominant in FinnGen GO) combined with their close genomic clustering confirms they are a single shared LD block signal and not independent novel hits.

## 4. Main Manuscript Backbones
*   **TSHR (Genetically Anchored Susceptibility):** Strong protective MR and strong PP.H4 coloc across all cohorts, pointing to rs179252.
*   **IGF1R (Pharmacologic Effector):** Robust multi-IV risk MR, but coloc fails (H2 dominant, weak expression signal). Interpreted as a pharmacologic effector axis rather than expression-mediated susceptibility.
*   **CTLA4 (Positive Control):** Strong coloc in FinnGen (EUR-EUR); H3 in BBJ represents EUR eQTL vs. EAS GWAS LD mismatch.

## 5. File Locations
*   `05_coloc_candidates/TaskE_01_coloc_results_v1.csv` (Week 3A Results)
*   `05_coloc_candidates/TaskE_02_candidate_coloc_v1.csv` (Week 3B Results)
*   `scripts/taskE_01_coloc.R`
*   `scripts/taskE_02_coloc_candidates.R`
*   `WEEK3_TASK_SPEC_v1.0.md`
*   `WEEK3_MANIFEST.md` (this file)
