# TED-TRAP v5 — Week 2 MANIFEST

**Status:** Week 2 COMPLETE (verified by Claude 2026-05-26)
**Working directory:** `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\`
**Purpose:** Single-page recovery reference for the druggable-genome-wide MR.

---

## 1. WEEK 2 IN ONE PARAGRAPH

Week 2 expanded TED-TRAP from a ~20-gene candidate study into a systematic druggable-genome-wide cis-eQTL MR. The Finan 2017 druggable genome (4,462 genes mapped to GRCh37) was screened against eQTLGen blood cis-eQTLs; 2,545 genes had a valid instrument (Bonferroni denominator; threshold 1.965e-5). MR was run against three local, MD5-verified outcomes: BBJ Graves disease (discovery, EAS), UKB hyperthyroidism (replication, EUR; LMM beta rescaled to log-odds), and FinnGen R12 Graves ophthalmopathy (TED sensitivity). Thirteen genes reached BBJ Bonferroni significance. After defensive filters (MHC spillover, chr16p11.2 LD cluster, single-IV unverifiability, FinnGen direction), **no defensible robust-novel druggable target remained.** The two genes passing all criteria — TSHR and CTLA4 — are already-known GD loci. **Decision: Option 1 defensive paper** anchored on TSHR (protective, upstream) vs IGF1R (risk, effector), with CTLA4 as positive control and the systematic screen resolving the narrow-scope criticism.

---

## 2. KEY NUMBERS (LOCKED)

| Item | Value |
|------|-------|
| Druggable genes (Finan 2017, GRCh37-mapped) | 4,462 |
| Genes with valid cis-eQTL instrument | **2,545** (Bonferroni denominator) |
| Single-IV / Multi-IV | 990 / 1,555 |
| Bonferroni threshold | **1.965e-5** (0.05/2,545) |
| MR result rows (3 outcomes) | 13,042 |
| Harmonization drops | 23 genes (0.9%) |
| BBJ Bonferroni hits | 13 |
| Robust-novel (defensible) | **0** |

---

## 3. THE THREE DEFENSIBLE GENES (paper backbone)

| Gene | Direction | BBJ | UKB | FinnGen | n_iv | Role |
|------|-----------|-----|-----|---------|------|------|
| **TSHR** | protective | −2.10 (1.1e-14) | −2.44 (8.8e-28) | −2.33 (2.8e-7) | 1 | upstream anchor |
| **IGF1R** | risk | +0.45 (0.021) | +0.30 (0.012) | +0.34 (0.18) | 4 | effector (teprotumumab) |
| **CTLA4** | protective | −1.74 (5.5e-15) | −1.57 (0.002) | −1.77 (0.012) | 1 | positive control |

Central narrative: TSHR (protective, upstream) vs IGF1R (risk, effector) = opposite-direction axes within a systematic druggable screen.

---

## 4. WHY NOVEL = 0 (do not lose this reasoning)

The 13 hits broke down as: 6 MHC spillover (chr6), 3 chr16p11.2 LD cluster (HSD3B7/VKORC1/PRSS36, <800 kb apart — likely one shared signal), 2 single-IV unverifiable (TNFSF14/IFNGR1), 1 discovery-only (MAPKAPK5, UKB replication fail), 2 known anchors (TSHR/CTLA4). Gemini's preliminary "3 novel" (TNFSF14/HSD3B7/PRSS36) all failed: every one is single-IV (no pleiotropy test) and TNFSF14/HSD3B7 reverse or null direction in FinnGen GO. chr16p11.2 trio is an LD-cluster artifact analogous to MHC.

---

## 5. FILE LOCATIONS

```
04_druggable_mr/
├── data/
│   ├── finan2017_druggable.csv         (4,479 input gene list)
│   ├── g1000_eur_freq.frq              (EUR MAF for beta/se derivation)
│   └── outcomes/
│       ├── GCST90018627_harmonised.tsv.gz  (BBJ Graves, MD5 325a1e...)
│       └── GCST90038636_harmonised.tsv.gz  (UKB hyperthyroid, MD5 afa33e...)
├── results/
│   ├── TaskD_01_druggable_gene_master_v1.csv      (4,462 genes + GRCh37 coords)
│   ├── TaskD_02a_eqtl_instruments_snp_level_v1.csv (6,136 instrument SNPs) ✅verified
│   ├── TaskD_02b_eqtl_instrument_gene_summary_v1.csv (4,462 gene rows) ✅verified
│   ├── TaskD_02_denominator_note.txt              (denominator=2,545)
│   ├── TaskD_03a_MR_BBJ_discovery_v1.csv          (3,348 rows)
│   ├── TaskD_03b_MR_UKB_replication_v1.csv        (4,834 rows)
│   ├── TaskD_03c_MR_FinnGen_GO_sensitivity_v1.csv (4,860 rows)
│   ├── TaskD_03d_MR_all_outcomes_merged_v1.csv    (13,042 rows) ✅verified
│   └── TaskD_03e_harmonization_qc_v1.csv          (QC)
05_coloc_candidates/  (Week 3 target)
WEEK2_D4_VERIFICATION.md                 (Claude-verified hit classification)
TaskD_04_hit_classification_VERIFIED.csv (13 hits, full flags)
scripts/ taskD_01 / taskD_02 / taskD_03 / taskD_04
```

FinnGen GO outcome: `c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz` (local, pre-existing).

---

## 6. CRITICAL CAVEATS FOR MANUSCRIPT

1. **NOT the first druggable/immune-gene GD MR.** Prior work: medRxiv 2025.01.02.25319932 (pQTL→GD); Sci Reports s41598-025-21754-4 (multi-omics MR, FGFRL1). Frame as contextualizing TSHR vs IGF1R, not "first."
2. **Cross-ancestry:** eQTLGen EUR instruments vs BBJ EAS discovery — document.
3. **Single-IV genes (incl. TSHR):** Wald only, no pleiotropy diagnostics.
4. **MHC + chr16p11.2 excluded from novel claims** (LD spillover).
5. **UKB beta rescaled** LMM→log-odds (β/se ÷ μ(1−μ), μ=3731/484598); P/direction preserved; OR now comparable to BBJ.

---

## 7. JOURNAL POSITION AFTER WEEK 2

No novel target → not a "new discovery" paper. Realistic tier: **IOVS / JCEM / EJE / Front Endocrinol.** Strength = clean TSHR-vs-IGF1R contrast + systematic scope + positive control + multi-modal (Week 1 fine-mapping + tissue RNA-seq still to integrate). Honest framing is the reviewer-proofing.

---

## 8. WEEK 3 PLAN (coloc)

Colocalization for: TSHR (confirm; v4.3 BBJ PP.H4=0.951), IGF1R (strengthen effector multi-IV), CTLA4 (positive control). Optionally chr16p11.2 to formally demonstrate LD-spillover. No novel-hit coloc needed (none claimed).

---

*Verified by Claude against uploaded D.3 merged + D.2 instrument files, 2026-05-26.*
*Next: Week 2 archive lock → Week 3 colocalization.*
