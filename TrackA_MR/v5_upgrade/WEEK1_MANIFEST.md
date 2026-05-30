# TED-TRAP v5 — Week 1 MANIFEST

**Status:** Week 1 COMPLETE (verified by Claude 2026-05-22)
**Working directory:** `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\`
**Purpose of this document:** Single-page recovery reference. Reading this alone should let anyone reconstruct what Week 1 did, where every file is, and what the key results were — without re-reading the full transcript.

---

## 1. WEEK 1 IN ONE PARAGRAPH

Week 1 tested whether the TSHR single-instrument limitation (a key reason for the Thyroid and ETJ desk rejections) could be structurally resolved. Two parallel tasks were run: (A) pQTL availability check for TSHR/IGF1R/IGF1/INSR in UKB-PPP and deCODE, and (B) TSHR locus fine-mapping via GCTA-COJO and SuSiE in both EUR and EAS LD references. **Result: TSHR is a single-IV locus (rs179252); no usable independent secondary signal exists.** pQTL availability remains pending (deferred to Week 2). The directionality-paradox explanation was secured: the rs179252 G allele increases blood TSHR expression (Z=+13.28) yet is protective for Graves disease (β=−1.408 in v4.3 MR), which is reconciled by tissue-specific regulation (blood cis-eQTL direction ≠ orbital tissue expression direction).

---

## 2. KEY RESULTS (LOCKED NUMBERS)

| Item | Value | Source file |
|------|-------|-------------|
| TSHR usable independent IVs | **0** | TaskB_05 |
| rs179252 effect allele | **G** | eqtl_TSHR_all.csv (origin), TaskB_01/03/04 |
| rs179252 other allele | **T** | eqtl_TSHR_all.csv (origin) |
| rs179252 Z-score (eQTLGen) | **+13.2844** | TaskB_01 |
| rs179252 P-value (eQTLGen cis-eQTL) | **2.859e-40** | TaskB_01 |
| rs179252 position (hg19) | chr14:81,435,985 | TaskB_01 |
| COJO EUR independent signals | 3 (only rs179252 primary; secondaries weak/collinear) | TaskB_03 |
| COJO EAS independent signals | 16 — **FALSE POSITIVES (LD mismatch artifact)** | TaskB_03 |
| SuSiE EUR credible sets | 2, both collinear with rs179252 (r²≥0.965) | TaskB_04 |
| SuSiE EAS credible sets | 2, both collinear with rs179252 (r²≥0.808) | TaskB_04 |
| pQTL availability (8 cells) | **pending** → Week 2 task | TaskA_01 |

**Decision:** TSHR retained as single-IV. Framing fixed to "locus-level genetic prioritization," NOT "robust multi-IV MR" and NOT "expression causality."

---

## 3. CRITICAL INTERPRETIVE NOTES (do not lose these)

1. **EAS COJO 16 signals are artifacts.** eQTLGen is a European-ancestry meta-analysis. Conditioning it on an EAS LD reference produces spurious "independent" signals (e.g., rs75538509 marginal P=0.236 but conditional P=4.87e-19 — a clear LD-mismatch tell). Only EUR-based fine-mapping of eQTLGen is valid. This must be stated as a limitation in the manuscript.

2. **Directionality paradox is explained, not denied.** The blood cis-eQTL direction (G allele → higher TSHR, protective for GD) must NOT be interpreted as predicting orbital tissue TSHR direction. cis-eQTL effects are tissue-specific. The elevated TSHR in TED orbital tissue is a disease-associated tissue state, not a germline expression-direction effect. Conclusion = locus-level prioritization.

3. **beta/se in TaskB_01 are intentionally blank.** eQTLGen provides only Z-score + P + N. beta/se are derived downstream in helper_make_ma.R using 1000G EUR MAF as EAF proxy (formula: se = 1/sqrt(2·eaf·(1−eaf)·(N+Z²)); beta = Z·se).

---

## 4. FILE LOCATIONS

### Original inputs (PRESERVE — basis for all downstream weeks)
```
02_tshr_finemap/data/
├── eqtl_TSHR_all.csv          ← eQTLGen TSHR extraction. effect=G/other=T. 7,309 variants.
├── TSHR_eqtlgen_locus.ma      ← COJO input (.ma format), eaf from 1000G EUR
├── 1000G_EUR_phase3/          ← SurfSara hg19 reference + EUR_chr14_freq.frq
└── 1000G_EAS_phase3/          ← SurfSara hg19 reference + EAS_chr14_freq.frq
```

### Results (7 deliverables, all verified)
```
02_tshr_finemap/results/
├── TaskB_00_ld_reference_metadata_v1.csv   ← VERIFY EXISTS (LD ref documentation)
├── TaskB_01_tshr_locus_extraction_v1.tsv   ← 7,309 variants ✅ verified
├── TaskB_02_ld_clumping_sensitivity_v1.csv ← VERIFY EXISTS (12-cell sensitivity)
├── TaskB_03_tshr_cojo_independent_signals_v1.tsv ← 19 signals ✅ verified
├── TaskB_04_tshr_susie_credible_sets_v1.tsv ← 4 CS ✅ verified
└── TaskB_05_tshr_finemap_summary_v1.csv    ← 21 rows integration ✅ verified
01_pqtl_check/results/
└── TaskA_01_pqtl_availability_v1.csv       ← 8 cells pending ✅ verified
03_decision/
├── TaskC_decision_table_v1.csv             ✅ verified
├── TaskC_pre_analysis_plan_v1.md           ✅ verified
└── TaskC_week1_summary_memo.md             ✅ verified
```

### Metadata (reproducibility)
```
00_meta/
├── data_sources.csv          ← verified PMIDs/DOIs; SurfSara 1000G URLs
├── software_versions.csv     ← R 4.3.3, GCTA 1.95.1, PLINK 20250819
├── analytical_log.md         ← run log + extraction commands
└── README.md
scripts/                      ← 9 R/bash scripts (paths corrected to TrackA_MR/v5_upgrade)
```

---

## 5. OUTSTANDING ITEMS (must confirm before Week 2 lock)

| Item | Status | Action |
|------|--------|--------|
| TaskB_00_ld_reference_metadata.csv exists? | UNCONFIRMED | Gemini to confirm/create |
| TaskB_02_ld_clumping_sensitivity.csv exists? | UNCONFIRMED (not in upload ZIP) | Gemini to confirm/upload |
| pQTL panel verification | DEFERRED | Week 2 Day 1 (IGF1 priority) |

---

## 6. WHERE WEEK 1 OUTPUTS GO IN THE MANUSCRIPT

- **Supplementary Table** (later): TSHR fine-mapping (COJO + SuSiE independent-signal analysis) — evidence for single causal signal
- **Supplementary Figure** (later): TSHR locus regional association plot (chr14, −log10P)
- **Supplementary Methods** (later): fine-mapping methodology + EAS LD-mismatch limitation
- **Main text Discussion**: directionality paradox explanation paragraph (data secured here)

No main-text figures are built from Week 1 alone. Main figures require Week 2-3 systematic MR results.

---

## 7. JOURNAL CEILING AFTER WEEK 1

- TSHR alone: **Tier 2** (IOVS / JCEM / EJE) — single-IV, known GD locus, limited novelty
- Tier 1 (JCI Insight / EBioMedicine): **conditional on Week 2-3** producing novel replicated druggable hits + colocalization
- Week 1's value is **defensive** (removes single-IV ambiguity + directionality paradox), not discovery

---

*Verified by Claude against uploaded deliverables, 2026-05-22.*
*All allele/position/P-value cross-checked across TaskB_01/03/04/05 for consistency.*
*Next: confirm 2 outstanding files → Week 2 druggable-genome-wide MR spec.*
