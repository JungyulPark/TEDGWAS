# TED-TRAP Manuscript — Major Upgrade Execution Plan

**Version**: v3.0 (Post-Review Revision)
**Target journal**: Thyroid (IF 5.4) — European Thyroid Journal as backup
**Current assessment**: Impact 5.5/10 → Target 7/10 after upgrade
**Execution environment**: VS Code + Antigravity, Windows
**Working directory**: `c:\ProjectTEDGWAS\` (as seen in screenshot)

---

## 🎯 UPGRADE OBJECTIVES

Based on the independent expert methodological review (6 CRITICAL + 9 MAJOR issues identified),
we are rebuilding the analytical pipeline with corrected statistics. This document
organizes the upgrade into **6 sequential phases**, each with explicit deliverables.

The prior analyses are preserved as **"v1 archive"**; all new analyses live in a **"v3" namespace**
to avoid file collision.

---

## 📋 PHASE STRUCTURE

| Phase | Title | Time | Deliverables |
|-------|-------|------|--------------|
| **P1** | MR re-analysis (scale fix + full sensitivity + all-gene coloc + GTEx parallel) | 5-7 days | new MR results Excel sheet |
| **P2** | RNA-seq re-analysis (DESeq2 + BH-FDR) | 2-3 days | DESeq2 output, new TPM table |
| **P3** | TED gene set cross-validation (Open Targets + DisGeNET API + PubTator + MAGMA) | 3-4 days | Concordance report, final curated list |
| **P4** | Pathway analysis upgrade (hypergeometric + network proximity) | 2-3 days | Statistical enrichment tables |
| **P5** | Table/Figure regeneration (v3 numbers) | 1 week | 6 Main/Supp tables + 3 main figures |
| **P6** | Manuscript drafting with corrected results | 2-3 weeks | Full manuscript + cover letter |

**Total timeline**: 4-6 weeks

---

## 📝 PHASE 1: MR RE-ANALYSIS

### Goals
1. Fix β-scale mismatch (primary vs. replication outcome)
2. Re-run MR with full sensitivity suite for ALL genes
3. Add colocalization for every candidate locus (not just TSHR)
4. Add GTEx v8 thyroid/adipose parallel MR
5. MR-PRESSO for heterogeneous loci (TNF)
6. MRlap for eQTLGen↔UKB sample overlap
7. Add BioBank Japan replication
8. Multivariable MR (TED | GD) to isolate TED-specific effects

### Critical fixes
**P1.1** β-scale: Primary log-odds vs. Replication linear — use β_logOR ≈ β_linear / [p(1−p)] with p=0.0077
  → OR switch replication to FinnGen R12 hyperthyroidism (log-odds SAIGE)

**P1.2** TSHR 2-IV: explicitly reframe as *cis*-MR / drug-target MR
  → Add coloc-SuSiE as formal consistency test
  → Report per-SNP Wald ratio as de facto leave-one-out
  → Cite Schmidt 2020, Burgess 2023, Zuber 2022

**P1.3** TNF Cochran Q = 0.028 → MR-PRESSO outlier correction

### Script deliverables (see /scripts_R)
- `01_token_refresh_and_setup.R` — OpenGWAS JWT renewal
- `02_eqtlgen_exposure_extraction_all_genes.R`
- `03_mr_all_genes_full_sensitivity.R` — IVW + Egger + WM + WMode + MR-PRESSO + CAUSE
- `04_mr_scale_rescale.R` — Convert linear betas to log-odds
- `05_gtex_thyroid_adipose_parallel.R`
- `06_coloc_all_loci.R` — coloc.abf + coloc.susie for every gene
- `07_mrlap_overlap_correction.R`
- `08_mvmr_ted_given_gd.R`
- `09_bbj_replication.R`

### Expected outputs
- `outputs/mr_v3_full_results.csv` (all genes × all methods × all outcomes)
- `outputs/coloc_v3_all_loci.csv`
- `outputs/mr_v3_mrlap_corrected.csv`
- `outputs/mvmr_v3.csv`

---

## 📝 PHASE 2: RNA-SEQ RE-ANALYSIS

### Goals
1. Replace Wilcoxon with DESeq2
2. Apply BH-FDR correction
3. Restrict testing to pre-specified MR-triangulated candidate set (~8 genes)
4. Direction-of-effect concordance focus (not p-value only)
5. Add public GEO dataset meta-validation (optional)

### Required inputs
- Raw count matrix from STAR/RSEM/salmon output
- OR: TPM + gene-level counts from prior analysis (user provides)

### Critical fix
**P2.1** Current Wilcoxon with n=4 vs n=1 is mathematically incoherent
  → Exact null for (4,1) has only 5 permutations; minimum two-sided p ≈ 0.4
  → Reported p=0.019 etc. cannot arise from exact Wilcoxon

**P2.2** Frame as "orthogonal tissue-localization confirmation" not genome-wide discovery

**P2.3** Verify INSR novelty against Kim 2021 IOVS and Kim 2024 JCI Insight

### Script deliverables (see /scripts_R)
- `10_deseq2_analysis.R` — Empirical Bayes DE with n=4 TED vs n=1 Ctrl
- `11_bh_fdr_correction.R`
- `12_gene_set_testing_candidate.R` — CAMERA / ROAST with MR-triangulated gene set
- `13_insr_novelty_check.R` — PubMed/GEO meta-analysis against Kim 2021/2024

### Expected outputs
- `outputs/deseq2_v3_results.csv`
- `outputs/candidate_gene_test_v3.csv`
- `outputs/insr_novelty_verification.md`

---

## 📝 PHASE 3: TED GENE SET CROSS-VALIDATION

### Goals
1. Replace manual 59-gene curation with reproducible multi-source intersection
2. Cross-validate with Open Targets + DisGeNET + PubTator + MAGMA

### Required sources
- **Open Targets Platform**: GraphQL API query for "thyroid-associated ophthalmopathy" EFO
- **DisGeNET v7+**: REST API with UMLS CUI C0018818 (Graves ophthalmopathy)
- **PubTator Central**: disease-gene co-occurrence from TED-indexed abstracts
- **MAGMA**: gene-set analysis on FinnGen R12 GO summary statistics

### Script deliverables (see /scripts_python)
- `14_opentargets_query.py`
- `15_disgenet_rest_api.py`
- `16_pubtator_mining.py`
- `17_magma_finngen_go.py` (R script for MAGMA integration)
- `18_ted_gene_set_concordance.py` — harmonize 4 sources

### Expected outputs
- `outputs/ted_genes_v3_final.csv` (multi-source validated)
- `outputs/ted_gene_concordance_report.md`

---

## 📝 PHASE 4: PATHWAY/NETWORK UPGRADE

### Goals
1. Replace Venn-only description with hypergeometric enrichment p-values
2. Add network-proximity z-scores (Cheng 2018 framework)
3. Remove circular KEGG hsa04024 tautology

### Critical fix
**P4.1** KEGG hsa04024 "TSHR present, IGF-1R absent" is circular (KEGG curates receptors
into canonical maps by definition). Remove from primary evidence; retain only as
descriptive background.

**P4.2** Compute hypergeometric p for intersection enrichment:
  - 15 genes in IGF1R ∩ TSHR ∩ TED vs. expected under null
  - 38 cochlear effectors expected ≈ 15.3 under null (IGF1R 359 × Cochlear 855 / 20000)
  - 51 insulin effectors expected ≈ 2.5 under null → 20-fold enrichment

**P4.3** Network proximity (Cheng 2018): compute drug-disease closeness metric
using STRING PPI, compare IGF1R-target-set vs TSHR-target-set proximity to TED disease genes.

### Script deliverables (see /scripts_python)
- `19_hypergeometric_enrichment.py`
- `20_network_proximity_zscore.py`

### Expected outputs
- `outputs/enrichment_v3_hypergeometric.csv`
- `outputs/network_proximity_v3.csv`

---

## 📝 PHASE 5: TABLE/FIGURE REGENERATION

Only proceed after Phase 1-4 complete and new results locked.

### Tables to regenerate
- Main Table 1: Study design — unchanged
- Main Table 2: MR+Coloc combined — UPDATED with all genes coloc + MR-PRESSO
- Supp S1: Full MR results — EXPANDED with all sensitivity methods
- Supp S2: TED genes — REPLACED with multi-source validated list
- Supp S3: Pathway genes — unchanged
- **Supp S4 (prior S5)**: DESeq2 RNA-seq results (RENUMBERED)
- **Supp S5 (prior S6)**: Off-target gene lists (RENUMBERED)
- **NEW Supp S6**: STROBE-MR checklist
- **NEW Supp S7**: Hypergeometric enrichment table
- **NEW Supp S8**: GTEx v8 parallel MR

### Figures to regenerate
- Figure 1: banner text — potentially update after Phase 1
- Figure 2: Venn gene placement unchanged, add hypergeometric p
- Figure 3 Panel A (triangulation matrix): UPDATED with new coloc+DESeq2 numbers
- Figure 3 Panels B/C: UPDATED with DESeq2 q-values
- Figure 3 Panel D: unchanged

### New figures
- **Supp Fig S1**: Coloc regional plots (LocusCompare) for TSHR + all candidates
- **Supp Fig S2**: MR forest plot with all sensitivity methods
- **Supp Fig S3**: Sample-level heatmap (z-scored, 15 genes × 5 samples)

---

## 📝 PHASE 6: MANUSCRIPT

Draft after all data/figures locked.

### Structure (Thyroid 4,000 words)
- Abstract (200w, unstructured) — emphasize triangulation framework as deliverable
- Introduction (800w) — TED/teprotumumab/IGF-1R–TSHR debate/computational gap
- Methods (1,000w) — refer to STROBE-MR supp; explicit *cis*-MR framing; DESeq2
- Results (1,400w) — Table 1 → Table 2 → Fig 2 → Fig 3 flow
- Discussion (600w) — Limitations: computational only, blood eQTL proxy,
  n=4 pilot-scale RNA-seq, TED-specific GWAS absence, no functional validation

### Cover letter
- Position as hypothesis-refining triangulation, not discovery
- Explicit restraint from therapy advocacy as strength
- Reusable triangulation-score framework as methodological contribution

---

## 🔧 DEPENDENCIES TO INSTALL

### R packages (current + NEW)
```r
# Existing
install.packages(c("TwoSampleMR", "ieugwasr", "coloc"))

# NEW — Critical methodology additions
install.packages(c(
  "MendelianRandomization",  # additional MR methods
  "MRPRESSO",                # outlier correction
  "MRlap",                   # sample overlap
  "MVMR",                    # multivariable MR
  "cause",                   # CAUSE
  "susieR"                   # SuSiE for coloc-SuSiE
))

# BiocManager for DESeq2
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "apeglm"))

# Optional
install.packages(c("pheatmap", "ggrepel", "locuscomparer"))
```

### Python packages
```bash
pip install requests pandas numpy scipy statsmodels scikit-learn networkx
pip install bioservices  # Open Targets / REST API wrappers
```

---

## 📊 VERSION TRACKING

**File naming convention**:
- `*_v1_*.ext` — original analyses (archive only, do not modify)
- `*_v3_*.ext` — upgraded analyses (use for final manuscript)

Prior `Final_Numbers_Frozen_20260421.xlsx` → retain as v1 archive.
New results → `Final_Numbers_Frozen_v3.xlsx`.

---

## ✅ DECISION POINTS FOR USER

Before proceeding, please confirm:

**D1. RNA-seq raw data availability**
  → Do you have the count matrix (gene × sample) or FASTQ files?
  → If only processed TPM is available, what was the original DE tool?
  → If Wilcoxon was manually computed, we will rebuild with DESeq2 from counts.

**D2. OpenGWAS token**
  → A new JWT token is needed (401 in screenshot).
  → Instructions in `01_token_refresh_and_setup.R`.

**D3. Computational budget**
  → All-gene coloc.susie will need ~30 min per locus × 8 loci
  → GTEx v8 download: ~50 MB per tissue × 3 tissues
  → All jobs can run locally in VS Code R Interactive

**D4. Approval to archive v1 Tables/Figures**
  → v1 tables in `/mnt/user-data/outputs/TED_TRAP_Supp_S*.docx` will be superseded.
  → Archive folder: `/Manuscript/archive_v1/`

---

## 🚦 START ORDER

Recommended execution sequence:

1. **Week 1**: Phase 1 (MR) + Phase 3 (Gene set) in parallel
2. **Week 2**: Phase 2 (RNA-seq) + Phase 4 (Pathway) in parallel
3. **Week 3**: Phase 5 (Tables/Figures)
4. **Week 4-6**: Phase 6 (Manuscript)

Phase 1 and Phase 3 are independent and can run concurrently.
Phase 2 depends on raw RNA-seq data.
Phase 4 depends on Phase 3 output (final TED gene list).
Phase 5 depends on Phases 1-4.
Phase 6 depends on everything.

---

## 📧 QUESTIONS TO RESOLVE

1. Where is the RNA-seq raw count matrix stored?
2. Preferred STRING PPI version for network proximity (v11.5 vs v12.0)?
3. Do you want multi-ancestry meta-analysis (European + Japanese) or European-only?
4. Should we archive all v1 files to `/Manuscript/archive_v1/` before regenerating?

These are needed before starting Phase 1. Once resolved, we begin coding.
