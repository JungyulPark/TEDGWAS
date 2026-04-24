# TED-TRAP v3 Upgrade — Installation & Execution Guide

**Environment**: Windows + VS Code + Antigravity + R Interactive
**Working directory**: `c:\ProjectTEDGWAS\`
**Generated date**: 2026-04-23

---

## 📦 Script Inventory (17 scripts)

### Phase 1 — MR Re-analysis (R, 8 scripts)
| # | File | Purpose | Critical fix |
|---|------|---------|--------------|
| 01 | `01_token_refresh_and_setup.R` | OpenGWAS token + package install | — |
| 02 | `02_eqtlgen_exposure_extraction_all_genes.R` | Extract IVs for 8 genes | — |
| 03 | `03_mr_all_genes_full_sensitivity.R` | Full MR with all sensitivity methods | **CRIT 1, 4, 5** |
| 04 | `04_mr_scale_rescale.R` | β-scale verification (Lloyd-Jones 2018) | **CRIT 1** |
| 05 | `05_gtex_thyroid_adipose_parallel.R` | GTEx tissue-specific MR | MAJOR 7 |
| 06 | `06_coloc_all_loci.R` | Coloc for all 8 loci + SuSiE | **CRIT 3**, MAJOR 11 |
| 07 | `07_mrlap_overlap_correction.R` | Sample-overlap correction | MAJOR 9 |
| 08 | `08_mvmr_ted_given_gd.R` | Multivariable MR (TED \| GD) | MAJOR 8 |
| 09 | `09_bbj_replication.R` | BBJ ancestry replication | MAJOR 8 |

### Phase 2 — RNA-seq Re-analysis (R, 4 scripts)
| # | File | Purpose | Critical fix |
|---|------|---------|--------------|
| 10 | `10_deseq2_analysis.R` | DESeq2 + apeglm + BH-FDR | **CRIT 2** |
| 11 | `11_bh_fdr_correction.R` | Unified BH-FDR across all tests | MAJOR 12 |
| 12 | `12_gene_set_testing_candidate.R` | CAMERA/ROAST for gene sets | MAJOR 12 |
| 13 | `13_insr_novelty_check.R` | INSR verification vs Kim 2021/2024 | **CRIT 6** |

### Phase 3 — Gene Set Cross-validation (Python, 5 scripts)
| # | File | Purpose | Critical fix |
|---|------|---------|--------------|
| 13b | `13b_kim2024_insr_scan.py` | Web scan of bioRxiv/JCI Insight | CRIT 6 aux |
| 14 | `14_opentargets_query.py` | Open Targets Platform GraphQL | MAJOR 10 |
| 15 | `15_disgenet_rest_api.py` | DisGeNET v7+ REST API | MAJOR 10 |
| 16 | `16_pubtator_mining.py` | PubTator3 literature mining | MAJOR 10 |
| 17 | `17_magma_finngen_go.R` | MAGMA on FinnGen R12 GO | MAJOR 10 |
| 18 | `18_ted_gene_set_concordance.py` | Multi-source harmonization | MAJOR 10 |

### Phase 4 — Pathway/Network (Python, 2 scripts)
| # | File | Purpose | Critical fix |
|---|------|---------|--------------|
| 19 | `19_hypergeometric_enrichment.py` | Enrichment tests ✅ **ALREADY RUN** | MAJOR 15 |
| 20 | `20_network_proximity_zscore.py` | Network proximity (Cheng 2018) | MAJOR 15 |

### Phase 5 — Figure Regeneration (Python, 1+ scripts)
| # | File | Purpose |
|---|------|---------|
| 21 | `21_fig3_panelA_triangulation_v3.py` | New triangulation matrix |

---

## 🚀 Execution Order (Sequential)

### Week 1 — Setup + Phase 1 (MR)
```powershell
# Day 1 — Environment
cd c:\ProjectTEDGWAS
Rscript scripts/01_token_refresh_and_setup.R
# → Install all packages, refresh JWT token

# Day 2 — Extract instruments
Rscript scripts/02_eqtlgen_exposure_extraction_all_genes.R
# → 8 gene IV files in TrackA_MR/data/instruments/

# Day 3 — Core MR
Rscript scripts/03_mr_all_genes_full_sensitivity.R
# → 03_mr_v3_full.csv (main), 03_mr_v3_sensitivity.csv, 03_mr_v3_presso.csv

# Day 3 — Scale verification
Rscript scripts/04_mr_scale_rescale.R
# → Diagnostic log with draft Methods text

# Day 4 — Coloc (the longest step)
Rscript scripts/06_coloc_all_loci.R
# → 06_coloc_v3_all_loci.csv
# Requires eQTLGen + FinnGen full summary stats (~500 MB total download first)

# Day 5 — GTEx, MRlap, MVMR, BBJ (parallel)
Rscript scripts/05_gtex_thyroid_adipose_parallel.R
Rscript scripts/07_mrlap_overlap_correction.R
Rscript scripts/08_mvmr_ted_given_gd.R
Rscript scripts/09_bbj_replication.R

# Day 6 — BH correction
Rscript scripts/11_bh_fdr_correction.R
# → 11_mr_v3_bh.csv, 11_coloc_v3_tiers.csv
```

### Week 2 — Phase 2 (RNA-seq) + Phase 3 (Gene set)

**IMPORTANT**: Phase 2 requires you to provide the RNA-seq count matrix first.
See Decision D1 below.

```powershell
# Day 7 — DESeq2 (requires count matrix)
Rscript scripts/10_deseq2_analysis.R
# → 10_deseq2_v3_results.csv, 10_deseq2_v3_candidate_set.csv

# Day 8 — Gene set testing
Rscript scripts/12_gene_set_testing_candidate.R
# → 12_camera_results.csv, 12_roast_results.csv

# Day 8 — INSR novelty
Rscript scripts/13_insr_novelty_check.R
# Optional: python scripts/13b_kim2024_insr_scan.py

# Day 9-10 — Gene set cross-validation (parallel)
python scripts/14_opentargets_query.py
python scripts/15_disgenet_rest_api.py   # requires DisGeNET API key
python scripts/16_pubtator_mining.py
Rscript scripts/17_magma_finngen_go.R    # requires MAGMA binary
python scripts/18_ted_gene_set_concordance.py
```

### Week 3 — Phase 4 (Pathway) + Phase 5 (Figures)
```powershell
# Phase 4
python scripts/19_hypergeometric_enrichment.py   # ALREADY RUN, results in hand
python scripts/20_network_proximity_zscore.py    # requires STRING PPI download

# Phase 5 — Figure regeneration
python scripts/21_fig3_panelA_triangulation_v3.py
# More scripts to come (Panels B, C, new Supp figures)
```

---

## 🔧 Prerequisite Downloads

### External data (user must download)

**1. OpenGWAS JWT token** (script 01 helps)
- https://api.opengwas.io/profile/ → Generate new token

**2. eQTLGen full summary stats** (for script 06 coloc)
- https://www.eqtlgen.org/cis-eqtls.html
- File: `2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` (~2 GB)
- Place in: `TrackA_MR/data/eqtlgen/`

**3. FinnGen R12 Graves ophthalmopathy** (for script 06)
- https://r12.finngen.fi/pheno/E4_GRAVES_OPHT
- Click "Download GWAS summary statistics"
- Place in: `TrackA_MR/data/finngen/`

**4. GTEx v8 tissue-specific eQTLs** (for script 05)
- https://gtexportal.org/home/downloads/adult-gtex/qtl
- Required files:
  - `Thyroid.v8.signif_variant_gene_pairs.txt.gz`
  - `Adipose_Visceral_Omentum.v8.signif_variant_gene_pairs.txt.gz`
- Place in: `TrackA_MR/data/gtex_v8/`

**5. LDSC reference files** (for script 07 MRlap)
- https://ibg.colorado.edu/mtag/
- Download `eur_w_ld_chr.tar.bz2`
- Extract to: `TrackA_MR/data/ldsc_ref/`

**6. STRING v12 human PPI** (for script 20 network proximity)
- https://string-db.org/cgi/download?species_text=Homo+sapiens
- Files:
  - `9606.protein.links.v12.0.txt.gz` (~90 MB)
  - `9606.protein.info.v12.0.txt.gz` (~10 MB)
- Place in: `TrackB_Network/data/`

**7. DisGeNET API key** (for script 15)
- Register: https://www.disgenet.com/
- Set environment variable: `setx DISGENET_API_KEY "your_key"` (Windows)

**8. MAGMA binaries** (for script 17)
- https://cncr.nl/research/magma/
- Download: `magma_win.zip` (Windows) or `magma_v1.10.zip`
- Also download: 1000G EUR reference + NCBI38.gene.loc
- Place in: `TrackA_MR/tools/magma/`

---

## 🎯 Decision Points (USER INPUT NEEDED)

### D1 — RNA-seq raw count matrix ❓
**Question**: Do you have integer count matrix for DESeq2?

**Options**:
- **A)** Yes, I have `counts_gene.tsv` (gene × sample, integer counts)
  - ✅ Provide path; Phase 2 runs as-is
- **B)** I have TPM only
  - ⚠️ DESeq2 requires counts. Options:
    - B.1) Re-run upstream quantification (STAR + featureCounts, salmon + tximport)
    - B.2) Use length-scaled TPM via tximport
    - B.3) Use limma-voom on log2(TPM+1) as approximation (disclose in Methods)
- **C)** I have only processed/normalized data
  - ❌ Need to re-do from FASTQ level, or disclose the limitation

**Action needed**: Tell me which option (A/B/C) applies.

### D2 — OpenGWAS token refresh 🔑
**Question**: Are you able to refresh the JWT token when prompted by script 01?
- Yes → run script 01 interactively
- No → paste token into `.Renviron` manually

### D3 — LD reference for coloc-SuSiE ❓
**Question**: Do you have a 1000G EUR LD matrix for chromosome 14 (TSHR locus)?
- Yes → script 06 can run coloc-SuSiE
- No → script 06 will skip coloc-SuSiE, run coloc.abf only (acceptable)

### D4 — Archive v1 tables/figures? 📁
**Question**: Approve archiving old v1 outputs?
- Yes → I create `/Manuscript/archive_v1/` and move old files
- No → keep v1 alongside v3 (can clutter but safer)

---

## 📊 First Results (Already Obtained)

### Phase 4 Hypergeometric Enrichment ✅
Results file: `TrackB_Network/results/19_hypergeometric_enrichment.csv`

| Test | Observed | Expected | Fold | p | BH q |
|------|----------|----------|------|---|------|
| IGF-1R pathway ∩ TED (59) | 12 | 0.11 | **112.99×** | 1.56e-22 | 3.90e-22 |
| TSHR pathway ∩ TED (59) | 11 | 0.10 | **112.99×** | 1.01e-20 | 1.68e-20 |
| Three-way intersection ∩ TED | 8 | 0.03 | **301.32×** | 3.14e-20 | 3.93e-20 |
| IGF-1R downstream ∩ Hearing loss | 38 | 15.35 | 2.48× | 2.73e-07 | 2.73e-07 |
| IGF-1R downstream ∩ KEGG insulin | 51 | 2.46 | **20.74×** | 8.30e-54 | 4.15e-53 |

**Draft Results text for manuscript**:
> *Pathway co-annotation with the TED disease-gene set was highly enriched for both IGF-1R (observed 12 vs expected 0.11, 113-fold, hypergeometric p=1.56e-22) and TSHR (observed 11 vs expected 0.10, 113-fold, p=1.01e-20) pathways. The three-way intersection of IGF-1R pathway, TSHR pathway, and TED genes was enriched 301-fold (p=3.14e-20), indicating convergent biology rather than independent signal. The expanded IGF-1R downstream effector set was 20.7-fold enriched within KEGG insulin signaling (hsa04910; p=8.30e-54), consistent with the pharmacological overlap that underlies the insulin-signaling off-target hypothesis for IGF-1R inhibitors.*

### INSR Novelty Verdict ✅
See: `TrackA_MR/results/13_insr_novelty_verdict.md`

**Verdict**: INSR is NOT a DEG in Kim 2024 JCI Insight main text; pathway-level insulin response IS reported. Use cautious "to our knowledge" language — full verdict and approved wording in the verdict document.

---

## 📅 Expected Timeline

| Week | Phase | Deliverable |
|------|-------|-------------|
| 1 | P1 (MR) | Complete MR with all sensitivity |
| 2 | P2 + P3 (RNA-seq, gene set) | DESeq2 + curated gene set |
| 3 | P4 + P5 start (pathway, figures) | Enrichment + new figures |
| 4 | P5 cont. + P6 start (tables, manuscript) | All tables v3, Results draft |
| 5-6 | P6 (manuscript) | Full manuscript + cover letter |

Total: **4-6 weeks** from start.

---

## 🛟 If anything breaks

Each script writes to `TrackA_MR/logs/` or `TrackB_Network/logs/` with the step number.
Check the log first, then report back the last lines of the log here.

Common issues + fixes:
- **OpenGWAS 401**: Token expired → run script 01 again
- **Package won't install on Windows**: Use `install.packages(..., type="binary")` or install Rtools
- **MR-PRESSO fails**: Requires ≥4 IVs; for TSHR (2 IVs) this is expected and handled
- **Coloc runs too slow**: Reduce locus window from ±500 kb to ±100 kb
- **MAGMA Windows build**: Use WSL or the precompiled Windows binary

---

## 📝 What's NOT YET covered by these scripts (future work)

These are handled in Phase 6 (Manuscript) — not scripted:

- Full Table 1 v3 regeneration (docx)
- Full Table 2 v3 (coloc ALL loci) regeneration (docx)
- Supplementary Tables S1-S8 v3 (docx)
- Figures 1, 2 banner updates (if results change)
- New Supp Figures: S1 coloc regional plots, S2 MR forest, S3 sample heatmap
- Methods section write-up
- Results section write-up
- Discussion with proper limitations
- Cover letter
- STROBE-MR checklist (supp)

These will be done AFTER Phase 1-5 results are in hand, because the numbers drive the text.
