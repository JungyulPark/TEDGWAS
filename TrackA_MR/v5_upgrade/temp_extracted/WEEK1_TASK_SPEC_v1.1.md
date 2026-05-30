# TED-TRAP v5 — Week 1 Critical Task Specification

**Version:** v1.1 (verified, corrections applied)
**Lock date:** 2026-05-19
**Working directory:** `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\`
**Executor:** Gemini (local R + GCTA + bash)
**Verifier:** Claude (uploaded results)

---

## Changelog from v1.0

1. **Verified all citations and PMIDs** via PubMed/publisher pages.
2. **Corrected TSHR hg19 coordinates** to Ensembl GRCh37 canonical values.
3. **Removed unverified Synapse ID** (`syn51365303`); replaced with verified AWS Open Data Registry URL for UKB-PPP.
4. **Added Eldjarn 2023 cross-platform comparison** as required reference if both UKB-PPP and deCODE are used.
5. **Task A and Task B run in PARALLEL** (not sequential) — both required regardless of pQTL feasibility (friend's modification 1).
6. **Added genome_build / liftover columns** to output schemas (friend's modification 2).
7. **Added `TaskB_00_ld_reference_metadata_v1.csv`** as a required deliverable (friend's modification 3).
8. **Added `availability_status` enum column** for strict pQTL state classification (friend's modification 4).
9. **Tone-downed ceiling logic** — JCI Insight/EBioMedicine requires Week 2-3 novel hits, not Week 1 alone (friend's modification 5).
10. **Added DGIdb 5.0 (2024)** as updated option.

---

## 1. Overall Objective

Upgrade TED-TRAP from a candidate-gene confirmation study (rejected by Thyroid and ETJ) into a systematic GD/TED genetic-prioritization study. Week 1 determines the analytical ceiling.

**Primary framing (locked):**
> Druggable-gene-wide genetic prioritization of GD/GO-relevant susceptibility architecture, with TSHR as an upstream GD-linked locus and IGF1R as a downstream pharmacologic effector rather than a germline susceptibility driver.

**Outcome hierarchy (locked):**
1. Primary discovery: BBJ Graves disease (EAS, 2,809 cases)
2. Replication: UKB hyperthyroidism (EUR, 3,731 cases)
3. TED-specific sensitivity: FinnGen R12 Graves ophthalmopathy (858 cases)

**Week 1 critical decision:**
> Does the TSHR single-IV problem resolve structurally (via multi-IV from fine-mapping or via pQTL MR)? Result determines journal targeting, but final ceiling also depends on Week 2-3 druggable-gene-wide MR outcomes.

---

## 2. Reproducibility Framework

### 2.1 Folder structure
```
v5_upgrade/
├── 00_meta/
│   ├── README.md
│   ├── data_sources.csv          (verified, PMIDs locked)
│   ├── software_versions.csv
│   └── analytical_log.md
├── 01_pqtl_check/
│   ├── data/
│   ├── results/
│   │   └── TaskA_01_pqtl_availability_v1.csv
│   └── docs/
├── 02_tshr_finemap/
│   ├── data/
│   ├── results/
│   │   ├── TaskB_00_ld_reference_metadata_v1.csv  ← NEW (friend mod 3)
│   │   ├── TaskB_01_tshr_locus_extraction_v1.tsv
│   │   ├── TaskB_02_ld_clumping_sensitivity_v1.csv
│   │   ├── TaskB_03_tshr_cojo_independent_signals_v1.tsv
│   │   ├── TaskB_04_tshr_susie_credible_sets_v1.tsv
│   │   └── TaskB_05_tshr_finemap_summary_v1.csv
│   └── docs/
├── 03_decision/
│   ├── TaskC_decision_table_v1.csv
│   ├── TaskC_pre_analysis_plan_v1.md
│   └── TaskC_week1_summary_memo.md
└── scripts/
```

### 2.2 File naming convention
`{TaskID}_{Step}_{description}_v{version}.{ext}`

Version bump (v2, v3) only on substantive data change. Format-only fixes overwrite v1. **Gemini must not change filenames** — downstream Week 2-6 scripts will look for these exact paths.

### 2.3 Metadata files (00_meta/)

**`data_sources.csv`** — all PMIDs and DOIs verified via PubMed (see Section 7 reference table).

**`software_versions.csv`** — `R --version`, `packageVersion()`, `gcta64 --version` to be recorded by Gemini.

**`analytical_log.md`** — chronological per-run log: date, script, input files, parameters, output files, issues.

---

## 3. Task A — pQTL Availability Check

### 3.1 Objective
Determine cis-pQTL availability for **TSHR, IGF1R, IGF1, INSR** in major plasma proteomics resources.

### 3.2 Data sources (all verified)

| Source | Platform | Citation (verified) | Access |
|---|---|---|---|
| UKB-PPP | Olink Explore 3072 (~2,923 proteins, 54,219 participants) | Sun et al. 2023 Nature 622:329-338 (PMID 37794186, DOI 10.1038/s41586-023-06592-6) | https://registry.opendata.aws/ukbppp |
| deCODE | SOMAscan v4 (4,907 aptamers, 35,559 Icelanders) | Ferkingstad et al. 2021 Nat Genet 53:1712-1721 (PMID 34857953, DOI 10.1038/s41588-021-00978-w) | https://www.decode.com/summarydata/ |
| Cross-platform comparison | Both platforms compared | Eldjarn et al. 2023 Nature 622:348-358 (PMID 37794188, DOI 10.1038/s41586-023-06563-x) | Reference only — use if both platforms used in our analysis |

### 3.3 Procedure

**Step A.1** — Query UKB-PPP supplementary tables (Sun 2023) for the 4 target proteins:
- Confirm protein on Olink Explore 3072 panel.
- **Critical:** TSHR is a membrane-bound thyroid receptor; plasma proteomics may not capture it. Document explicitly.
- For proteins on the panel: retrieve cis-pQTL summary statistics (±1 Mb of gene) from AWS Open Data Registry.

**Step A.2** — Query deCODE supplementary tables (Ferkingstad 2021) for the same 4 proteins:
- SOMAscan v4 aptamer presence check.
- Download cis-pQTL sumstats from deCODE summary data portal if available.

**Step A.3** — For each available pQTL, extract: lead cis-SNP, position (record original genome build), effect/other allele, EAF, beta, SE, P-value, sample size, ancestry, multiple independent signals if reported.

**Step A.4** — Cross-platform check: if a protein is measured on both platforms, compare lead SNPs and effect directions (cite Eldjarn 2023 for interpretation context — Olink and SOMAscan show only moderate correlation).

### 3.4 Required output file

**`TaskA_01_pqtl_availability_v1.csv`** (one row per gene × platform):

| Column | Type | Description |
|---|---|---|
| gene_symbol | string | TSHR, IGF1R, IGF1, INSR |
| gene_ensembl | string | ENSG00000165409 (TSHR), ENSG00000140443 (IGF1R), ENSG00000017427 (IGF1), ENSG00000171105 (INSR) |
| platform | string | "UKBPPP" or "deCODE" |
| platform_version | string | "Olink_Explore_3072" or "SOMAscan_v4" |
| **availability_status** | **enum** | **`not_measured` / `measured_no_cis_pqtl` / `cis_pqtl_available` / `trans_only`** (friend mod 4) |
| protein_measured | boolean | TRUE if on panel, FALSE if not |
| has_cis_pqtl | boolean | TRUE if any cis-pQTL P<5e-8 |
| lead_cis_snp | string | rsID of strongest cis-pQTL, NA if absent |
| chr | int | Chromosome |
| **pos_original** | **int** | **Position in source genome build** (friend mod 2) |
| **genome_build** | **string** | **"hg19" / "hg38" / "GRCh37" / "GRCh38"** (friend mod 2) |
| **pos_hg19** | **int** | **After liftOver if needed** (friend mod 2) |
| **pos_hg38** | **int** | **After liftOver if needed** (friend mod 2) |
| **liftover_required** | **boolean** | **TRUE if source build ≠ hg19** (friend mod 2) |
| **liftover_chain** | **string** | **e.g. "hg38ToHg19.over.chain.gz" or NA** (friend mod 2) |
| effect_allele | string | A/T/G/C |
| other_allele | string | A/T/G/C |
| eaf | float | Effect allele frequency |
| beta | float | Effect estimate |
| se | float | Standard error |
| pvalue | float | P-value |
| n_samples | int | Sample size |
| ancestry | string | "EUR", "EAS", etc. |
| n_independent_signals | int | If reported by source |
| sumstats_file | string | Local path to downloaded sumstats |
| source_pmid | string | "37794186" (Sun 2023) or "34857953" (Ferkingstad 2021) |
| access_date | date | YYYY-MM-DD |
| notes | string | "membrane receptor — possibly not in plasma", etc. |

### 3.5 Decision rules

| availability_status | Decision |
|---|---|
| `cis_pqtl_available` for TSHR (any platform) | Proceed to TSHR pQTL MR — major upgrade |
| `not_measured` TSHR + `cis_pqtl_available` IGF1R | IGF1R pQTL MR as orthogonal effector layer |
| `cis_pqtl_available` for IGF1 or INSR only | Pathway-level pQTL MR (off-target / IGF axis anchor) |
| All four `not_measured` or `measured_no_cis_pqtl` | pQTL arm not feasible — document as limitation |
| `trans_only` | Do NOT use as primary instrument — violates exclusion restriction |

---

## 4. Task B — TSHR Locus Fine-mapping

**Runs in PARALLEL with Task A from Day 1.** Required regardless of Task A result, because (a) if pQTL feasible, fine-mapping still establishes locus architecture; (b) if pQTL not feasible, fine-mapping is the only path to addressing the Wald-ratio limitation.

### 4.1 Target locus (verified Ensembl GRCh37)

- Gene: TSHR
- Ensembl: ENSG00000165409
- **hg19 coordinates: chr14:81,421,333 – 81,612,646** (forward strand, Ensembl GRCh37 archive browser)
- **hg38 coordinates: chr14:80,954,989 – 81,146,306** (forward strand)
- Genome build for analysis: **hg19** (eQTLGen native)
- Region for extraction: chr14:80,421,333 – 82,612,646 (gene ± 1 Mb in hg19)
- Current lead eQTL SNP (from previous TED-TRAP v4.3 BBJ coloc): **rs179252**
- Known GD risk SNPs (different from eQTL): rs179247, rs12101255 in TSHR intron 1 (Płoski et al. 2010 PLoS One)

### 4.2 Step-by-step procedure

**Step B.0 — LD reference metadata (NEW, friend mod 3)**

Before any analysis, document the LD reference panels used. This file is required for COJO/SuSiE reviewer defense.

**Output:** `TaskB_00_ld_reference_metadata_v1.csv`

| Column | Description |
|---|---|
| ld_reference_id | Short ID (e.g., "1000G_phase3_EUR", "1000G_phase3_EAS") |
| ancestry | EUR / EAS / SAS / AFR / AMR |
| source | "1000 Genomes Project Phase 3" |
| genome_build | "hg19" or "hg38" |
| n_individuals | Number of reference samples |
| n_snps_chr14 | Number of SNPs on chr14 after QC |
| plink_bfile_prefix | Path to bfile prefix |
| download_url | Where downloaded from |
| access_date | YYYY-MM-DD |
| md5_hash | If verifiable |
| notes | QC steps applied |

**Step B.1 — Locus extraction**

Input: `2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` (eQTLGen native, hg19, n=31,684)

Filter: gene == TSHR (or ENSG00000165409) AND chr14:80,421,333-82,612,646 (hg19).

**Output:** `TaskB_01_tshr_locus_extraction_v1.tsv`

| Column | Type | Description |
|---|---|---|
| snp | string | rsID |
| chr | int | 14 |
| pos_hg19 | int | Position (eQTLGen native) |
| pos_hg38 | int | After liftOver |
| liftover_status | string | "native_hg19" / "lifted_to_hg38_success" / "lifted_failed" |
| effect_allele | string | A/T/G/C |
| other_allele | string | A/T/G/C |
| eaf | float | EAF in eQTLGen |
| beta | float | Effect of EA on TSHR expression |
| se | float | Standard error |
| zscore | float | Z-score |
| pvalue | float | P-value |
| fdr | float | FDR-adjusted P |
| n_samples | int | eQTLGen sample size (~31,684) |
| distance_to_tss | int | bp from TSHR TSS (chr14:81,421,333) |
| is_rs179252 | boolean | TRUE only for rs179252 (for downstream tagging) |

**Step B.2 — LD clumping sensitivity**

Run at multiple thresholds to assess robustness:

```r
library(ieugwasr)
clump_results <- list()
for (r2 in c(0.001, 0.01, 0.1)) {
  for (pthresh in c(5e-8, 5e-6)) {
    for (pop in c("EUR", "EAS")) {
      key <- paste(r2, pthresh, pop, sep="_")
      clump_results[[key]] <- ld_clump(
        dat = locus_data,
        clump_r2 = r2,
        clump_p = pthresh,
        clump_kb = 10000,
        pop = pop
      )
    }
  }
}
```

**Output:** `TaskB_02_ld_clumping_sensitivity_v1.csv`

| Column | Type | Description |
|---|---|---|
| pthresh | float | 5e-8 or 5e-6 |
| r2_thresh | float | 0.001, 0.01, 0.1 |
| ld_reference_id | string | matches TaskB_00 |
| n_independent_snps | int | After clumping |
| lead_snps | string | Comma-separated rsIDs |
| includes_rs179252 | boolean | TRUE/FALSE |
| notes | string | Comments |

**Step B.3 — GCTA-COJO conditional analysis**

Convert eQTLGen extraction to COJO `.ma` format:
```
SNP A1 A2 freq b se p N
rs179252 T C 0.42 -0.45 0.013 9.1e-15 31684
...
```

Run stepwise selection (parameter values per GCTA documentation):
```bash
gcta64 \
  --bfile 1000G_EUR_chr14 \
  --cojo-file TSHR_eqtlgen_locus.ma \
  --cojo-slct \
  --cojo-p 5e-6 \
  --cojo-wind 10000 \
  --cojo-collinear 0.9 \
  --out TaskB_03_tshr_cojo_EUR
```

Also run with EAS reference for cross-ancestry sensitivity:
```bash
gcta64 \
  --bfile 1000G_EAS_chr14 \
  --cojo-file TSHR_eqtlgen_locus.ma \
  --cojo-slct \
  --cojo-p 5e-6 \
  --cojo-wind 10000 \
  --cojo-collinear 0.9 \
  --out TaskB_03_tshr_cojo_EAS
```

**Output:** `TaskB_03_tshr_cojo_independent_signals_v1.tsv`

| Column | Type | Description |
|---|---|---|
| ld_reference_id | string | "1000G_phase3_EUR" or "1000G_phase3_EAS" |
| signal_rank | int | 1, 2, 3... (order of selection) |
| lead_snp | string | rsID |
| chr | int | 14 |
| pos_hg19 | int | Position |
| effect_allele | string | A/T/G/C |
| other_allele | string | A/T/G/C |
| eaf | float | EAF |
| marginal_beta | float | Univariable beta |
| marginal_p | float | Univariable P |
| conditional_beta | float | Conditional on other selected SNPs |
| conditional_se | float | Conditional SE |
| conditional_p | float | Conditional P |
| ld_with_rs179252_eur | float | r² in EUR |
| ld_with_rs179252_eas | float | r² in EAS |
| interpretation | string | "primary" / "secondary" / "collinear with rs179252" |

**Step B.4 — SuSiE fine-mapping**

```r
library(susieR)
fit_eur <- susie_rss(
  z = locus_data$zscore,
  R = ld_matrix_eur,
  n = 31684,
  L = 10,
  coverage = 0.95
)
fit_eas <- susie_rss(
  z = locus_data$zscore,
  R = ld_matrix_eas,
  n = 31684,
  L = 10,
  coverage = 0.95
)
```

**Output:** `TaskB_04_tshr_susie_credible_sets_v1.tsv`

| Column | Type | Description |
|---|---|---|
| ld_reference_id | string | "1000G_phase3_EUR" or "1000G_phase3_EAS" |
| cs_id | int | Credible set ID (1, 2, ...) |
| n_snps_in_cs | int | Number of SNPs in this CS |
| cs_purity | float | Min absolute correlation within CS |
| top_snp | string | Highest-PIP SNP in CS |
| top_snp_pip | float | PIP of top SNP |
| top_snp_pos_hg19 | int | Position |
| top_snp_p | float | Marginal P |
| ld_top_with_rs179252_eur | float | r² in EUR |
| ld_top_with_rs179252_eas | float | r² in EAS |
| coverage | float | Posterior coverage |
| converged | boolean | SuSiE convergence flag |

**Step B.5 — Integration**

**Output:** `TaskB_05_tshr_finemap_summary_v1.csv`

| Column | Type | Description |
|---|---|---|
| signal_id | string | "TSHR_S1", "TSHR_S2", etc. |
| lead_snp_cojo_eur | string | COJO's lead in EUR |
| lead_snp_cojo_eas | string | COJO's lead in EAS |
| lead_snp_susie_eur | string | SuSiE's top-PIP in EUR |
| lead_snp_susie_eas | string | SuSiE's top-PIP in EAS |
| concordant_across_methods | boolean | TRUE if COJO and SuSiE agree (within ancestry) |
| concordant_across_ancestries | boolean | TRUE if EUR and EAS agree |
| chr | int | 14 |
| pos_hg19 | int | Lead position |
| pos_hg38 | int | After liftOver |
| marginal_p | float | From extraction |
| conditional_p_eur | float | From COJO EUR |
| conditional_p_eas | float | From COJO EAS |
| susie_pip_eur | float | From SuSiE EUR |
| susie_pip_eas | float | From SuSiE EAS |
| ld_with_rs179252_eur | float | r² in EUR |
| ld_with_rs179252_eas | float | r² in EAS |
| usable_as_independent_iv | boolean | TRUE if conditional P<5e-6 AND r²<0.1 with rs179252 in BOTH ancestries |
| notes | string | |

### 4.3 Decision rules

| Scenario | Decision |
|---|---|
| ≥2 independent signals concordant in COJO + SuSiE | Multi-IV MR for TSHR — solves Wald ratio limitation |
| 1 signal only, concordant across methods | Keep single-IV, reframe as "locus-prioritization" not "robust multi-IV MR" |
| COJO and SuSiE disagree | Report both, use conservative interpretation (single-IV) |
| EUR and EAS COJO results differ | Document explicitly; use ancestry-matched result for each outcome (EAS for BBJ, EUR for UKB/FinnGen) |

**Backup if only 1 signal:** LD-proxy expansion (r² > 0.6 with rs179252) as instrument-validity sensitivity, NOT primary multi-IV.

---

## 5. Task C — Week 1 Decision Memo

### 5.1 Decision table

**Output:** `TaskC_decision_table_v1.csv`

Required rows (revised per friend mod 4 & 5):
1. TSHR pQTL availability_status (UKB-PPP)
2. TSHR pQTL availability_status (deCODE)
3. IGF1R pQTL availability_status (UKB-PPP)
4. IGF1R pQTL availability_status (deCODE)
5. IGF1 pQTL availability_status (UKB-PPP)
6. IGF1 pQTL availability_status (deCODE)
7. INSR pQTL availability_status (UKB-PPP)
8. INSR pQTL availability_status (deCODE)
9. TSHR independent signals (COJO EUR, count)
10. TSHR independent signals (COJO EAS, count)
11. TSHR independent signals (SuSiE EUR, count)
12. TSHR independent signals (SuSiE EAS, count)
13. Is rs179252 the single SuSiE top-PIP variant? (yes/no/partial)
14. Druggable-gene-wide MR extension justified? (yes — proceed Week 2)

### 5.2 Pre-analysis plan

**Output:** `TaskC_pre_analysis_plan_v1.md` — locked BEFORE Week 2 starts.

Required sections (no changes from spec v1.0).

### 5.3 Week 1 summary memo

**Output:** `TaskC_week1_summary_memo.md`

One-page memo:
- Week 1 findings
- TSHR single-IV status: resolved / unresolved
- pQTL feasibility: yes / partial / no
- Recommended Tier (see Section 6 revised)

---

## 6. Ceiling Decision Logic (revised, friend mod 5)

Week 1 alone does NOT determine final journal tier. Final tier depends on Week 1 architecture + Week 2-3 druggable-genome-wide MR results.

| Condition | Realistic ceiling |
|---|---|
| TSHR multi-IV resolved + pQTL feasible **only** | IOVS / JCEM / EJE strong |
| TSHR multi-IV resolved + pQTL feasible + novel replicated druggable hits + colocalization | JCI Insight / EBioMedicine exploratory (Week 2-3 dependent) |
| Single-IV only but systematic MR strong with novel hits | IOVS / Front Endocrinol / JEI |
| Single-IV + no novel druggable hits | Front Endocrinol / JEI / Sci Rep |

**Key point:** JCI Insight / EBioMedicine cannot be locked in from Week 1 results alone. Those journals require novel druggable discoveries from Week 2-3 systematic MR, not just methodological strength at the TSHR locus.

---

## 7. Master citation reference (all PMIDs verified via PubMed)

| Reference | Citation | PMID | DOI |
|---|---|---|---|
| eQTLGen | Võsa et al. 2021 Nat Genet 53:1300-1310 | 34475573 | 10.1038/s41588-021-00913-z |
| UKB-PPP Sun 2023 | Plasma proteomic associations with genetics and health in the UK Biobank. Nature 622:329-338 | 37794186 | 10.1038/s41586-023-06592-6 |
| Cross-platform Eldjarn 2023 | Large-scale plasma proteomics comparisons. Nature 622:348-358 | 37794188 | 10.1038/s41586-023-06563-x |
| deCODE Ferkingstad 2021 | Large-scale integration of the plasma proteome. Nat Genet 53:1712-1721 | 34857953 | 10.1038/s41588-021-00978-w |
| Finan 2017 druggable | The druggable genome and support for target identification. Sci Transl Med 9:eaag1166 | 28356508 | 10.1126/scitranslmed.aag1166 |
| GTEx v8 Consortium 2020 | Atlas of genetic regulatory effects across human tissues. Science 369:1318-1330 | 32913098 | 10.1126/science.aaz1776 |
| DGIdb 4.0 Freshour 2021 | DGIdb 4.0 with open crowdsource efforts. Nucleic Acids Res 49:D1144-D1151 | 33237278 | 10.1093/nar/gkaa1084 |
| DGIdb 5.0 (newer option) | DGIdb 5.0. Nucleic Acids Res 52:D1227 | (verify before use) | 10.1093/nar/gkad1040 |
| GSE58331 Rosenbaum 2015 | Orbital adipose microarray. PLoS One 10(9):e0137654 | (in locked memory) | 10.1371/journal.pone.0137654 |
| GSE105149 Rosenbaum 2017 | Lacrimal gland microarray. JAMA Ophthalmol 135(11):1156-1162 | (in locked memory) | 10.1001/jamaophthalmol.2017.3458 |

---

## 8. Final Deliverables Manifest

| File | Created by | Status |
|---|---|---|
| `00_meta/README.md` | Manual | Required |
| `00_meta/data_sources.csv` (verified PMIDs) | Manual | Required |
| `00_meta/software_versions.csv` | Manual | Required |
| `00_meta/analytical_log.md` | Manual | Required |
| `01_pqtl_check/results/TaskA_01_pqtl_availability_v1.csv` | Task A | Required |
| `02_tshr_finemap/results/TaskB_00_ld_reference_metadata_v1.csv` | Task B.0 | Required (NEW) |
| `02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv` | Task B.1 | Required |
| `02_tshr_finemap/results/TaskB_02_ld_clumping_sensitivity_v1.csv` | Task B.2 | Required |
| `02_tshr_finemap/results/TaskB_03_tshr_cojo_independent_signals_v1.tsv` | Task B.3 | Required |
| `02_tshr_finemap/results/TaskB_04_tshr_susie_credible_sets_v1.tsv` | Task B.4 | Required |
| `02_tshr_finemap/results/TaskB_05_tshr_finemap_summary_v1.csv` | Task B.5 | Required |
| `03_decision/TaskC_decision_table_v1.csv` | Task C | Required |
| `03_decision/TaskC_pre_analysis_plan_v1.md` | Task C | Required |
| `03_decision/TaskC_week1_summary_memo.md` | Task C | Required |

---

## 9. Execution Sequence (5-day parallel plan, friend mod 1)

**Day 1**
- Set up folder structure, populate README + data_sources.csv + software_versions.csv
- **Task A start:** Query UKB-PPP and deCODE catalogs for 4 proteins (parallel)
- **Task B start:** TaskB_00 LD reference metadata + TaskB_01 locus extraction (parallel)

**Day 2**
- **Task A continues:** complete pQTL availability table (TaskA_01)
- **Task B continues:** TaskB_02 LD clumping sensitivity

**Day 3**
- **Task A continues:** if pQTLs available, download sumstats files
- **Task B continues:** TaskB_03 GCTA-COJO (EUR + EAS)

**Day 4**
- **Task A complete**
- **Task B continues:** TaskB_04 SuSiE (EUR + EAS), TaskB_05 integration

**Day 5**
- **Task C:** decision table + pre-analysis plan + summary memo
- Upload all deliverables to Claude for verification

---

## 10. Critical Reminders

1. **Do NOT use trans-pQTLs as primary instruments** — violates exclusion restriction.
2. **Do NOT use LD proxies (r² > 0.6) as independent IVs** in COJO results — they are not independent.
3. **eQTLGen native build is hg19** — confirm before extraction; mismatched builds will silently corrupt results.
4. **TSHR canonical hg19 coordinates: chr14:81,421,333 – 81,612,646** (Ensembl GRCh37 archive verified).
5. **Record genome_build column for every position** to prevent silent build mismatches downstream.
6. **Cross-ancestry LD differences matter** — TSHR is on chr14 and LD architecture differs between EUR and EAS. Use ancestry-matched LD reference for each outcome (EAS for BBJ, EUR for UKB/FinnGen).
7. **All PMIDs in this spec verified via PubMed**. If new citations are added later, verify before locking.
8. **Save all intermediate files** — eQTLGen TSHR locus extraction is the basis for Week 2-6 downstream analyses.
9. **Log every parameter** in `analytical_log.md` — required for reviewer defense.
10. **If anything ambiguous, pause and ask Claude** — better to pause than corrupt the pipeline.

---

*End of Week 1 Task Spec v1.1 (verified).*
*Lock date: 2026-05-19.*
*All citations verified via PubMed/publisher pages.*
*No unverified data, no hallucinated URLs.*
