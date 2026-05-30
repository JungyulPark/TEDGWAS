# TED-TRAP v5 Week 1 — Scripts

**Working directory on Gemini's machine:** `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\scripts\`

---

## Script inventory

| Script | Purpose | Output |
|---|---|---|
| `run_week1.sh` | Master orchestrator — runs all tasks in order with pre-flight checks | (orchestrator) |
| `taskA_01_pqtl_check.R` | Task A: query UKB-PPP + deCODE for 4 target proteins | `TaskA_01_pqtl_availability_v1.csv` |
| `taskB_01_locus_extract.R` | Task B.1: extract TSHR ±1Mb cis-eQTLs from eQTLGen | `TaskB_01_tshr_locus_extraction_v1.tsv` |
| `taskB_02_ld_clumping.R` | Task B.2: LD clumping sensitivity (2 P × 3 r² × 2 ancestries = 12 cells) | `TaskB_02_ld_clumping_sensitivity_v1.csv` |
| `taskB_03_gcta_cojo.sh` | Task B.3: GCTA-COJO bash driver (EUR + EAS) | `TaskB_03_tshr_cojo_*_raw.jma.cojo` |
| `taskB_03_helper_make_ma.R` | Helper: convert TSV → COJO .ma format | `TSHR_eqtlgen_locus.ma` |
| `taskB_03_helper_parse_cojo.R` | Helper: parse COJO raw output to spec format + add LD with rs179252 | `TaskB_03_tshr_cojo_independent_signals_v1.tsv` |
| `taskB_04_susie_finemap.R` | Task B.4: SuSiE-RSS fine-mapping with local LD matrices | `TaskB_04_tshr_susie_credible_sets_v1.tsv` |
| `taskB_05_integration.R` | Task B.5: integrate COJO + SuSiE → final independent-signal list | `TaskB_05_tshr_finemap_summary_v1.csv` |

---

## Dependency graph

```
run_week1.sh
├── taskA_01_pqtl_check.R                              [independent — Day 1-2]
└── taskB sequence                                     [Day 1-4]
    ├── taskB_01_locus_extract.R                       [Day 1]
    ├── taskB_02_ld_clumping.R                         [Day 3]
    └── taskB_03_gcta_cojo.sh                          [Day 3-4]
        ├── taskB_03_helper_make_ma.R                  (called by .sh)
        └── taskB_03_helper_parse_cojo.R               (called by .sh)
    └── taskB_04_susie_finemap.R                       [Day 4]
    └── taskB_05_integration.R                         [Day 4-5]
```

Task A runs entirely in parallel with Task B (friend modification 1).

---

## Prerequisites (Gemini to verify before running)

### Software
- R ≥ 4.2 with packages: `data.table`, `tibble`, `ieugwasr`, `susieR`, `Matrix`, `rtracklayer` (optional for liftOver)
- GCTA ≥ 1.94 (Windows binary: `gcta64.exe`)
- PLINK 1.9 or 2.0
- bash (Git Bash on Windows is sufficient)

### Data files at `c:\ProjectTEDGWAS\`
- `2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` (eQTLGen, already present per locked memory)

### LD references (download once, see `00_meta/data_sources.csv`)
- 1000G Phase 3 EUR bfile for chr14 (or whole-genome)
- 1000G Phase 3 EAS bfile for chr14 (or whole-genome)
- Allele frequency `.frq` files generated with `plink --freq`

### pQTL catalog files (download once for Task A)
- `01_pqtl_check/data/ukbpp_protein_panel.csv` — extracted from Sun 2023 Supp Table 1
- `01_pqtl_check/data/decode_protein_panel.csv` — extracted from Ferkingstad 2021 Supp tables

---

## Critical hardcoded values (verified, do NOT modify)

| Value | Source |
|---|---|
| TSHR Ensembl ID `ENSG00000165409` | Ensembl GRCh37 archive |
| TSHR hg19 TSS `81,421,333` | Ensembl GRCh37 |
| TSHR hg19 end `81,612,646` | Ensembl GRCh37 |
| TSHR hg38 coordinates `80,954,989-81,146,306` | Ensembl 109+ |
| Top eQTL SNP `rs179252` | Locked v4.3 BBJ coloc result |
| eQTLGen n=31,684 | Vosa 2021, PMID 34475573 |
| UKB-PPP n=54,219 | Sun 2023, PMID 37794186 |
| deCODE n=35,559 | Ferkingstad 2021, PMID 34857953 |

---

## Critical rules

1. **Do NOT change filenames.** Downstream Week 2-6 scripts will look for these exact paths.
2. **Do NOT change CSV column names or order.** Spec v1.1 is the contract.
3. **Genome build columns are required** for every position-bearing row.
4. **`availability_status` must be one of**: `not_measured`, `measured_no_cis_pqtl`, `cis_pqtl_available`, `trans_only`.
5. **Log every run** in `00_meta/analytical_log.md` (auto-handled by scripts).
6. **If a script errors out**, log it in `analytical_log.md` and pause — do not silently skip steps.

---

## Workflow for Gemini

1. Read this README and `WEEK1_TASK_SPEC_v1.1.md` in full.
2. Verify all paths in pre-flight section of `run_week1.sh` match local installation.
3. Update paths in each script's "Config" block if needed.
4. Download required data:
   - 1000G EUR + EAS bfiles
   - UKB-PPP and deCODE supplementary tables for protein panels
5. Run: `bash run_week1.sh`
6. If errors: pause, log to `analytical_log.md`, ask Claude.
7. On success: write Task C deliverables (`decision_table.csv`, `pre_analysis_plan_v1.md`, `week1_summary_memo.md`).
8. Upload all deliverables to Claude for verification.

---

*All citations and PMIDs verified via PubMed search 2026-05-19.*
*All hardcoded coordinates verified via Ensembl GRCh37/38 archive.*
*No unverified data, no hallucinated URLs.*
