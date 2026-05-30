# TED-TRAP v5 — Week 2 Task Specification: Druggable-Genome-Wide MR

**Version:** v1.0
**Lock date:** 2026-05-22
**Working directory:** `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\`
**Executor:** Gemini (local R)
**Verifier:** Claude
**Prerequisite:** Week 1 COMPLETE (TSHR single-IV confirmed, archived)

---

## 1. Objective

Transform TED-TRAP from a ~20-gene candidate study into a **systematic druggable-genome-wide genetic prioritization** of GD/TED-relevant targets. This directly answers the desk-rejection criticism that the candidate-gene scope was too narrow.

**One-sentence framing:**
> A druggable-genome-wide cis-eQTL Mendelian randomization screen across BBJ Graves disease (discovery), UKB hyperthyroidism (replication), and FinnGen Graves ophthalmopathy (TED-specific sensitivity), prioritizing actionable targets and positioning TSHR as a known upstream anchor while characterizing IGF1R as a downstream pharmacologic effector.

**Critical principle:** IGF1R and IGF1 are part of the druggable universe and enter the screen automatically. They are NOT pre-highlighted as candidates. Their rank within the systematic screen is reported objectively. This avoids the "TSHR but not IGF1R" framing that contributed to rejection.

---

## 2. Locked Decisions (do not deviate)

| Decision | Lock |
|---|---|
| Druggable gene list | **Finan 2017 (4,479 genes) as MR base.** DGIdb 5.0 used only for downstream druggability annotation of hits, NOT to expand the MR gene set. |
| eQTL instrument source | **eQTLGen blood (n=31,684) single source for discovery.** GTEx tissue only for top-hit sensitivity. |
| Bonferroni denominator | **Number of genes with ≥1 valid cis-eQTL instrument** (not all 4,479). Pre-specified here to prevent post-hoc adjustment. |
| Outcome hierarchy | BBJ Graves = discovery; UKB hyperthyroid = replication; FinnGen GO = TED-specific sensitivity (LOCKED from Week 1) |
| Novel hit definition | Bonferroni-significant in BBJ discovery + nominal P<0.05 same-direction in UKB replication + (colocalization PP.H4>0.8 assessed in Week 3) |
| IGF1R/IGF1 handling | Enter screen as ordinary druggable genes; reported by rank, not pre-flagged |

---

## 3. Reproducibility Framework (continues from Week 1)

```
v5_upgrade/
├── 00_meta/
│   ├── data_sources.csv          (append Week 2 sources)
│   ├── software_versions.csv     (append biomaRt, etc.)
│   ├── analytical_log.md         (append Week 2 runs)
│   └── Week2_file_manifest_v1.csv (NEW)
├── 01_pqtl_check/                (Week 2 Task A completes this)
├── 04_druggable_mr/              (NEW — Week 2 main)
│   ├── data/
│   ├── results/
│   └── docs/
├── 05_coloc_candidates/          (NEW — Week 2 produces candidate list; coloc in Week 3)
└── scripts/
```

File naming: `{TaskID}_{Step}_{description}_v{version}.{ext}`
Same rules as Week 1: no filename changes, no column changes, log every run.

---

## 4. Task A — pQTL Availability (Week 1 carryover, IGF1 priority)

Runs in **parallel** with Task D-1/D-2 (gene list build). Completes the 8 pending cells from Week 1.

### A.1 — Download protein panel lists
- UKB-PPP (Olink Explore 3072): protein list from Sun 2023 Nature 622:329-338 supplementary, or Olink official panel
- deCODE (SOMAscan v4): aptamer list from Ferkingstad 2021 Nat Genet 53:1712-1721 supplementary

### A.2 — Check 4 proteins against panels
For TSHR, IGF1R, IGF1, INSR: confirm presence/absence. **IGF1 is a secreted hormone — most likely present; verify first.**

### A.3 — For present proteins, retrieve cis-pQTL
Download cis-pQTL sumstats (±1 Mb), record lead SNP / allele / beta / SE / P / N.

### Output: `TaskA_01_pqtl_availability_v1.csv` (UPDATE existing, replace `pending` with confirmed status)
Same 26-column schema as Week 1 spec. `availability_status` enum: `not_measured` / `measured_no_cis_pqtl` / `cis_pqtl_available` / `trans_only`.

### Decision rule
- IGF1 cis-pQTL available → IGF1 pQTL MR feasible (orthogonal IGF-axis evidence) — high value
- IGF1R cis-pQTL available → IGF1R pQTL MR (effector-node characterization)
- TSHR cis-pQTL available → bonus orthogonal layer (low expectation; membrane receptor)

---

## 5. Task D — Druggable-Genome-Wide MR (Week 2 core)

### D.1 — Build druggable gene master list

**Input:** Finan et al. 2017 Sci Transl Med 9:eaag1166 (PMID 28356508) supplementary — the druggable genome (4,479 genes, tiers 1/2/3A/3B).

**Steps:**
1. Download Finan 2017 supplementary gene list.
2. Normalize to current HGNC symbols + Ensembl Gene IDs (GRCh37) via biomaRt:
   ```r
   library(biomaRt)
   grch37 <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
   coords <- getBM(
     attributes=c("hgnc_symbol","ensembl_gene_id","chromosome_name",
                  "start_position","end_position"),
     filters="hgnc_symbol", values=druggable_symbols, mart=grch37)
   ```
3. Annotate druggability tier from Finan + cross-reference DGIdb 5.0 (Cannon 2024, PMID 37953380) for current drug-gene interactions (annotation only).

**Output:** `TaskD_01_druggable_gene_master_v1.csv`

| Column | Description |
|---|---|
| gene_symbol | HGNC symbol |
| ensembl_gene_id | ENSG ID |
| chr | Chromosome |
| start_hg19 | Gene start (GRCh37) |
| end_hg19 | Gene end (GRCh37) |
| finan_tier | Finan 2017 druggability tier (1/2/3A/3B) |
| dgidb_n_interactions | DGIdb 5.0 drug-gene interaction count (annotation) |
| dgidb_drugs | Comma-separated approved/known drugs (annotation) |
| in_mr_screen | TRUE if gene proceeds to instrument extraction |
| notes | |

### D.2 — Extract cis-eQTL instruments for all druggable genes

**Input:** eQTLGen full cis-eQTL sumstats (hg19) + TaskD_01 gene coordinates.

**Per gene:**
- Window: cis = gene ±1 Mb (hg19)
- Primary instrument P-threshold: **5e-8**
- Sensitivity P-threshold: **5e-6**
- LD clumping: r²<0.001, EUR reference (eQTLGen is EUR-based — EUR only, per Week 1 lesson)
- Record genes with 0 instruments (needed for Bonferroni denominator)

**Output:** `TaskD_02_eqtl_instruments_v1.csv`

| Column | Description |
|---|---|
| gene_symbol | |
| ensembl_gene_id | |
| n_instruments_5e8 | IV count at P<5e-8 after clumping |
| n_instruments_5e6 | IV count at P<5e-6 after clumping |
| instrument_snps_5e8 | Comma-separated rsIDs |
| has_valid_instrument | TRUE if n_instruments_5e8 ≥ 1 |
| single_iv | TRUE if exactly 1 IV (Wald ratio only, like TSHR) |
| notes | |

**Bonferroni denominator** = count of `has_valid_instrument == TRUE`. Record this number explicitly in `TaskD_02_denominator_note.txt`.

### D.3 — Run MR across 3 outcomes

For each gene with ≥1 valid instrument, run MR against:
1. **BBJ Graves disease** (discovery, EAS, 2,809 cases) — Sakaue 2021 Nat Genet 53:1415-1424 (PMID 34594039)
2. **UKB hyperthyroidism** (replication, EUR, 3,731 cases) — same Sakaue 2021
3. **FinnGen R12 Graves ophthalmopathy** (TED-specific sensitivity, 858 cases) — local .gz

**Methods:**
- Multi-IV genes: IVW primary + MR-Egger + weighted median + weighted mode
- Single-IV genes: Wald ratio only (flag clearly, like TSHR)
- Steiger directionality filter
- Cochran Q heterogeneity where applicable

**CRITICAL ancestry note:** BBJ is EAS; eQTLGen instruments are EUR-derived. This is a known cross-ancestry limitation. Harmonize carefully on allele + flag. The discovery in BBJ with EUR instruments is a sensitivity consideration to document, not hide.

**Output:** `TaskD_03_mr_results_3outcomes_v1.csv`

| Column | Description |
|---|---|
| gene_symbol | |
| ensembl_gene_id | |
| outcome | "BBJ_Graves" / "UKB_hyperthyroid" / "FinnGen_GO" |
| method | "IVW" / "MR-Egger" / "WM" / "Wmode" / "Wald" |
| n_iv | Instrument count used |
| beta | MR estimate |
| se | |
| pvalue | |
| or | exp(beta) |
| or_ci_lower | |
| or_ci_upper | |
| egger_intercept | (multi-IV only) |
| egger_intercept_p | |
| cochran_q | |
| cochran_q_p | |
| steiger_direction | "correct" / "reverse" / NA |
| steiger_p | |
| notes | |

### D.4 — Hit classification

Apply pre-specified novel hit definition:

**Output:** `TaskD_04_hit_classification_v1.csv`

| Column | Description |
|---|---|
| gene_symbol | |
| bbj_ivw_p | Discovery P (IVW or Wald) |
| bbj_bonferroni_sig | TRUE if bbj_ivw_p < 0.05/denominator |
| ukb_ivw_p | Replication P |
| ukb_same_direction | TRUE if UKB beta same sign as BBJ + P<0.05 |
| finngen_go_p | TED-specific sensitivity P |
| finngen_go_same_direction | TRUE if same sign |
| hit_class | "robust_novel" / "discovery_only" / "known_anchor" / "not_significant" |
| is_known_gd_locus | TRUE for TSHR, CTLA4, etc. (annotate from literature) |
| druggable_tier | from TaskD_01 |
| coloc_candidate | TRUE if hit_class in (robust_novel, known_anchor) → goes to Week 3 coloc |
| notes | |

**hit_class definitions:**
- `robust_novel`: Bonferroni in BBJ + same-direction P<0.05 in UKB + NOT a known GD locus
- `known_anchor`: Bonferroni in BBJ + known GD locus (TSHR, CTLA4, HLA region, etc.)
- `discovery_only`: Bonferroni in BBJ but fails UKB replication
- `not_significant`: fails BBJ Bonferroni

### D.5 — Colocalization candidate list (Week 3 handoff)

**Output:** `TaskD_05_coloc_candidates_v1.csv`
- All genes with `coloc_candidate == TRUE`
- Their lead cis-eQTL SNP + locus coordinates
- Which outcome(s) to colocalize against
- This list is the input to Week 3 coloc; Week 2 does NOT run coloc itself.

---

## 6. Week 2 Deliverables Manifest

| File | Task | Status |
|---|---|---|
| `01_pqtl_check/results/TaskA_01_pqtl_availability_v1.csv` | A | UPDATE (resolve pending) |
| `04_druggable_mr/results/TaskD_01_druggable_gene_master_v1.csv` | D.1 | NEW |
| `04_druggable_mr/results/TaskD_02_eqtl_instruments_v1.csv` | D.2 | NEW |
| `04_druggable_mr/results/TaskD_02_denominator_note.txt` | D.2 | NEW |
| `04_druggable_mr/results/TaskD_03_mr_results_3outcomes_v1.csv` | D.3 | NEW |
| `04_druggable_mr/results/TaskD_04_hit_classification_v1.csv` | D.4 | NEW |
| `05_coloc_candidates/TaskD_05_coloc_candidates_v1.csv` | D.5 | NEW |
| `00_meta/Week2_pre_analysis_addendum_v1.md` | — | NEW (lock criteria before D.3 runs) |

---

## 7. Execution Sequence (suggested, parallel where possible)

**Day 1-2 (parallel):**
- Task A: download protein panels, check IGF1 first
- Task D.1: download Finan 2017, normalize gene symbols + coordinates

**Day 3-4:**
- Task D.2: extract cis-eQTL instruments for all druggable genes (heavy compute)
- Record Bonferroni denominator

**Day 5-7:**
- Task D.3: run 3-outcome MR (heaviest compute)
- Task A: finalize pQTL table

**Day 8-9:**
- Task D.4: hit classification
- Task D.5: coloc candidate list
- Write Week2_pre_analysis_addendum BEFORE looking at D.3 results in detail

**Day 10:**
- Upload all deliverables to Claude for verification
- Week 2 decision: how many robust_novel hits? → determines journal tier

---

## 8. Verified Citations (PubMed-confirmed)

| Reference | Citation | PMID | DOI |
|---|---|---|---|
| Finan druggable genome | Finan et al. 2017 Sci Transl Med 9:eaag1166 | 28356508 | 10.1126/scitranslmed.aag1166 |
| DGIdb 5.0 | Cannon et al. 2024 Nucleic Acids Res 52(D1):D1227-D1235 | 37953380 | 10.1093/nar/gkad1040 |
| eQTLGen | Vosa et al. 2021 Nat Genet 53:1300-1310 | 34475573 | 10.1038/s41588-021-00913-z |
| BBJ Graves + UKB hyperthyroid | Sakaue et al. 2021 Nat Genet 53:1415-1424 | 34594039 | 10.1038/s41588-021-00931-x |
| UKB-PPP | Sun et al. 2023 Nature 622:329-338 | 37794186 | 10.1038/s41586-023-06592-6 |
| deCODE | Ferkingstad et al. 2021 Nat Genet 53:1712-1721 | 34857953 | 10.1038/s41588-021-00978-w |
| 1000G Phase 3 | Auton et al. 2015 Nature 526:68-74 | 26432245 | 10.1038/nature15393 |
| TwoSampleMR | Hemani et al. 2018 eLife 7:e34408 | 29846171 | 10.7554/eLife.34408 |

---

## 9. Critical Reminders

1. **Bonferroni denominator = genes with valid instruments**, pre-specified BEFORE seeing results. Record the number in a separate note file.
2. **eQTLGen is EUR — use EUR LD reference only.** EAS clumping/conditioning produces artifacts (Week 1 lesson).
3. **BBJ is EAS ancestry; eQTLGen instruments are EUR.** Cross-ancestry MR is a documented limitation — harmonize on allele carefully, flag explicitly.
4. **IGF1R/IGF1 are NOT pre-highlighted.** Report their rank objectively within the systematic screen.
5. **Coloc is Week 3, not Week 2.** Week 2 produces only the candidate list.
6. **Write the pre-analysis addendum BEFORE examining D.3 results in detail** to prevent cherry-picking.
7. **Single-IV genes use Wald ratio only** — flag them exactly as TSHR was flagged.
8. **Log every run** in analytical_log.md.
9. **If anything ambiguous, pause and ask Claude.**

---

*End of Week 2 Task Spec v1.0.*
*All citations verified via PubMed 2026-05-22.*
*Builds on Week 1 (TSHR single-IV confirmed, locus-level prioritization framing locked).*
