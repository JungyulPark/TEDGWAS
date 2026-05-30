# TED-TRAP v5 — Week 3 Task Specification: Colocalization

**Version:** v1.0
**Lock date:** 2026-05-26
**Working directory:** `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\`
**Executor:** Gemini (local R) | **Verifier:** Claude
**Prerequisite:** Week 2 COMPLETE & ARCHIVED (13 hits, robust_novel=0, Option 1 defensive)

---

## 1. Objective

For the defensible backbone genes from Week 2, test whether the cis-eQTL signal
and the GD/TED GWAS signal share a single causal variant (colocalization).
This distinguishes a true shared signal from coincidental overlap / LD confounding,
and is the standard expectation reviewers have for any MR finding presented as
locus-level evidence.

**What coloc adds to the paper:**
- TSHR: confirm shared causal variant (strengthens the upstream-anchor claim)
- IGF1R: test whether the effector signal colocalizes (strengthens or qualifies it)
- CTLA4: positive-control coloc (pipeline validation)
- chr16p11.2 cluster (optional): formally demonstrate LD-spillover (NOT independent signals)

**Honest framing:** coloc is confirmatory here, not a new discovery. PP.H4 results
are reported as-is, including any that are weak/ambiguous (as in v4.3 where UKB TSHR
coloc was ambiguous, PP.H3=0.774).

---

## 2. Locked decisions

| Decision | Lock |
|---|---|
| Method | `coloc.abf` (single-causal-variant assumption); coloc.susie only if a locus clearly has multiple signals AND LD matrix available |
| eQTL dataset type | `type="quant"`, eQTLGen (N=31,684). sdY unknown → supply MAF + N (coloc estimates sdY) |
| GWAS dataset type | `type="cc"`, supply `s` = case proportion |
| Window | gene ±1 Mb (hg19), same as instrument extraction |
| PP.H4 threshold | > 0.80 = strong colocalization (pre-specified) |
| Primary coloc outcome | BBJ Graves (discovery). Also run FinnGen GO (TED-specific) |
| eQTLGen has no MAF column | use 1000G EUR .frq (already generated: g1000_eur_freq.frq) for MAF |

**Case proportions (s) for cc datasets:**
- BBJ Graves: s = 2809 / 175465 = 0.0160
- UKB hyperthyroid: s = 3731 / 484598 = 0.0077
- FinnGen GO: s = 858 / (858 + controls) → CONFIRM control count from FinnGen file

---

## 3. Genes to colocalize

| Gene | chr | Priority | Rationale | Outcome(s) |
|------|-----|----------|-----------|------------|
| TSHR | 14 | 1 | upstream anchor; v4.3 BBJ PP.H4=0.951 — reproduce | BBJ + FinnGen GO |
| IGF1R | 15 | 1 | effector, multi-IV; highest value to test | BBJ + FinnGen GO |
| CTLA4 | 2 | 2 | positive control | BBJ + FinnGen GO |
| HSD3B7 | 16 | 3 (optional) | chr16p11.2 cluster — show LD spillover | BBJ |
| PRSS36 | 16 | 3 (optional) | chr16p11.2 cluster — show LD spillover | BBJ |
| VKORC1 | 16 | 3 (optional) | chr16p11.2 cluster — show LD spillover | BBJ |

Priority 1-2 are essential. Priority 3 (chr16p11.2) only if reviewers likely to
probe; it would show these are NOT independent (high mutual PP.H4 or shared lead).

---

## 4. Inputs (all already local, verified)

| Input | Path | Status |
|---|---|---|
| eQTLGen full cis-eQTL | `c:/ProjectTEDGWAS/TrackA_MR/data/2019-12-11-cis-eQTLsFDR-...txt.gz` | rs179252 anchor verified |
| BBJ Graves harmonised | `04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz` | MD5 325a1e... |
| FinnGen GO | `c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz` | local |
| EUR MAF | `04_druggable_mr/data/g1000_eur_freq.frq` | generated Week 2 |

---

## 5. Procedure (per gene × outcome)

1. Extract eQTLGen cis-window (gene ±1 Mb) for the gene → SNP, P, Z, AssessedAllele, OtherAllele, N.
2. Extract outcome GWAS same window → SNP, P (or beta+se), effect/other allele.
3. Merge on rsID; harmonize alleles (drop SNPs that don't align; flag palindromic).
4. Attach MAF from 1000G EUR .frq.
5. Build coloc datasets:
   - eQTL: `list(snp=, pvalues=, type="quant", N=31684, MAF=)`
     (or beta/varbeta if using Z→beta with MAF; pvalues+MAF+N is acceptable)
   - GWAS: `list(snp=, pvalues=, type="cc", s=<case_prop>, N=<total>, MAF=)`
6. Run `coloc.abf(dataset1=eQTL, dataset2=GWAS)`.
7. Record PP.H0–H4, nsnps, top SNP (max SNP.PP.H4).

---

## 6. Output

`05_coloc_candidates/TaskE_01_coloc_results_v1.csv`

| Column | Description |
|---|---|
| gene_symbol | |
| outcome | BBJ_Graves / FinnGen_GO |
| n_snps_overlap | SNPs in both datasets after harmonization |
| pp_h0 | no association either |
| pp_h1 | assoc eQTL only |
| pp_h2 | assoc GWAS only |
| pp_h3 | both, distinct causal variants |
| pp_h4 | both, shared causal variant |
| top_snp | max SNP.PP.H4 |
| coloc_verdict | "strong" (H4>0.8) / "ambiguous" (H3≈H4) / "distinct" (H3>0.8) / "weak" |
| notes | |

Plus per-gene regional data saved for optional locus plots (Supp figures).

---

## 7. Interpretation guide (honest)

- **PP.H4 > 0.8** → strong shared-causal-variant evidence (e.g. expected for TSHR/BBJ)
- **PP.H3 > 0.8** → distinct causal variants (eQTL and GWAS signals are different)
- **PP.H3 ≈ PP.H4** (neither > 0.8) → ambiguous; possible LD-proxy tagging — report honestly, do NOT call it definitive coloc (this is what v4.3 UKB TSHR showed)
- **PP.H4 low + PP.H1/H2 high** → one trait associated, other not in this window
- **chr16p11.2 trio**: if they share the same top SNP / high mutual overlap → confirms LD-cluster artifact, supports excluding them as independent novel hits

---

## 8. Verified citations

| Reference | Citation | PMID/DOI |
|---|---|---|
| coloc method | Giambartolomei et al. 2014 PLoS Genet 10:e1004383 | 24830394 |
| coloc priors update | Wallace 2020 PLoS Genet 16:e1008720 | 32310995 |
| eQTLGen | Vosa et al. 2021 Nat Genet 53:1300 | 34475573 |
| (outcomes as Week 2) | Sakaue 2021; FinnGen R12 | 34594039 |

---

## 9. Critical reminders

1. eQTL = quant, GWAS = cc with correct `s`. Getting type/s wrong corrupts PP.
2. MAF from 1000G EUR (eQTLGen has no MAF column).
3. Confirm FinnGen GO control count for s before running FinnGen coloc.
4. Report ambiguous/weak results honestly — coloc that fails is still informative.
5. Single-causal-variant assumption: coloc.abf is fine for confirmation. Only use
   coloc.susie if a locus genuinely has multiple independent signals (TSHR does NOT
   — Week 1 confirmed single signal).
6. rsID position note (rs179252 eQTLGen 81.4 vs BBJ 81.0 region) — rsID match is
   what matters; document.
7. If any coloc fails to run (e.g. <50 overlapping SNPs), STOP and report, don't force.

---

*End of Week 3 Task Spec v1.0. Citations verified 2026-05-26.*
*Builds on Week 1 (TSHR single-IV) + Week 2 (TSHR/IGF1R/CTLA4 backbone, novel=0).*
