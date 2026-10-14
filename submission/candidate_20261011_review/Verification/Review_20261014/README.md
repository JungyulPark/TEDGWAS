# 2026-10-14 review revision: verification records

This folder records what was checked for the revision that resolved the
13 October review (`internal/REVIEW_codex_package_20261011.md`).

## derived_values.csv
Every number newly added to the manuscript in this revision, recomputed from
`Supplementary_Dataset.xlsx` (and the retained aggregate MR export) by
`reproducibility/derive_review_values_20261014.py`. The integrity audit checks
that each displayed string appears in the master.

## finngen_coloc_reproduction/
Independent reproduction of the new-endpoint results for *TSHR*, *IGF1R* and
*CTLA4* (FinnGen R12 E4_GRAVES_STRICT and E4_GRAVES_OPHT_STRICT), with the
original GRAVES_OPHT endpoint as an identity control.

- `A/` and `B/`: two separately written Python implementations of coloc.abf
  (disease = dataset 1, prior SD 0.2; expression = dataset 2, prior SD 0.15,
  sdY = 1; p1 = p2 = 1e-4; p12 = 1e-5, 5e-6, 1e-6) and of Wald/IVW MR.
  FinnGen regions were re-downloaded by tabix range request from the public R12
  bucket; eQTL regions are `TrackA_MR/data/eqtl_{TSHR,IGF1R,CTLA4}_all.csv`.
- `compare/`: reconciliation of A, B, the September values and the Codex values.

Outcome:
- Variant matching identical to the reported analysis (TSHR 6,708/6,675;
  IGF1R 7,894/7,850; CTLA4 5,023/4,993); same top SNP in all 27 settings.
- 1000 Genomes EUR was not reachable, so gnomAD v2.1.1 European frequencies were
  substituted. Identity control: max |ΔPP.H4| 4.1e-4, max |ΔPP| 2.8e-3.
- New endpoints: max |ΔPP.H4| vs reported 1.2e-3; no setting crosses 0.80.
  Restricting TSHR to the 6,579 variants in the 1000G EUR reference reproduces the
  reported values to ≤2.3e-6.
- IGF1R PP.H2/PP.H3 differ by up to 0.003 at the third decimal (same shift as the
  identity control: frequency source), so the printed IGF1R H2/H3 values are
  reproduced to within 0.003, not to the printed digit.
- MR for all three genes and both new endpoints reproduces to numerical precision
  on the original exposure scale.

Not included: FinnGen regional extracts, gnomAD extracts and per-SNP inputs
(source data are not redistributed). Machine paths in scripts are replaced by
`<scratch>`/`<home>` placeholders; the scripts need those inputs re-fetched.
