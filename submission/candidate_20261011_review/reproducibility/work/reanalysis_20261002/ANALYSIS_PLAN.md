# FinnGen endpoint extension — analysis plan

Recorded 2026-10-02 before computing associations for the two additional endpoints. This is a post hoc extension of the existing study, not a prospectively registered study.

## Questions and endpoint roles

1. Do the 13 previously identified BBJ Graves disease discovery genes show directionally concordant MR associations with FinnGen R12 **E4_GRAVES_STRICT** (3,962 cases; 496,386 controls)?
2. Which of the fixed 13 genes, plus the clinically selected comparator IGF1R, have shared-variant support in this GD endpoint?
3. Does the interpretation change under **E4_GRAVES_OPHT_STRICT** (753 cases; 499,595 controls), compared with the original GRAVES_OPHT (858 cases; 499,490 controls)?

FinnGen endpoints have overlapping participants. They are phenotype comparisons, not independent replications of each other. None provides a GD-only control group. No within-GD TED effect, treatment response, difference between endpoints, or pooled cross-phenotype causal estimate will be inferred.

## Frozen inputs and tests

- Reuse 6,135 selected SNP–gene records for 2,544 genes, with their original EUR clumping and selection. Do not select instruments using new outcome P values.
- Preserve BBJ, UKB and original FinnGen estimates; reproduce the original FinnGen MR with the current executable pipeline as an identity control before accepting new estimates.
- Report all 2,544 genes for each new endpoint, including non-estimable rows. Primary focused family: the existing 13 BBJ discoveries. FinnGen GD corroboration: same direction as BBJ and two-sided P < 0.05/13. Report nominal P < 0.05 and P < 0.05/2,544 separately; the full-screen extension is exploratory. IGF1R is a comparator and is not counted as a BBJ discovery.
- MR: TwoSampleMR harmonise_data action=2; one IV, Wald ratio; multiple IVs, multiplicative random-effects IVW with underdispersion correction. Report beta, SE, OR, 95% CI, P, instrument count, Cochran Q and Q P where estimable. No clinical covariate-adjusted regression is claimed.
- Frequency sensitivity: original EUR reference; same retained SNPs under reference vs eQTLGen frequencies; full eQTLGen re-harmonization. Confirm max absolute change in exposure Z < 1e-8 and same-SNP single-IV P invariance.
- Leave-one-out: all fixed 14 genes in new endpoints with >=2 retained instruments, on primary frequencies. Remaining single SNP uses a Wald ratio. Do not treat nonsignificance after removal as proof of biological absence.
- Colocalization: all available reported cis-eQTL statistics for the fixed 14 genes, not a significance-filtered subset. Join rsID and chromosome, verify compatible allele pairs. eQTL coordinates GRCh37, FinnGen GRCh38; do not join raw position across builds. Preserve existing region definitions.
- coloc.abf model, disease dataset1 and expression dataset2: H0 neither, H1 disease only, H2 expression only, H3 both with distinct variants, H4 shared variant. p1=p2=1e-4; p12=1e-5, 5e-6, 1e-6. Outcome prior SD=0.2; expression prior SD=0.15, sdY=1. Four frequency/availability scenarios matching the previous analysis.
- PP.H4 >=0.80 is a model-based evidence convention, not a P-value or proof of causality. No SuSiE or multi-signal fine mapping in this extension.
- Validate posteriors against installed R coloc, and independently check MR estimates from retained SNP-level inputs. Report coverage and exclusions. Stop any failing stage and diagnose before continuing.

## Provenance and output boundaries

Source data and the existing manuscript/submission package are read only. New scripts, logs, analysis tables, verification files and the interpretation report go under work/reanalysis_20261002 and outputs/TED_Reanalysis_20261002. Record input and script hashes, software versions and deviations. Do not silently rewrite previously completed analyses.

Runtime preflight: ordinary sandbox R did not see the user package library; an approved read with the user's library restored access. No analysis was run with missing/replacement packages.
