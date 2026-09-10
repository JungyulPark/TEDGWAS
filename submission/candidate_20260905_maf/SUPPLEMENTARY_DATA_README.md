# Supplementary data dictionary

These files contain public summary-statistic results, not individual participant records. Blank / NA entries mean unavailable or not computed, not zero.

## Supplementary Data 1 — Instruments

6,135 selected cis-eQTL instruments for 2,544 genes, before outcome-specific harmonization. Gene symbols and Ensembl identifiers identify the exposure. `chr` and `pos_hg19` use GRCh37/hg19. `snp` is the rsID; `effect_allele` is the assessed allele for `zscore`, and `other_allele` is the other allele. `pvalue` is the cis-eQTL association P value, `n_samples` is the SNP-level eQTL sample size, and `finan_tier` is the druggability classification. The primary selection uses P < 5×10⁻⁸ and European LD clumping r² < 0.001.

## Supplementary Data 2 — Primary MR

One primary estimate per estimable gene and outcome: 2,234 BBJ, 2,505 UKB and 2,480 FinnGen rows. `method` specifies Wald ratio or IVW; `n_iv` is the number of instruments used. `beta` and `se` are the MR estimate and standard error on the disease log-odds scale per unit of reconstructed exposure expression. `pvalue` is the primary association P value; `or`, `or_ci_lower` and `or_ci_upper` report exp(beta) and its 95% confidence interval. `egger_intercept`, `egger_intercept_p`, `cochran_q` and `cochran_q_p` are available diagnostics. The `steiger_*` fields are inherited empty fields: Steiger testing was not applied. Harmonization-count and scale columns retain the historical pipeline metadata. UKB estimates were rescaled from the original linear-model scale using the case fraction. All rows are restricted to genes in Supplementary Data 1.

## Supplementary Data 3 — Colocalization

81 rows: nine genes × three outcomes × three p12 priors. p1 = p2 = 10⁻⁴; p12 is the row-specific shared-association prior. Dataset 1 is the outcome GWAS and dataset 2 the eQTL exposure. `PP.H0` = neither associated; H1 = outcome only; H2 = eQTL only; H3 = both associated with distinct variants; H4 = shared causal variant under the single-causal-variant model. Low H4 is not proof of no disease association. `n_overlap` is the number of variants used; allele-mismatch and palindromic fields describe exclusions. `top_snp` and `top_SNP.PP.H4` summarize the highest conditional shared-variant SNP posterior. `cs95_*` fields are the coloc conditional H4 credible set, not SuSiE fine-mapping results. Minimum association P values are provided for each dataset. Exposure beta/variance use European reference frequencies with sdY = 1. These are the unchanged primary reference-frequency results; the completed cohort-frequency sensitivity is provided separately in Supplementary Data 4–6.

The manuscript and Tables S2/S9 define the phenotype hierarchy, model assumptions and prior sensitivity. Full-precision values are preserved in these CSVs; tables use rounded display values.


## Supplementary Data 4 — Frequency sensitivity of primary MR

Three scenarios for all 2,544 eligible genes and three outcomes (22,896 rows). `paired_reference` and `paired_eqtlgen` use identical retained variants to isolate frequency scaling. `reharmonized_eqtlgen` repeats action=2 harmonization. `n_iv=0` with blank estimates means no estimable association, not a null estimate. Methods are Wald ratio or IVW only; this frequency substitution did not rerun MR-Egger, weighted-median or weighted-mode estimators. The original estimators remain in Table S6. Outcome scales and thresholds are unchanged.

## Supplementary Data 5 — Frequency sensitivity of colocalization

324 rows = nine loci × three outcomes × three priors × four scenarios. `original_reference` reproduces Supplementary Data 3; `paired_reference` restricts reference results to cohort-frequency-covered variants; `paired_eqtlgen` substitutes frequencies on those identical variants; `available_eqtlgen` additionally includes newly frequency-covered variants. No original overlapping SNP lacked the cohort AF. `PP.H0`–`PP.H4`, priors, overlap counts and top SNPs use the definitions above.

## Supplementary Data 6 — Frequency sensitivity of MR detection thresholds

Four scenarios × three outcomes. Alpha remains 0.05/2,544 in discovery and 0.05 otherwise. `frac_OR1_5`, `frac_OR2`, `frac_OR3` are proportions of estimable genes whose normal-approximation detectable OR is no larger than the specified value; multiply by 100 for percentage. `or_q1`, `or_median`, `or_q3` exponentiate the corresponding detectable log-odds quantiles. These metrics are distinct from combined-filter power.
