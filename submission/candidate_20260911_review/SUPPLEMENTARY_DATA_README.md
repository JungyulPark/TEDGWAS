# Supplementary data guide

Five files retain the complete primary and frequency-sensitivity results without duplicating the original and sensitivity datasets as separate uploads. Blank values mean unavailable or not estimable, not zero. Clinical tables use rounded values; these files preserve full numerical precision.

## Supplementary Data 1 Instruments

6,135 selected cis-eQTL instrument records for 2,544 genes before outcome-specific harmonization. `snp` is the rsID, `chr` and `pos_hg19` use GRCh37, `effect_allele` is the assessed allele, and `zscore` and `n_samples` are the reported exposure association and variant-specific sample size. Gene symbols, Ensembl identifiers and druggability tier are retained. Selection used P < 5×10⁻⁸ and European LD clumping r² < 0.001.

## Supplementary Data 2 Mendelian randomization

30,115 rows combine 7,219 estimable original primary associations with 22,896 rows covering the 2,544 eligible genes, three outcomes and three frequency-sensitivity scenarios. Fields include gene, outcome, scenario, estimator, instrument count, beta, SE, P value, OR and 95% CI.

- `original_reference`: the primary analysis, using European reference frequencies. Estimable genes number 2,234 BBJ, 2,505 UKB and 2,480 FinnGen.
- `paired_reference`: the reference analysis on variants retained for the matched frequency comparison.
- `paired_eqtlgen`: the same variants with cohort frequencies substituted.
- `reharmonized_eqtlgen`: cohort frequencies followed by repeat harmonization. Estimable genes number 2,232, 2,506 and 2,481, respectively.

Sensitivity scenarios include `n_iv=0` rows with blank estimates, indicating no estimable association. These rows must not be counted as null associations. The primary estimators are Wald ratio and IVW. ORs and their CIs exponentiate beta and beta ± 1.96 SE, on a disease log-odds scale per reconstructed expression unit. UKB was rescaled before MR. This file does not claim that frequency substitution reran MR-Egger, weighted-median or weighted-mode analyses; the original supporting diagnostics are in Table S1 and complete estimator results in Supplementary Data 5.

## Supplementary Data 3 Colocalization

324 rows contain nine genes, three outcomes, three shared-association priors and four frequency scenarios. `original_reference` contains all 81 primary settings; `paired_reference` and `paired_eqtlgen` compare identical variants; `available_eqtlgen` also includes variants newly covered by cohort frequencies.

Dataset 1 is the disease outcome and dataset 2 is expression. H0 means neither associated, H1 outcome only, H2 expression only, H3 both associated at distinct variants, and H4 a shared variant under the single-causal-variant model. PP.H4 is a posterior probability, not a P value. p1 = p2 = 10⁻⁴; p12 is 10⁻⁵, 5×10⁻⁶ or 10⁻⁶. `n_overlap` gives the number of matched regional variants. `top_snp` is the highest-probability variant conditional on H4, not necessarily the GWAS lead SNP. For primary rows, `outcome_min_p` gives the smallest marginal outcome P value among the variants used. H2 dominance does not imply that this P value is non-significant. Low H4 does not prove that disease association is absent.

## Supplementary Data 4 Detection limits

Twelve rows retain all four frequency scenarios in each outcome. `or_median`, `or_q1` and `or_q3` summarize the minimum detectable association at 80% power under the normal approximation. `frac_OR1_5`, `frac_OR2` and `frac_OR3` are proportions of estimable genes able to detect the named OR at that power; multiply by 100 for percentages. Alpha is 0.05/2,544 in discovery and 0.05 in the other outcomes. These are MR association-detection estimates, not the power of colocalization or the combined filter.

## Supplementary Data 5 Complete MR estimators and diagnostics

13,039 aggregate results retain all five methods from the original reference-frequency analyses: 3,256 Wald ratios, 3,963 IVW, and 1,940 each for MR-Egger, weighted median and weighted mode. This supplies the effect sizes and uncertainty behind the compact P-value diagnostic table, including the rest of the screen.

`beta`, `se`, `pvalue`, instrument counts and diagnostic statistics are unchanged from the verified analysis. `test_distribution` and `test_df` make each test explicit. The 95% interval uses the matching normal distribution for Wald/IVW/median and Student t for Egger (n_iv−2 degrees of freedom) or mode (n_iv−1). `log_or_ci_lower/upper` and `or_ci_lower/upper` give its endpoints. These method-specific intervals differ from generic normal intervals in historical export files; no coefficient or P value was changed. The P values for all 13,039 estimates were independently verified against these distributions. If exponentiation overflows for an extremely wide interval, use the finite log-scale endpoints; this indicates imprecision, not a missing estimate.

Egger-intercept and Cochran-Q columns are gene–outcome diagnostics copied across estimator rows; they are not separate tests of each estimator. Unavailable values are NA, and the Steiger columns remain NA because directionality has not been established in this package. These analyses use the original reference frequencies; they are not a new frequency-sensitivity rerun.
