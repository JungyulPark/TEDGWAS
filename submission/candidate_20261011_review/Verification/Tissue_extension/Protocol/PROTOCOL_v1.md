# TSHR allele audit and prospective tissue-extension protocol, version 1

## Status and scope
This is a post hoc extension of an existing blood-eQTL manuscript, not preregistration of the original study. The existing blood MR results, earlier tissue work and published thymus hypothesis are already known. This document fixes the new tissue-analysis choices before new tissue-specific TSHR/IGF1R results are retrieved or examined in this extension. Earlier analyses must be labelled historical, not prospectively specified. Any deviations require a dated rationale and must retain the original analysis. No outcome-driven tissue selection is allowed.

Current execution: audit rs179252 (blood TSHR instrument) against rs179247 (published thymus variant), using raw exposure, disease data and signed reference LD. Tissue MR is a subsequent phase; this protocol does not imply that it has been performed. No SuSiE, MVMR or IGF-1 ligand MR is included.

## Part 1: allele audit
1. Extract both variants, where available, from original eQTLGen and the five existing disease datasets. Preserve reported alleles, coordinate build, beta/SE or Z, frequency and source identity. Do not substitute an illustrative task-specification row for data.
2. Verify that the eQTL effect allele really increases blood expression, align each disease effect to that allele, and independently reconstruct the single-SNP Wald estimate using the manuscript's original exposure scaling and UKB conversion. Compare full precision estimates with the current canonical results.
3. Recover both SNPs from the existing 1000 Genomes EUR and EAS reference genotypes. Report allele-labelled signed genotype correlation and r-squared; infer coupling only after allele orientation is explicit. If phased haplotypes are unavailable, identify phase as inferred from unphased genotype LD, not directly observed. EUR and EAS are mandatory; no best-population selection. Missingness and sample count are reported. JPT-only analysis, if later added, is an explicitly secondary check.
4. Compare the expression-increasing blood allele's linked rs179247 allele with the allele-specific thymus study. Classify as directionally consistent, directionally discordant, or unresolved. All classifications are reportable. High r-squared alone does not establish direction, causal identity, thymic mediation, or an orbital mechanism.

## Part 2: fixed tissue panel
Genes: TSHR and IGF1R, each examined in all five GTEx v8 tissues below, regardless of instrument availability or direction. GTEx v8 is selected to match the existing manuscript's reference version; newer versions are not substituted after viewing results.

| Tissue identifier | Role and biological question |
|---|---|
| Thyroid | Thyroid regulation, relevant to GD; not a thymus surrogate |
| Cells_Cultured_fibroblasts | Fibroblast regulation; cultured non-orbital cells, not orbital fibroblasts |
| Adipose_Subcutaneous | Adipose regulation; not orbital fat |
| Adipose_Visceral_Omentum | A second fixed adipose context; not orbital fat |
| Whole_Blood | GTEx blood comparator for the original eQTLGen exposure |

Thymus and orbital tissue are biologically important but not assumed to have usable GTEx eQTL summary data. Their availability is recorded explicitly. The in-house 4 TED/1 control RNA-seq cannot estimate a population eQTL or provide independent MR validation. A lack of instruments is not evidence of no biological effect.

ImmuNexUT: catalogue feasibility first for all 28 published immune-cell subsets and both genes. Do not pick subsets because of favourable gene associations. A separate dated amendment must fix the release, cell-subset list, access conditions, analysis family and correction before reading gene-specific association results or executing immune-cell MR. This availability work is not a completed ancestry-matched replication. Healthy-only and pooled autoimmune samples, if both available, must not be silently interchanged.

## Disease hierarchy and multiplicity
Primary family: 2 genes x 5 fixed GTEx tissues x FinnGen GD (E4_GRAVES_STRICT) = 10 planned tests. This is ancestry-aligned European exposure/outcome analysis with the directly relevant GD endpoint. Two-sided Bonferroni threshold: 0.05/10 = 0.005; retain the denominator even when some cells cannot be estimated.

Secondary family: the same 10 gene-tissue combinations x 4 endpoints (BBJ GD, UKB hyperthyroidism, original FinnGen GO, strict FinnGen GO) = 40 planned tests; threshold 0.05/40 = 0.00125. BBJ is a cross-ancestry portability check, UKB has a broader phenotype, and the two GO endpoints share population controls and overlap with FinnGen GD. Do not pool them as independent replications. Report raw P and family-adjusted P (min[1,m*P]); P<0.05 alone is nominal evidence, not family-corrected evidence. The endpoint labels/builds/source files remain those of the audited October 2 analysis.

## Instrument and data-quality criteria
- Use complete gene-region eQTL summary statistics with documented effect allele, beta and SE, sample size, genome build and normalization. Significant-pair files alone may support instrument selection but cannot support regional colocalization.
- Cis window: within 1 Mb of the source release's annotated TSS. Select P<5e-8 and F=(beta/SE)^2>10; require MAF>=0.01 where reliable frequencies exist. No relaxed-P rescue analysis in this protocol. Record a missing MAF rather than impute one silently.
- Clump by EUR reference LD r-squared<0.001 within 10 Mb; use the same fixed rule in all GTEx tissues. Selection uses exposure statistics only, never outcome P or MR direction. Keep one independently selected SNP if that is all that survives.
- Match genome builds before allele harmonization. Resolve palindromic SNPs only with compatible frequencies; exclude unresolved A/T or C/G SNPs (including intermediate frequency 0.42-0.58) and incompatible/multiallelic mappings. No automatic proxy replacement for missing outcome variants; report coverage and failures. If two or more selected SNPs map to one variant, resolve duplication before estimation.
- One SNP: Wald ratio with first-order SE. Two or more: IVW with multiplicative random-effects uncertainty bounded below by fixed-effect uncertainty; report Cochran Q, number of retained instruments, minimum/median F and leave-one-out when >=2 SNPs. MR-Egger is not a routine add-on for sparse cis instruments. Account for any residual LD or flag the combination as unsuitable rather than treating correlated SNPs as independent.
- Report OR, 95% CI and two-sided P for every estimable combination. OR units follow the actual tissue expression normalization, not an assumed universal SD. Compare direction across tissues; do not directly equate magnitudes if scales differ. Reusing the same outcome GWAS does not create independent replication of the original blood analysis.

## Ancillary colocalization and interpretation
When complete compatible regional data are available, run the existing single-causal-variant coloc model in the fixed TSS +/-1 Mb region using unfiltered overlapping variants, with p1=p2=1e-4 and p12=1e-5, plus p12=1e-6 and 1e-4 sensitivity. Report H0-H4, SNP count, coverage of exposure/outcome lead SNPs and exact window; do not impute absent lead signals or select a window to maximize H4. H4>=0.80 is supportive under that model, not a P value or proof of causality. Failure to reach it is not automatically evidence against association. Multi-signal uncertainty remains a limitation; no SuSiE is scheduled here.

Keep statistical strength and direction separate: (a) corrected evidence, (b) nominal evidence only, (c) imprecise/no nominal evidence with CI, (d) not estimable with reason. In each estimable category state concordant/opposite to the already known blood direction (TSHR lower odds; IGF1R higher odds), without treating significant versus nonsignificant as a between-tissue difference. Strong opposite-direction evidence must receive equal narrative prominence. With no estimable instruments, state unavailable, never negative.

## Required reporting, stopping and provenance
Publish a 50-row planned-result grid, including failures, with sources/releases, sample counts, normalization, coordinates, allele transformations, SNP count, F, beta/SE, OR/CI/P, adjusted P, direction, coloc status and exclusion reason. Report all ten primary cells and all forty secondary cells; show an availability matrix even if no MR is possible. Supplementary material retains detailed instrument rows and all sensitivity settings, not only favourable results.

Stop expansion after these fixed targets. Additional tissues/releases or thresholds require a dated prospective amendment and an exploratory label. Do not hold submission indefinitely while waiting for unavailable controlled-access data. This analysis cannot distinguish TED risk within GD without GD-only controls, establish thymic tolerance mediation, or predict treatment response.
