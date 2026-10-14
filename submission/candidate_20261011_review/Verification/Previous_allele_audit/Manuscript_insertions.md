# Manuscript integration text — allele-direction audit

Prepared against TED_Research_Contribution_20261003/Manuscript.md. Existing reference numbers 27 (Colobran 2011) and 28 (Marín-Sánchez 2019) remain correct. These are ready-to-integrate additions/replacements; the delivered manuscript DOCX/PDF has not been silently replaced. Merge the audit tables into the single detailed supplementary dataset at the next consolidated build, preserving Supplementary Tables S1–S3 and Figures S1–S2. A1–A4 here are audit-file identifiers, not new manuscript table numbers.

## Methods addition (post hoc audit, not a new primary hypothesis)

In a post hoc allele-direction audit, original eQTLGen and disease summary statistics were re-extracted for the TSHR instrument rs179252 and the thymic-expression-associated variant rs179247. Disease effects were aligned to the blood expression-increasing allele. Allele-labelled genotype correlations in 1000 Genomes European and East Asian reference samples were checked by direct genotype decoding and independently using PLINK. This analysis assessed directional consistency with published thymus findings, rather than testing thymic mediation.

## Results addition

Allele verification confirmed that rs179252-G increased blood TSHR expression and was associated with lower disease odds in all five outcomes. It was positively correlated with rs179247-G in European (r = 0.978; r² = 0.957; n = 503) and East Asian (r = 0.996; r² = 0.992; n = 504) reference samples. The latter allele also increased blood TSHR expression and corresponds to the protective allele preferentially transcribed in the published thymus study [28]. Recalculated MR estimates reproduced the existing results. These linked variants do not constitute independent genetic replications.

## Replacement for the current Discussion paragraph beginning “The lower-odds direction…”

The lower-odds direction for TSHR need not contradict the pathogenic role of receptor-stimulating antibodies. Blood transcript abundance, tissue-specific expression and receptor activation are different exposures. Allele verification showed that the blood expression-increasing rs179252-G allele was linked to rs179247-G in both European and East Asian reference samples. The latter is the GD-protective allele associated with greater thymic transcription in published studies [27, 28] and also increased blood expression in eQTLGen. This directional consistency strengthens the plausibility of a central-tolerance explanation, but does not demonstrate that thymic expression mediates the observed disease association. The variants remain correlated markers of a locus rather than proven causal regulators. The marked departure of the MR OR from 1 also reflects division by a small reconstructed expression coefficient. It should not be translated into the effect of experimentally increasing TSHR expression.

## Supplementary Methods addition

The audit used the original eQTLGen 2019-12-11 full cis-eQTL file and the same five disease GWAS sources as the main analysis. The GRCh37 positions of rs179252 and rs179247 were 14:81435985 and 14:81432546; their corresponding GRCh38 positions in disease data were 14:80969641 and 14:80966202. Both allele pairs were non-palindromic and matched between sources. UKB reported the opposite effect allele at both SNPs; its beta was reversed after the existing linear-to-log-odds conversion. No sign reversal was needed for BBJ or FinnGen. Reconstructed exposure coefficients and first-order Wald estimates were independently recalculated using the original reference-frequency scaling. Source GWAS P values and recalculated normal-approximation MR P values were retained separately.

European (n = 503) and East Asian (n = 504) reference samples had no missing genotypes at either variant. Correlations were explicitly oriented to G-allele counts at both SNPs, yielding r = 0.978322812 and 0.995803919, respectively. The corresponding genotype-based r² values were 0.957115524 and 0.991625445. PLINK's default A1 correlation in Europeans was negative because its A1 for rs179252 was T; reorientation to G resolved this apparent sign difference. Independent PLINK two-locus maximum-likelihood haplotype estimates supported G–G and T–A coupling. These are phase estimates from unphased reference genotypes, not directly measured phased haplotypes in the study participants. Haplotype-based and genotype-based r² were not interchanged. This audit does not resolve multiple causal signals, identify a causal variant, or establish thymic or orbital mediation.

## Integration checklist for the next manuscript build

- Add Methods/Results text and replace the relevant Discussion paragraph.
- Replace the earlier blanket statement that the matching direction is unestablished with the specific limitation that direction is consistent but mediation and causal identity are untested.
- Add A1–A4 audit tables as worksheets/sections in the existing single detailed supplementary dataset, with an updated data dictionary; do not proliferate separate numbered supplementary tables.
- No numeric changes to main Tables 1–3 or Figures 1–3 are required by this audit.
- Label the tissue protocol as a planned extension until the 50-cell grid is populated; do not write that replication has already been performed.
