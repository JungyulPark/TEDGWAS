**Running head:** TSHR and IGF1R genetics in Graves disease

# TSHR and IGF1R expression in Graves disease: Mendelian randomization and colocalization

Jungyul Park¹, Min-Seon Kim², Kyung-Hwa Shin³⁎, Suk-Woo Yang¹⁎

¹ Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
² Department of Ophthalmology, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
³ Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea  

\* These authors contributed equally as corresponding authors.

**Corresponding author:** Suk-Woo Yang, MD, PhD, Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea. E-mail: yswoph@catholic.ac.kr; Tel: +82-2-2258-2847.  
**Co-corresponding author:** Kyung-Hwa Shin, MD, Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea.

## Abstract

**Objective:** Graves disease (GD) is an autoimmune cause of hyperthyroidism; thyroid eye disease (TED) is an associated inflammatory orbital disorder. We compared genetic evidence linking blood *TSHR* and *IGF1R* expression to disease susceptibility. **Methods:** We screened 2,544 druggable genes using Wald-ratio or inverse-variance-weighted Mendelian randomization and colocalization, testing for shared expression–disease variants. Outcomes were Biobank Japan (BBJ) GD, UK Biobank (UKB) hyperthyroidism and FinnGen Graves ophthalmopathy. In-house orbital tissue RNA sequencing (four TED patients, one control) was summarized by transcript levels and fold changes without significance testing. **Results:** Higher genetically proxied *TSHR* expression was associated with lower GD odds; *IGF1R* findings were less conclusive. Thirteen genes met the multiple-testing threshold. Per approximately one standard deviation higher blood *TSHR* expression, the BBJ odds ratio (OR) was 0.12 (95% confidence interval [CI] 0.07–0.21; *P* = 1.09×10⁻¹⁴), with matching directions in UKB and FinnGen. Shared-variant support was strong in BBJ and FinnGen, but not UKB; BBJ support fell below threshold when shared variants were assumed less likely. *IGF1R* did not meet the discovery threshold or show strong shared-variant support. Associations reached *P* < 0.05 in BBJ (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB (*P* = 0.0117), but were inconclusive in FinnGen (*P* = 0.182). No additional gene had strong shared-variant support in BBJ and FinnGen. Using cohort allele frequencies preserved these findings. **Conclusions:** Evidence linking genetically proxied blood expression to GD susceptibility was stronger for *TSHR* than *IGF1R*. This neither establishes a TED-specific effect within GD nor predicts treatment response.

**Keywords:** Graves disease; Thyroid eye disease; *TSHR*; *IGF1R*; Mendelian randomization; Colocalization

## Introduction

Thyroid eye disease (TED) can cause proptosis, diplopia and visual impairment in patients with Graves disease (GD) [1, 2]. Teprotumumab, which blocks the insulin-like growth factor 1 receptor (IGF-1R), improves clinical outcomes in active TED [3]. Its therapeutic efficacy raises a related question: does inherited variation influencing *IGF1R* expression contribute to disease susceptibility in the same way as variation at *TSHR*, the established GD autoantigen and susceptibility locus?

TSHR and IGF-1R interact in experimental models of orbital disease [4, 5, 6]. However, a receptor can be an effective treatment target without the genetic evidence for its role in disease susceptibility being equally strong. Treatment response, inherited gene regulation and gene expression in inflamed tissue address different biological questions. Distinguishing them can help clinicians interpret claims arising from genetic studies of therapeutic targets.

Mendelian randomization (MR) uses genetic variants associated with an exposure to study its relationship with disease [7] — because the genetic variants a person inherits are assigned at conception in a way that mimics the random allocation of treatments in a clinical trial, this approach can support cause-and-effect inference without the confounding that affects ordinary observational studies. Colocalization is a complementary analysis that asks whether the exposure and disease associations at a locus are consistent with the same underlying genetic variant [8] — two signals that sit near each other in a region of the genome can look connected simply because they are inherited together, rather than because they reflect one shared causal variant, and colocalization checks which explanation the data actually support. This additional check matters because, without it, an exposure–disease association found by MR could reflect such a linked-but-distinct pair of signals rather than a truly shared one. The two analyses also report different kinds of numbers: MR results appear as odds ratios (ORs) with 95% confidence intervals and *P* values, in the same format as a standard clinical association study, while colocalization results appear instead as a posterior probability that the two signals share one causal variant — a different kind of number, and not a *P* value. Prior MR and multi-omics studies have investigated molecular contributors to GD [9, 10], but the comparison between *TSHR* and *IGF1R* remains clinically relevant.

Here, we used blood gene expression as the exposure and screened genes considered potentially druggable — that is, genes encoding proteins, such as receptors or enzymes, that are plausible targets for an existing or future drug — comparing *TSHR*, *IGF1R* and the known autoimmune gene *CTLA4* by both MR and colocalization. Our primary question was whether genetically proxied expression — a person's blood gene-expression level as predicted from their own inherited genetic variants, used as a stand-in for a measurement that was not taken directly — was associated with GD and supported by a shared expression–disease variant. We then assessed broader hyperthyroidism and a Graves ophthalmopathy outcome. A small in-house orbital RNA-sequencing dataset (four TED samples and one control) is reported alongside these analyses as descriptive context only; with a single control it supports no statistical comparison, and no genetic conclusion rests on it. The study was designed to inform interpretation of genetic target evidence; it did not test a treatment or directly compare GD patients with and without TED.

## Methods

### Study design and data sources

We conducted a two-sample MR study using public genetic association summary statistics and followed STROBE-MR reporting guidance [11]. Blood gene-expression data came from the eQTLGen Consortium, comprising up to 31,684 participants, predominantly of European ancestry [12]; the number contributing to each variant varied widely (median 28,092, range 232–31,684). Disease outcomes were Biobank Japan (BBJ) GD, with 2,809 cases [13]; UK Biobank (UKB) hyperthyroidism, with 3,731 cases [14]; and FinnGen R12 Graves ophthalmopathy, with 858 cases [15] (Table 1).

The defined outcome hierarchy used BBJ for discovery, UKB for comparison across ancestry and a broader phenotype, and FinnGen as a TED-enriched sensitivity outcome. FinnGen cases were compared with population controls. Consequently, that analysis includes susceptibility to GD and cannot isolate susceptibility to eye disease among people who already have GD. We selected *TSHR*, *IGF1R* and *CTLA4* on biological and therapeutic grounds and evaluated them across all outcomes, regardless of discovery significance.

### Genetic instruments and statistical analysis

We considered 4,462 druggable genes [16]. Eligible instruments were genetic variants within 1 Mb of a gene that were associated with its blood expression at *P* < 5×10⁻⁸. We selected approximately independent variants by linkage-disequilibrium clumping in PLINK [17] against the 1000 Genomes Phase 3 European reference panel [18], which matched the ancestry of the exposure data. This left 6,135 instruments for 2,544 genes. The minimum instrument-strength F statistic was 29.7, well above the conventional threshold of 10. Alleles were aligned between expression and disease datasets; ambiguous or unavailable variants were excluded. UKB estimates were converted from their reported linear-model scale to log-odds before analysis. Exposure effect sizes were reconstructed from the eQTLGen Z statistics, per-variant sample sizes and reference-panel allele frequencies [19]; the formulas, selection parameters and exclusions are provided in Supplementary Methods.

We used the Wald ratio for genes with one instrument and inverse-variance weighting for genes with multiple instruments [20]. Results are reported as odds ratios (ORs), 95% confidence intervals (CIs) and two-sided *P* values. ORs are per approximately one standard deviation higher genetically proxied blood expression, because the reconstruction assumes an expression variance of one; the genetic predictor itself is not standardized. The estimates are not per drug dose and should not be interpreted as expected treatment effects. Discovery significance required *P* < 1.965×10⁻⁵, correcting for 2,544 eligible genes. This Bonferroni correction is conservative when gene-level tests are correlated. Associations at *P* < 0.05 in the additional outcomes were considered nominal support only when their direction agreed with discovery. We report all estimates for the three selected genes, including non-significant results. Alternative MR estimators and tests for inconsistent instrument effects are presented in Table S1 [21, 22].

MR interpretation requires that instruments predict expression, are not associated with confounders, and do not affect disease through other pathways. Neither a small *P* value nor colocalization proves these assumptions. Single-instrument results, including *TSHR*, cannot undergo instrument-based tests of pleiotropy. We did not perform a formal genetic directionality test [23]. The reporting checklist records other analyses that were unavailable or not undertaken.

### Shared genetic signals and sensitivity analyses

For the three selected genes and six additional candidates, colocalization assessed whether gene expression and disease shared a genetic variant [8]. We considered a shared-variant posterior probability of at least 0.80 to provide strong support. This probability is a model-based measure of support, not a *P* value. We report it alongside the MR results and tested more conservative assumptions (Figure 3B; Supplementary Data 3) [24]. Colocalization used the regional expression and disease statistics without any significance filter, giving 1,479 to 7,740 variants per gene and outcome. Requiring support in both the discovery outcome and the TED-enriched outcome is a stringent cross-outcome filter: it selects candidates whose shared-variant evidence is not confined to one dataset, and it is not a test of Graves disease association on its own. The combined criterion was not prospectively registered, and the available analysis plan does not establish that this exact rule preceded inspection of results. The model assumes at most one causal variant per trait in the region; full model details and all probabilities are supplied in Supplementary Methods and Supplementary Data 3.

We checked instrument selection at *TSHR* and whether three neighbouring chromosome 16p11.2 hits represented independent findings. We also repeated the analyses using eQTLGen's own allele frequencies. Frequency substitution was assessed first with the same variants, then after repeating allele alignment. MR estimates and detectable effect sizes were compared with the primary analysis in Tables S2–S3, and all colocalization scenarios in Supplementary Data 3. Details of these checks, including their limitations, are in Supplementary Methods.

A small orbital RNA-seq dataset comprising four TED samples and one control was retained solely as descriptive context (Figure S1). It did not contribute to gene selection or genetic conclusions, and no differential-expression *P* values were calculated for this manuscript. Technical replicates were combined at the biological-sample level.

We additionally performed a post hoc leave-one-out analysis for the three selected genes, omitting each retained SNP in turn where at least two instruments were available (Figure S2; Supplementary Data 6). This assessed sensitivity to individual instruments rather than providing independent replication.

Analyses used R 4.3.3, TwoSampleMR 0.7.4, coloc 5.2.3 and PLINK 1.9. Python reproduced the colocalization calculations for verification and frequency sensitivity. The instruments and complete numerical results are supplied in Supplementary Data 1–6. No prospective study registration was available.

## Results

### Discovery across druggable genes

Of 2,544 genes with selected instruments, 2,234 had estimable MR results in BBJ; the corresponding numbers were 2,505 in UKB and 2,480 in FinnGen. The 310 genes without an estimable BBJ result were not treated as negative findings. Thirteen genes met the BBJ discovery threshold: seven were associated with lower odds and six with higher odds of disease (Table 3; Figure 1). They comprised *TSHR*, *CTLA4*, five genes in the major histocompatibility complex (MHC) and six additional candidates, so the discovery signal was not confined to the MHC region. The three genes selected on biological and therapeutic grounds are compared in Table 2 and Figure 2, which aligns each effect estimate with the corresponding shared-variant evidence.

### TSHR

Higher genetically proxied blood *TSHR* expression was associated with lower odds of BBJ GD (OR 0.12, 95% CI 0.07–0.21; *P* = 1.09×10⁻¹⁴). Estimates had the same direction in UKB (OR 0.09, 95% CI 0.06–0.14; *P* = 8.77×10⁻²⁸) and FinnGen (OR 0.10, 95% CI 0.04–0.24; *P* = 2.82×10⁻⁷).

Under the primary prior, expression and disease showed strong shared-variant support in BBJ and FinnGen (0.951 and 0.986), with rs179252 the highest-probability variant conditional on the shared-variant model in both. UKB did not, at 0.226, even though it carried the strongest MR association of the three outcomes; under the single-causal-variant model its posterior favoured distinct variants (PP.H3 = 0.774), which does not establish the number of independent regional signals (Figure 2). A highly significant MR result therefore did not by itself establish a shared genetic signal in every outcome. Substituting cohort allele frequencies did not resolve the discordance.

### IGF1R

Higher genetically proxied *IGF1R* expression was associated with higher odds of BBJ GD (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB hyperthyroidism (OR 1.35, 95% CI 1.07–1.70; *P* = 0.0117). These were nominal associations and did not meet the study-wide discovery threshold. The FinnGen estimate was in the same direction but imprecise (OR 1.41, 95% CI 0.85–2.33; *P* = 0.182).

Shared-variant support was below 0.80 in all three outcomes (BBJ 0.073, UKB 0.404 and FinnGen 0.032). Expression-only model support (PP.H2) was 0.690 in BBJ and 0.623 in FinnGen; this model preference does not establish an absence of marginal disease association. The UKB posterior was divided between an expression-only signal (PP.H2 = 0.400) and a shared variant (PP.H4 = 0.404), with little support for distinct variants (PP.H3 = 0.196), and should be regarded as unresolved rather than negative.

Alternative MR estimators agreed in direction but differed in statistical significance (Table S1): weighted-median estimates were nominally significant in BBJ and UKB, weighted-mode estimates were not, and neither provided nominal support in FinnGen. These methods used the same limited instrument sets and therefore do not constitute independent replication. All 11 *IGF1R* SNP-exclusion estimates retained the higher-odds direction, but the nominal associations were not robust to omission: in BBJ and UKB — the two outcomes that reached *P* < 0.05 — two of the four omissions removed nominal significance in each. Omitting rs2654980 yielded OR 1.49 (95% CI 0.96–2.30; *P* = 0.0751) in BBJ and OR 1.16 (0.78–1.72; *P* = 0.468) in UKB (Figure S2; Supplementary Data 6). rs2654980 also ranked first conditional on the shared-variant model, but because the MR and colocalization analyses reuse the same underlying association data, this is not independent corroboration. Taken together, these results provide weaker genetic support for *IGF1R* than for *TSHR* in BBJ Graves disease and FinnGen Graves ophthalmopathy; they do not show that *IGF1R* has no inherited contribution.

### CTLA4 and additional candidates

The known autoimmune locus *CTLA4* [25] was associated with lower odds of disease in BBJ (OR 0.18, 95% CI 0.11–0.27; *P* = 5.45×10⁻¹⁵), UKB (OR 0.21, 95% CI 0.08–0.57; *P* = 0.0022) and FinnGen (OR 0.17, 95% CI 0.04–0.68; *P* = 0.0118). Shared-variant support was strong in UKB and FinnGen (0.953 and 0.978), but not BBJ (0.201), where the posterior favoured distinct variants (PP.H3 = 0.799; Figure 2). In the two European outcomes, the MR association was largely supported by rs13030124: omitting it left imprecise estimates based on rs148849825 alone (UKB OR 1.02, 95% CI 0.37–2.82; FinnGen OR 1.72, 0.16–18.90). This agreed with the observed heterogeneity and limits interpretation of the two-SNP IVW estimates (Table S1; Figure S2). The screen therefore recovered a known autoimmune association — in BBJ on a single instrument — while *CTLA4* itself failed the combined BBJ-plus-FinnGen colocalization criterion.

None of the six additional candidates met the shared-variant criterion in both BBJ and FinnGen. *TNFSF14* and *IFNGR1* had strong support in BBJ only (0.994 and 0.989); the corresponding FinnGen probabilities were 0.017 and 0.020 (Figure 3A). *TNFSF14* nevertheless retained a direction-consistent UKB MR association (*P* = 0.0056); the absence of colocalization in that outcome was not an absence of association. The three neighbouring chromosome 16p11.2 hits were consistent with one regional signal rather than three independent findings (Supplementary Methods). Their BBJ shared-variant probabilities were below 0.80, and none gained strong support in UKB or FinnGen. *MAPKAPK5* also failed the shared-variant criterion despite passing the MR discovery threshold. Thus, the absence of an additional qualifying candidate resulted from examining the combined evidence, rather than from an absence of significant MR associations. All 13 discovery estimates and their classification are shown in Table 3; full comparisons are in Supplementary Data 2–3.

### Robustness and detectable effects

Using eQTLGen cohort frequencies preserved all 13 discovery hits. *IGF1R* remained nominally associated in BBJ and UKB (*P* = 0.022 and 0.016) and non-significant in FinnGen (*P* = 0.177; Table S2). Frequency substitution moved none of the 81 colocalization comparisons (nine genes × three outcomes × three priors) across the 0.80 threshold, so every classification was unchanged; the largest absolute change in shared-variant probability was 0.002142. The prior assumption was a different matter: under the most conservative prior, *TSHR* support in BBJ fell from 0.951 to 0.661, below the 0.80 threshold, while FinnGen remained above it at 0.875 (Figure 3B; Supplementary Data 3). Stability to the frequency substitution therefore does not imply stability to the prior.

The screen had limited ability to detect smaller associations. At 80% power, only 14.6% of estimable BBJ genes could detect an OR of 1.5 at the discovery threshold; this decreased to 12.5% with cohort frequencies. The corresponding proportions for an OR of 2.0 were 35.6% and 34.2% (Table S3). Gene-specific detection benchmarks are provided in Supplementary Methods. They describe the precision of MR association detection, not the probability of passing the complete colocalization filter, and cannot distinguish an absent effect from limited power; the ORs, 95% CIs and *P* values in Table 2 should guide interpretation of the observed associations. Failure to identify an additional candidate therefore does not exclude moderate or other undetected genetic effects.

## Discussion

This study found stronger genetic evidence linking blood *TSHR* expression to GD than linking *IGF1R* expression to the same outcomes. The *TSHR* association ran in the protective direction — higher genetically proxied expression corresponded to lower odds of disease — consistently across all three datasets, with shared expression–disease signals in BBJ GD and FinnGen Graves ophthalmopathy under the primary prior; in UKB neither gene reached the shared-variant threshold. Because the FinnGen outcome compares Graves ophthalmopathy cases with population controls, that recovery shows the signal extends to a TED-enriched population but does not isolate susceptibility to eye disease among people who already have GD. *IGF1R* showed nominal associations in BBJ and UKB, but its shared-variant evidence was unresolved in every outcome. No additional candidate passed the combined BBJ and FinnGen criterion. The comparison is therefore informative about the pattern and strength of genetic evidence, rather than a ranking of expected drug efficacy.

For clinicians, susceptibility and therapeutic response concern different stages of disease. Teprotumumab showed efficacy in the phase 2 and phase 3 randomized trials of active TED [3, 26]. OPTIC-X provided extension and re-treatment data [27], and a placebo-controlled trial in longstanding, low-activity TED also found improved proptosis [28]. These studies support treatment effects in established disease; none tests inherited blood-expression effects on disease onset. Two observations bear on how that gap might arise: the two receptors are physically scaffolded in experimental systems [5], and expression regulation is tissue-dependent [29]. A model compatible with our findings is that IGF-1R contributes to orbital effector processes even when inherited blood-expression variation provides weaker susceptibility evidence. This is a hypothesis rather than proof of an exclusively downstream role: the nominal BBJ and UKB associations and imprecise FinnGen estimate leave an inherited contribution unresolved.

The protective direction of the *TSHR* estimate is compatible with a central-tolerance hypothesis. Human thymic studies associated the GD-protective rs179247 allele with greater intrathymic *TSHR* transcription, a pattern that could improve immune recognition of TSHR as self [30, 31]. We calculated strong LD between our instrument rs179252 and rs179247 in European and East Asian reference samples (r² = 0.957 and 0.992, respectively; Supplementary Methods), which connects the instrument to a biologically relevant region: if the relevant allele increases thymic antigen expression, the protective association would be biologically coherent. That tissue-specific step remains untested here, since neither LD nor a blood eQTL demonstrates the same regulatory direction in thymus. At rs179252 itself the disease odds ratios were approximately 0.80, 0.77 and 0.78 across BBJ, UKB and FinnGen. The much larger MR estimates follow from the scale rather than from a larger effect: the variant shifts reconstructed expression by only about 0.11 SD, and dividing the per-allele effect by that shift expresses it per standard deviation (Table 2). Genetic target evidence [32] therefore requires tissue and exposure interpretation before translation into an intervention.

The UKB *TSHR* and BBJ *CTLA4* discordances illustrate the limits of the colocalization model as fitted. Their H3-dominant posteriors favour distinct variants under the single-causal-variant assumption; they do not demonstrate multiple independent signals. Multiple causal variants can nevertheless produce misleading single-variant colocalization summaries and are an important alternative explanation [33]. Both UKB and FinnGen are broadly European-ancestry outcomes, so ancestry mismatch alone cannot explain the UKB result. Phenotype definition, signal complexity, sample composition and variant coverage remain possible contributors. The reported top SNP is conditional on the shared-variant model and is not necessarily the lead disease-association SNP. Multi-signal colocalization was not performed, so the regional architecture remains unresolved.

The combined BBJ-plus-FinnGen filter favours evidence that is detectable in both outcomes and may miss valid susceptibility loci. Seven of nine evaluated FinnGen loci had H2-dominant posteriors, consistent with limited outcome evidence relative to expression evidence under this model. However, H2 is not a test proving that marginal disease association is absent. *CTLA4* illustrates a different failure mode: it passed in FinnGen (PP.H4 = 0.978) and failed in BBJ (0.201), where H3 was dominant. Failure of this known autoimmune locus cautions against interpreting filter failure as exclusion of susceptibility, although one comparator cannot quantify the filter's sensitivity. The criterion is therefore a stringent prioritization rule rather than a validated diagnostic test. Significant associations such as *TNFSF14* and *IFNGR1* remain relevant to follow-up even when they fail that rule [9, 10].

The sensitivity analyses improve confidence in specific aspects of the results without removing all uncertainty. Substituting eQTLGen cohort frequencies preserved the discovery-hit set and every tested shared-variant classification. It therefore addresses concern about the reference-frequency approximation used to reconstruct exposure effects. SNP-exclusion analyses further showed that nominal *IGF1R* significance was sensitive to instrument omission and that European *CTLA4* support was concentrated in one SNP. These post hoc checks refine interpretation of the same data; changes in significance reflect both effect estimates and precision and do not establish bias in a particular SNP. Prior sensitivity answers a separate question: *TSHR* support in BBJ fell below the chosen threshold under the most conservative shared-association prior, while FinnGen support remained above it. The BBJ finding should consequently be described as supported under the primary model and sensitive to the prior, rather than uniformly robust. Similarly, the screen's limited ability to detect moderate effects prevents the lack of additional qualifying candidates from excluding smaller inherited effects. Detection benchmarks quantify the design's precision; they do not determine whether a particular unresolved result reflects no effect or inadequate power.

The clinical question remains why some patients with GD develop eye disease. Established clinical risk factors and disease context discussed in the TED literature [1, 2] cannot be quantified from these summary statistics, and our analysis does not apportion genetic versus non-genetic risk. A direct GD-with-TED versus GD-without-TED comparison is needed to isolate eye-disease susceptibility. Similar effects in two population-control GWAS would not by themselves resolve this issue because the phenotypes and participants can overlap. Larger molecular datasets from relevant cell types could test whether the blood-derived signals transfer to thymus or orbital tissue. The present orbital observations—four TED samples and one control—remain descriptive and cannot establish that transfer or distinguish causes from consequences of established disease.

### Limitations

The main limitations concern phenotype definition, exposure relevance and available information. FinnGen compared 858 Graves ophthalmopathy cases with population controls rather than GD patients without eye disease, so its associations may reflect GD susceptibility and cannot isolate a TED-specific effect. UKB hyperthyroidism was broader than GD, making phenotype and ancestry differences partly confounded. The Wald ratio divides the outcome effect by the exposure effect; applying European eQTL weights to BBJ therefore assumes that the relevant expression effects transfer across ancestries. Differences in regulation, linkage patterns and variant availability may violate this assumption, and European blood expression may not represent thymic or orbital regulation. Participant overlap between eQTLGen and the European outcomes could not be quantified. Cohort-frequency substitution addressed one part of exposure reconstruction, but the published frequencies excluded Framingham participants and did not resolve these broader issues. The tissue series had only one control, was not externally replicated and cannot support differential-expression or causal inference.

Additional limitations concern the genetic models and detection range. *TSHR* was represented by one instrument under the European selection reference, precluding instrument-based pleiotropy tests; the other selected genes also had few instruments. Colocalization cannot exclude an effect through a neighbouring gene, and the single-causal-variant assumption may not capture more complex regional signals. The East Asian reference selected two *TSHR* variants at the primary threshold, so single-instrument status is not evidence for a single biological signal in every ancestry. Conclusions also depend on prior assumptions, particularly for *TSHR* in BBJ. Multi-variant colocalization, genetic directionality, formal colocalization power and the probability of passing the combined filter were not evaluated. In particular, the UKB TSHR and BBJ CTLA4 discordances remain unresolved by a multi-signal model. Together with the screen's limited MR detection range, these constraints leave smaller or otherwise undetected genetic effects unresolved.

### Conclusions

Within these limits, *TSHR* showed stronger expression–disease evidence than *IGF1R* in BBJ Graves disease and FinnGen Graves ophthalmopathy. The findings distinguish genetic susceptibility evidence from evidence of therapeutic efficacy, while preserving the uncertainty around *IGF1R* and the prior dependence of the BBJ *TSHR* result. They support further investigation in relevant tissues and a direct comparison of GD patients with and without TED; they do not establish treatment direction or expected clinical benefit.

## Declarations

**Funding.** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Conflict of interest.** The authors declare that they have no conflict of interest.

**Ethics approval.** This study used publicly available, de-identified summary statistics and an institutionally approved orbital tissue dataset. The in-house orbital transcriptomic component was approved by the Institutional Review Board of Pusan National University Hospital (approval number 2104-018-102) and was conducted in accordance with the Declaration of Helsinki.

**Informed consent.** Written informed consent was obtained from all individual participants included in the in-house orbital tissue study.

**Data availability.** Instruments, primary MR and colocalization results, and sensitivity results are supplied as Supplementary Data 1–6. Analysis code is available from the corresponding author on reasonable request. Public summary statistics analyzed in this study are available from their original repositories: blood *cis*-eQTL data from the eQTLGen Consortium, the Biobank Japan Graves disease and UK Biobank hyperthyroidism genome-wide association statistics through the GWAS Catalog, and the FinnGen Release 12 Graves ophthalmopathy statistics from FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions.

**Author contributions.** J.P. conceived and designed the study, performed the analyses, and drafted the manuscript. M.-S.K. contributed to data collection and interpretation. K.-H.S. and S.-W.Y. supervised the study and revised the manuscript.

**Acknowledgements.** We thank the eQTLGen Consortium, Biobank Japan, the UK Biobank, and the FinnGen study and its participants for making their summary statistics publicly available. We acknowledge the GWAS Catalog for hosting and distributing the genome-wide association summary statistics used as outcome data.

## References

1. Bahn RS. Graves' ophthalmopathy. N Engl J Med. 2010;362:726-738. doi:10.1056/nejmra0905750.

2. Smith TJ, Hegedüs L. Graves' Disease. N Engl J Med. 2016;375:1552-1565. doi:10.1056/nejmra1510030.

3. Douglas RS, Kahaly GJ, Patel A, Sile S, Thompson EHZ, Perdok R, Fleming JC, Fowler BT, Marcocci C, Marinò M, et al. Teprotumumab for the Treatment of Active Thyroid Eye Disease. N Engl J Med. 2020;382:341-352. doi:10.1056/nejmoa1910434.

4. Smith TJ. The insulin-like growth factor-I receptor and its role in thyroid-associated ophthalmopathy. Eye (Lond). 2019;33:200-205. doi:10.1038/s41433-018-0265-2.

5. Krieger CC, Boutin A, Jang D, Morgan SJ, Banga JP, Kahaly GJ, Klubo-Gwiezdzinska J, Neumann S, Gershengorn MC. Arrestin-β-1 Physically Scaffolds TSH and IGF1 Receptors to Enable Crosstalk. Endocrinology. 2019;160:1468-1479. doi:10.1210/en.2019-00055.

6. Smith TJ. Potential Roles of CD34+ Fibrocytes Masquerading as Orbital Fibroblasts in Thyroid-Associated Ophthalmopathy. J Clin Endocrinol Metab. 2019;104:581-594. doi:10.1210/jc.2018-01493.

7. Sanderson E, Glymour MM, Holmes MV, Kang H, Morrison J, Munafò MR, Palmer T, Schooling CM, Wallace C, Zhao Q, et al. Mendelian randomization. Nat Rev Methods Primers. 2022;2:6. doi:10.1038/s43586-021-00092-5.

8. Giambartolomei C, Vukcevic D, Schadt EE, Franke L, Hingorani AD, Wallace C, Plagnol V. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014;10:e1004383. doi:10.1371/journal.pgen.1004383.

9. Ji Q, Xu H, Chen H, Chen X, Wang S, Zou J. Pathogenic genes associated with immune-related genes in graves' disease: a multi-omics Mendelian randomization analysis. Sci Rep. 2025;15:37875. doi:10.1038/s41598-025-21754-4.

10. Li Y, Chen L, Lin S, An W, Miao L, Wan M, Zhang B. Inflammatory Proteins Mediate the Effect of Gut Microbiota on Graves' Ophthalmopathy: A Mendelian Randomization Study. Transl Vis Sci Technol. 2025;14:34. doi:10.1167/tvst.14.6.34.

11. Skrivankova VW, Richmond RC, Woolf BAR, Yarmolinsky J, Davies NM, Swanson SA, VanderWeele TJ, Higgins JPT, Timpson NJ, Dimou N, et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: The STROBE-MR Statement. JAMA. 2021;326:1614-1621. doi:10.1001/jama.2021.18236.

12. Võsa U, Claringbould A, Westra HJ, Bonder MJ, Deelen P, Zeng B, Kirsten H, Saha A, Kreuzhuber R, Yazar S, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. Nat Genet. 2021;53:1300-1310. doi:10.1038/s41588-021-00913-z.

13. Sakaue S, Kanai M, Tanigawa Y, Karjalainen J, Kurki M, Koshiba S, Narita A, Konuma T, Yamamoto K, Akiyama M, et al. A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet. 2021;53:1415-1424. doi:10.1038/s41588-021-00931-x.

14. Dönertaş HM, Fabian DK, Valenzuela MF, Partridge L, Thornton JM. Common genetic associations between age-related diseases. Nat Aging. 2021;1:400-412. doi:10.1038/s43587-021-00051-5.

15. Kurki MI, Karjalainen J, Palta P, Sipilä TP, Kristiansson K, Donner KM, Reeve MP, Laivuori H, Aavikko M, Kaunisto MA, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature. 2023;613:508-518. doi:10.1038/s41586-022-05473-8.

16. Finan C, Gaulton A, Kruger FA, Lumbers RT, Shah T, Engmann J, Galver L, Kelley R, Karlsson A, Santos R, et al. The druggable genome and support for target identification and validation in drug development. Sci Transl Med. 2017;9:eaag1166. doi:10.1126/scitranslmed.aag1166.

17. Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015;4:7. doi:10.1186/s13742-015-0047-8.

18. 1000 Genomes Project Consortium, Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, et al. A global reference for human genetic variation. Nature. 2015;526:68-74. doi:10.1038/nature15393.

19. Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, Montgomery GW, Goddard ME, Wray NR, Visscher PM, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016;48:481-487. doi:10.1038/ng.3538.

20. Hemani G, Zheng J, Elsworth B, Wade KH, Haberland V, Baird D, Laurin C, Burgess S, Bowden J, Langdon R, et al. The MR-Base platform supports systematic causal inference across the human phenome. Elife. 2018;7:e34408. doi:10.7554/elife.34408.

21. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. Int J Epidemiol. 2015;44:512-525. doi:10.1093/ije/dyv080.

22. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent Estimation in Mendelian Randomization with Some Invalid Instruments Using a Weighted Median Estimator. Genet Epidemiol. 2016;40:304-314. doi:10.1002/gepi.21965.

23. Hemani G, Tilling K, Davey Smith G. Orienting the causal relationship between imprecisely measured traits using GWAS summary data. PLoS Genet. 2017;13:e1007081. doi:10.1371/journal.pgen.1007081.

24. Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLoS Genet. 2020;16:e1008720. doi:10.1371/journal.pgen.1008720.

25. Ueda H, Howson JM, Esposito L, Heward J, Snook H, Chamberlain G, Rainbow DB, Hunter KM, Smith AN, Di Genova G, et al. Association of the T-cell regulatory gene CTLA4 with susceptibility to autoimmune disease. Nature. 2003;423:506-511. doi:10.1038/nature01621.

26. Smith TJ, Kahaly GJ, Ezra DG, Fleming JC, Dailey RA, Tang RA, Harris GJ, Antonelli A, Salvi M, Goldberg RA, et al. Teprotumumab for Thyroid-Associated Ophthalmopathy. N Engl J Med. 2017;376:1748-1761. doi:10.1056/nejmoa1614949.

27. Douglas RS, Kahaly GJ, Ugradar S, Elflein H, Ponto KA, Fowler BT, Dailey R, Harris GJ, Schiffman J, Tang R, et al. Teprotumumab Efficacy, Safety, and Durability in Longer-Duration Thyroid Eye Disease and Re-treatment: OPTIC-X Study. Ophthalmology. 2022;129:438-449. doi:10.1016/j.ophtha.2021.10.017.

28. Douglas RS, Couch S, Wester ST, Fowler BT, Liu CY, Subramanian PS, Tang R, Nguyen QT, Maamari RN, Ugradar S, et al. Efficacy and Safety of Teprotumumab in Patients With Thyroid Eye Disease of Long Duration and Low Disease Activity. J Clin Endocrinol Metab. 2024;109:25-35. doi:10.1210/clinem/dgad637.

29. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369:1318-1330. doi:10.1126/science.aaz1776.

30. Colobran R, Armengol Mdel P, Faner R, Gärtner M, Tykocinski LO, Lucas A, Ruiz M, Juan M, Kyewski B, Pujol-Borrell R. Association of an SNP with intrathymic transcription of TSHR and Graves' disease: a role for defective thymic tolerance. Hum Mol Genet. 2011;20:3415-3423. doi:10.1093/hmg/ddr247.

31. Marín-Sánchez A, Álvarez-Sierra D, González O, Lucas-Martin A, Sellés-Sánchez A, Rudilla F, Enrich E, Colobran R, Pujol-Borrell R. Regulation of TSHR Expression in the Thyroid and Thymus May Contribute to TSHR Tolerance Failure in Graves' Disease Patients via Two Distinct Mechanisms. Front Immunol. 2019;10:1695. doi:10.3389/fimmu.2019.01695.

32. Nelson MR, Tipney H, Painter JL, Shen J, Nicoletti P, Shen Y, Floratos A, Sham PC, Li MJ, Wang J, et al. The support of human genetic evidence for approved drug indications. Nat Genet. 2015;47:856-860. doi:10.1038/ng.3314.

33. Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genet. 2021;17:e1009440. doi:10.1371/journal.pgen.1009440.

## Figure Legends

**Figure 1. Study design and the druggable-gene association screen.** (A) Gene selection and discovery flow. Genes with independent blood-expression instruments were evaluated across the outcome hierarchy; genes without estimable BBJ results were not treated as negative. (B) Gene-level MR associations in BBJ Graves disease, plotted by the GRCh37 position of each gene's strongest selected expression instrument. Each point represents one gene, not a SNP association. The horizontal line marks the discovery threshold of *P* < 0.05/2,544; all 2,234 estimable genes are plotted and all 13 discovery genes are labelled. *TSHR*, *IGF1R* and *CTLA4* are highlighted as the three biologically selected genes. Their selection for comparison did not depend on discovery significance.

**Figure 2. Association estimates and regional genetic evidence for the three selected genes.** (A) Primary ORs and 95% CIs, with corresponding two-sided *P* values. Bold *P* values and asterisks (\*) indicate nominal *P* < 0.05; daggers (†) additionally mark BBJ discovery significance at *P* < 0.05/2,544. The OR axis is logarithmic, with a reference line at one. ORs refer to reconstructed genetically proxied blood expression, not treatment effects. All nine comparisons are included. (B) Regional colocalization probabilities under the primary prior, aligned with the same gene–outcome rows. Stacked bars show expression-only support (H2), support for distinct variants (H3), shared-variant support (H4), and residual H0/H1 support. Exact displayed H4 probabilities are given at right. These probabilities are not *P* values; ≥0.80 was the threshold for strong shared-variant support. The figure displays the discordance between association strength and shared-variant evidence without treating either as proof of a therapeutic effect.

**Figure 3. Shared-variant evidence across candidates and sensitivity to prior assumptions.** (A) Primary-prior PP.H4 values for all nine evaluated genes and three outcomes. Black cell borders indicate strong shared-variant support (PP.H4 ≥ 0.80), not *P*-value significance; probabilities below 0.001 are labelled <0.001. A gene passes the combined criterion only when both BBJ and FinnGen meet that threshold; UKB is an additional comparison. *TSHR*, *IGF1R* and *CTLA4* were biologically selected, while the other six genes were additional discovery candidates. (B) Prior sensitivity for the three selected genes. D, I and C denote default, intermediate and conservative shared-association priors (p12 = 10⁻⁵, 5×10⁻⁶ and 10⁻⁶), with p1 = p2 = 10⁻⁴ throughout. The dashed line marks 0.80. Prior changes can alter classification even when frequency substitution does not; Tables S2–S3 report the separate cohort-frequency sensitivity. Full H0–H4 probabilities and all 81 gene–outcome–prior settings are in Supplementary Data 3.

## Tables

**Table 1. Data sources and roles in the study**

| Dataset | Ancestry | Cases | Controls | Role |
|---|---|---|---|---|
| eQTLGen blood expression [12] | Predominantly European | — | — | Exposure; up to 31,684 participants |
| BBJ Graves disease [13] | East Asian | 2,809 | 172,656 | Discovery |
| UKB hyperthyroidism [14] | European | 3,731 | 480,867 | Broader-phenotype comparison |
| FinnGen R12 Graves ophthalmopathy [15] | European | 858 | 499,490 | TED-enriched sensitivity |

GWAS accessions are GCST90018627 (BBJ), GCST90038636 (UKB) and GRAVES_OPHT in FinnGen R12. eQTL sample size varies by variant. FinnGen uses population controls, so it cannot distinguish TED susceptibility from GD susceptibility. The small descriptive orbital dataset is detailed only in Supplementary Methods and Figure S1.

**Table 2. Primary associations for the three selected genes**

| Gene | Outcome | Instruments | OR (95% CI) | *P* value | Expression only (PP.H2) | Distinct variants (PP.H3) | Shared variant (PP.H4) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| *TSHR* | BBJ Graves disease | 1 | 0.12 (0.07–0.21) | 1.09×10⁻¹⁴ | <0.001 | 0.049 | 0.951 |
| *TSHR* | UKB hyperthyroidism | 1 | 0.09 (0.06–0.14) | 8.77×10⁻²⁸ | <0.001 | 0.774 | 0.226 |
| *TSHR* | FinnGen Graves ophthalmopathy | 1 | 0.10 (0.04–0.24) | 2.82×10⁻⁷ | <0.001 | 0.014 | 0.986 |
| *IGF1R* | BBJ Graves disease | 4 | 1.56 (1.07–2.28) | 0.0212 | 0.690 | 0.236 | 0.073 |
| *IGF1R* | UKB hyperthyroidism | 4 | 1.35 (1.07–1.70) | 0.0117 | 0.400 | 0.196 | 0.404 |
| *IGF1R* | FinnGen Graves ophthalmopathy | 3 | 1.41 (0.85–2.33) | 0.182 | 0.623 | 0.346 | 0.032 |
| *CTLA4* | BBJ Graves disease | 1 | 0.18 (0.11–0.27) | 5.45×10⁻¹⁵ | <0.001 | 0.799 | 0.201 |
| *CTLA4* | UKB hyperthyroidism | 2 | 0.21 (0.08–0.57) | 0.0022 | <0.001 | 0.047 | 0.953 |
| *CTLA4* | FinnGen Graves ophthalmopathy | 2 | 0.17 (0.04–0.68) | 0.0118 | <0.001 | 0.022 | 0.978 |

ORs are per approximately one SD higher reconstructed, genetically proxied blood expression. Wald ratio was used for one instrument and IVW for multiple instruments. *P* values are two-sided and shown to three significant digits; BBJ discovery significance requires *P* < 1.965×10⁻⁵, while *P* < 0.05 in other outcomes is nominal. H2, H3 and H4 denote expression only, distinct variants and a shared variant under the single-causal-variant model; these posterior probabilities are not *P* values. H0/H1 account for any remaining probability. H4 ≥ 0.80 denotes strong shared-variant support. H2 dominance does not prove absence of marginal disease association. For rs179252, the G-allele exposure coefficient was 0.10546; aligned disease log-ORs were −0.22103 (BBJ), −0.25691 (UKB, after scale conversion) and −0.24587 (FinnGen), corresponding to variant-level ORs of 0.802, 0.773 and 0.782. These variant effects and the scaled MR ORs answer different questions. European *CTLA4* MR support was concentrated in rs13030124; the other instrument alone gave imprecise estimates spanning OR = 1 (Figure S2; Supplementary Data 6). Full probabilities are in Supplementary Data 3; sensitivity results are in Tables S1–S3 and Figure S2.

**Table 3. All thirteen discovery associations**

| Gene | Instruments | BBJ OR (95% CI) | *P* value | Interpretation of the locus |
|---|---|---|---|---|
| *HLA-A* | 4 | 1.98 (1.73–2.27) | 2.6×10⁻²³ | MHC region |
| *HLA-DQA2* | 1 | 2.39 (1.95–2.93) | 8.3×10⁻¹⁷ | MHC region |
| *CTLA4* | 1 | 0.18 (0.11–0.27) | 5.45×10⁻¹⁵ | Known autoimmune locus |
| *TSHR* | 1 | 0.12 (0.07–0.21) | 1.09×10⁻¹⁴ | Known GD locus |
| *C4A* | 1 | 0.45 (0.35–0.57) | 2.23×10⁻¹⁰ | MHC region |
| *HSD3B7* | 1 | 0.29 (0.18–0.47) | 2.04×10⁻⁷ | Shared regional signal at 16p11.2 |
| *TUBB* | 2 | 1.49 (1.28–1.74) | 3.14×10⁻⁷ | MHC region |
| *VKORC1* | 1 | 0.13 (0.06–0.29) | 6.35×10⁻⁷ | Shared regional signal at 16p11.2 |
| *TNFSF14* | 1 | 0.63 (0.52–0.76) | 1.48×10⁻⁶ | Additional candidate |
| *PRSS36* | 1 | 6.00 (2.87–12.57) | 2.02×10⁻⁶ | Shared regional signal at 16p11.2 |
| *MAPKAPK5* | 4 | 0.977 (0.967–0.987) | 5.28×10⁻⁶ | Additional candidate |
| *PSMB8* | 1 | 1.25 (1.13–1.38) | 6.86×10⁻⁶ | MHC region |
| *IFNGR1* | 1 | 2.10 (1.51–2.91) | 9.41×10⁻⁶ | Additional candidate |

All genes met the same multiple-testing threshold and are listed by discovery *P* value. The *MAPKAPK5* estimate is shown to three decimals to preserve its small effect size. The five MHC genes were not treated as new independent targets. The six additional candidates were TNFSF14, IFNGR1, MAPKAPK5, HSD3B7, VKORC1 and PRSS36; none showed strong shared-variant support in both BBJ and FinnGen. Full outcome comparisons and model probabilities remain in Supplementary Data 2–3.

## Supplementary Material

## Supplementary Methods

### Expression instruments and allele alignment

Cis-eQTL instruments met *P* < 5×10⁻⁸ within ±1 Mb of the gene and were clumped at r² < 0.001 within 10 Mb using 1000 Genomes Phase 3 Europeans (n = 503). Exposure SE = 1/√[2f(1−f)(N + Z²)] and β = Z × SE, where f is allele frequency and N is the variant-specific sample size [19]. TwoSampleMR harmonise_data (action = 2) aligned alleles and excluded palindromic variants with MAF > 0.42; missing estimates were not imputed. BBJ/FinnGen retained their reported log-odds scales. UKB linear-model β and SE were divided by μ(1−μ), with μ = 3,731/484,598. Source GWAS covariate adjustments were retained; participant-level re-adjustment was unavailable.

### MR sensitivity and regional linkage

Single instruments used the Wald ratio; multiple instruments used multiplicative random-effects IVW with under-dispersion correction. Weighted-median/mode estimates, the MR-Egger intercept and Cochran's Q were assessed where estimable (Table S1; Supplementary Data 5) [21, 22]. Few instruments limit pleiotropy testing. Steiger directionality [23] and bidirectional MR were not performed. Post hoc leave-one-out analysis used the original reference-frequency harmonized sets, applying IVW when at least two SNPs remained and the Wald ratio for one. All 15 omissions (11 IGF1R, four CTLA4) and five comparators are in Figure S2 and Supplementary Data 6. TSHR in all outcomes and CTLA4 in BBJ had one original SNP and were ineligible. An independent base-R implementation reproduced all nine full-set estimates and verified the omissions. Two-sided normal tests and 95% CIs retained the primary SE conventions. These correlated checks are not independent replication. This analysis was not repeated under substituted eQTLGen frequencies.

TSHR clumping combined eQTL thresholds of 5×10⁻⁸/5×10⁻⁶, r² thresholds of 0.001/0.01/0.1 and European (n = 503)/East Asian (n = 504) panels. rs179252 remained an index variant in all 12 combinations. Primary settings selected one European and two East Asian instruments; relaxed settings selected up to 12. Post-review PLINK r² between rs179252 and rs179247 was 0.957116/0.991625 in Europeans/East Asians; for rs179252–rs12101255 it was 0.459227/0.855574. All three variants had complete genotypes. These reference-panel checks support instrument retention, not conditional association, a single-signal model or shared tissue-specific regulation.

At chromosome 16p11.2, rs4889606, rs34649473 and rs78924645 span approximately 143 kb (East Asian pairwise r² = 0.854, 0.761 and 0.863). GCTA-COJO selection at *P* < 5×10⁻⁶ retained rs8050588 (*P* = 1.15×10⁻⁸); conditioning gave *P* = 0.928, 0.839 and 0.806 for the three variants, respectively. This supports regional dependence without identifying a causal gene.

### Colocalization and prior assumptions

coloc.abf used all available gene-specific cis-eQTL statistics without significance filtering, matched to outcomes by rsID and allele pair; GRCh37 expression and GRCh38 outcome coordinates were not joined directly. Duplicate rsIDs and allele mismatches were excluded. Both datasets supplied β and variance: dataset 1 was the outcome, dataset 2 expression (sdY = 1). Outcome frequencies were not used in the Bayes factors. The top_snp field identifies the highest-probability variant conditional on H4, not necessarily the GWAS lead SNP. H0–H4 denote neither trait, outcome only, expression only, distinct variants and a shared variant. Priors were p1 = p2 = 10⁻⁴ and p12 = 10⁻⁵, with 5×10⁻⁶ and 10⁻⁶ sensitivity settings [24]. Figure 3B displays selected genes; Supplementary Data 3 preserves all nine genes × three outcomes × three priors. Low H4 does not establish absence of disease association. Multi-variant colocalization and formal colocalization power were not evaluated.

### Cohort-frequency sensitivity and detection limits

The eQTLGen frequency release covers 26,609 participants, excluding Framingham. All 50,245 requested SNPs, including 6,135 instrument records, had valid frequencies and matching GRCh37 positions. Assessed-allele frequencies replaced reference frequencies in β/SE reconstruction while Z and N were fixed (maximum |ΔZ| = 5.54×10⁻¹³). Identical-variant comparisons preceded repeat action = 2 harmonization. Of 27 changed inclusion decisions, 24 involved palindromic frequencies and three restored missing reference frequencies. No retained outcome estimate changed direction. Estimable genes became 2,232/2,506/2,481 in BBJ/UKB/FinnGen; the discovery denominator remained 2,544 (Table S2; Supplementary Data 2). Across 81 colocalization settings, maximum |ΔH4| was 0.000245 for identical variants and 0.002142 when newly covered variants were included; none crossed 0.80. Primary MR was independently checked. Substitution reran Wald/IVW, colocalization and detection limits, but not all alternative MR estimators.

The approximate 80% minimum detectable absolute log-OR was (z₁₋α/₂ + z₀.₈₀) × SE, with α = 0.05/2,544 in BBJ and 0.05 in UKB/FinnGen; exponentiation gives the detectable OR. Percentages use each scenario's estimable-gene denominator (Table S3; Supplementary Data 4). These describe MR association detection using observed SEs, not colocalization or combined-filter power. Gene-specific limits were 1.39 (TSHR) and 0.99 (IGF1R) in BBJ (α = 0.05/2,544), and 0.33/0.72 for IGF1R in UKB/FinnGen, respectively (α = 0.05). These benchmarks describe precision; comparing them with observed effects cannot distinguish an absent effect from limited power.

### Descriptive orbital expression

Orbital RNA-seq included four TED patients and one control. Technical replicates were combined by biological sample. Figure S1 shows log2[(mean TED TPM + 0.01)/(control TPM + 0.01)]. One control cannot establish control-group biological variance; no differential-expression *P* value or confirmatory test is reported. This unreplicated descriptive dataset contributes no inferential evidence to the genetic comparison.

**Table S1. Sensitivity tests for genes with multiple instruments**

| Gene | Outcome | Instruments | Primary *P* | Median *P* | Mode *P* | Egger intercept *P* | Heterogeneity *P* |
|---|---|---|---|---|---|---|---|
| *IGF1R* | BBJ Graves disease | 4 | 0.0212 | 0.0288 | 0.512 | 0.663 | 0.902 |
| *IGF1R* | UKB hyperthyroidism | 4 | 0.0117 | 0.00244 | 0.367 | 0.562 | 0.25 |
| *IGF1R* | FinnGen Graves ophthalmopathy | 3 | 0.182 | 0.124 | 0.633 | 0.433 | 0.309 |
| *CTLA4* | UKB hyperthyroidism | 2 | 0.0022 | NA | NA | NA | 0.00131 |
| *CTLA4* | FinnGen Graves ophthalmopathy | 2 | 0.0118 | NA | NA | NA | 0.0486 |

Median and mode denote weighted-median and weighted-mode MR. The Egger intercept assesses directional pleiotropy; heterogeneity uses Cochran’s Q. NA means not estimable, not a negative test. TSHR and CTLA4 in BBJ each used one instrument and are omitted from this diagnostic table. IGF1R estimates agreed in direction across estimators, but significance differed. Few instruments limit all these sensitivity checks. These are the original reference-frequency diagnostics; frequency substitution repeated primary Wald/IVW only. Supplementary Data 5 supplies all estimator coefficients, standard errors, ORs and 95% CIs alongside their original *P* values. Its intervals use the same reference distribution as each test: normal for Wald/IVW/weighted median and Student t for MR-Egger (instruments − 2 degrees of freedom) and weighted mode (instruments − 1).

**Table S2. Associations before and after cohort-frequency substitution**

| Gene | Outcome | Reference OR (95% CI) | Reference *P* | eQTLGen OR (95% CI) | eQTLGen *P* |
|---|---|---|---|---|---|
| *TSHR* | BBJ Graves disease | 0.12 (0.07–0.21) | 1.09×10⁻¹⁴ | 0.12 (0.07–0.21) | 1.09×10⁻¹⁴ |
| *TSHR* | UKB hyperthyroidism | 0.09 (0.06–0.14) | 8.77×10⁻²⁸ | 0.09 (0.06–0.14) | 8.77×10⁻²⁸ |
| *TSHR* | FinnGen Graves ophthalmopathy | 0.10 (0.04–0.24) | 2.82×10⁻⁷ | 0.10 (0.04–0.24) | 2.82×10⁻⁷ |
| *IGF1R* | BBJ Graves disease | 1.56 (1.07–2.28) | 0.0212 | 1.58 (1.07–2.34) | 0.022 |
| *IGF1R* | UKB hyperthyroidism | 1.35 (1.07–1.70) | 0.0117 | 1.37 (1.06–1.78) | 0.0157 |
| *IGF1R* | FinnGen Graves ophthalmopathy | 1.41 (0.85–2.33) | 0.182 | 1.44 (0.85–2.47) | 0.177 |
| *CTLA4* | BBJ Graves disease | 0.18 (0.11–0.27) | 5.45×10⁻¹⁵ | 0.18 (0.12–0.27) | 5.45×10⁻¹⁵ |
| *CTLA4* | UKB hyperthyroidism | 0.21 (0.08–0.57) | 0.0022 | 0.20 (0.08–0.49) | 3.86×10⁻⁴ |
| *CTLA4* | FinnGen Graves ophthalmopathy | 0.17 (0.04–0.68) | 0.0118 | 0.17 (0.05–0.57) | 0.00458 |

Reference estimates use 1000 Genomes European frequencies. eQTLGen estimates use cohort frequencies followed by repeat harmonization. For these nine comparisons, estimates using identical variants and repeat harmonization coincide; instrument counts are unchanged from Table 2. All 13 discovery genes remain significant. Results for every gene and each scenario are in Supplementary Data 2.

**Table S3. Association detection limits before and after frequency substitution**

| Outcome | Analysis | Estimable genes | Detectable OR median (IQR) | Detect OR 1.5 (%) | Detect OR 2.0 (%) |
|---|---|---|---|---|---|
| BBJ Graves disease | Reference | 2,234 | 2.55 (1.71–5.18) | 14.6 | 35.6 |
| BBJ Graves disease | Cohort frequencies | 2,232 | 2.61 (1.76–5.14) | 12.5 | 34.2 |
| UKB hyperthyroidism | Reference | 2,505 | 1.46 (1.26–1.91) | 53.7 | 77.4 |
| UKB hyperthyroidism | Cohort frequencies | 2,506 | 1.46 (1.26–1.94) | 53.1 | 76.4 |
| FinnGen Graves ophthalmopathy | Reference | 2,480 | 2.12 (1.58–3.81) | 19.4 | 45.5 |
| FinnGen Graves ophthalmopathy | Cohort frequencies | 2,481 | 2.17 (1.60–3.87) | 17.6 | 44.5 |

Percentages are the proportions of estimable genes with at least 80% power to detect the specified OR under the normal approximation. Alpha is 0.05/2,544 in BBJ and 0.05 in UKB/FinnGen. OR is per reconstructed expression unit. These calculations assess MR association detection, not colocalization or the combined filter. Full precision and identical-variant comparisons are retained in Supplementary Data 4.

**Figure S1. Descriptive orbital expression.** Transcript abundance for TSHR, IGF1R and CTLA4 in four TED samples and one control. Points are biological samples; the horizontal line is the TED mean. Log2 fold changes use a 0.01-TPM pseudocount. Axes differ by gene. With one control, group variance cannot be estimated reliably; no differential-expression test or *P* value is reported. These observations did not determine the genetic conclusions.


**Figure S2. Post hoc leave-one-out sensitivity of the selected-gene MR estimates.** Open circles show the 15 SNP-exclusion estimates; filled diamonds show the five corresponding all-SNP estimates. The SNP named on each row was omitted. Horizontal lines represent 95% CIs, and the dashed line marks OR = 1 on a logarithmic axis. ORs, CIs and two-sided *P* values are displayed for every row; bold *P* values denote nominal *P* < 0.05 for descriptive comparison. Original reference-frequency harmonized sets were used. If one SNP remained, its Wald ratio was calculated; otherwise IVW was used. TSHR in all outcomes and CTLA4 in BBJ are not plotted because a single original instrument precluded omission. These correlated post hoc analyses assess sensitivity to the included SNPs and do not constitute independent replication.
