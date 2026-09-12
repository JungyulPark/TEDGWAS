**Running head:** TSHR and IGF1R genetics in Graves disease

# Genetic evidence for TSHR and IGF1R in Graves disease and thyroid eye disease

Jungyul Park¹, Min-Seon Kim², Kyung-Hwa Shin³⁎, Suk-Woo Yang¹⁎

¹ Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
² Department of Ophthalmology, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
³ Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea  

\* These authors contributed equally as corresponding authors.

**Corresponding author:** Suk-Woo Yang, MD, PhD, Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea. E-mail: yswoph@catholic.ac.kr; Tel: +82-2-2258-2847.  
**Co-corresponding author:** Kyung-Hwa Shin, MD, Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea.

## Abstract

**Objective:** To compare genetic evidence linking blood expression of *TSHR* and *IGF1R* to Graves disease (GD) and thyroid eye disease (TED). **Methods:** We used genetic variants associated with blood gene expression to study 2,544 druggable genes by Mendelian randomization. Biobank Japan (BBJ) GD was the discovery outcome; UK Biobank (UKB) hyperthyroidism and FinnGen Graves ophthalmopathy provided additional comparisons. We reported odds ratios (ORs), 95% confidence intervals (CIs) and *P* values, and tested whether expression and disease associations shared a genetic variant. **Results:** Thirteen genes met the multiple-testing threshold. Higher genetically proxied *TSHR* expression was associated with lower odds of BBJ GD (OR 0.12, 95% CI 0.07–0.21; *P* = 1.09×10⁻¹⁴), with the same direction in UKB and FinnGen. Shared-variant support was strong in BBJ and FinnGen, but not UKB; the BBJ result was sensitive to a more conservative prior. *IGF1R* associations were nominal in BBJ (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB (*P* = 0.0117), and inconclusive in FinnGen (*P* = 0.182); shared-variant support remained below threshold. No additional candidate met the combined BBJ and FinnGen shared-variant criterion. Cohort-frequency sensitivity analyses preserved these findings. **Conclusions:** Genetic evidence was stronger for *TSHR* than *IGF1R* in this analysis. These results do not exclude an inherited *IGF1R* contribution or alter evidence supporting IGF-1R blockade. FinnGen comparisons with population controls cannot establish an effect specific to TED within GD.

**Keywords:** Graves disease; Thyroid eye disease; *TSHR*; *IGF1R*; Mendelian randomization; Colocalization

## Introduction

Thyroid eye disease (TED) can cause proptosis, diplopia and visual impairment in patients with Graves disease (GD) [1, 2]. Teprotumumab, which blocks the insulin-like growth factor 1 receptor (IGF-1R), improves clinical outcomes in active TED [3]. Its therapeutic efficacy raises a related question: does inherited variation influencing *IGF1R* expression contribute to disease susceptibility in the same way as variation at *TSHR*, the established GD autoantigen and susceptibility locus?

TSHR and IGF-1R interact in experimental models of orbital disease [4, 5, 6]. However, a receptor can be an effective treatment target without the genetic evidence for its role in disease susceptibility being equally strong. Treatment response, inherited gene regulation and gene expression in inflamed tissue address different biological questions. Distinguishing them can help clinicians interpret claims arising from genetic studies of therapeutic targets.

Mendelian randomization (MR) uses genetic variants associated with an exposure to study its relationship with disease [7]. Here, the exposure was blood gene expression. We also used colocalization, an analysis that asks whether the expression and disease associations are consistent with the same underlying genetic variant [8]. This additional check matters because neighbouring variants can be inherited together: an expression–disease association may otherwise reflect two linked signals rather than a shared signal. Prior MR and multi-omics studies have investigated molecular contributors to GD [9, 10], but the comparison between *TSHR* and *IGF1R* remains clinically relevant.

We screened genes considered potentially druggable and compared *TSHR*, *IGF1R* and the known autoimmune gene *CTLA4*. Our primary question was whether genetically proxied expression was associated with GD and supported by a shared expression–disease variant. We then assessed broader hyperthyroidism and a Graves ophthalmopathy outcome. The study was designed to inform interpretation of genetic target evidence; it did not test a treatment or directly compare GD patients with and without TED.

## Methods

### Study design and data sources

We conducted a two-sample MR study using public genetic association summary statistics and followed STROBE-MR reporting guidance [11]. Blood gene-expression data came from the eQTLGen Consortium, comprising up to 31,684 participants, predominantly of European ancestry [12]. Disease outcomes were Biobank Japan (BBJ) GD, with 2,809 cases [13]; UK Biobank (UKB) hyperthyroidism, with 3,731 cases [14]; and FinnGen R12 Graves ophthalmopathy, with 858 cases [15] (Table 1).

The defined outcome hierarchy used BBJ for discovery, UKB for comparison across ancestry and a broader phenotype, and FinnGen as a TED-enriched sensitivity outcome. FinnGen cases were compared with population controls. Consequently, that analysis includes susceptibility to GD and cannot isolate susceptibility to eye disease among people who already have GD. We selected *TSHR*, *IGF1R* and *CTLA4* on biological and therapeutic grounds and evaluated them across all outcomes, regardless of discovery significance.

### Genetic instruments and statistical analysis

We considered 4,462 druggable genes [16]. Eligible instruments were genetic variants within 1 Mb of a gene that were associated with its blood expression at *P* < 5×10⁻⁸. We selected approximately independent variants by linkage-disequilibrium clumping in PLINK [17] against the 1000 Genomes Phase 3 European reference panel [18], which matched the ancestry of the exposure data. This left 6,135 instruments for 2,544 genes. The minimum instrument-strength F statistic was 29.7. Alleles were aligned between expression and disease datasets; ambiguous or unavailable variants were excluded. UKB estimates were converted from their reported linear-model scale to log-odds before analysis. Exposure effect sizes were reconstructed from the eQTLGen Z statistics, per-variant sample sizes and reference-panel allele frequencies [19]; the formulas, selection parameters and exclusions are provided in Supplementary Methods.

We used the Wald ratio for genes with one instrument and inverse-variance weighting for genes with multiple instruments [20]. Results are reported as odds ratios (ORs), 95% confidence intervals (CIs) and two-sided *P* values. ORs represent genetically proxied effects per approximately one standard deviation higher blood gene expression on the standardized scale assumed by the reconstruction. This scale assumes an expression variance of one; it does not standardize the genetic predictor itself. The estimates are not per drug dose and should not be interpreted as expected treatment effects. Discovery significance required *P* < 1.965×10⁻⁵, correcting for 2,544 eligible genes. Associations at *P* < 0.05 in the additional outcomes were considered nominal support only when their direction agreed with discovery. We report all estimates for the three selected genes, including non-significant results. Alternative MR estimators and tests for inconsistent instrument effects are presented in Table S1 [21, 22].

MR interpretation requires that instruments predict expression, are not associated with confounders, and do not affect disease through other pathways. Neither a small *P* value nor colocalization proves these assumptions. Single-instrument results, including *TSHR*, cannot undergo instrument-based tests of pleiotropy. We did not perform a formal genetic directionality test [23]. The reporting checklist records other analyses that were unavailable or not undertaken.

### Shared genetic signals and sensitivity analyses

For the three selected genes and six additional candidates, colocalization assessed whether gene expression and disease shared a genetic variant [8]. We considered a shared-variant posterior probability of at least 0.80 to provide strong support. This probability is a model-based measure of support, not a *P* value. We report it alongside the MR results and tested more conservative assumptions (Table S2) [24]. Colocalization used the regional expression and disease statistics without any significance filter, giving 1,479 to 7,740 variants per gene and outcome. Requiring support in both the discovery outcome and the TED-enriched outcome is a stringent cross-outcome filter: it selects candidates whose shared-variant evidence is not confined to one dataset, and it is not a test of Graves disease association on its own. The model assumes at most one causal variant per trait in the region; full model details and all probabilities are supplied in Supplementary Methods and Supplementary Data 3.

We checked instrument selection at *TSHR* and whether three neighbouring chromosome 16p11.2 hits represented independent findings. We also repeated the analyses using eQTLGen's own allele frequencies. Frequency substitution was assessed first with the same variants, then after repeating allele alignment. MR, colocalization and detectable effect sizes were compared with the primary analysis (Tables S3–S4). Details of these checks, including their limitations, are in Supplementary Methods.

A small orbital RNA-seq dataset comprising four TED samples and one control was retained solely as descriptive context (Figure S1). It did not contribute to gene selection or genetic conclusions, and no differential-expression *P* values were calculated for this manuscript. Technical replicates were combined at the biological-sample level.

Analyses used R 4.3.3, TwoSampleMR 0.7.4, coloc 5.2.3 and PLINK 1.9. Python reproduced the colocalization calculations for verification and frequency sensitivity. The instruments and complete numerical results are supplied in Supplementary Data 1–4. No prospective study registration was available.

## Results

### Discovery across druggable genes

Of 2,544 genes with selected instruments, 2,234 had estimable MR results in BBJ; the corresponding numbers were 2,505 in UKB and 2,480 in FinnGen. Thirteen genes met the BBJ multiple-testing threshold: seven were associated with lower odds and six with higher odds of disease (Table 3; Figure 1). The 310 genes without estimable BBJ results were not classified as negative findings. The genomic distribution shows both the MHC cluster and signals outside that region; each plotted point represents a gene-level MR result, not a single-SNP GWAS association. The discovery genes included *TSHR*, *CTLA4*, five genes in the major histocompatibility complex (MHC), and six additional candidates. The three biologically selected genes are compared in Table 2 and Figure 2, which aligns effect estimates with the corresponding shared-variant evidence.

### TSHR

Higher genetically proxied blood *TSHR* expression was associated with lower odds of BBJ GD (OR 0.12, 95% CI 0.07–0.21; *P* = 1.09×10⁻¹⁴). Estimates had the same direction in UKB (OR 0.09, 95% CI 0.06–0.14; *P* = 8.77×10⁻²⁸) and FinnGen (OR 0.10, 95% CI 0.04–0.24; *P* = 2.82×10⁻⁷).

Expression and disease showed strong shared-variant support in BBJ and FinnGen, with probabilities of 0.951 and 0.986, respectively. The leading variant was rs179252 in both. UKB did not show the same support (0.226), despite a strong MR association. There the posterior favoured distinct causal variants (PP.H3 = 0.774), so the expression and disease signals were more consistent with two linked variants than with one shared variant. Thus, a highly significant MR result alone did not establish a shared genetic signal in every outcome. The contrast is visible in Figure 2: the MR estimates all favoured lower odds, whereas the UKB posterior favoured distinct expression and disease variants. This discordance was not resolved by substituting cohort allele frequencies.

### IGF1R

Higher genetically proxied *IGF1R* expression was associated with higher odds of BBJ GD (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB hyperthyroidism (OR 1.35, 95% CI 1.07–1.70; *P* = 0.0117). These were nominal associations and did not meet the study-wide discovery threshold. The FinnGen estimate was in the same direction but imprecise (OR 1.41, 95% CI 0.85–2.33; *P* = 0.182).

Shared-variant support was below 0.80 in all three outcomes (BBJ 0.073, UKB 0.404 and FinnGen 0.032). The UKB posterior was divided between an expression-only signal (PP.H2 = 0.400) and a shared variant (PP.H4 = 0.404), with little support for distinct variants (PP.H3 = 0.196), and should be regarded as unresolved rather than negative. Alternative MR estimators agreed in direction but differed in statistical significance (Table S1). Weighted-median estimates were nominally significant in BBJ and UKB, whereas weighted-mode estimates were not; neither provided nominal support in FinnGen. These methods used the same limited instrument sets and therefore do not constitute independent replication. Taken together, these results provide weaker genetic support for *IGF1R* than for *TSHR* in BBJ Graves disease and FinnGen Graves ophthalmopathy; they do not show that *IGF1R* has no inherited contribution.

### CTLA4 and additional candidates

The known autoimmune locus *CTLA4* [25] was associated with lower odds of disease in BBJ (OR 0.18, 95% CI 0.11–0.27; *P* = 5.45×10⁻¹⁵), UKB (OR 0.21, 95% CI 0.08–0.57; *P* = 0.0022) and FinnGen (OR 0.17, 95% CI 0.04–0.68; *P* = 0.0118). Shared-variant support was strong in UKB and FinnGen (0.953 and 0.978), but not BBJ (0.201), where the posterior favoured distinct variants (PP.H3 = 0.799; Figure 2). The screen therefore recovered a known autoimmune association, while *CTLA4* itself failed the combined BBJ-plus-FinnGen colocalization criterion.

None of the six additional candidates met the shared-variant criterion in both BBJ and FinnGen. *TNFSF14* and *IFNGR1* had strong support in BBJ only (0.994 and 0.989); the corresponding FinnGen probabilities were 0.017 and 0.020 (Figure 3A). *TNFSF14* nevertheless retained a direction-consistent UKB MR association (*P* = 0.0056); the absence of colocalization in that outcome was not an absence of association. The three neighbouring chromosome 16p11.2 hits were consistent with one regional signal rather than three independent findings. Their BBJ shared-variant probabilities were below 0.80, and none gained strong support in UKB or FinnGen. *MAPKAPK5* also failed the shared-variant criterion despite passing the MR discovery threshold. Thus, the absence of an additional qualifying candidate resulted from examining the combined evidence, rather than from an absence of significant MR associations. All 13 discovery estimates and their classification are shown in Table 3; full comparisons are in Supplementary Data 2–3.

### Robustness and detectable effects

Using eQTLGen cohort frequencies preserved all 13 discovery hits. *IGF1R* remained nominally associated in BBJ and UKB (*P* = 0.022 and 0.016) and non-significant in FinnGen (*P* = 0.177; Table S3). Frequency substitution moved none of the 81 colocalization comparisons (nine genes × three outcomes × three priors) across the 0.80 threshold, so every classification was unchanged; the largest absolute change in shared-variant probability was 0.002142. However, the original prior sensitivity remained relevant: under the most conservative prior, *TSHR* support in BBJ fell from 0.951 to 0.661, while FinnGen remained above threshold (Figure 3B; Table S2). This distinguished stability to the frequency substitution from sensitivity to the assumed prior: the former preserved classifications, while the latter changed the interpretation of the BBJ *TSHR* signal.

The screen had limited ability to detect smaller associations. At 80% power, only 14.6% of estimable BBJ genes could detect an OR of 1.5 at the discovery threshold; this decreased to 12.5% with cohort frequencies. The corresponding proportions for an OR of 2.0 were 35.6% and 34.2% (Table S4). Gene-specific detection benchmarks are provided in Supplementary Methods. They describe precision and cannot distinguish an absent effect from limited power; the ORs, 95% CIs and *P* values in Table 2 should guide interpretation of the observed associations. These calculations describe MR association detection, not the probability of passing the complete colocalization filter. Failure to identify an additional candidate does not exclude moderate or other undetected genetic effects.

## Discussion

This study found stronger genetic evidence linking blood *TSHR* expression to GD than linking *IGF1R* expression to the same outcomes. *TSHR* showed direction-consistent MR associations across all three datasets and shared expression–disease signals in BBJ GD and FinnGen Graves ophthalmopathy under the primary priors. *IGF1R* showed nominal associations in BBJ and UKB, but the shared-variant evidence was unresolved. No additional candidate passed the combined BBJ and FinnGen criterion. The comparison is therefore informative about the pattern and strength of genetic evidence, rather than a ranking of expected drug efficacy. Figures 1–3 make the screening results, the association–colocalization differences and the prior dependence visible together.

For clinicians, the distinction between disease susceptibility and therapeutic response is central. Teprotumumab improved outcomes in a randomized trial of active TED [3], while experimental work supports functional interaction between TSHR and IGF-1R in orbital fibroblasts [5]. These observations concern receptor activity in established disease. Our analysis concerns inherited variation associated with blood gene expression and disease odds. Tissue-specific regulation documented by GTEx further cautions against treating blood expression as a complete representation of orbital biology [26]. A receptor can remain relevant to treatment even when this particular genetic analysis provides limited evidence for an expression-mediated susceptibility effect. The *IGF1R* findings neither challenge the clinical evidence for IGF-1R blockade nor establish an exclusively downstream effector role.

The *TSHR* findings provide a useful susceptibility anchor, but they require equally careful translation. The established role of TSHR in GD makes an association at this locus biologically coherent [1, 2]. Nonetheless, an OR below one for genetically proxied expression does not imply that increasing TSHR expression or stimulating its receptor would prevent disease. Genetic evidence can inform target prioritization [27], but the measured exposure, relevant tissue, timing and direction of a pharmacologic intervention must still match the proposed mechanism. Here, reconstructed blood expression is neither receptor stimulation nor a drug dose. The strongest defensible interpretation is that *TSHR*-region expression and GD susceptibility share genetic support under the primary model in BBJ, with a corresponding signal recovered in the TED-enriched FinnGen outcome.

The UKB result shows why MR and colocalization should be read together. *TSHR* had a highly significant MR association in UKB, yet the regional evidence favoured distinct expression and disease variants. Colocalization addresses a different question from the association test: whether the two traits plausibly share the same underlying variant [8]. A small *P* value cannot answer that question by itself. The broader hyperthyroidism phenotype and differences in ancestry or available variants are possible contributors to the outcome contrast, but this study cannot identify which explanation is responsible. Conversely, the divided H2/H4 posterior for *IGF1R* in UKB should remain unresolved. Calling it negative would discard the shared-variant component, while calling it positive would overstate its support.

The *CTLA4* comparison also helps define what the combined filter can and cannot establish. As a known autoimmune susceptibility locus [25], *CTLA4* showed the expected disease-association signal, but its BBJ posterior favoured distinct variants and it did not pass the combined BBJ-plus-FinnGen criterion. This makes it a biologically informative comparator rather than a validating positive control for the full analysis. Among the additional candidates, *TNFSF14* and *IFNGR1* had strong shared-variant support only in BBJ, whereas the chromosome 16p11.2 genes were consistent with a shared regional signal. These results support separating association discovery from gene-level prioritization. They complement previous MR studies of GD-related molecular traits [9, 10] by showing how an explicit comparison across outcomes can change the interpretation of statistically significant findings.

The sensitivity analyses improve confidence in specific aspects of the results without removing all uncertainty. Substituting eQTLGen cohort frequencies preserved the discovery-hit set and every tested shared-variant classification. It therefore addresses concern about the reference-frequency approximation used to reconstruct exposure effects. Prior sensitivity answers a separate question: *TSHR* support in BBJ fell below the chosen threshold under the most conservative shared-association prior, while FinnGen support remained above it. The BBJ finding should consequently be described as supported under the primary model and sensitive to the prior, rather than uniformly robust. Similarly, the screen's limited ability to detect moderate effects prevents the lack of additional qualifying candidates from excluding smaller inherited effects. Detection benchmarks quantify the design's precision; they do not determine whether a particular unresolved result reflects no effect or inadequate power.

A clinically useful next step is to test whether these susceptibility patterns relate to eye-disease development among patients who already have GD. A direct GD-with-TED versus GD-without-TED comparison would address that question more closely than population controls. Larger, ancestry-matched molecular datasets from relevant cell types could then assess whether blood-derived signals transfer to orbital tissue. The present orbital observations cannot provide that confirmation: four TED samples and one control support descriptive context only, and altered expression in established disease may reflect its consequences as well as its causes. Accordingly, the genetic and tissue components should guide focused follow-up questions rather than be combined into a claim of demonstrated therapeutic mechanism.

### Limitations

The main limitations concern phenotype definition, exposure relevance and available information. FinnGen compared 858 Graves ophthalmopathy cases with population controls rather than GD patients without eye disease, so its associations may reflect GD susceptibility and cannot isolate a TED-specific effect. UKB hyperthyroidism was broader than GD, making phenotype and ancestry differences partly confounded. The predominantly European blood-expression data may not transfer fully to East Asian discovery participants or orbital cell types because regulation, linkage patterns and variant availability differ. Participant overlap between eQTLGen and the European outcomes could not be quantified. Cohort-frequency substitution addressed one part of exposure reconstruction, but the published frequencies excluded Framingham participants and did not resolve these broader issues. The tissue series had only one control, was not externally replicated and cannot support differential-expression or causal inference.

Additional limitations concern the genetic models and detection range. *TSHR* was represented by one instrument under the European selection reference, precluding instrument-based pleiotropy tests; the other selected genes also had few instruments. Colocalization cannot exclude an effect through a neighbouring gene, and the single-causal-variant assumption may not capture more complex regional signals. The East Asian reference selected two *TSHR* variants at the primary threshold, so single-instrument status is not evidence for a single biological signal in every ancestry. Conclusions also depend on prior assumptions, particularly for *TSHR* in BBJ. Multi-variant colocalization, formal colocalization power and the probability of passing the combined filter were not evaluated. Together with the screen's limited MR detection range, these constraints leave smaller or otherwise undetected genetic effects unresolved.

### Conclusions

Within these limits, *TSHR* showed stronger expression–disease evidence than *IGF1R* in BBJ Graves disease and FinnGen Graves ophthalmopathy. The findings distinguish genetic susceptibility evidence from evidence of therapeutic efficacy, while preserving the uncertainty around *IGF1R* and the prior dependence of the BBJ *TSHR* result. They support further investigation in relevant tissues and a direct comparison of GD patients with and without TED; they do not establish treatment direction or expected clinical benefit.

## Declarations

**Funding.** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Conflict of interest.** The authors declare that they have no conflict of interest.

**Ethics approval.** This study used publicly available, de-identified summary statistics and an institutionally approved orbital tissue dataset. The in-house orbital transcriptomic component was approved by the Institutional Review Board of Pusan National University Hospital (approval number 2104-018-102) and was conducted in accordance with the Declaration of Helsinki.

**Informed consent.** Written informed consent was obtained from all individual participants included in the in-house orbital tissue study.

**Data availability.** Instruments, primary MR and colocalization results, and frequency-sensitivity results are supplied as Supplementary Data 1–4. Analysis code is available from the corresponding author on reasonable request. Public summary statistics analyzed in this study are available from their original repositories: blood *cis*-eQTL data from the eQTLGen Consortium, the Biobank Japan Graves disease and UK Biobank hyperthyroidism genome-wide association statistics through the GWAS Catalog, and the FinnGen Release 12 Graves ophthalmopathy statistics from FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions.

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

26. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369:1318-1330. doi:10.1126/science.aaz1776.

27. Nelson MR, Tipney H, Painter JL, Shen J, Nicoletti P, Shen Y, Floratos A, Sham PC, Li MJ, Wang J, et al. The support of human genetic evidence for approved drug indications. Nat Genet. 2015;47:856-860. doi:10.1038/ng.3314.

## Figure Legends

**Figure 1. Study design and the druggable-gene association screen.** (A) Gene selection and discovery flow. Genes with independent blood-expression instruments were evaluated across the outcome hierarchy; genes without estimable BBJ results were not treated as negative. (B) Gene-level MR associations in BBJ Graves disease, plotted by the GRCh37 position of each gene's strongest selected expression instrument. Each point represents one gene, not a SNP association. The horizontal line marks the discovery threshold of *P* < 0.05/2,544; all 2,234 estimable genes are plotted and all 13 discovery genes are labelled. *TSHR*, *IGF1R* and *CTLA4* are highlighted as the three biologically selected genes. Their selection for comparison did not depend on discovery significance.

**Figure 2. Association estimates and regional genetic evidence for the three selected genes.** (A) Primary ORs and 95% CIs, with corresponding two-sided *P* values. The OR axis is logarithmic, with a reference line at one. ORs refer to reconstructed genetically proxied blood expression, not treatment effects. All nine comparisons are included. (B) Regional colocalization probabilities under the primary prior, aligned with the same gene–outcome rows. Stacked bars show expression-only support (H2), support for distinct variants (H3), shared-variant support (H4), and residual H0/H1 support. Exact displayed H4 probabilities are given at right. These probabilities are not *P* values; ≥0.80 was the threshold for strong shared-variant support. The figure displays the discordance between association strength and shared-variant evidence without treating either as proof of a therapeutic effect.

**Figure 3. Shared-variant evidence across candidates and sensitivity to prior assumptions.** (A) Primary-prior PP.H4 values for all nine evaluated genes and three outcomes. Outlined cells meet the ≥0.80 threshold; probabilities below 0.001 are labelled <0.001. A gene passes the combined criterion only when both BBJ and FinnGen meet that threshold; UKB is an additional comparison. *TSHR*, *IGF1R* and *CTLA4* were biologically selected, while the other six genes were additional discovery candidates. (B) Prior sensitivity for the three selected genes. D, I and C denote default, intermediate and conservative shared-association priors (p12 = 10⁻⁵, 5×10⁻⁶ and 10⁻⁶), with p1 = p2 = 10⁻⁴ throughout. The dashed line marks 0.80. Prior changes can alter classification even when frequency substitution does not; Tables S3–S4 report the separate cohort-frequency sensitivity. Full H0–H4 probabilities and all 81 gene–outcome–prior settings are in Supplementary Data 3.

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

| Gene | Outcome | Instruments | OR (95% CI) | *P* value | Distinct variants (PP.H3) | Shared variant (PP.H4) |
|---|---|---|---|---|---|---|
| *TSHR* | BBJ Graves disease | 1 | 0.12 (0.07–0.21) | 1.09×10⁻¹⁴ | 0.049 | 0.951 |
| *TSHR* | UKB hyperthyroidism | 1 | 0.09 (0.06–0.14) | 8.77×10⁻²⁸ | 0.774 | 0.226 |
| *TSHR* | FinnGen Graves ophthalmopathy | 1 | 0.10 (0.04–0.24) | 2.82×10⁻⁷ | 0.014 | 0.986 |
| *IGF1R* | BBJ Graves disease | 4 | 1.56 (1.07–2.28) | 0.0212 | 0.236 | 0.073 |
| *IGF1R* | UKB hyperthyroidism | 4 | 1.35 (1.07–1.70) | 0.0117 | 0.196 | 0.404 |
| *IGF1R* | FinnGen Graves ophthalmopathy | 3 | 1.41 (0.85–2.33) | 0.182 | 0.346 | 0.032 |
| *CTLA4* | BBJ Graves disease | 1 | 0.18 (0.11–0.27) | 5.45×10⁻¹⁵ | 0.799 | 0.201 |
| *CTLA4* | UKB hyperthyroidism | 2 | 0.21 (0.08–0.57) | 0.0022 | 0.047 | 0.953 |
| *CTLA4* | FinnGen Graves ophthalmopathy | 2 | 0.17 (0.04–0.68) | 0.0118 | 0.022 | 0.978 |

ORs are per unit increase in reconstructed, genetically proxied blood expression. Wald ratio was used for one instrument and IVW for multiple instruments. *P* values are two-sided and shown to three significant digits. BBJ discovery significance requires *P* < 1.965×10⁻⁵; *P* < 0.05 in other outcomes is nominal support. PP.H3 and PP.H4 are the colocalization posterior probabilities of two distinct causal variants and of one shared causal variant; they are model-based measures of support, not *P* values, and ≥0.80 was the criterion for strong support. The two do not sum to one: the remainder is mostly PP.H2, an expression association without a resolved disease signal, which is why the *IGF1R* results are reported as unresolved rather than as evidence of distinct variants. Full PP.H0–PP.H4 are in Supplementary Data 3. All nine estimates are shown. Sensitivity results are in Tables S1–S3.

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
| *MAPKAPK5* | 4 | 0.98 (0.97–0.99) | 5.28×10⁻⁶ | Additional candidate |
| *PSMB8* | 1 | 1.25 (1.13–1.38) | 6.86×10⁻⁶ | MHC region |
| *IFNGR1* | 1 | 2.10 (1.51–2.91) | 9.41×10⁻⁶ | Additional candidate |

All genes met the same multiple-testing threshold and are listed by discovery *P* value. The *MAPKAPK5* estimate is 0.977 (0.967–0.987) at three decimals: a precisely estimated but very small per-unit effect, not a weak one. The five MHC genes were not treated as new independent targets. The six additional candidates were TNFSF14, IFNGR1, MAPKAPK5, HSD3B7, VKORC1 and PRSS36; none showed strong shared-variant support in both BBJ and FinnGen. Full outcome comparisons and model probabilities remain in Supplementary Data 2–3.

## Supplementary Material

## Supplementary Methods

### Expression instruments and allele alignment

Instrument selection used cis-eQTL *P* < 5×10⁻⁸, a ±1-Mb gene window, and LD clumping at r² < 0.001 within 10 Mb using 1000 Genomes Phase 3 European samples (n = 503). Exposure standard errors were reconstructed as SE = 1/√[2f(1−f)(N + Z²)] and β = Z × SE, where f is allele frequency and N is the variant-specific eQTL sample size [19]. This defines the reconstructed expression unit used for reported ORs.

TwoSampleMR harmonise_data with action = 2 aligned exposure and outcome alleles. Palindromic variants with intermediate frequency (MAF > 0.42) were excluded. Missing estimates were not imputed. BBJ and FinnGen effect estimates were used on their reported log-odds scales. UKB linear-model β and SE were divided by μ(1−μ), where μ = 3,731/484,598, to approximate log-odds estimates. Original GWAS covariate adjustments were retained; participant-level re-adjustment was unavailable.

### Checks using multiple instruments and regional linkage

For multiple instruments, the primary IVW analysis used TwoSampleMR's multiplicative random-effects method with an under-dispersion correction. Weighted-median and weighted-mode estimates, the MR-Egger intercept and Cochran's Q were inspected where estimable (Table S1) [21, 22]. With three or four instruments, a non-significant intercept has limited ability to exclude pleiotropy. Steiger directionality [23], bidirectional MR and leave-one-out analyses were not reported. Single-instrument analyses were limited to the Wald ratio.

At *TSHR*, clumping was repeated using eQTL thresholds of 5×10⁻⁸ and 5×10⁻⁶, r² thresholds of 0.001, 0.01 and 0.1, and European (n = 503) or East Asian (n = 504) reference panels. rs179252 remained an index variant in all 12 combinations. The primary settings selected one instrument in Europeans and two in East Asians; relaxed settings yielded up to 12. These checks support retention of the selected instrument, not a single-signal model in every ancestry. The full clumping counts are preserved with the analysis records.

The chromosome 16p11.2 variants rs4889606, rs34649473 and rs78924645 span approximately 143 kb. Their East Asian pairwise r² values were 0.854, 0.761 and 0.863. GCTA-COJO conditional selection at *P* < 5×10⁻⁶ retained rs8050588 (*P* = 1.15×10⁻⁸). After conditioning, the three marginal associations had *P* values of 0.928, 0.839 and 0.806, respectively. This supports regional dependence and does not identify a particular causal gene.

### Colocalization and prior assumptions

coloc.abf was applied to available gene-specific cis-eQTL variants matched to each outcome by rsID and allele pair, avoiding comparisons of genomic windows across GRCh37 expression and GRCh38 outcome coordinates. Duplicate rsIDs and allele-pair mismatches were excluded. Both datasets used β and its variance. Dataset 1 was the outcome GWAS and dataset 2 the expression exposure; expression sdY was set to 1. Outcome allele frequencies were not used in the Bayes factors.

The five model probabilities describe neither trait associated (H0), outcome only (H1), expression only (H2), both traits with distinct variants (H3), and both traits with a shared variant (H4). The default priors were p1 = p2 = 10⁻⁴ and p12 = 10⁻⁵; p12 values of 10⁻⁶ and 5×10⁻⁶ tested more conservative shared-association assumptions [24]. Table S2 presents the selected genes; Supplementary Data 3 retains all nine genes, three outcomes and three priors. Low H4 does not demonstrate that disease association is absent. Neither multi-variant colocalization nor formal colocalization power was evaluated.

### Cohort-frequency sensitivity and detection limits

The eQTLGen frequency release represents 26,609 participants and excludes the Framingham Heart Study. All 50,245 requested SNPs, including all 6,135 instrument records, had valid frequencies and matching GRCh37 positions. Frequencies were aligned to the assessed allele, then β and SE were reconstructed while retaining the original Z and sample size. The maximum change in Z was 5.54×10⁻¹³. We compared identical retained variants and then repeated action = 2 harmonization. Across outcomes, 27 inclusion decisions changed: 24 involved palindromic frequency decisions and three restored missing reference frequencies. No retained outcome estimate changed direction. Estimable counts became 2,232 BBJ, 2,506 UKB and 2,481 FinnGen genes; the multiple-testing denominator remained 2,544.

The same frequency substitution was applied to all 81 colocalization settings. Maximum absolute H4 changes were 0.000245 on identical variants and 0.002142 after additionally including variants with newly available frequencies; none crossed 0.80. Primary MR was independently checked by a second implementation. Frequency substitution reran Wald/IVW estimates, colocalization and detection thresholds; it did not rerun every alternative MR estimator in Table S1.

For comparison with the primary analysis, the minimum detectable absolute MR effect at 80% power was approximated as (z₁₋α/₂ + z₀.₈₀) × SE. Alpha was 0.05/2,544 in BBJ and 0.05 in UKB/FinnGen. Exponentiation gives the detectable OR. Percentages use the estimable-gene denominator for each scenario (Table S4). These are normal-approximation association-detection metrics based on observed standard errors, not the power of colocalization or of the combined selection rule. They do not turn a non-significant association into evidence of absence.

At 80% power, the approximate minimum detectable absolute log-ORs were 1.39 for *TSHR* and 0.99 for *IGF1R* in BBJ (α = 0.05/2,544), and 0.33 and 0.72 for *IGF1R* in UKB and FinnGen, respectively (α = 0.05). These standard-error-based benchmarks describe precision; comparing them with observed effects cannot distinguish an absent effect from limited power. The ORs, 95% CIs and *P* values in Table 2 should guide interpretation of the observed associations.

### Descriptive orbital expression

Orbital RNA-seq comprised four TED patients and one control. Technical replicates were combined by biological sample, and abundance was expressed as transcripts per million. The plotted log2 fold change was log2[(mean TED TPM + 0.01)/(control TPM + 0.01)]. With one control, control-group biological variance cannot be estimated; no differential-expression *P* value or confirmatory test is reported. The dataset was not externally replicated and contributes no inferential evidence to the genetic comparison (Figure S1).

**Table S1. Sensitivity tests for genes with multiple instruments**

| Gene | Outcome | Instruments | Primary *P* | Median *P* | Mode *P* | Egger intercept *P* | Heterogeneity *P* |
|---|---|---|---|---|---|---|---|
| *IGF1R* | BBJ Graves disease | 4 | 0.0212 | 0.0288 | 0.512 | 0.663 | 0.902 |
| *IGF1R* | UKB hyperthyroidism | 4 | 0.0117 | 0.00244 | 0.367 | 0.562 | 0.25 |
| *IGF1R* | FinnGen Graves ophthalmopathy | 3 | 0.182 | 0.124 | 0.633 | 0.433 | 0.309 |
| *CTLA4* | UKB hyperthyroidism | 2 | 0.0022 | NA | NA | NA | 0.00131 |
| *CTLA4* | FinnGen Graves ophthalmopathy | 2 | 0.0118 | NA | NA | NA | 0.0486 |

Median and mode denote weighted-median and weighted-mode MR. The Egger intercept assesses directional pleiotropy; heterogeneity uses Cochran’s Q. NA means not estimable, not a negative test. TSHR and CTLA4 in BBJ each used one instrument and are omitted from this diagnostic table. IGF1R estimates agreed in direction across estimators, but significance differed. Few instruments limit all these sensitivity checks. These are the original reference-frequency diagnostics; frequency substitution repeated primary Wald/IVW only.

**Table S2. Shared-variant support under different prior assumptions**

| Gene | Outcome | Default prior | Intermediate prior | Conservative prior |
|---|---|---|---|---|
| *TSHR* | BBJ Graves disease | 0.951 | 0.907 | 0.661 |
| *TSHR* | UKB hyperthyroidism | 0.226 | 0.128 | 0.028 |
| *TSHR* | FinnGen Graves ophthalmopathy | 0.986 | 0.972 | 0.875 |
| *IGF1R* | BBJ Graves disease | 0.073 | 0.038 | 0.008 |
| *IGF1R* | UKB hyperthyroidism | 0.404 | 0.253 | 0.063 |
| *IGF1R* | FinnGen Graves ophthalmopathy | 0.032 | 0.016 | 0.003 |
| *CTLA4* | BBJ Graves disease | 0.201 | 0.112 | 0.025 |
| *CTLA4* | UKB hyperthyroidism | 0.953 | 0.911 | 0.672 |
| *CTLA4* | FinnGen Graves ophthalmopathy | 0.978 | 0.957 | 0.816 |

Entries are PP.H4, not *P* values. The shared-association prior p12 was 10⁻⁵ (default), 5×10⁻⁶ (intermediate) and 10⁻⁶ (conservative); p1 = p2 = 10⁻⁴. TSHR in BBJ and CTLA4 in UKB fall below 0.80 under the conservative prior. TSHR in FinnGen remains above it. Full H0–H4 probabilities, variant counts and all nine loci are in Supplementary Data 3.

**Table S3. Associations before and after cohort-frequency substitution**

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

**Table S4. Association detection limits before and after frequency substitution**

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
