# TSHR and IGF1R in Graves disease and thyroid eye disease: genetic evidence from blood expression

Jungyul Park¹, Min-Seon Kim², Kyung-Hwa Shin³\*, Suk-Woo Yang¹\*

¹ Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
² Department of Ophthalmology, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
³ Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea  

\* These authors contributed equally as corresponding authors.

**Corresponding author:** Suk-Woo Yang, MD, PhD, Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea. E-mail: yswoph@catholic.ac.kr; Tel: +82-2-2258-2847.  
**Co-corresponding author:** Kyung-Hwa Shin, MD, Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea.



## Abstract

**Objective:** Thyroid eye disease (TED) can accompany Graves disease (GD), an autoimmune cause of hyperthyroidism. We compared genetic evidence linking blood *TSHR* and *IGF1R* expression with susceptibility.

**Methods:** We screened 2,544 druggable genes using Mendelian randomization (MR). Disease outcomes were Biobank Japan (BBJ) GD, UK Biobank (UKB) hyperthyroidism and overlapping FinnGen GD, Graves ophthalmopathy (GO) and strict GO endpoints. Colocalization tested 24 genes for a shared causal variant. Tissue eQTL and orbital RNA-sequencing (four TED, one control) were supplementary.

**Results:** Thirteen genes passed BBJ discovery (*P* < 0.05/2,544). Higher genetically proxied *TSHR* expression was associated with lower GD odds in BBJ (odds ratio [OR] per approximately one standard deviation 0.12, 95% confidence interval [CI] 0.07–0.21; *P* = 1.09×10^−14^) and FinnGen (OR 0.08, 95% CI 0.05–0.12; *P* = 3.19×10^−32^), concordantly across all five datasets. Shared-variant probability exceeded 0.80 in FinnGen GO and, prior-sensitively, BBJ and strict GO, but UKB (0.226) and FinnGen GD (0.264) favored distinct variants. *IGF1R* ORs exceeded 1 but reached *P* < 0.05 only in BBJ (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB (*P* = 0.0117); the most precise estimate, FinnGen GD, was OR 1.07 (95% CI 0.86–1.33). Single-variant tissue estimates added no independent evidence. Mean orbital *TSHR* and *IGF1R* transcript levels exceeded the control.

**Conclusions:** Association evidence was stronger for *TSHR* than *IGF1R*, but absent colocalization in UKB and FinnGen GD leaves expression-mediated causality unproven. Population-controlled data cannot isolate TED susceptibility within GD; weaker *IGF1R* evidence does not argue against IGF-1R–targeted treatment.

**Keywords:** Graves disease; thyroid eye disease; *TSHR*; *IGF1R*; Mendelian randomization; colocalization

## Introduction

Thyroid eye disease (TED) can cause proptosis, diplopia and visual impairment in patients with Graves disease (GD) [1, 2]. Teprotumumab, which blocks the insulin-like growth factor 1 receptor (IGF-1R), improves clinical outcomes in active TED [3]. Its efficacy raises a related but distinct question: does inherited variation influencing *IGF1R* expression contribute to disease susceptibility in the same way as variation at *TSHR*, the established GD autoantigen and susceptibility locus?

TSHR and IGF-1R interact in experimental models of orbital disease [4, 5, 6]. However, treatment response, inherited gene regulation and gene expression in inflamed tissue address different biological questions. A receptor may be a useful treatment target without equally strong genetic evidence linking its expression to disease onset. Distinguishing these forms of evidence helps clinicians interpret genetic studies of therapeutic targets.

Mendelian randomization (MR) uses variants associated with an exposure to investigate its relationship with the disease phenotype being studied, termed the outcome [7]. Here, the exposure is blood gene expression and the outcomes are disease diagnoses, not responses to treatment. Because alleles are allocated at conception, MR is less prone than ordinary observational studies to confounding and reverse causation, but it supports causal interpretation only when its instrumental-variable assumptions hold. Colocalization asks whether expression and disease associations in a genomic region are consistent with the same causal variant [8]. Distinct nearby variants can appear connected because they are inherited together. Thus, an MR association may be highly significant even when the regional signals do not support a shared variant. MR reports odds ratios (ORs), confidence intervals (CIs) and *P* values; colocalization reports model-based posterior probabilities, which are not *P* values.

Previous genetic studies have already implicated *TSHR* in GD and Graves ophthalmopathy. Ji and colleagues integrated blood methylation, expression and protein associations with GD, reporting *TSHR* among the expression-associated genes also supported in a strict FinnGen ophthalmopathy dataset [9]. Other MR work has examined inflammatory proteins in ophthalmopathy [10], and recent integrative research has investigated *TSHR* across genetic and transcriptomic data in GD and ophthalmopathy [11]. These findings establish a basis for comparing clinically relevant receptors, but do not determine whether evidence for their blood expression is consistently supported across disease definitions under the same analysis framework. In particular, limited evidence for *IGF1R* could reflect imprecision or dependence on particular instruments, while a strong MR association for *TSHR* does not by itself establish that expression and disease share a causal variant. The remaining question is how robustly these two expression traits are associated with GD and related phenotypes, and how far those associations support biological interpretation.

To address this question, we compared genetically proxied blood expression of *TSHR* and *IGF1R* — expression predicted from inherited variants rather than measured directly — within a screen of genes encoding potentially druggable proteins, using *CTLA4* and other associated genes as context. We used Biobank Japan (BBJ) GD for discovery, compared the original findings with UK Biobank (UKB) hyperthyroidism and FinnGen ophthalmopathy, and subsequently added FinnGen GD and a stricter ophthalmopathy definition. We evaluated effect direction, precision and statistical support across these datasets, then examined dependence on individual instruments, allele frequencies and colocalization assumptions. A later fixed-panel tissue extension evaluated whether these patterns extended beyond the original blood instruments. Our aim was to characterize the consistency and limits of expression-related genetic evidence for the two receptors. The analyses distinguish GD from broader hyperthyroidism and population-controlled ophthalmopathy; they do not directly estimate TED susceptibility among patients who already have GD or predict treatment response. A small in-house orbital RNA-sequencing series (four TED samples and one control) is reported only as descriptive context.

## Methods

### Study design and disease datasets

This two-sample MR study used public association summary statistics and followed STROBE-MR reporting guidance [12]. The exposure used whole-blood expression quantitative trait locus (eQTL) statistics from eQTLGen, comprising up to 31,684 participants, predominantly of European ancestry [13]. The original outcome hierarchy comprised BBJ GD, UKB self-reported hyperthyroidism or thyrotoxicosis, and FinnGen R12 Graves ophthalmopathy [14, 15, 16]. We subsequently added FinnGen R12 E4_GRAVES_STRICT and E4_GRAVES_OPHT_STRICT. The five datasets and their case and control counts are shown in Table 1. The UKB source did not establish an exclusively European analysis sample, so ancestry was not inferred from the reference panel used in its source study.

BBJ remained the discovery dataset. FinnGen GD provided a same-disease comparison, UKB a broader-phenotype comparison, and the two ophthalmopathy endpoints sensitivity analyses of a TED-enriched phenotype. The three FinnGen endpoints contain overlapping participants and were not treated as independent replications or combined by meta-analysis. Their controls were population non-cases for the respective endpoint, not GD patients without TED. The FinnGen extensions, expanded candidate assessment and instrument diagnostics were post hoc; no prospective registration was available.

### Expression instruments and MR

From 4,462 potentially druggable genes [17], instruments were selected from significant blood cis-eQTLs (*P* < 5×10^−8^) within 1 Mb of the gene. PLINK clumping used a linkage disequilibrium (LD) threshold of r² < 0.001 within 10 Mb and 503 European participants from 1000 Genomes Phase 3 [18, 19]. This produced 6,135 single-nucleotide polymorphism (SNP)–gene records for 2,544 genes. The minimum F statistic was 29.7; the median variant-specific expression sample size was 28,092, with a range of 232–31,684. Exposure coefficients and standard errors were reconstructed from Z statistics, sample sizes and reference allele frequencies under unit expression variance [20]. Alleles were aligned with TwoSampleMR; ambiguous palindromic or unavailable variants were excluded. UKB linear-model coefficients and standard errors were divided by μ(1−μ), where μ = 3,731/484,598, to approximate a log-odds scale. BBJ and FinnGen coefficients already used log odds.

We used Wald ratios for single-instrument genes and multiplicative random-effects inverse-variance weighting (IVW), with correction for under-dispersion, for multiple instruments [21]. An instrumental variable (IV) is a genetic variant used to represent variation in expression. Results are reported as ORs per approximately one standard deviation higher blood expression, with 95% CIs and two-sided *P* values. This scale describes reconstructed expression effects, not a drug dose or an intervention. Primary discovery required *P* < 0.05/2,544 (1.965×10^−5^), retaining all eligible genes in the denominator even when an outcome estimate was unavailable. The 13 BBJ discoveries were assessed in FinnGen GD using a focused threshold of 0.05/13 and concordant effect direction. A separate full FinnGen screen used 0.05/2,544 and was exploratory. Nominal *P* < 0.05 was reported descriptively and was not equated with discovery significance.

MR interpretation assumes that instruments predict expression, are not associated with exposure–disease confounders, and influence disease only through the modeled exposure. A small *P* value or high PP.H4 cannot by itself establish these assumptions.

Where applicable, sensitivity assessment included Cochran's Q, alternative MR estimators and omission of individual instruments [22, 23]. Previously completed alternative-estimator results (MR-Egger, weighted median and weighted mode) for the original three outcomes were retained (Supplementary Dataset, MR_estimators_original); they were not represented as new analyses of every endpoint. Post hoc leave-one-out results assessed dependence on individual variants, not independent replication. Single-instrument *TSHR* results could not undergo these omission or multi-instrument pleiotropy tests.

In a post hoc allele-direction audit, we re-extracted rs179252 and rs179247 from the original expression and disease files. Disease effects were aligned to the blood expression-increasing allele, and G-allele correlations were calculated in European and East Asian reference genotypes and independently checked with PLINK. This audit assessed allele alignment and consistency with published thymus findings, rather than thymic mediation (Supplementary Dataset, TSHR_allele_audit).

### Colocalization and robustness assessment

We report MR associations as the main findings and use colocalization to assess their interpretation. In the original analysis, colocalization was reported alongside MR as a parallel evidence layer; the present emphasis was adopted during post hoc reappraisal and was not a prospectively registered hierarchy. Colocalization used complete available regional cis-eQTL statistics, including non-significant variants, matched to each disease dataset [8]. The expanded set included all 13 BBJ discoveries, *IGF1R*, nine additional FinnGen GD discoveries and *HLA-G*, which was selected because its MR significance changed after frequency substitution: 24 genes in total. We used coloc.abf with disease as dataset 1 and expression as dataset 2. Posterior probability of hypothesis 4 (PP.H4) measures support for a shared causal variant under the fitted model. PP.H3 instead supports distinct variants. We used PP.H4 ≥ 0.80 as a descriptive strong-support threshold, with priors p1 = p2 = 10^−4^ and p12 = 10^−5^, and repeated analyses at p12 = 5×10^−6^ and 10^−6^ [24]. All five hypotheses are defined in Supplementary Methods. Low PP.H4 alone does not establish no disease association, and H3 does not prove multiple independent signals.

We evaluated reference and eQTLGen cohort frequencies, distinguishing changes in effect scaling from changes in variant availability and harmonization. The 24-gene colocalization dataset contains 1,440 settings: 24 genes × five outcomes × four frequency/availability scenarios × three priors. For *TSHR*, we additionally examined narrower windows, a common five-outcome variant set, raw regional disease peaks and all available source cis-eQTLs. For newly evaluated candidates, we checked marginal signal coverage, shared instruments and variant aliases. *C4B* and *HLA-G* underwent instrument-level frequency, weighting, heterogeneity and omission checks. Full-precision results are supplied in one Supplementary Dataset workbook. The original BBJ–ophthalmopathy combined filter is retained as a historical descriptive comparison, not promoted to a prospectively specified validation rule or the primary conclusion of this reappraisal.

### Tissue-specific extension

In a post hoc extension, we assessed both genes in a fixed panel of five GTEx v8 tissues: thyroid, cultured fibroblasts, subcutaneous adipose, visceral adipose and whole blood [25]. The extension plan was fixed before retrieving these tissue results, after the blood analysis was known. GTEx participants were predominantly of European ancestry, not exclusively European. Instruments required cis-eQTL *P* < 5×10^−8^, F > 10 and European-reference LD r² < 0.001. Ten planned gene–tissue comparisons against FinnGen GD formed the primary family (*P* < 0.005); 40 comparisons against the other four datasets formed a secondary family (*P* < 0.00125). Unestimable combinations remained in the correction denominators. Effects used the original slopes and standard errors for inverse-normal-transformed expression. Complete regional summaries were not obtained for this extension, so tissue colocalization was not performed; significant-only pairs were used solely for instrument selection. All planned combinations are reported in the Supplementary Dataset.

### Orbital expression and implementation

An in-house orbital RNA-sequencing series of four patients with inactive TED and one control provided descriptive context only (Supplementary Methods). Technical replicates were combined by biological sample, and transcript abundance was displayed without differential-expression significance testing (Figure S1). These observations were not used to select genes or validate the genetic findings. Analyses used R 4.3.3, TwoSampleMR 0.7.4, coloc 5.2.3 and PLINK 1.9, with independent numerical verification in Python. No age- or sex-stratified patient analyses, individual-level multivariable regression, Steiger directionality test or multi-signal colocalization was performed.

## Results

### Discovery across druggable genes

Of the 2,544 instrumented genes, 2,234 were estimable in BBJ, 2,505 in UKB and 2,480 in each FinnGen endpoint. Thirteen genes passed the BBJ screen threshold (*P* < 0.05/2,544), including *TSHR* and *CTLA4* but not *IGF1R* (Figure 1). At this threshold, 14.6% and 35.6% of estimable BBJ genes had 80% power to detect ORs of 1.5 and 2.0, respectively (Supplementary Dataset, Detection_original).

### TSHR and IGF1R across disease datasets

Higher genetically proxied blood *TSHR* expression was associated with lower odds of GD in BBJ (OR 0.12, 95% CI 0.07–0.21; *P* = 1.09×10^−14^) and FinnGen (OR 0.08, 95% CI 0.05–0.12; *P* = 3.19×10^−32^). ORs were also below 1 for UKB hyperthyroidism and both ophthalmopathy endpoints, and all estimates except strict ophthalmopathy passed the screen threshold (Table 2; Figure 2). Each estimate used the single instrument rs179252. Its expression-increasing G allele had per-allele disease ORs of 0.77–0.82 across datasets and was correlated with rs179247-G (European r = 0.978; East Asian r = 0.996), which also increased blood *TSHR* expression. Allele-aligned recalculation reproduced the MR estimates (Supplementary Dataset, TSHR_allele_audit).

All five *IGF1R* ORs exceeded 1, and none passed the screen threshold. Nominal associations (uncorrected *P* < 0.05) were found only in BBJ GD (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB hyperthyroidism (OR 1.35, 95% CI 1.07–1.70; *P* = 0.0117). The most precise estimate, FinnGen GD, was OR 1.07 (95% CI 0.86–1.33; *P* = 0.539); this comparison had 80% power to detect an OR of 1.36 at *P* < 0.05 (Supplementary Methods). Weighted-median ORs were 1.56 (*P* = 0.0288) in BBJ, 1.42 (*P* = 0.00244) in UKB and 1.46 (*P* = 0.124) in FinnGen ophthalmopathy. In these three datasets, MR-Egger and weighted-mode estimates had *P* ≥ 0.32, and MR-Egger intercept *P* values were 0.66, 0.56 and 0.43, respectively (Supplementary Dataset, MR_estimators_original).

### Other associated genes

Four genes passed the screen threshold in both GD datasets: *TSHR*, *CTLA4*, *HLA-DQA2* and *TUBB* (Table 3). *CTLA4* was associated with lower odds, and *HLA-DQA2* and *TUBB*, both in the major histocompatibility complex (MHC) region, with higher odds. *CTLA4* and *TUBB* had direction-consistent nominal associations in all five datasets; *HLA-DQA2* did in four, but its two-instrument UKB estimate was non-significant and below 1 (Supplementary Dataset, MR_primary). At the focused FinnGen GD threshold (*P* < 0.05/13 = 0.00385), seven of the 13 BBJ discoveries had direction-consistent associations: these four genes plus *C4A*, *TNFSF14* and *IFNGR1*. *HLA-A*, which had the smallest BBJ *P* value, had a non-significant FinnGen GD estimate.

The exploratory full FinnGen GD screen identified 13 genes, nine of them outside the BBJ discovery set (Table S1). Seven of the nine were instrumented in the extended MHC region, which was not excluded from the screen, and two pairs each shared one instrument.

### Instrument and frequency sensitivity

Omitting rs13030124 left non-significant *CTLA4* estimates in UKB and all three FinnGen endpoints (*P* = 0.323–0.972; Figure S2). In BBJ and UKB, where *IGF1R* had four instruments, omitting each in turn gave *P* > 0.05 in two of the four estimates; omitting rs2654980 in UKB also lowered the OR from 1.35 to 1.16 (Figure S2). Among the BBJ discoveries, Cochran's Q test in FinnGen GD indicated heterogeneity among instruments for *C4A* (*P* = 1.27×10^−12^), *TUBB* (*P* = 3.67×10^−5^), *HLA-A* (*P* = 0.0387) and *CTLA4* (Table 2). Omitting rs149543464 (*TUBB*) or rs1269852 (*C4A*) raised the FinnGen GD *P* value to 0.0262 and 0.0471, respectively (Supplementary Dataset, LOO_new_fixed14).

With eQTLGen cohort allele frequencies, all 13 BBJ discoveries still passed the screen threshold, *TSHR* *P* values were unchanged and no *IGF1R* *P* value crossed 0.05 (Supplementary Dataset, MR_original_freq and MR_new_freq). Screen classification changed only in FinnGen GD: *C4B* no longer passed, and *HLA-G* newly passed; the IVW weight of rs150254595 in its estimate rose from 21.3% to 58.9% (Table S2; Supplementary Dataset, Instrument_weights). Omitting this variant left non-significant *HLA-G* estimates under either frequency source (Supplementary Dataset, LOO_frequency).

### Colocalization

*TSHR* PP.H4 (probability of one shared causal variant) was 0.951 in BBJ, 0.986 in FinnGen ophthalmopathy and 0.891 in strict ophthalmopathy, but 0.226 in UKB and 0.264 in FinnGen GD, where PP.H3 (distinct variants) was highest (Table 2; Figures 2 and 3). The UKB and FinnGen GD values were unchanged in narrower regions and in a variant set common to all five datasets, and the matched variants included the strongest regional disease variant in each dataset (Table S3). Lowering the prior probability of a shared variant (p12) from 10^−5^ to 10^−6^ gave PP.H4 of 0.875 in FinnGen ophthalmopathy, 0.661 in BBJ and 0.451 in strict ophthalmopathy.

*IGF1R* PP.H4 ranged from 0.010 to 0.404 (UKB); in the other four datasets, PP.H2 (expression signal only) was highest (0.571–0.731). *CTLA4* PP.H4 was ≥ 0.80 in UKB and both ophthalmopathy endpoints but not in BBJ or FinnGen GD. None of the 24 genes reached PP.H4 ≥ 0.80 in FinnGen GD, and *HLA-DQA2* and *TUBB* had PP.H4 < 0.001 in every dataset. Apart from *TSHR* and *CTLA4*, only *TNFSF14* and *IFNGR1* reached PP.H4 ≥ 0.80 in any dataset, both in BBJ only (Figure 3).

### Tissue-specific extension

Four of the ten gene–tissue combinations had a qualifying instrument, each a single SNP: *TSHR* in subcutaneous and visceral adipose and *IGF1R* in thyroid and cultured fibroblasts. Thus, 20 of the 50 planned comparisons were estimable and 30 were not (Table S4). All ten adipose *TSHR* ORs were below 1, and eight, including both primary FinnGen GD comparisons, passed the multiple-testing threshold of their test family. The adipose instruments, rs12101255 (subcutaneous) and rs74064796 (visceral), lie 50 bp apart in *TSHR* intron 1, about 15 kb from rs179252 (European r² 0.46 for rs12101255), and both are blood *TSHR* eQTLs (eQTLGen *P* = 4.1×10^−34^ and 6.5×10^−34^). With one instrument, each tissue *P* value equals that variant's disease-association z-test *P* value, so the adipose estimates re-express the regional *TSHR* association. No thyroid or cultured-fibroblast *IGF1R* estimate reached *P* < 0.05, and directions were inconsistent across datasets.

### Descriptive orbital expression

Mean orbital *TSHR*, *IGF1R* and *CTLA4* transcript levels in four patients with inactive TED exceeded those in the single control, whose sample came from a different fat depot; no significance testing was performed (Figure S1).

## Discussion

Higher genetically proxied blood *TSHR* expression was associated with lower GD odds in BBJ and FinnGen GD, with concordant directions for UKB hyperthyroidism and both ophthalmopathy endpoints. All *IGF1R* ORs exceeded 1, but only the BBJ and UKB estimates reached uncorrected *P* < 0.05, none met the discovery threshold, and the most precise estimate lay close to the null. Under the primary prior, colocalization supported a shared variant for *TSHR* in BBJ and both ophthalmopathy endpoints but favored distinct variants in UKB and FinnGen GD. These results do not establish a difference in biological importance between the receptors, and statistical significance in one gene but not another is not a formal test of different effects.

Our findings extend an existing body of genetic evidence. Ji and colleagues had already linked blood *TSHR* expression to GD and strict ophthalmopathy in an earlier FinnGen release [9]; recent integrative work has also investigated *TSHR* in GD and ophthalmopathy [11]. We add a focused comparison with *IGF1R*, a GD comparison spanning BBJ and FinnGen, and an assessment of both genes alongside the wider druggable-gene screen. Because Ji and colleagues also used eQTLGen and FinnGen, agreement is not wholly independent replication. Different MR methods, FinnGen releases and colocalization settings also preclude interpreting differences between reports as direct contradictions.

Among the other associated genes, the multi-instrument FinnGen and UKB *CTLA4* estimates relied heavily on rs13030124; with that caveat, *CTLA4*'s direction-consistent associations across all five outcomes support its use as an autoimmune comparator [26]. *HLA-DQA2* was discordant only in UKB, where its two instruments disagreed, while FinnGen GD *TUBB* and *C4A* estimates each depended largely on one variant. The strongest BBJ association, *HLA-A*, did not persist in FinnGen GD, so the smallest discovery *P* value need not identify the most transferable association. Shared instruments, complex MHC linkage and, for *HLA-DQA2* and *TUBB*, absent shared-variant support further prevent counting every associated expression trait as an independent causal target.

The lower-odds direction for *TSHR* need not contradict the pathogenic role of receptor-stimulating antibodies. Blood transcript abundance, tissue-specific expression and receptor activation are different exposures, and genetic regulatory effects can differ across tissues [25]. Our audit linked the blood expression-increasing rs179252-G allele to rs179247-G, the GD-protective allele associated with greater thymic transcription in published studies [27, 28]. This directional consistency is compatible with a central-tolerance explanation, but the thymic step was not tested here, and neither linked variant is a proven causal regulator. The marked departure of the MR ORs from 1 reflects scaling: rs179252 shifts reconstructed blood expression by only about 0.11 SD, and its per-allele disease ORs were 0.80 (BBJ), 0.77 (UKB), 0.78 (FinnGen ophthalmopathy), 0.77 (FinnGen GD) and 0.82 (strict ophthalmopathy). Rescaling these moderate per-allele effects to one SD of expression yields the per-SD ORs, which describe lifelong inherited differences, not the effect of raising *TSHR* expression or of TSHR-directed treatment. Because the adipose instruments tag the same intron-1 region as rs179252 and are themselves blood eQTLs, their concordant estimates re-express the regional association; they cannot resolve the distinct-variant pattern in UKB and FinnGen GD or establish an adipose-mediated mechanism, and neither depot is orbital fat.

Colocalization does not replace the MR associations, but it constrains causal wording. Strong *TSHR* associations in UKB and FinnGen GD coexisted with distinct-variant regional patterns, and BBJ and strict ophthalmopathy support weakened under a lower shared-variant prior. Thus, association strength and shared-variant evidence were not interchangeable, as also seen for the comparator *CTLA4* in BBJ and FinnGen GD. UKB and FinnGen GD had the most cases of the five datasets; their distinct-variant *TSHR* patterns could reflect finer resolution of a complex regional architecture, but this was not tested. Ji and colleagues chose a shared-variant prior intended to allow weaker signals to colocalize and accepted PP.H4 above 0.50 [9], whereas our primary settings used p12 = 10^−5^ and PP.H4 ≥ 0.80, so colocalization labels are not directly comparable between the studies. The regional checks excluded a simple missing-peak explanation for *TSHR* but did not establish the correct causal architecture.

Experimental and clinical studies provide context for relating TSHR-directed immunity to orbital disease. TSHR can undergo cleavage and release its extracellular A-subunit. Expression of this subunit in an adenoviral immunization model induced greater thyroid-stimulating autoimmunity than a minimally cleaving receptor [29]; plasmid immunization subsequently reproduced orbital adipogenesis, fibrosis and muscle pathology in mice [30]. Daser and colleagues further reported enhanced inflammatory and hyaluronan responses to TNFα in orbital fibroblasts from TSHR-immunized mice [31]. Graves ophthalmopathy has also been documented in primarily euthyroid and hypothyroid patients [32], supporting a distinction between thyroid functional status and orbital immune activity. These observations make the ophthalmopathy associations biologically plausible, but our instruments do not measure receptor cleavage or antibody activity. Because most patients with TED have GD and the ophthalmopathy endpoints used population controls, their MR estimates may largely reflect GD susceptibility rather than a TED-specific effect, and higher colocalization probabilities do not establish one.

For *IGF1R*, blood expression-related genetic evidence was limited: ORs above 1, two nominal associations, imprecise ophthalmopathy estimates and mostly weak regional disease signals do not establish a susceptibility effect. Some single-variant omissions removed nominal significance, although most left ORs close to the full estimates. The near-null FinnGen GD estimate argues against large, but not modest, blood-expression effects. Thyroid and cultured-fibroblast *IGF1R* estimates were non-significant and directions were inconsistent across datasets, but single instruments, wide intervals and absent orbital eQTLs leave tissue-specific effects open. Cultured fibroblasts are not orbital fibroblasts, and unestimable combinations are not negative findings. Weaker genetic evidence for *IGF1R* does not argue against IGF-1R blockade: teprotumumab efficacy rests on a randomized trial in active TED [3], and experimental TSHR–IGF-1R crosstalk [5] concerns receptor signaling rather than inherited susceptibility.

Frequency substitution preserved the principal findings and changed full-screen classifications only for exploratory *C4B* and *HLA-G*; like instrument omission, it reused the same data and is not independent replication. For research on TSHR and IGF-1R, blood expression-related associations are more consistently supported for *TSHR*, while its tissue mechanism and the extent of an inherited *IGF1R* contribution remain unresolved. Independent expression instruments and tissue-relevant regulatory data would help address these uncertainties. Comparing GD patients with and without TED, while accounting for follow-up and selection into the GD group, would address orbital susceptibility directly. The present results inform interpretation and follow-up research, but do not guide treatment selection.

### Limitations

The principal limitations involve the exposure and phenotype definitions. Predominantly European blood eQTL effects were applied to BBJ, requiring expression associations and LD patterns to transfer across ancestries. UKB self-reported hyperthyroidism is broader than GD, and an exclusively European source sample was not established. FinnGen's overlapping, population-controlled endpoints do not directly compare GD patients with and without TED; cohort differences and phenotype differences cannot be separated cleanly. Participant overlap with eQTLGen could not be quantified. Cohort-frequency replacement, using a release excluding Framingham participants, addresses only part of expression reconstruction and does not remove LD-reference limitations or independently validate rare-variant effects. The orbital series, from inactive TED, has a single control and no inferential role. The GTEx extension, also predominantly European, had four estimable combinations, each with one instrument. It reused the disease datasets and lacked tissue-specific regional colocalization.

The analytical limitations include post hoc selection, few instruments and a single-causal-variant colocalization model. *TSHR* has one retained European-reference instrument, precluding multi-instrument pleiotropy tests; several other genes have few variants or heterogeneous estimates. MHC analyses are particularly vulnerable to complex LD and correlated expression, and low PP.H4 can reflect model mismatch or limited signal coverage. Multi-signal colocalization [33], directionality testing and formal between-phenotype effect tests were not performed. Only 24 selected genes underwent colocalization; therefore its negative classifications cannot be generalized to the entire screen. Discovery and focused follow-up thresholds address specified families of tests, while the expanded scans and sensitivity analyses remain exploratory. At the discovery threshold, most BBJ genes lacked 80% power for an OR of 1.5, and neither a non-significant MR result nor a posterior below 0.80 establishes absence of an effect.

### Conclusions

Genetically proxied blood *TSHR* expression showed strong, directionally consistent associations with lower odds of GD and related phenotypes, whereas *IGF1R* evidence was weaker, with the most precise estimate near the null. Colocalization supported a shared *TSHR* variant in FinnGen ophthalmopathy and, depending on the prior, in BBJ and strict ophthalmopathy, but not in UKB or FinnGen GD, so these data establish association rather than expression-mediated causality. Weaker genetic evidence for *IGF1R* does not argue against the therapeutic relevance of IGF-1R blockade. The population-controlled ophthalmopathy analyses do not estimate TED susceptibility among patients with GD or predict treatment response.

## Declarations

**Funding.** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Conflict of interest.** The authors declare that they have no conflict of interest.

**Ethics approval.** This study used publicly available, de-identified summary statistics and an institutionally approved orbital tissue dataset. The in-house orbital transcriptomic component was approved by the Institutional Review Board of Pusan National University Hospital (approval number 2104-018-102) and was conducted in accordance with the Declaration of Helsinki.

**Informed consent.** Written informed consent was obtained from all individual participants included in the in-house orbital tissue study.

**Data availability.** Instruments, primary MR and colocalization results, and sensitivity results are supplied as Supplementary Dataset. Analysis code is available from the corresponding author on reasonable request. Public summary statistics analyzed in this study are available from their original repositories: blood *cis*-eQTL data from the eQTLGen Consortium, tissue eQTL data from GTEx v8, the Biobank Japan Graves disease and UK Biobank hyperthyroidism genome-wide association statistics through the GWAS Catalog, and the FinnGen Release 12 GD and ophthalmopathy statistics from FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions.

**Author contributions.** J.P. conceived and designed the study, performed the analyses, and drafted the manuscript. M.-S.K. contributed to data collection and interpretation. K.-H.S. and S.-W.Y. supervised the study and revised the manuscript.

**Acknowledgements.** We thank the eQTLGen Consortium, the GTEx Consortium, Biobank Japan, the UK Biobank, and the FinnGen study and its participants for making their summary statistics publicly available. We acknowledge the GWAS Catalog for hosting and distributing the genome-wide association summary statistics used as outcome data.

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

11. Zhang H, Jiang S, Zhu T, Liu Y, Li J, Fang S, Li Y, Sun J, He X, Tong C, et al. Thyroid-stimulating hormone receptor mediates peripheral-central neuroimmune crosstalk in autoimmune thyroid diseases. BMC Med. 2026;24:404. doi:10.1186/s12916-026-04957-y.

12. Skrivankova VW, Richmond RC, Woolf BAR, Yarmolinsky J, Davies NM, Swanson SA, VanderWeele TJ, Higgins JPT, Timpson NJ, Dimou N, et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: The STROBE-MR Statement. JAMA. 2021;326:1614-1621. doi:10.1001/jama.2021.18236.

13. Võsa U, Claringbould A, Westra HJ, Bonder MJ, Deelen P, Zeng B, Kirsten H, Saha A, Kreuzhuber R, Yazar S, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. Nat Genet. 2021;53:1300-1310. doi:10.1038/s41588-021-00913-z.

14. Sakaue S, Kanai M, Tanigawa Y, Karjalainen J, Kurki M, Koshiba S, Narita A, Konuma T, Yamamoto K, Akiyama M, et al. A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet. 2021;53:1415-1424. doi:10.1038/s41588-021-00931-x.

15. Dönertaş HM, Fabian DK, Valenzuela MF, Partridge L, Thornton JM. Common genetic associations between age-related diseases. Nat Aging. 2021;1:400-412. doi:10.1038/s43587-021-00051-5.

16. Kurki MI, Karjalainen J, Palta P, Sipilä TP, Kristiansson K, Donner KM, Reeve MP, Laivuori H, Aavikko M, Kaunisto MA, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature. 2023;613:508-518. doi:10.1038/s41586-022-05473-8.

17. Finan C, Gaulton A, Kruger FA, Lumbers RT, Shah T, Engmann J, Galver L, Kelley R, Karlsson A, Santos R, et al. The druggable genome and support for target identification and validation in drug development. Sci Transl Med. 2017;9:eaag1166. doi:10.1126/scitranslmed.aag1166.

18. Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015;4:7. doi:10.1186/s13742-015-0047-8.

19. 1000 Genomes Project Consortium, Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, et al. A global reference for human genetic variation. Nature. 2015;526:68-74. doi:10.1038/nature15393.

20. Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, Montgomery GW, Goddard ME, Wray NR, Visscher PM, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016;48:481-487. doi:10.1038/ng.3538.

21. Hemani G, Zheng J, Elsworth B, Wade KH, Haberland V, Baird D, Laurin C, Burgess S, Bowden J, Langdon R, et al. The MR-Base platform supports systematic causal inference across the human phenome. Elife. 2018;7:e34408. doi:10.7554/elife.34408.

22. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. Int J Epidemiol. 2015;44:512-525. doi:10.1093/ije/dyv080.

23. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent Estimation in Mendelian Randomization with Some Invalid Instruments Using a Weighted Median Estimator. Genet Epidemiol. 2016;40:304-314. doi:10.1002/gepi.21965.

24. Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLoS Genet. 2020;16:e1008720. doi:10.1371/journal.pgen.1008720.

25. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369:1318-1330. doi:10.1126/science.aaz1776.

26. Ueda H, Howson JM, Esposito L, Heward J, Snook H, Chamberlain G, Rainbow DB, Hunter KM, Smith AN, Di Genova G, et al. Association of the T-cell regulatory gene CTLA4 with susceptibility to autoimmune disease. Nature. 2003;423:506-511. doi:10.1038/nature01621.

27. Colobran R, Armengol Mdel P, Faner R, Gärtner M, Tykocinski LO, Lucas A, Ruiz M, Juan M, Kyewski B, Pujol-Borrell R. Association of an SNP with intrathymic transcription of TSHR and Graves' disease: a role for defective thymic tolerance. Hum Mol Genet. 2011;20:3415-3423. doi:10.1093/hmg/ddr247.

28. Marín-Sánchez A, Álvarez-Sierra D, González O, Lucas-Martin A, Sellés-Sánchez A, Rudilla F, Enrich E, Colobran R, Pujol-Borrell R. Regulation of TSHR Expression in the Thyroid and Thymus May Contribute to TSHR Tolerance Failure in Graves' Disease Patients via Two Distinct Mechanisms. Front Immunol. 2019;10:1695. doi:10.3389/fimmu.2019.01695.

29. Chen CR, Pichurin P, Nagayama Y, Latrofa F, Rapoport B, McLachlan SM. The thyrotropin receptor autoantigen in Graves disease is the culprit as well as the victim. J Clin Invest. 2003;111:1897-1904. doi:10.1172/JCI17069.

30. Berchner-Pfannschmidt U, Moshkelgosha S, Diaz-Cano S, Edelmann B, Görtz GE, Horstmann M, Noble A, Hansen W, Eckstein A, Banga JP. Comparative Assessment of Female Mouse Model of Graves' Orbitopathy Under Different Environments, Accompanied by Proinflammatory Cytokine and T-Cell Responses to Thyrotropin Hormone Receptor Antigen. Endocrinology. 2016;157:1673-1682. doi:10.1210/en.2015-1829.

31. Daser A, Gulbins A, Horstmann M, Kalisch J, Nolte I, Hashemi Arani F, Lang S, Mattheis S, Bechrakis NE, Banga JP, et al. TNFα drives disease-specific chemokine networks and orbital tissue remodeling in a preclinical model of thyroid eye disease. Eur Thyroid J. 2026;15:ETJ250407. doi:10.1530/ETJ-25-0407.

32. Eckstein AK, Lösch C, Glowacka D, Schott M, Mann K, Esser J, Morgenthaler NG. Euthyroid and primarily hypothyroid patients develop milder and significantly more asymmetrical Graves ophthalmopathy. Br J Ophthalmol. 2009;93:1052-1056. doi:10.1136/bjo.2007.137265.

33. Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genet. 2021;17:e1009440. doi:10.1371/journal.pgen.1009440.
