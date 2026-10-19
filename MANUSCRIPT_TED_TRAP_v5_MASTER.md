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

This two-sample MR study used public association summary statistics and followed STROBE-MR reporting guidance [12]. Exposure data were whole-blood expression quantitative trait locus (eQTL; a variant associated with gene expression) statistics from predominantly European eQTLGen participants [13]. The original outcomes were BBJ GD, UKB self-reported hyperthyroidism or thyrotoxicosis, and Graves ophthalmopathy in FinnGen release 12 (FinnGen ophthalmopathy) [14, 15, 16]. We later added two FinnGen endpoints: E4_GRAVES_STRICT (FinnGen GD) and E4_GRAVES_OPHT_STRICT (strict ophthalmopathy). Table 1 gives case and control counts. UKB ancestry was not established as exclusively European.

BBJ was the discovery dataset, FinnGen GD a same-disease comparison, UKB a broader-phenotype comparison, FinnGen ophthalmopathy a TED-enriched comparison and strict ophthalmopathy a definition sensitivity analysis. The three FinnGen endpoints share participants and were neither treated as independent replications nor meta-analyzed. Their controls were all other participants, not GD patients without TED. The FinnGen extensions, expanded candidate assessment and instrument diagnostics were post hoc; the study was not prospectively registered. Analyses used R 4.3.3, TwoSampleMR 0.7.4, coloc 5.2.3 and PLINK 1.9, with independent numerical verification in Python.

### Expression instruments and MR

For 4,462 potentially druggable genes [17], instruments (variants used to represent a gene's expression) were blood eQTLs with *P* < 5×10^−8^ within 1 Mb of the gene (cis-eQTLs). PLINK clumping kept only uncorrelated variants (linkage disequilibrium [LD] r² < 0.001 within 10 Mb, in 503 European participants from 1000 Genomes Phase 3) [18, 19]. This gave 6,135 single-nucleotide polymorphism (SNP)–gene pairs for 2,544 genes (minimum F statistic 29.7; F > 10 is the usual threshold for adequate instrument strength), and the extended major histocompatibility complex region was not excluded. Variant effects on expression, in standard deviation units, were derived from eQTLGen Z statistics, sample sizes and 1000 Genomes European allele frequencies [20]. TwoSampleMR aligned exposure and disease effects to the same allele, excluding unavailable and strand-ambiguous (palindromic) variants. UKB linear-model estimates were rescaled using the case proportion to approximate the log-odds scale of BBJ and FinnGen (Supplementary Methods).

For each gene and dataset, MR used the Wald ratio (disease effect divided by expression effect) for one instrument and, for several, multiplicative random-effects inverse-variance weighting (IVW), a precision-weighted average [21]. Estimates are ORs per approximately one standard deviation higher genetically proxied blood expression, with 95% CIs and two-sided *P* values. The BBJ screen threshold was Bonferroni-corrected for all 2,544 instrumented genes, estimable or not (*P* < 0.05/2,544 = 1.965×10^−5^); it was applied descriptively to the other datasets and in an exploratory full FinnGen GD screen. The 13 BBJ discoveries were assessed at a focused FinnGen GD threshold (*P* < 0.05/13 and the same direction). Nominal (uncorrected) *P* < 0.05 was reported descriptively. ORs detectable with 80% power were approximated from observed standard errors (Supplementary Methods).

MR assumes that instruments predict expression, are unrelated to confounders of expression and disease, and affect disease only through that gene's expression.

For estimates with several instruments, Cochran's Q test assessed heterogeneity among instruments and, for selected genes, post hoc omission of each instrument in turn (leave-one-out) assessed dependence on single variants. MR-Egger regression, whose intercept tests whether instruments affect disease through other pathways (directional pleiotropy), and weighted-median and weighted-mode estimators [22, 23] were applied only in the original three datasets (Supplementary Dataset, MR_estimators_original). Analyses were repeated with exposure effects reconstructed from eQTLGen cohort allele frequencies; for *C4B* and *HLA-G*, whose FinnGen GD screen result changed, instrument weights and omissions were compared between frequency sources (Supplementary Methods).

In a post hoc audit, rs179252 and rs179247 (the variant examined in published thymic *TSHR* studies) were re-extracted from the original expression and disease files. Disease effects were aligned to the blood expression-increasing allele to give per-allele ORs and recalculate the MR estimates, and the correlation between the two variants' G-allele counts was calculated in European and East Asian 1000 Genomes genotypes. The audit checked allele direction, not thymic mediation (Supplementary Dataset, TSHR_allele_audit).

### Colocalization

MR associations are the main findings, and colocalization was used to interpret them; the original analysis reported both in parallel, and this emphasis was adopted post hoc without prospective registration. We used coloc.abf, which assumes at most one causal variant per trait, with all available eQTLGen cis-eQTL statistics for each gene, significant or not, matched to each disease dataset [8]. It was applied to 24 genes: the 13 BBJ discoveries, *IGF1R*, the nine additional FinnGen GD screen discoveries and *HLA-G*, added because its MR result changed with the allele-frequency source. Posterior probabilities (PP.H0–PP.H4) were estimated for five hypotheses: neither trait associated (H0), disease only (H1), expression only (H2), distinct causal variants (H3) and one shared causal variant (H4). Prior probabilities that a variant is associated with disease (p1), expression (p2) or both (p12) were 10^−4^, 10^−4^ and 10^−5^, with p12 also set to 5×10^−6^ and 10^−6^ [24]. PP.H4 ≥ 0.80 was a descriptive reporting threshold.

For *TSHR*, colocalization was repeated in narrower regions and in a variant set common to all five datasets, and we checked whether the matched variants included each dataset's strongest regional disease variant (Supplementary Methods).

### Tissue-specific extension

In a post hoc extension, we assessed *TSHR* and *IGF1R* in a fixed panel of five Genotype-Tissue Expression (GTEx) v8 tissues: thyroid, cultured skin fibroblasts, subcutaneous adipose, visceral adipose and whole blood [25]. The plan was fixed after the blood results were known but before tissue results were retrieved. GTEx donors were predominantly European. Instruments required cis-eQTL *P* < 5×10^−8^, F > 10 and European-reference LD r² < 0.001. Ten planned gene–tissue comparisons against FinnGen GD formed the primary family (*P* < 0.005); 40 comparisons against the other four datasets formed a secondary family (*P* < 0.00125). Unestimable combinations remained in the correction denominators. Effects used GTEx slopes and standard errors for inverse-normal-transformed expression. Tissue colocalization was not performed because complete regional statistics were not obtained.

### Descriptive orbital expression

In-house RNA sequencing of orbital adipose tissue from four patients with inactive TED and of an upper-eyelid nasal fat pad from one control without thyroid disease was descriptive only (Supplementary Methods). Transcript abundance was displayed without significance testing (Figure S1) and was not used to select genes or validate the genetic findings.

## Results

### Discovery across druggable genes

Of the 2,544 instrumented genes, 2,234 were estimable in BBJ, 2,505 in UKB and 2,480 in each FinnGen endpoint. Thirteen genes passed the BBJ screen threshold (*P* < 0.05/2,544), including *TSHR* and *CTLA4* but not *IGF1R* (Figure 1). At this threshold, 14.6% and 35.6% of estimable BBJ genes had 80% power to detect ORs of 1.5 and 2.0, respectively (Supplementary Dataset, Detection_original).

### TSHR and IGF1R across disease datasets

Higher genetically proxied blood *TSHR* expression was associated with lower odds of GD in BBJ (OR 0.12, 95% CI 0.07–0.21; *P* = 1.09×10^−14^) and FinnGen (OR 0.08, 95% CI 0.05–0.12; *P* = 3.19×10^−32^). ORs were also below 1 for UKB hyperthyroidism and both ophthalmopathy endpoints, and all estimates except strict ophthalmopathy passed the screen threshold (Table 2; Figure 2). Each estimate used the single instrument rs179252. Its expression-increasing G allele had per-allele disease ORs of 0.77–0.82 across datasets and was correlated with rs179247-G (European r = 0.978; East Asian r = 0.996), which also increased blood *TSHR* expression. Allele-aligned recalculation reproduced the MR estimates (Supplementary Dataset, TSHR_allele_audit).

All five *IGF1R* ORs exceeded 1, and none passed the screen threshold. Nominal associations (uncorrected *P* < 0.05) were found only in BBJ GD (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB hyperthyroidism (OR 1.35, 95% CI 1.07–1.70; *P* = 0.0117). The most precise estimate, FinnGen GD, was OR 1.07 (95% CI 0.86–1.33; *P* = 0.539); this comparison had 80% power to detect an OR of 1.36 at *P* < 0.05 (Supplementary Methods). Weighted-median ORs were 1.56 (*P* = 0.0288) in BBJ, 1.42 (*P* = 0.00244) in UKB and 1.46 (*P* = 0.124) in FinnGen ophthalmopathy. In these three datasets, MR-Egger and weighted-mode estimates had *P* ≥ 0.32, and MR-Egger intercept *P* values were 0.66, 0.56 and 0.43, respectively (Supplementary Dataset, MR_estimators_original).

### Other associated genes

Four genes passed the screen threshold in both GD datasets: *TSHR*, *CTLA4*, *HLA-DQA2* and *TUBB* (Table 3). *CTLA4* was associated with lower odds, and *HLA-DQA2* and *TUBB*, both in the major histocompatibility complex (MHC) region, with higher odds. *CTLA4* and *TUBB* had direction-consistent nominal associations in all five datasets; *HLA-DQA2* did in four, but its two-instrument UKB estimate was non-significant and below 1 (Supplementary Dataset, MR_primary). At the focused FinnGen GD threshold (*P* < 0.05/13 = 0.00385), seven of the 13 BBJ discoveries had direction-consistent associations: these four genes plus *C4A*, *TNFSF14* and *IFNGR1*. *HLA-A*, which had the smallest BBJ *P* value, had a non-significant FinnGen GD estimate.

The exploratory full FinnGen GD screen identified 13 genes, nine of them outside the BBJ discovery set (Table S1). Seven of the nine were instrumented in the extended MHC region, and two pairs each shared one instrument.

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

Higher genetically proxied blood *TSHR* expression was associated with lower GD odds in BBJ and FinnGen GD, with concordant directions for UKB hyperthyroidism and both ophthalmopathy endpoints. All five *IGF1R* ORs exceeded 1, but only BBJ and UKB were nominally significant, none passed the screen threshold, and FinnGen GD, the most precise, lay close to the null. Under the primary prior, colocalization favored a shared *TSHR* variant except in UKB and FinnGen GD. Significance for one gene but not the other is not a formal test of different effects.

Ji and colleagues linked blood *TSHR* expression to GD and strict ophthalmopathy in an earlier FinnGen release [9], and integrative work also investigated *TSHR* [11]. We add an *IGF1R* comparison, GD in two populations and a druggable-gene screen. Both studies used eQTLGen and FinnGen, so agreement is not fully independent. Ji and colleagues accepted PP.H4 above 0.50 [9], whereas we used PP.H4 ≥ 0.80 with p12 = 10^−5^; with different methods and releases, the labels are not directly comparable and differences are not contradictions.

The multi-instrument FinnGen and UKB *CTLA4* estimates relied heavily on rs13030124; with that caveat, its direction-consistent associations across all five datasets support it as an autoimmune comparator [26]. *HLA-A*, the strongest BBJ association, was non-significant in FinnGen GD, so the smallest discovery *P* value need not transfer best. *HLA-DQA2* was discordant only in its one two-instrument estimate, UKB, whose instruments were strongly heterogeneous (Supplementary Dataset, MR_estimators_original). With shared instruments, MHC linkage, single-variant dependence of the FinnGen GD *TUBB* and *C4A* estimates, one regional association behind the three chromosome 16p11.2 genes, and no shared-variant support for *HLA-DQA2* or *TUBB*, these associations cannot each be counted as an independent causal target.

The lower-odds direction need not contradict stimulating-antibody pathogenesis: blood transcript abundance, tissue expression and receptor activation are different exposures with tissue-dependent regulatory effects [25]. The expression-increasing rs179252-G allele is linked to rs179247-G, the GD-protective allele associated with greater intrathymic *TSHR* transcription [27, 28]; this fits central tolerance, though the thymic step was untested and neither variant is a proven causal regulator. Each rs179252-G allele raises blood expression by 0.105 standard deviation, and its per-allele disease ORs were 0.80 (BBJ), 0.77 (UKB), 0.78 (FinnGen ophthalmopathy), 0.77 (FinnGen GD) and 0.82 (strict ophthalmopathy); rescaled to one standard deviation, these moderate effects give the reported odds ratios, which describe lifelong inherited differences, not raised *TSHR* expression or TSHR-directed treatment.

Very small *TSHR* MR *P* values in UKB and FinnGen GD coexisted with distinct-variant patterns, as for *CTLA4* in BBJ and FinnGen GD, so association strength and shared-variant evidence are not interchangeable: a non-zero MR estimate can arise from distinct but correlated variants [29]. These two datasets had the most cases, and their strongest regional disease variant was more strongly associated than rs179252 (Table S3); regional checks excluded a missing-peak explanation but did not identify the causal variant. Support in BBJ and strict ophthalmopathy fell below 0.80 only at the lowest prior tested (Figure 3B). The adipose instruments are themselves blood *TSHR* eQTLs, so their estimates re-express the regional association and cannot resolve the pattern or establish an adipose mechanism. These data establish association, not expression-mediated causality.

In mice, the shed TSHR A-subunit induced more thyroid-stimulating autoimmunity than a minimally cleaving receptor [30]; plasmid immunization reproduced orbital adipogenesis, fibrosis and muscle pathology [31], and orbital fibroblasts from immunized mice showed enhanced inflammatory and hyaluronan responses to TNFα [32]. Ophthalmopathy also occurs in euthyroid and hypothyroid patients [33]. Our instruments measure neither receptor cleavage nor antibody activity. Because most patients with TED have GD, the controls included people with GD but no recorded eye disease, and the ophthalmopathy per-allele ORs resembled those for GD (Table S3), these estimates may mostly reflect GD susceptibility; their higher PP.H4 values, from the two smallest, weakest-signal datasets, do not establish a TED-specific effect.

For *IGF1R*, blood expression-related evidence was limited: colocalization mostly favored an expression-only signal, consistent with weak regional disease signals, and some single-variant omissions removed nominal significance, though most left ORs close to the full estimates (Figure S2). The near-null FinnGen GD estimate argues against large, but not modest, effects; non-significance is not evidence of no inherited role. Single-instrument tissue estimates and absent orbital eQTL data leave tissue effects open. Weaker genetic evidence for *IGF1R* does not argue against IGF-1R blockade: teprotumumab reduced proptosis in a randomized trial in active TED [3] and in long-duration, low-activity TED [34], and TSHR–IGF-1R crosstalk [5] concerns receptor signaling, not inherited susceptibility.

These results inform research rather than treatment selection, leaving unresolved the tissue mediating the *TSHR* association and any inherited *IGF1R* contribution, which independent instruments, thymic and orbital regulatory data and multi-signal colocalization could address. Orbital susceptibility requires comparing GD patients with and without TED, accounting for follow-up and selection into GD.

### Limitations

BBJ analyses assumed that predominantly European eQTLGen effects and LD transfer to East Asians, in both the MR and colocalization. UKB self-reported hyperthyroidism is broader than GD and was not established as exclusively European. The three overlapping, population-controlled FinnGen endpoints are not independent replications, do not compare GD patients with and without TED, and cannot separate cohort from phenotype differences; participant overlap with eQTLGen was unquantified. Blood is neither thymus nor orbit: the predominantly European GTEx extension had four single-instrument gene–tissue combinations, reused the disease data and lacked colocalization, and adipose and cultured fibroblasts are not orbital tissue, with 30 unestimable comparisons that are not negative findings. The orbital series (inactive TED; one control from a different fat depot) is descriptive only.

The FinnGen extensions, candidate expansion and colocalization emphasis were post hoc, and only the screen, focused and tissue-family thresholds address specified test families. The single *TSHR* instrument precludes pleiotropy and omission tests and cannot exclude neighboring or correlated variants; other genes have few or heterogeneous instruments, MR-Egger has little power with three or four *IGF1R* instruments, and frequency substitution and omission reused the same data. Under the single-causal-variant model, a low PP.H4 can reflect model mismatch or complex MHC LD; multi-signal colocalization [35], directionality and formal between-phenotype or between-gene tests were not performed, and only 24 genes were assessed. At the screen threshold, most BBJ genes lacked 80% power for an OR of 1.5, so a non-significant MR result, a PP.H4 below 0.80 or an unestimable comparison does not establish absence of an effect.

### Conclusions

Higher genetically proxied blood *TSHR* expression was associated with lower odds of GD, hyperthyroidism and ophthalmopathy in every dataset. *IGF1R* odds ratios exceeded 1 in every dataset but nominally significant only in BBJ and UKB, and the most precise estimate, FinnGen GD, lay close to the null. Colocalization supported a shared *TSHR* variant in FinnGen ophthalmopathy and, except at the lowest prior tested, in BBJ and strict ophthalmopathy, but favored distinct variants in UKB and FinnGen GD. These analyses addressed disease susceptibility, not treatment response, and the *IGF1R* evidence does not argue against IGF-1R blockade. The ophthalmopathy analyses used population controls and did not estimate TED susceptibility among patients with GD.

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

29. Zuber V, Grinberg NF, Gill D, Manipur I, Slob EAW, Patel A, Wallace C, Burgess S. Combining evidence from Mendelian randomization and colocalization: Review and comparison of approaches. Am J Hum Genet. 2022;109:767-782. doi:10.1016/j.ajhg.2022.04.001.

30. Chen CR, Pichurin P, Nagayama Y, Latrofa F, Rapoport B, McLachlan SM. The thyrotropin receptor autoantigen in Graves disease is the culprit as well as the victim. J Clin Invest. 2003;111:1897-1904. doi:10.1172/JCI17069.

31. Berchner-Pfannschmidt U, Moshkelgosha S, Diaz-Cano S, Edelmann B, Görtz GE, Horstmann M, Noble A, Hansen W, Eckstein A, Banga JP. Comparative Assessment of Female Mouse Model of Graves' Orbitopathy Under Different Environments, Accompanied by Proinflammatory Cytokine and T-Cell Responses to Thyrotropin Hormone Receptor Antigen. Endocrinology. 2016;157:1673-1682. doi:10.1210/en.2015-1829.

32. Daser A, Gulbins A, Horstmann M, Kalisch J, Nolte I, Hashemi Arani F, Lang S, Mattheis S, Bechrakis NE, Banga JP, et al. TNFα drives disease-specific chemokine networks and orbital tissue remodeling in a preclinical model of thyroid eye disease. Eur Thyroid J. 2026;15:ETJ250407. doi:10.1530/ETJ-25-0407.

33. Eckstein AK, Lösch C, Glowacka D, Schott M, Mann K, Esser J, Morgenthaler NG. Euthyroid and primarily hypothyroid patients develop milder and significantly more asymmetrical Graves ophthalmopathy. Br J Ophthalmol. 2009;93:1052-1056. doi:10.1136/bjo.2007.137265.

34. Douglas RS, Couch S, Wester ST, Fowler BT, Liu CY, Subramanian PS, Tang R, Nguyen QT, Maamari RN, Ugradar S, et al. Efficacy and Safety of Teprotumumab in Patients With Thyroid Eye Disease of Long Duration and Low Disease Activity. J Clin Endocrinol Metab. 2024;109:25-35. doi:10.1210/clinem/dgad637.

35. Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genet. 2021;17:e1009440. doi:10.1371/journal.pgen.1009440.
