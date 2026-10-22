# TSHR and IGF1R in Graves disease and thyroid eye disease: genetic evidence from blood expression

Jungyul Park¹†, Min-Seon Kim²†, Suk-Woo Yang¹\*

¹ Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
² Department of Ophthalmology, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  

† These authors contributed equally to this work as first authors.  
\* Corresponding author.

**Corresponding author:** Suk-Woo Yang, MD, PhD, Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea. E-mail: yswoph@catholic.ac.kr; Tel: +82-2-2258-2847.  



## Abstract

**Objective:** Thyroid eye disease (TED) can accompany Graves disease (GD), an autoimmune cause of hyperthyroidism. We compared genetic evidence linking blood *TSHR* and *IGF1R* expression with susceptibility.

**Methods:** We screened 2,544 druggable genes using Mendelian randomization (MR). Disease outcomes were Biobank Japan (BBJ) GD, UK Biobank (UKB) hyperthyroidism and overlapping FinnGen GD, Graves ophthalmopathy (GO) and strict GO endpoints. Colocalization tested 24 genes for a shared causal variant. Tissue eQTL and orbital RNA-sequencing (four TED, one control) were supplementary.

**Results:** Thirteen genes passed BBJ discovery (*P* < 0.05/2,544). Higher genetically proxied *TSHR* expression was associated with lower GD odds in BBJ (odds ratio [OR] per approximately one standard deviation 0.12, 95% confidence interval [CI] 0.07–0.21; *P* = 1.09×10^−14^) and FinnGen (OR 0.08, 95% CI 0.05–0.12; *P* = 3.19×10^−32^), concordantly across all five datasets. Shared-variant probability exceeded 0.80 in FinnGen GO and, prior-sensitively, BBJ and strict GO, but UKB (0.226) and FinnGen GD (0.264) favored distinct variants. *IGF1R* ORs exceeded 1 but reached *P* < 0.05 only in BBJ (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB (*P* = 0.0117); the most precise estimate, FinnGen GD, was OR 1.07 (95% CI 0.86–1.33). Single-variant tissue estimates added no independent evidence. Mean orbital *TSHR* and *IGF1R* transcript levels exceeded the control.

**Conclusions:** Association evidence was stronger for *TSHR* than *IGF1R*, but absent colocalization in UKB and FinnGen GD leaves expression-mediated causality unproven. Population-controlled data cannot isolate TED susceptibility within GD; weaker *IGF1R* evidence does not argue against IGF-1R–targeted treatment.

**Keywords:** Graves disease; thyroid eye disease; *TSHR*; *IGF1R*; Mendelian randomization; colocalization

## Introduction

Thyroid eye disease (TED), or Graves ophthalmopathy, can cause proptosis, diplopia and visual loss in patients with Graves disease (GD) [1, 2]; a meta-analysis estimated that it affects about 40% of them [3]. In GD, autoantibodies that stimulate the thyroid-stimulating hormone receptor (TSHR) cause hyperthyroidism [2]. Common *TSHR* variants are established GD risk factors in European and East Asian populations [4, 5], but no TED-specific risk variant is known [6]. In orbital fibroblasts, TSHR forms a signaling complex with the insulin-like growth factor 1 receptor (IGF-1R) [7, 8], and the IGF-1R–blocking antibody teprotumumab reduced proptosis in randomized trials in active TED [9, 10]. Genetically supported drug targets succeed more often, but most approved drug mechanisms lack such support [11, 12], so an effective target need not be a susceptibility gene. Whether inherited variation in *IGF1R* expression, like *TSHR* variation, is associated with GD or TED is unclear.

Mendelian randomization (MR) treats inherited variants as a natural experiment: because alleles are allocated at conception, comparing people whose variants raise or lower a gene's expression is less prone to confounding and reverse causation than other observational comparisons, provided the variants affect disease only through that expression [13]. Their effects are lifelong, unlike a drug's, and can differ between tissues [14]. Here, genetically proxied expression means blood expression predicted from inherited variants rather than measured directly, and the outcomes are diagnoses, not treatment responses. MR reports odds ratios (ORs) with 95% confidence intervals (CIs). Colocalization asks whether expression and disease signals in a region share one causal variant and reports probabilities, not *P* values [15]. Because nearby variants are often inherited together, an MR association can be highly significant even when the two signals come from distinct variants [16, 17].

In FinnGen release 10, Ji and colleagues linked higher genetically proxied blood *TSHR* expression to lower odds of GD and strict ophthalmopathy [18]. Using release 12 endpoints also analyzed here, Zhang and colleagues reported shared-variant support for *TSHR* with the broader ophthalmopathy endpoint but weak support with GD [19]. A proteome-wide MR study of GD and ophthalmopathy prioritized other proteins [20]. None of these studies compared *IGF1R* with *TSHR*, and all relied mainly on FinnGen.

Within an MR screen of potentially druggable genes [21], we therefore compared genetic evidence linking blood *TSHR* and *IGF1R* expression with GD in Biobank Japan (BBJ; discovery) and FinnGen, hyperthyroidism in UK Biobank (UKB) and two ophthalmopathy endpoints in FinnGen. *CTLA4*, an established autoimmune susceptibility gene, served as comparator [5, 22]. The three FinnGen endpoints overlap and are not independent replications. We examined colocalization and sensitivity to single variants, allele frequencies and colocalization assumptions; FinnGen GD, strict ophthalmopathy, tissue-specific and several sensitivity analyses were post hoc. Because the ophthalmopathy endpoints used population controls, they cannot estimate TED risk among patients with GD, and no analysis predicts treatment response. Orbital RNA-sequencing data (four TED samples, one control) are descriptive only.

## Methods

### Study design and disease datasets

This two-sample MR study used public association summary statistics and followed STROBE-MR reporting guidance [23]. Exposure data were whole-blood expression quantitative trait locus (eQTL; a variant associated with gene expression) statistics from predominantly European eQTLGen participants [24]. The original outcomes were BBJ GD, UKB self-reported hyperthyroidism or thyrotoxicosis, and Graves ophthalmopathy in FinnGen release 12 (FinnGen ophthalmopathy) [25, 26, 27]. We later added two FinnGen endpoints: E4_GRAVES_STRICT (FinnGen GD) and E4_GRAVES_OPHT_STRICT (strict ophthalmopathy). Table 1 gives case and control counts. UKB ancestry was not established as exclusively European.

BBJ was the discovery dataset, FinnGen GD a same-disease comparison, UKB a broader-phenotype comparison, FinnGen ophthalmopathy a TED-enriched comparison and strict ophthalmopathy a definition sensitivity analysis. The three FinnGen endpoints share participants and were neither treated as independent replications nor meta-analyzed. Their controls were all other participants, not GD patients without TED. The FinnGen extensions, expanded candidate assessment and instrument diagnostics were post hoc; the study was not prospectively registered. Analyses used R 4.3.3, TwoSampleMR 0.7.4, coloc 5.2.3 and PLINK 1.9, with independent numerical verification in Python.

### Expression instruments and MR

For 4,462 potentially druggable genes [21], instruments (variants used to represent a gene's expression) were blood eQTLs with *P* < 5×10^−8^ within 1 Mb of the gene (cis-eQTLs). PLINK clumping kept only uncorrelated variants (linkage disequilibrium [LD] r² < 0.001 within 10 Mb, in 503 European participants from 1000 Genomes Phase 3) [28, 29]. This gave 6,135 single-nucleotide polymorphism (SNP)–gene pairs for 2,544 genes (minimum F statistic 29.7; F > 10 is the usual threshold for adequate instrument strength), and the extended major histocompatibility complex region was not excluded. Variant effects on expression, in standard deviation units, were derived from eQTLGen Z statistics, sample sizes and 1000 Genomes European allele frequencies [30]. TwoSampleMR aligned exposure and disease effects to the same allele, excluding unavailable and strand-ambiguous (palindromic) variants. UKB linear-model estimates were rescaled using the case proportion to approximate the log-odds scale of BBJ and FinnGen [31] (Supplementary Methods).

For each gene and dataset, MR used the Wald ratio (disease effect divided by expression effect) for one instrument and, for several, multiplicative random-effects inverse-variance weighting (IVW), a precision-weighted average [32]. Estimates are ORs per approximately one standard deviation higher genetically proxied blood expression, with 95% CIs and two-sided *P* values. The BBJ screen threshold was Bonferroni-corrected for all 2,544 instrumented genes, estimable or not (*P* < 0.05/2,544 = 1.965×10^−5^); it was applied descriptively to the other datasets and in an exploratory full FinnGen GD screen. The 13 BBJ discoveries were assessed at a focused FinnGen GD threshold (*P* < 0.05/13 and the same direction). Nominal (uncorrected) *P* < 0.05 was reported descriptively. ORs detectable with 80% power were approximated from observed standard errors (Supplementary Methods).

MR assumes that instruments predict expression, are unrelated to confounders of expression and disease, and affect disease only through that gene's expression.

For estimates with several instruments, Cochran's Q test assessed heterogeneity among instruments and, for selected genes, post hoc omission of each instrument in turn (leave-one-out) assessed dependence on single variants. MR-Egger regression, whose intercept tests whether instruments affect disease through other pathways (directional pleiotropy), and weighted-median and weighted-mode estimators [33, 34, 35] were applied only in the original three datasets (Supplementary Dataset, MR_estimators_original). Analyses were repeated with exposure effects reconstructed from eQTLGen cohort allele frequencies; for *C4B* and *HLA-G*, whose FinnGen GD screen result changed, instrument weights and omissions were compared between frequency sources (Supplementary Methods).

In a post hoc audit, rs179252 and rs179247 (the variant examined in published thymic *TSHR* studies [36, 37]) were re-extracted from the original expression and disease files. Disease effects were aligned to the blood expression-increasing allele to give per-allele ORs and recalculate the MR estimates, and the correlation between the two variants' G-allele counts was calculated in European and East Asian 1000 Genomes genotypes. The audit checked allele direction, not thymic mediation (Supplementary Dataset, TSHR_allele_audit).

### Colocalization

MR associations are the main findings, and colocalization was used to interpret them; the original analysis reported both in parallel, and this emphasis was adopted post hoc without prospective registration. We used coloc.abf, which assumes at most one causal variant per trait, with all available eQTLGen cis-eQTL statistics for each gene, significant or not, matched to each disease dataset [15]. It was applied to 24 genes: the 13 BBJ discoveries, *IGF1R*, the nine additional FinnGen GD screen discoveries and *HLA-G*, added because its MR result changed with the allele-frequency source. Posterior probabilities (PP.H0–PP.H4) were estimated for five hypotheses: neither trait associated (H0), disease only (H1), expression only (H2), distinct causal variants (H3) and one shared causal variant (H4). Prior probabilities that a variant is associated with disease (p1), expression (p2) or both (p12) were 10^−4^, 10^−4^ and 10^−5^, with p12 also set to 5×10^−6^ and 10^−6^ [38]. PP.H4 ≥ 0.80 was a descriptive reporting threshold.

For *TSHR*, colocalization was repeated in narrower regions and in a variant set common to all five datasets, and we checked whether the matched variants included each dataset's strongest regional disease variant (Supplementary Methods).

### Tissue-specific extension

In a post hoc extension, we assessed *TSHR* and *IGF1R* in a fixed panel of five Genotype-Tissue Expression (GTEx) v8 tissues: thyroid, cultured skin fibroblasts, subcutaneous adipose, visceral adipose and whole blood [39]. The plan was fixed after the blood results were known but before tissue results were retrieved. GTEx donors were predominantly European. Instruments required cis-eQTL *P* < 5×10^−8^, F > 10 and European-reference LD r² < 0.001. Ten planned gene–tissue comparisons against FinnGen GD formed the primary family (*P* < 0.005); 40 comparisons against the other four datasets formed a secondary family (*P* < 0.00125). Unestimable combinations remained in the correction denominators. Effects used GTEx slopes and standard errors for inverse-normal-transformed expression. Tissue colocalization was not performed because complete regional statistics were not obtained.

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

Ji and colleagues linked blood *TSHR* expression to GD and strict ophthalmopathy in an earlier FinnGen release [18], and integrative work also investigated *TSHR* [19]. We add an *IGF1R* comparison, GD in two populations and a druggable-gene screen. Both studies used eQTLGen and FinnGen, so agreement is not fully independent. Ji and colleagues accepted PP.H4 above 0.50 [18], whereas we used PP.H4 ≥ 0.80 with p12 = 10^−5^; with different methods and releases, the labels are not directly comparable and differences are not contradictions.

The multi-instrument FinnGen and UKB *CTLA4* estimates relied heavily on rs13030124; with that caveat, its direction-consistent associations across all five datasets support it as an autoimmune comparator [22]. *HLA-A*, the strongest BBJ association, was non-significant in FinnGen GD, so the smallest discovery *P* value need not transfer best. *HLA-DQA2* was discordant only in its one two-instrument estimate, UKB, whose instruments were strongly heterogeneous (Supplementary Dataset, MR_estimators_original). With shared instruments, MHC linkage, single-variant dependence of the FinnGen GD *TUBB* and *C4A* estimates, one regional association behind the three chromosome 16p11.2 genes, and no shared-variant support for *HLA-DQA2* or *TUBB*, these associations cannot each be counted as an independent causal target.

The lower-odds direction need not contradict stimulating-antibody pathogenesis: blood transcript abundance, tissue expression and receptor activation are different exposures with tissue-dependent regulatory effects [39]. The expression-increasing rs179252-G allele is linked to rs179247-G, the GD-protective allele associated with greater intrathymic *TSHR* transcription [4, 36, 37]; this fits central tolerance, though the thymic step was untested and neither variant is a proven causal regulator. Each rs179252-G allele raises blood expression by 0.105 standard deviation, and its per-allele disease ORs were 0.80 (BBJ), 0.77 (UKB), 0.78 (FinnGen ophthalmopathy), 0.77 (FinnGen GD) and 0.82 (strict ophthalmopathy); rescaled to one standard deviation, these moderate effects give the reported odds ratios, which describe lifelong inherited differences, not raised *TSHR* expression or TSHR-directed treatment.

Very small *TSHR* MR *P* values in UKB and FinnGen GD coexisted with distinct-variant patterns, as for *CTLA4* in BBJ and FinnGen GD, so association strength and shared-variant evidence are not interchangeable: a non-zero MR estimate can arise from distinct but correlated variants [16]. These two datasets had the most cases, and their strongest regional disease variant was more strongly associated than rs179252 (Table S3); regional checks excluded a missing-peak explanation but did not identify the causal variant. Support in BBJ and strict ophthalmopathy fell below 0.80 only at the lowest prior tested (Figure 3B). The adipose instruments are themselves blood *TSHR* eQTLs, so their estimates re-express the regional association and cannot resolve the pattern or establish an adipose mechanism. These data establish association, not expression-mediated causality.

In mice, the shed TSHR A-subunit induced more thyroid-stimulating autoimmunity than a minimally cleaving receptor [40]; plasmid immunization reproduced orbital adipogenesis, fibrosis and muscle pathology [41], and orbital fibroblasts from immunized mice showed enhanced inflammatory and hyaluronan responses to TNFα [42]. Ophthalmopathy also occurs in euthyroid and hypothyroid patients [43]. Our instruments measure neither receptor cleavage nor antibody activity. Because most patients with TED have GD, the controls included people with GD but no recorded eye disease, and the ophthalmopathy per-allele ORs resembled those for GD (Table S3), these estimates may mostly reflect GD susceptibility; their higher PP.H4 values, from the two smallest, weakest-signal datasets, do not establish a TED-specific effect.

For *IGF1R*, blood expression-related evidence was limited: colocalization mostly favored an expression-only signal, consistent with weak regional disease signals, and some single-variant omissions removed nominal significance, though most left ORs close to the full estimates (Figure S2). The near-null FinnGen GD estimate argues against large, but not modest, effects; non-significance is not evidence of no inherited role. Single-instrument tissue estimates and absent orbital eQTL data leave tissue effects open. Weaker genetic evidence for *IGF1R* does not argue against IGF-1R blockade: teprotumumab reduced proptosis in a randomized trial in active TED [10] and in long-duration, low-activity TED [44], and TSHR–IGF-1R crosstalk [8] concerns receptor signaling, not inherited susceptibility.

These results inform research rather than treatment selection, leaving unresolved the tissue mediating the *TSHR* association and any inherited *IGF1R* contribution, which independent instruments, thymic and orbital regulatory data and multi-signal colocalization could address. Orbital susceptibility requires comparing GD patients with and without TED, accounting for follow-up and selection into GD.

### Limitations

BBJ analyses assumed that predominantly European eQTLGen effects and LD transfer to East Asians, in both the MR and colocalization. UKB self-reported hyperthyroidism is broader than GD and was not established as exclusively European. The three overlapping, population-controlled FinnGen endpoints are not independent replications, do not compare GD patients with and without TED, and cannot separate cohort from phenotype differences; participant overlap with eQTLGen was unquantified. Blood is neither thymus nor orbit: the predominantly European GTEx extension had four single-instrument gene–tissue combinations, reused the disease data and lacked colocalization, and adipose and cultured fibroblasts are not orbital tissue, with 30 unestimable comparisons that are not negative findings. The orbital series (inactive TED; one control from a different fat depot) is descriptive only.

The FinnGen extensions, candidate expansion and colocalization emphasis were post hoc, and only the screen, focused and tissue-family thresholds address specified test families. The single *TSHR* instrument precludes pleiotropy and omission tests and cannot exclude neighboring or correlated variants; other genes have few or heterogeneous instruments, MR-Egger has little power with three or four *IGF1R* instruments, and frequency substitution and omission reused the same data. Under the single-causal-variant model, a low PP.H4 can reflect model mismatch or complex MHC LD; multi-signal colocalization [45], directionality and formal between-phenotype or between-gene tests were not performed, and only 24 genes were assessed. At the screen threshold, most BBJ genes lacked 80% power for an OR of 1.5, so a non-significant MR result, a PP.H4 below 0.80 or an unestimable comparison does not establish absence of an effect.

### Conclusions

Higher genetically proxied blood *TSHR* expression was associated with lower odds of GD, hyperthyroidism and ophthalmopathy in every dataset. *IGF1R* odds ratios exceeded 1 in every dataset but nominally significant only in BBJ and UKB, and the most precise estimate, FinnGen GD, lay close to the null. Colocalization supported a shared *TSHR* variant in FinnGen ophthalmopathy and, except at the lowest prior tested, in BBJ and strict ophthalmopathy, but favored distinct variants in UKB and FinnGen GD. These analyses addressed disease susceptibility, not treatment response, and the *IGF1R* evidence does not argue against IGF-1R blockade. The ophthalmopathy analyses used population controls and did not estimate TED susceptibility among patients with GD.

## Declarations

**Funding.** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Conflict of interest.** The authors declare that they have no conflict of interest.

**Ethics approval.** This study used publicly available, de-identified summary statistics and an institutionally approved orbital tissue dataset. The in-house orbital transcriptomic component was approved by the Institutional Review Board of Pusan National University Hospital (approval number 2104-018-102) and was conducted in accordance with the Declaration of Helsinki.

**Informed consent.** Written informed consent was obtained from all individual participants included in the in-house orbital tissue study.

**Data availability.** Instruments, primary MR and colocalization results, and sensitivity results are supplied as Supplementary Dataset. Analysis code is available from the corresponding author on reasonable request. Public summary statistics analyzed in this study are available from their original repositories: blood *cis*-eQTL data from the eQTLGen Consortium, tissue eQTL data from GTEx v8, the Biobank Japan Graves disease and UK Biobank hyperthyroidism genome-wide association statistics through the GWAS Catalog, and the FinnGen Release 12 GD and ophthalmopathy statistics from FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions.

**Author contributions.** J.P. and M.-S.K. contributed equally to this work; both performed the analyses and drafted the manuscript. J.P. conceived and designed the study. M.-S.K. contributed to data collection and interpretation. S.-W.Y. supervised the study and revised the manuscript.

**Acknowledgements.** We thank Kyung-Hwa Shin (Department of Laboratory Medicine, Pusan National University Hospital) for providing the orbital tissue RNA-sequencing data. We thank the eQTLGen Consortium, the GTEx Consortium, Biobank Japan, the UK Biobank, and the FinnGen study and its participants for making their summary statistics publicly available. We acknowledge the GWAS Catalog for hosting and distributing the genome-wide association summary statistics used as outcome data.

## References

1. Bahn RS. Graves' ophthalmopathy. *New England Journal of Medicine* 2010 **362** 726–738. (https://doi.org/10.1056/NEJMra0905750)

2. Smith TJ & Hegedüs L. Graves' Disease. *New England Journal of Medicine* 2016 **375** 1552–1565. (https://doi.org/10.1056/NEJMra1510030)

3. Chin YH, Ng CH, Lee MH, Koh JWH, Kiew J, Yang SP, Sundar G & Khoo CM. Prevalence of thyroid eye disease in Graves' disease: A meta-analysis and systematic review. *Clinical Endocrinology* 2020 **93** 363–374. (https://doi.org/10.1111/cen.14296)

4. Brand OJ, Barrett JC, Simmonds MJ, Newby PR, McCabe CJ, Bruce CK, Kysela B, Carr-Smith JD, Brix T, Hunt PJ, *et al.* Association of the thyroid stimulating hormone receptor gene (*TSHR*) with Graves' disease. *Human Molecular Genetics* 2009 **18** 1704–1713. (https://doi.org/10.1093/hmg/ddp087)

5. Chu X, Pan CM, Zhao SX, Liang J, Gao GQ, Zhang XM, Yuan GY, Li CG, Xue LQ, Shen M, *et al.* A genome-wide association study identifies two new risk loci for Graves' disease. *Nature Genetics* 2011 **43** 897–901. (https://doi.org/10.1038/ng.898)

6. Marinò M, Rotondo Dottore G, Menconi F, Comi S, Cosentino G, Rocchi R, Latrofa F, Figus M & Santini F. Role of genetics and epigenetics in Graves' orbitopathy. *European Thyroid Journal* 2024 **13** e240179. (https://doi.org/10.1530/ETJ-24-0179)

7. Tsui S, Naik V, Hoa N, Hwang CJ, Afifiyan NF, Sinha Hikim A, Gianoukakis AG, Douglas RS & Smith TJ. Evidence for an association between thyroid-stimulating hormone and insulin-like growth factor 1 receptors: a tale of two antigens implicated in Graves' disease. *Journal of Immunology* 2008 **181** 4397–4405. (https://doi.org/10.4049/jimmunol.181.6.4397)

8. Krieger CC, Boutin A, Jang D, Morgan SJ, Banga JP, Kahaly GJ, Klubo-Gwiezdzinska J, Neumann S & Gershengorn MC. Arrestin-β-1 Physically Scaffolds TSH and IGF1 Receptors to Enable Crosstalk. *Endocrinology* 2019 **160** 1468–1479. (https://doi.org/10.1210/en.2019-00055)

9. Smith TJ, Kahaly GJ, Ezra DG, Fleming JC, Dailey RA, Tang RA, Harris GJ, Antonelli A, Salvi M, Goldberg RA, *et al.* Teprotumumab for Thyroid-Associated Ophthalmopathy. *New England Journal of Medicine* 2017 **376** 1748–1761. (https://doi.org/10.1056/NEJMoa1614949)

10. Douglas RS, Kahaly GJ, Patel A, Sile S, Thompson EHZ, Perdok R, Fleming JC, Fowler BT, Marcocci C, Marinò M, *et al.* Teprotumumab for the Treatment of Active Thyroid Eye Disease. *New England Journal of Medicine* 2020 **382** 341–352. (https://doi.org/10.1056/NEJMoa1910434)

11. Nelson MR, Tipney H, Painter JL, Shen J, Nicoletti P, Shen Y, Floratos A, Sham PC, Li MJ, Wang J, *et al.* The support of human genetic evidence for approved drug indications. *Nature Genetics* 2015 **47** 856–860. (https://doi.org/10.1038/ng.3314)

12. Minikel EV, Painter JL, Dong CC & Nelson MR. Refining the impact of genetic evidence on clinical success. *Nature* 2024 **629** 624–629. (https://doi.org/10.1038/s41586-024-07316-0)

13. Sanderson E, Glymour MM, Holmes MV, Kang H, Morrison J, Munafò MR, Palmer T, Schooling CM, Wallace C, Zhao Q, *et al.* Mendelian randomization. *Nature Reviews Methods Primers* 2022 **2** 6. (https://doi.org/10.1038/s43586-021-00092-5)

14. Gill D, Georgakis MK, Walker VM, Schmidt AF, Gkatzionis A, Freitag DF, Finan C, Hingorani AD, Howson JMM, Burgess S, *et al.* Mendelian randomization for studying the effects of perturbing drug targets. *Wellcome Open Research* 2021 **6** 16. (https://doi.org/10.12688/wellcomeopenres.16544.2)

15. Giambartolomei C, Vukcevic D, Schadt EE, Franke L, Hingorani AD, Wallace C & Plagnol V. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genetics* 2014 **10** e1004383. (https://doi.org/10.1371/journal.pgen.1004383)

16. Zuber V, Grinberg NF, Gill D, Manipur I, Slob EAW, Patel A, Wallace C & Burgess S. Combining evidence from Mendelian randomization and colocalization: Review and comparison of approaches. *American Journal of Human Genetics* 2022 **109** 767–782. (https://doi.org/10.1016/j.ajhg.2022.04.001)

17. Zheng J, Haberland V, Baird D, Walker V, Haycock PC, Hurle MR, Gutteridge A, Erola P, Liu Y, Luo S, *et al.* Phenome-wide Mendelian randomization mapping the influence of the plasma proteome on complex diseases. *Nature Genetics* 2020 **52** 1122–1131. (https://doi.org/10.1038/s41588-020-0682-6)

18. Ji Q, Xu H, Chen H, Chen X, Wang S & Zou J. Pathogenic genes associated with immune-related genes in graves' disease: a multi-omics Mendelian randomization analysis. *Scientific Reports* 2025 **15** 37875. (https://doi.org/10.1038/s41598-025-21754-4)

19. Zhang H, Jiang S, Zhu T, Liu Y, Li J, Fang S, Li Y, Sun J, He X, Tong C, *et al.* Thyroid-stimulating hormone receptor mediates peripheral-central neuroimmune crosstalk in autoimmune thyroid diseases. *BMC Medicine* 2026 **24** 404. (https://doi.org/10.1186/s12916-026-04957-y)

20. Ke C, Yu Y, Li J, Yu Y, Sun Y, Wang Y, Wang B, Lu Y, Tang M, Wang N, *et al.* Genetic and Plasma Proteomic Approaches to Identify Therapeutic Targets for Graves' Disease and Graves' Ophthalmopathy. *ImmunoTargets and Therapy* 2025 **14** 87–98. (https://doi.org/10.2147/ITT.S494692)

21. Finan C, Gaulton A, Kruger FA, Lumbers RT, Shah T, Engmann J, Galver L, Kelley R, Karlsson A, Santos R, *et al.* The druggable genome and support for target identification and validation in drug development. *Science Translational Medicine* 2017 **9** eaag1166. (https://doi.org/10.1126/scitranslmed.aag1166)

22. Ueda H, Howson JM, Esposito L, Heward J, Snook H, Chamberlain G, Rainbow DB, Hunter KM, Smith AN, Di Genova G, *et al.* Association of the T-cell regulatory gene *CTLA4* with susceptibility to autoimmune disease. *Nature* 2003 **423** 506–511. (https://doi.org/10.1038/nature01621)

23. Skrivankova VW, Richmond RC, Woolf BAR, Yarmolinsky J, Davies NM, Swanson SA, VanderWeele TJ, Higgins JPT, Timpson NJ, Dimou N, *et al.* Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: The STROBE-MR Statement. *JAMA* 2021 **326** 1614–1621. (https://doi.org/10.1001/jama.2021.18236)

24. Võsa U, Claringbould A, Westra HJ, Bonder MJ, Deelen P, Zeng B, Kirsten H, Saha A, Kreuzhuber R, Yazar S, *et al.* Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. *Nature Genetics* 2021 **53** 1300–1310. (https://doi.org/10.1038/s41588-021-00913-z)

25. Sakaue S, Kanai M, Tanigawa Y, Karjalainen J, Kurki M, Koshiba S, Narita A, Konuma T, Yamamoto K, Akiyama M, *et al.* A cross-population atlas of genetic associations for 220 human phenotypes. *Nature Genetics* 2021 **53** 1415–1424. (https://doi.org/10.1038/s41588-021-00931-x)

26. Dönertaş HM, Fabian DK, Valenzuela MF, Partridge L & Thornton JM. Common genetic associations between age-related diseases. *Nature Aging* 2021 **1** 400–412. (https://doi.org/10.1038/s43587-021-00051-5)

27. Kurki MI, Karjalainen J, Palta P, Sipilä TP, Kristiansson K, Donner KM, Reeve MP, Laivuori H, Aavikko M, Kaunisto MA, *et al.* FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature* 2023 **613** 508–518. (https://doi.org/10.1038/s41586-022-05473-8)

28. Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM & Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 2015 **4** 7. (https://doi.org/10.1186/s13742-015-0047-8)

29. 1000 Genomes Project Consortium, Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, *et al.* A global reference for human genetic variation. *Nature* 2015 **526** 68–74. (https://doi.org/10.1038/nature15393)

30. Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, Montgomery GW, Goddard ME, Wray NR, Visscher PM, *et al.* Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. *Nature Genetics* 2016 **48** 481–487. (https://doi.org/10.1038/ng.3538)

31. Lloyd-Jones LR, Robinson MR, Yang J & Visscher PM. Transformation of Summary Statistics from Linear Mixed Model Association on All-or-None Traits to Odds Ratio. *Genetics* 2018 **208** 1397–1408. (https://doi.org/10.1534/genetics.117.300360)

32. Hemani G, Zheng J, Elsworth B, Wade KH, Haberland V, Baird D, Laurin C, Burgess S, Bowden J, Langdon R, *et al.* The MR-Base platform supports systematic causal inference across the human phenome. *eLife* 2018 **7** e34408. (https://doi.org/10.7554/eLife.34408)

33. Bowden J, Davey Smith G & Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. *International Journal of Epidemiology* 2015 **44** 512–525. (https://doi.org/10.1093/ije/dyv080)

34. Bowden J, Davey Smith G, Haycock PC & Burgess S. Consistent Estimation in Mendelian Randomization with Some Invalid Instruments Using a Weighted Median Estimator. *Genetic Epidemiology* 2016 **40** 304–314. (https://doi.org/10.1002/gepi.21965)

35. Hartwig FP, Davey Smith G & Bowden J. Robust inference in summary data Mendelian randomization via the zero modal pleiotropy assumption. *International Journal of Epidemiology* 2017 **46** 1985–1998. (https://doi.org/10.1093/ije/dyx102)

36. Colobran R, Armengol Mdel P, Faner R, Gärtner M, Tykocinski LO, Lucas A, Ruiz M, Juan M, Kyewski B & Pujol-Borrell R. Association of an SNP with intrathymic transcription of *TSHR* and Graves' disease: a role for defective thymic tolerance. *Human Molecular Genetics* 2011 **20** 3415–3423. (https://doi.org/10.1093/hmg/ddr247)

37. Marín-Sánchez A, Álvarez-Sierra D, González O, Lucas-Martin A, Sellés-Sánchez A, Rudilla F, Enrich E, Colobran R & Pujol-Borrell R. Regulation of *TSHR* Expression in the Thyroid and Thymus May Contribute to TSHR Tolerance Failure in Graves' Disease Patients via Two Distinct Mechanisms. *Frontiers in Immunology* 2019 **10** 1695. (https://doi.org/10.3389/fimmu.2019.01695)

38. Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. *PLoS Genetics* 2020 **16** e1008720. (https://doi.org/10.1371/journal.pgen.1008720)

39. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* 2020 **369** 1318–1330. (https://doi.org/10.1126/science.aaz1776)

40. Chen CR, Pichurin P, Nagayama Y, Latrofa F, Rapoport B & McLachlan SM. The thyrotropin receptor autoantigen in Graves disease is the culprit as well as the victim. *Journal of Clinical Investigation* 2003 **111** 1897–1904. (https://doi.org/10.1172/JCI17069)

41. Berchner-Pfannschmidt U, Moshkelgosha S, Diaz-Cano S, Edelmann B, Görtz GE, Horstmann M, Noble A, Hansen W, Eckstein A & Banga JP. Comparative Assessment of Female Mouse Model of Graves' Orbitopathy Under Different Environments, Accompanied by Proinflammatory Cytokine and T-Cell Responses to Thyrotropin Hormone Receptor Antigen. *Endocrinology* 2016 **157** 1673–1682. (https://doi.org/10.1210/en.2015-1829)

42. Daser A, Gulbins A, Horstmann M, Kalisch J, Nolte I, Hashemi Arani F, Lang S, Mattheis S, Bechrakis NE, Banga JP, *et al.* TNFα drives disease-specific chemokine networks and orbital tissue remodeling in a preclinical model of thyroid eye disease. *European Thyroid Journal* 2026 **15** ETJ250407. (https://doi.org/10.1530/ETJ-25-0407)

43. Eckstein AK, Lösch C, Glowacka D, Schott M, Mann K, Esser J & Morgenthaler NG. Euthyroid and primarily hypothyroid patients develop milder and significantly more asymmetrical Graves ophthalmopathy. *British Journal of Ophthalmology* 2009 **93** 1052–1056. (https://doi.org/10.1136/bjo.2007.137265)

44. Douglas RS, Couch S, Wester ST, Fowler BT, Liu CY, Subramanian PS, Tang R, Nguyen QT, Maamari RN, Ugradar S, *et al.* Efficacy and Safety of Teprotumumab in Patients With Thyroid Eye Disease of Long Duration and Low Disease Activity. *Journal of Clinical Endocrinology and Metabolism* 2024 **109** 25–35. (https://doi.org/10.1210/clinem/dgad637)

45. Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal variants. *PLoS Genetics* 2021 **17** e1009440. (https://doi.org/10.1371/journal.pgen.1009440)
