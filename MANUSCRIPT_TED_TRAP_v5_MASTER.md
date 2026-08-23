**Running head:** Diverging genetic architecture of TSHR and IGF1R in GD/TED

# Genetic susceptibility and therapeutic target biology diverge at TSHR and IGF1R in Graves disease and thyroid eye disease

Jungyul Park¹, Min-Seon Kim², Kyung-Hwa Shin³⁎, Suk-Woo Yang¹⁎

¹ Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
² Department of Ophthalmology, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
³ Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea  

\* These authors contributed equally as corresponding authors.

**Corresponding author:** Suk-Woo Yang, MD, PhD, Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea. E-mail: yswoph@catholic.ac.kr; Tel: +82-2-2258-2847.  
**Co-corresponding author:** Kyung-Hwa Shin, MD, Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea.

---

## Abstract

**Purpose:** Teprotumumab targets IGF-1R in thyroid eye disease (TED), but whether *IGF1R* represents a genetically anchored susceptibility locus or a pharmacologic effector axis distinct from the established autoantigen *TSHR* remains unclear. We applied a druggable-gene-wide Mendelian randomization (MR), colocalization, and orbital tissue framework to separate susceptibility from effector axes in Graves disease (GD) and TED. **Methods:** Blood *cis*-eQTL instruments from eQTLGen were used for two-sample MR across 2,545 MR-testable druggable genes. Outcomes were Biobank Japan (BBJ) Graves disease (discovery), UK Biobank (UKB) hyperthyroidism (replication), and FinnGen R12 Graves ophthalmopathy (TED-enriched sensitivity). Bonferroni-significant discovery hits underwent cross-outcome evaluation, fine-mapping, Bayesian colocalization, and candidate filtering. *TSHR*, *IGF1R*, and *CTLA4* were further integrated with in-house orbital RNA-seq. **Results:** Thirteen genes reached Bonferroni significance in discovery. *TSHR* showed consistent protective MR effects across all three outcomes, strong colocalization with the discovery and TED-enriched signals (PP.H4=0.951 and 0.986; shared variant rs179252), and exploratory upregulation in TED orbital tissue (log2FC=+2.33, adjusted P=0.032). *IGF1R* showed directionally consistent risk-direction MR estimates, reaching nominal significance in the BBJ and UKB outcomes but not in the smaller FinnGen outcome, and lacked expression colocalization or significant orbital tissue association. *CTLA4* served as a positive-control autoimmune locus. No non-known candidate survived robust cross-outcome colocalization filtering. **Conclusion:** This framework supports *TSHR* as an expression-colocalized susceptibility anchor for GD, whereas the *IGF1R* profile is more compatible with pharmacologic effector biology than with the *TSHR*-like expression-colocalized susceptibility architecture, while not excluding a more modest genetic contribution.

**Keywords:** Graves disease; Thyroid eye disease; *TSHR*; *IGF1R*; Mendelian randomization; Colocalization

---

## Introduction

Thyroid eye disease (TED) is the most common extrathyroidal manifestation of Graves disease (GD), an autoimmune disorder in which orbital inflammation and tissue remodeling can impair quality of life, ocular function, and vision [1, 2]. For decades, treatment relied largely on corticosteroids, orbital radiotherapy, and surgery, with limited disease-modifying options. The introduction of teprotumumab, a monoclonal antibody targeting the insulin-like growth factor 1 receptor (IGF-1R), changed this therapeutic landscape: randomized trials demonstrated improvement in proptosis, clinical activity, and diplopia in active TED [3], establishing IGF-1R blockade as a viable therapeutic strategy. However, clinical response to IGF-1R blockade does not by itself resolve whether *IGF1R* represents a genetically anchored susceptibility locus for GD/TED or a pharmacologic effector target that can be modulated without being a primary driver of disease risk.

*TSHR* and IGF-1R occupy different but biologically connected positions in GD/TED pathogenesis and therapy. *TSHR* is the central autoantigen of GD: loss of immune tolerance to the thyrotropin receptor gives rise to TSHR-directed stimulating autoantibodies, and *TSHR* is among the most consistently replicated genetic susceptibility loci for the disease [2]. IGF-1R, by contrast, has emerged as both a mechanistic candidate and a therapeutic target in TED [4]. Experimental studies have reported that IGF-1R can associate with TSHR, including through β-arrestin 1–dependent receptor scaffolding, and that the two receptors engage in functional crosstalk in orbital fibroblasts [5, 6], providing a biological rationale for IGF-1R blockade. These observations leave open a fundamental question of interpretation: whether *IGF1R* should be regarded as an independent genetic susceptibility locus, a signaling partner of TSHR, or a pharmacologic effector node that is therapeutically tractable without being a primary determinant of inherited disease risk. This distinction is difficult to resolve from tissue expression or therapeutic-response data alone, because both can reflect disease-state activation, cellular remodeling, or downstream pathway modulation rather than germline susceptibility architecture.

Human genetics offers a way to address this question that is less susceptible to confounding and reverse causation than observational expression studies. Mendelian randomization (MR) uses germline genetic variants as instruments to infer whether genetically proxied molecular exposures—here, *cis*-regulated gene expression—are associated with disease liability in a manner consistent with causality [7], while colocalization tests whether the genetic signals for expression and disease share an underlying causal variant rather than reflecting distinct linked variants [8]. Applied across a predefined druggable-gene set, this combination can prioritize loci systematically without restricting inference to a small number of preselected candidates. Prior MR and multi-omics studies have begun to nominate molecular contributors to GD, including immune-related and other molecular signals [9, 10], but they have not directly resolved whether *TSHR* and *IGF1R* occupy different positions as genetically anchored susceptibility versus pharmacologic effector axes. Integrating MR and colocalization with statistical fine-mapping and disease-tissue expression can further distinguish expression-colocalized susceptibility signals from non-colocalized or downstream effector biology.

Despite these advances, existing work has not systematically distinguished genetically anchored susceptibility loci from pharmacologic effector targets across GD- and TED-relevant outcomes. Candidate-focused analyses have examined a limited number of preselected loci, leaving the relative genetic positions of *TSHR* and *IGF1R* unresolved within a common framework; single-outcome designs can elevate locus-specific or linkage-driven signals that do not reproduce across related phenotypes, with TED outcomes often underpowered or treated only as secondary; and genetic association has seldom been integrated with colocalization, fine-mapping, and disease-relevant tissue expression in a way that separates expression-colocalized susceptibility from non-colocalized effector biology.

We therefore applied a druggable-gene-wide *cis*-eQTL MR and colocalization framework across a prespecified GD/TED outcome hierarchy: Biobank Japan Graves disease as the primary discovery outcome, UK Biobank hyperthyroidism as replication, and FinnGen Graves ophthalmopathy as a TED-enriched sensitivity outcome. We integrated this framework with *TSHR* fine-mapping, exploratory in-house orbital RNA-seq, and *CTLA4* as a positive-control autoimmune locus. Our aim was to distinguish expression-colocalized susceptibility anchors from pharmacologic effector axes across druggable genes relevant to GD/TED biology and therapeutic interpretation, rather than treating discovery MR hits alone as sufficient evidence for target nomination.

---

## Methods

### Study design and data sources

We performed a druggable-gene-wide Mendelian randomization (MR) study integrating MR, colocalization, fine-mapping, and orbital tissue evidence to distinguish genetically anchored susceptibility from pharmacologic effector biology in Graves disease (GD) and thyroid eye disease (TED) (Figure 1), following the STROBE-MR guideline [11] (Table S7). Exposures were blood *cis*-eQTL instruments from the eQTLGen Consortium (predominantly European; up to 31,684 participants) [12].

Disease outcomes followed a prespecified hierarchy of decreasing sample size but increasing TED specificity: Graves disease from Biobank Japan (BBJ; East Asian; 2,809 cases) as the primary discovery outcome, hyperthyroidism from UK Biobank (UKB; European; 3,731 cases) as cross-ancestry replication [13], and Graves ophthalmopathy from FinnGen R12 (European; 858 cases) as the TED-enriched sensitivity outcome [14] (Table 1). Because TED-specific GWAS remain substantially smaller than GD GWAS, this design tested the signal first in a better-powered discovery set and then in a TED-enriched outcome. An in-house orbital RNA-seq dataset (Korean; four TED, one control) provided exploratory tissue-level triangulation and was not used for instrument selection. Three genes were designated a priori as the interpretive backbone on biological and therapeutic grounds rather than on screen results: *TSHR* (the Graves disease autoantigen), *IGF1R* (the target of teprotumumab), and *CTLA4* (a positive-control autoimmune locus). All three were carried through every evidence layer irrespective of their discovery *P* values.

### Druggable gene set and instrument selection

Analyses were restricted to the druggable genome (Finan et al. [15]; 4,462 genes, GRCh37). For each gene, *cis*-eQTL variants within ±1 Mb reaching genome-wide significance (*P* < 5×10⁻⁸) were eligible as instruments; a relaxed threshold (*P* < 5×10⁻⁶) was used only for secondary or sensitivity analyses and did not define the primary denominator. Because eQTLGen reports Z-scores without allele frequencies, 1000 Genomes Phase 3 European frequencies (n = 503) were used to reconstruct exposure effect sizes [16]. Independent instruments were obtained by LD clumping (*r*² < 0.001, 10-Mb window) using local PLINK v1.9 [17] with the same reference panel. Instrument strength was quantified as *F* = *Z*² (minimum *F* = 14.2; no weak instruments).

After clumping, 2,545 of 4,462 genes retained at least one valid instrument and formed the MR-testable denominator (6,136 SNP-level instruments; Table S4), which defined the study-wide significance threshold (Bonferroni *P* < 1.965×10⁻⁵; 0.05/2,545).

### Outcome data and harmonization

Exposure and outcome data were harmonized with TwoSampleMR [18] (harmonise_data, action = 2), aligning alleles to a common effect allele and inferring strand for complementary alleles where possible; palindromic variants with intermediate frequency (MAF > 0.42) were removed to avoid strand ambiguity.

Outcome estimates were placed on a common log-odds scale: BBJ and FinnGen were used directly, while UKB hyperthyroidism (reported on a linear mixed-model scale) was rescaled to log-odds using case prevalence before MR. Given differences in ancestry, phenotype, and effect-scale reconstruction across cohorts, cross-outcome interpretation emphasized direction, replication, and colocalization rather than exact equality of effect magnitudes.

### Mendelian randomization

Two-sample MR was performed for each gene against all three outcomes. For single-instrument loci, the Wald ratio was used; for multi-instrument loci, the inverse-variance-weighted (IVW) method was the primary estimator, with MR-Egger [19], weighted median [20], and weighted mode as sensitivity estimators. Directional pleiotropy was assessed using the MR-Egger intercept, between-instrument heterogeneity using Cochran's Q, and exposure–outcome direction using the Steiger test [21] where estimable (Table S6). Discovery significance was defined at the study-wide Bonferroni threshold (*P* < 1.965×10⁻⁵; 0.05/2,545). Replication and TED-enriched support were evaluated using direction-consistent nominal association (*P* < 0.05).

Because a screen that advances no novel target is interpretable only in relation to the effect sizes it could have detected, we quantified the detectable effect range of the screen itself. For each gene and outcome, the minimum effect detectable with 80% power at the applicable significance threshold was computed from the empirical standard error of the primary MR estimate as |β|min = (*z*1−α/2 + *z*0.80) × SE, using α = 1.965×10⁻⁵ for the discovery outcome and α = 0.05 for the replication and TED-enriched outcomes. This design-based approach uses the standard errors actually obtained after harmonization rather than reconstructed variance components, and is reported for all genes and for the three backbone loci (Table S8).

### Fine-mapping and colocalization

To evaluate whether the single *TSHR* instrument reflected an independent locus-level signal rather than unresolved allelic heterogeneity, the *TSHR* locus was fine-mapped using SuSiE [22] (susieR; susie_rss) with summary statistics and 1000 Genomes Phase 3 European and East Asian LD reference panels. The primary signal was anchored on rs179252; secondary credible-set lead variants remained in substantial LD with rs179252 and did not provide clearly independent additional instruments. *TSHR* was therefore retained as a single-instrument, locus-level anchor rather than reconstructed as a multi-instrument exposure (Table S3; Figure S1).

For prioritized loci, we tested whether the *cis*-eQTL and outcome-association signals shared a causal variant using coloc.abf [8] (coloc) within ±1 Mb of the gene, under the single-causal-variant assumption and default priors (p1 = p2 = 10⁻⁴, p12 = 10⁻⁵) [23]. In this implementation dataset 1 was the outcome GWAS and dataset 2 the *cis*-eQTL, so the posterior probabilities distinguish association with the outcome only (H1), with the *cis*-eQTL only (H2), with both traits at distinct causal variants (H3), and with both traits at a shared causal variant (H4). Strong colocalization was defined as PP.H4 ≥ 0.80 and evaluated separately against the discovery and TED-enriched outcomes (Table S2; Figure S2).

### Orbital tissue transcriptomic analysis

Exploratory tissue-level support was assessed using an in-house bulk RNA-seq dataset of orbital tissue from a Korean cohort comprising four TED patients and one control. Technical replicates were collapsed to the biological-sample level by summing raw counts, and gene-level differential expression was estimated with DESeq2 [24] with Benjamini–Hochberg correction computed across genes. This study reports only the three backbone genes (*TSHR*, *IGF1R*, *CTLA4*) for hypothesis-driven tissue triangulation (Figure 3C; Table 2); no genome-wide transcriptomic catalogue, differential-expression list, or pathway-enrichment result is presented here. The dataset was used solely for exploratory tissue-level support of the genetically prioritized backbone genes, not for genetic instrument selection, and no external transcriptomic dataset was incorporated.

### Software, reporting, and data availability

Analyses used R 4.3.3 with TwoSampleMR (v0.7.4), ieugwasr, coloc (v5.2.3), susieR (v0.14.2), and DESeq2 (v1.42.1); LD operations used PLINK v1.9. Reporting follows the STROBE-MR guideline (Table S7). Public eQTL and GWAS summary statistics are available from their original repositories, including eQTLGen, the GWAS Catalog, and FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions. The study used de-identified summary statistics and an institutionally approved orbital tissue dataset (IRB approval, Pusan National University Hospital, 2104-018-102).

---

## Results

### Druggable-gene-wide MR discovery

The attrition of the screen and the resulting contrast between the two backbone receptors are summarized in **Figure 1**. The analysis followed the prespecified outcome hierarchy of BBJ Graves disease (discovery), UKB hyperthyroidism (replication), and FinnGen R12 Graves ophthalmopathy (TED-enriched sensitivity), with data sources in **Table 1**.

Of 4,462 druggable genes, 2,545 had at least one valid *cis*-eQTL instrument in eQTLGen and formed the MR-testable denominator; across these genes, 6,136 SNP-level instruments were retained, all exceeding conventional instrument-strength criteria (minimum F-statistic = 14.2, no weak instruments; **Table S4**). Primary MR estimates against the BBJ Graves disease discovery outcome were obtainable for 2,235 genes.

Applying a Bonferroni threshold corrected for the MR-testable denominator (*P* < 1.965×10⁻⁵; 0.05/2,545), 13 genes reached significance in discovery (**Figure 2**). These hits included the established Graves disease autoantigen *TSHR* and the known autoimmune locus *CTLA4*, both with strongly significant protective-direction estimates, together with additional loci that were carried forward for replication and colocalization-based filtering. The discovery volcano (**Figure 2**) shows the distribution of BBJ MR effect estimates and the 13 Bonferroni-significant hits; full per-hit statistics, cross-outcome colocalization, and classification are provided in **Table 3** (with extended detail in **Table S1**).

### Cross-outcome classification and colocalization filtering

The 13 discovery hits were not treated as putative targets on the strength of the discovery MR alone. Instead, each was classified according to prior biological status, genomic context, cross-outcome support, and colocalization with both the discovery and TED-enriched signals. Two hits—*TSHR* and *CTLA4*—were recognized as established Graves disease and autoimmune loci, respectively, and are examined in detail below. The remaining hits were evaluated as non-known candidates, including a cluster on chromosome 16p11.2 (*HSD3B7*, *VKORC1*, and *PRSS36*) whose top eQTL instruments were co-located within a narrow approximately 143-kb linkage-disequilibrium window (Table S5).

Candidate colocalization showed outcome-dependent patterns (Table S2; Figure S2). *TNFSF14* and *IFNGR1* showed strong colocalization with the BBJ discovery signal (PP.H4 = 0.994 and 0.989, respectively), but did not colocalize with the TED-enriched outcome (PP.H4 = 0.019 for both), where the posterior instead favored an eQTL-only pattern (PP.H2 = 0.63 and 0.64): the *cis*-eQTL signal was present, but no association with the TED outcome was detectable at these loci. *MAPKAPK5* showed a distinct-signal pattern in discovery (PP.H3 ≈ 1.0), whereas the chr16p11.2 candidates showed only weak or ambiguous colocalization in discovery (PP.H4 ≤ 0.62) and no shared colocalization support in the TED-enriched outcome, where they too showed an eQTL-only pattern (PP.H2 = 0.73–0.74).

Overall, no non-known candidate showed strong colocalization in both BBJ Graves disease and the TED-enriched outcome. Thus, after cross-outcome colocalization filtering, the screen yielded no robust novel target. This filtering result was retained as part of the final interpretation because it prevented overinterpretation of single-outcome or LD-driven discovery signals and strengthened the focus on the three backbone genes *TSHR*, *IGF1R*, and *CTLA4*.

The detectable range of the screen defines what this null does and does not exclude (**Table S8**). At the study-wide discovery threshold and 80% power, the median gene was powered to detect an odds ratio of approximately 2.55 per unit genetically proxied expression (interquartile range 1.71–5.19); 35.6% of genes were powered for OR ≥ 2.0 and 14.6% for OR ≥ 1.5. In the TED-enriched outcome at nominal significance, the corresponding median detectable OR was approximately 2.12. The screen was therefore well powered for large expression-mediated effects — consistent with its recovery of *TSHR* (β = −2.10) and *CTLA4* (β = −1.74), both of which exceeded their gene-specific detection thresholds — but not for modest effects. The absence of a robust novel target should accordingly be read as evidence against additional large, expression-colocalized druggable effects rather than against moderate ones, which larger GD- and TED-specific GWAS will be required to resolve.

### TSHR: an expression-colocalized susceptibility anchor

Among the established hits, *TSHR* showed the strongest convergence across the genetic, colocalization, fine-mapping, and tissue layers of the framework (Table 2; Figure 3). In MR, the *TSHR* *cis*-eQTL instrument was associated with protective-direction estimates across all three outcomes: BBJ Graves disease (β = −2.10, *P* = 1.1×10⁻¹⁴), UKB hyperthyroidism (β = −2.44, *P* = 8.8×10⁻²⁸), and the TED-enriched outcome (β = −2.33, *P* = 2.8×10⁻⁷). The consistency of direction and approximate magnitude across the independently ascertained Graves disease and hyperthyroidism outcomes argues against an outcome-specific artifact (Figure 3B). The concordant estimate in the TED-enriched outcome is interpreted more conservatively, because Graves ophthalmopathy cases are ascertained among patients with Graves disease and compared with population controls, so that outcome substantially re-measures shared Graves disease susceptibility rather than liability specific to orbital involvement (see Limitations).

Colocalization supported a shared causal variant between *TSHR* expression and disease association in both the discovery and TED-enriched outcomes (PP.H4 = 0.951 and 0.986, respectively; Figure 3A), with rs179252 prioritized as the shared lead variant (Table 2). Fine-mapping and ancestry-matched LD assessment supported retention of rs179252 as the single-IV/locus-level anchor: secondary credible-set lead variants remained in substantial LD with rs179252 and did not provide clearly independent additional instruments (Table S3; Figure S1).

At the tissue level, *TSHR* showed exploratory upregulation in TED orbital tissue relative to control tissue (log2 fold-change = +2.33, adjusted *P* = 0.032; Figure 3C); because the dataset contained a single control, this provides supportive but non-confirmatory tissue-level context rather than a replicated differential-expression result. Together, these results position *TSHR* as an expression-colocalized upstream susceptibility anchor for Graves disease — resolved to a shared causal variant and a single credible set rather than merely rediscovered as an association — whose signal is also recovered in a TED-enriched Graves ophthalmopathy outcome, although a TED-specific effect separable from Graves disease liability is not established by these data. The apparent contrast between the protective-direction germline MR estimate and disease-state orbital upregulation is addressed in the Discussion.

### IGF1R: nominal genetic association without expression colocalization

In contrast to *TSHR*, *IGF1R*, the therapeutic target of teprotumumab, showed a markedly different evidence profile (Table 2; Figure 3). *IGF1R* was instrumented by multiple independent *cis*-eQTL variants, with four instruments in the discovery and replication analyses. The inverse-variance-weighted MR estimate indicated a nominal positive-direction association with BBJ Graves disease (β = +0.45, *P* = 0.021) and UKB hyperthyroidism (β = +0.30, *P* = 0.012); the TED-enriched estimate was directionally consistent but non-significant (β = +0.34, *P* = 0.18). Sensitivity estimators were broadly concordant in direction without evidence of directional pleiotropy, although statistical significance was not uniform across estimators, consistent with a modest effect (Table S6). The magnitude of this effect lay below the gene-specific detection threshold of the discovery screen (minimum effect detectable with 80% power at the study-wide threshold, |β| = 0.99 for *IGF1R*; Table S8), so its failure to reach study-wide significance is expected for an effect of this size and does not by itself constitute evidence against a genetic contribution. This is one reason *IGF1R* was carried forward as a backbone gene evaluated across multiple evidence layers rather than judged on discovery significance alone.

Critically, this nominal MR signal was not accompanied by expression colocalization. At the *IGF1R* locus, the posterior favored an eQTL-only pattern rather than a shared causal variant in both the discovery and TED-enriched outcomes (PP.H2 = 0.79 and 0.68; PP.H4 = 0.043 and 0.036; Figure 3A), indicating a well-defined *cis*-eQTL at the locus without a detectable disease association, in contrast to the strong shared-variant evidence observed at *TSHR*. At the tissue level, *IGF1R* was not significantly differentially expressed in TED orbital tissue (log2 fold-change = +0.41, adjusted *P* > 0.99; Figure 3C).

Thus, *IGF1R* showed nominal positive-direction MR evidence that was not supported by eQTL–GWAS colocalization or by significant orbital tissue upregulation. This combination distinguishes the genetic profile of *IGF1R* from that of an expression-colocalized susceptibility locus such as *TSHR*. The interpretation of this profile—in particular its relationship to the established therapeutic role of IGF-1R blockade—is considered in the Discussion.

### CTLA4: a positive-control autoimmune locus

*CTLA4* (chromosome 2q33), a well-established autoimmune susceptibility locus [25], was included as a positive control to calibrate the framework against a gene of known relevance to Graves disease. Consistent with its role as a known autoimmune locus, *CTLA4* showed strong protective-direction MR estimates in the discovery analysis (BBJ Graves disease, β = −1.74, *P* = 5.5×10⁻¹⁵), with concordant protective-direction estimates in replication (UKB hyperthyroidism, β = −1.57, *P* = 0.002) and the TED-enriched outcome (β = −1.77, *P* = 0.012). The recovery of *CTLA4* as a strong autoimmune signal supports the sensitivity of the discovery screen.

Colocalization at *CTLA4* was outcome-dependent: the posterior favored distinct causal variants in the discovery outcome (PP.H3 = 0.79), but strong colocalization in the TED-enriched outcome (PP.H4 = 0.978; Figure 3A). This mixed pattern, together with heterogeneity across instruments in the replication analysis (Table S6), is consistent with locus and cross-ancestry complexity at *CTLA4* rather than a uniformly colocalized signal. We therefore interpret *CTLA4* as a positive-control autoimmune locus rather than a uniformly robust expression-colocalized anchor. At the tissue level, *CTLA4* was not significantly differentially expressed in TED orbital tissue (log2 fold-change = +1.27, adjusted *P* = 0.815; Figure 3C), consistent with an immune-system signal that was not accompanied by a significant local orbital expression change in this exploratory tissue dataset.

---

## Discussion

Using a druggable-gene-wide Mendelian randomization, colocalization, fine-mapping, and orbital tissue framework, we set out to distinguish genetically anchored susceptibility from pharmacologic effector biology in Graves disease (GD) and thyroid eye disease (TED). Two contrasting evidence profiles emerged. *TSHR*, the canonical thyroid autoantigen, showed convergent support across MR, expression colocalization, fine-mapping, and orbital tissue upregulation, supporting its role as an expression-colocalized susceptibility anchor. *IGF1R*, the target of teprotumumab, showed directionally consistent MR evidence that reached nominal significance in discovery and replication but lacked expression colocalization and significant orbital tissue upregulation — a profile more compatible with pharmacologic effector biology than with expression-colocalized susceptibility. Across 2,545 MR-testable druggable genes, no non-known candidate advanced as a robust novel target after cross-outcome colocalization, so the principal contribution of this work is not target discovery but the separation of susceptibility anchoring from effector biology.

*TSHR* is the central autoantigen of GD and a well-established genetic susceptibility locus [2]; a genetic association alone would therefore add little to the literature. The contribution of the present framework is not the rediscovery of *TSHR* but its resolution into a coherent expression-colocalized anchor: consistent protective-direction estimates across all three outcomes, a shared causal variant with disease association in both the discovery and TED-enriched outcomes at rs179252, and fine-mapping showing that secondary credible-set lead variants remained in substantial linkage disequilibrium with it and did not provide clearly independent additional instruments. The MR estimate therefore reflects a well-defined locus-level signal rather than an artifact of unmodeled allelic heterogeneity or arbitrary instrument selection. Its convergence with exploratory *TSHR* upregulation in TED orbital tissue connects the genetically prioritized locus to disease-state orbital expression and provides a reference profile against which other candidates, including *IGF1R*, can be compared, in keeping with a triangulation approach that integrates evidence from methods with differing assumptions [7].

A notable feature of the *TSHR* signal is an apparent directional mismatch: the protective-direction blood *cis*-eQTL MR estimate indicates that genetically proxied higher blood *TSHR* expression is associated with lower disease odds, whereas *TSHR* is upregulated in TED orbital tissue. These observations operate on different biological layers and are not necessarily contradictory. The MR estimate reflects inherited *cis*-regulatory variation captured in blood and indexing lifelong genetic liability, whereas the orbital finding reflects local disease-state expression in already-affected tissue, where *TSHR* upregulation may reflect orbital remodeling, inflammatory activation, or cellular-state shifts within the TED microenvironment. Because genetic regulatory effects on gene expression are frequently tissue-specific and shaped by local cell-type composition [26], a protective-direction blood *cis*-eQTL signal [12] and a disease-state increase in orbital expression are not expected to align in sign. Accordingly, we interpret the convergence across modalities as locus-level prioritization of *TSHR*, not as evidence that the direction of tissue expression change is itself causal for orbital disease, and we caution against reading tissue differential expression as directionally equivalent to a germline MR effect.

*IGF1R* presented a distinct evidence profile from *TSHR*. Although it showed nominal positive-direction MR evidence in the discovery and replication outcomes, this signal was not accompanied by expression colocalization—the posterior favored an eQTL-only pattern rather than a shared causal variant—and *IGF1R* was not significantly differentially expressed in orbital tissue. Importantly, the absence of expression colocalization is not evidence against the biological or therapeutic relevance of IGF-1R: colocalization tests a specific hypothesis—that a single shared causal variant links *cis*-regulated expression to disease liability—and a gene can be a mechanistically and therapeutically important effector without satisfying it, particularly where its disease relevance is mediated through receptor signaling or protein–protein interaction rather than through germline expression-level regulation. We therefore read the *IGF1R* profile as more compatible with pharmacologic effector biology than with a primary expression-colocalized susceptibility anchor, an interpretation these data support indirectly rather than establish. This interpretation is compatible with the clinical efficacy of teprotumumab [3], which targets IGF-1R signaling: an effective therapeutic target need not be a primary genetic determinant of disease, nor must its germline expression regulation colocalize with disease risk. The therapeutic rationale for IGF-1R blockade can be understood in relation to physical and functional crosstalk with TSHR in orbital fibroblasts [4, 5], without requiring *IGF1R* to be an independent germline susceptibility locus. Framing *IGF1R* as an effector rather than an anchor reconciles its therapeutic relevance with the comparatively weak germline expression-colocalization evidence observed in this study.

The contrast between *TSHR* and *IGF1R* illustrates a broader principle: genetic susceptibility loci and pharmacologic targets are overlapping but not identical categories, consistent with the observation that human genetic support for a target increases the likelihood of success in drug development even though not every effective target is a primary susceptibility locus [27]. That principle is already established; what this study contributes is a worked demonstration of it within a single disease area, in which the anchor and the candidate effector are resolved side by side under the same outcome hierarchy, analytic framework, and filters. Our findings do not test whether targeting an upstream anchor would offer clinical advantages over established effector-directed therapy — a question outside the scope of a genetic prioritization study.

A systematic screen of 2,545 MR-testable druggable genes did not advance any non-known candidate as a robust novel target, and we regard this as an informative result rather than a negative one. Several candidates were significant in the discovery outcome alone; for example, *TNFSF14* and *IFNGR1* colocalized strongly with the BBJ Graves disease signal but failed to show shared-variant colocalization in the TED-enriched outcome, while the chromosome 16p11.2 candidates showed only weak or ambiguous colocalization within a narrow linkage-disequilibrium window. Had we relied on discovery significance alone, any of these signals could have been overinterpreted as a putative druggable target. Requiring cross-outcome support and colocalization-based filtering guarded against the overinterpretation of single-outcome or LD-driven signals that can arise in large druggable-gene screens. This conservative outcome is compatible with the currently available genetic evidence, in which strong autoimmune and thyroid-specific loci are more consistently recovered than novel druggable candidates [9, 10]; it does not exclude additional contributors of smaller effect that larger TED-specific GWAS may eventually reveal.

Several methodological features support the robustness of these conclusions: the outcome hierarchy was prespecified, so interpretation rested on cross-outcome consistency rather than any single dataset; prioritization required both direction-consistent cross-outcome association and Bayesian colocalization [8] before a non-known signal was advanced as robust; the positive-control locus *CTLA4* [25] was recovered as a strong direction-consistent signal, indicating sensitivity to known autoimmune biology; and the *TSHR* instrument was fine-mapped rather than accepted at face value, supporting its interpretation as a locus-level anchor rather than an artifact of allelic heterogeneity. Together these reduce the likelihood that the central susceptibility–effector distinction is driven by single-outcome artifacts or linkage-driven signals.

This study has several limitations. First, the instruments were blood *cis*-eQTLs [12], which capture blood *cis*-regulatory variation rather than local expression in orbital or thyroid tissue, and genetic regulatory effects differ substantially across tissues [26]. Second, the eQTL data were predominantly European whereas the discovery outcome was East Asian [13]; concordant effects provide a form of cross-ancestry assessment, but allele-frequency and linkage-disequilibrium differences may attenuate or bias instrument performance, and exposure allele frequencies were approximated from a European reference panel. Third, the UKB replication phenotype was broader than Graves disease and may include heterogeneous causes of hyperthyroidism. Fourth, the *TSHR* signal rested on a single instrument; fine-mapping [22] supported this as a locus-level anchor, but single-instrument analyses preclude formal pleiotropy diagnostics such as the MR-Egger intercept. Fifth, the TED-enriched outcome was the smallest (858 cases), limiting power for colocalization and restricting the detectable effect range (Table S8); the eQTL-only posterior pattern seen for candidate loci in that outcome is therefore consistent with insufficient power to detect an outcome association and should not be read as positive evidence for distinct causal architecture. Sixth, and more fundamentally, that outcome ascertains Graves ophthalmopathy cases among individuals with Graves disease and compares them with population controls, so its associations substantially reflect shared Graves disease susceptibility rather than liability specific to orbital involvement; consistency of the *TSHR* estimate there should be read as reproducibility of the underlying Graves disease signal, not as evidence of a TED-specific causal effect. Separating the two would require a within-Graves contrast of patients with and without orbitopathy, which currently available GWAS do not provide at adequate power. Seventh, coloc.abf was applied under the single-causal-variant assumption without a multi-causal-variant extension, with one European reference allele-frequency vector used for both datasets at every locus including the East Asian discovery outcome, and no formal colocalization power analysis was performed; discovery-arm posteriors should therefore be regarded as approximate. Eighth, the framework included a single positive control, which did not behave uniformly: *CTLA4* colocalized in the TED-enriched outcome (PP.H4 = 0.978, lead variant rs1863800) but not in discovery (PP.H3 = 0.794, lead variant rs231811). Its recovery as a strong direction-consistent MR signal calibrates the MR component of the screen, but the sensitivity of the combined MR-plus-colocalization filter is not established by one control locus. Finally, the orbital tissue analysis was exploratory, rested on four TED samples and a single control, and was not externally replicated; it provides context rather than confirmation. These constraints should be considered when interpreting individual estimates, but the central contrast between *TSHR* and *IGF1R* was supported by multiple complementary layers.

In this druggable-gene-wide Mendelian randomization, colocalization, and exploratory orbital tissue framework, *TSHR* emerged as an expression-colocalized susceptibility anchor for GD — supported by consistent cross-outcome genetic effects and shared-variant colocalization at rs179252 — with its signal also recovered in a TED-enriched Graves ophthalmopathy outcome. *IGF1R* showed directionally consistent risk-direction estimates, nominally significant in BBJ and UKB but not FinnGen, without evidence of expression colocalization or a significant orbital tissue association. This profile is more compatible with pharmacologic effector biology than with the *TSHR*-like expression-colocalized susceptibility architecture, although it does not establish *IGF1R* as an effector axis or exclude a modest genetic contribution. No non-known druggable gene advanced as a robust novel target after cross-outcome colocalization. Distinguishing susceptibility anchoring from effector biology may help interpret existing therapies and guide the genetic prioritization of future targets in GD and TED.

---

## Declarations

**Funding.** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Conflict of interest.** The authors declare that they have no conflict of interest.

**Ethics approval.** This study used publicly available, de-identified summary statistics and an institutionally approved orbital tissue dataset. The in-house orbital transcriptomic component was approved by the Institutional Review Board of Pusan National University Hospital (approval number 2104-018-102) and was conducted in accordance with the Declaration of Helsinki.

**Informed consent.** Written informed consent was obtained from all individual participants included in the in-house orbital tissue study.

**Data availability.** Public summary statistics analyzed in this study are available from their original repositories: blood *cis*-eQTL data from the eQTLGen Consortium, the Biobank Japan Graves disease and UK Biobank hyperthyroidism genome-wide association statistics through the GWAS Catalog, and the FinnGen Release 12 Graves ophthalmopathy statistics from FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions.

**Author contributions.** J.P. conceived and designed the study, performed the analyses, and drafted the manuscript. M.-S.K. contributed to data collection and interpretation. K.-H.S. and S.-W.Y. supervised the study and revised the manuscript. All authors read and approved the final manuscript.

**Acknowledgements.** We thank the eQTLGen Consortium, Biobank Japan, the UK Biobank, and the FinnGen study and its participants for making their summary statistics publicly available. We acknowledge the GWAS Catalog for hosting and distributing the genome-wide association summary statistics used as outcome data.

---

## References

1. Bahn RS (2010) Graves' ophthalmopathy. N Engl J Med 362:726–738
2. Smith TJ, Hegedüs L (2016) Graves' disease. N Engl J Med 375:1552–1565
3. Douglas RS, Kahaly GJ, Patel A, et al (2020) Teprotumumab for the treatment of active thyroid eye disease. N Engl J Med 382:341–352
4. Smith TJ (2019) The insulin-like growth factor-I receptor and its role in thyroid-associated ophthalmopathy. Eye (Lond) 33:200–205
5. Krieger CC, Boutin A, Jang D, et al (2019) Arrestin-β-1 physically scaffolds TSH and IGF1 receptors to enable crosstalk. Endocrinology 160:1468–1479
6. Smith TJ (2019) Potential roles of CD34+ fibrocytes masquerading as orbital fibroblasts in thyroid-associated ophthalmopathy. J Clin Endocrinol Metab 104:581–594
7. Sanderson E, Glymour MM, Holmes MV, et al (2022) Mendelian randomization. Nat Rev Methods Primers 2:6
8. Giambartolomei C, Vukcevic D, Schadt EE, et al (2014) Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet 10:e1004383
9. Ji Q, Xu H, Chen H, et al (2025) Pathogenic genes associated with immune-related genes in Graves' disease: a multi-omics Mendelian randomization analysis. Sci Rep 15:37875
10. Li Y, Chen L, Lin S, et al (2025) Inflammatory proteins mediate the effect of gut microbiota on Graves' ophthalmopathy: a Mendelian randomization study. Transl Vis Sci Technol 14:34
11. Skrivankova VW, Richmond RC, Woolf BAR, et al (2021) Strengthening the reporting of observational studies in epidemiology using Mendelian randomization: the STROBE-MR statement. JAMA 326:1614–1621
12. Võsa U, Claringbould A, Westra HJ, et al (2021) Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and polygenic scores that regulate blood gene expression. Nat Genet 53:1300–1310
13. Sakaue S, Kanai M, Tanigawa Y, et al (2021) A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet 53:1415–1424
14. Kurki MI, Karjalainen J, Palta P, et al (2023) FinnGen provides genetic insights from a well-phenotyped isolated population. Nature 613:508–518
15. Finan C, Gaulton A, Kruger FA, et al (2017) The druggable genome and support for target identification and validation in drug development. Sci Transl Med 9:eaag1166
16. Wood AR, Perry JRB, Tanaka T, et al (2013) Imputation of variants from the 1000 Genomes Project modestly improves known associations and can identify low-frequency variant-phenotype associations undetected by direct genotyping. PLoS One 8:e64343
17. Chang CC, Chow CC, Tellier LC, et al (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience 4:7
18. Hemani G, Zheng J, Elsworth B, et al (2018) The MR-Base platform supports systematic causal inference across the human phenome. eLife 7:e34408
19. Bowden J, Davey Smith G, Burgess S (2015) Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. Int J Epidemiol 44:512–525
20. Bowden J, Davey Smith G, Haycock PC, Burgess S (2016) Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator. Genet Epidemiol 40:304–314
21. Hemani G, Tilling K, Davey Smith G (2017) Orienting the causal relationship between imprecisely measured traits using GWAS summary data. PLoS Genet 13:e1007081
22. Wang G, Sarkar A, Carbonetto P, Stephens M (2020) A simple new approach to variable selection in regression, with application to genetic fine mapping. J R Stat Soc Series B Stat Methodol 82:1273–1300
23. Wallace C (2020) Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLoS Genet 16:e1008720
24. Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15:550
25. Ueda H, Howson JMM, Esposito L, et al (2003) Association of the T-cell regulatory gene CTLA4 with susceptibility to autoimmune disease. Nature 423:506–511
26. GTEx Consortium (2020) The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369:1318–1330
27. Nelson MR, Tipney H, Painter JL, et al (2015) The support of human genetic evidence for approved drug indications. Nat Genet 47:856–860

---

## Figure Legends

**Figure 1. What survives the screen, and where *TSHR* and *IGF1R* diverge.** **(A)** Attrition through the druggable-gene-wide screen, labelled by the reason genes were removed. Of 4,462 druggable genes [15], 2,545 carried a valid eQTLGen *cis* instrument (6,136 instruments) and 2,235 yielded an estimable effect against the BBJ Graves disease discovery outcome; 13 reached the study-wide Bonferroni threshold (*P* < 1.965×10⁻⁵). Resolving those 13 (right) leaves no robust novel target: five are MHC-region signals, three fall in a single chromosome 16p11.2 linkage-disequilibrium block rather than representing three independent associations, two colocalized with the discovery outcome only (*TNFSF14* and *IFNGR1*, PP.H4 0.994 and 0.989 falling to 0.019 in the TED-enriched outcome), one reflects a distinct causal variant (*MAPKAPK5*), and two are established loci retained as the susceptibility anchor (*TSHR*) and the positive control (*CTLA4*). The screen was powered at 80% to detect a median odds ratio of 2.55 (Table S8), so this null excludes further large — but not moderate — expression-mediated effects. **(B)** The same three evidence layers were applied to *TSHR* and *IGF1R* using the same outcome hierarchy and analytic framework, with gene-specific *cis*-eQTL instruments. *TSHR* satisfies all three; *IGF1R* shows a nominal genetic association in discovery (β = +0.45, *P* = 0.021) but a *cis*-eQTL-only colocalization posterior (PP.H2 = 0.79 and 0.68) and no significant orbital tissue association — the therapeutic target does not display the susceptibility architecture seen at *TSHR*. Tissue evidence is exploratory (four TED, one control).

**Figure 2. Druggable-gene-wide BBJ discovery Mendelian randomization.** Volcano plot of *cis*-eQTL MR effects on Graves disease (BBJ discovery cohort, log-odds scale) for 2,235 druggable genes with estimable BBJ MR effects, among 2,545 MR-testable genes drawn from the druggable genome [15] with valid eQTLGen instruments. The horizontal dashed line marks the Bonferroni significance threshold (*P* = 1.965 × 10⁻⁵), computed against all 2,545 MR-testable genes. Thirteen genes exceed this threshold. *TSHR* and *CTLA4* (red), both established Graves disease loci, show strong protective-direction associations and serve as internal anchors/positive control. *HLA-A* and *HLA-DQA2* reflect MHC-region signal (grey). Three genes within the chr16p11.2 LD block (*HSD3B7*, *VKORC1*, *PRSS36*; light purple) represent a single linkage-disequilibrium signal rather than independent associations (Supplementary Table S5). *IGF1R* (blue), the therapeutic backbone gene, lies below the Bonferroni threshold with nominal association and is evaluated in the multi-layer backbone-gene analysis. No druggable gene outside the known anchors and LD/MHC artifacts achieved robust, colocalization-confirmed novel association (Supplementary Tables S1–S2).

**Figure 3. Multi-layer evidence integration for backbone genes.** Convergent genetic and tissue evidence for *TSHR*, *IGF1R*, and *CTLA4*. **(A)** Bayesian colocalization posterior probabilities (coloc.abf) for each gene against BBJ Graves disease and FinnGen Graves ophthalmopathy. Bars show PP.H2 (*cis*-eQTL association only), PP.H3 (distinct causal variants), and PP.H4 (shared causal variant); the dashed line marks PP.H4 = 0.8. *TSHR* shows strong colocalization in both outcomes; *IGF1R* is PP.H2-dominant (a *cis*-eQTL signal without a detectable outcome association); *CTLA4* shows a mixed pattern (PP.H3-dominant in BBJ, PP.H4 in FinnGen), consistent with locus complexity. **(B)** MR effect estimates (β, log-odds) with 95% confidence intervals across the three outcomes. *TSHR* and *CTLA4* show consistent protective-direction effects; *IGF1R* shows directionally consistent risk-direction estimates, with nominal significance in BBJ and UKB but not FinnGen. **(C)** Orbital tissue expression (TPM) from in-house RNA-seq; TED *n* = 4 and control *n* = 1. DESeq2 statistics (log2 fold-change, adjusted *P*) were computed at the biological-sample level after collapsing technical-replicate libraries and should be interpreted as exploratory tissue support given the single control. TPM panels use gene-specific y-axis scales. *TSHR* showed an exploratory differential-expression signal (log2FC = +2.33, adjusted *P* = 0.032), which is not confirmatory because the analysis included only one control sample; *IGF1R* and *CTLA4* were directionally up but non-significant.

---

## Tables

**Table 1. Study design and data sources.**

| Dataset | Role | Ancestry | Sample size / cases |
|---|---|---|---|
| eQTLGen blood *cis*-eQTL [12] | Exposure (instruments) | Predominantly European | up to 31,684 |
| Biobank Japan Graves disease (GCST90018627) [13] | Primary discovery outcome | East Asian | 2,809 cases |
| UK Biobank hyperthyroidism (GCST90038636) [13] | Cross-ancestry replication outcome | European | 3,731 cases |
| FinnGen R12 Graves ophthalmopathy [14] | TED-enriched sensitivity outcome | European | 858 cases |
| In-house orbital RNA-seq | Exploratory tissue support | Korean | 4 TED + 1 control |

*eQTLGen sample size indicates the maximum available blood cis-eQTL sample size; SNP-level sample sizes vary across instruments (e.g., rs179252 NrSamples = 31,566). BBJ Graves disease was used as the primary discovery phenotype, UKB hyperthyroidism as cross-ancestry/broader-phenotype replication, and FinnGen Graves ophthalmopathy as the TED-enriched sensitivity outcome. UKB binary-trait effect estimates were rescaled from a linear mixed-model scale to log-odds using case prevalence.*

**Table 2. Integrated genetic, statistical, and tissue evidence for backbone genes.**

| Gene | Outcome | N IV | OR (95% CI) | IVW/Wald *P* | Weighted median *P* | MR-Egger intercept *P* | Cochran's Q *P* | Coloc PP.H4 | Tissue log2FC (adj *P*) |
|---|---|---|---|---|---|---|---|---|---|
| ***TSHR*** | BBJ Graves | 1 | 0.12 (0.07–0.21) | 1.1×10⁻¹⁴ | NA (single IV) | NA (single IV) | NA (single IV) | 0.951 | +2.33 (0.032) |
| | UKB hyperthyroid | 1 | 0.09 (0.06–0.14) | 8.8×10⁻²⁸ | NA | NA | NA | — | |
| | FinnGen GO | 1 | 0.10 (0.04–0.24) | 2.8×10⁻⁷ | NA | NA | NA | 0.986 | |
| ***IGF1R*** | BBJ Graves | 4 | 1.56 (1.07–2.28) | 0.021 | 0.029 | 0.663 | 0.902 | 0.043 | +0.41 (>0.99, NS) |
| | UKB hyperthyroid | 4 | 1.35 (1.07–1.70) | 0.012 | 0.002 | 0.562 | 0.250 | — | |
| | FinnGen GO | 3 | 1.41 (0.85–2.33) | 0.182 | 0.124 | 0.433 | 0.309 | 0.036 | |
| ***CTLA4*** | BBJ Graves | 1 | 0.18 (0.11–0.27) | 5.5×10⁻¹⁵ | NA (single IV) | NA (single IV) | NA (single IV) | 0.206 | +1.27 (0.815, NS) |
| | UKB hyperthyroid | 2 | 0.21 (0.08–0.57) | 0.002 | — | NA (2 IV) | 0.001 | — | |
| | FinnGen GO | 2 | 0.17 (0.04–0.68) | 0.012 | — | NA (2 IV) | 0.049 | 0.978 | |

*MR effect estimates are odds ratios (OR) per unit increase in genetically proxied gene expression, by Wald ratio (single instrument) or inverse-variance weighted (IVW; multiple instruments). Sensitivity estimators (weighted median, weighted mode, MR-Egger intercept) and the Cochran's Q heterogeneity test are reported only where estimable; single-instrument loci are marked NA. The MR-Egger intercept P values for IGF1R (all > 0.43) indicate no detectable directional pleiotropy; Cochran's Q indicates significant between-instrument heterogeneity for CTLA4 in the UKB and FinnGen analyses (Q P = 0.001 and 0.049), consistent with cross-ancestry/phenotype locus complexity. Colocalization (coloc.abf): PP.H4 ≥ 0.80 indicates strong evidence for a shared causal variant; IGF1R is PP.H2-dominant (a cis-eQTL signal without a detectable outcome association). Tissue by DESeq2 (in-house orbital RNA-seq, n = 1 control — exploratory). NS, not significant; —, not applicable/not computed. All values verified against source analysis files.*

**Table 3. All thirteen Bonferroni-significant druggable-gene-wide BBJ discovery hits, with cross-outcome colocalization and classification.**

| Gene | N IV | OR (95% CI) | BBJ *P* | Coloc PP.H4 (BBJ) | Coloc PP.H4 (FinnGen) | Classification |
|---|---|---|---|---|---|---|
| *HLA-A* | 4 | 1.98 (1.73–2.27) | 2.6×10⁻²³ | — | — | MHC region |
| *HLA-DQA2* | 1 | 2.39 (1.95–2.93) | 8.3×10⁻¹⁷ | — | — | MHC region |
| ***CTLA4*** | 1 | 0.18 (0.11–0.27) | 5.5×10⁻¹⁵ | 0.206 | 0.978 | Established autoimmune locus (positive control) |
| ***TSHR*** | 1 | 0.12 (0.07–0.21) | 1.1×10⁻¹⁴ | 0.951 | 0.986 | Established GD locus (susceptibility anchor) |
| *C4A* | 1 | 0.45 (0.35–0.57) | 2.2×10⁻¹⁰ | — | — | MHC region |
| *HSD3B7* | 1 | 0.29 (0.18–0.47) | 2.0×10⁻⁷ | 0.616 | 0.026 | chr16p11.2 LD cluster |
| *TUBB* | 2 | 1.49 (1.28–1.74) | 3.1×10⁻⁷ | — | — | MHC-region-linked |
| *VKORC1* | 1 | 0.13 (0.06–0.29) | 6.4×10⁻⁷ | 0.360 | 0.029 | chr16p11.2 LD cluster |
| *TNFSF14* | 1 | 0.63 (0.52–0.76) | 1.5×10⁻⁶ | 0.994 | 0.019 | Candidate (single-outcome coloc) |
| *PRSS36* | 1 | 6.00 (2.87–12.57) | 2.0×10⁻⁶ | 0.157 | 0.038 | chr16p11.2 LD cluster |
| *MAPKAPK5* | 4 | 0.98 (0.97–0.99) | 5.3×10⁻⁶ | <0.001 | 0.030 | Candidate (distinct-variant signal) |
| *PSMB8* | 1 | 1.25 (1.13–1.38) | 6.9×10⁻⁶ | — | — | MHC region |
| *IFNGR1* | 1 | 2.10 (1.51–2.91) | 9.4×10⁻⁶ | 0.989 | 0.019 | Candidate (single-outcome coloc) |

*Primary discovery estimates against BBJ Graves disease (Wald ratio for single-instrument loci, IVW for multi-instrument loci), ordered by discovery P value. OR per unit increase in genetically proxied expression. Colocalization (coloc.abf) was evaluated against the BBJ discovery and FinnGen TED-enriched outcomes; "—" denotes loci in the MHC region or otherwise not carried into formal two-outcome colocalization. The two established loci (TSHR, CTLA4) and MHC-region signals are recognized a priori. Among non-known candidates, TNFSF14 and IFNGR1 colocalized strongly in BBJ (PP.H4 = 0.994 and 0.989) but not in the TED-enriched outcome (PP.H4 = 0.019), and the chr16p11.2 cluster reflects a single LD block (Supplementary Table S5); no non-known candidate showed shared-variant colocalization in both outcomes, yielding no robust novel target.*

---

## Supplementary Material

**Supplementary Tables**

- **Table S1.** Extended cross-outcome detail for the thirteen Bonferroni-significant BBJ discovery hits. Whereas Table 3 summarizes the discovery (BBJ) effect sizes and colocalization, this table provides the per-hit Mendelian randomization estimates across all three outcomes (BBJ Graves disease, UKB hyperthyroidism, FinnGen Graves ophthalmopathy) together with the colocalization lead variant, supporting the cross-outcome reproducibility assessment.

| Gene | BBJ β (*P*) | UKB β (*P*) | FinnGen β (*P*) | Coloc lead SNP | Classification |
|---|---|---|---|---|---|
| *HLA-A* | +0.68 (2.6×10⁻²³) | −0.17 (0.47) | +0.18 (0.43) | — | MHC region |
| *HLA-DQA2* | +0.87 (8.3×10⁻¹⁷) | −0.25 (0.45) | +0.83 (1.8×10⁻⁷) | — | MHC region |
| *CTLA4* | −1.74 (5.5×10⁻¹⁵) | −1.57 (2.2×10⁻³) | −1.77 (1.2×10⁻²) | rs231811 | Established autoimmune locus |
| *TSHR* | −2.10 (1.1×10⁻¹⁴) | −2.44 (8.8×10⁻²⁸) | −2.33 (2.8×10⁻⁷) | rs179252 | Established GD locus |
| *C4A* | −0.80 (2.2×10⁻¹⁰) | −0.14 (0.24) | −1.02 (1.4×10⁻¹¹) | — | MHC region |
| *HSD3B7* | −1.23 (2.0×10⁻⁷) | −0.28 (0.020) | +0.27 (0.27) | rs4889606 | chr16p11.2 LD cluster |
| *TUBB* | +0.40 (3.1×10⁻⁷) | +1.04 (3.1×10⁻⁶) | +1.24 (5.9×10⁻⁵) | — | MHC-region-linked |
| *VKORC1* | −2.02 (6.4×10⁻⁷) | −0.37 (0.073) | — | rs34649473 | chr16p11.2 LD cluster |
| *TNFSF14* | −0.47 (1.5×10⁻⁶) | −0.30 (5.6×10⁻³) | +0.03 (0.85) | rs2291668 | Candidate (single-outcome coloc) |
| *PRSS36* | +1.79 (2.0×10⁻⁶) | +0.45 (0.014) | +0.44 (0.23) | rs78924645 | chr16p11.2 LD cluster |
| *MAPKAPK5* | −0.02 (5.3×10⁻⁶) | −0.19 (0.078) | −0.10 (0.16) | rs79271898 | Candidate (distinct-variant signal) |
| *PSMB8* | +0.22 (6.9×10⁻⁶) | −0.04 (0.55) | −0.07 (0.57) | — | MHC region |
| *IFNGR1* | +0.74 (9.4×10⁻⁶) | +0.30 (0.070) | +0.14 (0.68) | rs11754268 | Candidate (single-outcome coloc) |

  *Effect estimates (β, log-odds) by Wald ratio (single instrument) or inverse-variance weighted (multiple instruments). "—" for FinnGen *VKORC1* indicates no estimable instrument in that outcome; "—" in the colocalization column denotes MHC-region loci not carried into formal colocalization. Several discovery hits do not reproduce in the broader-phenotype (UKB) or TED-enriched (FinnGen) outcomes. Among non-MHC established anchors, *TSHR* and *CTLA4* showed consistent protective-direction estimates across all three outcomes; MHC-region signals were interpreted separately because of complex regional linkage disequilibrium.*
- **Table S2.** Candidate colocalization results across the discovery (BBJ Graves disease) and TED-enriched (FinnGen Graves ophthalmopathy) outcomes for non-known candidate hits. Posterior probabilities are from coloc.abf (H0 no association; H1 GWAS only; H2 eQTL only; H3 distinct causal variants; H4 shared causal variant).

| Gene | Outcome | Overlapping SNPs | PP.H0 | PP.H1 | PP.H2 | PP.H3 | PP.H4 | Top SNP |
|---|---|---|---|---|---|---|---|---|
| *TNFSF14* | BBJ Graves | 4,788 | 0.000 | 0.000 | 0.001 | 0.005 | 0.994 | rs2291668 |
| *TNFSF14* | FinnGen GO | 7,696 | 0.000 | 0.000 | 0.626 | 0.355 | 0.019 | rs2291668 |
| *IFNGR1* | BBJ Graves | 3,958 | 0.000 | 0.000 | 0.005 | 0.006 | 0.989 | rs11754268 |
| *IFNGR1* | FinnGen GO | 5,859 | 0.000 | 0.000 | 0.642 | 0.339 | 0.019 | rs11754268 |
| *MAPKAPK5* | BBJ Graves | 1,909 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | rs79271898 |
| *MAPKAPK5* | FinnGen GO | 2,979 | 0.000 | 0.000 | 0.674 | 0.296 | 0.030 | rs79271898 |
| *HSD3B7* | BBJ Graves | 1,625 | 0.000 | 0.000 | 0.000 | 0.384 | 0.616 | rs4889606 |
| *HSD3B7* | FinnGen GO | 3,095 | 0.000 | 0.000 | 0.733 | 0.240 | 0.026 | rs4889606 |
| *VKORC1* | BBJ Graves | 1,512 | 0.000 | 0.000 | 0.000 | 0.639 | 0.360 | rs34649473 |
| *VKORC1* | FinnGen GO | 2,950 | 0.000 | 0.000 | 0.740 | 0.231 | 0.029 | rs34649473 |
| *PRSS36* | BBJ Graves | 1,479 | 0.000 | 0.000 | 0.000 | 0.843 | 0.157 | rs78924645 |
| *PRSS36* | FinnGen GO | 2,820 | 0.000 | 0.000 | 0.742 | 0.220 | 0.038 | rs78924645 |

  *No non-known candidate showed shared-variant colocalization (PP.H4 ≥ 0.80) in both outcomes. TNFSF14 and IFNGR1 colocalized strongly in the BBJ discovery outcome but collapsed to an eQTL-only pattern (high PP.H2) in the TED-enriched outcome; the chromosome 16p11.2 candidates (HSD3B7, VKORC1, PRSS36) showed only weak or ambiguous colocalization. coloc.abf priors: p1 = p2 = 1×10⁻⁴, p12 = 1×10⁻⁵.*
- **Table S3.** *TSHR* locus fine-mapping and LD-clumping sensitivity. Independent-instrument counts at the *TSHR* locus under varying clumping thresholds, using ancestry-matched 1000 Genomes Phase 3 reference panels (European, n = 503; East Asian, n = 504). rs179252 (the primary instrument) is retained as the lead variant under every threshold.

| eQTL P threshold | Clumping *r*² | LD reference | Independent SNPs | Includes rs179252 |
|---|---|---|---|---|
| 5×10⁻⁸ | 0.001 | EUR | 1 | Yes |
| 5×10⁻⁸ | 0.001 | EAS | 2 | Yes |
| 5×10⁻⁸ | 0.01 | EUR | 1 | Yes |
| 5×10⁻⁸ | 0.01 | EAS | 4 | Yes |
| 5×10⁻⁸ | 0.1 | EUR | 6 | Yes |
| 5×10⁻⁸ | 0.1 | EAS | 10 | Yes |
| 5×10⁻⁶ | 0.001 | EUR | 2 | Yes |
| 5×10⁻⁶ | 0.001 | EAS | 3 | Yes |

  *At the primary instrument-selection threshold (P < 5×10⁻⁸, r² < 0.001), the European reference yields a single independent instrument (rs179252); additional "independent" SNPs appear only at progressively relaxed r² thresholds and reflect the wider LD structure of the locus.*

SuSiE fine-mapping (susieR; 1000 Genomes Phase 3 European reference) resolved the *TSHR* *cis*-eQTL signal into a **single 95% credible set** (purity = 0.993) spanning a ~116-kb window on chromosome 14 (hg19). The credible-set variants and their posterior inclusion probabilities (PIP) are:

| Variant | Chr | Position (hg19) | PIP | *Z* |
|---|---|---|---|---|
| rs11603529 | 14 | 81,440,788 | 0.458 | −13.06 |
| rs179248 | 14 | 81,433,038 | 0.436 | −13.03 |
| rs2284720 | 14 | 81,489,065 | 0.052 | −12.45 |
| rs10137255 | 14 | 81,373,516 | 0.021 | −12.16 |
| rs2284722 | 14 | 81,442,129 | 0.020 | −12.15 |

  *A single credible set (log₁₀ Bayes factor = 88.2; purity = 0.993) was identified, with no second independent signal. The primary MR instrument rs179252 (chr14:81,435,985; marginal P = 2.86×10⁻⁴⁰) carries SuSiE PIP = 1.0 in the European panel and serves as the locus anchor. In a complementary per-variant conditional/LD assessment across both ancestries, every other lead variant at the locus remained in substantial LD with rs179252 (minimum r² = 0.808 in East Asian, 0.965 in European) and none qualified as an independent instrument (usable_as_independent_IV = FALSE for all). These analyses jointly support treatment of TSHR as a single-instrument, locus-level anchor rather than a multi-instrument exposure.*
- **Table S4.** Instrument selection and strength for the backbone genes and all thirteen discovery hits. Instruments were genome-wide significant *cis*-eQTLs (P < 5×10⁻⁸); the secondary threshold (P < 5×10⁻⁶) count is shown for reference. Across all 2,545 MR-testable genes (6,136 instruments), the minimum F-statistic was 14.2, with no weak instruments (all F > 10).

| Gene | Chr | Instruments (P < 5×10⁻⁸) | Instruments (P < 5×10⁻⁶) | Single-instrument | Druggability tier |
|---|---|---|---|---|---|
| *TSHR* | 14 | 1 | 710 | Yes | Tier 1 |
| *IGF1R* | 15 | 4 | 511 | No | Tier 1 |
| *CTLA4* | 2 | 2 | 333 | No | Tier 1 |
| *HLA-A* | 6 | 7 | 10,315 | No | Tier 1 |
| *HLA-DQA2* | 6 | 2 | 9,971 | No | Tier 3A |
| *C4A* | 6 | 5 | 8,175 | No | Tier 3A |
| *HSD3B7* | 16 | 2 | 736 | No | Tier 3B |
| *TUBB* | 6 | 3 | 1,476 | No | Tier 1 |
| *VKORC1* | 16 | 1 | 489 | Yes | Tier 1 |
| *TNFSF14* | 19 | 1 | 258 | Yes | Tier 1 |
| *PRSS36* | 16 | 1 | 333 | Yes | Tier 3A |
| *MAPKAPK5* | 12 | 4 | 1,349 | No | Tier 2 |
| *PSMB8* | 6 | 1 | 1,060 | Yes | Tier 1 |
| *IFNGR1* | 6 | 1 | 340 | Yes | Tier 1 |

  *Druggability tiers follow Finan et al. (Tier 1, approved drug targets or those in clinical development; Tier 2, targets of small molecules or biotherapeutics with druggable similarity; Tier 3, secreted or extracellular proteins and members of druggable gene families). Several MHC-region genes (e.g., HLA-A, HLA-DQA2, C4A) have very large secondary-threshold instrument counts reflecting the extended LD of the region; primary analyses used the stringent P < 5×10⁻⁸ instruments after clumping. The "Instruments (P < 5×10⁻⁸)" column reports the number of independent cis-eQTL instruments available after clumping; the number of instruments actually used in MR against a given outcome (N IV in Tables 2 and 3) may be smaller, reflecting outcome-specific SNP availability and harmonization (absence of a variant in the outcome GWAS, or removal of palindromic variants).*
- **Table S5.** Chromosome 16p11.2 candidate cluster. Three discovery hits (*HSD3B7*, *VKORC1*, *PRSS36*) map to a narrow (~143-kb) window on chromosome 16p11.2 and behave as a single linkage-disequilibrium signal rather than three independent associations.

| Gene | Chr | Position (hg19) | BBJ discovery *P* | Coloc PP.H4 (BBJ) | Coloc PP.H4 (FinnGen) | Top SNP |
|---|---|---|---|---|---|---|
| *HSD3B7* | 16 | 31,011,183 | 2.0×10⁻⁷ | 0.616 | 0.026 | rs4889606 |
| *VKORC1* | 16 | 31,066,538 | 6.4×10⁻⁷ | 0.360 | 0.029 | rs34649473 |
| *PRSS36* | 16 | 31,154,358 | 2.0×10⁻⁶ | 0.157 | 0.038 | rs78924645 |

  *All three lead variants lie within a 143.2-kb window on chromosome 16p11.2 (hg19 chr16:31,011,183–31,154,358) and none reaches strong colocalization (PP.H4 ≥ 0.80) in either outcome. Pairwise linkage disequilibrium in the 1000 Genomes Phase 3 European reference (n = 503) confirms that the cluster does not comprise three independent signals: rs4889606 (HSD3B7) and rs34649473 (VKORC1) are in strong LD (r² = 0.90), while rs78924645 (PRSS36) is in weaker LD with each (r² = 0.19 and 0.20, respectively). Their co-significance in discovery therefore reflects shared regional linkage disequilibrium across a single locus rather than three independent causal effects.*
- **Table S6.** Mendelian randomization sensitivity analyses for the backbone genes across all three outcomes. IVW, inverse-variance weighted; WM, weighted median; Wmode, weighted mode. MR-Egger intercept, Cochran's Q, and multi-estimator comparisons require multiple instruments and are not estimable (NA) at single-instrument loci.

| Gene | Outcome | N IV | IVW/Wald *P* | WM *P* | Wmode *P* | MR-Egger intercept *P* | Cochran's Q *P* |
|---|---|---|---|---|---|---|---|
| *TSHR* | BBJ Graves | 1 | 1.1×10⁻¹⁴ | NA | NA | NA | NA |
| *TSHR* | UKB hyperthyroid | 1 | 8.8×10⁻²⁸ | NA | NA | NA | NA |
| *TSHR* | FinnGen GO | 1 | 2.8×10⁻⁷ | NA | NA | NA | NA |
| *IGF1R* | BBJ Graves | 4 | 0.021 | 0.029 | 0.118 | 0.663 | 0.902 |
| *IGF1R* | UKB hyperthyroid | 4 | 0.012 | 0.002 | 0.040 | 0.562 | 0.250 |
| *IGF1R* | FinnGen GO | 3 | 0.182 | 0.124 | 0.241 | 0.433 | 0.309 |
| *CTLA4* | BBJ Graves | 1 | 5.5×10⁻¹⁵ | NA | NA | NA | NA |
| *CTLA4* | UKB hyperthyroid | 2 | 0.002 | NA | NA | NA | 0.001 |
| *CTLA4* | FinnGen GO | 2 | 0.012 | NA | NA | NA | 0.049 |

  *For IGF1R (the only multi-instrument backbone gene with a full estimator set), all sensitivity estimators are directionally concordant with the IVW estimate and the MR-Egger intercept is non-significant in every outcome (P > 0.43), indicating no detectable directional pleiotropy. CTLA4 showed significant between-instrument heterogeneity in the UKB and FinnGen analyses (Cochran's Q P = 0.001 and 0.049), consistent with cross-ancestry/phenotype locus complexity. Steiger directionality was not estimable for these summary-level single- or limited-instrument loci and was therefore not applied. Single-instrument loci (TSHR all outcomes; CTLA4 in BBJ) are supported by colocalization and fine-mapping rather than instrument-based pleiotropy diagnostics.*
- **Table S7.** STROBE-MR reporting checklist. Each item maps to where it is addressed in the manuscript (page/section references to be finalized at proof stage).

| # | STROBE-MR Item | Addressed in |
|---|---|---|
| 1 | MR stated in title/abstract; design described | Title, Abstract |
| 2 | Background; rationale for genetic instruments | Introduction |
| 3 | Explicit hypotheses / objectives | Introduction (end) |
| 4 | Study design & key elements | Methods; Table 1; Figure 1 |
| 5 | Data sources, settings, populations | Methods; Table 1 |
| 6a | Instrument selection criteria | Methods (eQTLGen; Finan druggable genome; cis P < 5×10⁻⁸) |
| 6b | MR assumptions (relevance, independence, exclusion) | Methods; Discussion (limitations) |
| 7 | Variable definition (exposure, outcome) | Methods |
| 8 | Effect measures & scale (log-odds; UKB rescaling) | Methods; Table 1 footnote |
| 9 | Statistical methods (IVW, MR-Egger, WM, Wmode; coloc.abf) | Methods |
| 10 | Sensitivity analyses (pleiotropy, heterogeneity, Steiger) | Methods; Table S6 |
| 11 | Software & package versions | Methods; below |
| 12 | Descriptive data (instrument counts, F-statistics, detectable effect range) | Results; Tables S4, S8 |
| 13 | Main MR estimates with CIs | Results; Tables 2, 3; Figures 2, 3 |
| 14 | Sensitivity / robustness results | Results; Tables S6, S8 |
| 15 | Colocalization results | Results; Table 2; Figure 3A; Tables S2, S5 |
| 16 | Additional analyses (tissue expression) | Results; Figure 3C |
| 17 | Key results in context | Discussion |
| 18 | Limitations (blood eQTL vs orbital tissue; European eQTL vs East Asian discovery; broader UKB phenotype; single-IV *TSHR*; limited power in the 858-case TED-enriched outcome; TED-enriched outcome is not a within-Graves TED-specific contrast; coloc single-causal-variant assumption and single European allele-frequency reference; no formal colocalization power analysis; non-uniform behaviour of the single positive control; exploratory n = 1-control tissue analysis) | Discussion (limitations paragraph) |
| 19 | Interpretation (susceptibility-anchor vs effector) | Discussion |
| 20 | Generalizability | Discussion |
| 21 | Funding | Declarations |
| 22 | Data availability | Declarations (public summary statistics: eQTLGen, GWAS Catalog, FinnGen) |

  *Software and package versions: R 4.3.3; TwoSampleMR v0.7.4; ieugwasr; coloc v5.2.3; susieR v0.14.2; DESeq2 v1.42.1; PLINK v1.9. Reporting follows the STROBE-MR statement (Skrivankova et al., JAMA 2021).*

- **Table S8.** Detectable effect range of the druggable-gene-wide screen. Minimum effect detectable with 80% power was computed per gene and outcome from the empirical standard error of the primary MR estimate as |β|min = (*z*1−α/2 + *z*0.80) × SE, with α = 1.965×10⁻⁵ (discovery) or α = 0.05 (replication, TED-enriched).

  **Panel A. Screen-wide detectable range.**

| Outcome | Significance threshold | Genes | Median \|β\|min (IQR) | Median detectable OR (IQR) | Genes powered for OR ≥ 1.5 | OR ≥ 2.0 | OR ≥ 3.0 |
|---|---|---|---|---|---|---|---|
| BBJ Graves (discovery) | *P* < 1.965×10⁻⁵ | 2,235 | 0.94 (0.54–1.65) | 2.55 (1.71–5.19) | 14.6% | 35.6% | 57.7% |
| UKB hyperthyroid (replication) | *P* < 0.05 | 2,506 | 0.38 (0.23–0.65) | 1.46 (1.26–1.92) | 53.7% | 77.4% | 93.1% |
| FinnGen Graves ophthalmopathy (TED-enriched) | *P* < 0.05 | 2,481 | 0.75 (0.46–1.34) | 2.12 (1.58–3.82) | 19.4% | 45.5% | 67.3% |

  **Panel B. Backbone genes: observed effect versus gene-specific detection threshold.**

| Gene | Outcome | Observed β | SE | \|β\|min (80% power) | Detectable OR | Observed effect detectable? |
|---|---|---|---|---|---|---|
| ***TSHR*** | BBJ Graves | −2.096 | 0.271 | 1.386 | 4.00 | Yes |
| | UKB hyperthyroid | −2.436 | 0.223 | 0.625 | 1.87 | Yes |
| | FinnGen GO | −2.331 | 0.454 | 1.272 | 3.57 | Yes |
| ***IGF1R*** | BBJ Graves | +0.446 | 0.194 | 0.989 | 2.69 | No |
| | UKB hyperthyroid | +0.299 | 0.119 | 0.333 | 1.39 | No |
| | FinnGen GO | +0.342 | 0.256 | 0.718 | 2.05 | No |
| ***CTLA4*** | BBJ Graves | −1.740 | 0.223 | 1.137 | 3.12 | Yes |
| | UKB hyperthyroid | −1.569 | 0.513 | 1.436 | 4.20 | Yes |
| | FinnGen GO | −1.768 | 0.702 | 1.968 | 7.15 | No |

  *The screen was well powered for large expression-mediated effects and recovered both established anchors (TSHR, CTLA4) at magnitudes exceeding their gene-specific thresholds. IGF1R effects fell below the corresponding thresholds in every outcome, so the absence of study-wide significance at IGF1R is expected for an effect of this size and is not evidence against a genetic contribution; IGF1R was therefore evaluated as a backbone gene across multiple evidence layers rather than on discovery significance alone. Conversely, the absence of a robust novel target constrains additional large, expression-colocalized druggable effects but not moderate ones. Values are computed from the standard errors of the locked primary MR estimates (reproduced by `scripts/22_power_analysis_null_screen.py`).*

**Supplementary Figures**

**Figure S1. *TSHR* fine-mapping and LD rationale for single-instrument selection.** SuSiE fine-mapping of the *TSHR* *cis*-eQTL signal (eQTLGen) using 1000 Genomes Phase 3 European and East Asian LD reference panels. **(A)** Genomic positions of credible-set lead SNPs at the *TSHR* locus (chr14, hg19); rs179252 (star) is the primary instrument used in all MR analyses, and all credible-set leads fall within a ~15-kb window. **(B)** Ancestry-matched LD (*r*²) of each secondary credible-set lead with rs179252; every secondary lead remains in substantial LD with rs179252 in its corresponding ancestry-matched panel (minimum *r*² = 0.808 in East Asian, 0.965 in European), exceeding the *r*² = 0.8 threshold (dashed line). Because no secondary signal is independent of rs179252, *TSHR* is retained as a single-instrument/locus-level anchor; formal multi-instrument pleiotropy diagnostics are therefore not applicable at this locus.

**Figure S2. Outcome-dependent colocalization of candidate discovery hits.** Stacked posterior probabilities (coloc.abf) for non-known candidate genes against BBJ Graves disease and FinnGen Graves ophthalmopathy, showing PP.H2 (*cis*-eQTL association only), PP.H3 (distinct causal variants), and PP.H4 (shared causal variant); the dashed line marks PP.H4 = 0.8. *TNFSF14* and *IFNGR1* exceed the shared-variant threshold in the BBJ discovery outcome (PP.H4 = 0.994 and 0.989) but not in the TED-enriched outcome (PP.H4 = 0.019 for both), and no candidate colocalized in both outcomes, yielding no robust novel target.
