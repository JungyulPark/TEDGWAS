**Running head:** TSHR and IGF1R genetics in Graves disease

# Genetic susceptibility at TSHR and IGF1R in Graves disease and thyroid eye disease

Jungyul Park¹, Min-Seon Kim², Kyung-Hwa Shin³⁎, Suk-Woo Yang¹⁎

¹ Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
² Department of Ophthalmology, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea  
³ Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea  

\* These authors contributed equally as corresponding authors.

**Corresponding author:** Suk-Woo Yang, MD, PhD, Department of Ophthalmology, Seoul St. Mary's Hospital, College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea. E-mail: yswoph@catholic.ac.kr; Tel: +82-2-2258-2847.  
**Co-corresponding author:** Kyung-Hwa Shin, MD, Department of Laboratory Medicine, Pusan National University Hospital, Busan, Republic of Korea.

---

## Abstract

Teprotumumab targets IGF-1R in thyroid eye disease (TED), but whether *IGF1R* carries the genetic architecture of a susceptibility locus, as the established autoantigen *TSHR* does, remains unclear. We applied a druggable-gene-wide Mendelian randomization (MR) and colocalization framework to separate susceptibility from effector axes in Graves disease (GD) and TED. eQTLGen blood *cis*-eQTL instruments were used for two-sample MR across 2,544 druggable genes. Outcomes were Biobank Japan (BBJ) Graves disease (discovery), UK Biobank (UKB) hyperthyroidism (replication), and FinnGen R12 Graves ophthalmopathy (TED-enriched sensitivity). Bonferroni-significant hits underwent cross-outcome evaluation and Bayesian colocalization against all three outcomes, with prior sensitivity analysis. Thirteen genes reached Bonferroni significance in discovery. *TSHR* showed protective MR effects across all three outcomes and strong colocalization in both Graves disease phenotypes (PP.H4=0.951 and 0.986; top-posterior variant rs179252), but not in the broader hyperthyroidism phenotype, where distinct causal variants were favoured (PP.H3=0.774). *IGF1R* showed directionally consistent risk-direction MR estimates, reaching nominal significance in BBJ and UKB but not in the smaller FinnGen outcome, and did not reach the colocalization threshold in any outcome (PP.H4=0.073, 0.404, 0.032). *CTLA4*, the positive control, colocalized in both European outcomes but not in discovery. No novel candidate satisfied the dual-phenotype colocalization criterion. These results support *TSHR* as an expression-colocalized susceptibility anchor for GD, while shared-variant evidence at *IGF1R* remains unresolved. The differing genetic profiles do not establish an exclusive effector role for *IGF1R* or exclude a modest inherited contribution.

**Keywords:** Graves disease; Thyroid eye disease; *TSHR*; *IGF1R*; Mendelian randomization; Colocalization

---

## Introduction

Thyroid eye disease (TED) is the most common extrathyroidal manifestation of Graves disease (GD), an autoimmune disorder in which orbital inflammation and tissue remodeling can impair quality of life, ocular function, and vision [1, 2]. For decades, treatment relied largely on corticosteroids, orbital radiotherapy, and surgery, with limited disease-modifying options. The introduction of teprotumumab, a monoclonal antibody targeting the insulin-like growth factor 1 receptor (IGF-1R), changed this therapeutic landscape: randomized trials demonstrated improvement in proptosis, clinical activity, and diplopia in active TED [3], establishing IGF-1R blockade as a viable therapeutic strategy. However, clinical response to IGF-1R blockade does not by itself resolve whether *IGF1R* represents a genetically anchored susceptibility locus for GD/TED or a pharmacologic effector target that can be modulated without being a primary driver of disease risk.

*TSHR* and IGF-1R occupy different but biologically connected positions in GD/TED pathogenesis and therapy. *TSHR* is the central autoantigen of GD: loss of immune tolerance to the thyrotropin receptor gives rise to TSHR-directed stimulating autoantibodies, and *TSHR* is among the most consistently replicated genetic susceptibility loci for the disease [2]. IGF-1R, by contrast, has emerged as both a mechanistic candidate and a therapeutic target in TED [4]. Experimental studies have reported that IGF-1R can associate with TSHR, including through β-arrestin 1–dependent receptor scaffolding, and that the two receptors engage in functional crosstalk in orbital fibroblasts [5, 6], providing a biological rationale for IGF-1R blockade. These observations leave open a fundamental question of interpretation: whether *IGF1R* should be regarded as an independent genetic susceptibility locus, a signaling partner of TSHR, or a pharmacologic effector node that is therapeutically tractable without being a primary determinant of inherited disease risk. This distinction is difficult to resolve from tissue expression or therapeutic-response data alone, because both can reflect disease-state activation, cellular remodeling, or downstream pathway modulation rather than germline susceptibility architecture.

Human genetics offers a way to address this question that is less susceptible to confounding and reverse causation than observational expression studies. Mendelian randomization (MR) uses germline genetic variants as instruments to infer whether genetically proxied molecular exposures—here, *cis*-regulated gene expression—are associated with disease liability in a manner consistent with causality [7], while colocalization tests whether the genetic signals for expression and disease share an underlying causal variant rather than reflecting distinct linked variants [8]. Applied across a predefined druggable-gene set, this combination can prioritize loci systematically without restricting inference to a small number of preselected candidates. Prior MR and multi-omics studies have begun to nominate molecular contributors to GD, including immune-related and other molecular signals [9, 10], but they have not directly resolved whether *TSHR* and *IGF1R* occupy different positions as genetically anchored susceptibility versus pharmacologic effector axes. Integrating MR with colocalization can further distinguish expression-colocalized susceptibility signals from non-colocalized or downstream effector biology.

We therefore applied a druggable-gene-wide *cis*-eQTL MR and colocalization framework across a defined GD/TED outcome hierarchy: Biobank Japan Graves disease as the primary discovery outcome, UK Biobank hyperthyroidism as replication, and FinnGen Graves ophthalmopathy as a TED-enriched sensitivity outcome. We evaluated *TSHR*, *IGF1R* and *CTLA4* — the autoantigen, the therapeutic target, and a positive-control autoimmune locus — across every outcome. Our aim was to distinguish expression-colocalized susceptibility anchors from pharmacologic effector axes across druggable genes relevant to GD/TED biology and therapeutic interpretation, rather than treating discovery MR hits alone as sufficient evidence for target nomination.

---

## Methods

### Study design and data sources

We performed a druggable-gene-wide Mendelian randomization (MR) study integrating MR and colocalization to distinguish genetically anchored susceptibility from pharmacologic effector biology in Graves disease (GD) and thyroid eye disease (TED) (Figure 1), following the STROBE-MR guideline [11] (Table S7). Exposures were blood *cis*-eQTL instruments from the eQTLGen Consortium (predominantly European; up to 31,684 participants) [12].

Disease outcomes followed a defined analytical hierarchy: Graves disease from Biobank Japan (BBJ; East Asian; 2,809 cases) served as the primary discovery outcome [13], hyperthyroidism from UK Biobank (UKB; European; 3,731 cases) as a cross-ancestry, broader-phenotype replication outcome [14], and Graves ophthalmopathy from FinnGen R12 (European; 858 cases) as the TED-enriched sensitivity outcome [15] (Table 1). Because adequately powered TED-specific GWAS are unavailable, this design prioritized discovery in Graves disease, assessed generalizability in a broader hyperthyroidism phenotype, and finally examined signal recovery in a smaller TED-enriched outcome. An in-house orbital RNA-seq dataset (Korean; four TED, one control) was included for descriptive context only and was not used for instrument selection. Three genes were selected on biological grounds as the interpretive backbone on biological and therapeutic grounds rather than on screen results: *TSHR* (the Graves disease autoantigen), *IGF1R* (the target of teprotumumab), and *CTLA4* (a positive-control autoimmune locus). All three were evaluated across both evidence layers — MR and colocalization — against every outcome, irrespective of their discovery *P* values, with an additional LD-clumping sensitivity assessment at the single-instrument *TSHR* locus.

### Druggable gene set and instrument selection

Analyses were restricted to the druggable genome (Finan et al. [16]; 4,462 genes, GRCh37). For each gene, *cis*-eQTL variants within ±1 Mb reaching genome-wide significance (*P* < 5×10⁻⁸) were eligible as instruments; a relaxed threshold (*P* < 5×10⁻⁶) was used only for secondary or sensitivity analyses and did not define the primary denominator. Exposure effect sizes were reconstructed from the reported Z-scores and per-SNP sample sizes using 1000 Genomes Phase 3 European allele frequencies (n = 503) [17]; the consequences of using the reference panel rather than eQTLGen's own allele-frequency file are considered in the Limitations. Independent instruments were obtained by LD clumping (*r*² < 0.001, 10-Mb window) using local PLINK v1.9 [18] with the same reference panel. Instrument strength was quantified as *F* = *Z*². Selection at *P* < 5×10⁻⁸ implies |*Z*| > 5.45 and therefore *F* > 29.7; the observed minimum was 29.7 (no weak instruments).

After clumping, 2,544 of 4,462 genes retained at least one valid instrument and formed the MR-testable denominator (6,135 SNP-level instruments; Table S4), which defined the study-wide significance threshold (Bonferroni *P* < 1.965×10⁻⁵; 0.05/2,544).

### Outcome data and harmonization

Exposure and outcome data were harmonized with TwoSampleMR [19] (harmonise_data, action = 2), aligning alleles to a common effect allele and inferring strand for complementary alleles where possible; palindromic variants with intermediate frequency (MAF > 0.42) were removed to avoid strand ambiguity.

Outcome estimates were placed on a common log-odds scale: BBJ and FinnGen were used directly, while UKB hyperthyroidism (reported on a linear mixed-model scale) was rescaled to log-odds using case prevalence before MR. Variants lacking usable outcome estimates or failing harmonization were omitted without imputation. Source-study GWAS covariate adjustments were retained; individual-level covariates were unavailable for common re-adjustment. Given differences in ancestry, phenotype, and effect-scale reconstruction across cohorts, cross-outcome interpretation emphasized direction, replication, and colocalization rather than exact equality of effect magnitudes.

### Mendelian randomization

Two-sample MR was performed for each gene against all three outcomes. Interpretation as a causal effect requires instrument relevance, independence from exposure–outcome confounders, and no effect on the outcome through pathways other than the proxied expression (exclusion restriction). These assumptions are not established by statistical significance or colocalization alone. For single-instrument loci, the Wald ratio was used; for multi-instrument loci, the inverse-variance-weighted (IVW) method was the primary estimator, with MR-Egger [20], weighted median [21], and weighted mode as sensitivity estimators. Directional pleiotropy was assessed using the MR-Egger intercept, between-instrument heterogeneity using Cochran's Q, while Steiger directionality testing [22] was not applied (Table S6). Discovery significance was defined at the study-wide Bonferroni threshold (*P* < 1.965×10⁻⁵; 0.05/2,544). Replication and TED-enriched support were evaluated using direction-consistent nominal association (*P* < 0.05).

Because a screen that advances no novel target is interpretable only in relation to the effect sizes it could have detected, we quantified the detectable effect range of the screen itself. For each gene and outcome, the minimum effect detectable with 80% power at the applicable significance threshold was computed from the empirical standard error of the primary MR estimate as |β|min = (*z*1−α/2 + *z*0.80) × SE, using α = 1.965×10⁻⁵ for the discovery outcome and α = 0.05 for the replication and TED-enriched outcomes. This design-based approach uses the standard errors actually obtained after harmonization rather than reconstructed variance components, and is reported for all genes and for the three backbone loci (Table S8).

### Instrument resolution at TSHR and colocalization

To assess the stability of instrument selection at *TSHR*, LD clumping at the locus was repeated across a grid of *P*-value and *r*² thresholds in both 1000 Genomes Phase 3 panels (European, *n* = 503; East Asian, *n* = 504). At the primary threshold (*P* < 5×10⁻⁸, *r*² < 0.001) the European panel — the one matching the eQTL data and therefore used for selection — yielded exactly one independent instrument, rs179252. rs179252 was among the index variants under all twelve thresholds tested, though not the only one at relaxed thresholds or in the East Asian panel (Table S3). *TSHR* was accordingly treated as a single-instrument, locus-level anchor rather than reconstructed as a multi-instrument exposure (Table S3).

For prioritized loci, we tested whether the *cis*-eQTL and outcome-association signals shared a causal variant using coloc.abf [8] (coloc) with the available gene-specific eQTLGen *cis*-variants, under the single-causal-variant assumption and default priors (p1 = p2 = 10⁻⁴, p12 = 10⁻⁵) [23]. In this implementation dataset 1 was the outcome GWAS and dataset 2 the *cis*-eQTL, so the posterior probabilities distinguish association with the outcome only (H1), with the *cis*-eQTL only (H2), with both traits at distinct causal variants (H3), and with both traits at a shared causal variant (H4). Colocalization used beta and varbeta for both datasets. Outcome effect estimates and standard errors came from each GWAS; exposure estimates were reconstructed from eQTLGen Z-scores, per-variant sample sizes and 1000 Genomes European frequencies, with expression variance set to one (sdY = 1). Variants were matched by rsID and allele pair, avoiding coordinate-window comparisons between GRCh37 eQTL and GRCh38 outcome files; duplicate rsIDs and allele-pair mismatches were removed. Outcome allele frequencies were not used in these Bayes factors. Strong colocalization was defined as PP.H4 ≥ 0.80, with prior sensitivity at p12 = 10⁻⁶ and 5×10⁻⁶ (Tables S2, S9; Figure S2).

### Chromosome 16p11.2 follow-up

The three chromosome 16p11.2 discovery hits were evaluated with pairwise LD in the 1000 Genomes East Asian panel (n = 504). GCTA-COJO applied stepwise conditional selection to the available BBJ regional summary statistics using the same reference and a selection threshold of *P* < 5×10⁻⁶; marginal and conditional associations were compared (Table S5). This follow-up assessed regional dependence rather than identifying a specific causal gene.

### Orbital tissue transcriptomic analysis

Exploratory tissue-level support was assessed using an in-house bulk RNA-seq dataset of orbital tissue from a Korean cohort comprising four TED patients and one control. Technical replicates were collapsed to the biological-sample level and abundance expressed as transcripts per million (TPM). Descriptive log2 fold changes were calculated as log2[(mean TED TPM + 0.01)/(control TPM + 0.01)]. With one control, control-side biological variance cannot be estimated, so differential-expression tests are not used for inference or reported in this manuscript; abundance is described for the three backbone genes only (Figure S1). No transcriptomic catalogue, differential-expression list or pathway-enrichment result is presented, the dataset was not used for instrument selection, and it contributes to no inference here.

### Software, reporting, and data availability

Analyses used R 4.3.3 with TwoSampleMR (v0.7.4), ieugwasr, and coloc (v5.2.3); LD operations used PLINK v1.9. Reporting is mapped to the STROBE-MR guideline, including limitations in reporting coverage (Table S7). No prospective registration identifier is available for this analysis. Public eQTL and GWAS summary statistics are available from their original repositories, including eQTLGen, the GWAS Catalog, and FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions. The study used de-identified summary statistics and an institutionally approved orbital tissue dataset (IRB approval, Pusan National University Hospital, 2104-018-102).

---

## Results

### Druggable-gene-wide MR discovery

The attrition of the screen and the resulting contrast between the two backbone receptors are summarized in **Figure 1**. The analysis followed the defined outcome hierarchy of BBJ Graves disease (discovery), UKB hyperthyroidism (replication), and FinnGen R12 Graves ophthalmopathy (TED-enriched sensitivity), with data sources in **Table 1**.

Of 4,462 druggable genes, 2,544 had at least one valid *cis*-eQTL instrument in eQTLGen and formed the MR-testable denominator; across these genes, 6,135 SNP-level instruments were retained, all exceeding conventional instrument-strength criteria (minimum F-statistic = 29.7, the value implied by the *P* < 5×10⁻⁸ selection threshold; **Table S4**). Primary MR estimates against the BBJ Graves disease discovery outcome were obtainable for 2,234 genes.

Applying a Bonferroni threshold corrected for the MR-testable denominator (*P* < 1.965×10⁻⁵; 0.05/2,544), 13 genes reached significance in discovery (**Figure 2**). These hits included the established Graves disease autoantigen *TSHR* and the known autoimmune locus *CTLA4*, both with strongly significant protective-direction estimates, together with additional loci that were carried forward for replication and colocalization-based filtering. The discovery volcano (**Figure 2**) shows the distribution of BBJ MR effect estimates and the 13 Bonferroni-significant hits; full per-hit statistics, cross-outcome colocalization, and classification are provided in **Table 3** (with extended detail in **Table S1**).

### Cross-outcome classification and colocalization filtering

The 13 discovery hits were not treated as putative targets on the strength of the discovery MR alone. Instead, each was classified according to prior biological status, genomic context, cross-outcome support, and colocalization with both the discovery and TED-enriched signals. Two hits—*TSHR* and *CTLA4*—were recognized as established Graves disease and autoimmune loci, respectively, and are examined in detail below. The remaining hits were evaluated as non-known candidates, including a cluster on chromosome 16p11.2 (*HSD3B7*, *VKORC1*, and *PRSS36*) whose top eQTL instruments were co-located within a narrow approximately 143-kb linkage-disequilibrium window (Table S5).

Candidate colocalization was outcome-dependent (Table S2; Figure S2). *TNFSF14* and *IFNGR1* colocalized strongly with the BBJ discovery signal (PP.H4 = 0.994 and 0.989) but reached the threshold in neither secondary outcome (UKB 0.237 and 0.052; TED-enriched 0.017 and 0.020), where the posteriors were dominated by the *cis*-eQTL-only and distinct-variant hypotheses. *TNFSF14* retains a direction-consistent nominal MR association in UKB (*P* = 0.0056), so its exclusion reflects the colocalization criterion, not an absence of association. *MAPKAPK5* was distinct-variant-dominant in discovery (PP.H3 = 1.000), and the chr16p11.2 candidates reached at most PP.H4 = 0.636 in discovery and no more than 0.115 in either secondary outcome.

Overall, no non-known candidate showed strong colocalization in both BBJ Graves disease and the TED-enriched outcome. Thus, after cross-outcome colocalization filtering, no novel candidate satisfied the dual-phenotype criterion. This filtering guarded against overinterpreting single-outcome or LD-driven discovery signals.

The MR detection range varied substantially across genes (**Table S8**). At the study-wide discovery threshold and 80% power, the median detectable odds ratio was approximately 2.55 per unit genetically proxied expression (interquartile range 1.71–5.18). Only 35.6% of genes had at least 80% power at OR 2.0 and 14.6% at OR 1.5. In the TED-enriched outcome at nominal significance, the median detectable OR was approximately 2.12. TSHR and CTLA4 exceeded their gene-specific discovery thresholds, but these calculations quantify MR association power, not power to pass colocalization or the combined filter. The absence of a qualifying novel candidate therefore cannot exclude other large or moderate effects.

### TSHR: an expression-colocalized susceptibility anchor

Among the established hits, *TSHR* showed the strongest convergence between the MR and colocalization layers of the framework (Table 2; Figure 3). In MR, the *TSHR* *cis*-eQTL instrument was associated with protective-direction estimates across all three outcomes: BBJ Graves disease (β = −2.10, *P* = 1.1×10⁻¹⁴), UKB hyperthyroidism (β = −2.44, *P* = 8.8×10⁻²⁸), and the TED-enriched outcome (β = −2.33, *P* = 2.8×10⁻⁷). The consistency of direction and approximate magnitude across the independently ascertained Graves disease and hyperthyroidism outcomes argues against an outcome-specific artifact (Figure 3B). The concordant estimate in the TED-enriched outcome is interpreted more conservatively, because Graves ophthalmopathy cases are ascertained among patients with Graves disease and compared with population controls, so that outcome substantially re-measures shared Graves disease susceptibility rather than liability specific to orbital involvement (see Limitations).

Colocalization supported a shared causal variant between *TSHR* expression and disease association in the Graves disease discovery outcome and in the TED-enriched outcome (PP.H4 = 0.951 and 0.986; rs179252 the top-posterior variant in both; Figure 3A), but not in the broader UKB hyperthyroidism phenotype, where the posterior instead favoured distinct causal variants (PP.H3 = 0.774; PP.H4 = 0.226; lead variant rs1023586). A weak outcome association alone does not explain this discordance: UKB carried a strong locus-wide association (minimum *P* = 1.4×10⁻⁴³), and 5,667 variants overlapped the eQTL data. Formal colocalization power was not evaluated. We therefore report *TSHR* colocalization as specific to the two Graves disease phenotypes rather than as a property of thyroid autoimmunity in general (Table 2).

In the in-house orbital dataset, *TSHR* transcript abundance was higher in the four TED samples than in the single control (log2 fold-change = +2.59; Figure S1). With one control, this dataset cannot support differential-expression inference; it is reported as a descriptive observation and no statistical claim is drawn from it. Together, the genetic results position *TSHR* as an expression-colocalized upstream susceptibility anchor for Graves disease whose signal is also recovered in a TED-enriched Graves ophthalmopathy outcome, although a TED-specific effect separable from Graves disease liability is not established by these data. The apparent contrast between the protective-direction germline MR estimate and higher orbital *TSHR* abundance in disease-state tissue is addressed in the Discussion.

### IGF1R: nominal genetic association without expression colocalization

In contrast to *TSHR*, *IGF1R*, the therapeutic target of teprotumumab, showed a markedly different evidence profile (Table 2; Figure 3). *IGF1R* was instrumented by multiple independent *cis*-eQTL variants, with four instruments in the discovery and replication analyses. The inverse-variance-weighted MR estimate indicated a nominal positive-direction association with BBJ Graves disease (β = +0.45, *P* = 0.021) and UKB hyperthyroidism (β = +0.30, *P* = 0.012); the TED-enriched estimate was directionally consistent but non-significant (β = +0.34, *P* = 0.18). Sensitivity estimators were broadly concordant in direction, and the MR-Egger intercept did not detect directional pleiotropy, although this test has limited power with three or four instruments. Statistical significance was not uniform across estimators (Table S6). The magnitude of this effect lay below the gene-specific detection threshold of the discovery screen (minimum effect detectable with 80% power at the study-wide threshold, |β| = 0.99 for *IGF1R*; Table S8), so its failure to reach study-wide significance is expected for an effect of this size and does not by itself constitute evidence against a genetic contribution. This is one reason *IGF1R* was carried forward as a backbone gene evaluated across multiple evidence layers rather than judged on discovery significance alone.

Critically, this nominal MR signal was not accompanied by expression colocalization in any outcome. At the *IGF1R* locus the posterior remained below the shared-variant threshold — PP.H4 = 0.073 (discovery), 0.404 (replication) and 0.032 (TED-enriched) — with an eQTL-only pattern dominating in the two Graves disease phenotypes (PP.H2 = 0.69 and 0.62; Figure 3A). In UKB, the posterior was divided between the eQTL-only and shared-variant hypotheses (PP.H2 = 0.400; PP.H4 = 0.404), so an eQTL-only interpretation does not apply uniformly. Descriptive orbital *IGF1R* abundance had a log2 fold-change of +0.80; no differential-expression inference is made.

Thus, *IGF1R* showed nominal positive-direction MR evidence without reaching the prespecified eQTL–GWAS colocalization threshold. This profile differs from the colocalized *TSHR* signal in the two Graves disease phenotypes. The interpretation of this profile—in particular its relationship to the established therapeutic role of IGF-1R blockade—is considered in the Discussion.

### CTLA4: a positive-control autoimmune locus

*CTLA4* (chromosome 2q33), a well-established autoimmune susceptibility locus [24], was included as a positive control to calibrate the framework against a gene of known relevance to Graves disease. Consistent with its role as a known autoimmune locus, *CTLA4* showed strong protective-direction MR estimates in the discovery analysis (BBJ Graves disease, β = −1.74, *P* = 5.5×10⁻¹⁵), with concordant protective-direction estimates in replication (UKB hyperthyroidism, β = −1.57, *P* = 0.002) and the TED-enriched outcome (β = −1.77, *P* = 0.012). The recovery of *CTLA4* shows that the discovery screen detected this known autoimmune signal.

Colocalization at *CTLA4* differed across outcomes: the posterior supported a shared causal variant in both European outcomes (UKB PP.H4 = 0.953, lead rs3087243; TED-enriched PP.H4 = 0.978, lead rs1863800) but favoured distinct causal variants in the East Asian discovery outcome (PP.H3 = 0.799; Figure 3A). Together with between-instrument heterogeneity in the replication analysis (Table S6), this pattern is consistent with cross-ancestry locus complexity at *CTLA4*. This locus illustrates detection of known autoimmune biology, but it does not pass the required discovery-plus-TED colocalization criterion and therefore cannot validate the sensitivity of that combined filter. Ancestry, phenotype and locus complexity cannot be separated by these three comparisons. Descriptive orbital *CTLA4* abundance had a log2 fold-change of +1.58; no differential-expression inference is made.

---

## Discussion

This study compared genetically proxied expression and disease colocalization across Graves disease (GD), broader hyperthyroidism and a TED-enriched Graves ophthalmopathy outcome. *TSHR* showed consistent protective-direction MR estimates and strong expression colocalization in the two Graves disease phenotypes at the default prior. *IGF1R* showed nominal MR associations in BBJ and UKB, but its colocalization posteriors remained below the prespecified threshold in all three outcomes. No novel candidate met the combined discovery-plus-TED colocalization criterion. The contribution is a within-disease comparison of these genetic profiles; the analysis does not establish a pharmacologic mechanism or prove that a non-colocalized gene lacks a genetic role.

*TSHR* is an established GD autoantigen and susceptibility locus [2]. Its colocalized signal, with rs179252 as the top-posterior variant in BBJ and FinnGen, provides a reference for interpreting other genes. The UKB result bounds that interpretation: despite a strong locus-wide disease association, the posterior favoured distinct variants (PP.H3 = 0.774). A weak outcome association alone therefore does not explain the discordance, although phenotype differences, LD structure and the single-causal-variant assumption remain possible contributors. Colocalization support is specific to the two Graves disease phenotypes evaluated here, and the posterior does not establish rs179252 as the uniquely causal variant.

The protective-direction blood-expression MR estimate and descriptively higher orbital *TSHR* abundance concern different biological settings. The former indexes inherited blood regulatory variation; the latter measures disease-state tissue from four affected individuals and one control. Tissue-specific regulation [25], cell composition, inflammation and remodeling can produce different directions across these settings. The tissue observation cannot establish a causal direction, confirm the MR result or resolve the mechanism linking the locus to disease.

For *IGF1R*, the nominal MR associations were not accompanied by PP.H4 ≥ 0.80. BBJ and FinnGen favoured the eQTL-only hypothesis, whereas UKB split support between eQTL-only and shared-variant hypotheses (PP.H2 = 0.400; PP.H4 = 0.404). This is an unresolved shared-variant architecture under the tested model, not evidence that disease association is absent. Blood expression may also be an incomplete proxy for tissue-specific or protein-level receptor biology. The results are compatible with the established therapeutic efficacy of IGF-1R blockade [3] and receptor crosstalk with TSHR [4, 5], but neither demonstrate that *IGF1R* is exclusively an effector nor exclude a modest inherited contribution.

Genetic susceptibility and therapeutic efficacy are related but distinct forms of evidence [26]. That general principle is established; our contribution is the comparison of *TSHR* and *IGF1R* using gene-specific instruments and a common analytic framework. A colocalized susceptibility signal does not show that manipulating expression will reproduce the effect of a drug, identify the appropriate tissue or treatment direction, or predict superiority over an established therapy.

The absence of a qualifying novel candidate is specific to the filters and power of this screen. *TNFSF14* and *IFNGR1* colocalized in BBJ but not in either secondary outcome; *TNFSF14* nevertheless retained a direction-consistent UKB MR association. The chromosome 16p11.2 hits were consistent with one regional signal under the East Asian LD and conditional analysis. These results caution against promoting discovery significance alone to a target claim. Power varied substantially across genes: only about 36% could detect an OR of 2.0 at 80% power in discovery, and the median detectable OR was about 2.55. These thresholds concern MR association, not the probability of passing colocalization or the combined filter. Failure to pass the complete filter therefore does not rule out other large or moderate druggable effects.

Several limitations determine the scope of inference. First, the exposure data were predominantly European blood *cis*-eQTLs [12], whereas BBJ discovery was East Asian [13]. Tissue context [25], LD differences and variant availability may affect instrument performance and colocalization. Participant overlap between eQTLGen and the European outcomes could not be quantified from the available summary data; overlapping samples can introduce bias, and its magnitude and direction were not assessed here. Second, exposure effect sizes were reconstructed using 1000 Genomes European frequencies rather than eQTLGen's own frequency file. The reconstruction rescales exposure β and SE together; for a fixed single-instrument Wald analysis, the MR effect magnitude and SE rescale together, preserving the effect direction and *P* value. Multi-instrument weights, reconstructed effect magnitudes, MR detection thresholds and beta/varbeta-based colocalization posteriors can change. A sensitivity analysis using eQTLGen's own allele frequencies has not been performed, so their stability under that substitution remains unverified.

Third, UKB hyperthyroidism is broader than GD. Fourth, *TSHR* used a single instrument selected with the European panel matching the eQTL data. The East Asian panel yields two index variants at the same primary threshold (Table S3); single-instrument status therefore follows from the selection reference and does not demonstrate a single signal in every ancestry. Single-instrument MR precludes instrument-based pleiotropy diagnostics, and colocalization cannot exclude effects through a neighbouring gene. MR-Egger also has limited power at the other backbone loci because instrument counts are small. Fifth, the FinnGen outcome includes only 858 cases and has limited locus-specific power; a low PP.H4 can reflect weak outcome evidence rather than distinct causal architecture.

Sixth, FinnGen compares Graves ophthalmopathy cases with population controls, not patients with GD without orbitopathy. Its associations therefore include GD susceptibility and do not isolate liability to orbital involvement. Agreement with BBJ should not be interpreted as a TED-specific causal effect. Seventh, coloc.abf assumes at most one causal variant per trait in the tested region, without a multi-variant extension or formal colocalization power analysis. Prior sensitivity was material: at p12 = 10⁻⁶, *TSHR* in BBJ and *CTLA4* in UKB fell below 0.80 (0.661 and 0.672), whereas *TSHR* in FinnGen remained above it (0.875; Table S9).

Eighth, the single positive control, *CTLA4* [24], colocalized in UKB and FinnGen but not BBJ (PP.H3 = 0.799). It demonstrates detection in individual outcomes, not the sensitivity of the combined discovery-plus-TED rule, which it fails. Ancestry and phenotype are confounded across these comparisons. Finally, the orbital dataset has only one control, no external replication and no inferential role. Its descriptive fold changes must not be treated as evidence for or against differential expression.

Within these limits, *TSHR* showed an expression-colocalized susceptibility profile in the two Graves disease phenotypes, while *IGF1R* showed nominal MR associations without meeting the shared-variant threshold. The contrast can inform interpretation of genetic target prioritization alongside established therapeutic evidence, while leaving the tissue-specific mechanism, modest genetic effects and TED-specific liability unresolved.

---

## Declarations

**Funding.** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

**Conflict of interest.** The authors declare that they have no conflict of interest.

**Ethics approval.** This study used publicly available, de-identified summary statistics and an institutionally approved orbital tissue dataset. The in-house orbital transcriptomic component was approved by the Institutional Review Board of Pusan National University Hospital (approval number 2104-018-102) and was conducted in accordance with the Declaration of Helsinki.

**Informed consent.** Written informed consent was obtained from all individual participants included in the in-house orbital tissue study.

**Data availability.** Instrument-level data, primary gene-level MR results and full colocalization posteriors are supplied as Supplementary Data 1–3. Analysis code is available from the corresponding author on reasonable request. Public summary statistics analyzed in this study are available from their original repositories: blood *cis*-eQTL data from the eQTLGen Consortium, the Biobank Japan Graves disease and UK Biobank hyperthyroidism genome-wide association statistics through the GWAS Catalog, and the FinnGen Release 12 Graves ophthalmopathy statistics from FinnGen. The in-house orbital RNA-seq data are available from the corresponding author on reasonable request, subject to institutional and ethical restrictions.

**Author contributions.** J.P. conceived and designed the study, performed the analyses, and drafted the manuscript. M.-S.K. contributed to data collection and interpretation. K.-H.S. and S.-W.Y. supervised the study and revised the manuscript. All authors read and approved the final manuscript.

**Use of AI-assisted tools.** OpenAI Codex assisted with manuscript editing, code preparation, numerical consistency checks, figure rebuilding and document formatting. Scientific interpretation and responsibility for the submitted work remain with the authors.

**Acknowledgements.** We thank the eQTLGen Consortium, Biobank Japan, the UK Biobank, and the FinnGen study and its participants for making their summary statistics publicly available. We acknowledge the GWAS Catalog for hosting and distributing the genome-wide association summary statistics used as outcome data.

---

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

17. Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, Montgomery GW, Goddard ME, Wray NR, Visscher PM, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016;48:481-487. doi:10.1038/ng.3538.

18. Chang CC, Chow CC, Tellier LC, Vattikuti S, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015;4:7. doi:10.1186/s13742-015-0047-8.

19. Hemani G, Zheng J, Elsworth B, Wade KH, Haberland V, Baird D, Laurin C, Burgess S, Bowden J, Langdon R, et al. The MR-Base platform supports systematic causal inference across the human phenome. Elife. 2018;7:e34408. doi:10.7554/elife.34408.

20. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. Int J Epidemiol. 2015;44:512-525. doi:10.1093/ije/dyv080.

21. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent Estimation in Mendelian Randomization with Some Invalid Instruments Using a Weighted Median Estimator. Genet Epidemiol. 2016;40:304-314. doi:10.1002/gepi.21965.

22. Hemani G, Tilling K, Davey Smith G. Orienting the causal relationship between imprecisely measured traits using GWAS summary data. PLoS Genet. 2017;13:e1007081. doi:10.1371/journal.pgen.1007081.

23. Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLoS Genet. 2020;16:e1008720. doi:10.1371/journal.pgen.1008720.

24. Ueda H, Howson JM, Esposito L, Heward J, Snook H, Chamberlain G, Rainbow DB, Hunter KM, Smith AN, Di Genova G, et al. Association of the T-cell regulatory gene CTLA4 with susceptibility to autoimmune disease. Nature. 2003;423:506-511. doi:10.1038/nature01621.

25. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369:1318-1330. doi:10.1126/science.aaz1776.

26. Nelson MR, Tipney H, Painter JL, Shen J, Nicoletti P, Shen Y, Floratos A, Sham PC, Li MJ, Wang J, et al. The support of human genetic evidence for approved drug indications. Nat Genet. 2015;47:856-860. doi:10.1038/ng.3314.

---

## Figure Legends

**Figure 1. Study flow and the *TSHR*–*IGF1R* comparison.** **(a)** Attrition from the druggable genome [16] to the thirteen genes reaching the study-wide Bonferroni threshold (*P* < 1.965 × 10⁻⁵), with exclusions at each step and the resolution of the thirteen at right. No novel candidate met the colocalization criterion. **(b)** Discovery MR estimates and colocalization posteriors for the two backbone receptors. *TSHR* shows shared-variant support in both Graves disease phenotypes; *IGF1R* does not reach PP.H4 ≥ 0.80 in any outcome.

**Figure 2. Discovery Mendelian randomization across the druggable genome.** **(A)** All estimable genes. **(B)** Enlarged effect range containing all thirteen discovery hits. *Cis*-eQTL MR effects on Biobank Japan Graves disease (log-odds) for the 2,234 genes with estimable estimates, among 2,544 MR-testable genes. The dashed line marks the Bonferroni threshold (*P* = 1.965 × 10⁻⁵) computed against all 2,544. Thirteen genes exceed it: *TSHR* and *CTLA4* (red), five MHC-region signals (grey), three chromosome 16p11.2 variants forming one conditional signal (light purple; Table S5), and three further candidates. *IGF1R* (blue) lies below the threshold.

**Figure 3. Colocalization and MR estimates for the backbone genes.** **(A)** Shared-variant posterior probabilities (PP.H4) against each outcome, computed from outcome beta/varbeta and exposure beta/varbeta reconstructed using European frequencies (sdY = 1). Numerical values are printed beside the points and the dashed line marks 0.80. *TSHR* colocalizes in both Graves disease phenotypes but not in UKB hyperthyroidism, where distinct variants are favoured. *IGF1R* reaches the threshold in no outcome; *CTLA4* colocalizes in both European outcomes. Default prior p12 = 10⁻⁵; sensitivity in Table S9. **(B)** MR estimates (β, log-odds) with 95% confidence intervals and *P* values.

---

## Tables

**Table 1. Study design and data sources.**

| Dataset | Accession / endpoint | Role | Ancestry | Cases | Controls | Total *N* | Case fraction |
|---|---|---|---|---|---|---|---|
| eQTLGen blood *cis*-eQTL [12] | — | Exposure (instruments) | Predominantly European | — | — | up to 31,684 | — |
| Biobank Japan Graves disease [13] | GCST90018627 | Primary discovery outcome | East Asian | 2,809 | 172,656 | 175,465 | 0.01601 |
| UK Biobank hyperthyroidism [14] | GCST90038636 | Cross-ancestry replication outcome | European | 3,731 | 480,867 | 484,598 | 0.00770 |
| FinnGen R12 Graves ophthalmopathy [15] | GRAVES_OPHT (R12) | TED-enriched sensitivity outcome | European | 858 | 499,490 | 500,348 | 0.00171 |
| In-house orbital RNA-seq | — | Descriptive tissue observation | Korean | 4 TED | 1 control | 5 | — |

eQTLGen sample size indicates the maximum available blood cis-eQTL sample size; SNP-level sample sizes vary across instruments (e.g., rs179252 NrSamples = 31,566). BBJ Graves disease was used as the primary discovery phenotype, UKB hyperthyroidism as cross-ancestry/broader-phenotype replication, and FinnGen Graves ophthalmopathy as the TED-enriched sensitivity outcome. UKB binary-trait effect estimates were rescaled from a linear mixed-model scale to log-odds using case prevalence.

**Table 2. Mendelian randomization and colocalization evidence for the backbone genes.**

| Gene | Outcome | N IV | OR (95% CI) | IVW/Wald *P* | Weighted median *P* | MR-Egger intercept *P* | Cochran's Q *P* | Coloc PP.H4 | Coloc PP.H3 | Coloc lead SNP |
|---|---|---|---|---|---|---|---|---|---|---|
| ***TSHR*** | BBJ Graves | 1 | 0.12 (0.07–0.21) | 1.1×10⁻¹⁴ | NA (single IV) | NA (single IV) | NA (single IV) | 0.951 | 0.049 | rs179252 |
| | UKB hyperthyroid | 1 | 0.09 (0.06–0.14) | 8.8×10⁻²⁸ | NA | NA | NA | 0.226 | 0.774 | rs1023586 |
| | FinnGen GO | 1 | 0.10 (0.04–0.24) | 2.8×10⁻⁷ | NA | NA | NA | 0.986 | 0.014 | rs179252 |
| ***IGF1R*** | BBJ Graves | 4 | 1.56 (1.07–2.28) | 0.021 | 0.029 | 0.663 | 0.902 | 0.073 | 0.236 | rs2654980 |
| | UKB hyperthyroid | 4 | 1.35 (1.07–1.70) | 0.012 | 0.002 | 0.562 | 0.250 | 0.404 | 0.196 | rs2654980 |
| | FinnGen GO | 3 | 1.41 (0.85–2.33) | 0.182 | 0.124 | 0.433 | 0.309 | 0.032 | 0.346 | rs2654980 |
| ***CTLA4*** | BBJ Graves | 1 | 0.18 (0.11–0.27) | 5.5×10⁻¹⁵ | NA (single IV) | NA (single IV) | NA (single IV) | 0.201 | 0.799 | rs231811 |
| | UKB hyperthyroid | 2 | 0.21 (0.08–0.57) | 0.002 | — | NA (2 IV) | 0.001 | 0.953 | 0.047 | rs3087243 |
| | FinnGen GO | 2 | 0.17 (0.04–0.68) | 0.012 | — | NA (2 IV) | 0.049 | 0.978 | 0.022 | rs1863800 |

MR effect estimates are odds ratios (OR) per unit increase in genetically proxied gene expression, by Wald ratio (single instrument) or inverse-variance weighted (IVW; multiple instruments). Sensitivity estimators (weighted median, weighted mode, MR-Egger intercept) and the Cochran's Q heterogeneity test are reported only where estimable; single-instrument loci are marked NA. MR-Egger intercept tests for IGF1R were non-significant (all P > 0.43), with limited power given the instrument counts; Cochran's Q indicates significant between-instrument heterogeneity for CTLA4 in the UKB and FinnGen analyses (Q P = 0.001 and 0.049), consistent with cross-ancestry/phenotype locus complexity. Colocalization (coloc.abf) used each outcome's beta and varbeta and exposure beta/varbeta reconstructed using European reference frequencies (sdY = 1), with dataset 1 = outcome and dataset 2 = eQTL; PP.H4 ≥ 0.80 indicates strong evidence for a shared causal variant and PP.H3 evidence for distinct causal variants. IGF1R is PP.H2-dominant in both Graves disease phenotypes (PP.H2 = 0.69 and 0.62), i.e. a cis-eQTL signal that does not resolve to a variant shared with the outcome; the five-hypothesis colocalization comparison does not negate the nominal MR association reported above. UKB instead splits support between H2 (0.400) and H4 (0.404). Posteriors shown are at the default prior p12 = 10⁻⁵; sensitivity at p12 = 10⁻⁶ and 5×10⁻⁶ is reported in Table S9 (backbone) and Table S2 (candidates). —, not applicable/not computed. All values verified against source analysis files.

**Table 3. All thirteen Bonferroni-significant druggable-gene-wide BBJ discovery hits, with cross-outcome colocalization and classification.**

| Gene | N IV | OR (95% CI) | BBJ *P* | Coloc PP.H4 (BBJ) | Coloc PP.H4 (UKB) | Coloc PP.H4 (FinnGen) | Classification |
|---|---|---|---|---|---|---|---|
| *HLA-A* | 4 | 1.98 (1.73–2.27) | 2.6×10⁻²³ | — | — | — | MHC region |
| *HLA-DQA2* | 1 | 2.39 (1.95–2.93) | 8.3×10⁻¹⁷ | — | — | — | MHC region |
| ***CTLA4*** | 1 | 0.18 (0.11–0.27) | 5.5×10⁻¹⁵ | 0.201 | 0.953 | 0.978 | Established autoimmune locus (positive control) |
| ***TSHR*** | 1 | 0.12 (0.07–0.21) | 1.1×10⁻¹⁴ | 0.951 | 0.226 | 0.986 | Established GD locus (susceptibility anchor) |
| *C4A* | 1 | 0.45 (0.35–0.57) | 2.2×10⁻¹⁰ | — | — | — | MHC region |
| *HSD3B7* | 1 | 0.29 (0.18–0.47) | 2.0×10⁻⁷ | 0.636 | 0.095 | 0.026 | chr16p11.2 LD cluster |
| *TUBB* | 2 | 1.49 (1.28–1.74) | 3.1×10⁻⁷ | — | — | — | MHC-region-linked |
| *VKORC1* | 1 | 0.13 (0.06–0.29) | 6.4×10⁻⁷ | 0.388 | 0.030 | 0.030 | chr16p11.2 LD cluster |
| *TNFSF14* | 1 | 0.63 (0.52–0.76) | 1.5×10⁻⁶ | 0.994 | 0.237 | 0.017 | Candidate (single-outcome coloc) |
| *PRSS36* | 1 | 6.00 (2.87–12.57) | 2.0×10⁻⁶ | 0.168 | 0.115 | 0.037 | chr16p11.2 LD cluster |
| *MAPKAPK5* | 4 | 0.98 (0.97–0.99) | 5.3×10⁻⁶ | <0.001 | 0.005 | 0.026 | Candidate (distinct-variant signal) |
| *PSMB8* | 1 | 1.25 (1.13–1.38) | 6.9×10⁻⁶ | — | — | — | MHC region |
| *IFNGR1* | 1 | 2.10 (1.51–2.91) | 9.4×10⁻⁶ | 0.989 | 0.052 | 0.020 | Candidate (single-outcome coloc) |

Primary discovery estimates against BBJ Graves disease (Wald ratio for single-instrument loci, IVW for multi-instrument loci), ordered by discovery P value. OR per unit increase in genetically proxied expression. Colocalization (coloc.abf) was evaluated against all three outcomes; "—" denotes loci in the MHC region or otherwise not carried into formal colocalization. The two established loci (TSHR, CTLA4) and MHC-region signals are classified by prior biological knowledge. Among non-known candidates, TNFSF14 and IFNGR1 colocalized strongly in BBJ (PP.H4 = 0.994 and 0.989) but in neither secondary outcome, and the chr16p11.2 cluster resolves to a single conditional signal (Supplementary Table S5); no non-known candidate showed shared-variant colocalization in more than one outcome.

---

## Supplementary Material

**Supplementary Tables**





**Table S1.** Extended cross-outcome detail for the thirteen Bonferroni-significant BBJ discovery hits. Whereas Table 3 summarizes the discovery (BBJ) effect sizes and colocalization, this table provides the per-hit Mendelian randomization estimates across all three outcomes (BBJ Graves disease, UKB hyperthyroidism, FinnGen Graves ophthalmopathy) together with the colocalization lead variant, supporting the cross-outcome reproducibility assessment.

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
| *IFNGR1* | +0.74 (9.4×10⁻⁶) | +0.30 (0.070) | +0.15 (0.68) | rs11754268 | Candidate (single-outcome coloc) |

  Effect estimates (β, log-odds) by Wald ratio (single instrument) or inverse-variance weighted (multiple instruments). "—" for FinnGen *VKORC1* indicates no estimable instrument in that outcome; "—" in the colocalization column denotes MHC-region loci not carried into formal colocalization. Several discovery hits do not reproduce in the broader-phenotype (UKB) or TED-enriched (FinnGen) outcomes. Among non-MHC established anchors, *TSHR* and *CTLA4* showed consistent protective-direction estimates across all three outcomes; MHC-region signals were interpreted separately because of complex regional linkage disequilibrium.




**Table S2.** Candidate colocalization across all three outcomes for non-known candidate hits. Posterior probabilities are from coloc.abf using outcome beta/varbeta and exposure beta/varbeta reconstructed from European reference frequencies (sdY = 1) (H0 no association; H1 outcome only; H2 *cis*-eQTL only; H3 distinct causal variants; H4 shared causal variant). The last two columns give PP.H4 under more conservative priors.

| Gene | Outcome | Overlapping SNPs | PP.H0 | PP.H1 | PP.H2 | PP.H3 | PP.H4 | Top SNP | PP.H4 (p12=10⁻⁶) | PP.H4 (p12=5×10⁻⁶) |
|---|---|---|---|---|---|---|---|---|---|---|
| *TNFSF14* | BBJ Graves | 4,790 | 0.000 | 0.000 | 0.001 | 0.005 | 0.994 | rs2291668 | 0.941 | 0.988 |
| *TNFSF14* | UKB hyperthyroid | 7,197 | 0.000 | 0.000 | 0.360 | 0.403 | 0.237 | rs2291668 | 0.030 | 0.134 |
| *TNFSF14* | FinnGen GO | 7,703 | 0.000 | 0.000 | 0.626 | 0.357 | 0.017 | rs2291668 | 0.002 | 0.009 |
| *IFNGR1* | BBJ Graves | 3,973 | 0.000 | 0.000 | 0.005 | 0.006 | 0.989 | rs11754268 | 0.899 | 0.978 |
| *IFNGR1* | UKB hyperthyroid | 5,145 | 0.000 | 0.000 | 0.736 | 0.211 | 0.052 | rs11754268 | 0.005 | 0.027 |
| *IFNGR1* | FinnGen GO | 5,891 | 0.000 | 0.000 | 0.641 | 0.339 | 0.020 | rs11754268 | 0.002 | 0.010 |
| *MAPKAPK5* | BBJ Graves | 2,077 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | rs79271898 | 0.000 | 0.000 |
| *MAPKAPK5* | UKB hyperthyroid | 2,763 | 0.000 | 0.000 | 0.052 | 0.944 | 0.005 | rs79271898 | 0.000 | 0.002 |
| *MAPKAPK5* | FinnGen GO | 3,294 | 0.000 | 0.000 | 0.663 | 0.311 | 0.026 | rs79271898 | 0.003 | 0.013 |
| *HSD3B7* | BBJ Graves | 1,630 | 0.000 | 0.000 | 0.000 | 0.364 | 0.636 | rs4889606 | 0.149 | 0.467 |
| *HSD3B7* | UKB hyperthyroid | 2,755 | 0.000 | 0.000 | 0.444 | 0.461 | 0.095 | rs4889606 | 0.010 | 0.050 |
| *HSD3B7* | FinnGen GO | 3,104 | 0.000 | 0.000 | 0.733 | 0.240 | 0.026 | rs4889606 | 0.003 | 0.013 |
| *VKORC1* | BBJ Graves | 1,516 | 0.000 | 0.000 | 0.000 | 0.612 | 0.388 | rs34649473 | 0.060 | 0.241 |
| *VKORC1* | UKB hyperthyroid | 2,644 | 0.000 | 0.000 | 0.477 | 0.493 | 0.030 | rs34649473 | 0.003 | 0.015 |
| *VKORC1* | FinnGen GO | 2,961 | 0.000 | 0.000 | 0.739 | 0.231 | 0.030 | rs2884737 | 0.003 | 0.015 |
| *PRSS36* | BBJ Graves | 1,479 | 0.000 | 0.000 | 0.000 | 0.832 | 0.168 | rs78924645 | 0.020 | 0.092 |
| *PRSS36* | UKB hyperthyroid | 2,577 | 0.000 | 0.000 | 0.436 | 0.449 | 0.115 | rs78924645 | 0.013 | 0.061 |
| *PRSS36* | FinnGen GO | 2,821 | 0.000 | 0.000 | 0.744 | 0.219 | 0.037 | rs78924645 | 0.004 | 0.019 |

  No non-known candidate reached shared-variant colocalization (PP.H4 ≥ 0.80) in more than one outcome. TNFSF14 and IFNGR1 colocalized strongly in the BBJ discovery outcome (PP.H4 = 0.994 and 0.989) but not in either replication outcome, and the chromosome 16p11.2 candidates showed only weak or ambiguous colocalization throughout. TNFSF14 retains a direction-consistent nominal MR association in UKB (P = 0.0056), so its exclusion reflects the colocalization criterion rather than absence of association; locus-specific power and phenotype heterogeneity across outcomes preclude a stronger interpretation. Default priors: p1 = p2 = 1×10⁻⁴, p12 = 1×10⁻⁵.




**Table S3.** *TSHR* locus LD-clumping sensitivity. Independent-instrument counts at the *TSHR* locus under varying clumping thresholds, using ancestry-matched 1000 Genomes Phase 3 reference panels (European, n = 503; East Asian, n = 504). rs179252 (the primary instrument) is retained among the index variants under every threshold.

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
| 5×10⁻⁶ | 0.01 | EUR | 3 | Yes |
| 5×10⁻⁶ | 0.01 | EAS | 7 | Yes |
| 5×10⁻⁶ | 0.1 | EUR | 12 | Yes |
| 5×10⁻⁶ | 0.1 | EAS | 12 | Yes |

  rs179252 (chr14:81,435,985 hg19; marginal cis-eQTL P = 2.86×10⁻⁴⁰) is among the index variants retained under all twelve threshold combinations tested; it is the sole index variant only at the primary threshold in the European panel. At the primary instrument-selection threshold (P < 5×10⁻⁸, r² < 0.001) the European reference yields a single independent instrument; further index variants appear there only as r² is relaxed. The East Asian reference yields two at the same threshold — rs179252 and rs72690955 (chr14:81,160,405; marginal P = 1.04×10⁻¹⁷; pairwise r² with rs179252 = 0.0002) — so the single-instrument treatment of TSHR follows from the European reference used for instrument selection, which matches the predominantly European eQTL panel, and not from an absence of independent signal in every ancestry. This is stated as a limitation. Because r² is invariant to allele coding, the assessment is unaffected by effect-allele orientation.





**Table S4.** Instrument selection and strength for the backbone genes and all thirteen discovery hits. Instruments were genome-wide significant *cis*-eQTLs (P < 5×10⁻⁸); a secondary-threshold (P < 5×10⁻⁶) variant count, before clumping, is shown for reference only. Across all 2,544 MR-testable genes (6,135 instruments), the minimum F-statistic was 29.7, consistent with the *P* < 5×10⁻⁸ selection threshold (which implies |*Z*| > 5.45 and therefore *F* = *Z*² > 29.7); no weak instruments were present.

| Gene | Chr | Instruments (P < 5×10⁻⁸, clumped) | Variants (P < 5×10⁻⁶, unclumped) | Single-instrument | Druggability tier |
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
| *IFNGR1* | 6 | 2 | 340 | No | Tier 1 |

  Druggability tiers follow Finan et al. (Tier 1, approved drug targets or those in clinical development; Tier 2, targets of small molecules or biotherapeutics with druggable similarity; Tier 3, secreted or extracellular proteins and members of druggable gene families). The column labelled "Variants (P < 5×10⁻⁶, unclumped)" counts variants passing the relaxed threshold before clumping and is therefore not an instrument count; it is shown only to indicate signal density, which is why MHC-region genes reach five figures. Only the clumped-instrument column contributed instruments to the primary analyses; the number of instruments actually used in MR against a given outcome (N IV in Tables 2 and 3) may be smaller, reflecting outcome-specific SNP availability and harmonization (absence of a variant in the outcome GWAS, or removal of palindromic variants).




**Table S5.** Chromosome 16p11.2 candidate cluster. Three discovery hits (*HSD3B7*, *VKORC1*, *PRSS36*) map to a narrow (~143-kb) window on chromosome 16p11.2 and behave as a single linkage-disequilibrium signal rather than three independent associations.

| Gene | Chr | Position (hg19) | BBJ discovery *P* | Coloc PP.H4 (BBJ) | Coloc PP.H4 (UKB) | Coloc PP.H4 (FinnGen) | Top SNP |
|---|---|---|---|---|---|---|---|
| *HSD3B7* | 16 | 31,011,183 | 2.0×10⁻⁷ | 0.636 | 0.095 | 0.026 | rs4889606 |
| *VKORC1* | 16 | 31,066,538 | 6.4×10⁻⁷ | 0.388 | 0.030 | 0.030 | rs34649473 |
| *PRSS36* | 16 | 31,154,358 | 2.0×10⁻⁶ | 0.168 | 0.115 | 0.037 | rs78924645 |

  All three lead variants lie within a 143.2-kb window on chromosome 16p11.2 (hg19 chr16:31,011,183–31,154,358) and none reaches strong colocalization (PP.H4 ≥ 0.80) in any outcome. Because these hits were discovered in the East Asian outcome, linkage disequilibrium was assessed in the ancestry-matched 1000 Genomes Phase 3 East Asian reference (n = 504), where all three lead variants are in substantial mutual LD (r² = 0.854, 0.761 and 0.863). Conditional analysis (GCTA-COJO, East Asian reference) selected exactly one independent signal at the locus, rs8050588 (P = 1.15×10⁻⁸); conditioning on it collapses all three associations (rs4889606, marginal P = 2.04×10⁻⁷ → conditional P = 0.928; rs34649473, 6.35×10⁻⁷ → 0.839; rs78924645, 2.02×10⁻⁶ → 0.806). Their co-significance in discovery therefore reflects a single regional signal rather than three independent causal effects. In the European reference the corresponding r² values are 0.896, 0.190 and 0.204, which is why the ancestry-matched panel is the appropriate one here.




**Table S6.** Mendelian randomization sensitivity analyses for the backbone genes across all three outcomes. IVW, inverse-variance weighted; WM, weighted median; Wmode, weighted mode. MR-Egger intercept, Cochran's Q, and multi-estimator comparisons require multiple instruments and are not estimable (NA) at single-instrument loci.

| Gene | Outcome | N IV | IVW/Wald *P* | WM *P* | Wmode *P* | MR-Egger intercept *P* | Cochran's Q *P* |
|---|---|---|---|---|---|---|---|
| *TSHR* | BBJ Graves | 1 | 1.1×10⁻¹⁴ | NA | NA | NA | NA |
| *TSHR* | UKB hyperthyroid | 1 | 8.8×10⁻²⁸ | NA | NA | NA | NA |
| *TSHR* | FinnGen GO | 1 | 2.8×10⁻⁷ | NA | NA | NA | NA |
| *IGF1R* | BBJ Graves | 4 | 0.021 | 0.029 | 0.512 | 0.663 | 0.902 |
| *IGF1R* | UKB hyperthyroid | 4 | 0.012 | 0.002 | 0.367 | 0.562 | 0.250 |
| *IGF1R* | FinnGen GO | 3 | 0.182 | 0.124 | 0.633 | 0.433 | 0.309 |
| *CTLA4* | BBJ Graves | 1 | 5.5×10⁻¹⁵ | NA | NA | NA | NA |
| *CTLA4* | UKB hyperthyroid | 2 | 0.002 | NA | NA | NA | 0.001 |
| *CTLA4* | FinnGen GO | 2 | 0.012 | NA | NA | NA | 0.049 |

  For IGF1R (the only multi-instrument backbone gene with a full estimator set), all sensitivity estimators are directionally concordant with the IVW estimate and the MR-Egger intercept is non-significant in every outcome (P > 0.43), with limited power to assess directional pleiotropy given the small instrument counts. CTLA4 showed significant between-instrument heterogeneity in the UKB and FinnGen analyses (Cochran's Q P = 0.001 and 0.049), consistent with cross-ancestry/phenotype locus complexity. Steiger directionality testing was not applied; no formal directionality result is available. Single-instrument loci are not amenable to instrument-based pleiotropy diagnostics. TSHR is instead supported by colocalization (PP.H4 = 0.951 in BBJ and 0.986 in the TED-enriched outcome) and by LD-clumping sensitivity analysis (Table S3); the single-instrument CTLA4 discovery arm is not, its discovery posterior favouring distinct causal variants (PP.H3 = 0.799, PP.H4 = 0.201), consistent with its role as a positive control rather than a colocalized anchor.




**Table S7. STROBE-MR reporting map.** Item numbering follows the official 20-item checklist; the descriptions below are adapted and explicitly identify unavailable analyses or information.

| Item | Reporting topic | Location and coverage |
|---|---|---|
| 1 | MR design | Abstract explicitly names two-sample Mendelian randomization. |
| 2 | Rationale | Introduction: blood expression, disease susceptibility and therapeutic interpretation. |
| 3 | Objectives | Introduction final paragraph; Methods defines the hierarchy and biological backbone. No prospective registration identifier is available. |
| 4a–b | Design, participants and sample size | Methods; Table 1; source reports [12–15]. Recruitment and eligibility details are those of the original studies; available sample sizes determined the analysis. Detection thresholds were calculated from observed MR standard errors, not for prospective recruitment. |
| 4c–d | Variants and phenotype definitions | Methods instrument selection and harmonization; Table 1 accessions; Table S4. Source studies provide their genotyping and phenotype procedures. |
| 4e | Ethics and consent | Declarations; in-house IRB approval and written consent. |
| 5 | Instrument assumptions | Methods states relevance, independence and exclusion restriction; Discussion discusses limits of verification. |
| 6a–b | Scales and variant weights | Methods: Z-score reconstruction, UKB log-odds rescaling, Wald/IVW and inverse-variance weighting. |
| 6c–e | Estimation, covariates, missingness and multiplicity | Methods: source GWAS adjustments retained; no individual-level re-adjustment; unavailable or excluded variants omitted; Bonferroni denominator 2,544. |
| 7 | Assumption assessment | Methods: instrument strength, harmonization, heterogeneity, MR-Egger and colocalization. These do not prove all instrument assumptions. |
| 8 | Additional analyses | Methods; Tables S2–S3, S5–S6, S8–S9. eQTLGen frequency substitution was not performed. |
| 9a–b | Software and registration | Methods lists R, TwoSampleMR, coloc and PLINK versions. No prospective registration identifier is available. |
| 10a–c | Descriptive data | Table 1 gives source sample counts; Figure 1 gives gene attrition; Table S4 gives instrument counts. Participant-level distributions and cohort-specific eQTL meta-analysis heterogeneity were not reanalysed. |
| 10d | Transportability and participant overlap | Discussion: European blood eQTL transfer to East Asian discovery is an assumption; overlap with European outcomes is unquantified. |
| 11a–d | Main results and uncertainty | Tables 2–3 and S1; Figures 2–3. ORs/CIs are per reconstructed expression unit. Variant-level instruments and primary MR results are supplied as Supplementary Data 1–2. Absolute risks were not estimated. |
| 12a–b | Results of assumption checks | Instrument strength in Table S4; heterogeneity and MR-Egger in Tables 2 and S6; clumping and colocalization sensitivity in Tables S2–S3 and S9. |
| 13a–e | Additional results | Tables S2–S3 and S5–S9; Figures S1–S2. Steiger and bidirectional MR were not applied; tissue data are descriptive; no leave-one-out result is reported. |
| 14 | Key findings | Discussion opening paragraph. |
| 15 | Limitations and bias | Discussion: tissue, ancestry, overlap, frequency reconstruction, phenotype, instrument counts, power, coloc model/priors and one tissue control; bias magnitudes were not quantified. |
| 16a–c | Interpretation and clinical meaning | Discussion: genetic expression effects do not establish a pharmacologic mechanism, treatment direction or intervention effect size. |
| 17 | Generalizability | Discussion: ancestry and tissue context, broad hyperthyroidism and population controls; inherited expression proxies do not estimate acute or dose-specific treatment effects. |
| 18 | Funding | Declarations states no specific grant; original data-source funding is reported in the cited source publications. |
| 19 | Data and code access | Declarations lists source repositories and access restrictions. Supplementary Data 1–3 provide instruments, MR and coloc summaries. Analysis code is available from the corresponding author on request. |
| 20 | Competing interests | Declarations. |

Adapted from the STROBE-MR checklist (EQUATOR Network, CC BY 3.0), https://www.strobe-mr.org/download/strobe-mr-checklist/. The checklist is a reporting map, not a claim that all study limitations have been resolved.





**Table S8.** Detectable effect range of the druggable-gene-wide screen. Minimum effect detectable with 80% power was computed per gene and outcome from the empirical standard error of the primary MR estimate as |β|min = (*z*1−α/2 + *z*0.80) × SE, with α = 1.965×10⁻⁵ (discovery) or α = 0.05 (replication, TED-enriched).

  **Panel A. Screen-wide detectable range.**

| Outcome | Significance threshold | Genes | Median absolute βmin (IQR) | Median detectable OR (IQR) | Genes with ≥80% power at OR 1.5 | At OR 2.0 | At OR 3.0 |
|---|---|---|---|---|---|---|---|
| BBJ Graves (discovery) | *P* < 1.965×10⁻⁵ | 2,234 | 0.94 (0.54–1.64) | 2.55 (1.71–5.18) | 14.6% | 35.6% | 57.7% |
| UKB hyperthyroid (replication) | *P* < 0.05 | 2,505 | 0.38 (0.23–0.65) | 1.46 (1.26–1.91) | 53.7% | 77.4% | 93.1% |
| FinnGen Graves ophthalmopathy (TED-enriched) | *P* < 0.05 | 2,480 | 0.75 (0.46–1.34) | 2.12 (1.58–3.81) | 19.4% | 45.5% | 67.3% |

  **Panel B. Backbone genes: observed effect versus gene-specific detection threshold.**

| Gene | Outcome | Observed β | SE | Absolute βmin (80% power) | Detectable OR | Observed magnitude ≥ threshold? |
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

  Detection thresholds vary across genes. TSHR and CTLA4 exceeded their gene-specific discovery thresholds, but the CTLA4 FinnGen estimate was below its 80% power threshold. IGF1R effects fell below the corresponding thresholds in every outcome, so the absence of study-wide significance at IGF1R is expected for an effect of this size and is not evidence against a genetic contribution; IGF1R was therefore evaluated as a backbone gene across multiple evidence layers rather than on discovery significance alone. These calculations describe MR association power and do not quantify power to pass colocalization or the combined filter. Values use the primary MR standard errors restricted to the 2,544 genes in the verified instrument manifest and a discovery threshold of 0.05/2,544.






**Table S9.** Full colocalization posteriors and prior sensitivity for the backbone genes, laid out as Table S2. Two rows fall below the PP.H4 ≥ 0.80 threshold under the most conservative prior: *TSHR* in the discovery outcome (0.951 → 0.661) and *CTLA4* in the replication outcome (0.953 → 0.672). The TED-enriched *TSHR* posterior stays above the threshold throughout (0.986 → 0.875).

| Gene | Outcome | Overlapping SNPs | PP.H0 | PP.H1 | PP.H2 | PP.H3 | PP.H4 | Top SNP | PP.H4 (p12=10⁻⁶) | PP.H4 (p12=5×10⁻⁶) |
|---|---|---|---|---|---|---|---|---|---|---|
| *TSHR* | BBJ Graves | 4,336 | 0.000 | 0.000 | 0.000 | 0.049 | 0.951 | rs179252 | 0.661 | 0.907 |
| *TSHR* | UKB hyperthyroid | 5,667 | 0.000 | 0.000 | 0.000 | 0.774 | 0.226 | rs1023586 | 0.028 | 0.128 |
| *TSHR* | FinnGen GO | 6,579 | 0.000 | 0.000 | 0.000 | 0.014 | 0.986 | rs179252 | 0.875 | 0.972 |
| *IGF1R* | BBJ Graves | 5,685 | 0.000 | 0.000 | 0.690 | 0.236 | 0.073 | rs2654980 | 0.008 | 0.038 |
| *IGF1R* | UKB hyperthyroid | 7,394 | 0.000 | 0.000 | 0.400 | 0.196 | 0.404 | rs2654980 | 0.063 | 0.253 |
| *IGF1R* | FinnGen GO | 7,740 | 0.000 | 0.000 | 0.623 | 0.346 | 0.032 | rs2654980 | 0.003 | 0.016 |
| *CTLA4* | BBJ Graves | 3,014 | 0.000 | 0.000 | 0.000 | 0.799 | 0.201 | rs231811 | 0.025 | 0.112 |
| *CTLA4* | UKB hyperthyroid | 4,283 | 0.000 | 0.000 | 0.000 | 0.047 | 0.953 | rs3087243 | 0.672 | 0.911 |
| *CTLA4* | FinnGen GO | 4,924 | 0.000 | 0.000 | 0.000 | 0.022 | 0.978 | rs1863800 | 0.816 | 0.957 |

  Default priors p1 = p2 = 1×10⁻⁴, p12 = 1×10⁻⁵. Posteriors use each GWAS's beta and varbeta, exposure beta/varbeta reconstructed with European reference frequencies and sdY = 1, and the sample definitions in Table 1.

**Supplementary Figures**

**Figure S1. Orbital tissue transcript abundance (descriptive).** Points represent individual samples; horizontal lines mark the mean TED abundance. Log2 fold changes use a 0.01-TPM pseudocount. Transcripts per million for the backbone genes in four TED samples and one control, with technical replicates collapsed to the biological-sample level. With one control, control-group biological variance cannot be estimated; no inferential test is reported and no conclusion depends on this panel. Axes are gene-specific.

**Figure S2. Candidate colocalization across the three outcomes.** Shared-variant posterior probabilities for the six non-known candidates, plotted as in Figure 3A; full posteriors are reported in Table S2. *TNFSF14* and *IFNGR1* colocalize in discovery (PP.H4 = 0.994 and 0.989) but in neither secondary outcome; *TNFSF14* retains a direction-consistent nominal MR association in UKB (*P* = 0.0056). The chromosome 16p11.2 candidates are weak or ambiguous throughout and *MAPKAPK5* is PP.H3-dominant in discovery.
