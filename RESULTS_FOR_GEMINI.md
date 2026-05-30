# RESULTS — full prose (for Gemini documentation)

All 6 sections FINAL & verified against locked ground truth. Gene symbols italic in manuscript.


## R1. Study overview
*Assets: Figure 1, Table 1*

We applied a druggable-gene-wide MR, colocalization, and orbital tissue-integration framework to distinguish genetically anchored susceptibility from pharmacologic effector biology in GD and TED (**Figure 1**). The analysis followed a prespecified outcome hierarchy using Biobank Japan (BBJ) Graves disease as the primary discovery outcome, UK Biobank (UKB) hyperthyroidism as replication, and FinnGen R12 Graves ophthalmopathy (a TED phenotype) as the TED-specific sensitivity outcome, with data sources summarized in **Table 1**.


## R2. Druggable-gene-wide MR discovery screen
*Assets: Figure 2, Table S1, Table S4*

Of 4,462 druggable genes, 2,545 had at least one valid *cis*-eQTL instrument in eQTLGen and formed the MR-testable denominator; across these genes, 6,136 SNP-level instruments were retained, all exceeding conventional instrument-strength criteria (minimum F-statistic = 14.2, no weak instruments; **Table S4**). Primary MR estimates against the BBJ Graves disease discovery outcome were obtainable for 2,235 genes.

Applying a Bonferroni threshold corrected for the MR-testable denominator (*P* < 1.965×10⁻⁵; 0.05/2,545), 13 genes reached significance in discovery (**Figure 2**). These hits included the established Graves disease autoantigen *TSHR* and the known autoimmune locus *CTLA4*, both with strongly significant protective-direction estimates, together with additional loci that were carried forward for replication and colocalization-based filtering. The discovery volcano (**Figure 2**) shows the distribution of BBJ MR effect estimates and the 13 Bonferroni-significant hits; full hit classification is provided in **Table S1**.


## R3. Cross-outcome classification and colocalization filtering
*Assets: Table S1, Table S2, Figure S2, Table S5*

The 13 discovery hits were not treated as putative targets on the strength of the discovery MR alone. Instead, each was classified according to prior biological status, genomic context, cross-outcome support, and colocalization with both the discovery and TED-specific signals. Two hits—*TSHR* and *CTLA4*—were recognized as established Graves disease and autoimmune loci, respectively, and are examined in detail below. The remaining hits were evaluated as non-known candidates, including a cluster on chromosome 16p11.2 (*HSD3B7*, *VKORC1*, and *PRSS36*) whose top eQTL instruments were co-located within a narrow approximately 143-kb linkage-disequilibrium window (Table S5).

Candidate colocalization showed outcome-dependent patterns (Table S2; Figure S2). *TNFSF14* and *IFNGR1* showed strong colocalization with the BBJ discovery signal (PP.H4 = 0.994 and 0.989, respectively), but did not colocalize with the TED-specific outcome (PP.H4 = 0.019 for both), where the posterior favored an outcome-only association pattern rather than shared eQTL–GWAS architecture. *MAPKAPK5* showed a distinct-signal pattern in discovery (PP.H3 ≈ 1.0), whereas the chr16p11.2 candidates showed only weak or ambiguous colocalization in discovery (PP.H4 ≤ 0.62) and no shared colocalization support in the TED-specific outcome.

Overall, no non-known candidate showed strong colocalization in both BBJ Graves disease and the TED-specific outcome. Thus, after cross-outcome colocalization filtering, the screen yielded no robust novel target. This filtering result was retained as part of the final interpretation because it prevented overinterpretation of single-outcome or LD-driven discovery signals and strengthened the focus on the prespecified backbone genes *TSHR*, *IGF1R*, and *CTLA4*.


## R4. TSHR as an expression-colocalized susceptibility anchor
*Assets: Table 2, Figure 3, Table S3, Figure S1*

Among the established hits, *TSHR* showed the strongest convergence across the genetic, colocalization, fine-mapping, and tissue layers of the framework (Table 2; Figure 3). In MR, the *TSHR* *cis*-eQTL instrument was associated with protective-direction estimates across all three outcomes: BBJ Graves disease (β = −2.10, *P* = 1.1×10⁻¹⁴), UKB hyperthyroidism (β = −2.44, *P* = 8.8×10⁻²⁸), and the TED-specific outcome (β = −2.33, *P* = 2.8×10⁻⁷). The consistency of direction and approximate magnitude across discovery, replication, and TED-specific analysis argues against an outcome-specific artifact (Figure 3B).

Colocalization supported a shared causal variant between *TSHR* expression and disease association in both the discovery and TED-specific outcomes (PP.H4 = 0.951 and 0.986, respectively; Figure 3A), with rs179252 prioritized as the shared lead variant (Table 2). Fine-mapping and ancestry-matched LD assessment supported retention of rs179252 as the single-IV/locus-level anchor: secondary credible-set lead variants remained in substantial LD with rs179252 and did not provide clearly independent additional instruments (Table S3; Figure S1).

At the tissue level, *TSHR* was significantly upregulated in TED orbital tissue relative to control tissue (log2 fold-change = +2.33, adjusted *P* = 0.032; Figure 3C), providing exploratory tissue-level support for the biological relevance of the prioritized locus. Together, these results position *TSHR* as an expression-colocalized upstream susceptibility anchor in GD/TED rather than merely a rediscovered Graves disease association. The apparent contrast between the protective-direction germline MR estimate and disease-state orbital upregulation is addressed in the Discussion.


## R5. IGF1R as a non-colocalized pharmacologic effector axis
*Assets: Table 2, Figure 3, Table S6*

In contrast to *TSHR*, *IGF1R*, the therapeutic target of teprotumumab, showed a markedly different evidence profile (Table 2; Figure 3). *IGF1R* was instrumented by multiple independent *cis*-eQTL variants, with four instruments in the discovery and replication analyses. The inverse-variance-weighted MR estimate indicated a nominal positive-direction association with BBJ Graves disease (β = +0.45, *P* = 0.021) and UKB hyperthyroidism (β = +0.30, *P* = 0.012); the TED-specific estimate was directionally consistent but non-significant (β = +0.34, *P* = 0.18). Sensitivity estimators were broadly concordant in direction without evidence of directional pleiotropy, although statistical significance was not uniform across estimators, consistent with a modest effect (Table S6).

Critically, this nominal MR signal was not accompanied by expression colocalization. At the *IGF1R* locus, the posterior favored an outcome-only association pattern rather than a shared causal variant with expression in both the discovery and TED-specific outcomes (PP.H4 = 0.043 and 0.036; H2-dominant; Figure 3A), in contrast to the strong shared-variant evidence observed at *TSHR*. At the tissue level, *IGF1R* was not significantly differentially expressed in TED orbital tissue (log2 fold-change = +0.41, adjusted *P* = 1.00; Figure 3C).

Thus, *IGF1R* showed nominal positive-direction MR evidence that was not supported by eQTL–GWAS colocalization or by significant orbital tissue upregulation. This profile is consistent with *IGF1R* acting as a pharmacologic effector axis, rather than an expression-colocalized germline susceptibility locus comparable to *TSHR*, and aligns with the clinical rationale for IGF-1R blockade in TED.


## R6. CTLA4 as a positive-control autoimmune locus
*Assets: Table 2, Figure 3, Table S6*

*CTLA4* (chromosome 2q33), a well-established autoimmune susceptibility locus, was included as a positive control to calibrate the framework against a gene of known relevance to Graves disease. Consistent with its role as a known autoimmune locus, *CTLA4* showed strong protective-direction MR estimates in the discovery analysis (BBJ Graves disease, β = −1.74, *P* = 5.5×10⁻¹⁵), with concordant protective-direction estimates in replication (UKB hyperthyroidism, β = −1.57, *P* = 0.002) and the TED-specific outcome (β = −1.77, *P* = 0.012). The recovery of *CTLA4* as a strong autoimmune signal supports the sensitivity of the discovery screen.

Colocalization at *CTLA4* was outcome-dependent: the posterior favored distinct causal variants in the discovery outcome (PP.H3 = 0.79), but strong colocalization in the TED-specific outcome (PP.H4 = 0.978; Figure 3A). This mixed pattern, together with heterogeneity across instruments in the replication analysis (Table S6), is consistent with locus and cross-ancestry complexity at *CTLA4* rather than a uniformly colocalized signal. We therefore interpret *CTLA4* as a positive-control autoimmune locus rather than a uniformly robust expression-colocalized anchor. At the tissue level, *CTLA4* was not significantly differentially expressed in TED orbital tissue (log2 fold-change = +1.27, adjusted *P* = 0.815; Figure 3C), consistent with an immune-system signal that was not accompanied by a significant local orbital expression change in this exploratory tissue dataset.
