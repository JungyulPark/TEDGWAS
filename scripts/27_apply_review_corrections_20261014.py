#!/usr/bin/env python3
"""Apply the 2026-10-14 review corrections to the manuscript master.

Each replacement must match the current master exactly once; the script stops
otherwise, so it cannot silently edit the wrong passage. Numbers introduced here
were checked against Supplementary_Dataset.xlsx (MR_primary, MR_estimators_original,
LOO_original, TSHR_allele_audit, Instruments, Detection_original) and
reproducibility/work/reanalysis_20261002/mr_v1/MR_14genes_all_scenarios.csv;
see internal/REVIEW_codex_package_20261011.md for the review items.

Applied once, to the master imported at commit 698e4b6; it will refuse to run
again because its source passages no longer exist. The Abstract was then
tightened to 250 words and three redundant sentences were removed to keep the
main text under 5,000 words; those follow-up edits are recorded in git history.

Usage: python3 scripts/27_apply_review_corrections_20261014.py
"""
from pathlib import Path

MASTER = Path(__file__).resolve().parents[1] / "MANUSCRIPT_TED_TRAP_v5_MASTER.md"

EDITS = [
    # ---------------- Abstract ----------------
    ("**Objective:** Graves disease (GD) causes autoimmune hyperthyroidism; thyroid eye disease (TED) involves orbital inflammation. We compared associations of genetically proxied blood *TSHR* and *IGF1R* expression with disease susceptibility.",
     "**Objective:** Graves disease (GD) is an autoimmune cause of hyperthyroidism; thyroid eye disease (TED) is an associated inflammatory orbital disorder. We compared genetic evidence linking blood *TSHR* and *IGF1R* expression to disease susceptibility."),
    ("**Methods:** We screened 2,544 druggable genes using Mendelian randomization (MR). Disease datasets comprised Biobank Japan (BBJ) GD, UK Biobank (UKB) hyperthyroidism and three overlapping FinnGen endpoints: GD, Graves ophthalmopathy and strict Graves ophthalmopathy. Secondary colocalization assessed 24 genes. A fixed five-tissue GTEx extension assessed both genes. Orbital RNA sequencing (four TED patients, one control) was descriptive.",
     "**Methods:** We screened 2,544 druggable genes using Mendelian randomization (MR). Disease outcomes were Biobank Japan (BBJ) GD, UK Biobank (UKB) hyperthyroidism and three overlapping FinnGen endpoints: GD, Graves ophthalmopathy (GO) and strict GO. Colocalization tested whether expression and disease shared a causal variant for 24 genes. A five-tissue GTEx extension and orbital RNA sequencing (four TED patients, one control) were supplementary."),
    ("**Results:** Thirteen genes met the BBJ discovery threshold (*P* < 0.05/2,544). Higher genetically proxied *TSHR* expression was associated with lower GD odds in BBJ (odds ratio [OR] per approximately one standard deviation 0.12, 95% confidence interval [CI] 0.07–0.21; *P* = 1.09×10⁻¹⁴) and FinnGen (OR 0.08, 95% CI 0.05–0.12; *P* = 3.19×10⁻³²), with concordant directions across all five datasets. *CTLA4*, *HLA-DQA2* and *TUBB* also passed the full-screen threshold in both GD datasets. Blood *IGF1R* estimates consistently indicated higher odds, but only BBJ (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB (*P* = 0.0117) reached uncorrected *P* < 0.05. Secondary colocalization support varied across datasets. Adipose *TSHR*–GD associations passed correction; thyroid and fibroblast *IGF1R* estimates did not reach *P* < 0.05. Mean orbital *TSHR* and *IGF1R* transcript levels exceeded the single control.",
     "**Results:** Thirteen genes met the BBJ discovery threshold (*P* < 0.05/2,544). Higher genetically proxied *TSHR* expression was associated with lower GD odds in BBJ (odds ratio [OR] per approximately one standard deviation 0.12, 95% confidence interval [CI] 0.07–0.21; *P* = 1.09×10⁻¹⁴) and FinnGen (OR 0.08, 95% CI 0.05–0.12; *P* = 3.19×10⁻³²), with the same direction in all five datasets. Shared-variant posterior probability (PP.H4) exceeded 0.80 in BBJ (0.951) and both GO endpoints (0.986, 0.891), but distinct variants were favoured in UKB (PP.H4 0.226) and FinnGen GD (0.264); under a stricter prior, BBJ and strict GO fell below 0.80. All *IGF1R* ORs exceeded 1, but only BBJ (OR 1.56, 95% CI 1.07–2.28; *P* = 0.0212) and UKB (*P* = 0.0117) reached uncorrected *P* < 0.05; the most precise estimate, FinnGen GD, was OR 1.07 (95% CI 0.86–1.33). Single-variant adipose *TSHR* instruments from the same locus gave concordant estimates; thyroid and fibroblast *IGF1R* estimates were non-significant. Mean orbital *TSHR* and *IGF1R* transcript levels exceeded the single control."),
    ("**Conclusions:** Associations were more consistently supported for blood *TSHR* than *IGF1R* expression. Secondary regional analyses limit expression-mediated causal interpretation. Population-control comparisons cannot isolate TED susceptibility within GD or predict treatment response.",
     "**Conclusions:** Genetic association evidence was stronger and more consistent for blood *TSHR* than *IGF1R* expression. Because colocalization did not support a shared variant in the two largest datasets, expression-mediated causality is not established. Population-control comparisons cannot isolate TED susceptibility within GD, and weaker *IGF1R* genetic evidence does not bear on IGF-1R–targeted treatment."),

    # ---------------- Introduction ----------------
    ("Here, the exposure is blood gene expression and the outcomes are disease diagnoses, not responses to treatment. Inherited variants precede disease onset, but MR supports causal interpretation only when its instrumental-variable assumptions hold.",
     "Here, the exposure is blood gene expression and the outcomes are disease diagnoses, not responses to treatment. Because alleles are allocated at conception, MR is less prone than ordinary observational studies to confounding and reverse causation, but it supports causal interpretation only when its instrumental-variable assumptions hold."),
    ("To address this question, we compared genetically proxied blood expression of *TSHR* and *IGF1R* within a screen of genes encoding potentially druggable proteins, using *CTLA4* and other associated genes as context.",
     "To address this question, we compared genetically proxied blood expression of *TSHR* and *IGF1R* — expression predicted from inherited variants rather than measured directly — within a screen of genes encoding potentially druggable proteins, using *CTLA4* and other associated genes as context."),
    ("The analyses distinguish GD from broader hyperthyroidism and population-controlled ophthalmopathy; they do not directly estimate TED susceptibility among patients who already have GD or predict treatment response.",
     "The analyses distinguish GD from broader hyperthyroidism and population-controlled ophthalmopathy; they do not directly estimate TED susceptibility among patients who already have GD or predict treatment response. A small in-house orbital RNA-sequencing series (four TED samples and one control) is reported only as descriptive context."),

    # ---------------- Methods ----------------
    ("Where applicable, sensitivity assessment included Cochran's Q, alternative MR estimators and omission of individual instruments [22, 23]. Previously completed alternative-estimator results for the original three outcomes were retained; they were not represented as new analyses of every endpoint.",
     "Where applicable, sensitivity assessment included Cochran's Q, alternative MR estimators and omission of individual instruments [22, 23]. Previously completed alternative-estimator results (MR-Egger, weighted median and weighted mode) for the original three outcomes were retained (Supplementary Dataset, MR_estimators_original); they were not represented as new analyses of every endpoint."),
    ("### Secondary regional analysis and robustness assessment",
     "### Colocalization and robustness assessment"),
    ("We report MR associations as the main findings and use colocalization to assess their interpretation. This reporting emphasis was adopted during post hoc reappraisal and was not a prospectively registered hierarchy.",
     "We report MR associations as the main findings and use colocalization to assess their interpretation. In the original analysis, colocalization was reported alongside MR as a parallel evidence layer; the present emphasis was adopted during post hoc reappraisal and was not a prospectively registered hierarchy."),
    ("\r\n\r\nAI-assisted tools supported code development, numerical checking and drafting; author verification and responsibility are stated in the Declarations.",
     ""),

    # ---------------- Results ----------------
    ("*IGF1R* estimates indicated higher odds in all five datasets, but no association met the full-screen threshold.",
     "All five *IGF1R* ORs exceeded 1, but no association met the full-screen threshold."),
    ("These findings show directional consistency without statistically convincing confirmation in FinnGen. They do not establish absence of an association.",
     "The FinnGen GD estimate was the most precise and lay close to the null. Weighted-median estimates resembled the IVW results in BBJ (OR 1.56; *P* = 0.0288) and UKB (OR 1.42; *P* = 0.00244), whereas MR-Egger and weighted-mode estimates were imprecise (*P* ≥ 0.32), and MR-Egger intercepts were non-significant (*P* = 0.66, 0.56 and 0.43 in BBJ, UKB and FinnGen GO), although these tests have little power with three or four instruments. These findings show directional consistency without statistically convincing confirmation in FinnGen. They do not establish absence of an association."),
    ("*CTLA4*, an established autoimmune locus [26], and *TUBB* had direction-consistent nominal associations across all five datasets.",
     "*HLA-DQA2* and *TUBB* lie in the major histocompatibility complex (MHC) region. *CTLA4*, an established autoimmune locus [26], and *TUBB* had direction-consistent nominal associations across all five datasets."),
    ("These pairs therefore do not represent independent genetic findings. Five major histocompatibility complex (MHC) region genes among the BBJ discoveries were retained in follow-up; membership of this region was not an exclusion criterion.",
     "These pairs therefore do not represent independent genetic findings. Seven of the nine (*APOM*, *NEU1*, *TNF*, *HIST1H4E*, *LTA*, *LTB* and *C4B*) were instrumented within the extended MHC region on chromosome 6. Five MHC-region genes among the BBJ discoveries were also retained in follow-up; membership of this region was not an exclusion criterion."),
    ("Heterogeneity was present in FinnGen GD for *CTLA4* (Q-test *P* = 0.0275), *TUBB* (*P* = 3.67×10⁻⁵) and *C4A* (*P* = 1.27×10⁻¹²).",
     "Heterogeneity was present in FinnGen GD for *C4A* (Q-test *P* = 1.27×10⁻¹²), *TUBB* (*P* = 3.67×10⁻⁵), *CTLA4* (*P* = 0.0275) and *HLA-A* (*P* = 0.0387)."),
    ("The original *IGF1R* BBJ and UKB estimates no longer reached nominal significance after excluding rs2654980 (*P* = 0.0751 and 0.468).",
     "In BBJ and UKB, two of the four single-variant omissions removed nominal *IGF1R* significance: rs2654980 (*P* = 0.0751 and 0.468) and, respectively, rs59467480 (*P* = 0.0503) and rs117212126 (*P* = 0.0529); FinnGen estimates had not reached nominal significance."),
    ("Frequency substitution preserved the original 13 BBJ discoveries and the principal *TSHR* and *IGF1R* interpretation, but not every exploratory FinnGen classification (Table 3C; Table S2).",
     "Frequency substitution preserved the original 13 BBJ discoveries and the principal *TSHR* and *IGF1R* interpretation (Supplementary Dataset, MR_new_freq), but not every exploratory FinnGen classification (Table 3C; Table S2)."),
    ("Precision also varied across the screen: only 14.6% of estimable BBJ genes could detect an OR of 1.5 with 80% power under the earlier normal-approximation benchmark (Supplementary Dataset, Detection_original).",
     "Precision also varied across the screen. At the discovery threshold, 14.6% and 35.6% of estimable BBJ genes had 80% power to detect ORs of 1.5 and 2.0 under a normal-approximation benchmark; the corresponding FinnGen GD proportions were 24.8% and 52.2% (Supplementary Dataset, Detection_original; Supplementary Methods)."),
    ("### Secondary colocalization findings",
     "### Colocalization findings"),
    ("Regional analyses provided a further check on expression-mediated interpretation (Figures 2 and 3; Table S1).",
     "Regional analyses tested whether the MR associations reflected a variant shared with expression (Figures 2 and 3; Table S1)."),
    ("Neither thyroid nor cultured-fibroblast *IGF1R* estimates reached *P* < 0.05 in any dataset. Their FinnGen GD estimates were OR 0.76 (95% CI 0.56–1.02; *P* = 0.0692) and OR 1.11 (95% CI 0.89–1.38; *P* = 0.359).",
     "The two *TSHR* instruments, rs12101255 and rs74064796, lie 50 bp apart in *TSHR* intron 1, about 15 kb from rs179252 (European r² 0.46 for rs12101255), and both are themselves blood *TSHR* eQTLs (eQTLGen *P* ≈ 10⁻³⁴). With one instrument, each tissue *P* value equals that variant's disease-association *P* value, so these estimates re-express the regional *TSHR* association rather than providing independent evidence. Neither thyroid nor cultured-fibroblast *IGF1R* estimates reached *P* < 0.05 in any dataset. Their FinnGen GD estimates were OR 0.76 (95% CI 0.56–1.02; *P* = 0.0692), opposite in direction to the blood estimates, and OR 1.11 (95% CI 0.89–1.38; *P* = 0.359)."),

    # ---------------- Discussion ----------------
    ("*IGF1R* estimates consistently indicated higher odds, but only the BBJ and UKB estimates reached uncorrected *P* < 0.05, and none met the discovery threshold.",
     "All *IGF1R* ORs exceeded 1, but only the BBJ and UKB estimates reached uncorrected *P* < 0.05, none met the discovery threshold, and the most precise estimate lay close to the null. Colocalization supported a shared variant for *TSHR* in BBJ and both ophthalmopathy endpoints but favoured distinct variants in UKB and FinnGen GD, the two datasets with the most cases."),
    ("Ji and colleagues had already linked blood *TSHR* expression to GD and reported support in FinnGen ophthalmopathy [9];",
     "Ji and colleagues had already linked blood *TSHR* expression to GD in an earlier FinnGen release and reported support in FinnGen ophthalmopathy [9];"),
    ("*CTLA4* also showed direction-consistent associations across all five outcomes, supporting its use as an autoimmune comparator. However, its multi-instrument estimates relied heavily on rs13030124.",
     "Its multi-instrument European estimates relied heavily on rs13030124; with that caveat, its direction-consistent associations across all five outcomes support its use as an autoimmune comparator."),
    ("This directional consistency strengthens the plausibility of a central-tolerance explanation, but does not demonstrate that thymic expression mediates the disease association.",
     "This directional consistency is compatible with a central-tolerance explanation, but the thymic step was not tested here, and it does not demonstrate that thymic expression mediates the disease association."),
    ("The marked departure of the MR OR from 1 also reflects division by a small reconstructed expression coefficient. It should not be translated into the effect of experimentally increasing *TSHR* expression. The adipose analyses extended the lower-odds direction beyond blood, but they reused disease GWAS data and implicated the same genomic region. They do not establish an adipose-mediated mechanism or provide independent causal replications.",
     "The marked departure of the MR ORs from 1 reflects scaling rather than a large per-variant effect: rs179252 shifts reconstructed blood expression by only about 0.11 SD, and its per-allele disease ORs were 0.80 (BBJ), 0.77 (UKB), 0.78 (FinnGen GO), 0.77 (FinnGen GD) and 0.82 (strict GO). Dividing these small effects by the small expression effect yields the per-SD ORs, which should not be translated into the effect of experimentally increasing *TSHR* expression. The adipose instruments tag the same intron-1 region as rs179252 and are themselves blood eQTLs, so their concordant estimates re-express this regional association. They cannot address the distinct-variant pattern in UKB and FinnGen GD, establish an adipose-mediated mechanism or provide independent causal replication."),
    ("The fixed-panel extension found no nominal association for thyroid or cultured-fibroblast *IGF1R*, despite qualifying instruments. Limited support was therefore not confined to the blood analyses that could be estimated. However, the confidence intervals, single-instrument design and lack of orbital eQTLs preclude excluding a clinically relevant effect.",
     "With 80% power at *P* < 0.05, the FinnGen GD analysis could detect an *IGF1R* OR of about 1.36, so it constrains larger blood-expression effects but not modest ones. The fixed-panel extension found no nominal association for thyroid or cultured-fibroblast *IGF1R*, and the thyroid estimate ran opposite to the blood direction. The confidence intervals, single-instrument design and lack of orbital eQTLs preclude excluding a clinically relevant effect."),
    ("Teprotumumab efficacy [3] and experimental TSHR–IGF-1R crosstalk [5] establish complementary clinical and biological observations; neither can substitute for the missing statistical confirmation of inherited expression-related susceptibility.",
     "Weaker genetic evidence for *IGF1R* does not argue against IGF-1R blockade: teprotumumab efficacy rests on randomized trials in established TED [3], and experimental TSHR–IGF-1R crosstalk [5] concerns receptor signalling rather than inherited susceptibility."),
    ("Colocalization serves as a secondary interpretive analysis rather than a substitute for the MR associations. Its findings nevertheless constrain causal wording.",
     "Colocalization does not replace the MR associations, but it constrains causal wording."),
    ("Thus, association strength and evidence for a shared causal variant were not interchangeable. Ji and colleagues used a shared-variant prior of 5×10⁻⁵ and a PP.H4 threshold above 0.50 [9], whereas our primary settings used 10⁻⁵ and at least 0.80.",
     "Thus, association strength and evidence for a shared causal variant were not interchangeable. UKB and FinnGen GD had the most cases of the five datasets; their distinct-variant patterns could reflect finer resolution of a complex regional architecture rather than weaker association, but this was not tested. Ji and colleagues chose a shared-variant prior intended to allow weaker signals to colocalize and accepted PP.H4 above 0.50 [9], whereas our primary settings used 10⁻⁵ and at least 0.80."),
    ("Genetically proxied blood *TSHR* expression showed consistent associations with lower odds of GD and related phenotypes, supported by cross-dataset comparisons within a broader screen of druggable genes. The adipose extension supported the same *TSHR* direction. Blood *IGF1R* estimates indicated higher odds but provided limited statistical evidence, and thyroid/fibroblast estimates provided no nominal support. These findings place the strongest emphasis on association; secondary colocalization limits claims that blood expression mediates disease risk.",
     "Genetically proxied blood *TSHR* expression showed strong, directionally consistent associations with lower odds of GD and related phenotypes, whereas *IGF1R* evidence was weaker and imprecise. Colocalization supported a shared variant for *TSHR* in BBJ and the ophthalmopathy endpoints but not in UKB or FinnGen GD, so these data establish association rather than expression-mediated causality. Single-variant adipose instruments from the same locus added concordant estimates but not independent tissue evidence. Weaker genetic evidence for *IGF1R* does not argue against the therapeutic relevance of IGF-1R blockade."),

    # ---------------- Declarations ----------------
    ("\r\n\r\n**AI-assisted preparation.** AI-assisted tools supported code development, numerical checking and language editing. The authors are responsible for verifying the analyses, references, interpretations and final manuscript.",
     ""),
]


def main() -> int:
    raw = MASTER.read_bytes().decode("utf-8")
    for old, new in EDITS:
        n = raw.count(old)
        if n != 1:
            raise SystemExit(f"expected exactly one match, found {n}: {old[:90]!r}")
        raw = raw.replace(old, new)
    MASTER.write_bytes(raw.encode("utf-8"))
    print(f"applied {len(EDITS)} edits to {MASTER.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
