# STROBE MR reporting checklist

TSHR and IGF1R expression in Graves disease: Mendelian randomization and colocalization

This adapted reporting map identifies where each item is addressed and explicitly records unavailable information. It does not certify that all methodological limitations are resolved.

| Item | Reporting topic | Location and coverage |
|---|---|---|
| 1 | MR design | Abstract identifies Mendelian randomization; Methods specifies the two-sample design. |
| 2 | Rationale | Introduction: blood expression, disease susceptibility and therapeutic interpretation. |
| 3 | Objectives | Introduction final paragraph; Methods defines the hierarchy and biological backbone. No prospective registration identifier is available. |
| 4a–b | Design, participants and sample size | Methods; Table 1; source reports [12, 13, 14, 15]. Recruitment and eligibility details are those of the original studies; available sample sizes determined the analysis. Detection thresholds were calculated from observed MR standard errors, not for prospective recruitment. |
| 4c–d | Variants and phenotype definitions | Methods and Supplementary Methods describe selection and harmonization; Table 1 gives accessions; Supplementary Data 1 lists instruments. Source studies provide their genotyping and phenotype procedures. |
| 4e | Ethics and consent | Declarations; in-house IRB approval and written consent. |
| 5 | Instrument assumptions | Methods states relevance, independence and exclusion restriction; Discussion discusses limits of verification. |
| 6a–b | Scales and variant weights | Methods: Z-score reconstruction, UKB log-odds rescaling, Wald/IVW and inverse-variance weighting. |
| 6c–e | Estimation, covariates, missingness and multiplicity | Methods: source GWAS adjustments retained; no individual-level re-adjustment; unavailable or excluded variants omitted; Bonferroni denominator 2,544. |
| 7 | Assumption assessment | Methods: instrument strength, harmonization, heterogeneity, MR-Egger and colocalization. These do not prove all instrument assumptions. |
| 8 | Additional analyses | Methods; Tables S1–S3. eQTLGen frequency substitution and re-harmonization were performed. |
| 9a–b | Software and registration | Methods lists R, TwoSampleMR, coloc and PLINK versions; Python was used for verification, frequency sensitivity and SNP exclusion, with an independent base-R check. No prospective registration identifier is available. |
| 10a–c | Descriptive data | Table 1 gives genetic-study sample counts; Supplementary Methods describes the orbital samples; Figure 1 gives gene attrition; Supplementary Data 1 gives instrument counts. Participant-level distributions and cohort-specific eQTL meta-analysis heterogeneity were not reanalysed. |
| 10d | Transportability and participant overlap | Discussion: European blood eQTL transfer to East Asian discovery is an assumption; overlap with European outcomes is unquantified. |
| 11a–d | Main results and uncertainty | Tables 2–3; Figure 2. ORs/CIs are per reconstructed expression unit. Variant-level instruments and primary MR results are supplied as Supplementary Data 1–2. Absolute risks were not estimated. |
| 12a–b | Results of assumption checks | Instrument strength in Methods; heterogeneity and MR-Egger in Table S1; clumping in Supplementary Methods; colocalization in Figure 3B and frequency sensitivity in Tables S2–S3. |
| 13a–e | Additional results | Tables S1–S3; Figure S1. Post hoc leave-one-out results are reported in Results, Figure S2, Supplementary Data 6 (15 omissions; five all-SNP comparators). Single-SNP starting sets were ineligible. Steiger and bidirectional MR were not applied; Figure S1 is descriptive. |
| 14 | Key findings | Discussion opening paragraph. |
| 15 | Limitations and bias | Discussion: tissue, ancestry, overlap, frequency reconstruction, phenotype, instrument counts, power, coloc model/priors and one tissue control; bias magnitudes were not quantified. |
| 16a–c | Interpretation and clinical meaning | Discussion: genetic expression effects do not establish a pharmacologic mechanism, treatment direction or intervention effect size. |
| 17 | Generalizability | Discussion: ancestry and tissue context, broad hyperthyroidism and population controls; inherited expression proxies do not estimate acute or dose-specific treatment effects. |
| 18 | Funding | Declarations states no specific grant; original data-source funding is reported in the cited source publications. |
| 19 | Data and code access | Declarations lists source repositories and access restrictions. Supplementary Data 1–6 provide instruments, primary results, frequency-sensitivity MR, coloc and power summaries, and post hoc leave-one-out estimates. Analysis code is available from the corresponding author on request. |
| 20 | Competing interests | Declarations. |

Adapted from the STROBE-MR checklist (EQUATOR Network, CC BY 3.0), https://www.strobe-mr.org/download/strobe-mr-checklist/. The checklist is a reporting map, not a claim that all study limitations have been resolved.