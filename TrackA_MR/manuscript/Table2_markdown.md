
## Table 2. Mendelian randomization and colocalization evidence for TSHR versus IGF1R causal effects on Graves disease and thyroid eye disease.

### Panel A — Primary cis-MR (Graves disease, Biobank Japan, East Asian)

|Gene  |Method                    | N_IV|Beta   |SE    |OR               |P_value  |Sig |
|:-----|:-------------------------|----:|:------|:-----|:----------------|:--------|:---|
|TSHR  |Wald ratio                |    1|-1.408 |0.182 |0.24 (0.17-0.35) |9.12e-15 |*** |
|IGF1R |Inverse variance weighted |   10|+0.200 |0.118 |1.22 (0.97-1.54) |0.089    |    |
|IGF1R |MR Egger                  |   10|+0.312 |0.235 |1.37 (0.86-2.17) |0.221    |    |
|IGF1R |Weighted median           |   10|+0.257 |0.157 |1.29 (0.95-1.76) |0.103    |    |
|IGF1R |Weighted mode             |   10|+0.257 |0.502 |1.29 (0.48-3.45) |0.621    |    |
|IGF1R |Simple mode               |   10|+0.287 |0.374 |1.33 (0.64-2.77) |0.463    |    |

### Panel B — Replication (UK Biobank, European; β rescaled to log-odds using Lloyd-Jones 2018 method)

|Gene  |Method                    | N_IV|Beta_rescaled |SE_rescaled |P_value  |Sig |
|:-----|:-------------------------|----:|:-------------|:-----------|:--------|:---|
|TSHR  |Wald ratio                |    1|-1.629        |0.149       |1.07e-27 |*** |
|IGF1R |Inverse variance weighted |   10|+0.152        |0.090       |0.094    |    |

### Panel C — Multivariable MR for TED-specific signal (outcome: FinnGen R12 Graves ophthalmopathy)

|Exposure                     |Method          | N_IV|Beta   |SE    |Conditional_F |P_value  |Sig |
|:----------------------------|:---------------|----:|:------|:-----|:-------------|:--------|:---|
|TSHR (adj. GD liability)     |MVMR-IVW        |   11|+0.127 |0.277 |25.4          |0.657    |    |
|GD liability (adj. TSHR)     |MVMR-IVW        |   11|+0.962 |0.074 |99.3          |3.67e-07 |*** |
|TSHR — Univariable reference |Univariable IVW |    1|-1.566 |0.305 |—             |2.82e-07 |*** |

### Panel D — Bayesian colocalization at TSHR locus (chr14:80.9-82.1 Mb)

|Outcome               |Cohort_N              |Common_SNPs          |PP_H3 |PP_H4        |Top_SNP   |Note                                                 |
|:---------------------|:---------------------|:--------------------|:-----|:------------|:---------|:----------------------------------------------------|
|Graves disease (BBJ)  |212,453 (2,809 cases) |2,707                |0.049 |0.951 ⭐⭐⭐ |rs179252  |Shared causal variant                                |
|Hyperthyroidism (UKB) |484,598 (3,731 cases) |3,376                |0.774 |0.227        |rs1023586 |Different lead tag; LD r²=0.85 in EAS with rs179252 |
|TED/GO (FinnGen R12)  |500,348 (858 cases)   |Technical limitation |NA    |NA           |—         |Technical limitation (attempted)                     |

### Footnotes

*** P < 2.08×10⁻³ (Bonferroni-adjusted for 8 genes × 3 outcomes = 24 tests).
* P < 0.05 (nominal).
β: log-odds effect per 1-SD increase in gene expression.
OR (95% CI): odds ratio computed as exp(β).
N_IV: number of independent instruments (r² < 0.001, clumping window 10 Mb).
Conditional F: Sanderson 2019 conditional F-statistic for instrument strength (>10 = strong).
PP.H4: posterior probability of shared causal variant (Giambartolomei 2014; default priors).
Panel B β values rescaled from UKB linear LMM to log-odds scale using factor 1/[p(1-p)] where p = case prevalence.

Primary outcome: Graves disease (Sakaue et al., 2021; GCST90018627).
Replication outcome: Hyperthyroidism (Dönertaş et al., 2021; GCST90038636).
Sensitivity/TED outcome: FinnGen Release 12, Graves ophthalmopathy endpoint (858 cases).
