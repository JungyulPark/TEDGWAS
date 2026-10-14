## Table S1. Additional FinnGen GD screen discoveries

| Gene | Chr | FinnGen GD IVs | FinnGen GD OR (95% CI) | FinnGen GD *P* | BBJ PP.H4 | FinnGen GD PP.H4 | FinnGen GO strict PP.H4 |
|----------|--------|--------|---------------|------------|---------|----------|-----------|
| *APOM* | 6 (MHC) | 1 | 0.05 (0.04–0.06) | 1.43×10^−108^ | <0.001 | <0.001 | 0.227 |
| *NEU1* | 6 (MHC) | 1 | 0.04 (0.03–0.06) | 1.43×10^−108^ | <0.001 | <0.001 | 0.061 |
| *IGLV2-18* | 22 | 1 | 0.56 (0.48–0.66) | 7.52×10^−12^ | 0.029 | 0.020 | 0.100 |
| *IGLV3-12* | 22 | 1 | 1.20 (1.14–1.26) | 7.52×10^−12^ | 0.019 | 0.022 | 0.101 |
| *TNF* | 6 (MHC) | 1 | 0.60 (0.51–0.71) | 4.65×10^−10^ | <0.001 | <0.001 | <0.001 |
| *HIST1H4E* | 6 (MHC) | 1 | 2.54 (1.88–3.43) | 1.51×10^−9^ | <0.001 | <0.001 | 0.007 |
| *LTA* | 6 (MHC) | 1 | 0.46 (0.35–0.61) | 3.75×10^−8^ | <0.001 | <0.001 | <0.001 |
| *LTB* | 6 (MHC) | 1 | 0.18 (0.10–0.36) | 4.89×10^−7^ | <0.001 | <0.001 | <0.001 |
| *C4B* | 6 (MHC) | 2 | 1.63 (1.31–2.03) | 1.10×10^−5^ | <0.001 | <0.001 | 0.311 |

Genes outside the BBJ discovery set that passed *P* < 0.05/2,544 in the post hoc full FinnGen GD screen (1000 Genomes European reference frequencies). These are not independent of the screen and seven of nine are in the MHC. *APOM*/*NEU1* share rs3101018 and *IGLV2-18*/*IGLV3-12* share rs5759410. The lead *IGLV* expression variant is absent in BBJ, weakening interpretation of its low PP.H4. None reached PP.H4 ≥0.80 (p12 = 10^−5^) in any dataset. Chr, chromosome; IV, instrumental variable; MHC, major histocompatibility complex; OR, odds ratio; PP.H4, shared-variant posterior probability; BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank.

## Table S2. Frequency-sensitive exploratory candidates

| Gene | Dataset | Frequency source | IVs | OR (95% CI) | *P* value | Q-test *P* |
|--------|---------------|-------------|-----|---------------|------------|------------|
| *C4B* | FinnGen GD | 1000G EUR | 2 | 1.63 (1.31–2.03) | 1.10×10^−5^ | 7.79×10^−7^ |
| *C4B* | FinnGen GD | eQTLGen cohort | 2 | 1.75 (1.31–2.33) | 1.48×10^−4^ | 1.36×10^−8^ |
| *C4B* | FinnGen GO strict | 1000G EUR | 2 | 1.68 (1.44–1.95) | 2.74×10^−11^ | 0.122 |
| *C4B* | FinnGen GO strict | eQTLGen cohort | 2 | 1.80 (1.46–2.23) | 6.17×10^−8^ | 0.0587 |
| *HLA-G* | FinnGen GD | 1000G EUR | 16 | 0.86 (0.72–1.03) | 0.107 | 9.98×10^−24^ |
| *HLA-G* | FinnGen GD | eQTLGen cohort | 16 | 0.60 (0.47–0.75) | 1.72×10^−5^ | 2.26×10^−10^ |
| *HLA-G* | FinnGen GO strict | 1000G EUR | 16 | 0.91 (0.76–1.08) | 0.267 | 0.0231 |
| *HLA-G* | FinnGen GO strict | eQTLGen cohort | 16 | 0.70 (0.53–0.92) | 0.0122 | 0.131 |

Exposure effects were reconstructed with 1000 Genomes European reference allele frequencies (primary) or eQTLGen cohort frequencies; variant sets and outcome effects were held fixed, so changes reflect rescaled expression effects and IVW weights. *HLA-G* was added because its eQTLGen-frequency *P* crossed the screen threshold (*P* < 0.05/2,544). Q tests heterogeneity, not directional pleiotropy; its *P* value is not a gene-association *P* value. All 72 omission estimates and instrument weights are in the workbook. CI, confidence interval; IV, instrumental variable; IVW, inverse-variance weighted; OR, odds ratio; BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank.

## Table S3. *TSHR* regional colocalization and prior sensitivity

| Dataset | Matched variants | Disease lead (*P*) | rs179252 *P* | rs179252-G OR (95% CI) | H4 top variant | PP.H4, 10^−5^ / 10^−6^ | PP.H4 range across windows |
|--------------|---------|--------------------|----------|-----------------|-----------|------------|------------|
| BBJ GD | 4,336 | rs28414437 (1.83×10^−15^) | 1.09×10^−14^ | 0.80 (0.76–0.85) | rs179252 | **0.951** / 0.661 | 0.951–0.951 |
| FinnGen GD | 6,579 | rs1023586 (1.58×10^−36^) | 3.19×10^−32^ | 0.77 (0.73–0.80) | rs1023586 | 0.264 / 0.035 | 0.264–0.264 |
| UKB hyperthyroidism | 5,667 | rs4903961 (1.17×10^−43^) | 8.77×10^−28^ | 0.77 (0.74–0.81) | rs1023586 | 0.226 / 0.028 | 0.226–0.226 |
| FinnGen GO | 6,579 | rs3783947 (1.23×10^−7^) | 2.82×10^−7^ | 0.78 (0.71–0.86) | rs179252 | **0.986** / **0.875** | 0.986–0.987 |
| FinnGen GO strict | 6,579 | rs2110695 (1.34×10^−5^) | 9.02×10^−5^ | 0.82 (0.74–0.90) | rs179252 | **0.891** / 0.451 | 0.891–0.905 |

Disease lead, the strongest disease-associated variant among the matched cis variants, with its *P* value; rs179252 *P*, z-test *P* of the aligned G-allele disease association (the UKB source BOLT-LMM *P* was 1.3×10^−27^). H4 top variant, the variant with the highest posterior under the shared-variant hypothesis; it is not necessarily the disease lead. PP.H4 under p12 = 10^−5^ and 10^−6^ (p1 = p2 = 10^−4^); bold, ≥0.80. Window range covers the full source cis-region, ±500 kb, ±250 kb and the five-outcome intersection at p12 = 10^−5^; it is not a confidence interval. All 240 window settings are in TSHR_windows. BBJ pairs European eQTL with East Asian disease statistics. CI, confidence interval; OR, odds ratio; PP.H4, shared-variant posterior probability; BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank.

## Table S4. *TSHR* and *IGF1R* tissue expression extension

### A. Availability

| Gene / tissue | GTEx N | Instrument (F) | Status |
|-------------------------|--------|--------------------|----------------------------------------|
| *TSHR* / Thyroid | 574 | — | No instrument (minimum eQTL *P* = 9.60×10^−8^) |
| *IGF1R* / Thyroid | 574 | rs1815009 (F = 61.1) | Estimable |
| *TSHR* / Cultured fibroblasts | 483 | — | Not tested by GTEx in this tissue |
| *IGF1R* / Cultured fibroblasts | 483 | rs6598538 (F = 33.8) | Estimable |
| *TSHR* / Subcutaneous adipose | 581 | rs12101255 (F = 38.5) | Estimable |
| *IGF1R* / Subcutaneous adipose | 581 | — | No instrument (minimum eQTL *P* = 3.56×10^−7^) |
| *TSHR* / Visceral adipose | 469 | rs74064796 (F = 32.1) | Estimable |
| *IGF1R* / Visceral adipose | 469 | — | No instrument (minimum eQTL *P* = 7.76×10^−6^) |
| *TSHR* / Whole blood | 670 | — | No instrument (minimum eQTL *P* = 2.46×10^−4^) |
| *IGF1R* / Whole blood | 670 | — | No instrument (minimum eQTL *P* = 5.42×10^−5^) |

### B. Estimates

| Gene / tissue | Dataset | Family | OR (95% CI) | *P* value | Below family threshold | Direction vs blood |
|-------------------------|----------------|---------|---------------|-----------|----------|----------|
| *TSHR* / Subcutaneous adipose | BBJ GD | Secondary | 0.46 (0.38–0.55) | 1.92×10^−15^ | Yes | Same |
| *TSHR* / Subcutaneous adipose | FinnGen GD | Primary | 0.36 (0.31–0.43) | 1.28×10^−34^ | Yes | Same |
| *TSHR* / Subcutaneous adipose | UKB hyperthyroidism | Secondary | 0.31 (0.27–0.37) | 2.86×10^−42^ | Yes | Same |
| *TSHR* / Subcutaneous adipose | FinnGen GO | Secondary | 0.45 (0.32–0.63) | 5.76×10^−6^ | Yes | Same |
| *TSHR* / Subcutaneous adipose | FinnGen GO strict | Secondary | 0.51 (0.35–0.74) | 3.67×10^−4^ | Yes | Same |
| *TSHR* / Visceral adipose | BBJ GD | Secondary | 0.56 (0.46–0.67) | 2.64×10^−10^ | Yes | Same |
| *TSHR* / Visceral adipose | FinnGen GD | Primary | 0.50 (0.43–0.58) | 4.09×10^−19^ | Yes | Same |
| *TSHR* / Visceral adipose | UKB hyperthyroidism | Secondary | 0.39 (0.34–0.45) | 9.19×10^−38^ | Yes | Same |
| *TSHR* / Visceral adipose | FinnGen GO | Secondary | 0.69 (0.49–0.98) | 0.0368 | No | Same |
| *TSHR* / Visceral adipose | FinnGen GO strict | Secondary | 0.78 (0.54–1.13) | 0.186 | No | Same |
| *IGF1R* / Thyroid | BBJ GD | Secondary | 1.23 (0.90–1.68) | 0.200 | No | Same |
| *IGF1R* / Thyroid | FinnGen GD | Primary | 0.76 (0.56–1.02) | 0.0692 | No | Opposite |
| *IGF1R* / Thyroid | UKB hyperthyroidism | Secondary | 1.33 (0.99–1.80) | 0.0594 | No | Same |
| *IGF1R* / Thyroid | FinnGen GO | Secondary | 1.02 (0.54–1.94) | 0.945 | No | Same |
| *IGF1R* / Thyroid | FinnGen GO strict | Secondary | 1.23 (0.62–2.44) | 0.549 | No | Same |
| *IGF1R* / Cultured fibroblasts | BBJ GD | Secondary | 0.91 (0.70–1.20) | 0.511 | No | Opposite |
| *IGF1R* / Cultured fibroblasts | FinnGen GD | Primary | 1.11 (0.89–1.38) | 0.359 | No | Same |
| *IGF1R* / Cultured fibroblasts | UKB hyperthyroidism | Secondary | 1.03 (0.82–1.29) | 0.815 | No | Same |
| *IGF1R* / Cultured fibroblasts | FinnGen GO | Secondary | 0.78 (0.49–1.25) | 0.301 | No | Opposite |
| *IGF1R* / Cultured fibroblasts | FinnGen GO strict | Secondary | 0.78 (0.47–1.30) | 0.341 | No | Opposite |

Post hoc fixed-panel extension using GTEx v8; tissue colocalization was not performed. A, availability of the ten gene–tissue combinations; F statistics are for the single retained instrument. 'Not tested' means the gene is not among the genes GTEx tested in that tissue. GTEx whole blood (n = 670) yielded no instrument at the fixed threshold, whereas the main analysis used eQTLGen (up to 31,684 participants). B, all 20 estimates. The primary family is the ten FinnGen GD comparisons (*P* < 0.005); the secondary family is the 40 comparisons with the other datasets (*P* < 0.00125); unestimable combinations stay in both denominators. ORs are per unit of inverse-normal-transformed tissue expression and are not comparable in size with the blood ORs. With one instrument, each *P* value equals that variant's disease-association *P* value. The adipose *TSHR* instruments (rs12101255, rs74064796; 50 bp apart in *TSHR* intron 1) are themselves blood *TSHR* eQTLs near rs179252 and re-express the regional association rather than providing independent evidence. Cultured fibroblasts and adipose tissue are not orbital tissue. Direction, compared with the blood estimate in the same dataset. CI, confidence interval; GTEx, Genotype-Tissue Expression; OR, odds ratio; BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank.
