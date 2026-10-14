# Main tables

## Table 1. Data sources and their roles

| Dataset / endpoint | Cases | Controls | Role |
|----------------------------|--------|---------|------------------------------------------------|
| eQTLGen blood cis-eQTL | — | — | Exposure; up to 31,684 participants (median 28,092 per variant); predominantly European |
| BBJ GD (GCST90018627) | 2,809 | 172,656 | Discovery; East Asian |
| FinnGen GD (E4_GRAVES_STRICT) | 3,962 | 496,386 | Post hoc same-disease comparison; Finnish |
| UKB hyperthyroidism (GCST90038636) | 3,731 | 480,867 | Broader-phenotype comparison; self-reported hyperthyroidism or thyrotoxicosis; ancestry not established as exclusively European |
| FinnGen GO (GRAVES_OPHT) | 858 | 499,490 | Original ophthalmopathy (TED-enriched) comparison; Finnish |
| FinnGen GO strict (E4_GRAVES_OPHT_STRICT) | 753 | 499,595 | Post hoc ophthalmopathy definition sensitivity; Finnish |
| GTEx v8 tissue eQTL | — | — | Post hoc tissue extension (thyroid, cultured fibroblasts, subcutaneous and visceral adipose, whole blood; n = 469–670; Table S4) |
| 1000 Genomes Phase 3 European reference (n = 503) | — | — | LD clumping of instruments and reference allele frequencies |

The three FinnGen R12 endpoints share participants and are not independent replications. FinnGen ophthalmopathy controls are all other participants, including people with GD but no recorded eye disease, so GO estimates do not estimate TED risk among patients with GD. UKB hyperthyroidism is broader than GD. 'Strict' is the FinnGen endpoint label, not proof of perfect specificity. The descriptive orbital RNA-sequencing series (four TED patients, one control) is separate from these datasets. eQTL, expression quantitative trait locus; GTEx, Genotype-Tissue Expression; LD, linkage disequilibrium; TED, thyroid eye disease; BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank.

## Table 2. MR and colocalization results for *TSHR*, *IGF1R* and *CTLA4* across five disease datasets

| Gene | Dataset | IVs | OR (95% CI) | *P* value | Q-test *P* | PP.H2 | PP.H3 | PP.H4 | PP.H4, p12 = 10^−6^ |
|-------|-----------------|----|---------------|-------------|---------|-------|-------|-------|-----------|
| *TSHR* | BBJ GD | 1 | 0.12 (0.07–0.21) | 1.09×10^−14^\*† | NA | <0.001 | 0.049 | **0.951** | 0.661 |
| *TSHR* | FinnGen GD | 1 | 0.08 (0.05–0.12) | 3.19×10^−32^\*† | NA | <0.001 | 0.736 | 0.264 | 0.035 |
| *TSHR* | UKB hyperthyroidism | 1 | 0.09 (0.06–0.14) | 8.77×10^−28^\*† | NA | <0.001 | 0.774 | 0.226 | 0.028 |
| *TSHR* | FinnGen GO | 1 | 0.10 (0.04–0.24) | 2.82×10^−7^\*† | NA | <0.001 | 0.014 | **0.986** | **0.875** |
| *TSHR* | FinnGen GO strict | 1 | 0.15 (0.06–0.39) | 9.02×10^−5^\* | NA | 0.027 | 0.082 | **0.891** | 0.451 |
| *IGF1R* | BBJ GD | 4 | 1.56 (1.07–2.28) | 0.0212\* | 0.902 | 0.690 | 0.236 | 0.073 | 0.008 |
| *IGF1R* | FinnGen GD | 3 | 1.07 (0.86–1.33) | 0.539 | 0.882 | 0.731 | 0.259 | 0.010 | 0.001 |
| *IGF1R* | UKB hyperthyroidism | 4 | 1.35 (1.07–1.70) | 0.0117\* | 0.250 | 0.400 | 0.196 | 0.404 | 0.063 |
| *IGF1R* | FinnGen GO | 3 | 1.41 (0.85–2.33) | 0.182 | 0.309 | 0.623 | 0.346 | 0.032 | 0.003 |
| *IGF1R* | FinnGen GO strict | 3 | 1.52 (0.92–2.49) | 0.0996 | 0.409 | 0.571 | 0.391 | 0.038 | 0.004 |
| *CTLA4* | BBJ GD | 1 | 0.18 (0.11–0.27) | 5.45×10^−15^\*† | NA | <0.001 | 0.799 | 0.201 | 0.025 |
| *CTLA4* | FinnGen GD | 2 | 0.17 (0.08–0.36) | 2.15×10^−6^\*† | 0.0275 | <0.001 | 0.746 | 0.254 | 0.033 |
| *CTLA4* | UKB hyperthyroidism | 2 | 0.21 (0.08–0.57) | 0.00220\* | 0.00131 | <0.001 | 0.047 | **0.953** | 0.672 |
| *CTLA4* | FinnGen GO | 2 | 0.17 (0.04–0.68) | 0.0118\* | 0.0486 | <0.001 | 0.022 | **0.978** | **0.816** |
| *CTLA4* | FinnGen GO strict | 2 | 0.17 (0.05–0.64) | 0.00888\* | 0.0755 | <0.001 | 0.021 | **0.979** | **0.820** |

MR estimates are ORs per approximately one SD higher genetically predicted blood expression: Wald ratio for one IV, multiplicative random-effects IVW for several. *TSHR* used one instrument, rs179252, in every dataset; its G allele raises blood expression by 0.105 SD per allele (F = 176), and its per-allele disease ORs were 0.77–0.82 (Table S3). The per-SD ORs therefore rescale a moderate per-allele association, and each single-IV *P* value equals the variant's disease-association z-test *P* value. *CTLA4* multi-IV estimates are concentrated in rs13030124, and three of the four have Q-test *P* < 0.05; in BBJ and UKB, two of four single-instrument omissions removed nominal *IGF1R* significance (Figure S2). \* Nominal *P* < 0.05; † *P* < 0.05/2,544, the screen threshold, applied descriptively outside BBJ. PP.H2–H4 are colocalization posterior probabilities (H2, expression signal only; H3, distinct variants; H4, one shared variant) with p1 = p2 = 10^−4^ and p12 = 10^−5^ unless stated; bold, ≥0.80. They are not *P* values. BBJ colocalization pairs European eQTL with East Asian disease statistics. FinnGen ophthalmopathy controls are all other participants, including people with GD but no recorded eye disease, so GO estimates do not estimate TED risk among patients with GD. CI, confidence interval; eQTL, expression quantitative trait locus; IV, instrumental variable; IVW, inverse-variance weighted; MR, Mendelian randomization; NA, not applicable (one IV); OR, odds ratio; PP, posterior probability; SD, standard deviation; TED, thyroid eye disease; BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank.

## Table 3. Thirteen BBJ discoveries and their FinnGen GD comparison

| BBJ discovery gene | Chr | IVs, BBJ/FinnGen GD | BBJ OR (95% CI) | BBJ *P* | FinnGen GD OR (95% CI) | FinnGen GD *P* | Focused support | FinnGen GD PP.H4 |
|----------|------------|-----------|---------------|----------|----------------|-----------|--------|---------|
| *HLA-A* | 6 (MHC) | 4/4 | 1.98 (1.73–2.27) | 2.60×10^−23^ | 1.08 (0.77–1.53) | 0.648 | No | <0.001 |
| *HLA-DQA2* | 6 (MHC) | 1/1 | 2.39 (1.95–2.93) | 8.30×10^−17^ | 1.90 (1.64–2.20) | 1.41×10^−17^ | Yes | <0.001 |
| *CTLA4* | 2 | 1/2 | 0.18 (0.11–0.27) | 5.45×10^−15^ | 0.17 (0.08–0.36) | 2.15×10^−6^ | Yes | 0.254 |
| *TSHR* | 14 | 1/1 | 0.12 (0.07–0.21) | 1.09×10^−14^ | 0.08 (0.05–0.12) | 3.19×10^−32^ | Yes | 0.264 |
| *C4A* | 6 (MHC) | 1/3 | 0.45 (0.35–0.57) | 2.23×10^−10^ | 0.41 (0.26–0.63) | 4.97×10^−5^ | Yes | <0.001 |
| *HSD3B7* | 16 (16p11.2 cluster) | 1/2 | 0.29 (0.18–0.47) | 2.04×10^−7^ | 0.90 (0.72–1.13) | 0.372 | No | 0.011 |
| *TUBB* | 6 (MHC) | 2/3 | 1.49 (1.28–1.74) | 3.14×10^−7^ | 3.28 (2.09–5.17) | 2.68×10^−7^ | Yes | <0.001 |
| *VKORC1* | 16 (16p11.2 cluster) | 1/0 | 0.13 (0.06–0.29) | 6.35×10^−7^ | NA | NA | NA | 0.011 |
| *TNFSF14* | 19 | 1/1 | 0.63 (0.52–0.76) | 1.48×10^−6^ | 0.77 (0.65–0.92) | 0.00361 | Yes | 0.303 |
| *PRSS36* | 16 (16p11.2 cluster) | 1/1 | 6.00 (2.87–12.57) | 2.02×10^−6^ | 1.40 (1.00–1.96) | 0.0490 | No | 0.063 |
| *MAPKAPK5* | 12 | 4/4 | 0.977 (0.967–0.987) | 5.28×10^−6^ | 0.990 (0.929–1.055) | 0.751 | No | 0.041 |
| *PSMB8* | 6 (MHC) | 1/1 | 1.25 (1.13–1.38) | 6.86×10^−6^ | 0.96 (0.85–1.08) | 0.512 | No | <0.001 |
| *IFNGR1* | 6 | 1/1 | 2.10 (1.51–2.91) | 9.41×10^−6^ | 1.62 (1.18–2.24) | 0.00330 | Yes | 0.313 |

Genes passed the BBJ screen threshold (*P* < 0.05/2,544) and are listed by BBJ *P* value. Focused support in FinnGen GD requires the same direction and *P* < 0.05/13; 7 of 13 qualify. IVs are counted separately for each dataset because instrument availability differs. The five MHC genes and the three chromosome 16p11.2 genes are not independent signals; the latter are compatible with one regional association. *MAPKAPK5* is shown to three decimals to preserve its small effect. NA, not estimable (no retained instrument), not a negative result; *VKORC1* colocalization does not require an MR instrument. PP.H4, shared-variant posterior probability under p12 = 10^−5^; bold, ≥0.80. Chr, chromosome of the strongest instrument; MHC, major histocompatibility complex; CI, confidence interval; IV, instrumental variable; MR, Mendelian randomization; OR, odds ratio; BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank.
