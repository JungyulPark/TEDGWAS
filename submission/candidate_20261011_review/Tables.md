# Main tables

## Table 1. Data sources

| Dataset / endpoint | Cases | Controls | Exposure or clinical role |
|---|---|---|---|
| eQTLGen blood expression | — | — | Exposure: up to 31,684; predominantly European |
| BBJ GD (GCST90018627) | 2,809 | 172,656 | Discovery; East Asian |
| FinnGen GD (E4_GRAVES_STRICT) | 3,962 | 496,386 | Post hoc GD comparison; Finnish |
| UKB hyperthyroidism (GCST90038636) | 3,731 | 480,867 | Self-report; ancestry not established as exclusively European |
| FinnGen GO (GRAVES_OPHT) | 858 | 499,490 | Original ophthalmopathy comparison; Finnish |
| FinnGen GO strict (E4_GRAVES_OPHT_STRICT) | 753 | 499,595 | Post hoc definition sensitivity; Finnish |

GD, Graves disease; GO, Graves ophthalmopathy; HT, hyperthyroidism/thyrotoxicosis. The three FinnGen R12 endpoints overlap and use population non-cases for their respective definitions. Strict is the source endpoint label, not proof of perfect specificity. The descriptive orbital series (four TED patients, one control) is separate from the genetic datasets.

## Table 2. TSHR, IGF1R and CTLA4

| Gene | Dataset | IVs | OR (95% CI) | P value | PP.H3 | PP.H4 |
|---|---|---|---|---|---|---|
| TSHR | BBJ GD | 1 | 0.12 (0.07–0.21) | 1.09×10⁻¹⁴ | 0.049 | 0.951 |
| TSHR | FinnGen GD | 1 | 0.08 (0.05–0.12) | 3.19×10⁻³² | 0.736 | 0.264 |
| TSHR | UKB HT | 1 | 0.09 (0.06–0.14) | 8.77×10⁻²⁸ | 0.774 | 0.226 |
| TSHR | FinnGen GO | 1 | 0.10 (0.04–0.24) | 2.82×10⁻⁷ | 0.014 | 0.986 |
| TSHR | FinnGen GO strict | 1 | 0.15 (0.06–0.39) | 9.02×10⁻⁵ | 0.082 | 0.891 |
| IGF1R | BBJ GD | 4 | 1.56 (1.07–2.28) | 0.0212 | 0.236 | 0.073 |
| IGF1R | FinnGen GD | 3 | 1.07 (0.86–1.33) | 0.539 | 0.259 | 0.010 |
| IGF1R | UKB HT | 4 | 1.35 (1.07–1.70) | 0.0117 | 0.196 | 0.404 |
| IGF1R | FinnGen GO | 3 | 1.41 (0.85–2.33) | 0.182 | 0.346 | 0.032 |
| IGF1R | FinnGen GO strict | 3 | 1.52 (0.92–2.49) | 0.0996 | 0.391 | 0.038 |
| CTLA4 | BBJ GD | 1 | 0.18 (0.11–0.27) | 5.45×10⁻¹⁵ | 0.799 | 0.201 |
| CTLA4 | FinnGen GD | 2 | 0.17 (0.08–0.36) | 2.15×10⁻⁶ | 0.746 | 0.254 |
| CTLA4 | UKB HT | 2 | 0.21 (0.08–0.57) | 0.0022 | 0.047 | 0.953 |
| CTLA4 | FinnGen GO | 2 | 0.17 (0.04–0.68) | 0.0118 | 0.022 | 0.978 |
| CTLA4 | FinnGen GO strict | 2 | 0.17 (0.05–0.64) | 0.00888 | 0.021 | 0.979 |

IVs, instrumental variables (retained genetic variants). Wald ratio for one IV; random-effects IVW for multiple IVs. ORs are per approximately one SD higher reconstructed blood expression. Bold P values denote nominal P < 0.05 only. Discovery uses 0.05/2,544. PP.H3 supports distinct variants; PP.H4 supports one shared variant. Bold PP.H4 denotes ≥0.80 under p12 = 10⁻⁵. These are probabilities, not P values. CTLA4 multi-IV results are concentrated in rs13030124; TSHR is single-IV in all datasets.

## Table 3. Screen and supporting analyses

### A. BBJ discoveries and FinnGen GD comparison

| BBJ discovery gene | BBJ OR (95% CI) | BBJ P | FinnGen GD OR (95% CI) | FinnGen GD P | Focused support | GD PP.H4 |
|---|---|---|---|---|---|---|
| HLA-A | 1.98 (1.73–2.27) | 2.60×10⁻²³ | 1.08 (0.77–1.53) | 0.648 | No | <0.001 |
| HLA-DQA2 | 2.39 (1.95–2.93) | 8.30×10⁻¹⁷ | 1.90 (1.64–2.20) | 1.41×10⁻¹⁷ | Yes | <0.001 |
| CTLA4 | 0.18 (0.11–0.27) | 5.45×10⁻¹⁵ | 0.17 (0.08–0.36) | 2.15×10⁻⁶ | Yes | 0.254 |
| TSHR | 0.12 (0.07–0.21) | 1.09×10⁻¹⁴ | 0.08 (0.05–0.12) | 3.19×10⁻³² | Yes | 0.264 |
| C4A | 0.45 (0.35–0.57) | 2.23×10⁻¹⁰ | 0.41 (0.26–0.63) | 4.97×10⁻⁵ | Yes | <0.001 |
| HSD3B7 | 0.29 (0.18–0.47) | 2.04×10⁻⁷ | 0.90 (0.72–1.13) | 0.372 | No | 0.011 |
| TUBB | 1.49 (1.28–1.74) | 3.14×10⁻⁷ | 3.28 (2.09–5.17) | 2.68×10⁻⁷ | Yes | <0.001 |
| VKORC1 | 0.13 (0.06–0.29) | 6.35×10⁻⁷ | NA | NA | NA | 0.011 |
| TNFSF14 | 0.63 (0.52–0.76) | 1.48×10⁻⁶ | 0.77 (0.65–0.92) | 0.00361 | Yes | 0.303 |
| PRSS36 | 6.00 (2.87–12.57) | 2.02×10⁻⁶ | 1.40 (1.00–1.96) | 0.049 | No | 0.063 |
| MAPKAPK5 | 0.977 (0.967–0.987) | 5.28×10⁻⁶ | 0.990 (0.929–1.055) | 0.751 | No | 0.041 |
| PSMB8 | 1.25 (1.13–1.38) | 6.86×10⁻⁶ | 0.96 (0.85–1.08) | 0.512 | No | <0.001 |
| IFNGR1 | 2.10 (1.51–2.91) | 9.41×10⁻⁶ | 1.62 (1.18–2.24) | 0.0033 | Yes | 0.313 |

### B. Additional FinnGen GD discoveries

| Additional gene | GD IVs | GD OR (95% CI) | GD P | BBJ PP.H4 | GD PP.H4 | GO strict PP.H4 |
|---|---|---|---|---|---|---|
| APOM | 1 | 0.05 (0.04–0.06) | 1.43×10⁻¹⁰⁸ | <0.001 | <0.001 | 0.227 |
| NEU1 | 1 | 0.04 (0.03–0.06) | 1.43×10⁻¹⁰⁸ | <0.001 | <0.001 | 0.061 |
| IGLV2-18 | 1 | 0.56 (0.48–0.66) | 7.52×10⁻¹² | 0.029 | 0.020 | 0.100 |
| IGLV3-12 | 1 | 1.20 (1.14–1.26) | 7.52×10⁻¹² | 0.019 | 0.022 | 0.101 |
| TNF | 1 | 0.60 (0.51–0.71) | 4.65×10⁻¹⁰ | <0.001 | <0.001 | <0.001 |
| HIST1H4E | 1 | 2.54 (1.88–3.43) | 1.51×10⁻⁹ | <0.001 | <0.001 | 0.007 |
| LTA | 1 | 0.46 (0.35–0.61) | 3.75×10⁻⁸ | <0.001 | <0.001 | <0.001 |
| LTB | 1 | 0.18 (0.10–0.36) | 4.89×10⁻⁷ | <0.001 | <0.001 | <0.001 |
| C4B | 2 | 1.63 (1.31–2.03) | 1.10×10⁻⁵ | <0.001 | <0.001 | 0.311 |

### C. Frequency-sensitive results

| Gene / GD analysis | Reference OR (95% CI) | Reference P | Cohort-frequency OR (95% CI) | Cohort P | Interpretation |
|---|---|---|---|---|---|
| C4B | 1.63 (1.31–2.03) | 1.10×10⁻⁵ | 1.75 (1.31–2.33) | 1.48×10⁻⁴ | Heterogeneous; loses screen significance |
| HLA-G | 0.86 (0.72–1.03) | 0.107 | 0.60 (0.47–0.75) | 1.72×10⁻⁵ | Weight-sensitive; single-variant dependence |

Panel A: focused support requires concordant direction and P < 0.05/13; 7 of 13 qualify. NA is unavailable, not a negative result. Panel B: nine additional reference-frequency FinnGen GD discoveries at P < 0.05/2,544; post hoc and not independent of the full screen. APOM/NEU1 share rs3101018; IGLV2-18/IGLV3-12 share rs5759410. The lead IGLV expression variant is absent in BBJ, weakening interpretation of its low PP.H4. Panel C: HLA-G was added because of frequency sensitivity, not primary discovery. Neither C4B nor HLA-G gains strong GD colocalization. All 24 genes and five datasets are available in the workbook.
In the formatted DOCX/PDF, bold P values denote nominal P < 0.05; screen and focused thresholds are stated separately.
