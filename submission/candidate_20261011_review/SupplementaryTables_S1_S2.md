## Table S1. *TSHR* regional and prior sensitivity

| Dataset | Full cis SNPs | Primary PP.H4 | Lower-prior PP.H4 | PP.H4 range across windows |
|---|---|---|---|---|
| BBJ GD | 4,336 | 0.951 | 0.661 | 0.951–0.951 |
| FinnGen GD | 6,579 | 0.264 | 0.035 | 0.264–0.264 |
| UKB HT | 5,667 | 0.226 | 0.028 | 0.226–0.226 |
| FinnGen GO | 6,579 | 0.986 | 0.875 | 0.986–0.987 |
| FinnGen GO strict | 6,579 | 0.891 | 0.451 | 0.891–0.905 |

Lower prior: p12 = 10⁻⁶; p1 = p2 = 10⁻⁴ throughout. Window range uses the full source cis-region, ±500 kb, ±250 kb and five-outcome intersection at the primary prior. It is not a confidence interval. Full regional disease lead variants were retained. All 240 settings are in TSHR_windows.

## Table S2. Frequency-sensitive exploratory candidates

| Gene / dataset | Frequency | IVs | OR (95% CI) | *P* value | Q-test *P* |
|---|---|---|---|---|---|
| *C4B* / FinnGen GD | Reference | 2 | 1.63 (1.31–2.03) | 1.10×10⁻⁵ | 7.79×10⁻⁷ |
| *C4B* / FinnGen GD | eQTLGen | 2 | 1.75 (1.31–2.33) | 1.48×10⁻⁴ | 1.36×10⁻⁸ |
| *C4B* / FinnGen GO strict | Reference | 2 | 1.68 (1.44–1.95) | 2.74×10⁻¹¹ | 0.122 |
| *C4B* / FinnGen GO strict | eQTLGen | 2 | 1.80 (1.46–2.23) | 6.17×10⁻⁸ | 0.0587 |
| *HLA-G* / FinnGen GD | Reference | 16 | 0.86 (0.72–1.03) | 0.107 | 9.98×10⁻²⁴ |
| *HLA-G* / FinnGen GD | eQTLGen | 16 | 0.60 (0.47–0.75) | 1.72×10⁻⁵ | 2.26×10⁻¹⁰ |
| *HLA-G* / FinnGen GO strict | Reference | 16 | 0.91 (0.76–1.08) | 0.267 | 0.0231 |
| *HLA-G* / FinnGen GO strict | eQTLGen | 16 | 0.70 (0.53–0.92) | 0.0122 | 0.131 |

*C4B* and *HLA-G* frequency sensitivity. Variant sets and outcome effects were held fixed. Frequency changes rescaled expression estimates and IVW weights. Q tests heterogeneity, not directional pleiotropy; its *P* value is not a gene-association *P* value. All 72 omission estimates and instrument weights are provided in the workbook.
