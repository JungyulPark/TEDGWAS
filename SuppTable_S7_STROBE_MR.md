# Supplementary Table S7 — STROBE-MR Checklist (TED-TRAP v5)

STROBE-MR (Skrivankova et al., JAMA 2021) reporting checklist for Mendelian
randomization studies. Each item maps to where it is addressed in the manuscript.
*(Page/section references to be filled once manuscript is paginated.)*

| # | STROBE-MR Item | Addressed in | Status |
|---|----------------|--------------|--------|
| **Title & Abstract** | | | |
| 1 | MR stated in title/abstract; design described | Title, Abstract | ✓ |
| **Introduction** | | | |
| 2 | Background; rationale for using genetic variants as instruments | Introduction | ✓ |
| 3 | Explicit hypotheses / objectives | Introduction (end) | ✓ |
| **Methods** | | | |
| 4 | Study design & key elements | Methods; Table 1; Figure 1 | ✓ |
| 5 | Data sources, settings, populations (exposure & outcome) | Methods; Table 1 | ✓ |
| 6a | Instrument selection criteria | Methods (eQTLGen; Finan druggable genome; cis P < 5×10⁻⁸) | ✓ |
| 6b | Assumptions (relevance, independence, exclusion restriction) | Methods; Discussion (limitations) | ✓ |
| 7 | Measurement / definition of variables (exposure, outcome, covariates) | Methods | ✓ |
| 8 | Effect measures & scale (log-odds; UKB rescaling) | Methods; Table 1 footnote | ✓ |
| 9 | Statistical methods (IVW, MR-Egger, WM, Wmode; coloc.abf) | Methods | ✓ |
| 10 | Sensitivity analyses (pleiotropy, heterogeneity, Steiger) | Methods; Supp Table S6 | ✓ |
| 11 | Software & packages (versions listed below) | Methods | ✓ |
| **Results** | | | |
| 12 | Descriptive data (instrument counts, F-statistics) | Results; Supp Table S4 | ✓ |
| 13 | Main MR estimates with CIs | Results; Table 2; Table 3; Figure 2/3 | ✓ |
| 14 | Sensitivity / robustness results | Results; Supp Table S6 | ✓ |
| 15 | Colocalization results | Results; Table 2; Figure 3A; Supp Tables S2/S5 | ✓ |
| 16 | Additional analyses (tissue expression) | Results; Figure 3C | ✓ |
| **Discussion** | | | |
| 17 | Key results in context | Discussion | ✓ |
| 18 | Limitations (cross-ancestry, single-IV, blood eQTL, n=1 control) | Discussion (limitations) | ✓ |
| 19 | Interpretation (susceptibility-anchor vs effector framing) | Discussion | ✓ |
| 20 | Generalizability | Discussion | ✓ |
| **Other** | | | |
| 21 | Funding sources | Funding statement | ✓ |
| 22 | Data availability | Data availability (public summary statistics: eQTLGen, GWAS Catalog, FinnGen) | ✓ |

## Software and package versions

R 4.3.3; TwoSampleMR v0.7.4; ieugwasr; coloc v5.2.3; susieR v0.14.2;
DESeq2 v1.42.1; PLINK v1.9.

## Key MR-assumption documentation (for reviewer reference)

**IV1 — Relevance:** All instruments are genome-wide significant cis-eQTLs (eQTLGen,
P < 5×10⁻⁸). F-statistics computed per instrument; across all 2,545 MR-testable genes
(6,136 SNP-level instruments) none were weak (all F > 10, minimum F ≈ 14.2), indicating
negligible weak-instrument bias (Supplementary Table S4).

**IV2 — Independence:** cis-eQTL instruments restricted to gene-local windows and
LD-clumped. Confounding by population structure mitigated by ancestry-matched LD
references where possible; the cross-ancestry relationship between predominantly
European eQTLGen instruments and the East Asian BBJ discovery outcome is acknowledged
in the Discussion.

**IV3 — Exclusion restriction:** For multi-instrument genes (IGF1R; CTLA4 in UKB/FinnGen),
the MR-Egger intercept was used to test for directional pleiotropy (IGF1R intercept
non-significant in all outcomes; Supplementary Table S6). For single-instrument genes
(TSHR in all outcomes; CTLA4 in BBJ), pleiotropy diagnostics are not estimable;
colocalization (coloc.abf) provides orthogonal evidence that the eQTL and GWAS signals
share a causal variant (TSHR PP.H4 = 0.951 in BBJ and 0.986 in FinnGen), supporting the
exclusion restriction for the key result.

**Directionality:** Exposure (cis-eQTL on gene expression) precedes outcome (disease
liability); reverse causation is implausible at the cis-eQTL level. Steiger filtering was
applied where estimable (Supplementary Table S6).
