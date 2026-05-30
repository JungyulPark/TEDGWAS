# TED-TRAP v5 — Main Figure & Table Captions (LOCKED text)

These captions go into the manuscript. They encode the reviewer-proofing language
(exploratory tissue support, denominator clarity, n=1 control disclosure).

---

## Table 1. Study design and data sources
*(footnote already in table)*

eQTLGen sample size indicates the maximum available blood cis-eQTL sample size;
SNP-level sample sizes vary across instruments (e.g., rs179252 NrSamples = 31,566).
BBJ Graves disease was used as the primary discovery phenotype, UKB hyperthyroidism
as cross-ancestry/broader-phenotype replication, and FinnGen Graves ophthalmopathy
as the TED-specific sensitivity outcome. UKB binary-trait effect estimates were
rescaled from linear scale to log-odds using case prevalence.

---

## Table 2. Integrated genetic and tissue evidence for backbone genes
*(footnote already in table)*

MR by inverse-variance weighted (multi-IV) or Wald ratio (single-IV). TSHR was
single-IV across all outcomes; CTLA4 included single- or limited-IV estimates
depending on outcome (BBJ single-IV; UKB/FinnGen 2-IV), limiting pleiotropy
diagnostics. Coloc by coloc.abf (eQTL quant, GWAS cc); PP.H4 > 0.8 = shared causal
variant. IGF1R coloc PP.H2-dominant (GWAS signal without shared eQTL). Tissue by
DESeq2 (in-house orbital, n=1 control — exploratory supporting layer). All values
verified against source files.

---

## Figure 2. Druggable-gene-wide BBJ discovery Mendelian randomization

Volcano plot of cis-eQTL MR effects on Graves disease (BBJ discovery cohort,
log-odds scale) for 2,235 druggable genes with estimable BBJ MR effects, among
2,545 MR-testable genes drawn from the druggable genome (Finan et al. 2017)
with valid eQTLGen instruments. The horizontal dashed line marks the Bonferroni
significance threshold (P = 1.965 × 10⁻⁵), computed against all 2,545 MR-testable
genes. Thirteen genes exceed this threshold. TSHR and CTLA4 (red), both established
Graves disease loci, show strong protective-direction associations and serve as
internal anchors/positive control. HLA-A and HLA-DQA2 reflect MHC-region signal
(grey). Three genes within the chr16p11.2 LD block (HSD3B7, VKORC1, PRSS36;
light purple) represent a single linkage-disequilibrium signal rather than
independent associations (Supplementary Table S5). IGF1R (blue) lies below the
Bonferroni threshold at nominal significance, consistent with a pharmacologic
effector rather than a primary susceptibility locus. No druggable gene outside the
known anchors and LD/MHC artifacts achieved robust, colocalization-confirmed
novel association (Supplementary Tables S1–S2).

---

## Figure 3. Multi-layer evidence integration for backbone genes

Convergent genetic and tissue evidence for TSHR, IGF1R, and CTLA4.
**(A)** Bayesian colocalization posterior probabilities (coloc.abf) for each gene
against BBJ Graves disease and FinnGen Graves ophthalmopathy. Bars show PP.H2
(GWAS association only), PP.H3 (distinct causal variants), and PP.H4 (shared causal
variant); the dashed line marks PP.H4 = 0.8. TSHR shows strong colocalization in
both outcomes; IGF1R is PP.H2-dominant (outcome association without a shared eQTL
signal); CTLA4 shows a mixed pattern (PP.H3-dominant in BBJ, PP.H4 in FinnGen),
consistent with locus complexity.
**(B)** MR effect estimates (β, log-odds) with 95% confidence intervals across the
three outcomes. TSHR and CTLA4 show consistent protective-direction effects;
IGF1R shows a consistent nominal risk-direction effect.
**(C)** Orbital tissue expression (TPM) from in-house RNA-seq; TED n = 4 and
control n = 1. DESeq2 statistics (log2 fold-change, adjusted P) were computed at
the biological-sample level after collapsing technical-replicate libraries, and
should be interpreted as exploratory tissue support given the single control. TPM
panels use gene-specific y-axis scales. Only TSHR reaches tissue-level significance
(log2FC = +2.33, padj = 0.032); IGF1R and CTLA4 are directionally up but
non-significant.

---

## Methods sentences (locked, for Methods section)

**Tissue / pseudoreplication:**
"Each of four TED patients and one control was sequenced as two technical-replicate
libraries. Read counts were summed across technical replicates to the
biological-sample level (4 TED, 1 control) prior to DESeq2 analysis, avoiding
pseudoreplication. Given the single control, tissue results are interpreted as
exploratory support rather than independent confirmation."

**Phenotype roles:**
"BBJ Graves disease was used as the primary discovery phenotype to maximize
autoimmune genetic power; UKB hyperthyroidism provided cross-ancestry,
broader-phenotype replication; and FinnGen Graves ophthalmopathy served as the
TED-specific sensitivity outcome."

**Bonferroni:**
"The druggable-gene-wide significance threshold was set at P = 0.05 / 2,545 =
1.965 × 10⁻⁵, corresponding to the number of MR-testable druggable genes with
valid instruments."

---

## SUPPLEMENTARY FIGURE CAPTIONS

### Figure S1. TSHR fine-mapping and LD rationale for single-instrument selection

SuSiE fine-mapping of the TSHR cis-eQTL signal (eQTLGen) using 1000 Genomes Phase 3
European and East Asian LD reference panels. **(A)** Genomic positions of credible-set
lead SNPs at the TSHR locus (chr14, hg19); rs179252 (star) is the primary instrument used
in all Mendelian randomization analyses. All credible-set leads fall within a ~15 kb window.
**(B)** Ancestry-matched LD (r²) of each secondary credible-set lead with rs179252. Every
secondary lead remains in substantial LD with rs179252 in its corresponding ancestry-matched
panel (minimum r² = 0.808 in East Asian, 0.965 in European), exceeding the r² = 0.8 threshold
(dashed line). Because no secondary signal is independent of rs179252, TSHR is retained as a
single-instrument / locus-level anchor rather than a multi-instrument exposure; formal
multi-instrument pleiotropy diagnostics are therefore not applicable at this locus. This figure
summarizes credible-set lead variants and their LD with the retained anchor SNP, rather than
showing a full regional association plot.

### Figure S2. Colocalization probability profiles for candidate genes

Stacked bars represent the colocalization posterior probabilities PP.H2 (GWAS association only, orange), PP.H3 (distinct causal variants, grey), and PP.H4 (shared causal variant, green) for six druggable candidate genes across BBJ Graves' disease and FinnGen Graves' ophthalmopathy. The horizontal dashed line marks PP.H4 = 0.8. Although candidate targets such as TNFSF14 and IFNGR1 show strong colocalization (PP.H4 > 0.8) in the BBJ discovery analysis, they fail to colocalize in FinnGen GO (PP.H4 <= 0.038), resulting in zero robust, colocalization-confirmed novel targets. No candidate showed strong colocalization (PP.H4 > 0.8) in both BBJ Graves' disease and FinnGen Graves' ophthalmopathy. Full H0–H4 posterior probabilities are provided in Table S2.


