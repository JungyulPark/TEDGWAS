# TED-TRAP Main Manuscript — Revised 5-Element Plan (v2)

**Target Journal**: Thyroid (Mary Ann Liebert)
**Constraint**: ≤ 5 figures + tables total (Main), unlimited Supplementary
**Word limit**: 4,000 words
**Abstract**: 200 words, unstructured

## User-Specified Adjustments (v2)
- Tables reduced to 2-3 (easier to handle than earlier 1-table plan); Figures reduced accordingly
- Figures will be created via **Nanobananapro** or **FigureLab** — need explicit natural-language prompts
- **No Cytoscape** (unavailable); pathway networks must be renderable by AI/diagram tools
- RNA-seq cohort: "2024 in-house cohort, IRB# placeholder" (no "submitted")

---

## Main Manuscript: Final 5-element Layout (v2)

| # | Type | Element | Purpose |
|---|------|---------|---------|
| 1 | **Figure 1** | Study design schematic | Workflow overview (AI-generated) |
| 2 | **Table 1** | Study design + data sources | Clean reference table |
| 3 | **Table 2** | MR results + Coloc summary (combined) | Genetic causal evidence |
| 4 | **Figure 2** | Differential pathway network | 3-zone drug coverage + KEGG |
| 5 | **Figure 3** | Triangulation + transcriptomic + off-target (composite) | Integrative central visual |

**Rationale for this re-balance**:
- 2 Tables + 3 Figures = easier proofing than 1 Table + 4 Figures
- Table 1 absorbs Study Design (freed Fig 1 to be purely schematic, cleaner)
- Table 2 combines MR + Coloc (condensed)
- Fig 3 is a composite that packages Triangulation + RNA-seq + off-target in one image

---

## FIGURE 1 — Study Design Schematic

**Purpose**: Single-image overview of the pipeline. First thing reviewers see.

### Prompt for Nanobananapro / FigureLab

```
Create a scientific workflow schematic for a computational biology study
on thyroid eye disease. Use a clean academic figure style similar to
Nature Methods or Thyroid journal illustrations. Portrait orientation,
approximately 7.5 × 10 inches.

The figure has 4 horizontal bands flowing top to bottom with connecting
arrows between bands.

=== BAND 1 (TOP): INPUT DATA — 3 vertical boxes side by side ===

Box 1 (left, blue theme #2E75B6):
Title: "Genetic Data"
Content:
  • eQTLGen cis-eQTL (n = 31,684, European)
  • GWAS Primary: Graves disease (ebi-a-GCST90018627; 2,809 cases)
  • GWAS Replication: Hyperthyroidism (ebi-a-GCST90038636; 3,731 cases)
  • GWAS Sensitivity: FinnGen R12 Graves ophthalmopathy (858 cases)
Icon: DNA double helix

Box 2 (center, orange theme #ED7D31):
Title: "Pathway Gene Sets"
Content:
  • IGF-1R pathway (n = 36 genes)
  • TSHR pathway (n = 33 genes)
  • TED disease genes (n = 59 genes)
  • Sources: CTD, STRING, DrugBank, literature
Icon: Network nodes

Box 3 (right, purple theme #7030A0):
Title: "Transcriptomic Data"
Content:
  • In-house bulk RNA-seq (2024)
  • 4 inactive TED + 1 control
  • Korean cohort, orbital adipose tissue
  • IRB# [INSERT]
Icon: RNA strand or orbital anatomy silhouette

=== BAND 2: ANALYTICAL METHODS — 4 parallel streams ===

From left to right, 4 rounded rectangles with labels:
  1. "Mendelian Randomization" (TwoSampleMR, IVW + sensitivity)
  2. "Bayesian Colocalization" (coloc R package)
  3. "Differential Pathway Network" (STRING PPI + KEGG hsa04024)
  4. "Off-target Pathway" (KEGG hsa04151 + hsa04910)

Arrows flow from Band 1 into these methods boxes.

=== BAND 3: INTERMEDIATE OUTPUTS ===

4 small result boxes aligned under each method:
  1. "TSHR genetic causality (P = 6.79e-17)"
  2. "Shared causal variant (PP.H4 = 0.985)"
  3. "15 intersection genes (4/3/8 zones)"
  4. "38 cochlear + 51 insulin shared effectors"

=== BAND 4 (BOTTOM): INTEGRATION ===

One large golden box (#BF9000) spanning full width:
Title: "Multi-modal Triangulation"
Content: "TSHR as convergent therapeutic target in TED — implications
for teprotumumab relapse and off-target effects"

=== STYLING ===
- Use flat, clean vector style (no 3D shading)
- Sans-serif font (Arial or Helvetica family) throughout
- Arrows: thin black with small arrowheads
- Avoid clipart or overly illustrative icons; keep icons minimal and schematic
- Leave space around edges for publication margins
- No watermarks, no logos
```

**Output**: PNG at ≥ 300 DPI, transparent or white background.

**Legend text (~120 words)** for the manuscript:
```
Figure 1. Study design and analytical workflow. Three classes of input
data (genetic, pathway, and transcriptomic) converge through four
parallel analytical streams: Mendelian randomization (MR), Bayesian
colocalization, differential pathway network analysis, and off-target
pathway interrogation. Exposure instruments were obtained from the
eQTLGen Consortium (blood cis-eQTL, n = 31,684). Outcome data comprised
three European-ancestry GWAS: Graves disease as the primary outcome
(GCST90018627), hyperthyroidism as the replication outcome
(GCST90038636), and FinnGen Release 12 Graves ophthalmopathy as a
TED-specific sensitivity analysis. Transcriptomic validation used an
in-house bulk RNA-seq cohort (four inactive TED patients and one
control, orbital adipose tissue, IRB#[INSERT]). Multi-modal results
were integrated through a triangulation framework to identify
convergent therapeutic targets.
```

---

## TABLE 1 — Study design and data sources

**Already built** in prior step (TED_TRAP_Tables_1_2_v2.docx, Table 1).

No changes needed. 5 rows: Primary / Replication / Sensitivity / Exposure / Transcriptomic.

**Adjustment**: RNA-seq row year → "2024"; Dataset column → "In-house bulk RNA-seq"; add IRB# to Methods text (not in table itself).

---

## TABLE 2 — MR Results + Colocalization (Combined)

**Change from prior plan**: Merge MR summary with Coloc summary into one compact table. Full MR details go to Supplementary.

**Layout**: Two panels in one table structure

### Panel A: Mendelian Randomization (condensed)

| Gene | Module | Primary β (P) | Replication β (P) | Sensitivity β (P) |
|------|--------|---------------|-------------------|-------------------|
| **TSHR** | TSHR pathway | **−1.394 (6.79e−17)** | **−0.012 (3.88e−29)** | **−1.453 (6.81e−07)** |
| IGF1R | IGF-1R pathway | 0.217 (0.053) | **0.001 (0.041)** | 0.233 (0.216) |
| TNF | Inflammatory | 0.881 (0.429) | **−0.004 (3.57e−05)†** | −0.543 (0.111) |
| PPARG | Adipogenesis | −0.079 (0.406) | −1.75e−5 (0.970) | — |
| ARRB1 | β-arrestin scaffold | −0.140 (0.239) | −6.70e−4 (0.130) | — |
| IRS1 | IGF-1R pathway | −0.138 (0.121) | 4.01e−6 (0.993) | — |
| AKT1 | IGF-1R pathway | 0.093 (0.431) | −2.40e−4 (0.599) | — |
| CTLA4 | Immune regulation | −0.655 (0.070) | −3.36e−4 (0.931) | — |

Bold = P < 0.05. Primary/Sensitivity on log-odds scale; Replication on LMM scale. Insufficient IVs: IL6, TGFB1, HAS2, IGF1, IGF2.
† Cochran Q P = 0.028; suggestive finding.

### Panel B: Colocalization (TSHR locus only)

| Region | chr14:80.5–81.7 Mb (TSHR) |
|--------|---------------------------|
| eQTL dataset | eQTLGen (n = 31,684, quantitative) |
| GWAS dataset | FinnGen R12 Graves ophthalmopathy (n = 858 cases, case-control s = 0.0017) |
| Overlapping SNPs | 4,046 |
| **PP.H4 (shared causal)** | **0.985** |
| PP.H3 (distinct causal) | 0.015 |
| Top GWAS SNP | rs3783947 |
| Top eQTL SNP | rs179252 |
| LD r² (EUR 1000G) | 0.737 |

**Visual treatment**:
- TSHR row in Panel A shaded yellow
- Panel B appears immediately below Panel A with clear heading separation

**Footnote** (~80 words):
```
MR was performed using two-sample inverse-variance weighted analysis
(TwoSampleMR). Effect sizes are in log-odds for Primary/Sensitivity
outcomes and in linear mixed model units for the Replication outcome;
values are therefore not directly comparable across columns. TSHR
analyses used 2 IVs, precluding formal MR-Egger pleiotropy assessment.
Colocalization was performed with coloc (default priors); a posterior
probability of shared causal variant (PP.H4) > 0.80 was considered
strong evidence. Full sensitivity statistics in Supplementary Table S1.
```

---

## FIGURE 2 — Differential Pathway Network

**Change from prior plan**: No Cytoscape. Two options.

### Prompt for Nanobananapro / FigureLab

```
Create a two-panel scientific figure showing pathway analysis for
thyroid eye disease drug targets. Academic publication style, portrait
orientation, approximately 7.5 × 9 inches.

=== PANEL A (TOP): Three-set Venn Diagram ===

Three large overlapping circles with clear transparency so overlaps
are visible.

Circle 1 (left, blue #2E75B6 at 40% opacity):
Label above: "IGF-1R pathway (n = 36)"

Circle 2 (right, red #C00000 at 40% opacity):
Label above: "TSHR pathway (n = 33)"

Circle 3 (bottom, gray #808080 at 40% opacity):
Label below: "TED disease genes (n = 59)"

In the 3-way overlap region (center), write:
  "Shared zone (n = 8):
   HAS1, HAS2, HAS3, TNF, IL6, ARRB1, PPARG, CEBPA"

In the 2-way overlap of IGF-1R ∩ TED (excluding TSHR), write:
  "IGF1R-only zone (n = 4):
   IGF1R, IGF1, IL1B, CXCL8"

In the 2-way overlap of TSHR ∩ TED (excluding IGF-1R), write:
  "TSHR-only zone (n = 3):
   TSHR, ADIPOQ, FABP4"

=== PANEL B (BOTTOM): KEGG cAMP Signaling Pathway Schematic ===

Horizontal flow diagram showing two parallel signaling routes
converging on the same downstream module.

Top route (solid arrows, red theme — "canonical, direct"):
  TSHR → GNAS → Adenylyl cyclase → cAMP → PKA → CREB

Bottom route (dashed arrows, blue theme — "indirect, no canonical
KEGG path to cAMP"):
  IGF-1R → IRS1/2 → PI3K → AKT → (dashed line crossing over to CREB?)

Right-side annotation box:
  "Extension Layer:
   KEGG hsa04024 (cAMP pathway)
   TSHR annotated as canonical input;
   IGF-1R not present in hsa04024"

Label at bottom:
  "TSHR directly accesses the cAMP/PKA module;
   IGF-1R has no canonical KEGG pathway to this node."

=== STYLING ===
- Vector style, flat colors
- Sans-serif font, readable at print size
- Solid arrows for canonical KEGG edges
- Dashed arrows for indirect/absent canonical edges
- White background, no 3D effects, no extra decorations
```

**Output**: PNG ≥ 300 DPI.

**Legend** (~150 words):
```
Figure 2. Differential pathway coverage between IGF-1R and TSHR
signaling in TED. (A) Three-set Venn diagram showing overlap between
the IGF-1R pathway (n = 36 genes), TSHR pathway (n = 33), and the
curated TED disease gene set (n = 59). Fifteen genes lie in the
triple- or two-way intersections, distributed across three zones:
IGF1R-only (4 genes), TSHR-only (3 genes), and Shared (8 genes).
(B) KEGG hsa04024 cAMP signaling pathway annotation indicates that
TSHR serves as a canonical input ligand directly connected through
Gs/adenylyl cyclase/PKA to CREB. In contrast, IGF-1R is not annotated
in hsa04024; its connection to the cAMP module is indirect and lacks
canonical KEGG support. This directional annotation replaces
undirected PPI hub-distance metrics, which are confounded by hub-node
bias (e.g., AKT1, PIK3CA).
```

---

## FIGURE 3 — Triangulation + Transcriptomic Validation + Off-target (Composite)

**Purpose**: Single integrative figure combining triangulation summary, in-tissue validation, and off-target mechanism. This replaces the earlier Fig 4.

### Prompt for Nanobananapro / FigureLab

```
Create a multi-panel scientific figure with 3 panels in a portrait
layout approximately 8 × 11 inches. Publication-quality style similar
to Thyroid journal or JCI Insight.

=== PANEL A (TOP, full width): Triangulation Heatmap/Table ===

A color-coded table showing 11 genes (rows) × 4 evidence columns:

Columns: "MR", "Colocalization", "RNA-seq", "Convergence score"

Rows (top to bottom, ordered by score):
  TSHR     — MR cell: dark green (value "6.79e-17"), Coloc: dark green
             (value "0.985"), RNA-seq: dark green (value "6.5x UP,
             P=0.019"), Score: "4/4" (gold)
  HAS2     — MR: gray (—), Coloc: gray, RNA-seq: light green
             (value "1.6x UP, P=0.048"), Score: "2/4"
  HAS3     — similar, RNA-seq: light green (2.1x UP, P=0.028), 2/4
  INSR*    — MR: gray, Coloc: gray, RNA-seq: light green
             (1.5x UP, P=0.027), Score: 2/4 (footnote: *off-target)
  IRS2*    — RNA-seq: light green (1.6x UP, P=0.036), 2/4
  TNF      — MR: light green (3.57e-5 rep only, ‡), RNA-seq: gray,
             Score: 1/4
  IGF1R    — MR: yellow (borderline 0.041/0.053), Score: 1/4
  IL6      — MR: gray, RNA-seq: yellow (11.8x UP but P=0.129), 0/4
  IL1B     — RNA-seq: yellow, 0/4
  CXCL8    — yellow, 0/4
  Others   — all gray, 0/4

Use a traffic-light color scheme: dark green = highly significant,
light green = P<0.05, yellow = borderline/suggestive, gray = no data
or non-significant.
Highlight the TSHR row with a bold yellow background.
Add annotation at right: "TSHR achieves 4/4 convergence. No other gene
exceeds 2/4."
Footnote: "* off-target pathway (not in original 15-gene intersection);
‡ heterogeneity detected (Cochran Q P=0.028)"

=== PANEL B (MIDDLE LEFT): TSHR expression bar plot ===

Bar plot showing TSHR TPM across 5 samples.
X-axis: Sample labels — "Control (n=1)", "T1", "T2", "T3", "T4"
Y-axis: TSHR TPM (range 0 to 1.2)
Bars:
  Control — 0.10 TPM (dark gray)
  T1 — 0.72 TPM (red)
  T2 — 0.43 TPM (red)
  T3 — 0.49 TPM (red)
  T4 — 0.96 TPM (red)
Annotation above bars: "6.5× upregulated in TED (P = 0.019)"
Mean line shown for TED group.

=== PANEL C (MIDDLE RIGHT): INSR and insulin signaling genes ===

Bar plot showing mean TPM (Control vs TED) for key off-target genes.
For each gene, two side-by-side bars: Control (gray) and TED (orange).
Genes (X-axis):
  INSR   — Control 5.73, TED 8.35 (P=0.027)
  IRS2   — Control 5.25, TED 8.22 (P=0.036)
  FOXO1  — Control 5.24, TED 10.64 (P=0.019)
  PIK3R1 — Control 40.87, TED 55.49 (P=0.025)
  IGF1R  — Control 2.89, TED 5.05 (n.s.)
Y-axis: TPM
Significance asterisks (*P<0.05) above significant comparisons.
Title: "Insulin receptor co-upregulation in TED orbital fat"
Annotation: "First report of INSR upregulation in TED orbital tissue"

=== PANEL D (BOTTOM, full width): Off-target mechanism schematic ===

Horizontal diagram showing teprotumumab vs TSHR-targeting therapy
and their tissue distribution.

Left half (IGF-1R blocker / teprotumumab):
  Central node: "IGF-1R" in red circle with red X through it
  Radiating outward to 4 tissue branches:
    1. Orbital fat — "Desired effect: ↓ PI3K/AKT, ↓ fibroblast proliferation"
    2. Cochlea — "Off-target: 38 shared effectors → hearing loss (~30%)"
    3. Pancreas + orbital INSR — "Off-target: 51 shared effectors →
       hyperglycemia (~15%)"
    4. CNS — "Rare cognitive effects"

Right half (TSHR-targeting therapy):
  Central node: "TSHR" in blue circle with blue X through it
  Single branch: "Orbital fat only — no IGF-1R off-target"

Label at bottom: "IGF-1R blockade affects 4 tissues; TSHR-targeting
therapies theoretically affect orbital tissue only."

=== STYLING ===
- Flat, vector style
- Sans-serif font
- Consistent color scheme across panels
- Panel labels (A, B, C, D) in bold at top-left of each panel
- No decorative elements
```

**Output**: PNG ≥ 300 DPI.

**Legend** (~200 words):
```
Figure 3. Multi-modal triangulation of TED therapeutic targets.
(A) Convergence scoring matrix across four evidence types (MR,
colocalization, RNA-seq, directional consistency) for intersection
genes and selected off-target candidates. TSHR achieves the maximum
convergence score (4/4), supported by strong MR causality (P = 6.79
e-17), colocalization (PP.H4 = 0.985), and 6.5-fold upregulation in
TED orbital fat (P = 0.019). No other gene exceeds 2/4. (B) TSHR
expression across individual samples; all four TED samples
consistently exceed the single control value. (C) Co-upregulation of
INSR and downstream insulin signaling effectors (IRS2, FOXO1, PIK3R1)
in TED orbital fat; INSR upregulation (P = 0.027) is reported here
for the first time in TED tissue and provides a tissue-level
mechanism for teprotumumab-associated hyperglycemia. (D) Schematic
summary: IGF-1R blockade affects at least four tissues with 38
shared cochlear effectors and 51 shared insulin-signaling effectors,
whereas TSHR-directed approaches theoretically affect orbital tissue
only.
```

---

## SUPPLEMENTARY MATERIALS

### Supplementary Tables
| # | Content | Source |
|---|---------|--------|
| S1 | Full MR results (19 rows × all methods + sensitivity) | Excel MR_Results + MR_Sensitivity |
| S2 | 59 TED disease genes with categories and literature refs | TED_disease_genes CSV |
| S3 | Full IGF-1R pathway (36 genes) and TSHR pathway (33 genes) with sources | Gene_Sets |
| S4 | Full KEGG + GO enrichment (3 zones × top 10 each) | Enrichment |
| S5 | Full 15-gene + INSR/IRS2 transcriptomic data (TPM per sample) | Transcriptome |
| S6 | Full cochlear (38) and insulin (51) off-target gene lists | Offtarget |
| S7 | CMap LINCS L1000 top compounds per zone | CMap |

### Supplementary Figures
| # | Content | Tool |
|---|---------|------|
| S1 | Sample-level consistency heatmap (15 genes × 5 samples, z-scored) | Python matplotlib |
| S2 | Coloc regional plot for TSHR locus (-log10 P across chr14:80.5–81.7 Mb) | R/coloc |
| S3 | MR scatter plots for each gene (instrument effect × outcome effect) | TwoSampleMR |

---

## MANUSCRIPT TEXT STRUCTURE (4,000 words target)

| Section | Words | Notes |
|---------|-------|-------|
| Title | — | Finalized after Results lock |
| Abstract | 200 | Unstructured; highlight Triangulation + INSR |
| Introduction | 800 | TED / teprotumumab / TSHR-IGF1R debate / gap |
| Methods | 1,000 | Refers to Fig 1; includes IRB# |
| Results | 1,400 | Flow: Table 1 → Table 2 → Fig 2 → Fig 3 |
| Discussion | 600 | TSHR primacy / INSR novelty / limitations |

---

## BUILD ORDER (Revised)

**Today (Day 1)**:
1. Build **Table 2** (MR + Coloc combined) — from Excel
2. Deliver **Figure 1 prompt** for Nanobananapro (done in this doc)

**Day 2**:
3. Build **Figure 3 source data** (Python matplotlib drafts of panels B and C as reference for AI tool)
4. Deliver **Figure 2 prompt** for Nanobananapro (done in this doc)

**Day 3**:
5. Deliver **Figure 3 prompt** for Nanobananapro (done in this doc)
6. All Supplementary Tables in one docx from Excel

**Day 4**:
7. Supplementary Figures (matplotlib, handled by Claude or Gemini)

**Day 5–7**:
8. Methods writing
9. Results writing
10. Discussion writing
11. Introduction + Abstract

---

## DECISIONS LOCKED

✅ Thyroid target journal
✅ 2 Tables + 3 Figures Main structure
✅ No Cytoscape — all figures via AI tools or matplotlib
✅ RNA-seq framed as "2024 in-house cohort, IRB#[INSERT]"
✅ 59/36/33 gene set numbers
✅ Triangulation inside Figure 3 (no separate Table 1 for it)
✅ Existing Table 1 (Study Design) and Table 2 draft reused

---

## NEXT IMMEDIATE STEP

Build **Table 2 (MR + Coloc combined)** from Excel, output as docx.

Upon approval, proceed with Figure 1 prompt delivery and Figure 3 matplotlib reference panels.

Awaiting your go-ahead.
