# TED-TRAP 논문 작성 PRD & Task Checklist

**Project**: TED-TRAP (Thyroid Eye Disease — TSHR vs IGF-1R Pathway Analysis)
**Target Journal**: Thyroid (IF 5.4), backup: JCEM (IF 5.8) → EJE (IF 5.5) → Front Endocrinol (IF 3.9)
**Current Phase**: Phase A (Number Freezing & Validation)
**Last Updated**: 2026-04-15

---

## 📋 Project Structure

```
ProjectTEDGWAS/
├── TrackA_MR/          # Mendelian Randomization
├── TrackB_Network/     # Differential Pathway Network
├── TrackC_Offtarget/   # Off-target Pathway
├── finn/               # FinnGen downloaded data
├── Literature/         # References
├── Manuscript/         # Final drafts
├── data.txt            # RNA-seq (46K genes × 5 samples)
├── RNAseqManuscript.docx   # Previously submitted RNA-seq paper
├── Validation_Full.txt     # Previous validation report
├── hsa04024_*.png/xml      # KEGG cAMP pathway reference
├── Figures.pdf         # Previous RNA-seq figures
└── Final_Numbers_Frozen.xlsx   # TO BE CREATED
```

---

## 🎯 Top-Level Goals

- [ ] Freeze all numbers before writing any text
- [ ] Create all Tables (6) and Figures (5) before body text
- [ ] Write Methods → Results → Discussion → Abstract → Introduction (in this order)
- [ ] Submit to Thyroid within 6-8 weeks

---

## PHASE A: Analysis Finalization & Number Freezing
**Duration**: 2-3 days | **Status**: IN PROGRESS

### A-0. Reproducibility Setup

- [ ] Create `/ProjectTEDGWAS/Manuscript/` subfolder structure:
    - [ ] `/drafts/`
    - [ ] `/tables/`
    - [ ] `/figures/`
    - [ ] `/supplementary/`
- [ ] Document R/Python package versions (`sessionInfo()`, `pip freeze`)
    - [ ] Save to `/ProjectTEDGWAS/Manuscript/session_info.txt`
- [ ] Set random seed = 42 for all analyses requiring randomization

**Completion Criterion**: Folder structure exists, session info saved.

---

### A-1. MR Results Reproducibility Check

**Goal**: Re-run MR to verify all previously reported numbers.

- [ ] Execute `TED_A_Step2_MR_primary.R` from scratch
    - [ ] Primary outcome: `ebi-a-GCST90018627` (Graves disease)
    - [ ] Confirm results match previous report for 8 genes
- [ ] Execute MR for Replication outcome: `ebi-a-GCST90038636` (Hyperthyroidism)
- [ ] Execute MR for Sensitivity outcome: FinnGen R12 Graves ophthalmopathy
- [ ] For each of 8 genes (TSHR, IGF1R, TNF, CTLA4, ARRB1, PPARG, IRS1, AKT1):
    - [ ] Record: beta, SE, P, nSNP for IVW
    - [ ] Record: MR-Egger beta, intercept, intercept P
    - [ ] Record: Weighted median beta, P
    - [ ] Record: Cochran Q P (heterogeneity)
    - [ ] Record: Steiger direction + P

**Critical Check**:
- [ ] TSHR IV count = 2 (confirmed in Validation_Full.txt)
- [ ] TSHR IVW P ≈ 6.79e-17 (reproducible?)
- [ ] IGF1R Primary P ≈ 0.053 (borderline, critical for interpretation)
- [ ] TNF Replication P ≈ 3.56e-05 (novel finding)

**Outputs**:
- [ ] `/TrackA_MR/results/MR_all_outcomes_final.csv` (24 rows: 8 genes × 3 outcomes)
- [ ] `/TrackA_MR/results/MR_sensitivity_final.csv`

**Completion Criterion**: All 24 MR result rows match prior reports or discrepancies documented.

---

### A-2. Colocalization Verification

- [ ] Re-run coloc for TSHR locus (chr14: 80,500,000–81,700,000)
    - [ ] eQTLGen TSHR cis-eQTL vs FinnGen R12 Graves ophthalmopathy
    - [ ] Confirm nsnps = 4,046
    - [ ] Confirm PP.H4 ≈ 0.985
- [ ] Record all 5 posterior probabilities (H0-H4)
- [ ] Verify LD r² between rs3783947 and rs179252
    - [ ] EUR 1000 Genomes: r² = 0.737 (confirm)
- [ ] **New check**: Minimum GWAS P in TSHR locus = 1.23e-07 (confirm)

**Outputs**:
- [ ] `/TrackA_MR/results/coloc_TSHR_final.rds`
- [ ] `/TrackA_MR/results/coloc_summary.csv`

**Completion Criterion**: Coloc numbers match Validation_Full.txt report exactly.

---

### A-3. Disease Gene Database Re-verification (CRITICAL)

**Rationale**: TED disease gene set of 44 genes is entirely manual curation per previous report. This is the biggest reviewer attack vector. Must re-verify database nulls.

- [ ] **CTD re-search**:
    - [ ] Query "Graves Ophthalmopathy" (D049970) → Genes tab
    - [ ] Record exact result count (expected: 0 or near-0)
    - [ ] Query "Graves Disease" (D006111) → Genes tab → Direct Evidence filter
    - [ ] Record exact result count
    - [ ] Export any hits to CSV
- [ ] **DisGeNET re-search**:
    - [ ] Query UMLS CUI C0342021 (Graves ophthalmopathy)
    - [ ] Query CUI C0018818 (Graves disease)
    - [ ] Record Score distribution
    - [ ] Export hits with Score > 0.1
- [ ] **OMIM check**:
    - [ ] Search "Graves ophthalmopathy" (MIM #300351)
    - [ ] Extract associated genes
- [ ] **Document source for EACH of the 44 manual genes**:
    - [ ] Gene | Category | Primary PMID reference
    - [ ] Save to `/TrackB_Network/data/disease_genes/TED_disease_genes_with_sources.csv`

**Completion Criterion**: Every one of the 44 genes has a PMID reference in the final CSV.

---

### A-4. Transcriptomic Validation (RNA-seq Integration)

**Input**: `/ProjectTEDGWAS/data.txt` (46,428 genes, 5 samples × 2 technical replicates)

- [ ] **Sample mapping verification**:
    - [ ] Read RNAseqManuscript.docx Methods section
    - [ ] Confirm: Sample 2 = Control (C), Samples 7/8/10/11 = TED (T1-T4)
    - [ ] If mapping differs from assumption, halt and flag
- [ ] **Normalization**:
    - [ ] Use TMM-normalized values (not raw TPM) for consistency with prior RNA-seq paper
    - [ ] If TMM column unavailable in data.txt, compute using edgeR::calcNormFactors
- [ ] **15-gene intersection analysis**:
    - [ ] For each of 15 intersection genes:
        - [ ] Gene | Zone | Control TPM/TMM | TED mean TPM/TMM | log2FC | P-value | Direction
    - [ ] Use one-sample t-test (4 TED samples vs Control value)
    - [ ] Save to `/TrackC_Offtarget/results/Transcriptome_15gene_validation.csv`
- [ ] **Extended analysis** (beyond 15 genes):
    - [ ] IGF-1R pathway (25 genes): mean log2FC, count upregulated
    - [ ] cAMP/PKA module (9 genes): mean log2FC, count upregulated
    - [ ] Pairwise comparison t-test
    - [ ] Save: `/TrackC_Offtarget/results/Pathway_level_comparison.csv`
- [ ] **INSR off-target validation**:
    - [ ] INSR, IRS1, IRS2, AKT2, PIK3CD expression
    - [ ] Document INSR upregulation (expected P ≈ 0.027)
    - [ ] Save: `/TrackC_Offtarget/results/Insulin_receptor_validation.csv`

**Completion Criterion**: TSHR upregulation (6.5×, P<0.05) reproduced; INSR upregulation confirmed; pathway-level comparison significant (P<0.05).

---

### A-5. INSR Novelty Check

- [ ] PubMed searches:
    - [ ] `"INSR" AND "thyroid eye disease"`
    - [ ] `"INSR" AND "Graves ophthalmopathy"`
    - [ ] `"insulin receptor" AND "orbital fibroblast"`
    - [ ] `"insulin receptor" AND "TED"`
    - [ ] `"IGF1R" AND "INSR" AND "orbital"`
- [ ] Document findings:
    - [ ] If prior report exists → record PMID, write "confirmation" framing
    - [ ] If no prior report → write "first report" framing
- [ ] Save: `/Literature/INSR_novelty_check.md`

**Completion Criterion**: Clear determination of whether INSR in TED orbital fat is novel finding.

---

### A-6. Sample-Level Consistency Analysis

- [ ] Create heatmap: 15 intersection genes × 5 samples (C, T1-T4)
    - [ ] Values: log2(TPM+1) z-score across samples
    - [ ] Hierarchical clustering on genes
    - [ ] Expected: TED samples cluster together, Control separates
- [ ] Individual gene plots for top findings:
    - [ ] TSHR across 5 samples (expect control=0.10, TED all >0.4)
    - [ ] INSR across 5 samples
    - [ ] HAS2, HAS3, IL6
- [ ] Save: `/TrackC_Offtarget/figures/SuppFig_Sample_Consistency.png`

**Completion Criterion**: Heatmap PNG + individual gene plots saved.

---

### A-7. Network Analysis Re-verification

- [ ] IGF-1R pathway gene set (25 genes):
    - [ ] List all genes with source (CTD/STRING/DrugBank/Manual)
- [ ] TSHR pathway gene set (24 genes):
    - [ ] List all genes with source
- [ ] Intersection results:
    - [ ] IGF1R-only: IGF1R, IGF1, IL1B, CXCL8 (confirm 4 genes)
    - [ ] TSHR-only: TSHR, ADIPOQ, FABP4 (confirm 3 genes)
    - [ ] Shared: HAS2, HAS1, TNF, ARRB1, IL6, HAS3, PPARG, CEBPA (confirm 8 genes)
- [ ] STRING PPI re-query:
    - [ ] Confidence ≥ 700
    - [ ] Record edge counts per zone
    - [ ] Record hub genes (degree ≥ threshold)

**Outputs**:
- [ ] `/TrackB_Network/data/gene_sets_final.csv` (all 3 sets with sources)
- [ ] `/TrackB_Network/results/intersection_genes_final.csv`
- [ ] `/TrackB_Network/results/ppi_statistics_final.csv`

**Completion Criterion**: All 15 intersection genes confirmed with sources documented.

---

### A-8. Enrichment Analysis (KEGG + GO)

- [ ] For each zone (IGF1R-only, TSHR-only, Shared):
    - [ ] KEGG enrichment via Enrichr → top 5 pathways with adjusted P
    - [ ] GO BP enrichment → top 5 terms
    - [ ] Note small gene set caveat (4/3/8 genes may have limited power)
- [ ] Off-target enrichment:
    - [ ] IGF1R expanded downstream (359 genes) KEGG enrichment
    - [ ] Confirm PI3K-Akt, Insulin signaling, Pathways in cancer appear
    
**Outputs**:
- [ ] `/TrackB_Network/results/enrichment_zones_final.csv`
- [ ] `/TrackC_Offtarget/results/enrichment_offtarget_final.csv`

**Completion Criterion**: Each zone has at least 3 significant (adj. P<0.05) KEGG/GO terms documented.

---

### A-9. CMap Drug Repurposing

- [ ] Submit each zone's gene list to Enrichr LINCS L1000 Chem Pert Down:
    - [ ] IGF1R-only → top 10 compounds
    - [ ] TSHR-only → top 10 compounds
    - [ ] Shared → top 10 compounds
- [ ] **M-LIGHT caveat**: Use ONLY "Enrichr LINCS L1000" (not L1000CDS2). Document clearly.
- [ ] Identify:
    - [ ] Any compounds already in TED clinical trials
    - [ ] Any known IGF-1R or TSHR targeting compounds (positive control)

**Output**: `/TrackB_Network/results/cmap_final.csv`

**Completion Criterion**: All three zones have top-10 compound lists with source = "Enrichr LINCS L1000" explicitly recorded.

---

### A-10. Off-target Intersection Verification

- [ ] IGF-1R expanded downstream: confirm 359 genes
    - [ ] Source: KEGG PI3K-Akt pathway (hsa04151) full gene list
- [ ] Hearing loss gene set: confirm 855 genes
    - [ ] Source: DisGeNET "sensorineural hearing loss" + OMIM "deafness"
- [ ] Insulin signaling gene set: confirm 137 genes
    - [ ] Source: KEGG Insulin signaling (hsa04910)
- [ ] Intersections:
    - [ ] Cochlear: 38 genes — list all
    - [ ] Insulin: 51 genes — list all
    - [ ] Sample members: COL1A2, FGFR1, BRAF, PIK3R1, EPHA2, NFKB1 (cochlear); IRS1, IRS2, AKT2, PIK3CD, PDPK1, INSR, GRB2 (insulin)

**Output**: `/TrackC_Offtarget/results/offtarget_intersections_final.csv`

**Completion Criterion**: All gene lists saved with source column indicating origin.

---

### A-11. Triangulation Table Construction

**This is the paper's central integrative analysis.**

- [ ] For each of 15 intersection genes + INSR, construct columns:
    - [ ] Gene | Zone | MR_Primary_P | MR_Rep_P | MR_FinnGen_P | Coloc_PP4 | RNAseq_log2FC | RNAseq_P | Direction_Consistent
- [ ] Convergence Score (out of 4):
    - [ ] +1 if any MR P < 0.05
    - [ ] +1 if Coloc PP4 > 0.80 (only TSHR qualifies)
    - [ ] +1 if RNAseq P < 0.05
    - [ ] +1 if direction consistent across all methods
- [ ] Expected distribution:
    - [ ] TSHR: 4/4 (star gene)
    - [ ] HAS2: 2/4 (Network + Transcriptome)
    - [ ] Others: 1-2/4

**Output**: `/Manuscript/tables/Table5_Triangulation.csv`

**Completion Criterion**: Triangulation table with scores for all 15+1 genes complete.

---

### A-12. Final Number Freezing (THE GATE)

**Goal**: Produce ONE master Excel file with ALL paper numbers.

- [ ] Create `/ProjectTEDGWAS/Final_Numbers_Frozen.xlsx` with sheets:
    - [ ] Sheet 1 — Study design: GWAS IDs, N, cases, populations
    - [ ] Sheet 2 — MR all results (8 genes × 3 outcomes × 4 methods)
    - [ ] Sheet 3 — MR sensitivity (Egger, Q, Steiger)
    - [ ] Sheet 4 — Coloc (H0-H4 + metadata)
    - [ ] Sheet 5 — Gene lists (IGF1R, TSHR, TED with sources)
    - [ ] Sheet 6 — Intersections (3 zones)
    - [ ] Sheet 7 — Enrichment (KEGG + GO for each zone)
    - [ ] Sheet 8 — CMap (top 10 per zone)
    - [ ] Sheet 9 — Off-target intersections (cochlear + insulin)
    - [ ] Sheet 10 — Transcriptome (15 genes + pathway-level + INSR)
    - [ ] Sheet 11 — Triangulation table
- [ ] Version control: date-stamp filename, never overwrite
- [ ] Distribute to selfm review: verify spot-check 10 random numbers

**Completion Criterion**: Final_Numbers_Frozen_20260415.xlsx created, all numbers cross-checked against source files.

**🔒 Once A-12 is complete, no more analysis. Phase B begins.**

---

## PHASE B: Table & Figure Production
**Duration**: 3-5 days | **Status**: PENDING

### B-1. Table 1 — Study Design & GWAS Summary

- [ ] Columns: Outcome | GWAS ID | N total | N cases | N controls | Population | Year | Source
- [ ] Rows:
    - [ ] Graves disease (ebi-a-GCST90018627): 175,465 / 2,809 / 172,656 / European / 2021 / GWAS Catalog
    - [ ] Hyperthyroidism (ebi-a-GCST90038636): 484,598 / 3,731 / 480,867 / European / 2021 / GWAS Catalog
    - [ ] Graves ophthalmopathy (FinnGen R12): 500,348 / 858 / 499,490 / Finnish / 2024 / FinnGen
- [ ] Additional rows: eQTLGen (31,684 European) | TED RNA-seq (4 TED, 1 control, inactive TED orbital fat)

**Format**: DOCX table, single column span.

---

### B-2. Table 2 — MR Primary Results

- [ ] Layout: 8 genes × columns (Method | Primary | Replication | FinnGen Sens)
- [ ] Each cell: beta (SE), P, nSNP
- [ ] Bold for significant (P<0.05) and bonferroni (P<6.25e-3)
- [ ] Footnote: explain different beta scales (log-OR vs linear mixed model)

---

### B-3. Table 3 — Intersection Gene Sources (Supplementary)

- [ ] Columns: Gene | Zone (IGF1R-only/TSHR-only/Shared) | IGF1R source | TSHR source | TED source | PMID
- [ ] 15 rows

---

### B-4. Table 4 — Enrichment Top Pathways

- [ ] 3 panels (one per zone)
- [ ] Each panel: KEGG top 5 + GO BP top 5
- [ ] Columns: Pathway | Adjusted P | Gene count | Overlapping genes

---

### B-5. Table 5 — Triangulation Table (CENTRAL TABLE)

- [ ] 16 rows (15 intersection genes + INSR)
- [ ] Columns: Gene | Zone | MR | Coloc | Network | RNA-seq | Score | Biological Role
- [ ] Visual: traffic-light coloring (green/yellow/red for strength)
- [ ] Sorted by Score descending (TSHR at top)

---

### B-6. Table 6 — Transcriptomic Validation

- [ ] 15 genes × (Zone, Control TPM, TED TPM, Fold change, log2FC, P-value, Direction)
- [ ] Include pathway-level summary row
- [ ] Supplementary: extended 359 IGF-1R downstream genes

---

### B-7. Figure 1 — Study Design Schematic

- [ ] Design elements:
    - [ ] Top: Input data (3 GWAS + eQTLGen + TED RNA-seq)
    - [ ] Middle: 4 analysis streams (MR → Coloc → Network → Off-target)
    - [ ] Bottom: Integration layer (Triangulation → Conclusions)
- [ ] Tool: BioRender or Illustrator, or PowerPoint (if no access)
- [ ] Size: full width, ~8" × 6"
- [ ] Color scheme: consistent with Thyroid journal aesthetics

---

### B-8. Figure 2 — MR Forest Plots

- [ ] Panel A: Primary (GCST90018627), log-OR scale, 8 genes
- [ ] Panel B: Replication (GCST90038636), linear scale, 8 genes
- [ ] Panel C: FinnGen Sensitivity, log-OR scale, 8 genes
- [ ] Each gene: point estimate ± 95% CI, P-value to right
- [ ] Color coding by module (blue=IGF-1R pathway, red=TSHR pathway, orange=inflammatory, gray=others)
- [ ] Vertical dashed line at β=0
- [ ] "Insufficient IV" marker for IL6, TGFB1, HAS2, IGF1

---

### B-9. Figure 3 — Differential Pathway Network (CENTRAL FIGURE)

- [ ] Upper panel: 3-zone Venn diagram
    - [ ] IGF-1R pathway (25) / TED disease (44) / TSHR pathway (24)
    - [ ] 3-zone overlap counts: 4 / 3 / 8
- [ ] Lower panel: PPI network with KEGG directed overlay
    - [ ] Nodes colored by zone
    - [ ] TSHR→GNAS→ADCY→cAMP→PKA→CREB as solid arrows (canonical)
    - [ ] IGF1R→IRS→PI3K→AKT as dashed (indirect to cAMP)
    - [ ] Extension Layer box for cAMP/PKA module
- [ ] Tool: Cytoscape for PPI, Illustrator for overlay

---

### B-10. Figure 4 — Off-target Pathway

- [ ] Center: IGF-1R node (large, highlighted)
- [ ] Radiating: 3 tissue destinations (Cochlea 38 genes, Pancreas 51 genes, CNS genes)
- [ ] Side: TSHR → orbital fat only (single destination)
- [ ] Key off-target genes labeled (IRS1, IRS2, INSR, AKT2, PIK3R1)
- [ ] Annotation: "IGF-1R blockade = 4-tissue impact; TSHR = orbital-specific"

---

### B-11. Figure 5 — Transcriptomic Triangulation

- [ ] Panel A: Heatmap 15 genes × 5 samples (log2 TPM z-score)
- [ ] Panel B: Bar plot of log2FC with significance markers
- [ ] Panel C: Pathway comparison (cAMP/PKA vs PI3K/AKT mean log2FC, with violin or box plot)
- [ ] Panel D: TSHR expression across 5 samples (individual points)

---

### Supplementary Figures

- [ ] S1: MR scatter plots (for each of 8 genes)
- [ ] S2: Venn diagrams (detailed)
- [ ] S3: Coloc regional plot (TSHR locus)
- [ ] S4: Sample consistency heatmap
- [ ] S5: Extended IGF-1R downstream transcriptome (all 359)
- [ ] S6: CMap top compounds per zone

---

## PHASE C: Manuscript Writing
**Duration**: 5-7 days | **Status**: PENDING

### Writing Order (Strict)

1. [ ] **Methods** (most objective, least interpretation) — 1-2 days
2. [ ] **Results** (follows Tables/Figures directly) — 1-2 days
3. [ ] **Discussion** (interpretation) — 2 days
4. [ ] **Introduction** (frames the paper based on what was found) — 1 day
5. [ ] **Abstract** (last, summarizes everything) — 0.5 day
6. [ ] **Title** (final) — 0.5 day

### C-1. Methods (~1,500 words)

- [ ] 2.1 Study design (brief, refers to Fig 1)
- [ ] 2.2 MR
    - [ ] 2.2.1 Exposure gene selection
    - [ ] 2.2.2 Outcome datasets (3 tiers)
    - [ ] 2.2.3 MR analysis
    - [ ] 2.2.4 Colocalization
- [ ] 2.3 Network analysis
    - [ ] 2.3.1 Drug pathway gene sets
    - [ ] 2.3.2 TED disease gene set (emphasize database nulls + manual curation)
    - [ ] 2.3.3 Intersection + PPI
    - [ ] 2.3.4 Enrichment
    - [ ] 2.3.5 Extension Layer (KEGG-based rationale)
    - [ ] 2.3.6 CMap
- [ ] 2.4 Off-target analysis
- [ ] 2.5 Transcriptomic validation (new section — RNA-seq integration)
- [ ] 2.6 Software & data availability
- [ ] 2.7 Ethics

**Reference baseline**: Previously drafted Methods at `/mnt/user-data/outputs/TED_TRAP_Methods_Draft.md`.

---

### C-2. Results (~2,000 words)

- [ ] 3.1 MR identifies TSHR as a causal gene for Graves disease
    - [ ] Refer Table 2, Figure 2
    - [ ] Acknowledge circularity, introduce coloc
- [ ] 3.2 Colocalization confirms shared causal variant at TSHR locus
    - [ ] PP.H4 = 0.985, LD r² = 0.737
    - [ ] Refer Supp Fig S3
- [ ] 3.3 Differential pathway network reveals distinct drug coverage
    - [ ] Refer Figure 3, Table 3
    - [ ] Emphasize 3-zone structure
- [ ] 3.4 Extension Layer: TSHR directly accesses cAMP/PKA via KEGG canonical pathway
    - [ ] Refer Figure 3 overlay
- [ ] 3.5 Off-target analysis: IGF-1R blockade impacts cochlear and insulin pathways
    - [ ] Refer Figure 4
    - [ ] "Shared effectors, not teprotumumab-specific" framing
- [ ] 3.6 Transcriptomic validation in TED orbital fat
    - [ ] Refer Figure 5, Table 5, Table 6
    - [ ] TSHR 6.5× upregulation (P<0.05) as headline
    - [ ] INSR co-upregulation as novel finding
- [ ] 3.7 Triangulation identifies TSHR as the strongest convergent target
    - [ ] Refer Table 5 (central)

---

### C-3. Discussion (~1,500 words)

- [ ] Paragraph 1: Summary of key findings
- [ ] Paragraph 2: TSHR as fundamental target (MR + Coloc + Transcriptome convergence)
- [ ] Paragraph 3: 30% relapse explanation (IGF-1R-independent TSHR→cAMP pathway)
- [ ] Paragraph 4: Off-target mechanism (IRS1/2/INSR shared signaling)
- [ ] Paragraph 5: Implications for TSHR-targeting therapies (K-1-70, cyclic peptides, decoy receptors — mention in general terms, not specifically TSHR-ATrap)
- [ ] Paragraph 6: Comparison to Kim et al. 2024 JCI Insight (constructive, not adversarial)
- [ ] Paragraph 7: Limitations
    - [ ] TED-specific GWAS absent
    - [ ] Blood eQTL (not orbital tissue)
    - [ ] Computational only
    - [ ] Small transcriptome sample size
    - [ ] TED disease gene manual curation
    - [ ] TSHR MR uses only 2 IVs
    - [ ] TNF heterogeneity (Q P=0.028)
- [ ] Paragraph 8: Conclusion

---

### C-4. Introduction (~800 words)

- [ ] Para 1: TED epidemiology, clinical burden, unmet needs
- [ ] Para 2: Teprotumumab breakthrough (2020 FDA), mechanism (IGF-1R)
- [ ] Para 3: Remaining problems (30% relapse, hearing loss ~30%, hyperglycemia ~15%)
- [ ] Para 4: TSHR vs IGF-1R debate, β-arrestin scaffold, TSHR-targeting emerging
- [ ] Para 5: Gap — no systematic computational comparison of the two targets; study aim

---

### C-5. Abstract (250 words)

- [ ] Background (3 sentences)
- [ ] Methods (4 sentences)
- [ ] Results (6 sentences, include key numbers: PP.H4=0.985, TSHR 6.5× up, 4-tissue off-target)
- [ ] Conclusion (2 sentences)

---

### C-6. Title Options

- [ ] Option A: "Mendelian randomization and differential pathway network analysis reveal TSHR as a convergent therapeutic target in thyroid eye disease"
- [ ] Option B: "Multi-modal computational target validation identifies distinct therapeutic advantages of TSHR over IGF-1R inhibition in thyroid eye disease"
- [ ] Option C: "Genetic and network-based evidence supports TSHR as a fundamental target for thyroid eye disease: implications for teprotumumab relapse and adverse effects"

Pick based on final message strength.

---

## PHASE D: Pre-submission Polish
**Duration**: 3-5 days | **Status**: PENDING

### D-1. Internal QC

- [ ] Every number in text cross-checked against Final_Numbers_Frozen.xlsx
- [ ] Every figure citation in text verified
- [ ] Every table citation verified
- [ ] References completeness (check all 15+ in-text citations)
- [ ] STROBE-MR checklist completed
- [ ] English polish / proofreading

### D-2. Reviewer-Simulation

- [ ] Read as hostile reviewer — list potential attacks
- [ ] Preempt each attack in text or limitations
- [ ] Likely attacks:
    - [ ] "TED-specific GWAS 없음" → addressed with FinnGen sensitivity
    - [ ] "TSHR IV 2개뿐" → addressed with coloc
    - [ ] "Manual curation bias" → addressed with CTD/DisGeNET null documentation
    - [ ] "No experimental validation" → framed as hypothesis-generating + RNA-seq support
    - [ ] "Why not test TSHR-ATrap here" → framed as future work

### D-3. Cover Letter

- [ ] Hook: TED unmet need + TSHR vs IGF-1R debate
- [ ] Novelty claims: first MR+Coloc for TSHR in TED, INSR off-target, triangulation framework
- [ ] Why Thyroid: alignment with journal's TED focus
- [ ] Author contributions, funding, conflict of interest

### D-4. Supplementary Files

- [ ] Supp Methods (expanded)
- [ ] All Supp Tables (S1-S5)
- [ ] All Supp Figures (S1-S6)
- [ ] Data availability statement
- [ ] Code availability (GitHub repo link if applicable)

### D-5. Submission Checklist

- [ ] Thyroid author guidelines read (word limits, figure format, ref style)
- [ ] Format manuscript per Thyroid requirements
- [ ] All co-author approvals (if any)
- [ ] ORCID IDs
- [ ] Ethics statement
- [ ] Upload and submit

---

## Risk Register

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|-----------|
| Number reproduction inconsistency | Medium | High | A-1 to A-11 enforce re-running |
| Kim et al. 2024 is reviewer | Medium | High | Constructive comparison framing |
| Manual curation cherry-picking attack | High | Medium | A-3 provides PMID for each gene |
| TSHR circularity attack despite coloc | Low | High | Defend with PP.H4=0.985 + RNA-seq support |
| INSR finding already published | Low | Low | A-5 checks; reframe as confirmation if needed |
| RNA-seq n=5 too small | High | Medium | Frame as validation cohort, emphasize direction consistency |
| Thyroid desk reject | Low-Medium | High | Cover letter preempts concerns |

---

## Progress Tracking

**Phase A Progress**: __ / 12 tasks complete
**Phase B Progress**: __ / 11 tables+figures complete
**Phase C Progress**: __ / 6 sections complete
**Phase D Progress**: __ / 5 QC tasks complete

**Overall**: __ / 34 major milestones

---

## Current Blocker / Next Action

**Current Task**: A-1 MR reproducibility check (Gemini in progress)

**Next Action**: Gemini to execute A-1 to A-12 in sequence, produce Final_Numbers_Frozen.xlsx.

**Gate to Phase B**: A-12 completion (all numbers frozen in single Excel).

---

*End of PRD*
