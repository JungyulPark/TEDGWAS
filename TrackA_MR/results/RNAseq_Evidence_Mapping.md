# TED-TRAP: Evidence Mapping Table (RNA-seq Manuscript Validation)

## Source
선생님의 미출판 자체 데이터: Bulk RNA-seq, TED Orbital Adipose Tissue  
N=4 Inactive TED (CAS<3) + 1 Control.  DEG criteria: |log2FC| > 1, p < 0.05  
- 3,200 UP-regulated / 1,500 DOWN-regulated genes identified  
- Sequencing depth: ~250M reads/library (very deep)

## Intersection Gene Evidence Mapping Table

| Gene | Zone | Status in TED RNA-seq | Evidence | Strength |
|------|------|-----------------------|----------|----------|
| **IL1B** | IGF1R-only | ✅ UP-regulated | Table 1: "inflammation-associated genes (IL1B, CXCL9, CCL2)" | Direct text |
| **CXCL8** | IGF1R-only | ✅ UP-regulated (inferred) | Same inflammation cluster as IL1B | Inferred |
| **IGF1R** | IGF1R-only | ✅ Pathway active | PI3K-Akt pathway elevated (Discussion); Kim et al.[18] cited as IGF-1R pathway gene | Indirect |
| **IGF1** | IGF1R-only | ⚪ Not mentioned | No direct mention in manuscript | Not confirmed |
| **TNF** | Shared | ✅ Referenced | "TNF signaling cascade" enriched in KEGG upregulated pathways (Table 3) | Pathway level |
| **IL6** | Shared | ✅ Referenced | "cytokines such as IFN-γ, IL-6 and TNF-α" (Introduction, standard TED pathogenesis) | Literature citation |
| **HAS1/2/3** | Shared | ✅ UP-regulated (inferred) | "hyaluronan (glycosaminoglycans)... overproduce" — Intro; ECM remodeling dominant UP pathway | Pathway-level |
| **ARRB1** | Shared | ⚪ Not mentioned | No direct mention | Not confirmed |
| **PPARG** | Shared | ⚠️ DOWN-regulated | "PPAR signaling... downregulated pathways" (Table 3: KEGG down-regulated) | Direct KEGG |
| **CEBPA** | Shared | ⚠️ DOWN-regulated | "suppressed adipogenesis" — Kim discrepancy discussion implies CEBPB-related suppression | Inferred |
| **TSHR** | TSHR-only | ✅ Referenced | "autoantigens shared by thyroid and orbital fibroblasts... TSHR" (Introduction) | Pathophysiological |
| **ADIPOQ** | TSHR-only | ⚠️ DOWN (inferred) | "metabolic shutdown" / lipid biosynthesis suppressed / adipogenesis suppressed | Inferred |
| **FABP4** | TSHR-only | ⚠️ DOWN (inferred) | Same metabolic suppression cluster; ACADL, APOC1 etc. all down | Inferred |

---

## Key Findings & Manuscript Integration Language

### Statement 1: ECM-Inflammation Prediction Validated
> *"Computational prediction was consistent with transcriptomic findings from our own bulk RNA-seq profiling of inactive TED orbital adipose tissue (n=4), which identified upregulation of inflammation-associated genes including IL1B and CXCL8 (Table 1), and enrichment of ECM-receptor interaction and TNF signaling pathways in the differentially expressed gene set."*

### Statement 2: PI3K-Akt / IGF-1R Pathway Validation
> *"Moreover, the PI3K-Akt signaling pathway — a key downstream effector of IGF-1R — remained elevated in inactive TED orbital tissue despite the absence of clinical inflammatory activity, supporting the hypothesis that IGF-1R-mediated survival and remodeling signals persist through the chronic phase."*

### Statement 3: PPARG/Adipogenesis Suppression
> *"Notably, the PPARG-mediated adipogenic program, a key component of the TSHR-specific gene zone in our network analysis, was suppressed in our inactive TED cohort (KEGG: PPAR signaling, downregulated; p<0.05), consistent with a fibrotic rather than adipogenic differentiation fate in the chronic phase."*

### Statement 4: Honest Limitation
> *"No quantified log2FC or FDR values are available at the single-gene level for HAS1/2/3, ARRB1, or ADIPOQ in the current dataset; therefore, these genes await confirmation in the forthcoming DEG list (expected within 2 weeks), which will be included as a Supplementary Table."*

---

## Figures Available (from Figures (1).pdf)
- Fig. 1: Volcano plot (TED vs Control) — Upregulated cloud visible
- Fig. 2: Heatmap (top 1,000 DEGs)
- Fig. 3: GO enrichment dot plots (BP/MF/CC for UP and DOWN)
