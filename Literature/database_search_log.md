# Database Search Log for TED Disease Gene Curation
**Date**: 2026-04-21
**Analyst**: Gemini AI + Manual verification

## 1. CTD (Comparative Toxicogenomics Database)
- **URL**: https://ctdbase.org/
- **Query**: Disease ID = MESH:D049970 ("Graves Ophthalmopathy")
- **Search type**: Gene-Disease Associations
- **Results**: 2,255 gene entries
  - **Direct Evidence (marker/mechanism)**: 3 genes
    - PTGS2 (COX-2): 7 references, score 8.48
    - SCD: 2 references, score 3.33
    - TSHR: 1 reference, no inference score
  - **Inference-based**: 2,252 genes (via Rosiglitazone, Methylprednisolone, Prednisone, Botulinum Toxin networks)
- **Decision**: Only Direct Evidence genes considered. Inference-based entries excluded per CTD best practices.
- **Raw data**: `/ProjectTEDGWAS/CTD_D049970_genes_20260421243921.csv`

## 2. OMIM (Online Mendelian Inheritance in Man)
- **URL**: https://omim.org/entry/275000
- **Query**: OMIM #275000 ("Graves Disease; GRD")
- **Gene-Phenotype Relationships**: 1 entry only
  - Location: 14q31, Phenotype: {Graves disease, susceptibility to, 1}
  - Mapping key: 2 (linked by gene linkage analysis)
- **Narrative text mentions**: TSHR, CTLA4, HLA-DRB1, VDR, ICAM1, IFIH1, GC, SLC26A4 (in literature review sections, NOT curated gene table)
- **Note**: OMIM does not have a separate entry for "Graves Ophthalmopathy" / TED

## 3. DisGeNET
- **URL**: https://www.disgenet.org/
- **Query**: CUI = C0342021 ("Graves ophthalmopathy")
- **Access**: Subscription required at time of analysis (April 2026)
- **Status**: Not completed

## Cross-reference: Our 59 curated genes vs CTD
- 36/59 genes found in CTD (mostly inferred)
- 23/59 genes NOT in CTD (manual curation unique)
- TSHR is the ONLY gene with Direct Evidence in both CTD and OMIM

## Conclusion
Structured databases provide sparse, predominantly inference-based coverage of TED-associated genes. Manual literature curation across 10 functional categories was the only viable approach to construct a comprehensive disease gene set for pathway intersection analysis.
