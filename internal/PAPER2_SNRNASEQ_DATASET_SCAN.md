# Paper 2 GATE A — public TED orbital single-cell dataset scan (2026-06)

Background-agent scan for a PUBLIC + ACTIVE + orbital-fat + single-cell TED dataset
to check INSR / insulin-cassette / IGF1R direction. **Honest bottom line: no dataset
satisfies all four criteria.** Verify CAS/activity and n from PDFs before any claim.

## Ranked candidates
| Rank | Dataset / paper | Accession | Tissue | Activity | Modality | Open? | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | Li et al. 2022, Cell Rep Med (TAO scRNA-seq) PMID 35896115 | **GSA-Human HRA000870** (not GEO) | orbital connective/adipose | unverified (CAS not confirmed) | 10x scRNA-seq (true single-cell) | **controlled-access** (data request) | best single-cell atlas; has adipogenic lipofibroblast compartment → INSR/IGF1R queryable IF access granted |
| 2 | SPARC/IGF1R orbital-fat scRNA-seq (medRxiv 2026 10.64898/2026.02.24.26346524; SfE-BES 2025) | **unverified** | orbital fat | decompression (likely stable) | 10x scRNA-seq | likely not deposited (preprint) | explicitly IGF1R-focused; contact authors for deposition |
| 3 | Kim et al. 2024 JCI Insight (orbital-fat snRNA-seq) PMID 39704170 | **none public** | orbital fat (confirmed) | **inactive, CAS<4** (disqualifier) | 10x snRNA-seq | **NO** (IRB on request) | wrong activity + not public |
| — | Sci Rep 2025/26 "orbital tissue" 10.1038/s41598-025-30716-9 | **GSE308553** (verified to this paper) | orbital adipose | some active analyses | **likely BULK** (confirm) | **GEO open** | not single-cell → no cell-type resolution, BUT see note |
| — | TED lacrimal-gland single-cell | unverified | lacrimal gland | — | — | — | **DISQUALIFIED** (lacrimal, not orbital fat) |

## Implication for GATE A (honest)
- The **ideal gate** (open + active + orbital-fat + single-cell) is **not achievable** right now.
- Cleanest tractable paths:
  1. **GSE308553 bulk orbital-adipose direction check (recommended first).** It is OPEN
     and orbital adipose, and—crucially—**bulk matches the in-house bulk RNA-seq modality**,
     making a direction comparison MORE apples-to-apples than snRNA-seq. Confirm (a) modality
     is bulk, (b) activity, (c) INSR/IGF1R present. If INSR/insulin direction concordant →
     a legitimate external direction check (not single-cell, but honest and modality-matched).
  2. **HRA000870 controlled-access request** (Li 2022) for true single-cell resolution —
     weeks of lead time + must confirm activity; pursue in parallel only if (1) is encouraging.
  3. SPARC preprint authors — contact for deposition; preprint-stage, low immediacy.
- If neither (1) nor (2) yields a concordant, defensible signal → Paper 2 stays a
  hypothesis/correspondence or is shelved (a valid non-failure outcome).

## Revised gate decision
Replace "public snRNA-seq" with **"public orbital-adipose external direction check,
GSE308553-first (bulk, modality-matched); single-cell (HRA000870) only if access
obtained."** Confirm modality/activity from the source PDFs before running.

## Sources
- Li 2022 CRM: https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(22)00235-X ; GSA: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000870
- Kim 2024 JCI Insight: https://insight.jci.org/articles/view/182352 (data on request)
- SPARC scRNA-seq: https://www.medrxiv.org/content/10.64898/2026.02.24.26346524v1.full
- Sci Rep orbital tissue / GSE308553: https://www.nature.com/articles/s41598-025-30716-9 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308553
