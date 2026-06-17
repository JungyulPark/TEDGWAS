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

## Round 2 (deeper sweep) — verified-OPEN leads + a key activity-matching correction

**Decisive context:** the in-house TED orbital RNA-seq is **inactive-stage** (companion
RNA-seq paper: "clinically stable, inactive TED, CAS < 3, orbital decompression").
So the correct external comparator is **inactive orbital fat**, not active — which is
exactly what is openly available. This makes GATE A materially MORE tractable.

| Dataset | Accession | Tissue / modality | Activity | Open? | Note |
|---|---|---|---|---|---|
| **Rosenberg-group orbital fat** | **GSE174139** (bulk **+ snRNA-seq**, intraconal fat/differentiating orbital fibroblasts) & **GSE158464** (bulk) | orbital fat | inactive | **GEO OPEN** ✅ | **best open lead** — modality-matched (bulk) to in-house + a snRNA layer; verify snRNA contents in GEO record |
| Cheng et al. 2025 Adv Sci | accession **unverified** (likely China GSA/OMIX) | orbital fat SVF, scRNA (58k cells; fibroblasts/preadipocytes) | inactive (fibrotic) | unverified | strong single-cell fit if matrix is public |
| iScience 2024 (CD169+ monocyte GO) | NGDC **HRA005673** | orbit + blood, 10x | — | controlled (HRA) | partial |
| medRxiv 2026 SPARC / bioRxiv 2026 sorafenib | none confirmed | orbital fat/fibroblast scRNA | mixed | preprint, no accession | low immediacy |

**Corrections:** Hua 2025 *Sci Rep* **GSE308553 is confirmed BULK** (not single-cell); the
"31,353-cell / 27 TAO+21 control pooled atlas" description is **unverified / a likely
conflation** — do not rely on it. Li 2022 (GSA HRA000870) and Kim 2024 (non-public)
unchanged.

### GATE A — revised again (final)
- In-house is **inactive** → compare against **open inactive orbital-fat GEO sets**.
- **Run order:** **GSE174139 + GSE158464 (open, orbital fat, bulk + snRNA) first**, then
  GSE308553 (open bulk) as a second open check; Cheng/HRA005673 only if matrices prove public.
- Check **INSR / IGF1R / insulin-cassette direction** (bulk, modality-matched to in-house).
  Concordant → legitimate external direction support; discordant/absent → hypothesis/shelve.
- Verify each GEO record's modality, group labels, and gene presence **before** running
  (audit-first; do not assume snRNA contents of GSE174139 without checking).

## Sources
- Li 2022 CRM: https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(22)00235-X ; GSA: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000870
- Kim 2024 JCI Insight: https://insight.jci.org/articles/view/182352 (data on request)
- SPARC scRNA-seq: https://www.medrxiv.org/content/10.64898/2026.02.24.26346524v1.full
- Sci Rep orbital tissue / GSE308553: https://www.nature.com/articles/s41598-025-30716-9 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308553
- Cheng 2025 Adv Sci orbital-fat SVF scRNA: https://doi.org/10.1002/advs.202511404
- iScience 2024 (orbit+blood, NGDC HRA005673): https://doi.org/10.1016/j.isci.2024.109213
- Rosenberg orbital-fat GEO (OPEN): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174139 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158464
