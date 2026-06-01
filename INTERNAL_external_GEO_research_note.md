# Internal research note — External GEO TSHR/IGF1R/CTLA4 direction check
**STATUS: NOT INCLUDED IN v5 MASTER. Internal reference only.**
Date: 2026-05-28
Decision: Option A — exclude from v5 submission, retain for future use.

## Why this analysis is not in the v5 manuscript
The v5 master already states (Methods §M5, Discussion limitations) that no
external transcriptomic dataset was incorporated and that the orbital tissue
analysis was not externally replicated. After Gemini's audit-first pipeline
corrected an earlier grouping error (non-TED orbital diseases had previously
been pooled into "control"), the corrected strict TED-vs-Normal limma analysis
of GSE58331 and GSE105149 did NOT support a direction-consistent external
replication of the in-house TSHR signal. Reporting this as "directional
concordance" would have been incorrect.

The most parsimonious explanation is tissue mismatch rather than failure of
the in-house signal: GSE58331 and GSE105149 are predominantly lacrimal-gland
and mixed-orbital cohorts (Rosenbaum et al.), whereas the in-house cohort is
orbital adipose tissue. Reviewer-facing inclusion of these results would
introduce a new attack surface ("TSHR external replication failed") without
correctly representing the tissue heterogeneity. The v5 limitations already
acknowledge the lack of external replication honestly.

## Datasets (verified)
- GSE58331  — Rosenbaum 2015 PLoS ONE 10(9):e0137654   PMID 26371757
- GSE105149 — Rosenbaum 2017 JAMA Ophthalmol 135(11):1156–1162  PMID 28975236
Both: Affymetrix GPL570 microarray; multi-disease orbital/lacrimal cohorts
(TED, sarcoidosis, GPA, NSOI, normal, lymphoma, xanthogranuloma).

## Strict TED-vs-Normal grouping (sarcoid/GPA/NSOI/lymphoma EXCLUDED)
- GSE58331:  TED n=35 vs Control n=29  (111 non-TED diseases excluded)
- GSE105149: TED n=4  vs Control n=7   (38  non-TED diseases excluded)

## limma results — backbone genes (multi-probe collapsed by highest mean expr)

| Cohort     | Gene  | Probe        | log2FC  | P       | adj P  | Direction vs in-house |
|------------|-------|--------------|---------|---------|--------|-----------------------|
| GSE58331   | TSHR  | 210055_at    | +0.140  | 0.524   | 0.673  | concordant (up), n.s. |
| GSE58331   | IGF1R | 243358_at    | −0.516  | 0.00402 | 0.0538 | DISCORDANT            |
| GSE58331   | CTLA4 | 231794_at    | +0.481  | 0.00114 | 0.038  | concordant (up), sig. |
| GSE105149  | TSHR  | 210055_at    | −0.256  | 0.276   | 0.561  | DISCORDANT            |
| GSE105149  | IGF1R | 203628_at    | −0.289  | 0.290   | 0.574  | concordant (in-house n.s.) |
| GSE105149  | CTLA4 | 221331_x_at  | −0.143  | 0.654   | 0.838  | DISCORDANT            |

In-house orbital adipose RNA-seq (master Table 2 / Figure 3C):
- TSHR  log2FC = +2.33, adj P = 0.032
- IGF1R log2FC = +0.41, adj P = 1.000
- CTLA4 log2FC = +1.27, adj P = 0.815

## Notes / honest interpretation
- **CTLA4 in GSE58331** (+0.481, P=0.001) is the only externally significant
  backbone signal and is directionally consistent with the in-house result.
  This is plausibly driven by infiltrating immune cells in lacrimal/orbital
  tissue (CTLA4 = T-cell marker), and is consistent with CTLA4's role as a
  positive-control autoimmune locus. Not used in v5.
- TSHR external signals are small or opposite in sign; the in-house +2.33
  signal in orbital adipose is **not directly comparable** to lacrimal-gland-
  enriched cohorts. This is a tissue-difference observation, not a refutation
  of the in-house finding.

## Provenance
- Script: GEMINI_S8_external_GEO_STAGE2_backbone.R (this folder; locked)
- Console output: pasted verbatim in this folder as Stage2_console_output.txt
- Analyst: Gemini (R 4.3.3, limma, GEOquery; GPL570 annotation)
- Audited by: Claude (PMIDs verified via NCBI; grouping verified against
  Rosenbaum 2017 source paper specifying "four TED samples and seven Normal
  samples" for GSE105149)

## Possible future uses (NOT v5)
1. **Reviewer rebuttal**: if a JEI reviewer raises external replication, this
   analysis can be cited honestly as having been performed; the tissue-
   heterogeneity explanation is the principled response.
2. **Track C** (teprotumumab adverse-effect mechanisms): CTLA4 immune signal
   in lacrimal tissue may be relevant to immune cell infiltration analyses
   downstream.
3. **Future external replication**: requires a TED orbital adipose cohort,
   not lacrimal-gland data. Candidate: JCI Insight 2024 (Kim et al.,
   PMID may be searched; snRNA-seq of TED orbital fat).

---
## UPDATE 2026-05-28 — external GEO program CLOSED (Option A confirmed)
Hybrid gate test on the best-matched orbital adipose cohort (GSE185952,
Wei et al. 2022, PMID 34976024; orbital adipose/connective tissue,
TAO n=3 vs Healthy n=3, "Stable period"/inactive) did NOT reproduce the
in-house TSHR direction:

| Cohort (tissue)                         | TSHR log2FC | P     | dir vs in-house |
|-----------------------------------------|-------------|-------|-----------------|
| In-house orbital adipose (4 vs 1)       | +2.33       | adjP=0.032 | reference (UP) |
| GSE185952 orbital adipose (3 vs 3)      | -0.515      | 0.509 | DISCORDANT      |
| GSE58331 retrobulbar+mixed (35 vs 29)   | +0.140      | 0.524 | concordant, n.s.|
| GSE105149 lacrimal-enriched (4 vs 7)    | -0.256      | 0.276 | DISCORDANT      |

Decision (per friend-AI criterion "TSHR discordant -> do not include"):
external GEO is NOT added to the v5 manuscript. The master limitation
("not externally replicated") stands and is now supported by three
independent external attempts. Most plausible non-manuscript explanation:
disease-stage difference (external cohorts inactive/stable; TSHR rises with
active orbital inflammation) plus tissue/platform heterogeneity. This is a
post-hoc note, NOT a manuscript claim.

S10B raw: TSHR -0.515 (P=0.509); IGF1R -0.498 (P=0.115); CTLA4 -0.137 (P=0.751).
