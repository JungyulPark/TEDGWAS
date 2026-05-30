# Supplementary — CONSOLIDATED FINAL fixes before lock (Claude-verified)

Combines friend's S3 edits + renumbering + the still-pending S4 title / S5 footnote.
After these, supplementary package is fully locked → move to Supp Figures.

---

## EDIT 1 — Table S3: column rename + footnote precision

**1a. First column rename:** "Ancestry Locus" → "LD reference panel"
(eQTLGen is a single predominantly-European blood eQTL dataset; EUR/EAS denote the
fine-mapping LD reference panel, NOT separate-ancestry eQTLs. Avoids misreading.)
Row values: "European (1000G)" / "East Asian (1000G)".

**1b. Footnote — replace entirely with (verified against S3 cross-ancestry r² values):**
"TSHR cis-eQTL fine-mapping was performed using SuSiE with 1000 Genomes Phase 3 European
and East Asian LD reference panels. Credible sets resolved to single lead SNPs with
PIP=1.0. Secondary credible-set lead variants showed substantial LD with rs179252 in the
corresponding ancestry-matched LD reference panel, with minimum ancestry-matched r²=0.808
in EAS and r²=0.965 in EUR, supporting retention of TSHR as a single-IV/locus-level anchor
rather than a multi-IV exposure."

Rationale (verified): cross-ancestry r² is low (EAS variants: r²(EUR)=0.207–0.229), so
"both ancestries high LD" was imprecise. "Ancestry-matched panel" is correct.

Rebuild TableS3, update MD5.

---

## EDIT 2 — Table S4: title "Table 4" → "Table S4" (still pending)

In build_supp_s4.R tab_header, fix title. Content unchanged (6,136 / minF=14.2 / 0 weak).
Rebuild, update MD5.

---

## EDIT 3 — Table S5: soften footnote final clause (still pending)

"...cannot be statistically distinguished using blood cis-eQTL instruments alone."
→ "...suggesting limited ability to distinguish independent gene-specific causal signals
using blood cis-eQTL instruments alone."
Rebuild, update MD5.

---

## EDIT 4 — Supplementary RENUMBERING (S6 excluded → close the gap)

External GEO (old S6) excluded. Renumber so numbering is continuous:

| Old | New | Content |
|-----|-----|---------|
| S1 | S1 | 13 hits |
| S2 | S2 | candidate coloc |
| S3 | S3 | TSHR fine-mapping |
| S4 | S4 | instruments + F-stat |
| S5 | S5 | chr16p11.2 cluster |
| S6 | — | external GEO (EXCLUDED) |
| **S7** | **S6** | MR sensitivity |
| **S8** | **S7** | STROBE-MR |

Apply renumbering to: file names, table titles inside each file, and
MANUSCRIPT_ASSETS_LOCKED.md. Manuscript text (not yet written) will reference the NEW
numbers (S6 = MR sensitivity, S7 = STROBE-MR). Do this now to avoid cross-reference errors later.

---

## EDIT 5 — S8(→S7) file cleanup

Ensure the corrected STROBE-MR (genetically prioritize / Supported / no off-target /
GitHub URL) is the single canonical file; delete/overwrite any stale copy still containing
the old overclaim phrases or local path c:\ProjectTEDGWAS.

---

## After all edits — FINAL supplementary package
- Tables: S1 (13 hits), S2 (coloc), S3 (fine-mapping), S4 (instruments),
  S5 (chr16 cluster), S6 (MR sensitivity), S7 (STROBE-MR). [7 tables]
- External GEO: excluded.
- All values Claude-verified, MD5-locked.
→ Then build Supplementary Figures: SF1 (TSHR locus), SF2 (tissue heatmap), SF3 (MR scatter).
