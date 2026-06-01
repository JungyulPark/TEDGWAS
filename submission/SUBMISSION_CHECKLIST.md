# JEI submission package — TED-TRAP v5

**Manuscript:** Druggable-gene-wide Mendelian randomization distinguishes
TSHR-anchored susceptibility from IGF1R effector biology in Graves disease and
thyroid eye disease.

Master source of truth: `MANUSCRIPT_TED_TRAP_v5_MASTER.md`
(md5 `42199cab697a926e507ec3f65b1aadda`, placeholders = 0).

## Files to upload

| Item | File | Status |
|---|---|---|
| Main manuscript (Word) | `MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx` | ✅ final (compressed Discussion, softened conclusion, no placeholders) |
| Cover letter | `COVER_LETTER_JEI.docx` | ✅ aligned to Abstract framing |
| Figure 1 | `figures/Figure1.png` | ✅ study design schematic |
| Figure 2 | `figures/Figure2.png` | ✅ BBJ MR volcano (2,235 genes, 13 hits, P=1.965×10⁻⁵) |
| Figure 3 | `figures/Figure3.png` | ✅ composite A coloc / B MR forest / C tissue |
| Figure S1 | `figures/FigureS1.png` | ✅ TSHR fine-mapping (rs179252 anchor) |
| Figure S2 | `figures/FigureS2.png` | ✅ candidate coloc (NOT the dropped insulin scatter) |

Tables (Main 1–3) and Supplementary (S1–S7) are embedded in the main `.docx`.

## Content inventory (verified)

- References: 27; in-text citations [1]→[27], ascending.
- Main tables: 3 (Table 1 data sources, Table 2 backbone evidence, Table 3 thirteen hits).
- Figure legends: 5 (Fig 1–3, S1–S2).
- Supplementary tables: S1–S7.
- Figure data verified against locked master — see `../FIGURE_VERIFICATION.md` (all 5 match).
- External-GEO decision: Option A (not included); honest "not externally replicated"
  limitation retained. Rationale archived in `../INTERNAL_external_GEO_research_note.md`
  (internal only — do NOT submit).

## Author still-to-do before clicking submit

1. **Word count.** Main text ≈ 4,637 words (JEI "preferably 4,500"). Within tolerance
   if Methods is excluded (≈3,673); trim Discussion further only if an editor requests.
2. **Abstract length.** ≈238–262 words depending on count method — confirm against the
   JEI 250-word abstract limit before upload.
3. **Figures.** Eyeball each PNG once at full size (data verified; pixel/style not auto-checked).
   Ensure ≥300 DPI export meets JEI figure spec.
4. **ORCID / author details / ethics & data-availability statements** per JEI portal fields.

## Data NOT in this repo (by design)

In-house RNA-seq raw counts (`data.txt`), FinnGen/eQTLGen licensed sumstats, and any
patent material are excluded via `.gitignore` and were purged from git history.
