# Submission package — TED-TRAP Paper 1 (v5, final)

**Manuscript:** Genetic susceptibility and therapeutic target biology diverge at TSHR and
IGF1R in Graves disease and thyroid eye disease.

**Authors:** Jungyul Park¹; Min-Seon Kim²; Kyung-Hwa Shin³\*; Suk-Woo Yang¹\* (\*corresponding)

Master source of truth: `MANUSCRIPT_TED_TRAP_v5_MASTER.md`
(md5 `ade746389fd0c8f8bc3847db7bc4bfec`, placeholders = 0, integrity audit ALL PASS).

**Target venue:** Endocrine Connections (1st choice) → Endocrine → JEI → Frontiers-Endo/BMC.
See `../internal/SUBMISSION_VENUE_CASCADE.md`.

## Files to upload

| Item | File | Status |
|---|---|---|
| Main manuscript (Word) | `MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx` | ✅ final — 13 tables embedded (Tables 1–3 + S1–S8), LibreOffice-validated |
| Cover letter (1st choice) | `COVER_LETTER_EndocrineConnections.docx` | ✅ **Endocrine Connections version** — sound-science framing, power/limitations, dual-use disclosure |
| Cover letter (fallback) | `COVER_LETTER_JEI.docx` | ⚠️ JEI version — **author list still says Yae-Eun Kang; fix before use** |
| Figure 1 | `figures/Figure1.png` | ✅ study design schematic |
| Figure 2 | `figures/Figure2.png` | ✅ BBJ MR volcano (2,235 genes, 13 hits, P = 1.965×10⁻⁵) |
| Figure 3 | `figures/Figure3.png` | ✅ composite A/B/C — **regenerated 2026-08** with corrected coloc legend |
| Figure S1 | `figures/FigureS1.png` | ✅ TSHR fine-mapping (rs179252 anchor) |
| Figure S2 | `figures/FigureS2.png` | ✅ candidate coloc — **regenerated 2026-08** with corrected legend |
| Conference abstract (not for journal) | `KSO_oral_abstract_KR_EN.docx` | ✅ KSO oral, KR + EN |

Tables 1–3 and Supplementary S1–S8 are embedded in the main `.docx` (no separate table files).

## Content inventory (verified)
- References 27; in-text citations [1]→[27], all cited.
- Main tables 3; Supplementary tables **S1–S8**; figure legends 5.
- OR = exp(β) verified for all 9 backbone estimates.
- Abstract ≈ **226 words**. Main text ≈ **5,195 words** (Intro 704 / Methods 1,050 / Results 1,540 / Discussion 1,901).

## Substantive upgrades applied in the final round (2026-08)
1. **Table S8 — detectable effect range.** Design-based minimum detectable effect per gene/outcome
   (|β|min = (z₁₋α/₂ + z₀.₈₀) × SE) from the locked standard errors. Makes `robust_novel = 0`
   interpretable (median detectable OR 2.55 at the discovery threshold) and shows IGF1R's effect fell
   **below** its detection threshold in every outcome, so its sub-threshold status is expected rather
   than evidence against a genetic contribution. Script: `../scripts/22_power_analysis_null_screen.py`.
2. **Colocalization H1/H2 labels corrected (text + figures).** The pipeline passes dataset1 = outcome
   GWAS and dataset2 = *cis*-eQTL, so PP.H2 = "eQTL association only"; the manuscript and both figure
   builders had it inverted. 10 text occurrences fixed; Figures 3 and S2 rebuilt
   (`../scripts/23_rebuild_figures_coloc_label_fix.py`). No plotted value changed. See
   `../FIGURE_VERIFICATION.md`.
3. **FinnGen GO ascertainment overclaim removed.** GO cases are ascertained among Graves disease
   patients and compared with population controls, so that outcome re-measures GD susceptibility.
   TSHR is now framed as an anchor **for GD** whose signal is recovered in a TED-enriched case series,
   with a TED-specific effect explicitly not established (new Limitation "Sixth").
4. **New Limitations "Seventh" and "Eighth":** colocalization method limits (single-causal-variant
   assumption; one European allele-frequency vector used for both datasets including the East Asian
   discovery outcome; no colocalization power analysis) and the positive-control operating
   characteristic (CTLA4 colocalized only in the TED outcome, PP.H4 = 0.978 / rs1863800, and not in
   discovery, PP.H3 = 0.794 / rs231811 — one control does not establish filter sensitivity).
5. **Conceptual positioning made explicit:** the general "drug target ≠ susceptibility locus"
   principle is acknowledged as already established; the contribution is the worked, prespecified
   within-disease demonstration.

## Author to-do before clicking submit
1. **Cover letter** — use `COVER_LETTER_EndocrineConnections.docx` for Endocrine Connections. The JEI fallback letter still carries the superseded author name and must be corrected before reuse.
2. **Word count** — ~5,195 words of main text. Fine for Endocrine Connections; **over the JEI
   "preferably 4,500"** guidance. If submitting to JEI, trim Discussion (the newest limitation
   items are the least compressible; the "broader principle" paragraph is the most).
3. **Figures** — inspect each PNG at full size and confirm the journal's DPI/format requirement.
4. **Portal fields** — ORCID, author affiliations, ethics + informed-consent + data-availability
   statements (all already written in Declarations).
5. **GitHub security** — old commit `88eabcb` still resolvable server-side; recreate the repo or
   contact Support before the repo is shared. See `../internal/GITHUB_SECURITY_CLEANUP.md`.

## Not submitted (internal only)
`../internal/` — external-GEO decision, MVMR forensic verdict, literature landscape, venue cascade,
loop-engineering plan, security cleanup. Do not upload.
