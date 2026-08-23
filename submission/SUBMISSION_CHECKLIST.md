# Submission package — TED-TRAP Paper 1 (v5, final)

**Manuscript:** Genetic susceptibility and therapeutic target biology diverge at TSHR and
IGF1R in Graves disease and thyroid eye disease.

**Authors:** Jungyul Park¹; Min-Seon Kim²; Kyung-Hwa Shin³\*; Suk-Woo Yang¹\* (\*corresponding)

Master source of truth: `MANUSCRIPT_TED_TRAP_v5_MASTER.md`
(md5 `7f0786f3d8e6939bd39bf701fb46053b`, placeholders = 0, integrity audit ALL PASS).

**Target venue:** Endocrine Connections (1st choice) → Endocrine → JEI → Frontiers-Endo/BMC.
See `../internal/SUBMISSION_VENUE_CASCADE.md`.

## Files to upload

| Item | File | Status |
|---|---|---|
| Main manuscript (Word) | `MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx` | ✅ final — 13 tables embedded (Tables 1–3 + S1–S8), LibreOffice-validated |
| Cover letter (1st choice) | `COVER_LETTER_EndocrineConnections.docx` | ✅ **Endocrine Connections version** — sound-science framing, power/limitations, dual-use disclosure |
| Cover letter (fallback) | `COVER_LETTER_JEI.docx` | ⚠️ JEI version — **author list still says Yae-Eun Kang; fix before use** |
| Figure 1 | `figures/Figure1.png` | ✅ **claim figure** — screening funnel with reasons (A) + TSHR/IGF1R divergence (B); rebuilt 2026-08 |
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
- Abstract **241 words, single paragraph** with inline structured headings (Endocrine Connections format).
- Main text **4,886 words** (Intro 655 / Methods 1,062 / Results 1,495 / Discussion 1,674) — under the 5,000-word
  Endocrine Connections limit. Recount with `python3 ../scripts/26_wordcount_main_text.py`.

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
   principle is acknowledged as already established; the contribution is the worked
   within-disease demonstration.

## External review round 2 — applied 2026-08

| # | Item | Applied |
|---|---|---|
| 1 | Main text ≤ 5,000 words | ✅ 5,155 → **4,886** (Intro gaps paragraph, Results Methods-restatement, Discussion general-principle and TSHR paragraphs, Limitations duplication) |
| 2 | Abstract as a single paragraph | ✅ one paragraph, inline `Purpose:/Methods:/Results:/Conclusion:` headings, 241 words |
| 3 | IGF1R "nominal" precision | ✅ Abstract, Figure 1 and Figure 3 legends now say nominal in **BBJ and UKB but not FinnGen** (P = 0.18) |
| 4 | Soften the final Conclusion | ✅ Conclusion and the IGF1R Discussion paragraph now read "more compatible with … than with …", explicitly not establishing effector status nor excluding a modest genetic contribution |
| 5 | Figure 1 "identical instruments" error | ✅ **factual fix** — TSHR and IGF1R do not share instruments. Legend and the rendered Panel B subtitle now read "the same outcome hierarchy and analytic framework, with gene-specific *cis*-eQTL instruments"; Figure 1 regenerated |
| 6 | Complete author / corresponding-author info | ⛔ **not applied — requires the author's own data** (see to-do 1 above) |
| 7 | Justify "prespecified" | ✅ retained only for the outcome hierarchy (documented in `TrackA_MR/v5_upgrade/03_decision/TaskC_pre_analysis_plan_v1.md`); backbone genes now described as designated *a priori on biological and therapeutic grounds*, with a new Methods sentence stating that basis |
| 8 | Soften single-control tissue "significance" | ✅ Results and Figure 3 legend now report an *exploratory differential-expression signal* that "is not confirmatory because the analysis included only one control sample" |
| 9 | Stop calling FinnGen a "case series" | ✅ all 3 occurrences → "TED-enriched Graves ophthalmopathy outcome" |
| 10 | Update STROBE-MR Table S7 limitation row | ✅ item 18 now lists all ten limitations carried in the Discussion |

## Author to-do before clicking submit
1. **Author and corresponding-author details — BLOCKING, needs the author's own information.**
   The title page is still incomplete:
   - Kyung-Hwa Shin (co-corresponding): **no e-mail, no telephone**, and the degree is written as
     "MD" — confirm whether it should be "MD, PhD". Target format, matching the first
     corresponding author: `Kyung-Hwa Shin, MD[, PhD], Department of Laboratory Medicine, Pusan
     National University Hospital, Busan, Republic of Korea. E-mail: …; Tel: …`
   - Min-Seon Kim's affiliation ² reads "Department of Ophthalmology, College of Medicine, The
     Catholic University of Korea" with **no hospital named** — confirm the actual hospital.
   - ORCIDs for all four authors; confirm the portal accepts two corresponding authors; confirm all
     co-authors approve the author order.
   - Confirm IRB 2104-018-102 covers the orbital RNA-seq component and this investigator team.
2. **Cover letter** — use `COVER_LETTER_EndocrineConnections.docx` for Endocrine Connections. The JEI fallback letter still carries the superseded author name and must be corrected before reuse.
3. **Figures** — inspect each PNG at full size and confirm the journal's DPI/format requirement.
4. **Portal fields** — ORCID, author affiliations, ethics + informed-consent + data-availability
   statements (all already written in Declarations).
5. **Proof check** — after the portal builds the PDF, eyeball tables, equations, figures, and line
   numbers.
6. **GitHub security** — old commit `88eabcb` still resolvable server-side; recreate the repo or
   contact Support before the repo is shared. See `../internal/GITHUB_SECURITY_CLEANUP.md`.

## Not submitted (internal only)
`../internal/` — external-GEO decision, MVMR forensic verdict, literature landscape, venue cascade,
loop-engineering plan, security cleanup. Do not upload.
