# Figure verification against locked master (v5)

Verified that the input data underlying each rendered figure matches the locked
values in `MANUSCRIPT_TED_TRAP_v5_MASTER.md` (Tables 2, 3, S2, S3). Verification
is at the **data level** (the values each build script plots), confirmed against
the source result tables in `TrackA_MR/v5_upgrade/`.

Canonical rendered figures (tracked, on remote):
`TrackA_MR/v5_upgrade/07_manuscript/figures/`

| Figure | Item | Source table | Master value | Verified value | OK |
|---|---|---|---|---|---|
| **Fig 2** | genes plotted (BBJ IVW/Wald) | TaskD_03d | 2,235 | 2,235 | ✅ |
| Fig 2 | Bonferroni threshold | 0.05 / 2,545 | 1.965×10⁻⁵ | 1.965×10⁻⁵ | ✅ |
| Fig 2 | Bonferroni hits | TaskD_03d | 13 | 13 (C4A, CTLA4, HLA-A, HLA-DQA2, HSD3B7, IFNGR1, MAPKAPK5, PRSS36, PSMB8, TNFSF14, TSHR, TUBB, VKORC1) | ✅ |
| **Fig 3A** | TSHR PP.H4 (BBJ/FinnGen) | TaskE_01 | 0.951 / 0.986 | 0.951 / 0.986 (top=rs179252) | ✅ |
| Fig 3A | IGF1R PP.H4 (BBJ/FinnGen) | TaskE_01 | 0.043 / 0.036 | 0.043 / 0.036 (weak) | ✅ |
| Fig 3A | CTLA4 (BBJ H3 / FinnGen H4) | TaskE_01 | 0.79 / 0.978 | 0.794 / 0.978 | ✅ |
| **Fig 3B** | TSHR β (BBJ/UKB/FinnGen) | TaskD_03d | −2.10 / −2.44 / −2.33 | −2.096 / −2.436 / −2.331 | ✅ |
| Fig 3B | IGF1R β | TaskD_03d | +0.45 / +0.30 / +0.34 | +0.446 / +0.299 / +0.342 | ✅ |
| Fig 3B | CTLA4 β | TaskD_03d | −1.74 / −1.57 / −1.77 | −1.740 / −1.569 / −1.768 | ✅ |
| **Fig 3C** | TSHR tissue (log2FC, padj) | TaskF_01 | +2.33, 0.032 | +2.328, 0.0319 | ✅ |
| Fig 3C | IGF1R tissue | TaskF_01 | +0.41, NS | +0.411, padj 0.9999 | ✅ |
| Fig 3C | CTLA4 tissue | TaskF_01 | +1.27, NS | +1.273, padj 0.815 | ✅ |
| **Fig S1** | TSHR anchor SNP | TaskB_05 | rs179252 | rs179252 (concordant) | ✅ |
| Fig S1 | anchor r² (EAS / EUR) | TaskB_05 | 0.808 / 0.965 | 0.808 / 0.965 present | ✅ |
| **Fig S2** | TNFSF14 (BBJ / FinnGen) | TaskE_02 | 0.994 / 0.019 | 0.9939 / 0.0185 | ✅ |
| Fig S2 | IFNGR1 (BBJ / FinnGen) | TaskE_02 | 0.989 / 0.019 | 0.9893 / 0.0190 | ✅ |

**Result:** all 5 figures are consistent with the final master. No rebuild required.

## Notes / caveats
- Verification is data-level. Pixel-level rendering of each PNG was not re-checked
  (cannot view images here); the canonical build scripts read the verified source
  tables, and Fig 3B axis scale was independently confirmed earlier.
- **Fig S2** is the correct *candidate-colocalization* panel
  (`FigureS2_candidate_coloc_failure.png`) — NOT the superseded insulin-cassette
  scatter that was dropped from v5.
- **Stale scripts (do not use for v5):** `make_fig3_panelA.py` and
  `make_fig3_panelC.py` at the repo root contain the insulin cassette
  (INSR/IRS2/…), which belongs to Track C and was removed from v5. The canonical
  v5 composite is built by `build_figure3_composite.R` (uses TaskE_01 backbone).

---

## 2026-08 correction: colocalization H2 label was inverted

**Defect.** `taskE_01_coloc.R` and `taskE_02_coloc_candidates.R` call `coloc.abf` with
`dataset1 = outcome GWAS` and `dataset2 = cis-eQTL`. In coloc.abf, H1 = trait-1-only and
H2 = trait-2-only, so **PP.H2 = "cis-eQTL association only"**. The manuscript text (Methods,
Results, Table 2 footnote, Table S2 footnote, Figure 3 and S2 legends) and the ggplot builders
(`build_figure3_composite.R` L37, `TrackA_MR/v5_upgrade/build_figureS2_coloc.R` L32/L40) all
described PP.H2 as "GWAS/outcome association only" — inverted. Only the Table S2 inline header
was correct, so the paper contradicted itself.

**Fix.** 10 text occurrences corrected in the master; Figure 3 and Figure S2 regenerated with the
label **"H2 (eQTL assoc. only)"** by `scripts/23_rebuild_figures_coloc_label_fix.py` (matplotlib;
R unavailable in this environment). All plotted values are read from the locked result tables —
no plotted number changed. Figure 3C now shows IGF1R `padj > 0.99` instead of `1.0e+00`, matching
the text.

**Interpretation impact (no conclusion reversed).**
- IGF1R PP.H2 = 0.79 (BBJ) / 0.68 (FinnGen) now correctly reads as *a cis-eQTL at the locus with no
  detectable disease association* — a cleaner statement of "not a susceptibility locus" than the
  previous wording, and consistent with its sub-threshold MR effect (Table S8).
- TNFSF14/IFNGR1 and the chr16p11.2 cluster "fail" in the TED outcome with PP.H2 ≈ 0.63–0.74, i.e.
  eQTL present but no TED-GWAS association — consistent with limited power in an 858-case GWAS
  rather than established distinct causal architecture.
- TSHR (PP.H4 0.951/0.986) and CTLA4 are unaffected.

**Verified values re-plotted (unchanged):** TSHR H4 0.951/0.986; IGF1R H2 0.793/0.675, H4 0.043/0.036;
CTLA4 H3 0.794 (BBJ) / H4 0.978 (FinnGen); TNFSF14 H4 0.994→0.019; IFNGR1 0.989→0.019.

---

## 2026-08: Figure 1 rebuilt as a claim figure

**Why.** The original Figure 1 was a data-sources → analysis → integration workflow diagram: it
listed what was done rather than what was found, so the paper's argument was not readable from the
figure alone. It also had a rendering defect (the FinnGen outcome text overflowed its box and was
cut by the border) and carried the superseded "TED-specific" label.

**Replacement** (`scripts/25_rebuild_figure1_claim.py`, matplotlib):
- **Panel A — screening funnel with reasons.** 4,462 druggable → 2,545 MR-testable (6,136 IVs) →
  2,235 estimable in BBJ → 13 Bonferroni hits → **0 robust novel targets**, with the 13 hits resolved
  by *why* they were removed: 5 MHC, 3 chr16p11.2 LD block, 2 discovery-only colocalization
  (TNFSF14/IFNGR1 0.99 → 0.02), 1 distinct causal variant (MAPKAPK5), 2 established loci retained.
  The verdict strip states the null together with its detectable range (median OR 2.55, Table S8).
- **Panel B — the divergence.** The same three evidence layers (MR / colocalization / orbital tissue)
  scored side by side for TSHR and IGF1R, so the title's claim is visible at a glance.

Counts are parsed from Table 3 of the locked master and asserted to sum to 13; backbone values are
read from TaskD_03d, TaskE_01 and TaskF_01. Nothing is hard-coded. The legend was rewritten to match
and every number in it was checked against the source tables.
