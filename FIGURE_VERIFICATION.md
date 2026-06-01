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
