# Supplementary FINAL fixes — last 2 edits before full lock (Claude-verified)

Friend's final review: S1/S3/S8 PASS. Only S4 title + S5 footnote remain.

---

## FIX 1 — Table S4: title "Table 4" → "Table S4"

The rendered PNG title still reads "Table 4. cis-eQTL genetic instruments..."
Must be "Table S4. cis-eQTL genetic instruments and F-statistics for key candidate genes".
In build_supp_s4.R, the gt tab_header title — change "Table 4" → "Table S4".
(Everything else correct: 6,136 instruments, minF=14.2, 0 weak, F=Z² noted.)
Rebuild S4 CSV+HTML+PNG, update MD5.

---

## FIX 2 — Table S5: soften footnote final clause

Current ending: "...indicating their causal signals cannot be statistically distinguished
using blood cis-eQTL instruments alone."

Replace with (reviewer-safer):
"...suggesting limited ability to distinguish independent gene-specific causal signals
using blood cis-eQTL instruments alone."

Keep everything else (top eQTL SNP cols, 143 kb span, 800 kb TSS note — all verified).
Rebuild S5, update MD5.

---

## VERIFIED-PASS (no edit needed)

- **S1**: 4 columns (UKB/FinnGen × same-dir/nominal), consistent definitions, footnote clear. LOCK.
- **S3**: rs179252 lead, eQTL P=2.86e-40, pos 81435985, all CS PIP=1.0, 0 independent
  secondary IVs (r²≥0.808 EAS / ≥0.965 EUR). Z-score (13.2844) and G/T allele are
  reported in S4 (instruments table) — intentional division of labor, no duplication. LOCK.
- **S8**: Latest version is clean — Item 2 "genetically prioritize", Item 17 "Supported...",
  Item 19 no "off-target", Item 21 GitHub URL. The forbidden phrases are ABSENT from the
  current file. LOCK.

---

## FILE-VERSION NOTE for Gemini

An OLD copy of TableS8_STROBE_MR.md (with the old overclaim phrases) appears to still
exist in one location and was uploaded by mistake. Please ensure the CORRECTED S8
(genetically prioritize / Supported / no off-target / GitHub) is the single canonical
file in 07_manuscript/, and overwrite/delete any stale copy so there is no version
confusion in the final package.

---

## After these 2 edits + file cleanup
All supplementary tables LOCKED: S1, S2, S3, S4, S5, S7, S8 (S6 excluded).
Then proceed to Supplementary Figures (S1 locus plot, S2 tissue heatmap, S3 MR scatter).
