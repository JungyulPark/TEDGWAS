# Supplementary Final Fixes — Spec for Gemini (Claude-verified)

Friend review identified 5 issues. All verified against source data. Apply and re-lock.

---

## FIX 1 — Table S1: split ambiguous "Same Dir" into 4 explicit columns

**Problem (verified):** Current "UKB Same Dir" = (same direction AND P<0.05), but
"FinnGen Same Dir" = (pure direction only). Same column name, different definitions →
reviewer confusion (e.g., PRSS36 FinnGen P=0.226 shows "Yes"; VKORC1 UKB P=0.073 shows "No").

**Fix:** Replace 2 columns (UKB Same Dir, FinnGen Same Dir) with 4 columns, explicit defs:
- `UKB same dir` = sign(beta_UKB) == sign(beta_BBJ)
- `UKB nominal` = same dir AND P_UKB < 0.05
- `FinnGen same dir` = sign(beta_FinnGen) == sign(beta_BBJ)
- `FinnGen nominal` = same dir AND P_FinnGen < 0.05

**Verified correct values for all 13 hits (use these exactly):**

| Gene | UKB same dir | UKB nominal | FinnGen same dir | FinnGen nominal |
|------|------|------|------|------|
| HLA-A | No | No | Yes | No |
| HLA-DQA2 | No | No | Yes | Yes |
| CTLA4 | Yes | Yes | Yes | Yes |
| TSHR | Yes | Yes | Yes | Yes |
| C4A | Yes | No | Yes | Yes |
| HSD3B7 | Yes | Yes | No | No |
| TUBB | Yes | Yes | Yes | Yes |
| VKORC1 | Yes | No | NA | NA |
| TNFSF14 | Yes | Yes | No | No |
| PRSS36 | Yes | Yes | Yes | No |
| MAPKAPK5 | Yes | No | Yes | No |
| PSMB8 | No | No | No | No |
| IFNGR1 | Yes | No | Yes | No |

**Implementation:** Recompute in build_supp_tables.R from TaskD_03d primary-method betas/P
(don't hardcode — compute sign-match and P<0.05 so it's reproducible). Footnote must state
both definitions explicitly:
"Same direction: outcome effect sign matches BBJ discovery. Nominal: same direction and P<0.05."

---

## FIX 2 — Table S4: title "Table 4" → "Table S4"

Title currently reads "Table 4..."; must be "Table S4. cis-eQTL genetic instruments and
F-statistics for key candidate genes". (Keep 6,136 / minF=14.2 / 0 weak — all correct.)

---

## FIX 3 — Table S5: confirm top-eQTL-SNP position columns present

The 143-kb chr16p11.2 cluster claim must rest on **top eQTL SNP co-location**, not gene
position. S5 must include columns: Gene, Chr, Gene position (Mb), **Top eQTL SNP**,
**Top SNP position (Mb)**, BBJ β/P, UKB β/P, FinnGen β/P, Classification.
If top-SNP-position columns are missing, add them (from TaskD_02a pos_hg19 of each gene's
lead instrument). Footnote: cluster defined by top-SNP co-location within ~143 kb.

---

## FIX 4 — Table S8 STROBE-MR: remove overclaim language

- Item 2: "resolve causal therapeutic targets" →
  "genetically prioritize druggable susceptibility and effector axes"
- Item 17: "Confirmed TSHR...; verified IGF1R..." →
  "Supported TSHR as an expression-colocalized susceptibility anchor and IGF1R as a
  non-colocalized pharmacologic effector axis."
- Item 19: remove "off-target pathways" (Part C excluded from main) →
  "Discussed clinical relevance to TSHR/IGF1R biology and limitations of interpreting
  cis-eQTL MR as therapeutic perturbation."
- Item 21 (data availability): keep GitHub but phrase fully:
  "Analytical code and derived summary tables are available at
  https://github.com/JungyulPark/TEDGWAS. Public GWAS/eQTL data are available from their
  original sources; in-house orbital RNA-seq data are available from the corresponding
  author upon reasonable request, subject to institutional restrictions."

---

## FIX 5 — Table S7: APPROVED (no change)

UKB outcome labels fixed, pleiotropy spelling fixed. Lock as-is.
**Body-text caution (for manuscript, not the table):** describe IGF1R as "broadly
concordant positive-direction across IVW and weighted median, no significant heterogeneity
or directional pleiotropy, though MR-Egger estimates were imprecise." Do NOT call CTLA4
"robust" (it has UKB Q heterogeneity P=0.001).

---

## After all fixes
Rebuild S1/S4/S5 CSV+HTML+PNG, rebuild S8 markdown, recompute MD5, update
MANUSCRIPT_ASSETS_LOCKED.md. Report values to Claude for final verification before lock.
