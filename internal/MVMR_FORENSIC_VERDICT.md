# MVMR forensic audit — VERDICT: must-not-report (2026-08)

An unreported multivariable MR (MVMR) sits in the repo. It appears to show that TSHR's effect on
Graves ophthalmopathy vanishes after conditioning on Graves disease liability. **It must not be reported.**
Recorded here so it is never re-imported.

## The artifacts
- `TrackA_MR/results/04v3_mvmr_finngen_summary.csv` (from `TrackA_MR/scripts/03b_04v3_finngen_local.R`)
  - UVMR_TSHR_only: β −1.566, SE 0.305, P 2.82e-07, 1 SNP
  - MVMR_TSHR_adj: β +0.127, SE 0.277, P 0.657, 11 SNPs, "Cond_F" 25.35
  - MVMR_GD_liability: β +0.962, SE 0.0735, P 3.67e-07, "Cond_F" 99.28
- `TrackA_MR/results/04_mvmr_summary.csv` — superseded; its log shows the outcome silently fell through to
  `finn-b-E4_THYROID` (a thyroid-disorder phenotype, **not** TED), with an EAS BBJ exposure.
- `TrackA_MR/manuscript/Table2_PanelC_MVMR.csv` — v3-era draft table, superseded by the v5 master.
- `scripts/08_mvmr_ted_given_gd.R` — **never ran** (no log, no output).
- MVMR appears nowhere in `TrackA_MR/v5_upgrade/` and the word "MVMR"/"multivariable" appears **zero** times
  in the v5 master or the submission .docx. It was dropped, not adopted.

## Why it is unreportable (fatal defects)
1. **100% sample overlap with the outcome nested inside the exposure.** Exposure-2 (FinnGen GD) and the
   outcome (FinnGen GO) are the same N = 500,348 FinnGen R12 participants; GO cases (858) are drawn from
   within the GD case pool, and GD-without-GO patients sit in the outcome's control group. β_GD ≈ 0.96 with
   model R² = 0.958 — the signature of regressing a GWAS on a nested sub-GWAS of the same trait in the same
   people. Conditioning TSHR on this is conditioning on the outcome's own case definition: over-adjustment
   that drives *any* GD-susceptibility gene to null by construction.
2. **Correlated instruments.** rs179252 and rs1023586 are 26 kb apart at the same TSHR locus, EUR r² = 0.490
   (project standard is r² < 0.001). They were never jointly clumped. Weighted cor(β_TSHR, β_GD) = −0.477
   with them in, +0.026 with them out — the model's entire "conditional" structure is this one LD pair.
3. **Nine of eleven "TSHR instruments" are not cis** (chr1/2/6/6/11/15/16/20/21), violating the manuscript's
   ±1 Mb cis rule; all have |Z| ≤ 2.3 for TSHR expression.
4. **An MHC variant is the highest-leverage point** (rs1617322, chr6:32,706,116) — a region the manuscript
   excludes from causal interpretation. Dropping it flips β_TSHR from +0.127 to −0.061.
5. **Sign instability.** Leave-one-out β_TSHR: +0.127 → −0.041 / −0.061 / +0.196. (P stays non-significant
   throughout, but on a model that cannot identify the parameter.)
6. **"Cond_F" is not the Sanderson–Windmeijer statistic** (no genetic-covariance term; wrong denominator),
   and 69% of it comes from the two correlated SNPs above.
7. **Off-pipeline sourcing**: TSHR exposure from the OpenGWAS API (raw betas), GD/TED from local FinnGen —
   not the v5 harmonised, 1000G-EUR-reconstructed instrument set. The v5 rescale is SNP-specific, so a
   re-run on the v5 scale would change the coefficient non-proportionally.

Cleared concerns: the reported 04v3 run is **EUR-only** (no BBJ/EAS — that was the superseded run), and its
instrument rs179252 clears P < 5×10⁻⁸ by a wide margin (eQTLGen P = 2.86×10⁻⁴⁰). Heterogeneity was fine
(weighted RSS 12.02, 9 df, P = 0.21). None of this rescues defects 1–7.

Also note: the producing run **partially failed** (the PART 1 coloc errored and
`03b_coloc_FinnGen_GRAVES_OPHT.csv` was never written) yet the log's closing "Saved:" banner is a hardcoded
`cat`. Do not trust that log's success messages.

## What WAS adopted into the manuscript (and why)
The audit surfaced a real problem that stands **independently of the MVMR**, provable from study design alone:
FinnGen Graves ophthalmopathy ascertains cases among Graves disease patients and compares them with
population controls, so its associations substantially re-measure GD susceptibility rather than TED-specific
liability. The v5 master previously used cross-outcome consistency *including* FinnGen GO as evidence against
an outcome-specific artifact, and called TSHR an anchor "for GD and TED". That was an overclaim.

Applied to the master (framing only — **no number changed**; β = −2.33, P = 2.8e-7, PP.H4 = 0.986 all stand):
- Results/TSHR: consistency argument now rests on the independently ascertained GD and hyperthyroidism
  outcomes; the TED-specific estimate is explicitly interpreted conservatively.
- Results/TSHR: "anchor in GD/TED" → anchor for GD "whose signal is also recovered in a TED-enriched case
  series, although a TED-specific effect separable from Graves disease liability is not established".
- Conclusion: "anchor for GD and TED" → same softened form.
- Limitations: new "Sixth" item stating the ascertainment problem and that separating the two needs a
  within-Graves (orbitopathy vs no orbitopathy) contrast, unavailable at adequate power.

## If anyone ever wants a real answer to "TSHR → TED beyond GD?"
Not this design. Required: outcome = GO cases vs **GD-without-GO controls**; exposure-2 = a GD GWAS from a
**non-overlapping** sample; joint clumping at r² < 0.001 (EUR); cis-only TSHR instruments, MHC excluded;
palindromes resolved by allele frequency; v5-scale exposure betas; genuine Sanderson–Windmeijer conditional F
(`MVMR::strength_mvmr` with `gencov`) plus Q_A. With an 858-case outcome it will almost certainly be
underpowered — which is itself the honest answer.

## Not verified
FinnGen endpoint definitions could not be read from source (`.gz` are local-only/gitignored): the 100%
*sample* overlap is verified directly from the script's FG_META (identical N), while the near-total
*phenotype* nesting is a strong inference from β ≈ 1 and R² = 0.958 across independent loci. Analysis
chronology rests on in-file log timestamps, not file mtimes. The audit was reproduced in Python (no R
available) by parsing the stored RDS at byte level; the refit reproduces the published coefficients exactly.
