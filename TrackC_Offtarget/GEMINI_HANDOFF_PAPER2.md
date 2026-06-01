# Antigravity (Gemini R) handoff — Paper 2 (Track C), batch 1

Role: **run R + report only.** Do not edit the v5 master, Paper 1, or any
manuscript text. Do not touch `submission/`. Do not invent numbers — if a file
is missing, print MISSING and report back.

## Run order
1. **`TrackC_Offtarget/scripts/P2_00_data_audit.R`** — audit FIRST.
   - Edit the `WORK` path at top if your local root differs from `c:/ProjectTEDGWAS`.
   - It prints what is FOUND/MISSING (pQTL sources, outcome GWAS, EUR LD, snRNA-seq)
     and pulls the existing INSR/IGF1R blood eQTL-MR from Paper 1.
   - **Paste the full console output back to Claude.** Do not proceed to step 2
     until Claude confirms grouping/paths.
2. **`TrackC_Offtarget/scripts/P2_01_pQTL_MR.R`** — pQTL-MR + coloc (Aim 2).
   - Run ONLY if audit shows pQTL + 3 outcomes + EUR LD all FOUND.
   - You must standardize the pQTL and outcome column names to the TwoSampleMR
     set noted in the script comments (SNP/beta/se/pval/effect_allele/other_allele/
     eaf/chr/pos) before MR — do not guess allele coding; if unclear, stop and report.
   - Returns `P2_01_pQTL_MR_results.csv` (+ coloc when region frames are supplied).

## What I need back (report verbatim)
- P2_00: full console output + `P2_00_audit_manifest.csv`.
- P2_01: full console + `P2_01_pQTL_MR_results.csv` (gene, outcome, method, nsnp,
  b, se, pval, OR) and min-F per gene.

## Data you need locally (public)
- **pQTL**: UKB-PPP (Sun 2023) and/or deCODE (Ferkingstad 2021) cis-pQTL summary
  stats for INSR, IGF1R (and INS/IGF1 if available). Put per-protein files under
  `data/pqtl/UKB-PPP/` or `data/pqtl/deCODE/`.
- **Outcomes**: same as Paper 1 — BBJ Graves GCST90018627, UKB hyperthyroid
  GCST90038636, FinnGen R12 Graves ophthalmopathy. Put under `data/outcomes/`.
- **LD**: 1000G EUR phase3 panel already at
  `TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/g1000_eur` (EUR-only locked rule).

## Honest expectation (so a null isn't mistaken for failure)
INSR was already NULL in blood eQTL-MR across all three outcomes (Paper 1:
P=0.93/0.56/0.46). pQTL-MR is **expected to also be null** for INSR — that is the
*supportive* result for "INSR is a pharmacologic off-target, not a susceptibility
locus," consistent with Paper 1. We are confirming a dissociation, not hunting a
new risk gene. If something IS positive, flag it clearly for re-check.

## Next (after batch 1 lands)
- P2_02: public TED orbital snRNA-seq direction check for INSR/insulin cassette
  (Aim 1 — the real strength). Script will follow once batch 1 paths are confirmed.
