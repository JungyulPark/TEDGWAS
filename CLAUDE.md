> **15 September 2026 compact supplementary review:** Leave-one-out is now reported in Methods/Results/Discussion, Figure S2, Data 6 and STROBE item 13. Fifteen omissions comprise 11 IGF1R and 4 CTLA4 results; all nine full-set estimates and every omission were independently checked in base R. The European CTLA4 dependence on rs13030124 is disclosed in Table 2. No new multi-signal colocalization was undertaken: its absence and the single-causal-variant assumption remain explicit limitations. Formal directionality and participant overlap remain unverified. Read the current REVIEW_REPORT_KO.md for review status and audit counts.

> **2026-09-11 reviewed revision:** The current master and `submission/candidate_20260911_review/` supersede the historical numerical claims, fixed scientific assertions, journal status and build instructions below. See the current revision report before editing. Outdated SuSiE/tissue P-value claims must not be restored; IGF1R UKB remains split H2/H4. The requested eQTLGen frequency sensitivity is completed (Tables S2–S3; Supplementary Data 2–4). This is a final author-review candidate, with author details/declarations and live journal checks pending. The non-negotiable data-management rules below remain in force. The 2026-09-05 author decision to omit the manuscript AI-tool declaration is retained; journal disclosure requirements remain an explicit final-review item. "Prespecified" must not describe unregistered thresholds. The old root package remains archived in `archives/submission_pre20260905/`.

# CLAUDE.md — TED-TRAP project rules

This file orients any AI assistant working in this repository. Read it before
making changes. The guiding value of this project is **정확하고 진실한게 생명**
— accuracy and truth above all. No hallucination; verification first.

## What this project is
Druggable-gene-wide Mendelian randomization (MR) + colocalization +
orbital transcriptomics distinguishing **TSHR-anchored genetic susceptibility**
from **IGF1R pharmacologic effector biology** in Graves disease (GD) and
thyroid eye disease (TED). Current state: **v5 leave-one-out integrated revision (14 Sep 2026)**,
prepared for **Endocrine Connections**. The one live package is
`submission/candidate_20260911_review/`; everything else is archived.

## Non-negotiable data rules
1. **NEVER commit IRB raw data.** The in-house orbital RNA-seq (`data.txt`,
   FASTQ/BAM/counts) is IRB-protected (Pusan National University Hospital
   2104-018-102) and is "available on request" per the manuscript. It must never
   enter git. `.gitignore` enforces this; do not weaken it. (Note: `data.txt`
   was committed once early and has since been **purged from all history** — see
   `internal/GITHUB_SECURITY_CLEANUP.md`.)
2. **NEVER commit TSHR-ATrap patent material** (sequences, structures,
   `.pdb`/`.fasta`, AlphaFold/HDOCK outputs).
3. **Do not redistribute license-bound data** (eQTLGen full sumstats, FinnGen
   `.gz`, 1000G panels). Commit the *download script*, not the data.
4. Commit only **scripts** and **aggregate results** (summary tables, figure
   code, final figures).

## Locked scientific ground truths (do not silently change)
- **EUR-only LD reference is a locked rule** for clumping/coloc (eQTLGen is
  EUR-based; EAS LD mismatch produced artifactual COJO signals).
- **TSHR** — single-instrument locus (rs179252, chr14:81,435,985 hg19) under the
  European selection reference; the East Asian panel gives two. OR 0.12 (BBJ,
  *P* = 1.09×10⁻¹⁴) / 0.09 (UKB, 8.77×10⁻²⁸) / 0.10 (FinnGen, 2.82×10⁻⁷).
  PP.H4 0.951 BBJ / **0.226 UKB** / 0.986 FinnGen — the UKB non-colocalization is
  a headline finding and must never be dropped. The SuSiE fine-mapping layer was
  **withdrawn as invalid**; do not restore credible sets, PIPs or SuSiE r².
- **IGF1R** — 4 instruments (3 in FinnGen). OR 1.56 (BBJ, *P* = 0.0212) / 1.35
  (UKB, 0.0117) / 1.41 (FinnGen, **0.182, not significant**). PP.H4 0.073 / 0.404
  / 0.032 — UKB splits H2 0.400 vs H4 0.404, i.e. **unresolved**, not "eQTL only".
  Write it comparatively (weaker genetic support than *TSHR*), never as "IGF1R is
  not a susceptibility locus" and never as proof of an exclusive effector role.
- **CTLA4** — biologically selected comparator, not a validating positive control:
  PP.H4 0.201 BBJ (H3 0.799) / 0.953 UKB / 0.978 FinnGen, so it **fails** the
  combined BBJ-plus-FinnGen criterion. OR 0.18 / 0.21 / 0.17.
- **eQTLGen per-variant sample size**: max **31,684**, median **28,092**, min **232** across the 6,135
  instruments; 27.4% have n < 20,000. Methods reports the median and range, not just "up to 31,684".
  Instrument strength is unaffected (P < 5e-8 forces F >= 29.7 at any n), and every instrument behind
  the reported backbone genes has n 20,515-31,567.
- **eQTLGen allele-frequency sensitivity is DONE** (2026-09-11). All 13 discovery
  hits retained, no PP.H4 crossed 0.80 (max |Δ| 0.002142), *IGF1R* stayed nominal
  in BBJ/UKB and non-significant in FinnGen. Detection fell: BBJ OR 1.5 power
  14.6% → 12.5%. Never describe this analysis as unperformed again.
- **OPEN: the [9, 10] citation on "*TNFSF14* and *IFNGR1* remain relevant to follow-up"** (Discussion,
  filter paragraph). Ji et al. [9] does appear to support *IFNGR1* — published summaries list it among
  genes positively associated with GD risk and among putative causal genes, alongside *TSHR* — but
  lists **TNFSF13, not TNFSF14**. Li et al. [10] is a gut-microbiota → inflammatory-protein mediation
  study of GO (FinnGen 691 cases, a different release from our 858) and does not appear to report
  either gene. Checked from search summaries only; nature.com and PMC are egress-blocked here, so the
  papers were not read. **Authors should confirm against the PDFs** and consider either dropping [10]
  from that sentence or, if Ji et al. really does report *IFNGR1*, saying so explicitly — independent
  corroboration of a candidate our filter rejected would strengthen that paragraph's own argument.
- **Why the *TSHR* MR OR looks so large.** rs179252 shifts reconstructed expression by only **0.10546
  SD**, and the MR estimate is the per-allele effect divided by that shift: variant log OR / 0.10546
  reproduces the MR β exactly (−0.22103 → −2.09595 BBJ; −0.25691 → −2.43623 UKB; −0.24587 → −2.33149
  FinnGen). **Per-allele disease ORs are 0.80 / 0.77 / 0.78.** Source:
  `provenance/major_review_evidence.json` → `TSHR_scale_check`. Keep this explanation in the Discussion;
  without it OR 0.12 reads as implausible.
- **In UKB, *IGF1R*'s PP.H4 (0.404) is HIGHER than *TSHR*'s (0.226)** — both below 0.80. So "stronger
  *TSHR* evidence" holds unscoped only on the MR axis (TSHR beats IGF1R in all three outcomes:
  1.09×10⁻¹⁴ vs 0.0212; 8.77×10⁻²⁸ vs 0.0117; 2.82×10⁻⁷ vs 0.182). On colocalization it holds in BBJ
  and FinnGen only. Results, Discussion and Conclusions must all scope it, or say that neither gene
  reached the threshold in UKB.
- **TSHR PP.H4 across the three priors** (p12 = 10⁻⁵ / 5×10⁻⁶ / 10⁻⁶): BBJ **0.951 / 0.907 / 0.661**,
  FinnGen **0.986 / 0.972 / 0.875**, UKB 0.226 / 0.128 / 0.028. Only BBJ crosses 0.80; FinnGen keeps a
  margin. **Frequency stability and prior stability are different claims** — never let one stand in for
  the other.
- **Where the frequency-substitution results live**: Table S2 = MR estimates before/after, Table S3 =
  detection limits. **Colocalization under the four frequency scenarios is only in Supplementary
  Data 3** (324 rows) — no supplementary table carries it, so never cite Tables S2–S3 for the coloc
  comparison.
- **MHC composition of the 13 discovery genes**: exactly **5 inside** chr6:25-34 Mb (HLA-A, HLA-DQA2,
  C4A, TUBB, PSMB8) and **8 outside**. *IFNGR1* is on chr6 but at 137 Mb, so it is NOT an MHC gene —
  do not recount it as one. Results may therefore say the signal is not confined to the MHC region.
- **Prior sensitivity is already stated six times** (Abstract, Results/Robustness, Discussion opening,
  Discussion sensitivity paragraph, Limitations, Conclusions). Do not add a seventh — the TSHR Results
  paragraph labels its posteriors "under the primary prior" and stops there, which is specification,
  not another hedge.
- **UKB carries the strongest TSHR MR association** of the three outcomes (*P* = 8.77×10⁻²⁸ vs BBJ
  1.09×10⁻¹⁴ and FinnGen 2.82×10⁻⁷). That is precisely why its PP.H4 of 0.226 is worth reporting;
  keep the contrast explicit.
- **Results must not restate the figure legends.** The Figure 1 legend already carries "Each point
  represents one gene, not a SNP association" and the not-estimable-is-not-negative caveat; Results
  had duplicated both and no longer does.
- **The "screen recovered a known autoimmune locus" claim rests on BBJ**, where *CTLA4* has a single
  instrument at *P* = 5.45×10⁻¹⁵. The rs13030124 dependence is a UKB/FinnGen two-SNP issue and does not
  touch that claim — keep the two separate, and keep the takeaway after the caveat, not before it.
- **IGF1R leave-one-out detail**: BBJ (full-set *P* = 0.0212) loses nominal significance under
  rs2654980 (0.0751) and rs59467480 (0.0503); UKB (0.0117) under rs2654980 (0.468) and rs117212126
  (0.0529) — **two of four omissions in each**. **FinnGen never reached nominal significance**
  (full-set *P* = 0.182), so it had none to lose; never write "each outcome lost significance".
- **Leave-one-out is reported** (2026-09-14): 15 omissions = 11 IGF1R + 4 CTLA4; Figure S2, Data 6 and STROBE item 13. European CTLA4 support is concentrated in rs13030124; imprecision of the other SNP is not proof of no effect. New multi-signal coloc is outside this revision; retain model limitations. Do not restore the superseded "not reported" statement.
- **robust_novel = 0** after MHC + chr16p11.2 LD-spillover + cross-outcome coloc
  filtering. This is an *informative* result, not a negative one — do not reframe
  it as a discovery.
- **Tissue n=4 TED + 1 control**, biological-sample level (technical replicates
  collapsed). **Descriptive only** — with one control, control-side variance is not
  estimable, so no *P* value, no DESeq2 result and no significance language may
  appear anywhere. Earlier padj values (0.032, 6.6e-5, 0.006) are all withdrawn.
- **No "first systematic" claim** — druggable GD MR is already published.
- **FinnGen Graves ophthalmopathy is NOT a TED-specific contrast.** Cases are ascertained among Graves
  disease patients and compared with population controls, so its associations substantially re-measure GD
  susceptibility. TSHR is an anchor **for GD**, "whose signal is also recovered in a TED-enriched
  Graves ophthalmopathy outcome" — never claim a TED-specific effect separable from GD liability. See
  `internal/MVMR_FORENSIC_VERDICT.md`.
- **MVMR results in the repo are UNREPORTABLE** (100% sample overlap with the outcome nested in the
  exposure, correlated/trans/MHC instruments, sign instability). Never import
  `04v3_mvmr_finngen_summary.csv`, `04_mvmr_summary.csv`, or `Table2_PanelC_MVMR.csv` into any manuscript.

## Framing rules (reviewer-proof, locked)
- "Replicated" → "directionally reproduced".
- IGF1R: the genetic evidence is **weaker and unresolved**, not absent. Absence of
  colocalization is *not* evidence against its therapeutic role, and the manuscript
  must never imply that the teprotumumab evidence is challenged by these data.
- **IGF1R significance is outcome-specific**: nominal in BBJ (P=0.021) and UKB (P=0.012), **not** in
  FinnGen (P=0.182). Never write "nominal across outcomes" — say "directionally consistent, nominally
  significant in BBJ and UKB but not FinnGen".
- **Never call the FinnGen outcome a "case series"** — it is a case–control GWAS with population
  controls. Use "TED-enriched Graves ophthalmopathy outcome".
- **"Prespecified" applies ONLY to the outcome hierarchy** (documented in
  `TrackA_MR/v5_upgrade/03_decision/TaskC_pre_analysis_plan_v1.md`, which covers the BBJ/UKB/FinnGen
  hierarchy and TSHR — *not* IGF1R or CTLA4). The backbone genes were designated **a priori on
  biological and therapeutic grounds**; write it that way, never "prespecified backbone genes".
- **TSHR and IGF1R do NOT share instruments** (TSHR = 1 IV, IGF1R = 4 IVs). What is shared is the
  outcome hierarchy and the analytic framework — never write "identical/same instruments".
- **PP.H2 is NOT "no disease association".** IGF1R has a nominal MR association (BBJ P=0.0212, UKB
  P=0.0117). A high H2 says the *cis*-eQTL signal does not resolve to a variant shared with the
  outcome. Never write "no detectable disease/outcome association" for IGF1R — it contradicts the
  paper's own Table 2 row. Say "does not resolve to a variant shared with the outcome", and for UKB
  say the posterior is **split** between H2 (0.400) and H4 (0.404).
- **The null CONSTRAINS, it does not EXCLUDE.** Only 35.6% of genes were powered for OR≥2.0 and
  14.6% for OR≥1.5. Write "constrains additional large expression-mediated effects, particularly
  among well-powered genes, but does not exclude moderate effects" — never "excludes" / "evidence
  against" / "rules out".
- **Fine-mapping is withdrawn entirely.** The SuSiE run had an allele-harmonisation defect; no
  credible set, PIP, purity or SuSiE-derived r² may appear. The evidence layers are MR and
  colocalization; the orbital tissue is descriptive context only.
- **One master, one package.** `submission/` must not hold a second copy of the manuscript markdown
  (a mirror went stale once and re-introduced fixed errors), and must hold exactly **one**
  `candidate_*` directory — three accumulated once, each with a different manuscript, cover letter
  and figure set. Superseded packages go to `archives/submission_candidates/`.
  `scripts/audit_paper1_integrity.py` fails on either.
- Tissue evidence is *exploratory* (single control), never confirmatory.
- **Call the FinnGen outcome "TED-enriched", never "TED-specific"** (it is GO cases vs population
  controls). "TED-specific" survives ONLY where it denotes the *concept* of a TED-specific effect or
  the field's future TED GWAS — never as a label for our outcome.
- **IGF1R wording is comparative, never categorical**: "weaker genetic support for *IGF1R* than for
  *TSHR* in the two Graves disease outcomes; they do not show that *IGF1R* has no inherited
  contribution" — do not write that IGF1R "is not a susceptibility locus".
- External GEO is **not** included (Option A): three external cohorts did not
  reproduce the TSHR tissue direction (lacrimal-enriched / inactive TED). The
  "not externally replicated" limitation is honest and stays. Record:
  `internal/INTERNAL_external_GEO_research_note.md`.

## Working style
- **One task at a time (하나씩).** Verify before moving on.
- When in doubt about a number, check it against the source analysis file under
  `TrackA_MR/v5_upgrade/` — do not fill from memory.
- Figures are built by the R/Python scripts in the repo; **verify figure inputs
  against the locked master** before trusting a render (`FIGURE_VERIFICATION.md`).
  The current set is Figures 1-3 plus Figures S1-S2; Figure 2's OR axis MUST match Table 2. All figure text must remain black; asterisks denote nominal P<0.05, with an additional dagger only for BBJ discovery P<0.05/2,544. Figure3 borders denote PP.H4>=0.80, not P-value significance.
- **Run `python3 scripts/audit_paper1_integrity.py` after every manuscript edit.** It runs the
  framing guards (each one encodes a defect that was actually removed) and then the candidate's
  numeric audit, which compares every displayed value with its full-precision source.
- Master integrity is tracked by MD5. The master is stored with CRLF line endings
  (`.gitattributes -text`); `scripts/audit_paper1_integrity.py` prints the LF-normalised
  hash. Current master `MANUSCRIPT_TED_TRAP_v5_MASTER.md`: raw `9f6a59e51e78a4aa7fd0c687e4b162b8`,
  normalised `28c30a9c4c9d8f9ef45ed59f6f8836cd` (placeholders = 0). **33 references** since the 1000 Genomes panel got its own.
- **Gene names are italic; protein names and headings are not.** The bare `TSHR` in "TSHR and IGF-1R
  interact" and "recognition of TSHR as self" is the **protein** and stays bare, as do the `### TSHR`
  / `### IGF1R` / `### CTLA4` subsection headings. Everything else referring to the gene takes
  asterisks. Do not "fix" the protein or heading usages.
- **Em-dashes in the main text are spaced** (` — `), 15 of them; no unspaced ones. Keep it that way.
- **UKB and FinnGen participants do not overlap** (different national biobanks). The reason two
  population-control outcomes cannot settle TED specificity is that both compare affected people with
  the general population, and a GO case is by definition a GD case — never write that their
  *participants* overlap. (eQTLGen-vs-European-outcome overlap is a separate, real, unquantified issue.)
- ***P* is always italic** — in "*P* value" and "*P* <" too, not only before a number.
- **The Abstract carries the tissue direction and the numeric thresholds** (2026-09-16). It reports
  PP.H4 > 0.80, the 10⁻⁵ → 10⁻⁶ prior drop and "Mean orbital *TSHR* and *IGF1R* transcript levels
  exceeded the control" (verified: TED mean 0.651 vs 0.0998 TPM and 5.051 vs 2.894;
  `provenance/tissue_descriptive.csv`) — direction only, never a test. A parallel Abstract lives on
  `origin/codex/abstract-clinical-refinement-20260915` (3468bf3) and is **not** merged; its useful
  parts are now in the master, so do not re-merge that branch wholesale.
- **Length budget (Endocrine Connections):** main text ≤ 5,000 words — currently **3,921**
  (Introduction–Discussion; 3,953 with sub-headings, which is what Word reports). Abstract is a
  **single paragraph** with inline `Objective:/Methods:/Results:/Conclusions:` labels, **250** words
  (≤250). `scripts/26_wordcount_main_text.py` is the number of record.
  Recount with `python3 scripts/26_wordcount_main_text.py` after any edit.
- **Introduction is written for a first-time MR/colocalization reader** (medical-student/resident level,
  2026-09-15): MR and colocalization are each glossed with their core intuition on first use,
  "genetically proxied expression" and "druggable genes" are defined at first use, and one sentence
  previews that MR results are OR/95% CI/*P* while colocalization results are a posterior probability,
  not a *P* value. **Paragraph 3 stays impersonal** — it defines the two methods generically, with no
  "we used X" statement and only a generic "exposure"; every first-person statement of what this study
  did belongs in the final aims paragraph, which opens "Here, we...". That aims paragraph also
  discloses the in-house orbital RNA-seq as descriptive context only — the Abstract announces the
  dataset, so the Introduction must not stay silent about it, but it is never framed as an evidence
  layer, a validation or a result.
  Keep this register if the Introduction is edited again; do not revert to the terser
  original phrasing.

## Repo layout (actual)
```
submission/candidate_20260911_review/   # THE package: docx, cover letter,
                   #   STROBE-MR checklist, Figures 1-3 + S1-S2, Supplementary Data 1-6,
                   #   provenance/ (incl. maf/) and reproducibility/
archives/submission_candidates/         # superseded packages -- never submit from here
TrackA_MR/         # v5 core: MR, coloc, tissue (fine-mapping withdrawn)
  v5_upgrade/      #   final analysis results + 07_manuscript/figures (canonical PNGs)
TrackB_Network/    # network / pathway analysis
TrackC_Offtarget/  # off-target + insulin cassette (SEPARATE paper — not v5)
scripts/           # numbered R analysis + download scripts (no data)
finn/  Literature/ Manuscript/  archives/
internal/          # notes NOT for submission (external-GEO, venue, security)
MANUSCRIPT_TED_TRAP_v5_MASTER.md   # markdown source of truth
*.docx             # final submission deliverables (only these are tracked)
```
Local-only (gitignored) raw/reference data lives under each track's `data/`.
