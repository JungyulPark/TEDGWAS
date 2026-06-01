# TED-TRAP — Project master plan (paper portfolio)

_Last updated: 2026-06-01. Guiding principle: 정확하고 진실한게 생명 — accuracy
and truth first. One task at a time. Verify before locking._

This project yields a **portfolio of separate papers**, not one. Each paper has a
single thesis and a matched evidence type. Keeping them apart is a deliberate
decision: mixing a rigorous genetics paper with a fragile transcriptomic one
lowers the average rigor of both. A **scope firewall** (below) prevents bleed.

---

## Paper 1 — v5 (READY TO SUBMIT) ✅

**Title:** Druggable-gene-wide Mendelian randomization distinguishes
TSHR-anchored susceptibility from IGF1R effector biology in Graves disease and
thyroid eye disease.

- **Thesis:** Genetic susceptibility ≠ drug target. *TSHR* = expression-colocalized
  susceptibility anchor; *IGF1R* = pharmacologic effector axis (nominal MR, no
  coloc); *CTLA4* = positive control. Rigorous druggable-genome screen yields
  **no robust novel target** (informative, conservative result).
- **Evidence type:** human genetics — MR + colocalization + fine-mapping +
  cross-ancestry/cross-outcome, with exploratory in-house tissue support.
- **Strength profile:** high rigor, low discovery-novelty. Value = methodological
  rigor + conceptual framework + honest consolidation (see strength assessment).
- **Status:** master `MANUSCRIPT_TED_TRAP_v5_MASTER.md` (md5 `42199cab…`),
  placeholders = 0, 12 cross-document audits passed. Package in `submission/`.
- **IN scope:** TSHR / IGF1R / CTLA4 backbone; 2,545-gene MR; coloc; fine-mapping;
  13 Bonferroni hits; exploratory backbone tissue (3 genes only).
- **OUT of scope (firewalled to Paper 2):** INSR / insulin cassette; off-target
  enrichment (hearing loss, hyperglycemia); whole-transcriptome DE beyond the 3
  backbone genes; external GEO (Option A — not included).
- **Venue:** see `SUBMISSION_VENUE_CASCADE.md` (JEI → Endocrine Connections → Genes/Frontiers).
- **Author to-dos before submit:** word count vs JEI 4,500 (currently ~4,637;
  trim Discussion only if asked); abstract ≤250; eyeball figures at full size;
  portal fields (ORCID, ethics, data-availability).
- **Optional (parked, not yet applied):** one **interpretation-only** Discussion
  sentence linking IGF1R-as-effector to insulin-receptor–shared signaling
  (metabolic/otologic adverse effects), *with no new data*. Decision deferred;
  default is leave v5 untouched.

## Paper 2 — Track C: IGF-1R antagonist off-target axis (ROADMAP) 🚧

**Working title:** Insulin-receptor and otologic off-target signatures of IGF-1R
antagonism in thyroid eye disease orbital tissue.

- **Thesis:** IGF-1R blockade (teprotumumab/linsitinib) engages INSR / insulin
  signaling and hearing-loss–associated pathways, a transcriptomic basis for its
  metabolic and otologic adverse effects.
- **Current evidence (in repo):** in-house orbital RNA-seq INSR upregulation;
  insulin-pathway cassette co-upregulation; enrichment IGF-1R-downstream ∩ KEGG
  insulin (51 genes) and ∩ hearing-loss genes (38 genes). Files: `TrackC_Offtarget/`,
  `13_insr_novelty_verdict.md`, `19_hypergeometric_enrichment.csv`.
- **Honest weaknesses:** (1) n=4 TED vs n=1 control → fragile DE; (2) enrichment
  p-values partly reflect known IGF-1R/INSR downstream sharing (semi-tautological);
  (3) INSR novelty is "UPHOLD with caveat" — not "first report" (Kim 2024 reports
  GO-level insulin response; linsitinib INSR off-target is known pharmacology).
- **Roadmap to publishable (do NOT bolt onto v5):**
  1. Genetically anchor INSR via **pQTL-MR** (UKB-PPP / deCODE) + coloc.
  2. Replicate tissue direction in **public TED orbital snRNA-seq** (active disease,
     correct tissue) — also addresses v5's n=1 critique if a reviewer asks.
  3. Verify Kim 2021 IOVS Supp Table 2 for INSR before any novelty wording.
  4. Frame as a translational adverse-effect-mechanism paper, not a discovery.

## Future bet — pQTL-MR novelty study (OPTIONAL) 🔭

The only structurally sound route to a genuine **novel target** (v5 cannot
produce one by design). Repeat the v5 framework with plasma protein QTLs as
exposures. Same pipeline; public data; conservative coloc retained. A target that
is pQTL-MR-significant + colocalized + cross-outcome-consistent = real novelty —
but a null is a real possibility (run before claiming). May feed Paper 2 (INSR) or
stand alone.

---

## Scope firewall (prevents version divergence / thesis bleed)
1. v5 touches only TSHR/IGF1R/CTLA4 + the genome-wide MR screen. No INSR, no
   insulin cassette, no whole-transcriptome DE, no external GEO.
2. New evidence layers (pQTL, snRNA-seq, tissue DE) belong to Paper 2/3, never
   retrofitted into v5 unless a v5 reviewer explicitly requests it.
3. One master per paper. v5 master md5 is the integrity anchor; never edit a stale
   fork. Track C keeps its own master when drafted.
4. Figures: verify inputs against the owning paper's locked values before trusting
   a render. The insulin-scatter FigureS2 was a v5 contamination — never reuse it.

## Decision log (locked)
- **External GEO = Option A (excluded).** 3 cohorts did not reproduce TSHR tissue
  direction (lacrimal/inactive). "Not externally replicated" limitation stays.
- **Insulin/INSR cassette = Track C (Paper 2),** not v5.
- **EUR-only LD reference** for clumping/coloc (EAS mismatch → COJO artifacts).
- **No "first systematic"/"first report" claims.**
- **Tissue = exploratory** (n=1 control), never confirmatory.

## Data governance
IRB raw (`data.txt`), patent (TSHR-ATrap), and licensed sumstats (eQTLGen/FinnGen)
are never committed (`.gitignore`); `data.txt` was purged from history
(`GITHUB_SECURITY_CLEANUP.md`). Commit scripts + aggregate results only.
