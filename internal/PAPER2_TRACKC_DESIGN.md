# Paper 2 (Track C) — design document

_Status: ROADMAP / not started. Separate from Paper 1 (v5). Do not bolt any of
this onto v5. Same guiding principle: verify before claiming; a null is a real
possible outcome._

## Working title
Insulin-receptor and otologic off-target signatures of IGF-1R antagonism in
thyroid eye disease orbital tissue: a transcriptomic and genetic framework.

## One-line thesis
IGF-1R blockade (teprotumumab/linsitinib) engages INSR / insulin signaling and
hearing-loss–associated pathways — a transcriptomic and target-network basis for
its **metabolic (hyperglycemia) and otologic (hearing loss) adverse effects**.

## Clinical motivation (the hook)
Teprotumumab improves proptosis/diplopia in active TED but carries notable
**hyperglycemia and hearing-loss** adverse effects. IGF-1R and the insulin
receptor (INSR) share downstream signaling; IGF-1R inhibitors are known to
cross-react with INSR pharmacologically. A tissue- and genetics-level account of
this off-target liability is clinically useful and under-reported.

## Central hypothesis
The adverse-effect profile of IGF-1R antagonism is reflected in (i) orbital-tissue
co-upregulation of INSR/insulin-axis genes and (ii) the off-target gene network of
IGF-1R downstream signaling overlapping insulin and otologic pathways — distinct
from genetic *susceptibility* (which Paper 1 shows IGF1R lacks).

## Specific aims
1. **Tissue signature + external concordance.** Characterize the INSR/insulin
   cassette signature in in-house TED orbital RNA-seq and test direction
   concordance in a **public orbital-adipose dataset**. ⚠️ Per the 2026 dataset scan
   (`PAPER2_SNRNASEQ_DATASET_SCAN.md`), no open + active + orbital-fat + single-cell
   dataset exists; **GATE A is therefore GSE308553-first (open bulk orbital adipose,
   modality-matched to the in-house bulk RNA-seq)**, with controlled-access single-cell
   (Li 2022, GSA HRA000870) only if access is obtained. *Aims to remedy the n=1-control weakness.*
2. **Genetic anchoring (is INSR a signal or a pure off-target?).** Test INSR and
   insulin-axis genes with **pQTL-MR (UKB-PPP, deCODE)** and cis-eQTL MR +
   colocalization on the same GD/TED outcome hierarchy used in Paper 1. Expected/
   honest prior: INSR is an effector/off-target, not a susceptibility locus —
   confirming this strengthens the "off-target, not driver" message.
3. **Off-target liability map.** Integrate IGF-1R-downstream target sets with KEGG
   insulin pathway and hearing-loss gene sets (DisGeNET) + adverse-event pharmacology
   to build a mechanistic hypothesis for teprotumumab's metabolic/otologic AEs.

## Data sources
| Data | Role | Access |
|---|---|---|
| In-house orbital RNA-seq (4 TED + 1 ctrl) | Tissue signature (Aim 1) | IRB-restricted (on request) |
| Public TED orbital snRNA-seq (e.g. Kim 2024 JCI Insight) | External concordance (Aim 1) | Public |
| UKB-PPP / deCODE plasma pQTLs | Genetic anchoring (Aim 2) | Public |
| eQTLGen / GTEx thyroid cis-eQTL | Genetic anchoring (Aim 2) | Public |
| BBJ Graves / UKB hyperthyroid / FinnGen GO | Outcomes (Aim 2) | Public (as Paper 1) |
| KEGG hsa04910 insulin, DisGeNET hearing-loss, IGF-1R downstream | Network (Aim 3) | Public |

## Methods (per aim, reusing Paper-1 pipeline where possible)
- Aim 1: DESeq2 at biological-sample level (collapse technical replicates);
  report as exploratory; snRNA-seq pseudobulk/cell-type direction check.
- Aim 2: TwoSampleMR (pQTL & eQTL exposures) + coloc.abf; EUR-only LD reference
  (locked rule); conservative coloc filter (PP.H4 ≥ 0.8); positive/negative controls.
- Aim 3: hypergeometric enrichment with an honest universe; **report fold + q**,
  and explicitly state that IGF-1R∩INSR overlap is partly expected (shared
  downstream) — not a discovery (avoid over-reading p≈1e-54).

## Planned display items
- Fig 1: off-target concept (IGF-1R/INSR shared signaling → AE phenotypes).
- Fig 2: INSR/insulin cassette expression (in-house) + snRNA-seq concordance.
- Fig 3: pQTL/eQTL MR + coloc for INSR & insulin-axis genes across GD/TED outcomes.
- Fig 4: off-target liability network (insulin + hearing-loss overlap).
- Tables: DE summary; MR/coloc estimates; enrichment with q-values.

## Novelty positioning (honest)
- **Genuinely incremental, not paradigm:** INSR off-target → metabolic AE is known
  pharmacology; GO-level insulin response already reported (Kim 2024). What is new:
  **tissue-level transcript evidence in TED orbital fat + a genetics-anchored
  off-target framework tied to specific AEs.** Do **not** use "first report" language
  (see `13_insr_novelty_verdict.md`); verify Kim 2021 IOVS Supp Table 2 first.

## Risk register / go–no-go gates
- **Gate A (Aim 1):** if public snRNA-seq does NOT show concordant INSR/insulin
  direction → demote to "in-house exploratory only"; reconsider whether Paper 2 is viable.
- **Gate B (Aim 2):** pQTL-MR may be null for INSR. A null still supports the
  "off-target, not susceptibility" framing (consistent with Paper 1) — usable, but
  it caps the paper's novelty. Decide venue accordingly.
- **Overall:** if Aim 1 and Aim 2 are both weak, fold the off-target map into a
  short hypothesis/perspective piece rather than a full original article.

## Target venue (tentative)
Translational/ophthalmology or clinical-pharmacology leaning: e.g. *Investigative
Ophthalmology & Visual Science (IOVS)*, *Eye*, *Frontiers in Endocrinology/Pharmacology*,
or a drug-safety venue — depending on how Aims 1–2 land.

## Firewall vs Paper 1 (v5)
Paper 1 stays purely genetics (TSHR/IGF1R/CTLA4) and ships now. No INSR/insulin
content, n=1 DE, or off-target enrichment enters v5 unless a v5 reviewer explicitly
asks. Paper 2 gets its own master and integrity tracking when drafting begins.
