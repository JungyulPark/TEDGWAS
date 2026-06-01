# Strengthening v5 & establishing novelty — honest options (internal)

Honest framing up front: the external-GEO attempt already taught us not to
promise a result before the data is in. Everything below is an **option with
real uncertainty**, not a guarantee. Ratings: Effort (low/med/high),
Feasibility (chance it produces usable evidence), ROI (impact on acceptance).

---

## Part A — Can the v5 weaknesses be patched *within this paper*?

| # | Weakness | Patch option | Effort | Feasibility | ROI | Verdict |
|---|---|---|---|---|---|---|
| 1 | Tissue **n=1 control** / not externally replicated | (a) Reanalyze a **public single-nucleus/single-cell TED orbital dataset** (active disease, correct tissue = orbital fat) for TSHR direction; (b) GTEx thyroid TSHR as supportive context | med | medium | **high** | **Best single patch.** Bulk GEO failed because it was lacrimal/inactive; snRNA-seq of *active orbital fat* is the right tissue. Could meaningfully soften n=1 if direction holds. Not guaranteed. |
| 2 | IGF1R "effector" = interpretation, not result | Already mitigated in text (absence-of-coloc defense; "We interpret…"). Optional: add one sentence explicitly stating the data show *absence of susceptibility architecture*, and effector status rests on prior pharmacology, not these data | low | high | med | Cheap, fully defensible. Do if a reviewer presses. |
| 3 | "expression-colocalized" label vs blood eQTL (TSHR low in blood) | Soften label to "blood-eQTL–colocalized" in 2–3 places; or add thyroid-eQTL coloc (GTEx) as confirmation | low–med | medium | med | Low cost; thyroid-eQTL coloc would actually *strengthen* the anchor claim if it holds. |
| 4 | TED oversold vs GD-heavy data | Re-balance title/abstract to "GD, with TED-specific sensitivity"; already partly done | low | high | med | Pure framing; safe. |
| 5 | Low novelty / no novel target / not "first" | **Cannot be patched within v5** — the design filters to robust=0 by construction. Honest "framework/cautionary" positioning is the only in-paper lever | — | — | — | Accept as identity of the paper. |

**Bottom line for v5:** weaknesses #1–#4 are *mitigable* (top pick: public
snRNA-seq orbital tissue for #1; thyroid-eQTL coloc helps #1/#3). #5 is
structural and cannot be retrofitted — submit v5 as the rigorous framework paper.

---

## Part B — Can we establish *genuine novelty*? (new analysis, next study)

v5 by design cannot yield a novel target. Genuine novelty needs a different
input layer or outcome. Ranked by promise:

| Path | Idea | Effort | Feasibility of NEW target | Why it's novel | Notes |
|---|---|---|---|---|---|
| **B1. pQTL-MR** ⭐ | Repeat the framework using **plasma protein QTLs** (UKB-PPP ~2,900 proteins, deCODE ~4,900) as exposures instead of blood eQTLs, MR + coloc on the same 3 GD/TED outcomes | med–high | **medium–high** | Protein-level instruments capture druggable biology blood-eQTL misses; pQTL-MR in GD/TED is far less worked than eQTL-MR | **Strongest path.** Public pQTL sumstats; same pipeline/tools. This is a *new paper*, not a v5 patch. Directly addresses "no novel target." |
| **B2. Disease-tissue eQTL** | Swap blood eQTLGen for **thyroid (GTEx) and, if available, orbital** cis-eQTLs as exposures | med | medium | Tissue-correct instruments can reveal effects masked in blood (TSHR itself is low in blood) | Smaller eQTL n in thyroid → power limited; still a cleaner biological story. |
| **B3. Larger / newer GWAS** | Re-run with the newest, larger Graves/TED GWAS (more cases, esp. TED-specific) as they appear | low–med | medium | More power → loci that were sub-threshold may pass + colocalize | Depends on availability of a better TED GWAS; FinnGen GO (n=858) is the current bottleneck. |
| **B4. snRNA-seq–anchored target discovery** | Use public TED orbital snRNA-seq to nominate cell-type–specific druggable genes, then MR/coloc-test them | high | medium | Cell-type resolution is genuinely new for TED target ID | Overlaps with A#1; could be its own translational paper. |
| **B5. Drug-target / pathway enrichment** | Formal enrichment of MR-prioritized genes against approved-drug-target sets, immune pathways | low | low–med (confirmatory, not new targets) | Adds interpretive novelty, not new targets | Cheap supporting analysis for any of the above. |

**Most realistic route to a defensible *novel* finding:** **B1 (pQTL-MR)**,
optionally combined with B2 (thyroid-eQTL) for cross-layer triangulation. Same
methods we already trust; public data; conservative coloc filter retained. A
target that is (i) significant in pQTL-MR, (ii) colocalized, and (iii) consistent
across outcomes would be a genuine novelty claim — but, as with the GEO episode,
**we will not know until the analysis runs**, and a null is a real possibility.

---

## Suggested sequencing
1. Submit **v5 as-is** (it is clean and self-consistent) — do not delay it for novelty work.
2. In parallel, start **B1 (pQTL-MR)** as the next study (the real novelty bet).
3. Keep **A#1 (public snRNA-seq orbital)** ready as (a) rebuttal material if a
   v5 reviewer attacks the n=1 tissue, and (b) a seed for B4.
4. Only expand v5's own scope if a reviewer explicitly requests it — otherwise
   new layers belong in the next paper, not bolted onto v5.
