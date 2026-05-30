# Task F.3 (REVISED): Directionality Paradox Resolution Note

**Locus Focus:** TSHR (anchored by rs179252)
**Date:** 2026-05-26
**Status:** REVISED by Claude — unverified mechanism claims removed; paradox-resolution logic retained
**Use:** Methods/Discussion backbone

> **Revision rationale:** The original draft embedded specific *mechanistic* claims
> (e.g., "increases central tolerance," "prevents autoantibody generation") that we did
> NOT measure. Those are hypotheses, not findings, and would invite reviewer attack
> ("where is the tolerance/T-cell evidence?"). The paradox is resolved WITHOUT them.
> This version states only what the data support.

---

## 1. The Apparent Paradox

In integrated analysis of TSHR (rs179252):

1. **Genetic susceptibility (blood cis-eQTL + MR + coloc):**
   - rs179252 G allele → increased TSHR expression in whole blood (eQTLGen Z=+13.28).
   - Genetically predicted higher TSHR expression is protective against GD
     (BBJ β=−2.10, P=1.1e-14; UKB β=−2.44, P=8.8e-28).
   - Colocalization confirms a shared causal variant (PP.H4>0.95 in BBJ and FinnGen).

2. **Target-tissue pathobiology (orbital RNA-seq):**
   - TSHR expression is increased in TED orbital tissue vs control
     (in-house log2FC=+2.31, padj=6.6e-5; directionally consistent in GSE58331
     [+0.25, P=0.06] and GSE105149 [+0.24]).

**Paradox:** genetically predicted higher TSHR expression is protective (genetics),
yet TSHR expression is elevated in diseased orbital tissue (tissue).

---

## 2. Resolution: Three Analytical Dimensions

The discrepancy is not a biological contradiction; it reflects differences in tissue,
time, and the causal question being asked. We do NOT invoke any specific unverified
mechanism (e.g., tolerance, T-cell selection) to resolve it.

### Pillar 1 — Tissue specificity (blood ≠ orbit)
The cis-eQTL is derived from whole blood (eQTLGen); the regulatory landscape of blood
differs from that of orbital fibroblasts/pre-adipocytes. A blood-derived eQTL direction
cannot be directly extrapolated to baseline expression regulation in the orbit. The
genetic instrument indexes a blood-expression-associated variant, not an orbital one.

### Pillar 2 — Germline tendency vs acquired disease-state
MR uses germline variants (randomized at meiosis) to index lifelong baseline expression
tendency. Orbital RNA-seq is a snapshot of the transcriptome in established disease. The
elevated orbital TSHR in TED is an acquired disease-state observation (in tissue already
undergoing inflammation, adipogenesis, and autoantibody-mediated receptor activation) —
a different quantity from the germline tendency that indexes susceptibility.

### Pillar 3 — Etiological cause vs downstream consequence
MR/coloc test whether germline-indexed expression is an etiological *cause* of disease
onset (cause → effect). Orbital RNA-seq measures expression *after* onset (a marker/
consequence). The two measure expression in opposite positions along the causal chain.
A variant indexing protective susceptibility at the germline/onset level is therefore
not expected to predict the direction of expression change in already-diseased tissue.

*(We deliberately do not assert WHY the germline direction is protective — that
mechanism is not measured here and is left as an open question for functional studies.)*

---

## 3. Locus-Level Prioritization Framing (LOCKED)

> "We do not claim that orbital TSHR overexpression is protective. We demonstrate that
> the TSHR locus is causally linked to GD/TED susceptibility (genetic evidence: MR +
> colocalization, shared causal variant rs179252) and is transcriptionally upregulated
> in the primary disease tissue (orbital RNA-seq). Both lines of evidence converge to
> prioritize the TSHR locus as the central susceptibility anchor and therapeutic target.
> The directionality divergence reflects the distinction between germline susceptibility
> and acquired disease-state pathobiology — measured in different tissues, at different
> points along the causal chain — and is not a contradiction."

This framing: (a) makes a locus-level claim, not an expression-direction claim;
(b) acknowledges the blood-vs-orbit tissue limitation; (c) does not over-interpret the
blood eQTL direction; (d) invokes no unmeasured mechanism. It is the reviewer-proof
position.

---

## 4. What we can and cannot say (quick reference)

| CAN say (data-supported) | CANNOT say (not measured) |
|--------------------------|---------------------------|
| TSHR locus causally linked to GD/TED (MR+coloc) | TSHR↑ enhances central tolerance |
| TSHR expression elevated in TED orbit | TSHR↑ prevents autoantibody generation |
| Blood eQTL direction ≠ orbital direction | Specific immune-cell mechanism of protection |
| Germline (cause) vs disease-state (consequence) differ | Why the germline direction is protective |
| Both converge to prioritize the TSHR locus | Causal mechanism of the protective effect |

---

*Revised by Claude 2026-05-26. Original logic retained; unverified mechanism claims
(tolerance, T-cell selection, autoantibody prevention) removed to avoid reviewer attack.
Mechanism of the protective germline effect is explicitly left as an open question.*
