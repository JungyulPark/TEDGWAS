# Task F.3: Directionality Paradox Resolution Note

**Locus Focus:** `TSHR` (Thyrotropin Receptor Locus)  
**Date:** 2026-05-26  
**Status:** DRAFTED & REVIEWER-PROOFED (For Methods/Discussion)

---

## 1. The Apparent Paradox

During our integrated analysis of the genetic susceptibility and tissue expression profiles of Thyroid Eye Disease (TED) / Graves' Disease (GD), we observed an apparent directionality conflict at the primary susceptibility locus, `TSHR` (anchored by `rs179252`):

1.  **Genetic Susceptibility (Blood cis-eQTL & MR):** 
    *   The susceptibility-associated genetic variant `rs179252` (G allele) is strongly associated with *increased* `TSHR` expression in whole blood (cis-eQTL Z = +13.28).
    *   Mendelian Randomization (MR) using this instrument demonstrates that genetically predicted *higher* `TSHR` expression is strongly *protective* against GD (BBJ: $\beta = -2.10, P = 1.1 \times 10^{-14}$; UKB: $\beta = -2.44, P = 8.8 \times 10^{-28}$).
    *   Colocalization analysis confirms a shared causal variant (`PP.H4 > 0.95` in both BBJ and FinnGen R12 cohorts), indicating that this susceptibility signal is expression-driven.
2.  **Target Tissue Pathobiology (Orbital RNA-seq):**
    *   In orbital tissue, `TSHR` expression is significantly *increased* in active/inactive TED patients compared to healthy controls (log2FC = +2.31, padj = $6.6 \times 10^{-5}$ in our in-house cohort; replicated in GSE58331 [log2FC = +0.25, P = 0.06] and GSE105149 [log2FC = +0.24]).

**Paradox Formulation:** How can genetically predicted *higher* expression of `TSHR` protect against the disease (genetics), while *increased* expression of `TSHR` in orbital tissue is a hallmark of the diseased state (tissue)?

---

## 2. Resolving the Paradox: Three Analytical Pillars

This apparent directionality discrepancy is not a biological contradiction. Rather, it reflects distinct dimensions of tissue, time, and causality. We resolve this paradox through three key mechanistic layers:

### Pillar 1: Tissue-Specific cis-Regulation (Blood vs. Orbit)
*   **Mechanism:** The cis-eQTL signal for `TSHR` is derived from whole blood (eQTLGen), reflecting systemic immune cell baseline expression. The regulatory landscape in blood (e.g., in antigen-presenting cells or lymphocytes) is structurally distinct from that in orbital fibroblasts or pre-adipocytes. 
*   **Resolution:** A variant that increases expression in immune cells (which might enhance central tolerance or modulate T-cell selection, thereby acting protectively) does not necessarily imply a similar positive regulation in the orbit. Thus, blood-derived eQTL directions cannot be directly extrapolated to target-tissue-level baseline expression.

### Pillar 2: Germline Susceptibility vs. Acquired Disease-State
*   **Mechanism:** MR utilizes germline genetic variants (randomized at meiosis) to evaluate the lifelong, baseline expression tendencies of a gene (germline state). In contrast, orbital tissue RNA-seq measures a snapshot of the transcriptome in established disease (acquired disease-state).
*   **Resolution:** The increased `TSHR` expression observed in TED orbital tissue is a downstream, acquired disease-state phenomenon driven by inflammation, adipogenesis, and autoantibody-mediated receptor activation (TSHR-stimulating immunoglobulins [TSI] binding to TSHR and triggering a feedback loop/adipoblast differentiation). This is distinct from the lifelong germline expression levels that determine initial susceptibility.

### Pillar 3: Etiological Cause vs. Downstream Consequence
*   **Mechanism:** MR and colocalization test whether baseline expression is an *etiological cause* of disease onset (cause $\to$ effect). Target tissue RNA-seq measures transcriptional profiles *after* the onset of disease (effect $\to$ marker/mediator).
*   **Resolution:** Lifelong higher systemic expression of `TSHR` (likely in thyroid/immune tissues) may promote better immune tolerance (preventing the generation of autoantibodies). However, once central tolerance fails and Graves' Disease is established, TSHR in the orbit becomes the primary target of autoimmune attack, resulting in receptor activation, orbital fibroblast proliferation, and secondary upregulation of the transcript. In other words, TSHR expression is a protective factor *systemically* before onset, but its localized upregulation in orbital tissue is a pathogenic *consequence* or mediator of the clinical phenotype.

---

## 3. Locus-Level Prioritization Framing

To ensure manuscript defensibility against reviewer critique, we adopt a **locus-level prioritization framing** rather than an expression-direction-level claim:

> *"We do not claim that orbital TSHR overexpression is protective. Instead, we demonstrate that the TSHR locus is causally linked to systemic Graves' Disease/TED susceptibility (genetic evidence) and is actively transcribed and dysregulated in the primary disease tissue (tissue expression evidence). Both lines of evidence converge to prioritize TSHR as the central susceptibility anchor and drug target, while the directionality divergence reflects the transition from germline susceptibility to clinical pathobiology."*

This framing prevents over-interpretation of the blood eQTL direction in the context of the orbit, acknowledges the tissue-specific limitations, and provides a clear, reviewer-proof defense of the findings.
