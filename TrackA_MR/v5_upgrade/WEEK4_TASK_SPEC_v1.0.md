# TED-TRAP v5 — Week 4 Task Specification: Tissue RNA-seq Integration

**Version:** v1.0
**Lock date:** 2026-05-26
**Working directory:** `c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\`
**Executor:** Gemini (local R) | **Verifier:** Claude
**Prerequisite:** Weeks 1-3 COMPLETE & ARCHIVED (genetic evidence locked)

---

## 1. Objective & scope

Week 4 does NOT re-discover or re-analyze the orbital transcriptome from scratch.
It (a) re-confirms the already-validated v4.x tissue results for the **backbone genes only**,
and (b) integrates the blood-genetic evidence (MR + coloc) with the orbital-tissue
evidence into a single coherent layer, with the directionality paradox handled explicitly.

**Scope discipline (lesson from prior rejections):** focus on backbone (TSHR, IGF1R,
CTLA4). Do NOT foreground the n=5 RNA-seq or expand into the full candidate-gene list —
that diluted the message before. Tissue is a *supporting* layer under the genetic backbone.

**IMPORTANT — methods distinguishability:** the in-house RNA-seq raw data is shared with a
separate descriptive paper. Analytical methods here MUST be fully distinguishable
(this paper = backbone-gene confirmation in service of genetic triangulation, NOT a
transcriptome-wide discovery). State this in Methods.

---

## 2. What is REUSED vs RE-CONFIRMED

| Asset | Value (v4.x, locked) | Week 4 action |
|-------|---------------------|---------------|
| In-house orbital TSHR | log2FC=+2.33, padj=0.006 (TPM Ctrl=0.10, TED mean=0.65) | re-confirm once from data.txt |
| In-house orbital IGF1R | TPM Ctrl=2.89, TED=5.05, padj=0.710 (NS) | re-confirm |
| In-house orbital CTLA4 | (check — immune gene, may be low in adipose) | compute |
| GSE58331 (Rosenbaum 2015 PLoS One 10:e0137654, PMID 26371757; 23 TED orbital + 20 normal) | TSHR +0.249 P=0.061 | reuse (no re-analysis) |
| GSE105149 (Rosenbaum 2017 JAMA Ophthalmol, PMID 28910444; 4 lacrimal + 7 normal) | TSHR +0.24 | reuse |

In-house pipeline = DESeq2; GEO = limma. (Do NOT mix.)

---

## 3. Tasks

### F.1 — Backbone tissue re-confirmation (in-house)
- Input: `c:\ProjectTEDGWAS\data.txt` (46,427 × 20, FPKM+TPM)
- Samples: Sample#2 = Control; Samples 7/8/10/11 = TED (4 inactive TED + 1 control)
- For TSHR, IGF1R, CTLA4: extract TPM per sample, compute TED-vs-Control comparison.
- Re-confirm TSHR log2FC ≈ +2.33, padj ≈ 0.006 (validation that pipeline reproduces v4.x).
- Output per gene: TPM per sample, mean Ctrl, mean TED, log2FC, p, padj.
- If TSHR does NOT reproduce ≈+2.33/0.006 → STOP and report (pipeline/data drift).

**Note on n=1 control:** with a single control, formal DESeq2 padj for individual genes
is fragile. The v4.x padj=0.006 came from the full-model dispersion estimate. Re-confirm
using the same approach; report honestly that n=1 control limits per-gene inference and
that tissue is a *supporting/descriptive* layer, not a powered DE analysis.

### F.2 — Genetic × Tissue integration table
Build one master table, backbone genes as rows, evidence layers as columns:

| Gene | MR β (P), BBJ | MR replication (UKB) | Coloc PP.H4 (BBJ/FinnGen) | Coloc verdict | Tissue log2FC (padj) | Tissue direction | Integrated role |
|------|---------------|---------------------|---------------------------|---------------|---------------------|------------------|-----------------|
| TSHR | −2.10 (1.1e-14) | −2.44 (8.8e-28) | 0.951 / 0.986 | shared causal | +2.33 (0.006) | UP in TED | upstream anchor, expression-driven |
| IGF1R | +0.45 (0.021) | +0.30 (0.012) | 0.043 / 0.036 (H2) | no coloc | +0.41* (0.710 NS) | UP (NS) | effector, not expression-mediated |
| CTLA4 | −1.74 (5.5e-15) | −1.57 (0.002) | 0.206 / 0.978 | pos control | (compute) | (compute) | known autoimmune locus |

*IGF1R tissue log2FC = log2(5.05/2.89) ≈ +0.81; confirm exact value from F.1.

### F.3 — Directionality paradox documentation
Produce a clear written explanation (for Methods/Discussion) of the apparent paradox
and its resolution (see Section 4 below). This is the intellectual core of the integration.

---

## 4. DIRECTIONALITY PARADOX — resolution logic (LOCKED)

**Observed:**
- Blood cis-eQTL: rs179252 G allele → TSHR expression ↑ (Z=+13.28)
- MR: higher TSHR expression → lower GD risk (β=−2.10, protective)
- Coloc: blood eQTL and GD GWAS share causal variant rs179252 (H4>0.95)
- Orbital tissue: TSHR expression ↑ in TED (log2FC=+2.33)

**Apparent paradox:** blood "TSHR↑ = protective" vs orbital "TSHR↑ = disease state."

**Resolution (three layers):**
1. **Tissue specificity** — blood eQTL ≠ orbital eQTL; the same variant has different
   expression consequences across tissues.
2. **Germline vs disease-state** — blood eQTL reflects germline (lifelong) expression
   tendency; orbital tissue reflects acquired disease-state expression.
3. **Causal direction** — MR/coloc test germline expression as a *cause* of disease
   onset; tissue measures expression as a *consequence/marker* in established disease.

**Framing (locus-level prioritization):**
We do NOT claim orbital TSHR expression is protective. We claim the TSHR *locus* is
causally linked to GD/TED susceptibility (genetic) AND is transcriptionally active in
orbital tissue (tissue). The germline-vs-disease-state direction difference is not a
contradiction but reflects tissue/temporal context; both point to the TSHR locus as a
prioritized target. This is consistent with the locus-level-prioritization framing
established in Week 1.

**Reviewer-proofing:** never say "orbital TSHR expression is protective"; never present
the directions as contradictory without the resolution; keep the genetic claim at
locus/causal level and the tissue claim at expression-marker level.

---

## 5. Outputs

```
06_tissue_integration/
├── TaskF_01_backbone_tissue_inhouse_v1.csv   (TPM per sample + log2FC/padj, 3 genes)
├── TaskF_02_genetic_tissue_integration_v1.csv (master integration table)
└── TaskF_03_directionality_note.md            (paradox resolution text)
scripts/taskF_01_tissue_confirm.R
```

---

## 6. Interpretation guide

- **TSHR**: 3-layer convergence (MR + coloc + tissue) — the strongest integrated signal.
  Tissue UP confirms locus is active in target tissue; paradox resolved as above.
- **IGF1R**: MR risk + NO coloc + tissue UP-but-NS → consistent effector picture
  (active/upregulated in tissue, drug target, but not expression-driven susceptibility).
- **CTLA4**: positive control; tissue likely low/NS in adipose (immune gene) — report as-is.

---

## 7. Critical reminders

1. Backbone genes ONLY. No transcriptome-wide expansion, no candidate-gene tissue dump.
2. n=1 control → tissue is supporting/descriptive, not powered DE. State honestly.
3. In-house = DESeq2; GEO = limma. Never mix.
4. Methods must be distinguishable from the separate descriptive RNA-seq paper.
5. If TSHR log2FC does not reproduce ≈+2.33 → STOP, report (no silent acceptance).
6. Directionality paradox: locus-level framing only; never claim orbital expression protective.
7. GEO results (GSE58331/105149) are REUSED, not re-run. Cite Rosenbaum 2015/2017 correctly.

---

## 8. Verified citations

| Reference | Citation | PMID |
|-----------|----------|------|
| GSE105149 | GEO dataset, lacrimal gland (4 TED + 7 control) | **verify primary citation at GEO page before manuscript** |
| DESeq2 | Love et al. 2014 Genome Biol 15:550 | 25516281 |
| limma | Ritchie et al. 2015 NAR 43:e47 | 25605792 |

**Citation note:** GSE58331 PMID 26371757 (Rosenbaum 2015 PLoS One) is confirmed. GSE105149
(lacrimal gland, 4 TED + 7 control) primary-source citation must be verified directly from
its GEO record before manuscript submission — multiple downstream papers reuse it, so the
original-source attribution needs confirmation (do NOT cite from memory as 28910444 without
checking GEO).

---

*End of Week 4 Task Spec v1.0.*
*Integrates Weeks 1-3 genetic backbone with orbital tissue; locus-level prioritization framing.*
