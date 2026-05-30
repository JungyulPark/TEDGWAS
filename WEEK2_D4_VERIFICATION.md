# TED-TRAP v5 — Week 2 D.4 Hit Classification (VERIFIED by Claude)

**Date:** 2026-05-26
**Verified by:** Claude (direct analysis of uploaded TaskD_03d/02a/02b)
**Decision:** Option 1 — defensive paper (no overclaimed novel targets)

---

## 1. SUMMARY VERDICT

Of 2,545 instrumented druggable genes, **13 reached BBJ Bonferroni significance** (P < 1.965e-5). After applying defensive filters (MHC spillover, chr16p11.2 LD cluster, single-IV unverifiability, FinnGen GO direction consistency), **the number of defensible robust-novel druggable targets is effectively zero.** This is an honest result and does not weaken the paper: the core finding is the contrast between TSHR (protective, upstream) and IGF1R (risk, effector), supported within a systematic druggable-genome-wide screen.

---

## 2. THE 13 BONFERRONI HITS — FULL CLASSIFICATION

| Gene | chr:pos(Mb) | n_iv | BBJ β (P) | UKB rep | FinnGen dir | Class |
|------|-------------|------|-----------|---------|-------------|-------|
| HLA-A | 6:29.6 | 4 | +0.685 (2.6e-23) | no | — | MHC spillover |
| HLA-DQA2 | 6:32.5 | 1 | +0.870 (8.3e-17) | no | — | MHC spillover |
| **TSHR** | 14:81.4 | 1 | **−2.096 (1.1e-14)** | **yes** | **same** | **known_anchor ★** |
| **CTLA4** | 2:204.6 | 1 | **−1.740 (5.5e-15)** | **yes** | **same** | **known_anchor (pos. control) ★** |
| C4A | 6:32.0 | 1 | −0.803 (2.2e-10) | no | — | MHC spillover |
| HSD3B7 | 16:30.4 | 1 | −1.227 (2.0e-7) | yes | **opposite** | chr16p11.2 LD cluster |
| TUBB | 6:30.1 | 2 | +0.399 (3.1e-7) | yes | same | MHC spillover |
| VKORC1 | 16:31.1 | 1 | −2.017 (6.4e-7) | no (P=0.073) | NA | chr16p11.2 LD cluster |
| TNFSF14 | 19:6.7 | 1 | −0.468 (1.5e-6) | yes | **opposite (≈0)** | single-IV unverifiable |
| PRSS36 | 16:31.2 | 1 | +1.792 (2.0e-6) | yes | same | chr16p11.2 LD cluster |
| MAPKAPK5 | 12:112.0 | 4 | −0.023 (5.3e-6) | no | same | discovery_only |
| PSMB8 | 6:32.8 | 1 | +0.222 (6.9e-6) | no | no | MHC spillover |
| IFNGR1 | 6:137.2 | 1 | +0.741 (9.4e-6) | no (P=0.070) | same | single-IV unverifiable |

★ = the two genes that pass all criteria are both already-known GD loci.

---

## 3. WHY THE THREE "NOVEL" CANDIDATES FAILED

Gemini's preliminary run flagged TNFSF14, HSD3B7, PRSS36 as robust_novel. Direct verification rejects all three:

- **TNFSF14** (chr19): single-IV (Wald only, no pleiotropy test); FinnGen GO β = +0.034 (essentially null, opposite sign to BBJ) — fails TED-specific direction. Also an immune gene in a space already covered by published GD MR (TNF-superfamily); novelty claim weak.
- **HSD3B7** (chr16:30.4): single-IV; FinnGen GO β = +0.272 (opposite to BBJ −1.227) — direction inconsistent. Sits in chr16p11.2 with VKORC1/PRSS36 (~800 kb) — likely shared LD signal, not independent.
- **PRSS36** (chr16:31.2): single-IV; co-located with HSD3B7/VKORC1 in chr16p11.2 — same LD-cluster concern. Cannot be called independent without colocalization.

**chr16p11.2 cluster (HSD3B7/VKORC1/PRSS36)** behaves like a second MHC-type spillover: 3 genes in <800 kb sharing a likely common signal. Treating them as 3 discoveries would be a colocalization error.

---

## 4. WHAT IS DEFENSIBLE (the paper's real backbone)

| Gene | Role | Evidence | Why it holds |
|------|------|----------|--------------|
| **TSHR** | Upstream susceptibility, protective direction | BBJ −2.10, UKB −2.44, FinnGen −2.33; all genome-wide/strong; 3-outcome consistent | Most consistent signal in the entire screen |
| **IGF1R** | Effector, risk direction | BBJ +0.45 (P=0.021, multi-IV n=4), UKB +0.30 (P=0.012) | Cleanest replicated multi-IV signal; matches teprotumumab mechanism |
| **CTLA4** | Known autoimmune locus | 3-outcome consistent negative | Positive control — proves pipeline works |

**Central narrative:** Within a systematic druggable-genome-wide screen, TSHR (protective, upstream) and IGF1R (risk, effector) emerge as opposite-direction axes of GD/TED genetic architecture, consistent with TSHR-directed vs IGF1R-directed therapeutic strategies.

---

## 5. CRITICAL CAVEATS (must be in manuscript)

1. **Not the first druggable/immune-gene GD MR.** Prior work exists: medRxiv 2025.01.02.25319932 (pQTL→GD via immune cells); Sci Reports s41598-025-21754-4 (multi-omics MR, FGFRL1). Frame as "systematic druggable-genome-wide screen contextualizing TSHR vs IGF1R," NOT "first systematic GD MR."
2. **Cross-ancestry limitation.** eQTLGen instruments are EUR; BBJ discovery is EAS. Document explicitly.
3. **Single-IV genes (incl. TSHR) cannot undergo pleiotropy diagnostics.** Wald ratio only.
4. **MHC + chr16p11.2 hits excluded from novel claims** due to LD spillover; require conditional/coloc analysis (Week 3) before any independent-signal claim.
5. **rs179252 position note:** eQTLGen hg19 81.4 Mb vs BBJ harmonised 81.0 Mb region — rsID-matched; harmonised liftover (hm_code=10). Documented.

---

## 6. WEEK 3 HANDOFF (coloc candidates)

Colocalization in Week 3 should target:
- **TSHR** (confirm shared causal variant; v4.3 BBJ PP.H4=0.951 already strong)
- **IGF1R** (multi-IV; coloc would strengthen effector claim)
- **CTLA4** (positive control coloc)
- Optionally chr16p11.2 cluster (to formally show LD-spillover, i.e. NOT independent) — only if reviewers likely to ask.

robust_novel coloc: none required (none claimed).

---

*Verified by Claude against uploaded D.3 merged results + D.2 instrument files, 2026-05-26.*
*All 13 Bonferroni hits classified with chr/pos, n_iv, 3-outcome direction, MHC/LD flags.*
*Decision: Option 1 defensive paper. No novel-target overclaim.*
