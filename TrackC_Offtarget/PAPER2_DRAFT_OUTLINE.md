# Paper 2 (Track C) — draft outline / working asset

_Status: analysis starting. Separate master from Paper 1; do NOT edit the v5
master. Honest principle: a null is a real and acceptable result._

## Title (working)
Insulin-receptor and otologic off-target signatures of IGF-1R antagonism in
thyroid eye disease orbital tissue: a transcriptomic and genetic framework.

## Thesis (one line)
IGF-1R blockade (teprotumumab/linsitinib) engages INSR/insulin signaling and
hearing-loss–associated pathways — a transcriptomic and target-network basis for
its metabolic (hyperglycemia) and otologic (hearing-loss) adverse effects, and
**distinct from genetic susceptibility** (which INSR lacks).

## Key grounding facts (already in repo — no new analysis needed)
- **INSR is NOT a GD/TED genetic susceptibility signal.** Paper-1 blood eQTL-MR:
  INSR β = +0.025 (BBJ, P=0.93), +0.207 (UKB, P=0.56), +0.364 (FinnGen, P=0.46) —
  null in all three outcomes. → supports "off-target, not susceptibility".
- In-house orbital RNA-seq: INSR/insulin-cassette up-regulation (Track C results).
- Enrichment: IGF-1R-downstream ∩ KEGG insulin (51 genes) and ∩ hearing-loss
  genes (38) — but partly reflects known downstream sharing (report honestly).

## Aims & evidence status
| Aim | Question | Evidence now | New analysis needed (Antigravity R) |
|---|---|---|---|
| 1 | INSR/insulin tissue signature + external concordance | in-house DE (n=1 ctrl, fragile) | **P2_02**: public TED orbital snRNA-seq direction check |
| 2 | Genetic anchoring (is INSR a signal or pure off-target?) | INSR eQTL-MR null (Paper 1) | **P2_01**: pQTL-MR (UKB-PPP/deCODE) + coloc for INSR/IGF1R/insulin axis |
| 3 | Off-target liability map (insulin + otologic) | enrichment done | formalize universe + q (light, P2_03) |

## Manuscript skeleton
1. **Introduction** — teprotumumab efficacy + its metabolic/otologic AEs; IGF-1R/INSR
   shared signaling; gap: tissue- and genetics-level account of the off-target axis.
2. **Methods** — in-house RNA-seq DE (exploratory, biological-sample level); public
   snRNA-seq concordance; pQTL/eQTL-MR + colocalization (EUR-only LD, F=Z², coloc.abf);
   enrichment (hypergeometric, honest universe, BH-q).
3. **Results** —
   3.1 INSR/insulin-cassette up-regulation in TED orbital tissue (+ snRNA-seq concordance).
   3.2 Genetic dissociation: INSR shows no eQTL/pQTL-MR association or colocalization
       with GD/TED → off-target, not susceptibility.
   3.3 Off-target network: insulin and hearing-loss pathway overlap with IGF-1R targets.
4. **Discussion** — mechanistic hypothesis for hyperglycemia/hearing-loss AEs; honest
   limitations (n=1 control; enrichment tautology; INSR novelty "provisional"); relation
   to Paper 1 (effector axis → off-target consequence).

## Honest expectations (set before running)
- P2_01 (pQTL-MR) will **likely be null** for INSR — that is *consistent and supportive*
  (off-target, not driver), not a failure. It does NOT add discovery novelty.
- The paper's real strength = Aim 1 (tissue + snRNA-seq concordance) + Aim 3 (AE
  mechanism). If snRNA-seq concordance fails (Gate A), demote to in-house-exploratory
  or a hypothesis/perspective piece.

## Display items (planned)
- Fig 1: IGF-1R/INSR shared-signaling → AE phenotype concept.
- Fig 2: INSR/insulin cassette expression (in-house) + snRNA-seq concordance.
- Fig 3: eQTL + pQTL MR/coloc for INSR vs IGF1R across GD/TED (the dissociation).
- Fig 4: off-target liability network (insulin + hearing-loss overlap).

## Firewall
No content here enters the v5 master. Paper 2 gets its own master when drafting.
