# Loop engineering plan — TED-TRAP to publication

_Goal-driven, gated, iterative plan to take the portfolio to publication with
minimum wasted cycles. "Loop" = a closed feedback cycle with a defined goal,
inputs, gate (exit criterion), owner, and artifact. Guiding value: 정확하고
진실한게 생명 — never advance a loop on unverified data._

## North-star goal
Publish the portfolio honestly and efficiently:
- **G1:** Paper 1 (genetics framework) published in an Endocrine-tier journal.
- **G2:** Paper 2 (INSR off-target) published in a lower/niche tier — *only if its
  snRNA-seq gate passes*; otherwise downgraded or shelved (a legitimate outcome).
- **G3 (optional):** Paper 3 (GD↔TED genetic relationship, or pQTL novelty) — only
  if a clean question survives scoping.

## Owners / roles (who closes each step)
- **PI (human):** anything outward-facing or judgement-bound — journal submission,
  ORCID/portal, IRB, co-author sign-off, GitHub repo recreate, final go/no-go.
- **Claude Code (this repo):** text edits, consistency audits, docx build, figure
  verification, literature checks, plan/issue tracking, PR/CI monitoring, drafting.
- **Antigravity/Gemini (R execution):** runs R that needs the local licensed data
  (pQTL-MR, snRNA-seq DE) and reports console output back. Never edits manuscripts.

---

## LOOP 1 — Paper 1 submission cycle (ACTIVE)
**Goal:** acceptance at an Endocrine-tier journal. **State machine:**

```
READY ─submit→ UNDER REVIEW ─decision→ { ACCEPT → done
                                         MINOR/MAJOR → REVISE loop → resubmit
                                         REJECT → next venue (cascade) }
```

| Step | Owner | Gate / exit criterion | Artifact |
|---|---|---|---|
| Pre-flight | Claude | md5 stable, 0 placeholders, audits pass, docx valid | this is **done** (md5 a112b854) |
| Author tasks | PI | ORCID, portal fields, figures eyeballed ≥300 DPI, **repo recreated (security)** | submit-ready |
| Submit | PI | Endocrine Connections (1st), else Endocrine/JEI | submission ID |
| Watch | Claude | (after a PR/preprint exists) monitor; for journal, PI relays decision | status log |
| Triage decision | Claude+PI | classify reviewer points: accept-as-is / quick-fix / needs-analysis | `internal/REVIEW_RESPONSE_*.md` |
| Revise | Claude | every reviewer point answered; new numbers verified before insert | point-by-point + revised master |
| Re-audit | Claude | full consistency audit re-run; md5 re-locked | updated package |

**Cascade rule (no idle time between rejections):** Endocrine Connections → Endocrine
→ JEI → Frontiers-Endo/BMC. Do **not** re-shop to a thyroid-specialist journal that
already desk-rejected. (See `SUBMISSION_VENUE_CASCADE.md`.)

**Pre-staged rebuttal ammunition (already in repo):**
- "external replication?" → `INTERNAL_external_GEO_research_note.md` (Option A rationale).
- "how is this different from 2025 GD/GO MR papers?" → `LITERATURE_LANDSCAPE_CHECK.md`
  + (optional) the differentiation sentence parked for Discussion.
- "n=1 tissue?" → exploratory framing + (if pushed) public snRNA-seq from Loop 2.

---

## LOOP 2 — Paper 2 build-then-submit cycle (GATED, not started)
**Goal:** a defensible off-target paper, OR an honest decision to shelve.
**Hard gate first — do not write the paper before the gate.**

```
P2_00 audit ─→ P2_01 pQTL-MR (expect null, supportive) 
            └→ P2_02 public snRNA-seq direction check  ── GATE A ──┐
GATE A: INSR/insulin concordant in active orbital snRNA-seq?       │
   PASS → draft Paper 2 (Frontiers/Sci Rep) ───────────────────────┘
   WEAK/NULL → downgrade to hypothesis/correspondence OR fold into the
               (submitted) RNA-seq paper's revision OR shelve.
```

| Step | Owner | Gate / exit | Artifact |
|---|---|---|---|
| Data audit | Antigravity | `P2_00` reports pQTL+outcomes+LD+snRNA-seq availability | `P2_00_audit_manifest.csv` |
| pQTL-MR | Antigravity | `P2_01` returns INSR/IGF1R MR+coloc (null = supportive) | `P2_01_*.csv` |
| **GATE A** snRNA-seq | Antigravity+Claude | INSR direction concordant in active TED orbital snRNA-seq | `P2_02` result |
| Draft (if PASS) | Claude | own master, own md5, anti-salami vs RNA-seq paper | Paper 2 master |
| Verify Kim 2021 IOVS S2 | PI/Claude | INSR novelty wording justified before any claim | note |

**Honest prior:** pQTL-MR for GD/GO is already published; INSR off-target is known
pharmacology → Paper 2 is *incremental*, lower-tier, gate-dependent. Treat shelving
as a valid, non-failure outcome.

---

## LOOP 3 — Optional novelty bet (only if a clean question survives)
- **3a GD↔TED genetic relationship:** LDSC genetic correlation + bidirectional MR
  (public GWAS only). Answers "is TED genetically independent of GD or shared?" —
  a genuinely open question this project is positioned to ask.
- **3b pQTL novelty:** structurally the only route to a *new target*, but the
  GD/GO pQTL-MR space is now crowded → require a differentiated angle before
  investing. Run before claiming; a null is real.
- Gate: scope it in `internal/` first; do not start analysis until the question is
  shown to be unoccupied and tractable.

---

## Cadence & automation hooks (available in this environment)
- **`/loop`** (skill): schedule a recurring check (e.g. weekly: "re-audit master
  integrity + summarize open loops + flag any stale md5"). Good for keeping the
  package submission-ready without manual nudges.
- **`subscribe_pr_activity`** / **Monitor:** once any GitHub PR or CI exists, watch
  it and auto-triage events. (Not applicable until a PR is opened.)
- **Background Agent:** fan-out literature re-checks or data audits without blocking.
- Not available here: `/goal`, `ultracode`, a "Workflow" engine — so loops are run
  as the above primitives + this document, not a dedicated orchestrator.

## Definition of done (per loop)
- Loop 1: Paper 1 accepted (or cascade exhausted with a documented reason).
- Loop 2: Paper 2 accepted, OR a recorded decision to downgrade/shelve after GATE A.
- Loop 3: a scoping note concluding go or no-go.

## Cross-loop invariants (never violated)
- One master per paper; md5 is the integrity anchor; never edit a stale fork.
- No unverified number enters any manuscript (Gemini-data lesson).
- Scope firewall between Paper 1 / Paper 2 / RNA-seq paper (anti-salami).
- IRB/patent/licensed data never committed.
