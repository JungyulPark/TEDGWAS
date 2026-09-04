# SUBMISSION HOLD — provenance audit of external review round 5

_Date: 2026-08. Status: **제출 보류.** Every claim below was checked against the
source files in this repository, not against memory. Where I could not verify
something, that is stated explicitly._

An external reviewer raised eight objections and judged the manuscript
submission-stopping. **Five are confirmed true against the source data, one is
true and worse than diagnosed, one is partly misattributed, and one could not be
verified from this container.** The manuscript must not be submitted in its
current state.

---

## Confirmed 1 — one instrument violates the stated selection threshold

**Reviewer's argument.** Methods says all primary instruments have
*P* < 5×10⁻⁸. With *F* = *Z*², that forces |*Z*| > 5.45 and therefore
*F* > 29.7. A minimum *F* of 14.2 is arithmetically impossible.

**Verified.** The argument is correct, and the cause is a single contaminating row.

Checked `TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv`
(6,136 rows, 2,545 genes):

| check | result |
|---|---|
| instruments with *P* ≥ 5×10⁻⁸ | **1** |
| instruments with \|*Z*\| < 5.4513 | **1** |
| the offending row | `GAN` / `rs310019`, *Z* = −3.7684, *P* = 1.643×10⁻⁴, *N* = 22,151 |
| its *F* = *Z*² | **14.2008** ← this is the reported minimum |
| next smallest *F* | 29.724 (`PPIB` / `rs895885`, *P* = 4.986×10⁻⁸) |

`GAN` has exactly one instrument, and it is that SNP. So `GAN` should never have
entered the MR-testable set at all. It fails the relaxed *P* < 5×10⁻⁶ threshold too.

**Corrected figures.**

| reported | correct |
|---|---|
| 6,136 SNP-level instruments | **6,135** |
| 2,545 MR-testable genes | **2,544** |
| minimum *F* = 14.2 | **29.72** |
| 2,235 BBJ-estimable genes | **2,234** (`GAN` is in the 2,235) |
| Bonferroni 0.05/2,545 = 1.9646×10⁻⁵ | 0.05/2,544 = 1.9654×10⁻⁵ |

**Impact on conclusions.** Both thresholds print as 1.965×10⁻⁵ at four significant
figures, and `GAN` is neither a discovery hit nor a backbone gene, so **no hit
changes and no backbone estimate changes.** But every count in the funnel, in
Figure 1 Panel A, in Table S4's footnote and in the Results paragraph is off by
one, and the reported minimum *F* is wrong by a factor of two. The screening
denominator is load-bearing for the Bonferroni claim, so this must be corrected,
not waved off as a typo.

**Required action.** Rebuild the SNP-level manifest with an assertion that every
primary instrument satisfies *P* < 5×10⁻⁸ **and** *F* = *Z*² ≥ 29.7, then
propagate 6,135 / 2,544 / 2,234 / 29.72 through the manuscript, Table S4,
Table S8 (the power calculation uses the per-gene set) and Figure 1.

---

## Confirmed 2 — Table S3's fine-mapping block is not reproducible (most serious)

**Reviewer's argument.** Table S3 lists a single credible set whose five PIPs sum
to 0.987, none of them rs179252, while the text below claims rs179252 has PIP 1.0.
Those cannot both hold.

**Verified — and the problem is larger than an internal inconsistency.**

Source of record, `TrackA_MR/v5_upgrade/07_manuscript/supp_tables/TableS3_finemapping.csv`:

```
European (1000G),CS_1,1,1,rs179252,1,81435985,2.86e-40,1,1
European (1000G),CS_2,1,1,rs179248,1,81433038,1.64e-39,0.965,0.992
East Asian (1000G),CS_1,1,1,rs3783950,1,81448282,7.75e-23,0.229,0.808
East Asian (1000G),CS_2,1,1,rs3783949,1,81448382,2.38e-20,0.207,0.808
```

i.e. SuSiE returned **two credible sets in the European panel**, each containing
**one SNP**, each with **purity = 1**. `TaskB_05_tshr_finemap_summary_v1.csv`
agrees (rows TSHR_S1 and TSHR_S19: `susie_pip_eur = 1` for both rs179252 and
rs179248). `WEEK1_MANIFEST.md` states "SuSiE EUR credible sets | **2**, both
collinear with rs179252 (r² ≥ 0.965)", and `TaskC_decision_table_v1.csv` records
"TSHR independent signals (SuSiE EUR count), **2**".

The manuscript instead asserts:

| manuscript | source |
|---|---|
| "a **single** 95% credible set" | **two** credible sets (EUR) |
| "purity = 0.993" | purity = **1** |
| "log₁₀ Bayes factor = 88.2" | **not present in any file in the repository** |
| five variants, PIP 0.458 / 0.436 / 0.052 / 0.021 / 0.020 | each CS contains **one** SNP with PIP **1.0** |
| lead variant `rs11603529` | **appears in no data file in this repository** |
| `rs10137255` | **appears in no data file in this repository** |
| "~116-kb window" | derived from that unverified variant list |

Greps for `rs11603529`, `rs10137255`, `0.993` and `88.2` across every `.csv`,
`.tsv` and `.txt` in the repo return no fine-mapping source; `0.993` and `88.2`
occur only in unrelated `GEO_validation/` limma output.

The **clumping half** of Table S3 is correct — it matches
`TaskB_02_ld_clumping_sensitivity_v1.csv` row for row.

**What survives.** The *conclusion* is still defensible on the real data: both
European credible sets are single-SNP, and rs179248 sits at r² = 0.965 with
rs179252, so neither yields an independent instrument
(`usable_as_independent_iv = FALSE` for every row of TaskB_05). The single-IV,
locus-level treatment of *TSHR* stands. **The reported numbers do not.**

**What must be withdrawn now.** "Fine-mapping resolved the *TSHR* signal into a
single 95% credible set (purity 0.993)" and the five-variant PIP table. Also the
Methods sentence and Figure S1 legend must be reconciled — Methods says
"secondary credible-set lead variants" (plural, consistent with two CSs) while
Table S3 says "single credible set, no second signal". They contradict each other.

**Required action.** Re-run SuSiE and report, from the actual run: input SNP
count, LD reference and build, credible-set membership, every PIP, coverage and
purity, and rs179252's PIP — plus a `susie_rss` LD-vs-summary-statistics
consistency diagnostic, since reference-LD mismatch destabilises `susie_rss`.

---

## Confirmed 3 — the colocalization input is p-value based, with one European MAF vector for both outcomes

**Verified from `TrackA_MR/v5_upgrade/scripts/taskE_01_coloc.R`:**

```r
coloc.abf(
  dataset1 = list(snp=m$SNP, pvalues=m$p_gwas,  type="cc",    s=s_val, N=n_val),
  dataset2 = list(snp=m$SNP, pvalues=m$p_eqtl,  type="quant", N=N_EQTL),
  MAF = m$MAF)
```

Three consequences, all as the reviewer described:

1. **No beta/varbeta.** ABF is computed from *p*-values, so **MAF, N and the case
   fraction `s` are load-bearing inputs that directly move the posterior.**
   They are not incidental metadata.
2. **A single MAF vector** — `04_druggable_mr/data/g1000_eur_freq.frq`, 1000G
   **European** — is applied to both datasets and to **both outcomes, including
   the East Asian BBJ discovery GWAS.** PP.H4 = 0.951 for *TSHR* in BBJ therefore
   rests on European allele frequencies applied to an East Asian outcome. The
   manuscript acknowledges this in Limitations, but reports 0.951 in the abstract,
   Table 2 and Figure 3 as a settled value.
3. **`SNP.PP.H4` was computed** (it selects `top_snp`) but never reported. It is
   recoverable and should be.

**Note on the reviewer's eQTLGen remark.** They wrote that the manuscript claims
"eQTLGen reports no allele frequencies". **That sentence is not in our
manuscript** — ours says exposure allele frequencies "were approximated from a
European reference panel". The underlying criticism nevertheless holds: eQTLGen
does publish its own allele frequencies (n = 26,609), so the approximation was
avoidable.

---

## Confirmed 4 — UK Biobank colocalization was never run

`TaskE_01_coloc_results_v1.csv` contains exactly two outcomes: `BBJ_Graves` and
`FinnGen_GO`. The `OUTCOMES` list in `taskE_01_coloc.R` has only those two.

This is the reviewer's strongest structural point. UKB is the outcome whose
ancestry best matches eQTLGen, and it carries *IGF1R*'s strongest MR signal
(*P* = 0.012). Omitting it while reporting colocalization for the East Asian
discovery arm and the 858-case sensitivity arm is exactly the pattern a reviewer
reads as selective analysis.

---

## Confirmed 5 — RESOLVED by Task 2: the case count is right, the denominator is wrong

_Status updated 2026-09-04 after the local session's Task 2 (commit `a7e5cbc`).
Every claim below was re-verified in this container against the repository._

**The reviewer was wrong on this point, and the manuscript's "858 cases" stands.**

### What Task 2 established

The FinnGen file carries `af_alt`, `af_alt_cases` and `af_alt_controls`, so the
case fraction can be recovered from the file itself with no external dependency:

```
af_alt = s·af_cases + (1−s)·af_controls   ⇒   s = (af_alt − af_ctrl)/(af_case − af_ctrl)
```

Re-checked here: the identity recovers `s` to within 5×10⁻¹⁷. Applied to the three
FinnGen files it yields **858.01 / 752.98 / 3962.02** cases at N = 500,348 — integers.

Those exact numbers were **already in the repository**, in
`TrackA_MR/scripts/03b_04v3_finngen_local.R:50-54`, verified here:

```r
TED_main    = list(N = 500348, cases = 858,  controls = 499490, label = "GRAVES_OPHT"),
GD_exposure = list(N = 500348, cases = 3962, controls = 496386, label = "E4_GRAVES_STRICT"),
TED_strict  = list(N = 500348, cases = 753,  controls = 499595, label = "E4_GRAVES_OPHT_STRICT")
```

The local session additionally confirmed against FinnGen R12 PheWeb (rendered in a
browser; PheWeb is a JS app and does not serve content to a plain fetch):
`GRAVES_OPHT` = **858 cases / 499,490 controls**. That web check is theirs, not mine;
the three internal lines of evidence above are what I verified directly.

### The actual defect

| | value |
|---|---|
| file used by v5 | `finngen_R12_GRAVES_OPHT.gz` (**broader** endpoint), md5 `1b26a77e1e4eb451e80ead00e4ec6f28` |
| endpoint / release | `GRAVES_OPHT` / R12 |
| cases / controls / N | 858 / 499,490 / **500,348** |
| `taskE_01_coloc.R` used | `858 / **520,387**` ← wrong denominator |
| provenance of 520,387 | **none.** `analytical_log.md:142` records "Total N used: 520,387" with no source. The placeholder was never replaced. |

**The correct value was in the repository the whole time; `taskE_01_coloc.R` wrote a
placeholder instead of referencing it.** That is the defect — a provenance failure,
not an arithmetic one.

### Numerical impact: negligible, and here is why

`coloc:::Var.data.cc = 1/(2·N·f·(1−f)·s·(1−s))`, so N and s enter **only as the
product** `N·s·(1−s)`. Both parameterisations encode the same numerator
(`N·s = 858` exactly, since s was defined as 858/N in each), so only the `(1−s)`
factor differs. Re-computed here:

```
correct : N=500,348  s=0.001714806  N·s(1−s) = 856.5287
used    : N=520,387  s=0.001648773  N·s(1−s) = 856.5854   → 0.0066% apart
lABF difference at z = 2 / 5 / 8 : 2.2e-05 / 2.5e-05 / 1.1e-04
```

**PP.H4 = 0.986 does not move from this defect.** It says nothing, however, about
the *other* Task 4 defects (p-value-based ABF, a European MAF vector on the East
Asian outcome, and the missing UKB arm), which remain open.

### Why the reviewer's 786 / 894 were wrong

They are **Risteys** figures, which count a different stage from the **PheWeb**
summary statistics that were actually analysed:

| endpoint | Risteys R12 | PheWeb R12 = the analysed GWAS |
|---|---|---|
| `GRAVES_OPHT` | 894 | **858 cases** |
| `E4_GRAVES_OPHT_STRICT` | 786 | **753 cases** |

The β, SE and p-values in the file come from the 858 analysis set. **No change is
needed to Table 1 or to any "858 cases" statement in the manuscript.**

### New defect found during Task 2

`TrackA_MR/v5_upgrade/00_meta/data_sources.csv` row `FinnGen_R12_GO_v1` points at the
**broader** file while labelling it with the **strict** endpoint name
(`E4_GRAVES_OPHT_STRICT`) and giving controls as **376,419** (actual: 499,490); it
also still carries a literal `TBD_by_Gemini` field. Verified here. No evidence it fed
any calculation, but it is a reproducibility document and must be corrected.

### Remaining action

- Set `FINNGEN_S <- 858/500348`, `FINNGEN_N <- 500348` in the Task 4 rewrite.
- Fix the `data_sources.csv` row.
- Delete the stale `520,387` from `analytical_log.md` or annotate it as superseded.

## Confirmed 6 — chr16p11.2 is not established as a single LD signal

Table S5's own footnote reports rs4889606 (*HSD3B7*) and rs34649473 (*VKORC1*) at
r² = 0.90 — one signal — but rs78924645 (*PRSS36*) at **r² = 0.19 and 0.20** with
those two. The footnote then concludes the r² values "confirm that the cluster
does not comprise three independent signals". **That does not follow from
r² = 0.19.** Physical proximity within 143 kb is not evidence of a shared signal.

The defensible statement is two signals, not one: *HSD3B7*/*VKORC1* are collinear;
*PRSS36* is a separate, weakly correlated variant that happens to lie nearby.
Establishing otherwise requires conditional analysis or conditional/multi-variant
colocalization.

This changes the Figure 1 attrition tally (currently "3 — chr16p11.2, one LD block").

---

## Confirmed 7 — "no robust novel target" over-reaches

The manuscript concedes FinnGen is underpowered, that a high PP.H2 there is
consistent with absent outcome signal rather than distinct architecture, and that
the phenotype is not a within-Graves TED contrast — then uses failure in that same
outcome as the disqualifying filter.

Checked in `TaskD_03d_MR_all_outcomes_merged_v1.csv`:

| gene | BBJ | UKB | FinnGen | BBJ coloc PP.H4 |
|---|---|---|---|---|
| *TNFSF14* | β = −0.468, *P* = 1.5×10⁻⁶ | β = −0.298, ***P* = 0.0056** | β = +0.034, *P* = 0.855 | 0.994 |
| *IFNGR1* | β = +0.741, *P* = 9.4×10⁻⁶ | β = +0.301, *P* = 0.070 | β = +0.145, *P* = 0.677 | 0.989 |

*TNFSF14* has discovery colocalization at 0.994 **and** direction-consistent
nominal replication in UKB. Calling that "no robust novel target" discards a
legitimate **Graves disease candidate requiring replication**.

**Defensible wording:** "No novel candidate satisfied our stringent dual-phenotype
colocalization criterion." Not "no robust novel target".

---

## Confirmed 8 — the single-control RNA-seq should leave the inferential chain

Four TED vs **one** control cannot estimate control-side biological variance;
DESeq2's own documentation declines inference without biological replication. The
manuscript already labels the analysis exploratory and non-confirmatory, but still
carries adjusted *P* = 0.032 into the abstract and lists tissue as a third
evidence layer for *TSHR*.

**Recommendation:** demote to descriptive illustration, remove the *P* value from
abstract and conclusion, or replicate in an external orbital dataset.

---

## Also flagged, not yet resolved

- Single instrument plus colocalization does not establish that *TSHR* **expression**
  is causal; the same variant may act through a neighbouring gene or another
  mechanism. Current wording is stronger than the design supports.
- *IGF1R*'s MR-Egger intercept is estimated from 3–4 instruments. "No evidence of
  directional pleiotropy" should read "no detectable evidence, with limited power".
- The UKB linear-mixed-model → log-odds rescaling is not specified: formula, total
  *N*, case prevalence, and whether both β and SE were transformed.
- Table S4's secondary-threshold counts are not labelled clumped vs raw.
- "Prespecified" is supported for the outcome hierarchy only by
  `TaskC_pre_analysis_plan_v1.md`, which carries no registration or timestamp.
  Safer: "defined analytical hierarchy".

---

## Verdict

| element | status |
|---|---|
| Research question and the *TSHR* vs *IGF1R* contrast | **worth keeping** |
| Instrument counts, minimum *F*, screening denominator | **wrong — one invalid instrument** |
| Table S3 fine-mapping numbers | **not reproducible from the repository** |
| "rs179252 resolved as the shared causal variant" | **withdraw** — eQTL-only fine-mapping cannot establish it |
| FinnGen endpoint, cases, `s` | **resolved (Task 2)** — 858 cases correct; denominator 520,387 wrong, should be 500,348; numerically negligible |
| BBJ PP.H4 = 0.951 | **conditional** on a European MAF vector applied to an East Asian GWAS |
| UKB colocalization | **never run** |
| "no robust novel target" | **over-stated** — *TNFSF14* replicates nominally in UKB |
| chr16p11.2 as one LD block | **not established** (r² = 0.19/0.20) |
| Orbital RNA-seq as an evidence layer | **remove from inference** |
| *TSHR* locus MR + colocalization, broad direction | **likely survives re-analysis** |

## Order of work — do not edit prose first

1. **Rebuild the instrument manifest** with automated assertions
   (`P < 5e-8` and `F = Z² ≥ 29.7` for every primary instrument). Repropagate
   6,135 / 2,544 / 2,234 / 29.72.
2. ~~**Pin the FinnGen file**~~ — **done** (Task 2, commit `a7e5cbc`):
   `GRAVES_OPHT` / R12 / 858 cases / 499,490 controls / N = 500,348 / s = 858/500348.
3. **Re-run SuSiE** and report the real credible sets, PIPs, coverage, purity and
   an LD-consistency diagnostic.
4. **Re-run colocalization for all three outcomes** with beta/varbeta where
   available or ancestry-matched MAF otherwise; report H0–H4, `SNP.PP.H4`, the
   outcome minimum *P* per locus, and prior sensitivity at
   p12 = 10⁻⁶, 5×10⁻⁶, 10⁻⁵.

Only after 1–4 pass should tables, figures and the abstract be regenerated.

**Blocking constraint:** steps 2–4 need eQTLGen, FinnGen, BBJ and 1000G files that
are licence-bound and not in this repository. The v5 scripts point at
`c:/ProjectTEDGWAS/`. These re-runs must happen on the machine that holds those
files; this container can prepare and verify the scripts but cannot execute them.

---

# Task 3 addendum — remote verification (2026-09-04)

The local session's Task 3 forensics were re-checked here. Three findings verified
independently, plus one analytical result that bounds how far the new bug reaches.

## Verified: the fabricated block was written at packaging time

`git log -S"rs11603529"` returns exactly one commit, and its diff is unambiguous.
Commit **`e72f84c` "Finalize JEI submission package"** replaced this —

```
[Per-credible-set PIP values from the SuSiE output to be tabulated in the final
supplementary file.]
```

— with the five-variant PIP table, "single 95% credible set (purity = 0.993)",
"log₁₀ Bayes factor = 88.2" and the "~116-kb window". Its own commit message states
the goal: *"S3 PIP table + S5 chr16 coords/r2 filled (placeholders now 0)"*.

**A placeholder was closed by inventing its contents rather than by reading the
SuSiE output.** This is the same failure mode as Task 2's FinnGen denominator, one
step worse: Task 2 produced a wrong number, this produced a whole results table.

The immediately preceding version said the credible-set leads fell within a
**~15-kb** window — which matches the real lead span the local session measured
(15.3 kb). The "~116-kb" figure is arithmetic on the two fabricated coordinates.

## Verified: the same commit's other placeholder was filled correctly

`e72f84c` also filled Table S5's chr16 coordinates. Those **are** genuine — all
three match `TaskD_02a_eqtl_instruments_snp_level_v1.csv` exactly:

```
HSD3B7  rs4889606   31,011,183      VKORC1  rs34649473  31,066,538
PRSS36  rs78924645  31,154,358
```

So the commit is not uniformly untrustworthy. The difference is that S5's values sat
in an accessible CSV while S3's PIPs required opening the SuSiE output.
**Still unverified from that commit: the S5 pairwise r² values (0.90 / 0.19 / 0.20).
Task 5 must confirm them from an LD computation, not assume them.**

## Verified: there is no allele harmonization in the SuSiE script

`TrackA_MR/v5_upgrade/scripts/taskB_04_susie_finemap.R`, checked line by line:

- **L62** `plink --r square --keep-allele-order` — produces **signed** correlations,
  oriented to the 1000G `.bim` A1.
- **L106** `sub <- sub[match(retained, snp)]` — reorders the z-vector to the LD
  matrix **by rsID only**.
- Nothing anywhere compares eQTLGen's `effect_allele` to the `.bim` A1. Greps for
  `allele`, `flip`, `harmon`, `strand`, `sign` return only those two lines.

So z is on eQTLGen's orientation and R is on the `.bim` orientation, with no
reconciliation. The local session's diagnosis is correct.

## New: how far the bug actually reaches

Sign errors in R do not propagate everywhere, and the boundary matters for what has
to be re-run.

| consumer | uses | affected? |
|---|---|---|
| PLINK **clumping** (Table S3 clumping half, `TaskB_02`) | **r²** — sign-invariant | **No.** This is why that half verified row-for-row against source. |
| `ld_with_rs179252_eur/eas` in `TaskB_05` (0.965, 0.808) | **r²** — sign-invariant | **No.** |
| `susie_rss` (credible sets, PIPs) | **signed R** | **Yes.** |
| `coloc.abf` under the single-causal-variant assumption | **no LD matrix at all** — only per-SNP p-values, MAF, N, s | **No.** PP.H4 = 0.951 / 0.986 are not contaminated by this bug. |
| `coloc.susie` (proposed as a Task 4 option) | **signed LD** | **Yes — would inherit it.** Harmonize before running, or do not run it. |

Two consequences:

1. **Task 4's `coloc.abf` re-run is not blocked by Task 3.** The colocalization
   defects (p-value-based ABF, European MAF on the East Asian outcome, missing UKB
   arm) are independent of the allele bug and can proceed in parallel.
2. **The single-IV decision survives either way.** It rests on r² — every other lead
   at the locus sits at r² ≥ 0.808 (EAS) / 0.965 (EUR) with rs179252 — and r² is
   sign-invariant. So whether harmonization merges the two European credible sets
   into one or leaves them at two, *TSHR* remains a single-instrument, locus-level
   anchor. What changes is the **stated** fine-mapping result, not the instrument
   choice built on it.

## Still open, pending the local session's harmonization test

- Whether the two European credible sets collapse to one after harmonization.
  The local session found rs179252 (mismatched) and rs179248 (matched) sit in
  *opposite* orientation groups, so the signed LD between the two credible-set leads
  is the pair most likely to be wrong. Neither outcome is being claimed until the
  test returns.
- `rs179252` carrying **EAS PIP = 0** while the discovery outcome is East Asian.
  This needs an explanation in the manuscript regardless of how the harmonization
  test resolves, and the EAS run did not converge (`converged = FALSE`,
  `estimate_s_rss` 0.868 vs 0.763 in EUR).
