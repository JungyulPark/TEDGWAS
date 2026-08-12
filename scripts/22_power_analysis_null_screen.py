#!/usr/bin/env python3
"""
Power analysis for the druggable-gene-wide MR screen (TED-TRAP v5).

Purpose: make the informative null (robust_novel = 0) INTERPRETABLE by quantifying
what effect sizes the screen was powered to detect.

Method (design-based, uses only the locked v5 results — no new data, no assumptions
about case/control split): for a Wald/IVW estimate with standard error SE, the
minimum |beta| detectable with power 1-gamma at two-sided level alpha is
    beta_min = (z_{1-alpha/2} + z_{1-gamma}) * SE
This is exact for the normal-approximation test actually used, and is preferable to
reconstructing R^2/K because every SE is an empirical output of the harmonised analysis.

Input : TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv
Output: TrackA_MR/v5_upgrade/08_power/TaskH_01_power_summary_v1.csv
        TrackA_MR/v5_upgrade/08_power/TaskH_02_power_per_gene_v1.csv
"""
import csv, math, os, statistics as st
from scipy.stats import norm

SRC = "TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv"
OUT = "TrackA_MR/v5_upgrade/08_power"
os.makedirs(OUT, exist_ok=True)

BONF = 0.05 / 2545          # study-wide discovery threshold = 1.965e-5
NOMINAL = 0.05              # replication / TED-specific support threshold
POWER = 0.80
z_pow = norm.ppf(POWER)     # 0.8416

rows = [r for r in csv.DictReader(open(SRC))
        if r["method"] in ("Inverse variance weighted", "Wald ratio")
        and r.get("se") not in (None, "", "NA")]

def zcrit(alpha):
    return norm.ppf(1 - alpha / 2)

def summarize(outcome, alpha, label):
    sub = [r for r in rows if r["outcome"] == outcome]
    zc = zcrit(alpha)
    mdb = []           # minimum detectable |beta| per gene
    for r in sub:
        se = float(r["se"])
        if se > 0 and math.isfinite(se):
            mdb.append((r["gene_symbol"], (zc + z_pow) * se))
    vals = sorted(v for _, v in mdb)
    n = len(vals)
    q = st.quantiles(vals, n=4)
    med = st.median(vals)
    # fraction of genes powered to detect given ORs (two-sided, |beta| = |ln OR|)
    fr = {}
    for orr in (1.2, 1.5, 2.0, 3.0):
        b = abs(math.log(orr))
        fr[orr] = sum(1 for v in vals if v <= b) / n
    return dict(outcome=outcome, threshold=label, alpha=alpha, n_genes=n,
                mdb_median=med, mdb_q1=q[0], mdb_q3=q[2],
                or_median=math.exp(med), or_q1=math.exp(q[0]), or_q3=math.exp(q[2]),
                frac_pow_OR1_2=fr[1.2], frac_pow_OR1_5=fr[1.5],
                frac_pow_OR2=fr[2.0], frac_pow_OR3=fr[3.0]), mdb

summaries, pergene = [], {}
for outcome, alpha, label in [
    ("BBJ_Graves", BONF, "study-wide Bonferroni (P<1.965e-5)"),
    ("UKB_hyperthyroid", NOMINAL, "nominal (P<0.05)"),
    ("FinnGen_GO", NOMINAL, "nominal (P<0.05)"),
]:
    s, mdb = summarize(outcome, alpha, label)
    summaries.append(s)
    pergene[outcome] = dict(mdb)

with open(f"{OUT}/TaskH_01_power_summary_v1.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summaries[0].keys()))
    w.writeheader(); w.writerows(summaries)

# per-gene detail for the backbone + the 13 discovery hits
HITS = ["TSHR","CTLA4","IGF1R","HLA-A","HLA-DQA2","C4A","HSD3B7","TUBB","VKORC1",
        "TNFSF14","PRSS36","MAPKAPK5","PSMB8","IFNGR1"]
with open(f"{OUT}/TaskH_02_power_per_gene_v1.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene_symbol","outcome","observed_beta","se",
                "min_detectable_beta_80pct","min_detectable_OR_80pct","powered_for_observed"])
    for r in rows:
        g = r["gene_symbol"]
        if g not in HITS: continue
        alpha = BONF if r["outcome"] == "BBJ_Graves" else NOMINAL
        se = float(r["se"]); b = float(r["beta"])
        m = (zcrit(alpha) + z_pow) * se
        w.writerow([g, r["outcome"], f"{b:+.3f}", f"{se:.4f}", f"{m:.3f}",
                    f"{math.exp(m):.2f}", "yes" if abs(b) >= m else "no"])

print("=== POWER SUMMARY (80% power) ===")
for s in summaries:
    print(f"\n{s['outcome']}  [{s['threshold']}]  n={s['n_genes']} genes")
    print(f"  min detectable |beta| : median {s['mdb_median']:.3f}  (IQR {s['mdb_q1']:.3f}-{s['mdb_q3']:.3f})")
    print(f"  equivalent OR         : median {s['or_median']:.2f}  (IQR {s['or_q1']:.2f}-{s['or_q3']:.2f})")
    print(f"  genes powered for OR>=1.2 : {s['frac_pow_OR1_2']*100:.1f}%")
    print(f"  genes powered for OR>=1.5 : {s['frac_pow_OR1_5']*100:.1f}%")
    print(f"  genes powered for OR>=2.0 : {s['frac_pow_OR2']*100:.1f}%")
    print(f"  genes powered for OR>=3.0 : {s['frac_pow_OR3']*100:.1f}%")

print("\n=== BACKBONE genes: was the screen powered for the observed effect? ===")
for g in ("TSHR","IGF1R","CTLA4"):
    for r in rows:
        if r["gene_symbol"]==g and r["outcome"]=="BBJ_Graves":
            alpha=BONF; se=float(r["se"]); b=float(r["beta"])
            m=(zcrit(alpha)+z_pow)*se
            print(f"  {g:6} BBJ  beta={b:+.3f}  min-detectable={m:.3f} (OR {math.exp(m):.2f})  -> powered: {'YES' if abs(b)>=m else 'NO'}")
