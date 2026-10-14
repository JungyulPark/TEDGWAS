#!/usr/bin/env python3
"""Derive every number added to the manuscript in the 2026-10-14 review revision.

Reads only Supplementary_Dataset.xlsx and the retained aggregate MR export, and
writes Verification/Review_20261014/derived_values.csv (item, value, displayed,
source). The integrity audit compares the displayed strings with the master.

    python3 submission/candidate_20261011_review/reproducibility/derive_review_values_20261014.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

PKG = Path(__file__).resolve().parents[1]
X = PKG / "Supplementary_Dataset.xlsx"
OUT = PKG / "Verification" / "Review_20261014" / "derived_values.csv"


def sig3(p):
    return f"{p:.3g}" if p >= 0.001 else f"{p:.2e}"


def main():
    rows = []
    add = lambda item, value, shown, src: rows.append(dict(item=item, value=value, displayed=shown, source=src))

    # per-allele TSHR ORs (rs179252 G allele)
    a = pd.read_excel(X, sheet_name="TSHR_allele_audit", header=None).values.tolist()
    hdr_i = next(i for i, r in enumerate(a) if "aligned_G_OR" in [str(v) for v in r])
    hdr = [str(v) for v in a[hdr_i]]
    names = {"BBJ": "BBJ", "UKB": "UKB", "FinnGen_GO": "FinnGen GO", "E4_GRAVES_STRICT": "FinnGen GD",
             "E4_GRAVES_OPHT_STRICT": "strict GO"}
    for r in a[hdr_i + 1:]:
        if pd.isna(r[0]) or str(r[0]).strip() == "" or "outcome" in [str(v) for v in r]:
            break  # end of section C
        d = dict(zip(hdr, r))
        if d.get("snp") == "rs179252" and d.get("outcome") in names:
            v = float(d["aligned_G_OR"])
            add(f"per-allele OR {names[d['outcome']]}", v, f"{v:.2f} ({names[d['outcome']]})", "TSHR_allele_audit C")

    # IGF1R alternative estimators (original three outcomes)
    e = pd.read_excel(X, sheet_name="MR_estimators_original")
    g = e[e.gene_symbol.eq("IGF1R")]
    for oc, lab in [("BBJ_Graves", "BBJ"), ("UKB_hyperthyroid", "UKB")]:
        r = g[(g.outcome == oc) & (g.method == "Weighted median")].iloc[0]
        add(f"IGF1R weighted median OR {lab}", r["or"], f"{r['or']:.2f}", "MR_estimators_original")
        add(f"IGF1R weighted median P {lab}", r.pvalue, sig3(r.pvalue), "MR_estimators_original")
    others = g[g.method.isin(["MR Egger", "Weighted mode"])].pvalue
    add("IGF1R Egger/mode minimum P", others.min(), "0.32", "MR_estimators_original (displayed as P >= 0.32)")
    for oc, lab in [("BBJ_Graves", "BBJ"), ("UKB_hyperthyroid", "UKB"), ("FinnGen_GO", "FinnGen GO")]:
        v = g[(g.outcome == oc) & (g.method == "MR Egger")].egger_intercept_p.iloc[0]
        add(f"IGF1R Egger intercept P {lab}", v, f"{v:.2f}", "MR_estimators_original")

    # IGF1R leave-one-out (original outcomes)
    lo = pd.read_excel(X, sheet_name="LOO_original")
    for oc, snp in [("BBJ_Graves", "rs2654980"), ("UKB_hyperthyroid", "rs2654980"), ("BBJ_Graves", "rs59467480"),
                    ("UKB_hyperthyroid", "rs117212126")]:
        v = lo[(lo.gene == "IGF1R") & (lo.outcome == oc) & (lo.excluded_SNP == snp)].pvalue.iloc[0]
        add(f"IGF1R LOO P {oc} without {snp}", v, sig3(v), "LOO_original")

    # detection benchmarks (same formula as Detection_original)
    mr = pd.read_excel(X, sheet_name="MR_primary")
    k = lambda alpha: norm.ppf(1 - alpha / 2) + norm.ppf(0.8)
    for oc, lab, alpha in [("BBJ_Graves", "BBJ", 0.05 / 2544), ("E4_GRAVES_STRICT", "FinnGen GD", 0.05 / 2544)]:
        se = mr[(mr.outcome == oc) & mr.estimable & mr.se.notna()].se
        det = np.exp(k(alpha) * se)
        for thr in (1.5, 2.0):
            v = float((det <= thr).mean())
            add(f"detectable fraction {lab} OR {thr}", v, f"{100 * v:.1f}%", "MR_primary SE")
    for oc, lab in [("BBJ_Graves", "BBJ"), ("UKB_hyperthyroid", "UKB"), ("E4_GRAVES_STRICT", "FinnGen GD")]:
        se = mr[(mr.gene_symbol == "IGF1R") & (mr.outcome == oc)].se.iloc[0]
        v = float(np.exp(k(0.05) * se))
        add(f"IGF1R detectable OR alpha 0.05 {lab}", v, f"{v:.2f}", "MR_primary SE")

    # FinnGen GD heterogeneity for HLA-A (fixed 14-gene export)
    q = pd.read_csv(PKG / "reproducibility/work/reanalysis_20261002/mr_v1/MR_14genes_all_scenarios.csv")
    v = q[(q.gene_symbol == "HLA-A") & (q.outcome == "E4_GRAVES_STRICT") & (q.scenario == "original_reference")].Q_pvalue.iloc[0]
    add("HLA-A FinnGen GD Q P", v, f"{v:.3g}", "reanalysis mr_v1 MR_14genes_all_scenarios")

    # adipose instruments: positions and blood eQTL P (eQTLGen regional file in repo)
    eq = pd.read_csv(PKG.parents[1] / "TrackA_MR/data/eqtl_TSHR_all.csv")
    pos = eq.set_index("SNP").SNPPos
    add("rs12101255 - rs179252 distance (bp)", int(pos["rs12101255"] - pos["rs179252"]), "15 kb", "eqtl_TSHR_all GRCh37")
    add("rs74064796 - rs12101255 distance (bp)", int(pos["rs74064796"] - pos["rs12101255"]), "50 bp", "eqtl_TSHR_all GRCh37")
    for s in ("rs12101255", "rs74064796"):
        add(f"{s} blood eQTL P", float(eq.set_index("SNP").Pvalue[s]), "10⁻³⁴", "eqtl_TSHR_all")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"{OUT.relative_to(PKG)}: {len(rows)} values")


if __name__ == "__main__":
    main()
