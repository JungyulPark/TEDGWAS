#!/usr/bin/env python3
"""Step 1 of the provenance repair — verify and rebuild the instrument manifest.

Why this exists
---------------
Methods states that every primary instrument is a genome-wide significant
*cis*-eQTL (P < 5e-8), and Table S4 reports F = Z^2 with a minimum of 14.2.
Those two statements are arithmetically incompatible: P < 5e-8 implies
|Z| > 5.4513 and therefore F > 29.72. An external reviewer caught this.

The cause is a single contaminating row (GAN / rs310019, P = 1.643e-04), which
is also the *only* instrument that gene has. This script proves that claim,
removes the violation, and recomputes every count that depends on the manifest.

What it checks
--------------
1. F as published equals Z^2 (i.e. the F formula itself is not the bug).
2. Every primary instrument satisfies P < P_THRESHOLD.
3. Every primary instrument satisfies F = Z^2 >= F_FLOOR, where F_FLOOR is
   derived from P_THRESHOLD rather than hard-coded, so the two can never drift.
4. The published Table S4 contains exactly the same rows as the source manifest.
5. After removing violations, the discovery hit list is unchanged.

Outputs (written only with --write)
-----------------------------------
  TaskD_02a_eqtl_instruments_snp_level_v2_verified.csv
  TableS4_instruments_v2_verified.csv

Exit code 0 if the corrected manifest passes every assertion, 1 otherwise.

Usage:
  python3 scripts/27_verify_instrument_manifest.py           # report only
  python3 scripts/27_verify_instrument_manifest.py --write   # emit corrected files
"""
from __future__ import annotations

import csv
import math
import os
import sys
from statistics import NormalDist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SRC_MANIFEST = f"{ROOT}/TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v1.csv"
SRC_TABLES4 = f"{ROOT}/TrackA_MR/v5_upgrade/07_manuscript/supp_tables/TableS4_instruments.csv"
SRC_MR = f"{ROOT}/TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv"

OUT_MANIFEST = f"{ROOT}/TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v2_verified.csv"
OUT_TABLES4 = f"{ROOT}/TrackA_MR/v5_upgrade/07_manuscript/supp_tables/TableS4_instruments_v2_verified.csv"

P_THRESHOLD = 5e-8
BACKBONE = ("TSHR", "IGF1R", "CTLA4")
PRIMARY_METHODS = ("Inverse variance weighted", "Wald ratio")

# F floor implied by the P threshold, derived so the two cannot drift apart.
Z_FLOOR = NormalDist().inv_cdf(1 - P_THRESHOLD / 2)
F_FLOOR = Z_FLOOR**2

fails: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(f"  {'ok ' if ok else 'FAIL'}  {msg}")
    if not ok:
        fails.append(msg)


def main() -> int:
    write = "--write" in sys.argv

    print("== Step 1: instrument manifest verification ==")
    print(f"P threshold      : {P_THRESHOLD:.0e}")
    print(f"implied |Z| floor: {Z_FLOOR:.4f}")
    print(f"implied F floor  : {F_FLOOR:.2f}   (F = Z^2)\n")

    src = list(csv.DictReader(open(SRC_MANIFEST, encoding="utf-8")))
    pub = list(csv.DictReader(open(SRC_TABLES4, encoding="utf-8")))

    # --- 1. the published F column is genuinely Z^2 -------------------------
    f_mismatch = [
        r for r in pub
        if abs(float(r["F-statistic"]) - float(r["Z-score"]) ** 2)
        > 0.02 * max(1.0, float(r["Z-score"]) ** 2)
    ]
    check(not f_mismatch, f"published F equals Z^2 (mismatches: {len(f_mismatch)})")

    # --- 2. source and published table describe the same rows --------------
    key_src = {(r["gene_symbol"], r["snp"]) for r in src}
    key_pub = {(r["Gene Symbol"], r["SNP"]) for r in pub}
    check(key_src == key_pub,
          f"Table S4 matches the source manifest row for row "
          f"(source-only {len(key_src - key_pub)}, table-only {len(key_pub - key_src)})")

    # --- 3. threshold violations -------------------------------------------
    viol = []
    for r in src:
        p = float(r["pvalue"])
        z = float(r["zscore"])
        f = z * z
        if p >= P_THRESHOLD or f < F_FLOOR:
            viol.append((r, p, z, f))

    print()
    if viol:
        print(f"  -- {len(viol)} instrument(s) violate the stated selection criterion --")
        for r, p, z, f in viol:
            print(f"     {r['gene_symbol']:<10} {r['snp']:<13} "
                  f"Z={z:+8.4f}  P={p:.3e}  F={f:.3f}  N={r['n_samples']}")
    else:
        print("  no threshold violations")

    clean = [r for r in src if float(r["pvalue"]) < P_THRESHOLD
             and float(r["zscore"]) ** 2 >= F_FLOOR]

    genes_before = {r["gene_symbol"] for r in src}
    genes_after = {r["gene_symbol"] for r in clean}
    dropped_genes = sorted(genes_before - genes_after)

    minF_before = min(float(r["zscore"]) ** 2 for r in src)
    minF_after = min(float(r["zscore"]) ** 2 for r in clean)

    # --- 4. downstream counts ----------------------------------------------
    mr = list(csv.DictReader(open(SRC_MR, encoding="utf-8")))
    bbj = [r for r in mr if r["outcome"] == "BBJ_Graves" and r["method"] in PRIMARY_METHODS]
    est_before = {r["gene_symbol"] for r in bbj}
    est_after = est_before - set(dropped_genes)

    thr_before = 0.05 / len(genes_before)
    thr_after = 0.05 / len(genes_after)
    hits_before = sorted(r["gene_symbol"] for r in bbj if float(r["pvalue"]) < thr_before)
    hits_after = sorted(r["gene_symbol"] for r in bbj
                        if r["gene_symbol"] in genes_after and float(r["pvalue"]) < thr_after)

    print(f"""
  -- accounting --
                             reported      corrected
    SNP-level instruments    {len(src):>8,}      {len(clean):>8,}
    MR-testable genes        {len(genes_before):>8,}      {len(genes_after):>8,}
    BBJ-estimable genes      {len(est_before):>8,}      {len(est_after):>8,}
    minimum F                {minF_before:>8.2f}      {minF_after:>8.2f}
    Bonferroni threshold     {thr_before:.4e}   {thr_after:.4e}
    both print as            {thr_before:.3e}      {thr_after:.3e}
    genes dropped            {dropped_genes}
""")

    # --- 5. conclusions must not move --------------------------------------
    print("  -- conclusions --")
    check(hits_before == hits_after,
          f"discovery hit list unchanged ({len(hits_after)} hits)")
    check(all(g not in dropped_genes for g in BACKBONE),
          "no backbone gene affected")
    check(f"{thr_before:.3e}" == f"{thr_after:.3e}",
          "Bonferroni threshold unchanged at 4 significant figures")

    # --- 6. the corrected manifest itself passes ---------------------------
    print("\n  -- corrected manifest self-check --")
    check(all(float(r["pvalue"]) < P_THRESHOLD for r in clean),
          "every instrument satisfies P < 5e-8")
    check(all(float(r["zscore"]) ** 2 >= F_FLOOR for r in clean),
          f"every instrument satisfies F = Z^2 >= {F_FLOOR:.2f}")
    check(all(len([x for x in clean if x["gene_symbol"] == g]) >= 1 for g in genes_after),
          "every retained gene keeps at least one instrument")

    if write:
        keep = {(r["gene_symbol"], r["snp"]) for r in clean}
        with open(OUT_MANIFEST, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(src[0].keys()))
            w.writeheader()
            w.writerows(clean)
        pub_clean = [r for r in pub if (r["Gene Symbol"], r["SNP"]) in keep]
        with open(OUT_TABLES4, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(pub[0].keys()))
            w.writeheader()
            w.writerows(pub_clean)
        print(f"\n  wrote {os.path.relpath(OUT_MANIFEST, ROOT)}  ({len(clean)} rows)")
        print(f"  wrote {os.path.relpath(OUT_TABLES4, ROOT)}  ({len(pub_clean)} rows)")
    else:
        print("\n  (report only — pass --write to emit the corrected files)")

    print("\nRESULT:", "ALL PASS" if not fails else f"{len(fails)} FAIL")
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
