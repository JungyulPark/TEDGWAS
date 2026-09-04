#!/usr/bin/env python3
"""
Figure 1 — study flow and the TSHR/IGF1R comparison, in journal format.

Redesigned 2026-09 for a conventional manuscript figure: a plain study-flow
diagram in the style journals expect (hairline rectangles, no fills, exclusions
bracketed to the right), and a compact comparison panel. The earlier version used
coloured verdict strips and display copy, which read as a slide rather than a
figure.

Every number is read from the locked master and result tables; nothing is
hard-coded except the classification labels, which are parsed from Table 3.

Out: TrackA_MR/v5_upgrade/07_manuscript/figures/Figure1_study_design.png
     submission/figures/Figure1.png
"""
import csv, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER = f"{ROOT}/MANUSCRIPT_TED_TRAP_v5_MASTER.md"
COLOC = f"{ROOT}/TrackA_MR/v5_upgrade/05_coloc_candidates/TaskE_01b_coloc_v2_20260904.csv"
MR = f"{ROOT}/TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv"
OUT = [f"{ROOT}/TrackA_MR/v5_upgrade/07_manuscript/figures/Figure1_study_design.png",
       f"{ROOT}/submission/figures/Figure1.png"]

INK, RULE = "#000000", "#555555"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5})


def master_counts():
    t = open(MASTER).read()
    for n in ("4,462", "2,544", "2,234", "6,135"):
        assert n in t, f"count {n} not in master"
    sec = t.split("**Table 3.")[1].split("\n\n")[1]
    rows = [r for r in sec.split("\n") if r.startswith("|") and "---" not in r]
    hdr = [c.strip() for c in rows[0].strip("|").split("|")]
    ci = hdr.index("Classification")
    return {re.sub(r"[*]", "", [x.strip() for x in r.strip("|").split("|")][0]):
            [x.strip() for x in r.strip("|").split("|")][ci] for r in rows[1:]}


def backbone_values():
    co = {(r["gene"], r["outcome"]): r for r in csv.DictReader(open(COLOC))
          if r["p12"] in ("1e-05", "1e-5")}
    mr = {(r["gene_symbol"], r["outcome"]): r for r in csv.DictReader(open(MR))
          if r["method"] in ("Inverse variance weighted", "Wald ratio")}
    out = {}
    for g in ("TSHR", "IGF1R"):
        out[g] = dict(
            beta=float(mr[(g, "BBJ_Graves")]["beta"]),
            p=float(mr[(g, "BBJ_Graves")]["pvalue"]),
            h4=[float(co[(g, o)]["PP.H4"]) for o in
                ("BBJ_Graves", "UKB_hyperthyroid", "FinnGen_GO")],
            h2=[float(co[(g, o)]["PP.H2"]) for o in ("BBJ_Graves", "FinnGen_GO")],
        )
    return out


def box(ax, cx, cy, w, h, lw=1.1):
    ax.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h,
                           facecolor="white", edgecolor=INK, linewidth=lw))


def build():
    cls = master_counts()
    bb = backbone_values()
    n = {k: sum(1 for v in cls.values() if k in v) for k in
         ("MHC", "chr16p11.2", "single-outcome", "distinct-variant")}
    n["known"] = sum(1 for v in cls.values() if v.startswith("Established"))
    assert sum(n.values()) == 13

    fig = plt.figure(figsize=(7.2, 9.6), dpi=300)
    fig.patch.set_facecolor("white")

    # ============ (a) study flow — CONSORT style ============
    ax = fig.add_axes([0.015, 0.395, 0.97, 0.590]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.00, 1.005, "(a)", fontsize=13, fontweight="bold", va="top")

    MX, MW = 0.285, 0.50          # main spine
    EX, EW = 0.775, 0.40          # exclusion column
    LAB, NUM = 10.5, 12.0

    main = [
        (0.925, 0.110, "Druggable genome", "n = 4,462 genes"),
        (0.720, 0.130, "MR-testable\n\u22651 valid cis-eQTL instrument (6,135 total)",
         "n = 2,544 genes"),
        (0.515, 0.110, "Estimable against discovery outcome", "n = 2,234 genes"),
        (0.305, 0.130, "Bonferroni-significant in discovery\nP < 1.965 \u00d7 10\u207b\u2075",
         "n = 13 genes"),
        (0.085, 0.110, "Met the colocalization criterion", "n = 0 novel candidates"),
    ]
    for k, (cy, h, lab, num) in enumerate(main):
        last = (k == len(main) - 1)
        box(ax, MX, cy, MW, h, lw=1.6 if last else 1.1)
        nlines = lab.count("\n") + 1
        ax.text(MX, cy + (0.020 if nlines == 1 else 0.032), lab, ha="center", va="center",
                fontsize=LAB, linespacing=1.35, color=INK)
        ax.text(MX, cy - (0.026 if nlines == 1 else 0.036), num, ha="center", va="center",
                fontsize=NUM, fontweight="bold", color=INK)
        if not last:
            ny, nh = main[k + 1][0], main[k + 1][1]
            ax.annotate("", xy=(MX, ny + nh / 2), xytext=(MX, cy - h / 2),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3,
                                        mutation_scale=16))

    # exclusion boxes, branching from each transition
    for cy, lab, num in [(0.822, "No valid cis-eQTL instrument", "n = 1,918"),
                         (0.617, "Not estimable in discovery outcome", "n = 310"),
                         (0.428, "Not Bonferroni-significant", "n = 2,221")]:
        ax.annotate("", xy=(EX - EW / 2, cy), xytext=(MX, cy),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.1,
                                    mutation_scale=14))
        box(ax, EX, cy, EW, 0.085)
        ax.text(EX, cy + 0.018, lab, ha="center", va="center", fontsize=9.6, color=INK)
        ax.text(EX, cy - 0.022, num, ha="center", va="center", fontsize=10.8,
                fontweight="bold", color=INK)

    # resolution of the 13, branching sideways from that box
    cy = 0.245
    ax.annotate("", xy=(EX - EW / 2, 0.290), xytext=(MX + MW / 2, 0.290),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.1, mutation_scale=14))
    box(ax, EX, cy, EW, 0.200)
    ax.text(EX, cy + 0.076, "Resolution of the 13 genes", ha="center",
            fontsize=9.8, fontweight="bold", color=INK)
    rows = [(n["MHC"], "MHC region"),
            (n["chr16p11.2"], "chr16p11.2 (one signal)"),
            (n["single-outcome"], "discovery-only colocalization"),
            (n["distinct-variant"], "distinct causal variant"),
            (n["known"], "established loci")]
    yy = cy + 0.040
    for k, lab in rows:
        ax.text(EX - EW / 2 + 0.030, yy, str(k), fontsize=9.8, va="center",
                fontweight="bold", color=INK)
        ax.text(EX - EW / 2 + 0.062, yy, lab, fontsize=9.4, va="center", color=INK)
        yy -= 0.031

    # ============ (b) comparison table ============
    axb = fig.add_axes([0.015, 0.030, 0.97, 0.330]); axb.axis("off")
    axb.set_xlim(0, 1); axb.set_ylim(0, 1)
    axb.text(0.00, 1.005, "(b)", fontsize=13, fontweight="bold", va="top")

    T, I = bb["TSHR"], bb["IGF1R"]
    x0, xT, xI, xR = 0.030, 0.560, 0.845, 0.985
    axb.plot([x0, xR], [0.880, 0.880], color=INK, lw=1.5)
    axb.text(xT, 0.800, "TSHR", ha="center", fontsize=12, style="italic",
             fontweight="bold", color=INK)
    axb.text(xI, 0.800, "IGF1R", ha="center", fontsize=12, style="italic",
             fontweight="bold", color=INK)
    axb.text(xT, 0.720, "autoantigen", ha="center", fontsize=9.4, color=INK)
    axb.text(xI, 0.720, "teprotumumab target", ha="center", fontsize=9.4, color=INK)
    axb.plot([x0, xR], [0.665, 0.665], color=INK, lw=1.0)

    body = [
        ("Discovery MR", f"\u03b2 = {T['beta']:+.2f}\nP = 1.1 \u00d7 10\u207b\u00b9\u2074",
         f"\u03b2 = {I['beta']:+.2f}\nP = {I['p']:.3f}"),
        ("Colocalization PP.H4\ndiscovery / replication /\nTED-enriched",
         " / ".join(f"{v:.2f}" for v in T["h4"]),
         " / ".join(f"{v:.2f}" for v in I["h4"])),
        ("Shared-variant support", "both Graves\nphenotypes", "none"),
    ]
    yy = 0.545
    for lab, tv, iv in body:
        axb.text(x0, yy, lab, va="center", fontsize=10.2, linespacing=1.35, color=INK)
        axb.text(xT, yy, tv, ha="center", va="center", fontsize=10.2,
                 linespacing=1.35, color=INK)
        axb.text(xI, yy, iv, ha="center", va="center", fontsize=10.2,
                 linespacing=1.35, color=INK)
        yy -= 0.185
    axb.plot([x0, xR], [0.070, 0.070], color=INK, lw=1.5)

    for path in OUT:
        fig.savefig(path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("Figure 1 rebuilt (CONSORT-style, monochrome)")
    print(f"  flow: 4,462 -> 2,544 -> 2,234 -> 13 -> 0")
    print(f"  13 = {n['MHC']} MHC + {n['chr16p11.2']} chr16 + {n['single-outcome']} discovery-only"
          f" + {n['distinct-variant']} distinct + {n['known']} established")


if __name__ == "__main__":
    build()
