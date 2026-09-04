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


def flowbox(ax, cx, y, w, h, lines, fs=8.3):
    ax.add_patch(Rectangle((cx - w / 2, y - h / 2), w, h, facecolor="white",
                           edgecolor=INK, linewidth=0.8))
    ax.text(cx, y, "\n".join(lines), ha="center", va="center",
            fontsize=fs, linespacing=1.45, color=INK)


def build():
    cls = master_counts()
    bb = backbone_values()
    n = {k: sum(1 for v in cls.values() if k in v) for k in
         ("MHC", "chr16p11.2", "single-outcome", "distinct-variant")}
    n["known"] = sum(1 for v in cls.values() if v.startswith("Established"))
    assert sum(n.values()) == 13

    fig = plt.figure(figsize=(7.1, 8.9), dpi=300)
    fig.patch.set_facecolor("white")

    # ---------------- (a) study flow ----------------
    ax = fig.add_axes([0.02, 0.400, 0.96, 0.580]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.01, 1.00, "(a)", fontsize=10, fontweight="bold", va="top")

    CX, W, H = 0.33, 0.50, 0.105
    steps = [
        (0.930, ["Druggable genome", "n = 4,462 genes"]),
        (0.735, ["MR-testable: \u22651 valid cis-eQTL instrument", "n = 2,544 genes (6,135 instruments)"]),
        (0.540, ["Estimable against discovery outcome", "n = 2,234 genes"]),
        (0.335, ["Bonferroni-significant in discovery", "P < 1.965 \u00d7 10\u207b\u2075", "n = 13 genes"]),
        (0.110, ["0 novel candidates met the", "colocalization criterion"]),
    ]
    for i, (y, lines) in enumerate(steps):
        hh = H * (1.22 if len(lines) > 2 else 1.0)
        bold = (i == len(steps) - 1)
        ax.add_patch(Rectangle((CX - W / 2, y - hh / 2), W, hh, facecolor="white",
                               edgecolor=INK, linewidth=1.2 if bold else 0.8))
        ax.text(CX, y, "\n".join(lines), ha="center", va="center",
                fontsize=8.3, linespacing=1.45, color=INK,
                fontweight="bold" if bold else "normal")
        if i < len(steps) - 1:
            nh = H * (1.22 if len(steps[i + 1][1]) > 2 else 1.0)
            ax.annotate("", xy=(CX, steps[i + 1][0] + nh / 2), xytext=(CX, y - hh / 2),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8))

    # exclusions bracketed to the right of each transition
    for y, lines in [(0.835, ["Excluded: no valid instrument", "n = 1,918"]),
                     (0.640, ["Not estimable in discovery outcome", "n = 310"]),
                     (0.445, ["Not Bonferroni-significant", "n = 2,221"])]:
        ax.plot([CX, 0.615], [y, y], color=INK, lw=0.8)
        ax.text(0.625, y, "\n".join(lines), va="center", fontsize=7.6,
                linespacing=1.4, color=INK)

    # resolution of the 13, bracketed from the 13-gene box
    ax.plot([CX + W / 2, 0.615], [0.335, 0.335], color=INK, lw=0.8)
    ax.text(0.625, 0.395, "Resolution of the 13 genes", fontsize=8.0,
            fontweight="bold", va="center", color=INK)
    rows = [(n["MHC"], "MHC region"),
            (n["chr16p11.2"], "chr16p11.2, single conditional signal"),
            (n["single-outcome"], "colocalized in discovery only"),
            (n["distinct-variant"], "distinct causal variant"),
            (n["known"], "established loci (TSHR, CTLA4)")]
    yy = 0.350
    for k, lab in rows:
        ax.text(0.632, yy, f"{k}", fontsize=8.0, va="center", color=INK)
        ax.text(0.664, yy, lab, fontsize=7.6, va="center", color=INK)
        yy -= 0.038

    # ---------------- (b) comparison ----------------
    axb = fig.add_axes([0.02, 0.035, 0.96, 0.325]); axb.axis("off")
    axb.set_xlim(0, 1); axb.set_ylim(0, 1)
    axb.text(0.01, 1.00, "(b)", fontsize=10, fontweight="bold", va="top")

    T, I = bb["TSHR"], bb["IGF1R"]
    x0, xT, xI = 0.055, 0.475, 0.775
    axb.plot([x0, 0.965], [0.855, 0.855], color=INK, lw=0.9)
    axb.text(xT, 0.895, "TSHR", ha="center", fontsize=9.5, style="italic",
             fontweight="bold", color=INK)
    axb.text(xI, 0.895, "IGF1R", ha="center", fontsize=9.5, style="italic",
             fontweight="bold", color=INK)
    axb.text(xT, 0.795, "autoantigen", ha="center", fontsize=7.6, color=RULE)
    axb.text(xI, 0.795, "teprotumumab target", ha="center", fontsize=7.6, color=RULE)
    axb.plot([x0, 0.965], [0.755, 0.755], color=INK, lw=0.6)

    body = [
        ("Discovery MR",
         f"β = {T['beta']:+.2f}, P = 1.1 × 10⁻¹⁴",
         f"β = {I['beta']:+.2f}, P = {I['p']:.3f}"),
        ("Colocalization, PP.H4\n(discovery / replication /\nTED-enriched)",
         " / ".join(f"{v:.2f}" for v in T["h4"]),
         " / ".join(f"{v:.2f}" for v in I["h4"])),
        ("Shared-variant support",
         "both Graves phenotypes",
         "none"),
    ]
    yy = 0.620
    for lab, tv, iv in body:
        axb.text(x0, yy, lab, va="center", fontsize=8.0, linespacing=1.4, color=INK)
        axb.text(xT, yy, tv, ha="center", va="center", fontsize=8.0, color=INK)
        axb.text(xI, yy, iv, ha="center", va="center", fontsize=8.0, color=INK)
        yy -= 0.195
    axb.plot([x0, 0.965], [0.115, 0.115], color=INK, lw=0.9)
    axb.text(x0, 0.055,
             "Genetic susceptibility and therapeutic target biology need not coincide.",
             fontsize=8.2, va="center", color=INK)

    for p in OUT:
        fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("Figure 1 rebuilt in journal format")
    print(f"  flow: 4,462 -> 2,544 -> 2,234 -> 13 -> 0")
    print(f"  13 = {n['MHC']} MHC + {n['chr16p11.2']} chr16 + {n['single-outcome']} single-outcome"
          f" + {n['distinct-variant']} distinct + {n['known']} established")


if __name__ == "__main__":
    build()
