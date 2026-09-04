#!/usr/bin/env python3
"""Rebuild Figures 3, S1 and S2 on the three-outcome colocalization (v2).

What changed and why
--------------------
The v1 figures were built from a two-outcome colocalization run (BBJ and FinnGen)
that used p-value-based ABF with a single European allele-frequency vector, and
they carried an orbital-tissue panel that was presented as an evidence layer.
Both are gone:

  * colocalization was re-run for all three outcomes with beta/varbeta and
    ancestry-matched allele frequencies (TaskE_01b_coloc_v2), so Figure 3A and
    Figure S2 now show three bars per gene rather than two;
  * the orbital tissue data cannot support inference with a single control, so
    it leaves Figure 3 entirely and becomes Figure S1, drawn from TPM (not the
    DESeq2 shrunken estimate, since DESeq2 is no longer used), without p-values
    and labelled as a descriptive observation.

Every plotted number is read from the locked result tables. Nothing is hard-coded.

Out: TrackA_MR/v5_upgrade/07_manuscript/figures/{Figure3,FigureS1,FigureS2}*.png
     submission/figures/{Figure3,FigureS1,FigureS2}.png
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLOC = f"{ROOT}/TrackA_MR/v5_upgrade/05_coloc_candidates/TaskE_01b_coloc_v2_20260904.csv"
MR = f"{ROOT}/TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv"
TISSUE = f"{ROOT}/TrackA_MR/v5_upgrade/06_tissue_integration/TaskF_01_backbone_tissue_inhouse_v1.csv"

FIGDIR = f"{ROOT}/TrackA_MR/v5_upgrade/07_manuscript/figures"
SUBDIR = f"{ROOT}/submission/figures"

INK, MUTED, RULE = "#1a1d23", "#5c6673", "#c9d0d9"
ANCHOR, EFFECTOR = "#a4262c", "#2e75b6"
H2C, H3C, H4C = "#b9c0c9", "#e0a458", "#3f6f3f"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})

BACKBONE = ["TSHR", "IGF1R", "CTLA4"]
CANDIDATES = ["TNFSF14", "IFNGR1", "MAPKAPK5", "HSD3B7", "VKORC1", "PRSS36"]
OUTCOMES = [("BBJ_Graves", "BBJ\nGraves"),
            ("UKB_hyperthyroid", "UKB\nhyperthyroid"),
            ("FinnGen_GO", "FinnGen\nGO")]
PRIMARY = ("Inverse variance weighted", "Wald ratio")


def load_coloc():
    """p12 = 1e-5 rows only — the default prior reported in the tables."""
    out = {}
    with open(COLOC) as fh:
        for r in csv.DictReader(fh):
            if r["p12"] in ("1e-05", "1e-5"):
                out[(r["gene"], r["outcome"])] = r
    return out


def load_mr():
    out = {}
    with open(MR) as fh:
        for r in csv.DictReader(fh):
            if r["method"] in PRIMARY:
                out[(r["gene_symbol"], r["outcome"])] = r
    return out


def load_tissue():
    with open(TISSUE) as fh:
        return {r["gene"]: r for r in csv.DictReader(fh)}


def save(fig, stem, also_submission=True):
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(SUBDIR, exist_ok=True)
    paths = [f"{FIGDIR}/{stem}.png"]
    if also_submission:
        paths.append(f"{SUBDIR}/{stem}.png")
    for p in paths:
        fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return paths


def stacked_coloc(ax, genes, co):
    """H2 / H3 / H4 stacked bars, three outcomes per gene."""
    labels, h2, h3, h4, marks = [], [], [], [], []
    for g in genes:
        for key, lab in OUTCOMES:
            r = co.get((g, key))
            if r is None:
                continue
            labels.append((g, lab))
            h2.append(float(r["PP.H2"]))
            h3.append(float(r["PP.H3"]))
            h4.append(float(r["PP.H4"]))
            marks.append(r["top_snp"])
    x = range(len(labels))
    ax.bar(x, h2, color=H2C, label="PP.H2 (cis-eQTL only)")
    ax.bar(x, h3, bottom=h2, color=H3C, label="PP.H3 (distinct variants)")
    ax.bar(x, h4, bottom=[a + b for a, b in zip(h2, h3)], color=H4C,
           label="PP.H4 (shared variant)")
    for xi, v in enumerate(h4):
        ax.text(xi, 1.02, f"{v:.2f}", ha="center", fontsize=7.2, color=INK)
    ax.axhline(0.8, color=ANCHOR, ls="--", lw=1.1)
    ax.text(len(labels) - 0.4, 0.815, "PP.H4 = 0.8", fontsize=7.5,
            color=ANCHOR, ha="right")
    ax.set_xticks(list(x))
    ax.set_xticklabels([l for _, l in labels], fontsize=7.6)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("posterior probability", fontsize=9)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # gene bands under the outcome labels
    step = len(OUTCOMES)
    for i, g in enumerate(genes):
        c = ANCHOR if g == "TSHR" else EFFECTOR if g == "IGF1R" else INK
        ax.text(i * step + (step - 1) / 2, -0.20, g, ha="center",
                fontsize=11, fontweight="bold", style="italic",
                color=c, transform=ax.get_xaxis_transform())
        if i:
            ax.axvline(i * step - 0.5, color=RULE, lw=0.8)
    return labels, marks


def figure3(co, mr):
    fig = plt.figure(figsize=(11.0, 9.4), dpi=300)
    fig.patch.set_facecolor("white")

    # ---------------- Panel A ----------------
    fig.text(0.055, 0.972, "A", fontsize=15, fontweight="bold", va="top")
    fig.text(0.085, 0.972, "Colocalization across all three outcomes",
             fontsize=12.5, fontweight="bold", va="top")
    fig.text(0.085, 0.949, "TSHR colocalizes in both Graves disease phenotypes but not in the broader "
                           "hyperthyroidism outcome, where distinct variants are favoured.",
             fontsize=8.8, color=MUTED, va="top")

    axA = fig.add_axes([0.085, 0.605, 0.875, 0.275])
    stacked_coloc(axA, BACKBONE, co)
    axA.legend(loc="lower center", fontsize=8.0, frameon=False, ncol=3,
               bbox_to_anchor=(0.5, 1.02))

    # ---------------- Panel B ----------------
    fig.text(0.055, 0.470, "B", fontsize=15, fontweight="bold", va="top")
    fig.text(0.085, 0.470, "MR effect estimates with 95% confidence intervals",
             fontsize=12.5, fontweight="bold", va="top")
    fig.text(0.085, 0.447, "TSHR and CTLA4 protective in every outcome; IGF1R risk-direction, "
                           "nominally significant in BBJ and UKB only.",
             fontsize=8.8, color=MUTED, va="top")

    axB = fig.add_axes([0.235, 0.075, 0.615, 0.335])
    ys, lab, col = [], [], []
    y = 0
    for g in BACKBONE:
        gy = []
        for key, l in OUTCOMES:
            r = mr.get((g, key))
            if r is None:
                continue
            b, se = float(r["beta"]), float(r["se"])
            c = ANCHOR if g == "TSHR" else EFFECTOR if g == "IGF1R" else INK
            axB.errorbar(b, y, xerr=1.96 * se, fmt="o", ms=5.5, color=c,
                         ecolor=c, elinewidth=1.5, capsize=3)
            pv = float(r["pvalue"])
            ptxt = (f"P = {pv:.3f}" if pv >= 1e-3
                    else f"P = {pv:.1e}".replace("e-0", " \u00d7 10\u207b").replace("e-", " \u00d7 10\u207b"))
            axB.text(1.015, y, ptxt, transform=axB.get_yaxis_transform(),
                     va="center", ha="left", fontsize=7.8, color=c)
            lab.append(l.replace("\n", " "))
            col.append(c)
            ys.append(y); gy.append(y)
            y -= 1
        c = ANCHOR if g == "TSHR" else EFFECTOR if g == "IGF1R" else INK
        axB.text(-0.355, sum(gy) / len(gy), g, transform=axB.get_yaxis_transform(),
                 ha="left", va="center", fontsize=11, fontweight="bold",
                 style="italic", color=c)
        y -= 0.8
    axB.axvline(0, color=MUTED, lw=1)
    axB.set_yticks(ys)
    axB.set_yticklabels(lab, fontsize=8.4)
    for t, c in zip(axB.get_yticklabels(), col):
        t.set_color(c)
    axB.set_xlabel("MR effect on disease (β, log-odds per unit genetically proxied expression)",
                   fontsize=9.2)
    for sp in ("top", "right", "left"):
        axB.spines[sp].set_visible(False)
    axB.tick_params(axis="y", length=0)
    return save(fig, "Figure3")


def figure_s1(ti):
    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.5), dpi=300)
    fig.patch.set_facecolor("white")
    for ax, g in zip(axes, BACKBONE):
        r = ti[g]
        ctrl = float(r["ctrl_tpm"])
        ted = [float(x) for x in r["ted_per_sample"].split(";")]
        ax.scatter([0] * len(ted), ted, s=46, color=EFFECTOR, zorder=3,
                   label="TED (n = 4)")
        ax.scatter([1], [ctrl], s=64, marker="s", color=INK, zorder=3,
                   label="control (n = 1)")
        ax.set_xlim(-0.6, 1.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["TED", "control"], fontsize=9)
        ax.set_ylim(0, max(max(ted), ctrl) * 1.35)
        ax.set_title(g, fontsize=11.5, fontweight="bold", style="italic",
                     color=ANCHOR if g == "TSHR" else EFFECTOR if g == "IGF1R" else INK)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if ax is axes[0]:
            ax.set_ylabel("transcript abundance (TPM)", fontsize=9)
            ax.legend(fontsize=7.4, frameon=False, loc="upper left")
    return save(fig, "FigureS1")


def figure_s2(co):
    fig = plt.figure(figsize=(12.0, 5.4), dpi=300)
    fig.patch.set_facecolor("white")
    fig.text(0.055, 0.965, "Candidate colocalization across all three outcomes",
             fontsize=12.5, fontweight="bold", va="top")
    fig.text(0.055, 0.930, "TNFSF14 and IFNGR1 colocalize in the discovery outcome only; neither "
                           "reproduces in UKB, the largest and ancestry-matched outcome.",
             fontsize=8.8, color=MUTED, va="top")
    ax = fig.add_axes([0.062, 0.215, 0.925, 0.615])
    stacked_coloc(ax, CANDIDATES, co)
    ax.legend(loc="lower center", fontsize=7.8, frameon=False, ncol=3,
              bbox_to_anchor=(0.5, 1.02))
    return save(fig, "FigureS2")


def main():
    co, mr, ti = load_coloc(), load_mr(), load_tissue()
    for g in BACKBONE + CANDIDATES:
        for key, _ in OUTCOMES:
            assert (g, key) in co, f"missing coloc: {g} x {key}"
    print("Rebuilding figures on the three-outcome colocalization (v2)")
    for fn, args in ((figure3, (co, mr)), (figure_s1, (ti,)), (figure_s2, (co,))):
        for p in fn(*args):
            print("  wrote", os.path.relpath(p, ROOT))
    print(f"  TSHR PP.H4  BBJ {float(co[('TSHR','BBJ_Graves')]['PP.H4']):.3f}"
          f" / UKB {float(co[('TSHR','UKB_hyperthyroid')]['PP.H4']):.3f}"
          f" / FinnGen {float(co[('TSHR','FinnGen_GO')]['PP.H4']):.3f}")
    print(f"  IGF1R PP.H2 BBJ {float(co[('IGF1R','BBJ_Graves')]['PP.H2']):.3f}"
          f" / FinnGen {float(co[('IGF1R','FinnGen_GO')]['PP.H2']):.3f}")


if __name__ == "__main__":
    main()
