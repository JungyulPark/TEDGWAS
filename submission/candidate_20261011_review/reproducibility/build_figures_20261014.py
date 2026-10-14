#!/usr/bin/env python3
"""Redraw Figures 1-3 and S2 from Supplementary_Dataset.xlsx (2026-10-14).

Portable replacement for the October local figure code, which read local CSV
exports. Same data, layout and number formats, with these corrections:
in-image figure titles removed (the legends carry them); gene symbols and P in
italics; Figure 1 chromosome labels no longer collide; Figure 2 marks the
screen threshold (P < 0.05/2,544) with a dagger in addition to the nominal
asterisk. Figure S1 (orbital RNA-seq) needs the IRB-protected sample data and is
not redrawn here.

    python3 submission/candidate_20261011_review/reproducibility/build_figures_20261014.py
"""
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

PKG = Path(__file__).resolve().parents[1]
X = PKG / "Supplementary_Dataset.xlsx"
OUT = PKG / "Figures"

ORDER = ["BBJ_Graves", "E4_GRAVES_STRICT", "UKB_hyperthyroid", "FinnGen_GO", "E4_GRAVES_OPHT_STRICT"]
SHORT = {"BBJ_Graves": "BBJ GD", "E4_GRAVES_STRICT": "FinnGen GD", "UKB_hyperthyroid": "UKB HT",
         "FinnGen_GO": "FinnGen GO", "E4_GRAVES_OPHT_STRICT": "FinnGen GO strict"}
COLORS = ["#28658b", "#43866b", "#777777", "#896348", "#865290"]
GENES = ["TSHR", "IGF1R", "CTLA4"]
NEWGENES = ["APOM", "NEU1", "IGLV2-18", "IGLV3-12", "TNF", "HIST1H4E", "LTA", "LTB", "C4B"]
SCREEN = 0.05 / 2544
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": "black",
                     "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                     "savefig.dpi": 250, "mathtext.default": "regular"})


def fmtp(x):
    if pd.isna(x):
        return "NA"
    if x < 0.001:
        a, b = f"{x:.2e}".split("e")
        return a + "×10" + str(int(b)).translate(SUP)
    return f"{x:.3g}"


def pp(x):
    return "<0.001" if x < 0.001 else f"{x:.3f}"


def orci(r):
    return f"{r.OR:.2f} ({r.CI_lower:.2f}–{r.CI_upper:.2f})"


def it(s):
    """Italic gene symbol for mathtext (symbols used here contain no hyphen)."""
    assert "-" not in s
    return "$\\mathit{" + s + "}$"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load():
    mr = pd.read_excel(X, sheet_name="MR_primary")
    co = pd.read_excel(X, sheet_name="Colocalization")
    ins = pd.read_excel(X, sheet_name="Instruments")
    loo = pd.read_excel(X, sheet_name="LOO_anchors_plot")
    prim = co[(co.scenario == "original_reference") & np.isclose(co.p12, 1e-5)]
    return mr, co, prim, ins, loo


def row(mr, prim, g, o):
    a = mr[(mr.gene_symbol == g) & (mr.outcome == o)].iloc[0]
    b = prim[(prim.gene == g) & (prim.outcome == o)].iloc[0]
    return pd.concat([a, b[["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]]])


def figure1(mr, ins, hits):
    fig = plt.figure(figsize=(12, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.8, 2.6], hspace=0.22)
    flow = fig.add_subplot(gs[0])
    flow.axis("off")
    flow.set_title("A  Study design", loc="left", fontweight="bold", pad=6)
    texts = ["Blood expression instruments\n4,462 druggable genes considered\n2,544 genes; 6,135 SNP–gene records",
             "BBJ GD discovery\n2,234 estimable genes; 13 discoveries\nThreshold $\\mathit{P}$ < 0.05/2,544",
             "Disease-definition comparisons\nFinnGen GD, UKB hyperthyroidism\nTwo overlapping GO definitions"]
    for i, txt in enumerate(texts):
        x = 0.01 + i * 0.334
        flow.add_patch(FancyBboxPatch((x, 0.30), 0.31, 0.62, boxstyle="round,pad=0.008", facecolor="#f1f4f5",
                                      edgecolor="#71838b", transform=flow.transAxes))
        flow.text(x + 0.015, 0.61, txt, ha="left", va="center", transform=flow.transAxes, fontsize=10, linespacing=1.7)
        if i < 2:
            flow.annotate("", xy=(x + 0.331, 0.61), xytext=(x + 0.312, 0.61), xycoords="axes fraction",
                          arrowprops={"arrowstyle": "->", "color": "black"})
    flow.text(0.01, 0.06, "Regional follow-up: " + it("TSHR") + " and " + it("IGF1R") + " at the center; "
              + it("CTLA4") + " comparator; 24 selected genes in total.", transform=flow.transAxes, fontsize=10.5)

    ax = fig.add_subplot(gs[1])
    lead = ins.sort_values("pvalue").drop_duplicates("gene_symbol")
    d = mr[(mr.outcome == "BBJ_Graves") & mr.estimable].merge(lead[["gene_symbol", "chr", "pos_hg19"]], on="gene_symbol")
    cmax = d.groupby("chr").pos_hg19.max()
    offset, total, ticks = {}, 0.0, []
    for c in range(1, 23):
        offset[c] = total
        size = float(cmax.get(c, 1e8))
        ticks.append(total + size / 2)
        total += size + 1e7
    d["x"] = [offset[int(c)] + p for c, p in zip(d.chr, d.pos_hg19)]
    d["y"] = -np.log10(d.pvalue)
    assert len(d) == 2234, len(d)
    for c, g in d.groupby("chr"):
        ax.scatter(g.x, g.y, s=8, color="#bac7ce" if c % 2 else "#92a6b2", alpha=0.7, linewidths=0)
    thr = -np.log10(SCREEN)
    ax.axhline(thr, color="#555555", ls="--", lw=1.2)
    ax.text(offset[7] + 3e7, thr + 0.45, "Bonferroni threshold: $\\mathit{P}$ < 1.965×10⁻⁵", fontsize=9)
    offs = {"HLA-A": (0, 12), "HLA-DQA2": (-30, 14), "C4A": (-26, 23), "TUBB": (-42, 18), "PSMB8": (-38, -24),
            "IFNGR1": (26, 25), "HSD3B7": (0, 56), "VKORC1": (-24, 36), "PRSS36": (24, 14), "TNFSF14": (10, 22),
            "MAPKAPK5": (-12, 28), "CTLA4": (0, 13), "TSHR": (0, 13), "IGF1R": (0, 20)}
    special = {"TSHR": "#146891", "IGF1R": "#934912", "CTLA4": "#73558a"}
    for g in hits + ["IGF1R"]:
        r = d[d.gene_symbol == g].iloc[0]
        ax.scatter([r.x], [r.y], s=60 if g == "IGF1R" else 30, color=special.get(g, "#344f5d"), zorder=5)
        dx, dy = offs[g]
        ax.annotate(g, (r.x, r.y), xytext=(dx, dy), textcoords="offset points", ha="center", fontsize=9,
                    fontstyle="italic", color="black", arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "black"})
    labels = [str(c) if c <= 18 or c in (20, 22) else "" for c in range(1, 23)]
    ax.set_xticks(ticks, labels)
    ax.set_ylim(-0.3, 25.5)
    ax.set_xlim(-1e7, total)
    ax.set_ylabel("Association strength ($-\\log_{10}\\mathit{P}$)")
    ax.set_xlabel("Chromosome of strongest selected expression instrument (GRCh37)")
    ax.set_title("B  BBJ gene-level MR associations", loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=0.18)
    save(fig, "Figure1")
    return len(d)


def figure2(mr, prim):
    fig = plt.figure(figsize=(14, 8.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[3, 1.7, 2.7], wspace=0.07)
    ax = fig.add_subplot(gs[0])
    tx = fig.add_subplot(gs[1], sharey=ax)
    bar = fig.add_subplot(gs[2], sharey=ax)
    rows = [(g, o) for g in GENES for o in ORDER]
    shown = []
    for i, (g, o) in enumerate(rows):
        r = row(mr, prim, g, o)
        c = COLORS[ORDER.index(o)]
        ax.errorbar(r.OR, i, xerr=[[r.OR - r.CI_lower], [r.CI_upper - r.OR]], fmt="o", color=c, ecolor=c,
                    capsize=3, markersize=5)
        ax.text(-0.53, i, g, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=10,
                fontstyle="italic")
        ax.text(-0.33, i, SHORT[o], transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=10)
        mark = (" *" if r.pvalue < 0.05 else "") + ("†" if r.pvalue < SCREEN else "")
        tx.text(0.02, i, orci(r), ha="left", va="center", fontsize=9)
        tx.text(0.73, i, fmtp(r.pvalue) + mark, ha="left", va="center", fontsize=9,
                fontweight="bold" if r.pvalue < 0.05 else "normal")
        vals = [r["PP.H0"] + r["PP.H1"], r["PP.H2"], r["PP.H3"], r["PP.H4"]]
        left = 0
        for v, col in zip(vals, ["#eeeeee", "#c3cbd0", "#dcab80", "#72a9bf"]):
            bar.barh(i, v, left=left, height=0.6, color=col, edgecolor="white", lw=0.5)
            left += v
        bar.text(1.025, i, pp(r["PP.H4"]), va="center", fontsize=10, fontweight="bold" if r["PP.H4"] >= 0.8 else "normal")
        shown.append((g, o, orci(r), fmtp(r.pvalue), pp(r["PP.H2"]), pp(r["PP.H3"]), pp(r["PP.H4"])))
    for a in (ax, tx, bar):
        a.set_ylim(14.8, -0.8)
        for y in (4.5, 9.5):
            a.axhline(y, color="#cccccc", lw=0.7)
    ax.set_xscale("log")
    ax.set_xlim(0.025, 4)
    ax.set_xticks([0.03, 0.1, 0.3, 1, 3], ["0.03", "0.1", "0.3", "1", "3"])
    ax.axvline(1, color="black", ls="--", lw=0.8)
    ax.set_yticks(range(15), [""] * 15)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("OR per approximately 1 SD higher expression")
    ax.set_title("A  MR estimate", loc="left", fontweight="bold")
    tx.axis("off")
    tx.set_xlim(0, 1.32)
    tx.set_title("OR (95% CI)             $\\mathit{P}$ value", loc="left", fontsize=10)
    bar.tick_params(axis="y", left=False, labelleft=False)
    bar.set_xlim(0, 1.18)
    bar.set_xticks([0, 0.5, 1])
    bar.axvline(0.8, color="black", ls=":", lw=0.8)
    bar.set_xlabel("Posterior probability")
    bar.set_title("B  Colocalization                         H4", loc="left", fontweight="bold")
    fig.legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(
        ["#eeeeee", "#c3cbd0", "#dcab80", "#72a9bf"],
        ["H0/H1", "H2: expression only", "H3: distinct variants", "H4: shared variant"])],
        loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.57, 0.015))
    fig.subplots_adjust(left=0.24, right=0.95, top=0.94, bottom=0.12)
    save(fig, "Figure2")
    return shown


def figure3(co, prim, hits):
    fig = plt.figure(figsize=(12, 11.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1], wspace=0.34)
    hm = fig.add_subplot(gs[0])
    pr = fig.add_subplot(gs[1])
    order = GENES + [g for g in hits if g not in GENES] + NEWGENES + ["HLA-G"]
    assert len(order) == 24
    m = np.array([[prim[(prim.gene == g) & (prim.outcome == o)]["PP.H4"].iloc[0] for o in ORDER] for g in order])
    hm.imshow(m, cmap="Blues", vmin=0, vmax=1, aspect="auto", alpha=0.48)
    for i in range(24):
        for j in range(5):
            hm.text(j, i, pp(m[i, j]), ha="center", va="center", fontsize=9, color="black")
            if m[i, j] >= 0.8:
                hm.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="black", linewidth=1.2))
    hm.set_xticks(range(5), ["BBJ\nGD", "FinnGen\nGD", "UKB\nHT", "FinnGen\nGO", "FinnGen\nGO strict"], fontsize=9)
    hm.set_yticks(range(24), order, fontstyle="italic")
    hm.tick_params(length=0)
    hm.set_title("A  Shared-variant probability (PP.H4)", loc="left", fontweight="bold", pad=14)
    for y in (2.5, 13.5, 22.5):
        hm.axhline(y, color="black", lw=1)
    for o, c in zip(ORDER, COLORS):
        z = co[(co.gene == "TSHR") & (co.outcome == o) & (co.scenario == "original_reference")].sort_values("p12", ascending=False)
        pr.plot(range(3), z["PP.H4"], marker="o", lw=1.7, color=c, label=SHORT[o])
    pr.axhline(0.8, color="black", ls="--", lw=0.8)
    pr.set_ylim(0, 1.05)
    pr.set_xticks(range(3), ["10⁻⁵", "5×10⁻⁶", "10⁻⁶"])
    pr.set_xlabel("Shared-variant prior (p12)")
    pr.set_ylabel("PP.H4")
    pr.set_title("B  $\\mathbfit{TSHR}$ prior sensitivity", loc="left", fontweight="bold", pad=14)
    pr.legend(frameon=False, loc="center left", bbox_to_anchor=(0, 0.4), fontsize=10)
    pr.grid(alpha=0.15)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.96, bottom=0.07)
    save(fig, "Figure3")
    return m, order


def figureS2(loo):
    fig, axes = plt.subplots(1, 2, figsize=(14, 9.4), gridspec_kw={"width_ratios": [1, 1]})
    fig.subplots_adjust(left=0.17, right=0.86, wspace=0.9, top=0.93, bottom=0.08)
    assert len(loo) == 25
    for ax, g in zip(axes, ["IGF1R", "CTLA4"]):
        z = loo[loo.gene == g].copy()
        z["sort"] = z.outcome.map({o: i for i, o in enumerate(ORDER)})
        z = z.sort_values(["sort", "removed_snp"])
        labels = []
        for i, (_, r) in enumerate(z.iterrows()):
            ax.errorbar(r.OR, i, xerr=[[r.OR - r.CI_lower], [r.CI_upper - r.OR]], fmt="o",
                        color=COLORS[ORDER.index(r.outcome)], capsize=3, markersize=4)
            labels.append(SHORT[r.outcome] + "\n− " + r.removed_snp)
            ax.text(1.04, i, fmtp(r.pvalue), transform=ax.get_yaxis_transform(), va="center", fontsize=9,
                    fontweight="bold" if r.pvalue < 0.05 else "normal")
        ax.set_yticks(range(len(z)), labels, fontsize=8.5)
        ax.set_ylim(len(z) - 0.25, -0.75)
        ax.set_xscale("log")
        ticks = [0.5, 1, 2, 4] if g == "IGF1R" else [0.1, 1, 10]
        ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.axvline(1, color="black", ls="--", lw=0.8)
        ax.set_xlabel("OR after omission (95% CI)")
        ax.set_title("$\\mathbfit{" + g + "}$", loc="left", fontsize=11)
        ax.text(1.04, 1.005, "$\\mathit{P}$ value", transform=ax.transAxes, fontsize=10, va="bottom")
        ax.grid(axis="x", alpha=0.15)
    save(fig, "FigureS2")


def main():
    mr, co, prim, ins, loo = load()
    hits = mr[(mr.outcome == "BBJ_Graves") & mr.P_lt_0_05_over_2544].sort_values("pvalue").gene_symbol.tolist()
    assert len(hits) == 13
    n = figure1(mr, ins, hits)
    shown = figure2(mr, prim)
    figure3(co, prim, hits)
    figureS2(loo)
    pd.DataFrame(shown, columns=["gene", "outcome", "OR_CI", "P", "PP.H2", "PP.H3", "PP.H4"]).to_csv(
        PKG / "Verification" / "Review_20261014" / "figure2_displayed_values.csv", index=False)
    print(f"Figure 1: {n} genes; Figures 1-3 and S2 written to {OUT.relative_to(PKG.parents[1])}")


if __name__ == "__main__":
    main()
