#!/usr/bin/env python3
"""Draw Figures 1-3 and S2 from Supplementary_Dataset.xlsx at final print size.

Revised 2026-10-14 after the three-reviewer check of the displays:
- drawn at the final full-page width (6.85 in, 174 mm) with text of at least 7 pt;
- no in-image titles (the legends carry them); gene symbols and P in italics;
- Figure 1 marks discoveries by shape, not by colours that mean datasets elsewhere;
- Figure 2 shows estimates only (numbers are in Table 2) and stacks H4 from zero, so
  the dotted line at 0.80 is the PP.H4 reporting threshold;
- Figure 3 has a colour bar, marks cells that stay >= 0.80 at p12 = 1e-6, and plots
  the three priors as separate points;
- Figure S2 shows the all-instrument estimate above each dataset's omissions.
Figure S1 needs IRB-protected sample data and is not redrawn here.

    python3 submission/candidate_20261011_review/reproducibility/build_figures_20261014.py
"""
from pathlib import Path

import matplotlib
import matplotlib.colors
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

PKG = Path(__file__).resolve().parents[1]
X = PKG / "Supplementary_Dataset.xlsx"
OUT = PKG / "Figures"
W = 6.85  # final full-page width, inches

ORDER = ["BBJ_Graves", "E4_GRAVES_STRICT", "UKB_hyperthyroid", "FinnGen_GO", "E4_GRAVES_OPHT_STRICT"]
SHORT = {"BBJ_Graves": "BBJ GD", "E4_GRAVES_STRICT": "FinnGen GD", "UKB_hyperthyroid": "UKB hyperthyroidism",
         "FinnGen_GO": "FinnGen GO", "E4_GRAVES_OPHT_STRICT": "FinnGen GO strict"}
COLORS = ["#28658b", "#43866b", "#777777", "#896348", "#865290"]
SEG = {"PP.H4": "#5f9fba", "PP.H3": "#d9a273", "PP.H2": "#b9c3c8", "H01": "#e9ecee"}
GENES = ["TSHR", "IGF1R", "CTLA4"]
NEWGENES = ["APOM", "NEU1", "IGLV2-18", "IGLV3-12", "TNF", "HIST1H4E", "LTA", "LTB", "C4B"]
SCREEN = 0.05 / 2544
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
INK = "#1a1f22"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "text.color": "black",
                     "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42, "savefig.dpi": 400,
                     "mathtext.default": "regular"})


def fmtp(x):
    """3 significant figures, trailing zeros kept; scientific below 0.001."""
    if x < 0.001:
        a, b = f"{x:.2e}".split("e")
        return a + "×10" + str(int(b)).translate(SUP)
    d = 2 - int(np.floor(np.log10(x)))
    return f"{x:.{d}f}"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", facecolor="white")
    plt.close(fig)


def load():
    mr = pd.read_excel(X, sheet_name="MR_primary")
    co = pd.read_excel(X, sheet_name="Colocalization")
    ins = pd.read_excel(X, sheet_name="Instruments")
    loo = pd.read_excel(X, sheet_name="LOO_anchors_plot")
    prim = co[(co.scenario == "original_reference") & np.isclose(co.p12, 1e-5)]
    low = co[(co.scenario == "original_reference") & np.isclose(co.p12, 1e-6)]
    return mr, co, prim, low, ins, loo


def get(df, g, o, gcol="gene"):
    return df[(df[gcol] == g) & (df.outcome == o)].iloc[0]


def figure1(mr, ins, hits):
    fig = plt.figure(figsize=(W, 5.3))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.95, 2.7], hspace=0.17, left=0.08, right=0.985, top=0.975, bottom=0.085)
    flow = fig.add_subplot(gs[0])
    flow.axis("off")
    flow.text(-0.075, 1.0, "A", fontweight="bold", fontsize=9, transform=flow.transAxes, va="top")
    texts = ["$\\bf{Blood\\ expression\\ instruments}$\n4,462 potentially druggable genes\n2,544 genes with instruments\n"
             "6,135 SNP–gene records",
             "$\\bf{BBJ\\ GD\\ discovery}$\n2,234 estimable genes\n13 discoveries\nThreshold $\\mathit{P}$ < 0.05/2,544",
             "$\\bf{Comparison\\ datasets}$\nUKB hyperthyroidism, FinnGen GO\nPost hoc: FinnGen GD,\n"
             "FinnGen GO strict, GTEx tissues"]
    bw, gap = 0.300, 0.04
    for i, txt in enumerate(texts):
        x = 0.008 + i * (bw + gap)
        flow.add_patch(FancyBboxPatch((x, 0.22), bw, 0.74, boxstyle="round,pad=0.006", facecolor="#f1f4f5",
                                      edgecolor="#71838b", lw=0.6, transform=flow.transAxes, clip_on=False))
        flow.text(x + 0.012, 0.59, txt, ha="left", va="center", transform=flow.transAxes, fontsize=7.2, linespacing=1.35)
        if i < 2:
            flow.annotate("", xy=(x + bw + gap - 0.006, 0.59), xytext=(x + bw + 0.008, 0.59), xycoords="axes fraction",
                          arrowprops={"arrowstyle": "->", "color": "black", "lw": 0.7})
    flow.text(0.008, 0.02, "Colocalization: 24 selected genes, centered on $\\mathit{TSHR}$ and $\\mathit{IGF1R}$, "
              "with $\\mathit{CTLA4}$ as comparator", transform=flow.transAxes, fontsize=7.5)

    ax = fig.add_subplot(gs[1])
    lead = ins.sort_values("pvalue").drop_duplicates("gene_symbol")
    d = mr[(mr.outcome == "BBJ_Graves") & mr.estimable].merge(lead[["gene_symbol", "chr", "pos_hg19"]], on="gene_symbol")
    assert len(d) == 2234
    cmax = d.groupby("chr").pos_hg19.max()
    offset, total, ticks = {}, 0.0, []
    for c in range(1, 23):
        offset[c] = total
        size = float(cmax.get(c, 1e8))
        ticks.append(total + size / 2)
        total += size + 1e7
    d["x"] = [offset[int(c)] + p for c, p in zip(d.chr, d.pos_hg19)]
    d["y"] = -np.log10(d.pvalue)
    for c, g in d.groupby("chr"):
        ax.scatter(g.x, g.y, s=3, color="#bac7ce" if c % 2 else "#92a6b2", alpha=0.75, linewidths=0, rasterized=True)
    thr = -np.log10(SCREEN)
    ax.axhline(thr, color="#555555", ls="--", lw=0.8)
    ax.text(offset[8], thr + 0.35, "$\\mathit{P}$ = 0.05/2,544", fontsize=7)
    offs = {"HLA-A": (0, 8), "HLA-DQA2": (-24, 9), "C4A": (-22, 14), "TUBB": (-30, 10), "PSMB8": (-28, -15),
            "IFNGR1": (20, 15), "HSD3B7": (0, 36), "VKORC1": (-22, 22), "PRSS36": (20, 9), "TNFSF14": (8, 14),
            "MAPKAPK5": (-10, 18), "CTLA4": (0, 8), "TSHR": (0, 8), "IGF1R": (0, 12)}
    for g in hits:
        r = d[d.gene_symbol == g].iloc[0]
        ax.scatter([r.x], [r.y], s=16, color=INK, zorder=5, linewidths=0)
    r = d[d.gene_symbol == "IGF1R"].iloc[0]
    ax.scatter([r.x], [r.y], s=26, facecolor="white", edgecolor=INK, linewidths=1.0, zorder=6)
    for g in hits + ["IGF1R"]:
        r = d[d.gene_symbol == g].iloc[0]
        dx, dy = offs[g]
        ax.annotate(g, (r.x, r.y), xytext=(dx, dy), textcoords="offset points", ha="center", fontsize=7,
                    fontstyle="italic", color="black", arrowprops={"arrowstyle": "-", "lw": 0.45, "color": "black"})
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=INK, markersize=4, label="BBJ discovery (13 genes)"),
                       Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=INK, markersize=5,
                              label="$\\mathit{IGF1R}$ (focal gene)")],
              loc="upper right", frameon=False, handletextpad=0.3)
    labels = [str(c) if c <= 18 or c in (20, 22) else "" for c in range(1, 23)]
    ax.set_xticks(ticks, labels)
    ax.set_ylim(-0.3, 25.5)
    ax.set_xlim(-1e7, total)
    ax.set_ylabel("$-\\log_{10}\\mathit{P}$, BBJ GD")
    ax.set_xlabel("Chromosome of strongest selected expression instrument (GRCh37)")
    ax.text(-0.075, 1.0, "B", fontweight="bold", fontsize=9, transform=ax.transAxes, va="top")
    ax.grid(axis="y", alpha=0.18, lw=0.5)
    save(fig, "Figure1")


def figure2(mr, prim):
    fig = plt.figure(figsize=(W, 4.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=0.08, left=0.255, right=0.985, top=0.93, bottom=0.17)
    ax = fig.add_subplot(gs[0])
    bar = fig.add_subplot(gs[1], sharey=ax)
    rows = [(g, o) for g in GENES for o in ORDER]
    for i, (g, o) in enumerate(rows):
        m, c = get(mr, g, o, "gene_symbol"), get(prim, g, o)
        col = COLORS[ORDER.index(o)]
        ax.errorbar(m.OR, i, xerr=[[m.OR - m.CI_lower], [m.CI_upper - m.OR]], fmt="o", color=col, ecolor=col,
                    elinewidth=0.9, capsize=2, markersize=3.5)
        if o == ORDER[0]:
            ax.text(-0.585, i, g, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7.5, fontstyle="italic")
        ax.text(-0.03, i, SHORT[o], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7.5)
        left = 0
        for key in ("PP.H4", "PP.H3", "PP.H2"):
            bar.barh(i, c[key], left=left, height=0.62, color=SEG[key], edgecolor="white", lw=0.4)
            left += c[key]
        bar.barh(i, c["PP.H0"] + c["PP.H1"], left=left, height=0.62, color=SEG["H01"], edgecolor="white", lw=0.4)
    for a in (ax, bar):
        a.set_ylim(14.7, -0.7)
        for y in (4.5, 9.5):
            a.axhline(y, color="#c9c9c9", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.025, 4)
    ax.set_xticks([0.03, 0.1, 0.3, 1, 3], ["0.03", "0.1", "0.3", "1", "3"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.axvline(1, color="black", ls="--", lw=0.6)
    ax.set_yticks(range(15), [""] * 15)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("OR (95% CI) per approximately 1 SD\nhigher genetically predicted expression")
    ax.set_title("A  MR estimate", loc="left", fontweight="bold")
    bar.tick_params(axis="y", left=False, labelleft=False)
    bar.set_xlim(0, 1)
    bar.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    bar.axvline(0.8, color="black", ls=":", lw=0.8)
    bar.set_xlabel("Posterior probability")
    bar.set_title("B  Colocalization", loc="left", fontweight="bold")
    fig.legend(handles=[Patch(facecolor=SEG[k], edgecolor="#9aa4a8", lw=0.4, label=l) for k, l in [
        ("PP.H4", "H4: shared variant"), ("PP.H3", "H3: distinct variants"), ("PP.H2", "H2: expression only"),
        ("H01", "H0/H1")]], loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.62, 0.0),
        handlelength=1.2, columnspacing=1.2)
    save(fig, "Figure2")


def figure3(co, prim, low, hits):
    fig = plt.figure(figsize=(W, 7.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], wspace=0.42, left=0.115, right=0.985, top=0.955, bottom=0.075)
    hm = fig.add_subplot(gs[0])
    pr = fig.add_subplot(gs[1])
    order = GENES + [g for g in hits if g not in GENES] + NEWGENES + ["HLA-G"]
    assert len(order) == 24
    m = np.array([[get(prim, g, o)["PP.H4"] for o in ORDER] for g in order])
    m6 = np.array([[get(low, g, o)["PP.H4"] for o in ORDER] for g in order])
    # Truncated Blues keeps every cell light enough for black lettering.
    cmap = matplotlib.colors.ListedColormap(plt.cm.Blues(np.linspace(0, 0.55, 256)))
    im = hm.imshow(m, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i in range(24):
        for j in range(5):
            v = m[i, j]
            hm.text(j, i, "<0.001" if v < 0.001 else f"{v:.3f}", ha="center", va="center", fontsize=7, color="black")
            if v >= 0.8:
                hm.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="black", linewidth=1.0))
            if m6[i, j] >= 0.8:
                hm.plot(j + 0.38, i - 0.30, marker="o", markersize=2.6, color="black")
    hm.set_xticks(range(5), ["BBJ\nGD", "FinnGen\nGD", "UKB\nhyperthyroidism", "FinnGen\nGO", "FinnGen\nGO strict"], fontsize=7)
    hm.set_yticks(range(24), order, fontstyle="italic", fontsize=7.5)
    hm.tick_params(length=0)
    for s in hm.spines.values():
        s.set_visible(False)
    hm.set_title("A  Shared-variant probability (PP.H4)", loc="left", fontweight="bold", pad=8)
    for y in (2.5, 13.5, 22.5):
        hm.axhline(y, color="black", lw=0.8)
    cb = fig.colorbar(im, ax=hm, orientation="horizontal", fraction=0.035, pad=0.065, aspect=30)
    cb.set_label("PP.H4 (p12 = 10⁻⁵)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    priors = [1e-5, 5e-6, 1e-6]
    for k, (o, c) in enumerate(zip(ORDER, COLORS)):
        z = co[(co.gene == "TSHR") & (co.outcome == o) & (co.scenario == "original_reference")]
        vals = [z[np.isclose(z.p12, p)]["PP.H4"].iloc[0] for p in priors]
        pr.scatter(np.arange(3) + (k - 2) * 0.07, vals, s=16, color=c, label=SHORT[o], zorder=3)
    pr.axhline(0.8, color="black", ls="--", lw=0.7)
    pr.set_ylim(0, 1.05)
    pr.set_xlim(-0.4, 2.4)
    pr.set_xticks(range(3), ["10⁻⁵", "5×10⁻⁶", "10⁻⁶"])
    pr.set_xlabel("Shared-variant prior (p12)")
    pr.set_ylabel("PP.H4")
    pr.set_title("B  $\\mathbfit{TSHR}$ prior sensitivity", loc="left", fontweight="bold", pad=8)
    pr.legend(frameon=False, loc="lower left", fontsize=7, handletextpad=0.2, borderaxespad=0.2)
    pr.grid(axis="y", alpha=0.15, lw=0.5)
    save(fig, "Figure3")
    return [(order[i], ORDER[j]) for i in range(24) for j in range(5) if m6[i, j] >= 0.8]


def figureS2(mr, loo):
    groups = {"IGF1R": ORDER, "CTLA4": [o for o in ORDER if o != "BBJ_Graves"]}
    nrows = {g: sum(1 + loo[(loo.gene == g) & (loo.outcome == o)].shape[0] for o in groups[g]) for g in groups}
    assert len(loo) == 25 and nrows == {"IGF1R": 22, "CTLA4": 12}
    fig, axes = plt.subplots(2, 1, figsize=(W, 7.4), gridspec_kw={"height_ratios": [nrows["IGF1R"] + 2, nrows["CTLA4"] + 2]})
    fig.subplots_adjust(left=0.235, right=0.86, hspace=0.30, top=0.965, bottom=0.07)
    for ax, g, letter in zip(axes, ["IGF1R", "CTLA4"], ["A", "B"]):
        labels, bold, y = [], [], 0
        for o in groups[g]:
            full = get(mr, g, o, "gene_symbol")
            col = COLORS[ORDER.index(o)]
            ax.errorbar(full.OR, y, xerr=[[full.OR - full.CI_lower], [full.CI_upper - full.OR]], fmt="D", color=col,
                        markerfacecolor="white", markersize=3.8, capsize=2, elinewidth=0.9)
            ax.text(1.02, y, fmtp(full.pvalue), transform=ax.get_yaxis_transform(), va="center", fontsize=7,
                    fontweight="bold" if full.pvalue < 0.05 else "normal")
            labels.append(f"{SHORT[o]}, all IVs")
            bold.append(True)
            y += 1
            for _, r in loo[(loo.gene == g) & (loo.outcome == o)].sort_values("removed_snp").iterrows():
                ax.errorbar(r.OR, y, xerr=[[r.OR - r.CI_lower], [r.CI_upper - r.OR]], fmt="o", color=col,
                            markersize=3, capsize=2, elinewidth=0.8)
                ax.text(1.02, y, fmtp(r.pvalue), transform=ax.get_yaxis_transform(), va="center", fontsize=7,
                        fontweight="bold" if r.pvalue < 0.05 else "normal")
                labels.append(f"omit {r.removed_snp}")
                bold.append(False)
                y += 1
        ax.set_yticks(range(len(labels)), labels, fontsize=7.2)
        ax.tick_params(axis="y", length=0, pad=4)
        for k, b in enumerate(bold):
            if b and k:
                ax.axhline(k - 0.5, color="#d0d0d0", lw=0.5)
        ax.set_ylim(len(labels) - 0.4, -0.6)
        ax.set_xscale("log")
        ticks = [0.5, 1, 2, 4] if g == "IGF1R" else [0.03, 0.1, 0.3, 1, 3, 10]
        ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.axvline(1, color="black", ls="--", lw=0.6)
        ax.set_xlabel("OR (95% CI) per approximately 1 SD higher genetically predicted expression")
        ax.set_title(f"{letter}  $\\mathbfit{{{g}}}$", loc="left", fontweight="bold")
        ax.text(1.02, 1.005, "$\\mathit{P}$", transform=ax.transAxes, fontsize=7.5, va="bottom")
        ax.grid(axis="x", alpha=0.15, lw=0.5)
    save(fig, "FigureS2")


def main():
    mr, co, prim, low, ins, loo = load()
    hits = mr[(mr.outcome == "BBJ_Graves") & mr.P_lt_0_05_over_2544].sort_values("pvalue").gene_symbol.tolist()
    assert len(hits) == 13
    figure1(mr, ins, hits)
    figure2(mr, prim)
    robust = figure3(co, prim, low, hits)
    figureS2(mr, loo)
    print("Figures 1-3 and S2 written; cells >=0.80 at p12 = 1e-6:", robust)


if __name__ == "__main__":
    main()
