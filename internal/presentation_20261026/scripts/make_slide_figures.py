"""Slide versions (16:9 deck, 10 x 5.625 in) of Figures 1B, 2, 3A and S2, drawn from the workbook.

Same data and logic as reproducibility/build_figures_20261014.py, re-laid out at the size they are
placed on the slide so text is 10-13 pt. Carlito is metric-compatible with the deck's Calibri.
"""
import importlib.util
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # importing the package script must not leave __pycache__ in the delivery

import matplotlib
import matplotlib.colors
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

BF = Path("/home/user/TEDGWAS/submission/candidate_20261011_review/reproducibility/build_figures_20261014.py")
spec = importlib.util.spec_from_file_location("bf", BF)
bf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bf)
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "Carlito", "font.size": 12, "axes.titlesize": 12.5, "axes.labelsize": 12,
                     "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11, "savefig.dpi": 300,
                     "mathtext.fontset": "custom", "mathtext.rm": "Carlito", "mathtext.it": "Carlito:italic",
                     "mathtext.bf": "Carlito:bold", "mathtext.bfit": "Carlito:bold:italic"})
ORDER, SHORT, COLORS, SEG, GENES, INK = bf.ORDER, bf.SHORT, bf.COLORS, bf.SEG, bf.GENES, bf.INK
mr, co, prim, low, ins, loo = bf.load()
hits = mr[(mr.outcome == "BBJ_Graves") & mr.P_lt_0_05_over_2544].sort_values("pvalue").gene_symbol.tolist()
assert len(hits) == 13


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", facecolor="white")
    plt.close(fig)


def manhattan():
    fig = plt.figure(figsize=(6.0, 4.3))
    ax = fig.add_axes([0.105, 0.135, 0.885, 0.845])
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
        ax.scatter(g.x, g.y, s=5, color="#bac7ce" if c % 2 else "#92a6b2", alpha=0.8, linewidths=0, rasterized=True)
    thr = -np.log10(bf.SCREEN)
    ax.axhline(thr, color="#555555", ls="--", lw=0.9)
    ax.text(offset[8], thr + 0.45, "$\\mathit{P}$ = 0.05/2,544", fontsize=11)
    offs = {"HLA-A": (0, 9), "HLA-DQA2": (-30, 10), "C4A": (-26, 14), "TUBB": (-34, 10), "PSMB8": (-30, -18),
            "IFNGR1": (26, 16), "HSD3B7": (0, 42), "VKORC1": (-30, 24), "PRSS36": (-28, -5), "TNFSF14": (6, 30),
            "MAPKAPK5": (-12, 20), "CTLA4": (0, 9), "TSHR": (0, 9), "IGF1R": (0, 14)}
    for g in hits:
        r = d[d.gene_symbol == g].iloc[0]
        ax.scatter([r.x], [r.y], s=26, color=INK, zorder=5, linewidths=0)
    r = d[d.gene_symbol == "IGF1R"].iloc[0]
    ax.scatter([r.x], [r.y], s=40, facecolor="white", edgecolor=INK, linewidths=1.2, zorder=6)
    for g in hits + ["IGF1R"]:
        r = d[d.gene_symbol == g].iloc[0]
        dx, dy = offs[g]
        big = g in ("TSHR", "IGF1R", "CTLA4")
        ax.annotate(g, (r.x, r.y), xytext=(dx, dy), textcoords="offset points", ha="center",
                    fontsize=12 if big else 10.5, fontstyle="italic", fontweight="bold" if big else "normal",
                    color="black", arrowprops={"arrowstyle": "-", "lw": 0.5, "color": "black"})
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=INK, markersize=5.5, label="BBJ discovery (13 genes)"),
                       Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=INK, markersize=6.5,
                              label="$\\mathit{IGF1R}$")],
              loc="upper right", frameon=False, handletextpad=0.3, fontsize=11)
    labels = [str(c) if c <= 12 or c in (14, 16, 18, 20, 22) else "" for c in range(1, 23)]
    ax.set_xticks(ticks, labels)
    ax.set_ylim(-0.3, 25.5)
    ax.set_xlim(-1e7, total)
    ax.set_ylabel("$-\\log_{10}\\mathit{P}$, BBJ GD")
    ax.set_xlabel("Chromosome of strongest expression instrument")
    ax.grid(axis="y", alpha=0.18, lw=0.5)
    save(fig, "slide_fig1_manhattan")


def overview():
    fig = plt.figure(figsize=(8.8, 4.35))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.06, left=0.235, right=0.985, top=0.93, bottom=0.2)
    ax = fig.add_subplot(gs[0])
    bar = fig.add_subplot(gs[1], sharey=ax)
    rows = [(g, o) for g in GENES for o in ORDER]
    for i, (g, o) in enumerate(rows):
        m, c = bf.get(mr, g, o, "gene_symbol"), bf.get(prim, g, o)
        col = COLORS[ORDER.index(o)]
        ax.errorbar(m.OR, i, xerr=[[m.OR - m.CI_lower], [m.CI_upper - m.OR]], fmt="o", color=col, ecolor=col,
                    elinewidth=1.2, capsize=2.5, markersize=5)
        if o == ORDER[0]:
            ax.text(-0.545, i, g, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=12.5,
                    fontstyle="italic", fontweight="bold")
        ax.text(-0.025, i, SHORT[o], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=11)
        left = 0
        for key in ("PP.H4", "PP.H3", "PP.H2"):
            bar.barh(i, c[key], left=left, height=0.66, color=SEG[key], edgecolor="white", lw=0.4)
            left += c[key]
        bar.barh(i, c["PP.H0"] + c["PP.H1"], left=left, height=0.66, color=SEG["H01"], edgecolor="white", lw=0.4)
    for a in (ax, bar):
        a.set_ylim(14.6, -0.6)
        for y in (4.5, 9.5):
            a.axhline(y, color="#b5b5b5", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.025, 4)
    ax.set_xticks([0.03, 0.1, 0.3, 1, 3], ["0.03", "0.1", "0.3", "1", "3"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.axvline(1, color="black", ls="--", lw=0.8)
    ax.set_yticks(range(15), [""] * 15)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("OR (95% CI) per ~1 SD higher expression")
    ax.set_title("A  MR estimate", loc="left", fontweight="bold")
    bar.tick_params(axis="y", left=False, labelleft=False)
    bar.set_xlim(0, 1)
    bar.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0", "0.2", "0.4", "0.6", "0.8", "1"])
    bar.axvline(0.8, color="black", ls=":", lw=1.0)
    bar.set_xlabel("Posterior probability")
    bar.set_title("B  Colocalization", loc="left", fontweight="bold")
    fig.legend(handles=[Patch(facecolor=SEG[k], edgecolor="#9aa4a8", lw=0.4, label=l) for k, l in [
        ("PP.H4", "H4: shared variant"), ("PP.H3", "H3: distinct variants"), ("PP.H2", "H2: expression only"),
        ("H01", "H0/H1")]], loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.61, -0.01),
        handlelength=1.2, columnspacing=1.4, fontsize=11)
    save(fig, "slide_fig2_overview")


def heatmap():
    order = GENES + [g for g in hits if g not in GENES] + bf.NEWGENES + ["HLA-G"]
    assert len(order) == 24
    m = np.array([[bf.get(prim, g, o)["PP.H4"] for g in order] for o in ORDER])
    m6 = np.array([[bf.get(low, g, o)["PP.H4"] for g in order] for o in ORDER])
    fig = plt.figure(figsize=(8.9, 3.75))
    hm = fig.add_axes([0.165, 0.25, 0.83, 0.67])
    cmap = matplotlib.colors.ListedColormap(plt.cm.Blues(np.linspace(0, 0.55, 256)))
    hm.imshow(m, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i in range(5):
        for j in range(24):
            v = m[i, j]
            hm.text(j, i, "" if v < 0.001 else (f"{v:.3f}" if v < 0.01 else f"{v:.2f}"), ha="center", va="center", fontsize=8.2,
                    color="black")
            if v >= 0.8:
                hm.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="black", linewidth=1.3))
            if m6[i, j] >= 0.8:
                hm.plot(j + 0.36, i - 0.30, marker="o", markersize=2.8, color="black")
    hm.set_yticks(range(5), [SHORT[o] for o in ORDER], fontsize=11)
    hm.set_xticks(range(24), order, fontstyle="italic", fontsize=10.5, rotation=55, ha="right", rotation_mode="anchor")
    hm.tick_params(length=0)
    for s in hm.spines.values():
        s.set_visible(False)
    for x in (2.5, 13.5, 22.5):
        hm.axvline(x, color="black", lw=1.0)
    for x, t in [(1, "Focal and comparator"), (8, "Other BBJ discoveries"), (18, "Additional FinnGen GD screen"), (23, "")]:
        hm.text(x, -0.62, t, ha="center", va="bottom", fontsize=10.5, fontweight="bold")
    save(fig, "slide_fig3_heatmap")
    return order, m


def loo_plot():
    groups = {"IGF1R": ORDER, "CTLA4": [o for o in ORDER if o != "BBJ_Graves"]}
    fig, axes = plt.subplots(1, 2, figsize=(8.9, 4.35), gridspec_kw={"wspace": 1.2})
    fig.subplots_adjust(left=0.245, right=0.925, top=0.93, bottom=0.12)
    for ax, g, letter in zip(axes, ["IGF1R", "CTLA4"], ["A", "B"]):
        labels, bold, y = [], [], 0
        for o in groups[g]:
            full = bf.get(mr, g, o, "gene_symbol")
            col = COLORS[ORDER.index(o)]
            ax.errorbar(full.OR, y, xerr=[[full.OR - full.CI_lower], [full.CI_upper - full.OR]], fmt="D", color=col,
                        markerfacecolor="white", markersize=4.5, capsize=2, elinewidth=1.0)
            ax.text(1.03, y, bf.fmtp(full.pvalue), transform=ax.get_yaxis_transform(), va="center", fontsize=10,
                    fontweight="bold" if full.pvalue < 0.05 else "normal")
            labels.append(f"{SHORT[o]} (all IVs)")
            bold.append(True)
            y += 1
            for _, r in loo[(loo.gene == g) & (loo.outcome == o)].sort_values("removed_snp").iterrows():
                ax.errorbar(r.OR, y, xerr=[[r.OR - r.CI_lower], [r.CI_upper - r.OR]], fmt="o", color=col,
                            markersize=3.5, capsize=2, elinewidth=0.9)
                ax.text(1.03, y, bf.fmtp(r.pvalue), transform=ax.get_yaxis_transform(), va="center", fontsize=10,
                        fontweight="bold" if r.pvalue < 0.05 else "normal")
                labels.append(f"– {r.removed_snp}")
                bold.append(False)
                y += 1
        ax.set_yticks(range(len(labels)), labels, fontsize=10)
        for t, b in zip(ax.get_yticklabels(), bold):
            t.set_fontweight("bold" if b else "normal")
        ax.tick_params(axis="y", length=0, pad=3)
        for k, b in enumerate(bold):
            if b and k:
                ax.axhline(k - 0.5, color="#c8c8c8", lw=0.6)
        ax.set_ylim(len(labels) - 0.4, -0.6)
        ax.set_xscale("log")
        ticks = [0.5, 1, 2, 4] if g == "IGF1R" else [0.03, 0.1, 0.3, 1, 3, 10]
        ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.axvline(1, color="black", ls="--", lw=0.7)
        ax.set_xlabel("OR (95% CI) per ~1 SD", fontsize=11)
        ax.set_title(f"{letter}  $\\mathbfit{{{g}}}$", loc="left", fontweight="bold", fontsize=12.5)
        ax.text(1.03, 1.005, "$\\mathit{P}$", transform=ax.transAxes, fontsize=10, va="bottom")
        ax.grid(axis="x", alpha=0.15, lw=0.5)
    save(fig, "slide_figS2_loo")


manhattan()
overview()
order, m = heatmap()
loo_plot()
print("written to", OUT)
