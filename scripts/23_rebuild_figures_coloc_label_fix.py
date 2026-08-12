#!/usr/bin/env python3
"""
Rebuild Figure 3 (composite) and Figure S2 with the CORRECTED colocalization legend.

Why: taskE_01_coloc.R / taskE_02_coloc_candidates.R call coloc.abf with
dataset1 = outcome GWAS and dataset2 = cis-eQTL. In coloc.abf, H1 = trait-1-only and
H2 = trait-2-only, so PP.H2 = "cis-eQTL association only". The original ggplot builders
labelled it "H2 (GWAS assoc. only)", which is inverted. This rebuild fixes the label
(and shows IGF1R's tissue padj as >0.99 rather than 1.0e+00, matching the text).

All plotted values are read from the locked v5 result tables; nothing is hard-coded.
Outputs overwrite the canonical PNGs and the submission copies.
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLOC1 = f"{ROOT}/TrackA_MR/v5_upgrade/05_coloc_candidates/TaskE_01_coloc_results_v1.csv"
COLOC2 = f"{ROOT}/TrackA_MR/v5_upgrade/05_coloc_candidates/TaskE_02_candidate_coloc_v1.csv"
MR     = f"{ROOT}/TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv"
TISSUE = f"{ROOT}/TrackA_MR/v5_upgrade/06_tissue_integration/TaskF_01_backbone_tissue_inhouse_v1.csv"
FIGDIR = f"{ROOT}/TrackA_MR/v5_upgrade/07_manuscript/figures"
SUBDIR = f"{ROOT}/submission/figures"

C_H2, C_H3, C_H4 = "#F4B183", "#BFBFBF", "#70AD47"      # keep original palette
GENE_COL = {"TSHR": "#C00000", "IGF1R": "#2E75B6", "CTLA4": "#548235"}
STRIP = "#D9D9D9"
LBL_H2 = "H2 (eQTL assoc. only)"                          # <- CORRECTED
LBL_H3, LBL_H4 = "H3 (distinct)", "H4 (shared)"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.edgecolor": "#404040", "axes.linewidth": 0.8})

def rows(p): return list(csv.DictReader(open(p)))

def strip_title(ax, text):
    ax.set_title(text, fontsize=9, fontstyle="italic", fontweight="bold", pad=7,
                 bbox=dict(facecolor=STRIP, edgecolor="#404040", linewidth=0.8,
                           boxstyle="square,pad=0.32"))

def stacked_coloc(ax, h2, h3, h4, xlabels, ylab=None, ymax=1.0):
    x = range(len(xlabels))
    ax.bar(x, h2, color=C_H2, width=0.62, edgecolor="none")
    ax.bar(x, h3, bottom=h2, color=C_H3, width=0.62, edgecolor="none")
    ax.bar(x, h4, bottom=[a + b for a, b in zip(h2, h3)], color=C_H4, width=0.62, edgecolor="none")
    ax.axhline(0.8, ls="--", lw=1.0, color="#404040")
    ax.set_xticks(list(x)); ax.set_xticklabels(xlabels)
    ax.set_ylim(0, ymax); ax.set_yticks([0, .25, .5, .75, 1.0])
    ax.set_yticklabels([f"{v:.2f}" for v in [0, .25, .5, .75, 1.0]])
    ax.grid(axis="y", color="#EDEDED", lw=0.7); ax.set_axisbelow(True)
    if ylab: ax.set_ylabel(ylab)
    else: ax.set_yticklabels([])

# ---------------------------------------------------------------- Figure 3
def figure3():
    co = {(r["gene_symbol"], r["outcome"]): r for r in rows(COLOC1)}
    mr = {(r["gene_symbol"], r["outcome"]): r for r in rows(MR)
          if r["method"] in ("Inverse variance weighted", "Wald ratio")}
    ti = {r["gene"]: r for r in rows(TISSUE)}
    genes = ["TSHR", "IGF1R", "CTLA4"]
    oc2 = [("BBJ_Graves", "BBJ"), ("FinnGen_GO", "FinnGen GO")]
    oc3 = [("BBJ_Graves", "BBJ Graves"), ("UKB_hyperthyroid", "UKB hyperthyroid"),
           ("FinnGen_GO", "FinnGen GO")]

    fig = plt.figure(figsize=(9.2, 11.6), dpi=300)
    L, R = 0.145, 0.975
    gsA = fig.add_gridspec(1, 3, left=L, right=R, top=0.905, bottom=0.722, wspace=0.16)
    gsB = fig.add_gridspec(1, 3, left=L, right=R, top=0.590, bottom=0.365, wspace=0.16)
    gsC = fig.add_gridspec(1, 3, left=L, right=R, top=0.255, bottom=0.055, wspace=0.30)

    # ---- Panel A
    fig.text(0.015, 0.955, "A", fontsize=15, fontweight="bold")
    fig.text(0.050, 0.955, "Colocalization", fontsize=13, fontweight="bold")
    for j, g in enumerate(genes):
        ax = fig.add_subplot(gsA[0, j])
        h2 = [float(co[(g, o)]["pp_h2"]) for o, _ in oc2]
        h3 = [float(co[(g, o)]["pp_h3"]) for o, _ in oc2]
        h4 = [float(co[(g, o)]["pp_h4"]) for o, _ in oc2]
        stacked_coloc(ax, h2, h3, h4, [lab for _, lab in oc2],
                      ylab="Posterior probability" if j == 0 else None)
        strip_title(ax, g)
    fig.legend(handles=[Patch(facecolor=C_H2, label=LBL_H2),
                        Patch(facecolor=C_H3, label=LBL_H3),
                        Patch(facecolor=C_H4, label=LBL_H4)],
               loc="center", bbox_to_anchor=(0.55, 0.655), ncol=3,
               frameon=False, title="Coloc PP", title_fontproperties={"weight": "bold"})

    # ---- Panel B
    fig.text(0.015, 0.612, "B", fontsize=15, fontweight="bold")
    fig.text(0.050, 0.612, "MR effect across outcomes", fontsize=13, fontweight="bold")
    for j, g in enumerate(genes):
        ax = fig.add_subplot(gsB[0, j])
        ys = list(range(len(oc3)))[::-1]
        for y, (o, _) in zip(ys, oc3):
            r = mr[(g, o)]; b, se = float(r["beta"]), float(r["se"])
            ax.errorbar(b, y, xerr=1.96 * se, fmt="o", ms=5.5, color=GENE_COL[g],
                        ecolor=GENE_COL[g], elinewidth=1.4, capsize=3.2, capthick=1.4)
        ax.axvline(0, ls=":", lw=1.1, color="#808080")
        ax.set_yticks(ys); ax.set_yticklabels([lab for _, lab in oc3] if j == 0 else [])
        ax.set_ylim(-0.6, len(oc3) - 0.4)
        ax.grid(axis="x", color="#EDEDED", lw=0.7); ax.set_axisbelow(True)
        strip_title(ax, g)
        if j == 1: ax.set_xlabel("MR β (log-odds, 95% CI)")
    # ---- Panel C
    fig.text(0.015, 0.292, "C", fontsize=15, fontweight="bold")
    fig.text(0.050, 0.292, "Orbital tissue expression (in-house RNA-seq)",
             fontsize=13, fontweight="bold")
    for j, g in enumerate(genes):
        ax = fig.add_subplot(gsC[0, j])
        t = ti[g]
        ted = [float(v) for v in t["ted_per_sample"].split(";")]
        ctrl = [float(t["ctrl_tpm"])]
        bp = ax.boxplot([ted], positions=[0], widths=0.5, patch_artist=True,
                        medianprops=dict(color=GENE_COL[g], lw=1.8),
                        boxprops=dict(facecolor="none", edgecolor=GENE_COL[g], lw=1.4),
                        whiskerprops=dict(color=GENE_COL[g], lw=1.4),
                        capprops=dict(color=GENE_COL[g], lw=1.4), showfliers=False)
        ax.scatter([0] * len(ted), ted, color=GENE_COL[g], s=26, zorder=3)
        ax.scatter([1] * len(ctrl), ctrl, color=GENE_COL[g], s=26, zorder=3)
        padj = float(t["deseq2_padj"])
        ptxt = "padj > 0.99" if padj > 0.99 else f"padj = {padj:.2g}"
        ax.text(0.97, 0.97, f"log2FC = {float(t['deseq2_log2fc']):+.2f}\n{ptxt}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["TED", "Control"])
        ax.set_xlim(-0.6, 1.6)
        ax.grid(axis="y", color="#EDEDED", lw=0.7); ax.set_axisbelow(True)
        if j == 0: ax.set_ylabel("TPM")
        strip_title(ax, g)

    for p in (f"{FIGDIR}/Figure3_TSHR_vs_IGF1R_composite.png", f"{SUBDIR}/Figure3.png"):
        fig.savefig(p, dpi=300, facecolor="white")
    plt.close(fig); print("Figure 3 rebuilt (H2 label corrected; IGF1R padj >0.99)")

# ---------------------------------------------------------------- Figure S2
def figureS2():
    rr = rows(COLOC2)
    genes, seen = [], set()
    for r in rr:
        if r["gene_symbol"] not in seen:
            seen.add(r["gene_symbol"]); genes.append(r["gene_symbol"])
    co = {(r["gene_symbol"], r["outcome"]): r for r in rr}
    oc2 = [("BBJ_Graves", "BBJ"), ("FinnGen_GO", "FinnGen GO")]

    fig, axes = plt.subplots(1, len(genes), figsize=(2.05 * len(genes), 3.5), dpi=300)
    for j, (g, ax) in enumerate(zip(genes, axes)):
        h2 = [float(co[(g, o)]["pp_h2"]) for o, _ in oc2]
        h3 = [float(co[(g, o)]["pp_h3"]) for o, _ in oc2]
        h4 = [float(co[(g, o)]["pp_h4"]) for o, _ in oc2]
        stacked_coloc(ax, h2, h3, h4, [lab for _, lab in oc2],
                      ylab="Posterior probability" if j == 0 else None)
        ax.tick_params(axis="x", labelrotation=30)
        strip_title(ax, g)
    fig.legend(handles=[Patch(facecolor=C_H2, label=LBL_H2),
                        Patch(facecolor=C_H3, label=LBL_H3),
                        Patch(facecolor=C_H4, label=LBL_H4)],
               loc="lower center", ncol=3, frameon=False, title="Coloc PP",
               title_fontproperties={"weight": "bold"}, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.14, 1, 1])
    for p in (f"{FIGDIR}/FigureS2_candidate_coloc_failure.png", f"{SUBDIR}/FigureS2.png"):
        fig.savefig(p, dpi=300, facecolor="white")
    plt.close(fig); print("Figure S2 rebuilt (H2 label corrected)")

if __name__ == "__main__":
    figure3(); figureS2()
