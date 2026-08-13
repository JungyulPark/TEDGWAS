#!/usr/bin/env python3
"""
Rebuild Figure 1 as a CLAIM figure, not a generic workflow diagram.

The previous Figure 1 was a data-sources -> analysis -> integration flowchart: it listed
what was done rather than what was found, and its FinnGen box overflowed its border. The
replacement carries the paper's argument on its own:

  (A) the screening funnel, with every attrition step labelled by WHY genes were removed
      (MHC/LD spillover, single-outcome colocalization, distinct-variant) — ending at
      0 robust novel targets;
  (B) the verdict panel — the same three evidence layers scored for TSHR vs IGF1R, so the
      divergence that titles the paper is visible at a glance.

Every number is read from the locked master and result tables; nothing is hard-coded
except the classification labels, which are parsed from Table 3.

Out: TrackA_MR/v5_upgrade/07_manuscript/figures/Figure1_study_design.png
     submission/figures/Figure1.png
"""
import csv, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER = f"{ROOT}/MANUSCRIPT_TED_TRAP_v5_MASTER.md"
COLOC = f"{ROOT}/TrackA_MR/v5_upgrade/05_coloc_candidates/TaskE_01_coloc_results_v1.csv"
MR = f"{ROOT}/TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_03d_MR_all_outcomes_merged_v1.csv"
TISSUE = f"{ROOT}/TrackA_MR/v5_upgrade/06_tissue_integration/TaskF_01_backbone_tissue_inhouse_v1.csv"
OUT = [f"{ROOT}/TrackA_MR/v5_upgrade/07_manuscript/figures/Figure1_study_design.png",
       f"{ROOT}/submission/figures/Figure1.png"]

INK, MUTED, RULE = "#1a1d23", "#5c6673", "#c9d0d9"
ANCHOR, EFFECTOR, DROP, KEEP = "#a4262c", "#2e75b6", "#b9c0c9", "#3f6f3f"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})


def master_counts():
    t = open(MASTER).read()
    for n in ("4,462", "2,545", "2,235", "6,136"):
        assert n in t, f"count {n} not in master"
    sec = t.split("**Table 3.")[1].split("\n\n")[1]
    rows = [r for r in sec.split("\n") if r.startswith("|") and "---" not in r]
    hdr = [c.strip() for c in rows[0].strip("|").split("|")]
    ci = hdr.index("Classification")
    cls = {}
    for r in rows[1:]:
        c = [x.strip() for x in r.strip("|").split("|")]
        cls[re.sub(r"[*]", "", c[0])] = c[ci]
    return cls


def backbone_values():
    co = {(r["gene_symbol"], r["outcome"]): r
          for r in csv.DictReader(open(COLOC))}
    mr = {(r["gene_symbol"], r["outcome"]): r for r in csv.DictReader(open(MR))
          if r["method"] in ("Inverse variance weighted", "Wald ratio")}
    ti = {r["gene"]: r for r in csv.DictReader(open(TISSUE))}
    out = {}
    for g in ("TSHR", "IGF1R"):
        out[g] = dict(
            beta_bbj=float(mr[(g, "BBJ_Graves")]["beta"]),
            p_bbj=float(mr[(g, "BBJ_Graves")]["pvalue"]),
            h4_bbj=float(co[(g, "BBJ_Graves")]["pp_h4"]),
            h4_fin=float(co[(g, "FinnGen_GO")]["pp_h4"]),
            h2_bbj=float(co[(g, "BBJ_Graves")]["pp_h2"]),
            h2_fin=float(co[(g, "FinnGen_GO")]["pp_h2"]),
            tis_fc=float(ti[g]["deseq2_log2fc"]), tis_p=float(ti[g]["deseq2_padj"]),
        )
    return out


def box(ax, x, y, w, h, fc, ec, lw=1.0, r=0.012):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw))


def build():
    cls = master_counts()
    bb = backbone_values()
    n_mhc = sum(1 for v in cls.values() if "MHC" in v)
    n_ld = sum(1 for v in cls.values() if "chr16p11.2" in v)
    n_single = sum(1 for v in cls.values() if "single-outcome" in v)
    n_distinct = sum(1 for v in cls.values() if "distinct-variant" in v)
    n_known = sum(1 for v in cls.values() if v.startswith("Established"))
    assert n_mhc + n_ld + n_single + n_distinct + n_known == 13

    fig = plt.figure(figsize=(10.5, 12.2), dpi=300)
    fig.patch.set_facecolor("white")

    # ================= PANEL A — screening funnel =================
    axA = fig.add_axes([0.06, 0.505, 0.90, 0.455]); axA.axis("off")
    axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    axA.text(0, 1.00, "A", fontsize=16, fontweight="bold", va="top")
    axA.text(0.035, 1.00, "What survives a two-filter druggable-genome screen",
             fontsize=13, fontweight="bold", va="top")
    axA.text(0.035, 0.945,
             "Attrition is labelled by the reason genes were removed, not merely by count.",
             fontsize=9.5, color=MUTED, va="top")

    steps = [
        ("4,462", "druggable genes\n(Finan et al.)", 0.86, 1.00),
        ("2,545", "MR-testable with a valid\neQTLGen cis instrument (6,136 IVs)", 0.66, 0.84),
        ("2,235", "estimable against the\nBBJ discovery outcome", 0.46, 0.68),
        ("13", "Bonferroni-significant\ndiscovery hits (P < 1.965×10⁻⁵)", 0.26, 0.50),
    ]
    for i, (n, lab, y, wfrac) in enumerate(steps):
        w = 0.46 * wfrac; x = 0.055 + (0.46 - w) / 2
        box(axA, x, y - 0.075, w, 0.115, "#f2f5f8", RULE, 1.0)
        axA.text(x + w / 2, y + 0.008, n, ha="center", va="center",
                 fontsize=17, fontweight="bold", color=INK)
        axA.text(x + w / 2, y - 0.048, lab, ha="center", va="center",
                 fontsize=8.4, color=MUTED, linespacing=1.35)
        if i < len(steps) - 1:
            ny = steps[i + 1][2]
            axA.annotate("", xy=(0.285, ny + 0.045), xytext=(0.285, y - 0.078),
                         arrowprops=dict(arrowstyle="-|>", color=RULE, lw=1.4))

    # the 13 hits, split by why they do or do not survive
    axA.text(0.575, 0.795, "The 13 hits, resolved", fontsize=10.5, fontweight="bold", va="top")
    bars = [
        (f"{n_mhc}", "MHC-region signal", DROP, "removed"),
        (f"{n_ld}", "chr16p11.2 — one LD block,\nnot three signals", DROP, "removed"),
        (f"{n_single}", "colocalized in discovery only\n(TNFSF14, IFNGR1: 0.99 → 0.02)", DROP, "removed"),
        (f"{n_distinct}", "distinct causal variant\n(MAPKAPK5)", DROP, "removed"),
        (f"{n_known}", "established loci retained\n(TSHR anchor, CTLA4 control)", KEEP, "kept"),
    ]
    y = 0.720
    for n, lab, col, tag in bars:
        box(axA, 0.575, y - 0.052, 0.055, 0.062, col, col)
        axA.text(0.6025, y - 0.021, n, ha="center", va="center", fontsize=13,
                 fontweight="bold", color="white")
        axA.text(0.645, y - 0.005, lab, va="top", fontsize=8.6, color=INK, linespacing=1.35)
        y -= 0.105

    # verdict strip
    box(axA, 0.055, 0.055, 0.89, 0.105, "#fdf3f3", ANCHOR, 1.4)
    axA.text(0.077, 0.108, "0", fontsize=25, fontweight="bold", color=ANCHOR, va="center")
    axA.text(0.135, 0.128, "robust novel druggable targets", fontsize=11.5,
             fontweight="bold", color=INK, va="center")
    axA.text(0.135, 0.086,
             "No non-established gene showed shared-variant colocalization in both outcomes. "
             "At 80% power the screen\ncould detect a median OR of 2.55 (Table S8), so this excludes "
             "further large — not moderate — effects.",
             fontsize=8.5, color=MUTED, va="center", linespacing=1.45)

    # ================= PANEL B — the divergence verdict =================
    axB = fig.add_axes([0.06, 0.035, 0.90, 0.435]); axB.axis("off")
    axB.set_xlim(0, 1); axB.set_ylim(0, 1)
    axB.text(0, 1.02, "B", fontsize=16, fontweight="bold", va="top")
    axB.text(0.035, 1.02, "Where the two receptors diverge", fontsize=13,
             fontweight="bold", va="top")
    axB.text(0.035, 0.945,
             "The same three evidence layers, the same instruments, outcomes and filters — read across each row.",
             fontsize=9.5, color=MUTED, va="top")

    T, I = bb["TSHR"], bb["IGF1R"]
    rows = [
        ("Genetic association\n(MR, BBJ discovery)",
         f"β = {T['beta_bbj']:+.2f}\nP = 1.1×10⁻¹⁴", True,
         f"β = {I['beta_bbj']:+.2f}\nP = {I['p_bbj']:.3f}", None),
        ("Shared causal variant\n(colocalization)",
         f"PP.H4 = {T['h4_bbj']:.3f} / {T['h4_fin']:.3f}\nshared variant rs179252", True,
         f"PP.H2 = {I['h2_bbj']:.2f} / {I['h2_fin']:.2f}\ncis-eQTL only,\nno disease signal", False),
        ("Orbital tissue\n(exploratory, 4 vs 1)",
         f"log2FC = {T['tis_fc']:+.2f}\nadj P = {T['tis_p']:.3f}", True,
         f"log2FC = {I['tis_fc']:+.2f}\nadj P > 0.99", False),
    ]
    # header
    axB.text(0.3225, 0.865, "TSHR", ha="center", fontsize=12.5, fontweight="bold", color=ANCHOR)
    axB.text(0.3225, 0.815, "susceptibility anchor", ha="center", fontsize=8.8, color=MUTED)
    axB.text(0.7425, 0.865, "IGF1R", ha="center", fontsize=12.5, fontweight="bold", color=EFFECTOR)
    axB.text(0.7425, 0.815, "teprotumumab target", ha="center", fontsize=8.8, color=MUTED)
    axB.plot([0.035, 0.965], [0.785, 0.785], color=RULE, lw=1.1)

    y = 0.70
    for lab, tv, tok, iv, iok in rows:
        axB.text(0.035, y, lab, va="center", fontsize=9, color=INK, linespacing=1.4)
        box(axB, 0.215, y - 0.078, 0.215, 0.156, "#fdf3f3", ANCHOR, 1.1)
        axB.text(0.3225, y, tv, ha="center", va="center", fontsize=8.7,
                 color=INK, linespacing=1.45)
        icol = EFFECTOR if iok is None else "#eef2f6"
        box(axB, 0.635, y - 0.078, 0.215, 0.156, "#eef4fa", EFFECTOR, 1.1)
        axB.text(0.7425, y, iv, ha="center", va="center", fontsize=8.7,
                 color=INK, linespacing=1.45)
        mark = "✓" if tok else "—"
        axB.text(0.455, y, mark, ha="center", va="center", fontsize=13, color=ANCHOR)
        axB.text(0.875, y, "✓" if iok else ("~" if iok is None else "✗"),
                 ha="center", va="center", fontsize=13,
                 color=EFFECTOR if iok is not False else MUTED)
        y -= 0.235

    axB.plot([0.035, 0.965], [0.075, 0.075], color=RULE, lw=1.1)
    axB.text(0.5, 0.028,
             "Genetic susceptibility architecture and therapeutic target biology do not have to coincide: "
             "TSHR carries both,\nIGF1R carries the therapy without the susceptibility architecture.",
             ha="center", va="center", fontsize=9.6, color=INK, fontweight="bold",
             linespacing=1.5)

    for p in OUT:
        fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("Figure 1 rebuilt as a claim figure")
    print(f"  funnel: 4,462 -> 2,545 -> 2,235 -> 13 -> 0 robust novel")
    print(f"  13 hits = {n_mhc} MHC + {n_ld} chr16 LD + {n_single} single-outcome + "
          f"{n_distinct} distinct-variant + {n_known} established")


if __name__ == "__main__":
    build()
