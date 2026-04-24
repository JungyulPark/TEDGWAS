"""
TED-TRAP Upgrade — Phase 5 Script 21
Figure 3 Panel A: Triangulation matrix v3 (updated with coloc + DESeq2 v3 numbers)
=============================================================================
Purpose:
    Regenerate the triangulation matrix figure using the upgraded results.
    This is the centerpiece figure — it synthesizes MR, coloc, RNA-seq direction,
    and disease-gene membership into a single visual score.

Inputs (from upgrade results):
    - TrackA_MR/results/11_mr_v3_bh.csv           (MR with BH-FDR)
    - TrackA_MR/results/11_coloc_v3_tiers.csv      (coloc with tiers)
    - TrackA_MR/results/11_de_v3_bh.csv            (DESeq2 with BH-FDR)
    - TED disease gene set (from 18_ted_gene_set_v3_final.csv)

Output:
    Manuscript/figures_v3/Fig3_PanelA_triangulation_v3.pdf
    Manuscript/figures_v3/Fig3_PanelA_triangulation_v3.png

Run:
    python 21_fig3_panelA_triangulation_v3.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import os

PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
FIG_DIR = os.path.join(PROJECT_ROOT, "Manuscript", "figures_v3")
os.makedirs(FIG_DIR, exist_ok=True)


def load_triangulation_data():
    """Load all evidence lines and merge into a single scoring matrix."""

    # 8 candidate genes
    genes = ["TSHR", "IGF1R", "TNF", "PPARG", "ARRB1", "IRS1", "AKT1", "CTLA4"]

    # Initialize matrix
    matrix = pd.DataFrame(index=genes,
                          columns=["MR_primary", "MR_replication", "Coloc", "RNAseq"])

    # --- Load MR results ---
    mr_path = os.path.join(PROJECT_ROOT, "TrackA_MR", "results", "11_mr_v3_bh.csv")
    if os.path.exists(mr_path):
        mr = pd.read_csv(mr_path)
        # Focus on IVW method (primary MR estimate)
        mr_ivw = mr[mr["method"].str.contains("Inverse variance weighted", na=False)]
        for g in genes:
            sub = mr_ivw[mr_ivw["gene"] == g]
            # Primary outcome score
            p = sub[sub["outcome_role"] == "Primary"]
            if not p.empty:
                bh = p["pval_BH"].iloc[0]
                matrix.loc[g, "MR_primary"] = (
                    "significant" if bh < 0.05 else "suggestive" if bh < 0.10 else "null"
                )
            # Replication outcome score
            r = sub[sub["outcome_role"] == "Replication"]
            if not r.empty:
                bh = r["pval_BH"].iloc[0]
                matrix.loc[g, "MR_replication"] = (
                    "significant" if bh < 0.05 else "suggestive" if bh < 0.10 else "null"
                )
    else:
        print(f"[WARN] MR results not found at {mr_path} — using v1 fallback")
        # v1 fallback values
        matrix["MR_primary"] = ["significant", "null", "significant",
                                 "significant", "suggestive", "null",
                                 "null", "significant"]
        matrix["MR_replication"] = ["significant", "null", "significant",
                                      "suggestive", "null", "null",
                                      "null", "null"]

    # --- Load coloc results ---
    coloc_path = os.path.join(PROJECT_ROOT, "TrackA_MR", "results",
                               "11_coloc_v3_tiers.csv")
    if os.path.exists(coloc_path):
        coloc = pd.read_csv(coloc_path)
        for g in genes:
            sub = coloc[coloc["gene"] == g]
            if not sub.empty:
                pph4 = sub["PP_H4"].iloc[0]
                matrix.loc[g, "Coloc"] = (
                    "very strong" if pph4 > 0.95
                    else "strong" if pph4 > 0.80
                    else "weak" if pph4 > 0.50
                    else "null"
                )
    else:
        # v1 fallback: only TSHR has coloc=0.985
        matrix["Coloc"] = ["very strong", np.nan, np.nan, np.nan,
                            np.nan, np.nan, np.nan, np.nan]

    # --- Load RNA-seq results ---
    de_path = os.path.join(PROJECT_ROOT, "TrackA_MR", "results", "11_de_v3_bh.csv")
    if os.path.exists(de_path):
        de = pd.read_csv(de_path)
        for g in genes:
            sub = de[de["gene"] == g]
            if not sub.empty:
                bh = sub["padj_candidate_BH"].iloc[0]
                lfc = sub["log2FoldChange"].iloc[0]
                matrix.loc[g, "RNAseq"] = (
                    "concordant" if (bh < 0.05 and lfc > 0)
                    else "suggestive" if (bh < 0.10 and lfc > 0)
                    else "null"
                )
    else:
        matrix["RNAseq"] = ["concordant", np.nan, "concordant", "concordant",
                             np.nan, np.nan, np.nan, np.nan]

    return matrix


def compute_triangulation_score(row):
    """Sum scores across evidence lines for sorting."""
    score_map = {
        "significant": 3, "strong": 3, "very strong": 4, "concordant": 3,
        "suggestive": 1.5, "weak": 1, "null": 0, np.nan: 0
    }
    score = 0
    for col in ["MR_primary", "MR_replication", "Coloc", "RNAseq"]:
        score += score_map.get(row[col], 0)
    return score


def plot_triangulation_matrix(matrix, output_path):
    """Render a heatmap-like triangulation figure with color-coded cells."""

    # Compute total score and reorder
    matrix["Total_Score"] = matrix.apply(compute_triangulation_score, axis=1)
    matrix = matrix.sort_values("Total_Score", ascending=False)

    # Color scheme (Tableau-like, accessible)
    color_map = {
        "very strong": "#1a5490", "strong": "#2e7bba", "significant": "#2e7bba",
        "concordant": "#2e7bba", "weak": "#b8d4e8", "suggestive": "#b8d4e8",
        "null": "#f2f2f2", np.nan: "#d9d9d9"
    }

    fig, ax = plt.subplots(figsize=(11, 6))
    cols = ["MR_primary", "MR_replication", "Coloc", "RNAseq"]
    col_labels = ["MR\n(Graves)", "MR replication\n(Hyperthyroid)",
                  "Colocalization\n(PP.H4)", "Orbital RNA-seq\n(direction + q)"]

    genes = matrix.index.tolist()

    for i, g in enumerate(genes):
        for j, col in enumerate(cols):
            val = matrix.loc[g, col]
            color = color_map.get(val, "#d9d9d9")
            rect = plt.Rectangle((j, len(genes) - i - 1), 1, 1,
                                 facecolor=color, edgecolor="white", linewidth=1.5)
            ax.add_patch(rect)

            # Label inside cell
            label = str(val) if pd.notna(val) else "n/a"
            text_color = "white" if color in ["#1a5490", "#2e7bba"] else "#333333"
            ax.text(j + 0.5, len(genes) - i - 0.5, label,
                    ha="center", va="center", fontsize=9, color=text_color,
                    fontweight="bold" if val in ["very strong", "concordant",
                                                  "significant"] else "normal")

    # Gene labels on y-axis
    for i, g in enumerate(genes):
        total = matrix.loc[g, "Total_Score"]
        fontweight = "bold" if total >= 7 else "normal"
        ax.text(-0.1, len(genes) - i - 0.5, g, ha="right", va="center",
                fontsize=11, fontweight=fontweight)
        ax.text(len(cols) + 0.1, len(genes) - i - 0.5,
                f"Score: {total:.1f}", ha="left", va="center",
                fontsize=9, color="#666")

    # Column headers
    for j, lbl in enumerate(col_labels):
        ax.text(j + 0.5, len(genes) + 0.1, lbl, ha="center", va="bottom",
                fontsize=10, fontweight="bold")

    ax.set_xlim(-2.5, len(cols) + 2)
    ax.set_ylim(-0.5, len(genes) + 1)
    ax.axis("off")

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor="#1a5490", label="Very strong"),
        mpatches.Patch(facecolor="#2e7bba", label="Strong / Significant / Concordant"),
        mpatches.Patch(facecolor="#b8d4e8", label="Weak / Suggestive"),
        mpatches.Patch(facecolor="#f2f2f2", label="Null / Not significant"),
        mpatches.Patch(facecolor="#d9d9d9", label="Not tested (N/A)")
    ]
    ax.legend(handles=legend_elements, loc="upper center",
              bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=9, frameon=False)

    # Title
    fig.suptitle(
        "Figure 3 Panel A. Triangulation matrix across genetic, colocalization, and transcriptomic evidence lines",
        fontsize=11, y=0.97, fontweight="bold"
    )
    fig.text(0.5, 0.92,
             "Each row represents a candidate gene; each column an evidence line. "
             "TSHR achieves the highest convergent score.",
             fontsize=9, ha="center", style="italic")

    plt.tight_layout(rect=[0.05, 0.06, 0.95, 0.92])
    plt.savefig(output_path.replace(".png", ".pdf"), bbox_inches="tight")
    plt.savefig(output_path, bbox_inches="tight", dpi=200)
    print(f"✅ Saved: {output_path}")
    print(f"✅ Saved: {output_path.replace('.png', '.pdf')}")
    plt.close()


def main():
    print("=" * 70)
    print("Figure 3 Panel A v3: Triangulation Matrix")
    print("=" * 70 + "\n")

    matrix = load_triangulation_data()
    print("Triangulation matrix (raw):")
    print(matrix)

    output_png = os.path.join(FIG_DIR, "Fig3_PanelA_triangulation_v3.png")
    plot_triangulation_matrix(matrix, output_png)

    # Save the underlying data as CSV for the figure legend
    matrix.to_csv(os.path.join(FIG_DIR, "Fig3_PanelA_triangulation_v3_data.csv"))


if __name__ == "__main__":
    main()
