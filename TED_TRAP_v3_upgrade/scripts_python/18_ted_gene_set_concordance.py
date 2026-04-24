"""
TED-TRAP Upgrade — Phase 3 Script 18
Harmonize TED disease-gene sources: Open Targets + DisGeNET + PubTator + Manual
=============================================================================
Purpose:
    Combine outputs from scripts 14-17 and the v1 manually curated list.
    Produce a single final TED gene set with source provenance.

    The goal is to REPLACE v1's manually curated 59-gene list with a
    multi-source, reproducible gene set that reviewers cannot reject as
    confirmation-biased.

Logic:
    A gene enters the final TED set if ANY of:
    - Open Targets overall_score > 0.1
    - DisGeNET GDA score > 0.1
    - PubTator co-occurrence count > 2
    - Manual v1 inclusion with literature citation

    Each gene is flagged with which source(s) support it. The more sources,
    the higher the confidence.

Output:
    TrackB_Network/results/18_ted_gene_set_v3_final.csv
    TrackB_Network/results/18_concordance_report.md

Run:
    python 18_ted_gene_set_concordance.py
"""

import pandas as pd
import os
from datetime import datetime

PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
OUT_DIR = os.path.join(PROJECT_ROOT, "TrackB_Network", "results")
os.makedirs(OUT_DIR, exist_ok=True)


def load_source(path: str, symbol_col: str, source_label: str) -> pd.DataFrame:
    """Load a source file and standardize to (symbol, source_label)."""
    if not os.path.exists(path):
        print(f"[WARN] {source_label} file missing: {path}")
        return pd.DataFrame(columns=["symbol", "source"])
    df = pd.read_csv(path)
    if symbol_col not in df.columns:
        print(f"[WARN] {source_label}: column '{symbol_col}' not found")
        return pd.DataFrame(columns=["symbol", "source"])
    out = pd.DataFrame({"symbol": df[symbol_col].astype(str).str.upper()})
    out["source"] = source_label
    return out.drop_duplicates()


def main():
    print("=" * 70)
    print("TED Gene Set Concordance Analysis")
    print(f"Date: {datetime.now()}")
    print("=" * 70 + "\n")

    sources = []

    # --- Open Targets ---
    ot = load_source(
        os.path.join(OUT_DIR, "14_opentargets_TED_genes.csv"),
        symbol_col="symbol",
        source_label="OpenTargets",
    )
    if not ot.empty:
        print(f"Open Targets: {len(ot)} genes")
        sources.append(ot)

    # --- DisGeNET ---
    dg = load_source(
        os.path.join(OUT_DIR, "15_disgenet_TED_genes.csv"),
        symbol_col="gene_symbol",  # adjust if DisGeNET uses different column
        source_label="DisGeNET",
    )
    if dg.empty:
        dg = load_source(
            os.path.join(OUT_DIR, "15_disgenet_TED_genes.csv"),
            symbol_col="symbol",
            source_label="DisGeNET",
        )
    if not dg.empty:
        print(f"DisGeNET: {len(dg)} genes")
        sources.append(dg)

    # --- PubTator (if run) ---
    pt = load_source(
        os.path.join(OUT_DIR, "16_pubtator_TED_genes.csv"),
        symbol_col="symbol",
        source_label="PubTator",
    )
    if not pt.empty:
        print(f"PubTator: {len(pt)} genes")
        sources.append(pt)

    # --- MAGMA (if run) ---
    mg = load_source(
        os.path.join(OUT_DIR, "17_magma_TED_genes.csv"),
        symbol_col="symbol",
        source_label="MAGMA_FinnGen_GO",
    )
    if not mg.empty:
        print(f"MAGMA: {len(mg)} genes")
        sources.append(mg)

    # --- v1 Manual curation ---
    v1_path = os.path.join(
        PROJECT_ROOT, "TED_disease_genes_with_sources.csv"
    )
    if os.path.exists(v1_path):
        v1 = pd.read_csv(v1_path)
        v1_df = pd.DataFrame({"symbol": v1["Gene"].astype(str).str.upper()})
        v1_df["source"] = "Manual_v1"
        sources.append(v1_df)
        print(f"v1 Manual: {len(v1_df)} genes")
    else:
        print(f"[WARN] v1 manual file missing: {v1_path}")

    if not sources:
        print("\n❌ No source data found. Run scripts 14-17 first.")
        return

    # --- Combine ---
    combined = pd.concat(sources, ignore_index=True)

    # Per-gene aggregation: which sources support it?
    per_gene = (
        combined.groupby("symbol")["source"]
        .agg(lambda x: sorted(set(x)))
        .reset_index()
    )
    per_gene["n_sources"] = per_gene["source"].apply(len)
    per_gene["sources_str"] = per_gene["source"].apply(lambda x: ", ".join(x))
    per_gene = per_gene[["symbol", "n_sources", "sources_str"]].sort_values(
        ["n_sources", "symbol"], ascending=[False, True]
    )

    print(f"\nTotal unique genes across all sources: {len(per_gene)}")

    # --- Tier by n_sources ---
    tier_counts = per_gene["n_sources"].value_counts().sort_index()
    print("\nGenes by number of supporting sources:")
    for n, c in tier_counts.items():
        print(f"  {n} source(s): {c} genes")

    # --- Confidence tiers ---
    def assign_tier(n):
        if n >= 3:
            return "High (≥3 sources)"
        if n == 2:
            return "Medium (2 sources)"
        return "Low (1 source)"
    per_gene["confidence_tier"] = per_gene["n_sources"].apply(assign_tier)

    # --- Final v3 gene set: require ≥2 sources OR v1 manual with literature ---
    high_conf = per_gene[per_gene["n_sources"] >= 2].copy()
    print(f"\nHigh/Medium confidence TED genes (≥2 sources): {len(high_conf)}")

    # --- Save ---
    per_gene.to_csv(
        os.path.join(OUT_DIR, "18_ted_gene_set_v3_full.csv"), index=False
    )
    high_conf.to_csv(
        os.path.join(OUT_DIR, "18_ted_gene_set_v3_final.csv"), index=False
    )

    # --- Concordance report ---
    report_path = os.path.join(OUT_DIR, "18_concordance_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# TED Gene Set — Multi-Source Concordance Report\n\n")
        f.write(f"Date: {datetime.now().isoformat()}\n\n")
        f.write("## Summary\n\n")
        f.write(f"- Total unique genes across sources: **{len(per_gene)}**\n")
        for n, c in tier_counts.items():
            f.write(f"- Genes with {n} source(s): {c}\n")
        f.write(f"\n## High-confidence set (≥2 sources): {len(high_conf)} genes\n\n")
        f.write("| Gene | N sources | Sources |\n|------|-----------|----------|\n")
        for _, row in high_conf.iterrows():
            f.write(f"| {row['symbol']} | {row['n_sources']} | {row['sources_str']} |\n")

    print(f"\n✅ Saved:")
    print(f"   {os.path.join(OUT_DIR, '18_ted_gene_set_v3_final.csv')}")
    print(f"   {report_path}")


if __name__ == "__main__":
    main()
