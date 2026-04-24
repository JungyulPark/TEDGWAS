"""
TED-TRAP Upgrade — Phase 3 Script 15
DisGeNET REST API query for TED/Graves ophthalmopathy
=============================================================================
Purpose:
    Query DisGeNET v7+ for gene-disease associations using UMLS CUI codes
    for TED. This replaces the v1 manual curation that returned zero TED
    entries (likely because the v1 query used the free-text search rather
    than the UMLS CUI API endpoint).

UMLS CUI codes:
    C0018818  — Graves disease (broad)
    C0042373  — Graves ophthalmopathy (endocrine exophthalmos)
    C0013936  — Endocrine exophthalmos

Output:
    TrackB_Network/results/15_disgenet_TED_genes.csv

Authentication:
    DisGeNET requires a free registration at https://www.disgenet.org/signup
    After registration, set env var DISGENET_API_KEY or paste interactively.

Run:
    python 15_disgenet_rest_api.py
"""

import requests
import pandas as pd
import os
import sys
from datetime import datetime

PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
os.makedirs(os.path.join(PROJECT_ROOT, "TrackB_Network", "results"), exist_ok=True)

DISGENET_BASE = "https://www.disgenet.org/api"
TED_CUIS = {
    "C0018818": "Graves disease",
    "C0042373": "Graves ophthalmopathy",
    "C0013936": "Endocrine exophthalmos",
}


def get_api_key():
    key = os.environ.get("DISGENET_API_KEY", "").strip()
    if not key:
        print("DisGeNET API key not found in environment.")
        print("Get one at: https://www.disgenet.com/")
        print("Then either:")
        print("  a) Set env var: setx DISGENET_API_KEY 'your_key' (Windows)")
        print("  b) Enter interactively below")
        key = input("\nEnter API key (or press ENTER to skip): ").strip()
    return key


def fetch_cui_associations(cui: str, api_key: str) -> pd.DataFrame:
    """Fetch gene-disease associations for a CUI."""
    url = f"{DISGENET_BASE}/gda/disease/{cui}"
    headers = {"Authorization": f"Bearer {api_key}"} if api_key else {}

    r = requests.get(url, headers=headers, timeout=60)
    if r.status_code == 401:
        print(f"  [401] Authentication failed. Check API key.")
        return pd.DataFrame()
    if r.status_code == 404:
        print(f"  [404] CUI {cui} not found.")
        return pd.DataFrame()
    r.raise_for_status()

    data = r.json()
    if isinstance(data, dict) and "error" in data:
        print(f"  API error: {data['error']}")
        return pd.DataFrame()

    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)


def main():
    print("=" * 70)
    print("DisGeNET REST API — TED gene-disease associations")
    print(f"Date: {datetime.now()}")
    print("=" * 70 + "\n")

    api_key = get_api_key()
    if not api_key:
        print("❌ No API key provided. Cannot proceed.")
        sys.exit(1)

    all_rows = []
    for cui, label in TED_CUIS.items():
        print(f"\nCUI {cui} ({label}):")
        df = fetch_cui_associations(cui, api_key)
        if df.empty:
            print("  No associations retrieved.")
            continue
        df["query_cui"] = cui
        df["query_label"] = label
        print(f"  Retrieved {len(df)} associations")
        all_rows.append(df)

    if not all_rows:
        print("\n❌ No data retrieved.")
        return

    combined = pd.concat(all_rows, ignore_index=True)
    print(f"\nTotal rows across all CUIs: {len(combined)}")

    # DisGeNET columns typically include: gene_symbol, geneid, score,
    # ei (evidence index), yearInitial, yearFinal, pmid (when available)
    # Print columns for user clarity
    print(f"\nColumns: {list(combined.columns)[:15]}")

    # Score filter: DisGeNET GDA score > 0.1 is common threshold
    score_col = [c for c in combined.columns if "score" in c.lower() and c != "ei_score"]
    score_col = score_col[0] if score_col else None

    if score_col:
        plausible = combined[combined[score_col] > 0.05].copy()
        print(f"\nAssociations with score > 0.05: {len(plausible)}")
    else:
        plausible = combined
        print("(No 'score' column detected; saving all)")

    # Unique genes
    gene_col = None
    for c in ["gene_symbol", "symbol", "geneSymbol"]:
        if c in plausible.columns:
            gene_col = c
            break

    if gene_col:
        unique_genes = plausible[gene_col].dropna().unique()
        print(f"Unique gene symbols: {len(unique_genes)}")
        print(f"Top 20: {sorted(unique_genes)[:20]}")

    # --- Save ---
    out_dir = os.path.join(PROJECT_ROOT, "TrackB_Network", "results")
    combined.to_csv(os.path.join(out_dir, "15_disgenet_raw.csv"), index=False)
    plausible.to_csv(
        os.path.join(out_dir, "15_disgenet_TED_genes.csv"), index=False
    )
    print(f"\n✅ Saved: {os.path.join(out_dir, '15_disgenet_TED_genes.csv')}")


if __name__ == "__main__":
    main()
