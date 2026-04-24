"""
TED-TRAP Upgrade — Phase 3 Script 14
Open Targets Platform GraphQL query for Graves ophthalmopathy / TED
=============================================================================
Purpose:
    Resolve MAJOR issue #10: Manual TED gene curation lacks reproducibility.
    This script queries Open Targets Platform (Ochoa et al. 2023, NAR 51:D1353)
    for disease-associated targets using the EFO identifier for
    "thyroid-associated ophthalmopathy" (EFO_0005606 or MONDO_0007303).

Output:
    TrackB_Network/results/14_opentargets_TED_genes.csv

Run:
    python 14_opentargets_query.py
"""

import requests
import pandas as pd
import json
import os
from datetime import datetime

# Working directory
PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
os.makedirs(os.path.join(PROJECT_ROOT, "TrackB_Network", "results"), exist_ok=True)
os.makedirs(os.path.join(PROJECT_ROOT, "TrackB_Network", "logs"), exist_ok=True)

# Open Targets GraphQL endpoint
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Target disease EFO IDs to try (TED synonyms)
# Open Targets uses EFO; Graves ophthalmopathy is sometimes MONDO
TED_EFO_IDS = [
    "EFO_0005606",  # Graves ophthalmopathy
    "EFO_0004237",  # Thyroid eye disease
    "MONDO_0007303",  # Graves disease (broader)
    "Orphanet_466551",  # Possibly TED
]

# GraphQL query
QUERY = """
query AssociatedTargets($efoId: String!, $pageSize: Int!, $cursor: String) {
  disease(efoId: $efoId) {
    id
    name
    description
    associatedTargets(page: {size: $pageSize, index: 0}) {
      count
      rows {
        target {
          approvedSymbol
          id
          biotype
          approvedName
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"""


def query_disease(efo_id: str, page_size: int = 500) -> dict:
    """Query Open Targets for targets associated with a disease."""
    payload = {
        "query": QUERY,
        "variables": {"efoId": efo_id, "pageSize": page_size},
    }
    resp = requests.post(OT_URL, json=payload, timeout=60)
    resp.raise_for_status()
    return resp.json()


def parse_targets(data: dict, efo_id: str) -> pd.DataFrame:
    """Extract target info into a dataframe."""
    if "errors" in data:
        print(f"  GraphQL errors: {data['errors']}")
        return pd.DataFrame()
    disease = data.get("data", {}).get("disease")
    if not disease:
        print(f"  No disease found for EFO: {efo_id}")
        return pd.DataFrame()
    print(f"  Disease: {disease['name']} ({disease['id']})")
    print(f"  Total associated targets: {disease['associatedTargets']['count']}")

    rows = []
    for row in disease["associatedTargets"]["rows"]:
        t = row["target"]
        record = {
            "efo_id": efo_id,
            "disease_name": disease["name"],
            "target_id": t["id"],
            "symbol": t.get("approvedSymbol", ""),
            "approved_name": t.get("approvedName", ""),
            "biotype": t.get("biotype", ""),
            "overall_score": row["score"],
        }
        # Datatype breakdown
        for dt_score in row.get("datatypeScores", []):
            record[f"score_{dt_score['id']}"] = dt_score["score"]
        rows.append(record)
    return pd.DataFrame(rows)


def main():
    print("=" * 70)
    print("Open Targets Platform — TED target query")
    print(f"Date: {datetime.now()}")
    print("=" * 70 + "\n")

    all_results = []
    for efo in TED_EFO_IDS:
        print(f"\nQuerying EFO: {efo}")
        try:
            data = query_disease(efo)
            df = parse_targets(data, efo)
            if not df.empty:
                all_results.append(df)
        except Exception as e:
            print(f"  [ERROR] {e}")

    if not all_results:
        print("\n❌ No results retrieved from any EFO ID.")
        print("   Check API access and EFO ID validity.")
        return

    # Combine
    combined = pd.concat(all_results, ignore_index=True)
    print(f"\nTotal rows across all EFOs: {len(combined)}")

    # Keep unique (gene, efo_id) pairs
    unique_per_efo = combined.drop_duplicates(subset=["efo_id", "symbol"])
    print(f"Unique (target × disease): {len(unique_per_efo)}")

    # Aggregate: best score per gene across EFOs
    best_per_gene = (
        combined.groupby("symbol")
        .agg(
            best_efo=("efo_id", "first"),
            best_disease=("disease_name", "first"),
            max_overall_score=("overall_score", "max"),
            n_efos_hit=("efo_id", "nunique"),
        )
        .reset_index()
        .sort_values("max_overall_score", ascending=False)
    )

    # Filter to plausible disease genes (score > 0.1 is OT's common threshold)
    plausible = best_per_gene[best_per_gene["max_overall_score"] > 0.1].copy()
    print(f"\nGenes with association score > 0.1: {len(plausible)}")
    print(f"Top 20:")
    print(plausible.head(20).to_string(index=False))

    # --- Save ---
    out_dir = os.path.join(PROJECT_ROOT, "TrackB_Network", "results")
    combined.to_csv(os.path.join(out_dir, "14_opentargets_raw.csv"), index=False)
    plausible.to_csv(
        os.path.join(out_dir, "14_opentargets_TED_genes.csv"), index=False
    )

    # --- Sanity check with known TED genes ---
    print("\n=== Concordance with v1 manually curated genes ===")
    v1_known = [
        "TSHR", "IGF1R", "IGF1", "HLA-DRB1", "CTLA4", "CD40", "PTPN22",
        "FOXP3", "PPARG", "HAS2", "IL6", "TNF", "TGFB1",
    ]
    v1_found = plausible[plausible["symbol"].isin(v1_known)]
    missing = [g for g in v1_known if g not in plausible["symbol"].values]
    print(f"v1 genes confirmed by OT: {len(v1_found)} / {len(v1_known)}")
    print(f"  Found: {v1_found['symbol'].tolist()}")
    print(f"  Missing (possibly manual-only): {missing}")

    print(f"\n✅ Saved: {os.path.join(out_dir, '14_opentargets_TED_genes.csv')}")


if __name__ == "__main__":
    main()
