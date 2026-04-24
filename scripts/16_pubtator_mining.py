"""
TED-TRAP Upgrade — Phase 3 Script 16
PubTator Central mining: TED disease-gene co-occurrence from literature
=============================================================================
Purpose:
    Query NCBI PubTator3 for gene-disease co-occurrences in articles indexed
    for thyroid eye disease / Graves ophthalmopathy. This provides a
    literature-derived gene set that is reproducible and independent of
    manual curation.

Reference:
    Wei CH, Allot A, Leaman R, Lu Z (2024). PubTator 3.0: an AI-powered
    Literature Resource for Unlocking Biomedical Knowledge.
    Nucleic Acids Res doi:10.1093/nar/gkae235.

API:
    https://www.ncbi.nlm.nih.gov/research/pubtator3/api/v1/publications

Strategy:
    1. Search PubMed for TED-related articles via eutils
    2. For each PMID, fetch entity annotations from PubTator
    3. Count gene mentions across all TED-indexed papers
    4. Filter by co-occurrence count threshold

Output:
    TrackB_Network/results/16_pubtator_TED_genes.csv

Run:
    python 16_pubtator_mining.py
"""

import requests
import pandas as pd
import os
import sys
import time
from datetime import datetime
from collections import Counter

PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
OUT_DIR = os.path.join(PROJECT_ROOT, "TrackB_Network", "results")
os.makedirs(OUT_DIR, exist_ok=True)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUBTATOR_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3/api/v1"

# Search terms for TED literature
TED_QUERIES = [
    '"thyroid eye disease"[Title/Abstract]',
    '"Graves ophthalmopathy"[Title/Abstract]',
    '"Graves\' orbitopathy"[Title/Abstract]',
    '"thyroid-associated ophthalmopathy"[Title/Abstract]',
    '"thyroid-associated orbitopathy"[Title/Abstract]',
]


def search_pubmed(query: str, retmax: int = 1000) -> list:
    """Search PubMed and return list of PMIDs."""
    url = f"{EUTILS_BASE}/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": retmax,
        "retmode": "json",
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    return data.get("esearchresult", {}).get("idlist", [])


def fetch_pubtator_annotations(pmids_batch: list) -> dict:
    """Fetch PubTator annotations for a batch of PMIDs."""
    # PubTator3 accepts batches of IDs
    url = f"{PUBTATOR_BASE}/publications/export/biocjson"
    params = {"pmids": ",".join(pmids_batch[:100])}  # 100 per request
    r = requests.get(url, params=params, timeout=60)
    if r.status_code != 200:
        return {}
    # PubTator returns BioC JSON
    return r.json() if r.text else {}


def extract_gene_entities(bioc_json: dict) -> list:
    """Extract (PMID, gene_symbol) pairs from BioC JSON."""
    pairs = []
    if not bioc_json:
        return pairs

    # BioC JSON structure: PubTator3Document with documents list
    documents = bioc_json.get("documents", [])
    for doc in documents:
        pmid = doc.get("id", "")
        for passage in doc.get("passages", []):
            for annot in passage.get("annotations", []):
                if annot.get("infons", {}).get("type") == "Gene":
                    gene_id = annot.get("infons", {}).get("identifier", "")
                    text = annot.get("text", "")
                    if gene_id and text:
                        pairs.append((pmid, gene_id, text))
    return pairs


def main():
    print("=" * 70)
    print("PubTator Central — TED literature gene mining")
    print(f"Date: {datetime.now()}")
    print("=" * 70 + "\n")

    # Step 1: Collect PMIDs
    all_pmids = set()
    for q in TED_QUERIES:
        print(f"Searching PubMed: {q}")
        try:
            ids = search_pubmed(q, retmax=2000)
            print(f"  Found {len(ids)} PMIDs")
            all_pmids.update(ids)
            time.sleep(0.5)
        except Exception as e:
            print(f"  [ERROR] {e}")

    print(f"\nTotal unique TED PMIDs: {len(all_pmids)}")
    if len(all_pmids) == 0:
        print("❌ No PMIDs retrieved. Check network connection.")
        return

    # Step 2: Fetch annotations in batches
    pmid_list = sorted(all_pmids)
    all_pairs = []
    batch_size = 50
    for i in range(0, len(pmid_list), batch_size):
        batch = pmid_list[i : i + batch_size]
        print(f"Fetching PubTator annotations for batch {i // batch_size + 1} "
              f"({len(batch)} PMIDs)...")
        try:
            data = fetch_pubtator_annotations(batch)
            pairs = extract_gene_entities(data)
            all_pairs.extend(pairs)
        except Exception as e:
            print(f"  [ERROR] {e}")
        time.sleep(1.0)   # rate-limit

    print(f"\nTotal (PMID, gene) co-occurrences: {len(all_pairs)}")

    if not all_pairs:
        print("❌ No gene entities extracted.")
        print("   Possible causes: PubTator API limit; JSON format change.")
        return

    # Step 3: Aggregate by gene
    df_pairs = pd.DataFrame(all_pairs, columns=["pmid", "gene_id", "symbol"])
    df_pairs["symbol"] = df_pairs["symbol"].str.upper()

    gene_counts = (
        df_pairs.groupby("symbol")
        .agg(
            n_papers=("pmid", "nunique"),
            gene_id=("gene_id", "first"),
            example_pmids=("pmid", lambda x: ",".join(list(set(x))[:5])),
        )
        .reset_index()
        .sort_values("n_papers", ascending=False)
    )

    # Step 4: Threshold: genes appearing in ≥3 TED papers
    threshold = 3
    filtered = gene_counts[gene_counts["n_papers"] >= threshold].copy()
    print(f"\nGenes with ≥{threshold} TED paper co-occurrences: {len(filtered)}")
    print(f"Top 30:")
    print(filtered.head(30).to_string(index=False))

    # --- Save ---
    out_full = os.path.join(OUT_DIR, "16_pubtator_all_mentions.csv")
    out_filt = os.path.join(OUT_DIR, "16_pubtator_TED_genes.csv")
    gene_counts.to_csv(out_full, index=False)
    filtered.to_csv(out_filt, index=False)

    print(f"\n✅ Saved:")
    print(f"   {out_full}   (all mentions)")
    print(f"   {out_filt}   (filtered, ≥{threshold} papers)")


if __name__ == "__main__":
    main()
