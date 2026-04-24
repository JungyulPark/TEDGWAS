"""
TED-TRAP Upgrade — Phase 4 Script 20
Network proximity z-scores for drug target sets vs disease gene set
=============================================================================
Purpose:
    Resolve MAJOR issue #15: Replace Venn-only description with
    quantitative network-proximity analysis using the Cheng et al. 2018
    framework (Nat Commun 9:2691).

    Network proximity measures how close a drug's target set is to a
    disease's gene set in a protein-protein interaction network. It
    corrects for network topology via a degree-preserving randomization.

Method:
    1. Build STRING human PPI network (confidence ≥ 0.7)
    2. Compute closest distance d_c(S, T) = mean over targets of min
       distance to any disease gene
    3. Random permutation: resample target set matched by degree
    4. Z-score = (observed d_c - mean_random) / sd_random

Two comparisons:
    A) IGF-1R target set (36) vs TED gene set (59)
    B) TSHR target set (33) vs TED gene set (59)

    Smaller d_c and more negative z → closer proximity → stronger target-disease link.

Output:
    TrackB_Network/results/20_network_proximity.csv

Requires:
    - STRING PPI file (protein.links.v12.0.txt.gz from string-db.org)
    - networkx

Run:
    python 20_network_proximity_zscore.py
"""

import networkx as nx
import pandas as pd
import numpy as np
import os
import gzip
import random
from collections import defaultdict
from datetime import datetime

PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
OUT_DIR = os.path.join(PROJECT_ROOT, "TrackB_Network", "results")
os.makedirs(OUT_DIR, exist_ok=True)


# ----- STRING PPI loading -----
STRING_FILE = os.path.join(
    PROJECT_ROOT, "TrackB_Network", "data",
    "9606.protein.links.v12.0.txt.gz"
)
STRING_INFO_FILE = os.path.join(
    PROJECT_ROOT, "TrackB_Network", "data",
    "9606.protein.info.v12.0.txt.gz"
)
# Download instructions:
#   https://string-db.org/cgi/download?species_text=Homo+sapiens
#   - 9606.protein.links.v12.0.txt.gz
#   - 9606.protein.info.v12.0.txt.gz


def load_string_network(confidence_threshold: int = 700):
    """Load STRING human PPI, filter to high-confidence edges, convert to symbols."""
    print(f"Loading STRING network from {STRING_FILE}...")
    if not os.path.exists(STRING_FILE) or not os.path.exists(STRING_INFO_FILE):
        print("❌ STRING files missing. Download from:")
        print("   https://string-db.org/cgi/download?species_text=Homo+sapiens")
        return None, None

    # Protein ID → symbol map
    info = pd.read_csv(STRING_INFO_FILE, sep="\t")
    sp_to_symbol = dict(zip(info["#string_protein_id"], info["preferred_name"]))

    # Build graph
    G = nx.Graph()
    n_added = 0
    with gzip.open(STRING_FILE, "rt") as f:
        header = f.readline()
        for line in f:
            p1, p2, score = line.strip().split()
            score = int(score)
            if score >= confidence_threshold:
                s1 = sp_to_symbol.get(p1)
                s2 = sp_to_symbol.get(p2)
                if s1 and s2:
                    G.add_edge(s1, s2, weight=score / 1000.0)
                    n_added += 1
    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges "
          f"(threshold={confidence_threshold})")
    return G, sp_to_symbol


def closest_distance(G: nx.Graph, source_set: set, target_set: set):
    """Mean of minimum shortest-path distances from source nodes to target set."""
    distances = []
    source_in_net = [s for s in source_set if s in G.nodes]
    target_in_net = [t for t in target_set if t in G.nodes]

    if not source_in_net or not target_in_net:
        return np.nan

    for s in source_in_net:
        min_d = np.inf
        for t in target_in_net:
            if s == t:
                min_d = 0
                break
            try:
                d = nx.shortest_path_length(G, s, t)
                if d < min_d:
                    min_d = d
            except nx.NetworkXNoPath:
                continue
        if min_d != np.inf:
            distances.append(min_d)
    return np.mean(distances) if distances else np.nan


def degree_binning(G: nx.Graph, bin_size: int = 100):
    """Bin nodes by degree, for degree-preserving sampling."""
    degrees = dict(G.degree())
    nodes_by_degree = defaultdict(list)
    for node, deg in degrees.items():
        nodes_by_degree[deg].append(node)

    # Merge bins to have at least bin_size nodes each
    sorted_degrees = sorted(nodes_by_degree.keys())
    bins = []
    current_bin = []
    for deg in sorted_degrees:
        current_bin.extend(nodes_by_degree[deg])
        if len(current_bin) >= bin_size:
            bins.append(current_bin)
            current_bin = []
    if current_bin:
        if bins:
            bins[-1].extend(current_bin)
        else:
            bins.append(current_bin)

    node_to_bin = {}
    for i, b in enumerate(bins):
        for n in b:
            node_to_bin[n] = i

    return bins, node_to_bin


def proximity_zscore(G: nx.Graph, source_set: set, target_set: set,
                     n_random: int = 1000, seed: int = 42):
    """Compute network proximity z-score with degree-preserving randomization."""
    d_observed = closest_distance(G, source_set, target_set)
    if np.isnan(d_observed):
        return {"observed": np.nan, "mean_random": np.nan,
                "sd_random": np.nan, "z": np.nan, "p_empirical": np.nan}

    print(f"  Observed d_c = {d_observed:.3f}")
    print(f"  Generating {n_random} degree-matched random sets...")

    # Degree binning for source set
    bins, node_to_bin = degree_binning(G, bin_size=100)
    rng = random.Random(seed)

    # Extract source nodes by bin
    source_in_net = [s for s in source_set if s in G.nodes]
    bin_assignment = [node_to_bin[s] for s in source_in_net]

    random_distances = []
    for i in range(n_random):
        sampled = set()
        for bin_idx in bin_assignment:
            candidate = rng.choice(bins[bin_idx])
            sampled.add(candidate)
        d_r = closest_distance(G, sampled, target_set)
        if not np.isnan(d_r):
            random_distances.append(d_r)
        if (i + 1) % 100 == 0:
            print(f"    {i+1}/{n_random} done")

    random_distances = np.array(random_distances)
    mean_r = np.mean(random_distances)
    sd_r = np.std(random_distances)
    z = (d_observed - mean_r) / sd_r if sd_r > 0 else np.nan
    p_emp = np.mean(random_distances <= d_observed)

    return {
        "observed": d_observed,
        "mean_random": mean_r,
        "sd_random": sd_r,
        "z": z,
        "p_empirical": p_emp,
    }


def main():
    print("=" * 70)
    print("Network Proximity Analysis (Cheng et al. 2018 framework)")
    print(f"Date: {datetime.now()}")
    print("=" * 70 + "\n")

    # Load network
    G, _ = load_string_network(confidence_threshold=700)
    if G is None:
        return

    # Define target sets (same as script 19)
    igf1r_pathway = {
        "IGF1R","IGF1","IGF2","IRS1","IRS2","SHC1","GRB2","SOS1",
        "PIK3CA","PIK3CB","PIK3R1","AKT1","AKT2","MTOR","RPS6KB1",
        "KRAS","BRAF","MAP2K1","MAP2K2","MAPK1","MAPK3",
        "JAK1","JAK2","STAT3","STAT5A",
        "ARRB1","ARRB2",
        "HAS1","HAS2","HAS3","PPARG","CEBPA",
        "IL6","CXCL8","TNF","IL1B"
    }
    tshr_pathway = {
        "TSHR","GNAS","GNAQ","GNA11",
        "ADCY1","ADCY2","ADCY3","ADCY4","ADCY5",
        "ADCY6","ADCY7","ADCY8","ADCY9",
        "PRKACA","PRKACB","PRKAR1A","CREB1","CREB3",
        "ARRB1","ARRB2","GRK2","GRK5",
        "FOXO1","FOXO3",
        "HAS1","HAS2","HAS3","IL6","TNF",
        "PPARG","CEBPA","FABP4","ADIPOQ"
    }

    # TED disease gene set
    v3_ted_path = os.path.join(OUT_DIR, "18_ted_gene_set_v3_final.csv")
    v1_ted_path = os.path.join(PROJECT_ROOT, "TED_disease_genes_with_sources.csv")
    if os.path.exists(v3_ted_path):
        ted_df = pd.read_csv(v3_ted_path)
        ted_set = set(ted_df["symbol"].str.upper())
        ted_source = "v3"
    else:
        ted_df = pd.read_csv(v1_ted_path)
        ted_set = set(ted_df["Gene"].str.upper())
        ted_source = "v1"
    print(f"TED gene set: {len(ted_set)} genes ({ted_source})\n")

    # Compute proximity for each drug target set
    results = []
    for name, set_src in [("IGF-1R_pathway", igf1r_pathway),
                          ("TSHR_pathway", tshr_pathway)]:
        print(f"\n--- {name} vs TED gene set ---")
        res = proximity_zscore(G, set_src, ted_set, n_random=500, seed=42)
        res["source_set"] = name
        res["source_size"] = len(set_src)
        res["target_size"] = len(ted_set)
        results.append(res)

    df = pd.DataFrame(results)
    # Reorder columns
    df = df[["source_set", "source_size", "target_size",
             "observed", "mean_random", "sd_random", "z", "p_empirical"]]

    print("\n=== Network Proximity Results ===")
    print(df.to_string(index=False))

    out = os.path.join(OUT_DIR, "20_network_proximity.csv")
    df.to_csv(out, index=False)
    print(f"\n✅ Saved: {out}")

    # Interpretation
    print("\n=== INTERPRETATION ===")
    print("Lower observed d_c and more negative z-score indicate CLOSER network")
    print("proximity between a drug's target set and the TED gene set.")
    print("Both IGF-1R and TSHR sets should show significant proximity (z < -2)")
    print("given they overlap with known TED genes; the comparison is about")
    print("which set is CLOSER to the disease neighborhood.")


if __name__ == "__main__":
    main()
