"""
TED-TRAP Upgrade — Phase 4 Script 19
Hypergeometric enrichment tests for pathway intersections
=============================================================================
Purpose:
    Resolve MAJOR issue #15: "Venn-based differential coverage needs statistical
    enrichment testing."

    For each intersection (e.g., IGF-1R pathway ∩ TED gene set), compute:
      - Observed overlap (k)
      - Expected overlap under null hypothesis
      - Hypergeometric p-value
      - Fold-enrichment
      - BH-corrected q-value (if multiple tests)

    Uses scipy.stats.hypergeom (test = P(X >= k)).

    This upgrades the Venn diagram from descriptive to inferential.

Key tests:
    1. IGF-1R pathway ∩ TED gene set
    2. TSHR pathway ∩ TED gene set
    3. (IGF-1R ∩ TSHR) ∩ TED gene set  (three-way)
    4. IGF-1R expanded (359) ∩ Hearing loss gene set (855)
    5. IGF-1R expanded (359) ∩ KEGG insulin signaling (137)

Output:
    TrackB_Network/results/19_hypergeometric_enrichment.csv

Run:
    python 19_hypergeometric_enrichment.py
"""

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import os
from datetime import datetime

PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
OUT_DIR = os.path.join(PROJECT_ROOT, "TrackB_Network", "results")
os.makedirs(OUT_DIR, exist_ok=True)

# Human protein-coding gene universe (Ensembl GRCh38 ≈ 19,962; use 20000 round)
GENE_UNIVERSE_SIZE = 20000


def hypergeom_test(set_A: set, set_B: set, universe_size: int, test_name: str):
    """Compute hypergeometric enrichment p-value for set A ∩ set B."""
    M = universe_size            # total pop
    n = len(set_A)               # 'successes' in pop
    N = len(set_B)               # sample size
    k = len(set_A & set_B)       # observed successes in sample
    expected = n * N / M
    fold_enrichment = k / expected if expected > 0 else np.inf

    # P(X >= k) — enrichment test
    p_value = hypergeom.sf(k - 1, M, n, N)

    return {
        "test": test_name,
        "set_A_size": n,
        "set_B_size": N,
        "universe": M,
        "observed_overlap": k,
        "expected_overlap": round(expected, 2),
        "fold_enrichment": round(fold_enrichment, 2),
        "hypergeom_p": p_value,
        "overlap_genes": sorted(set_A & set_B),
    }


def main():
    print("=" * 70)
    print("Hypergeometric Enrichment Tests")
    print(f"Date: {datetime.now()}")
    print("=" * 70 + "\n")

    # ===== Define gene sets =====

    # IGF-1R pathway (n=36 from v1)
    igf1r_pathway = {
        "IGF1R","IGF1","IGF2","IRS1","IRS2","SHC1","GRB2","SOS1",
        "PIK3CA","PIK3CB","PIK3R1","AKT1","AKT2","MTOR","RPS6KB1",
        "KRAS","BRAF","MAP2K1","MAP2K2","MAPK1","MAPK3",
        "JAK1","JAK2","STAT3","STAT5A",
        "ARRB1","ARRB2",
        "HAS1","HAS2","HAS3","PPARG","CEBPA",
        "IL6","CXCL8","TNF","IL1B"
    }

    # TSHR pathway (n=33)
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
    # v1 used 59 manually curated; v3 will use Phase 3 output
    # For now, load whichever is available
    v3_ted_path = os.path.join(OUT_DIR, "18_ted_gene_set_v3_final.csv")
    # Search v1 file in multiple likely locations
    v1_candidates = [
        os.path.join(PROJECT_ROOT, "TED_disease_genes_with_sources.csv"),
        os.path.join(PROJECT_ROOT, "data", "TED_disease_genes_with_sources.csv"),
        os.path.join(PROJECT_ROOT, "TrackB_Network", "data", "TED_disease_genes_with_sources.csv"),
    ]
    v1_ted_path = next((p for p in v1_candidates if os.path.exists(p)), v1_candidates[0])

    if os.path.exists(v3_ted_path):
        ted_df = pd.read_csv(v3_ted_path)
        ted_genes = set(ted_df["symbol"].str.upper())
        ted_source = "v3 multi-source"
    elif os.path.exists(v1_ted_path):
        ted_df = pd.read_csv(v1_ted_path)
        ted_genes = set(ted_df["Gene"].str.upper())
        ted_source = "v1 manual"
    else:
        print("❌ TED gene set not found.")
        return

    print(f"TED gene set source: {ted_source} ({len(ted_genes)} genes)\n")

    # IGF-1R expanded downstream (359 genes from Excel off-target sheet)
    # These would be loaded from file; for now simulate with KEGG PI3K-Akt
    # Load from Phase 1 output if available
    # TODO: replace with actual IGF-1R 359-gene list load

    # Cochlear / hearing loss gene set (855 from DisGeNET)
    # Load from file; for now, use placeholder size
    n_cochlear = 855

    # Insulin signaling gene set (137 from KEGG hsa04910)
    n_insulin = 137

    # IGF-1R expanded ∩ cochlear = 38 (observed)
    # IGF-1R expanded ∩ insulin = 51 (observed)

    # ===== Run tests =====
    tests = []

    # Test 1: IGF-1R pathway ∩ TED
    tests.append(hypergeom_test(
        igf1r_pathway, ted_genes, GENE_UNIVERSE_SIZE,
        "IGF-1R pathway (n=36) ∩ TED genes"
    ))

    # Test 2: TSHR pathway ∩ TED
    tests.append(hypergeom_test(
        tshr_pathway, ted_genes, GENE_UNIVERSE_SIZE,
        "TSHR pathway (n=33) ∩ TED genes"
    ))

    # Test 3: three-way intersection
    shared = igf1r_pathway & tshr_pathway
    tests.append(hypergeom_test(
        shared, ted_genes, GENE_UNIVERSE_SIZE,
        "(IGF-1R ∩ TSHR) pathways ∩ TED genes"
    ))

    # Test 4 & 5: off-target (using observed counts and set sizes;
    # no actual gene list needed for hypergeometric with known overlap)
    # Here we compute directly from counts:
    def hypergeom_direct(n_A: int, n_B: int, k_obs: int, M: int, name: str):
        expected = n_A * n_B / M
        fold = k_obs / expected if expected > 0 else np.inf
        p = hypergeom.sf(k_obs - 1, M, n_A, n_B)
        return {
            "test": name,
            "set_A_size": n_A,
            "set_B_size": n_B,
            "universe": M,
            "observed_overlap": k_obs,
            "expected_overlap": round(expected, 2),
            "fold_enrichment": round(fold, 2),
            "hypergeom_p": p,
            "overlap_genes": "(not listed here — see Supp Table S5)",
        }

    tests.append(hypergeom_direct(
        359, n_cochlear, 38, GENE_UNIVERSE_SIZE,
        "IGF-1R downstream (n=359) ∩ Hearing loss (n=855)"
    ))
    tests.append(hypergeom_direct(
        359, n_insulin, 51, GENE_UNIVERSE_SIZE,
        "IGF-1R downstream (n=359) ∩ KEGG hsa04910 Insulin (n=137)"
    ))

    # ===== BH correction =====
    pvals = [t["hypergeom_p"] for t in tests]
    _, q_bh, _, _ = multipletests(pvals, method="fdr_bh")
    for t, q in zip(tests, q_bh):
        t["BH_q"] = q

    # ===== Output =====
    df = pd.DataFrame(tests)
    df["overlap_genes"] = df["overlap_genes"].apply(
        lambda x: ", ".join(x) if isinstance(x, list) else x
    )

    print("=== Enrichment Results ===\n")
    display_cols = [
        "test", "set_A_size", "set_B_size",
        "observed_overlap", "expected_overlap",
        "fold_enrichment", "hypergeom_p", "BH_q"
    ]
    # Nicely formatted print
    for t in tests:
        print(f"▸ {t['test']}")
        print(f"    Observed: {t['observed_overlap']}   Expected: {t['expected_overlap']}")
        print(f"    Fold enrichment: {t['fold_enrichment']}×   p = {t['hypergeom_p']:.3e}   BH q = {t['BH_q']:.3e}")
        if isinstance(t.get('overlap_genes'), list) and len(t['overlap_genes']) <= 15:
            print(f"    Overlap genes: {', '.join(t['overlap_genes'])}")
        print()

    out = os.path.join(OUT_DIR, "19_hypergeometric_enrichment.csv")
    df.to_csv(out, index=False)
    print(f"✅ Saved: {out}")

    # --- Draft text for Results section ---
    print("\n=== DRAFT TEXT FOR RESULTS SECTION ===\n")
    for t in tests:
        if t["fold_enrichment"] > 2 and t["hypergeom_p"] < 0.05:
            p_str = f"{t['hypergeom_p']:.2e}" if t["hypergeom_p"] < 1e-3 else f"{t['hypergeom_p']:.3f}"
            print(f"The {t['test']} was significantly enriched "
                  f"(observed={t['observed_overlap']} vs expected={t['expected_overlap']}, "
                  f"{t['fold_enrichment']}-fold, hypergeometric p={p_str}).")


if __name__ == "__main__":
    main()
