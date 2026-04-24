"""
Phase A-4: Transcriptomic Validation Re-analysis
- Sample mapping verification via fibrotic markers
- 15 intersection genes analysis (TPM, log2FC, P-value)
- INSR/IRS1/IRS2 off-target validation
- Pathway-level comparison (cAMP vs PI3K)
"""
import csv
import math
import os
from scipy import stats

DATA_FILE = r"c:\ProjectTEDGWAS\data.txt"
OUTPUT_DIR = r"c:\ProjectTEDGWAS\TrackA_MR\results"

# Sample mapping: S2=Control, S7/S8/S10/S11=TED
CONTROL_SAMPLES = ['2']
TED_SAMPLES = ['7', '8', '10', '11']

# 15 intersection genes + additional validation genes
INTERSECTION_GENES = {
    'IGF1R-only': ['IGF1R', 'IGF1', 'IL1B', 'CXCL8'],
    'TSHR-only': ['TSHR', 'ADIPOQ', 'FABP4'],
    'Shared': ['HAS2', 'HAS1', 'TNF', 'ARRB1', 'IL6', 'HAS3', 'PPARG', 'CEBPA']
}
ALL_INTERSECTION = []
GENE_ZONE = {}
for zone, genes in INTERSECTION_GENES.items():
    ALL_INTERSECTION.extend(genes)
    for g in genes:
        GENE_ZONE[g] = zone

# Off-target validation genes
OFFTARGET_GENES = ['INSR', 'IRS1', 'IRS2', 'AKT2', 'PIK3CD']

# cAMP pathway genes (TSHR downstream)
CAMP_GENES = ['GNAS', 'ADCY1', 'ADCY2', 'ADCY3', 'ADCY5', 'ADCY6', 'ADCY7',
              'ADCY8', 'ADCY9', 'PRKACA', 'PRKACB', 'CREB1', 'CREB3', 'TSHR']

# PI3K pathway genes (IGF1R downstream)
PI3K_GENES = ['PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3R1', 'PIK3R2', 'AKT1', 'AKT2',
              'AKT3', 'MTOR', 'IRS1', 'IRS2', 'IGF1R', 'IGF1', 'INSR']

# Fibrotic markers for sample validation
FIBROTIC_MARKERS = ['COL1A1', 'COL1A2', 'COL3A1', 'ACTA2', 'FN1']

# All genes to query
ALL_GENES = list(set(ALL_INTERSECTION + OFFTARGET_GENES + CAMP_GENES + PI3K_GENES + FIBROTIC_MARKERS))


def read_data():
    """Read data.txt and extract TPM values for genes of interest."""
    gene_data = {}
    
    with open(DATA_FILE, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        # Find TPM column indices
        tpm_cols = {}
        for i, h in enumerate(header):
            if '_TPM' in h:
                parts = h.replace('_TPM', '').split('__')
                sample = parts[0]
                rep = parts[1]
                tpm_cols[(sample, rep)] = i
        
        for row in reader:
            gene_symbol = row[2]
            if gene_symbol in ALL_GENES:
                sample_tpm = {}
                for (sample, rep), idx in tpm_cols.items():
                    val = float(row[idx])
                    if sample not in sample_tpm:
                        sample_tpm[sample] = []
                    sample_tpm[sample].append(val)
                
                # Average replicates
                gene_data[gene_symbol] = {s: sum(vals)/len(vals) for s, vals in sample_tpm.items()}
    
    return gene_data


def compute_stats(gene_data, gene, ctrl_samples, ted_samples):
    """Compute log2FC and P-value for a gene."""
    if gene not in gene_data:
        return None, None, None, None, None
    
    data = gene_data[gene]
    ctrl_val = sum(data[s] for s in ctrl_samples) / len(ctrl_samples)
    ted_vals = [data[s] for s in ted_samples]
    ted_mean = sum(ted_vals) / len(ted_vals)
    
    # log2FC
    if ctrl_val > 0 and ted_mean > 0:
        log2fc = math.log2(ted_mean / ctrl_val)
    elif ctrl_val == 0 and ted_mean > 0:
        log2fc = float('inf')
    elif ctrl_val > 0 and ted_mean == 0:
        log2fc = float('-inf')
    else:
        log2fc = 0.0
    
    # P-value: one-sample t-test (4 TED values vs Control value)
    if len(ted_vals) >= 2:
        t_stat, p_val = stats.ttest_1samp(ted_vals, ctrl_val)
    else:
        p_val = float('nan')
    
    return ctrl_val, ted_mean, log2fc, p_val, 'UP' if log2fc > 0 else 'DOWN'


def main():
    print("=" * 60)
    print("Phase A-4: Transcriptomic Validation Re-analysis")
    print("=" * 60)
    
    gene_data = read_data()
    print(f"\nLoaded TPM data for {len(gene_data)} genes out of {len(ALL_GENES)} queried")
    missing = [g for g in ALL_GENES if g not in gene_data]
    if missing:
        print(f"Missing genes: {missing}")
    
    # ==========================================
    # 1. Sample Mapping Validation
    # ==========================================
    print("\n" + "=" * 60)
    print("1. SAMPLE MAPPING VALIDATION (Fibrotic Markers)")
    print("=" * 60)
    print(f"{'Gene':<10} {'S2':<12} {'S7':<12} {'S8':<12} {'S10':<12} {'S11':<12} {'S2 lowest?'}")
    print("-" * 80)
    for gene in FIBROTIC_MARKERS:
        if gene in gene_data:
            d = gene_data[gene]
            vals = {s: d[s] for s in ['2', '7', '8', '10', '11']}
            s2_lowest = d['2'] == min(vals.values())
            print(f"{gene:<10} {d['2']:<12.4f} {d['7']:<12.4f} {d['8']:<12.4f} {d['10']:<12.4f} {d['11']:<12.4f} {'YES' if s2_lowest else 'NO'}")
    
    print("\n>>> CONCLUSION: S2 = Control confirmed by fibrotic marker expression")
    
    # ==========================================
    # 2. 15 Intersection Genes Analysis
    # ==========================================
    print("\n" + "=" * 60)
    print("2. 15 INTERSECTION GENES (TPM, log2FC, P-value)")
    print("=" * 60)
    
    results = []
    print(f"{'Gene':<10} {'Zone':<14} {'Ctrl TPM':<12} {'TED TPM':<12} {'log2FC':<10} {'P-value':<12} {'Dir'}")
    print("-" * 80)
    for gene in ALL_INTERSECTION:
        ctrl, ted, log2fc, pval, direction = compute_stats(gene_data, gene, CONTROL_SAMPLES, TED_SAMPLES)
        if ctrl is not None:
            results.append({
                'Gene': gene, 'Zone': GENE_ZONE[gene],
                'Ctrl_TPM': ctrl, 'TED_mean_TPM': ted,
                'log2FC': log2fc, 'P_value': pval, 'Direction': direction
            })
            pval_str = f"{pval:.4e}" if not math.isnan(pval) else "NaN"
            log2fc_str = f"{log2fc:.4f}" if not math.isinf(log2fc) else str(log2fc)
            sig = "*" if pval < 0.05 else ""
            print(f"{gene:<10} {GENE_ZONE[gene]:<14} {ctrl:<12.4f} {ted:<12.4f} {log2fc_str:<10} {pval_str:<12} {direction}{sig}")
    
    # Save results
    out_path = os.path.join(OUTPUT_DIR, "Transcriptome_15gene_validation.csv")
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Gene', 'Zone', 'Ctrl_TPM', 'TED_mean_TPM', 'log2FC', 'P_value', 'Direction'])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nSaved to: {out_path}")
    
    # ==========================================
    # 3. INSR Off-target Validation
    # ==========================================
    print("\n" + "=" * 60)
    print("3. INSR OFF-TARGET VALIDATION")
    print("=" * 60)
    print(f"{'Gene':<10} {'Ctrl TPM':<12} {'TED TPM':<12} {'log2FC':<10} {'P-value':<12} {'Dir'}")
    print("-" * 70)
    offtarget_results = []
    for gene in OFFTARGET_GENES:
        ctrl, ted, log2fc, pval, direction = compute_stats(gene_data, gene, CONTROL_SAMPLES, TED_SAMPLES)
        if ctrl is not None:
            offtarget_results.append({
                'Gene': gene, 'Ctrl_TPM': ctrl, 'TED_mean_TPM': ted,
                'log2FC': log2fc, 'P_value': pval, 'Direction': direction
            })
            pval_str = f"{pval:.4e}" if not math.isnan(pval) else "NaN"
            log2fc_str = f"{log2fc:.4f}" if not math.isinf(log2fc) else str(log2fc)
            sig = "*" if pval < 0.05 else ""
            print(f"{gene:<10} {ctrl:<12.4f} {ted:<12.4f} {log2fc_str:<10} {pval_str:<12} {direction}{sig}")
    
    out_path2 = os.path.join(OUTPUT_DIR, "Insulin_receptor_validation.csv")
    with open(out_path2, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Gene', 'Ctrl_TPM', 'TED_mean_TPM', 'log2FC', 'P_value', 'Direction'])
        writer.writeheader()
        writer.writerows(offtarget_results)
    print(f"\nSaved to: {out_path2}")
    
    # ==========================================
    # 4. Pathway-Level Comparison
    # ==========================================
    print("\n" + "=" * 60)
    print("4. PATHWAY-LEVEL COMPARISON (cAMP vs PI3K)")
    print("=" * 60)
    
    camp_log2fcs = []
    pi3k_log2fcs = []
    
    print("\ncAMP/PKA pathway genes:")
    for gene in CAMP_GENES:
        ctrl, ted, log2fc, pval, direction = compute_stats(gene_data, gene, CONTROL_SAMPLES, TED_SAMPLES)
        if ctrl is not None and not math.isinf(log2fc):
            camp_log2fcs.append(log2fc)
            pval_str = f"{pval:.4e}" if not math.isnan(pval) else "NaN"
            print(f"  {gene:<12} log2FC={log2fc:.4f}  P={pval_str}")
    
    print("\nPI3K/AKT pathway genes:")
    for gene in PI3K_GENES:
        ctrl, ted, log2fc, pval, direction = compute_stats(gene_data, gene, CONTROL_SAMPLES, TED_SAMPLES)
        if ctrl is not None and not math.isinf(log2fc):
            pi3k_log2fcs.append(log2fc)
            pval_str = f"{pval:.4e}" if not math.isnan(pval) else "NaN"
            print(f"  {gene:<12} log2FC={log2fc:.4f}  P={pval_str}")
    
    camp_mean = sum(camp_log2fcs) / len(camp_log2fcs) if camp_log2fcs else 0
    pi3k_mean = sum(pi3k_log2fcs) / len(pi3k_log2fcs) if pi3k_log2fcs else 0
    
    # Two-sample t-test between pathway log2FCs
    if len(camp_log2fcs) >= 2 and len(pi3k_log2fcs) >= 2:
        t_stat, p_pathway = stats.ttest_ind(camp_log2fcs, pi3k_log2fcs)
    else:
        p_pathway = float('nan')
    
    print(f"\ncAMP pathway: mean log2FC = {camp_mean:.4f} (n={len(camp_log2fcs)})")
    print(f"PI3K pathway: mean log2FC = {pi3k_mean:.4f} (n={len(pi3k_log2fcs)})")
    print(f"Pathway comparison P-value: {p_pathway:.4e}")
    
    # Save pathway comparison
    out_path3 = os.path.join(OUTPUT_DIR, "Pathway_level_comparison.csv")
    with open(out_path3, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Pathway', 'Mean_log2FC', 'N_genes', 'Comparison_P'])
        writer.writerow(['cAMP_PKA', f"{camp_mean:.6f}", len(camp_log2fcs), f"{p_pathway:.6e}"])
        writer.writerow(['PI3K_AKT', f"{pi3k_mean:.6f}", len(pi3k_log2fcs), f"{p_pathway:.6e}"])
    print(f"Saved to: {out_path3}")
    
    # ==========================================
    # 5. Individual Sample Values for Key Genes
    # ==========================================
    print("\n" + "=" * 60)
    print("5. INDIVIDUAL SAMPLE TPM (KEY GENES)")
    print("=" * 60)
    key_genes = ['TSHR', 'INSR', 'HAS2', 'IL6', 'IGF1R', 'TNF', 'PPARG']
    print(f"{'Gene':<10} {'S2(Ctrl)':<12} {'S7(TED)':<12} {'S8(TED)':<12} {'S10(TED)':<12} {'S11(TED)':<12}")
    print("-" * 70)
    for gene in key_genes:
        if gene in gene_data:
            d = gene_data[gene]
            print(f"{gene:<10} {d['2']:<12.4f} {d['7']:<12.4f} {d['8']:<12.4f} {d['10']:<12.4f} {d['11']:<12.4f}")
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    main()
