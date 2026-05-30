"""
Figure 3 Panel C — Insulin receptor pathway co-upregulation (standalone)
Journal-quality grouped bar plot with INSR novelty callout.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ===== Colors =====
C_BLUE = '#2E75B6'
C_RED = '#C00000'
C_GOLD = '#BF9000'
C_LIGHT_GOLD = '#FFF2CC'

# ===== Data (from Excel Transcriptome sheet) =====
xl = 'c:/ProjectTEDGWAS/Final_Numbers_Frozen_20260422_v2.xlsx'
trans = pd.read_excel(xl, sheet_name='Transcriptome')

genes_c = ['INSR', 'IRS2', 'IRS1', 'AKT2', 'PIK3CD', 'IGF1R']
plot_data = []
for g in genes_c:
    row = trans[trans['Gene'] == g].iloc[0]
    plot_data.append({
        'gene': g,
        'ctrl': row['Ctrl_TPM'],
        'ted': row['TED_mean_TPM'],
        'p': row['P_value'],
        'sig': row['Significant']
    })

# ===== Figure =====
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

n = len(plot_data)
x = np.arange(n)
w = 0.38
ctrl_vals = [d['ctrl'] for d in plot_data]
ted_vals = [d['ted'] for d in plot_data]

bars1 = ax.bar(x - w/2, ctrl_vals, w, label='Control (n = 1)',
               color='#8C8C8C', edgecolor='black', linewidth=1.0)
bars2 = ax.bar(x + w/2, ted_vals, w, label='TED (n = 4, mean)',
               color=C_BLUE, edgecolor='black', linewidth=1.0)

# Significance asterisks
ymax = max(ted_vals)
for i, d in enumerate(plot_data):
    if d['sig'] == 'YES':
        top = max(d['ctrl'], d['ted'])
        ax.text(i, top * 1.10, '*', ha='center', va='bottom',
                fontsize=22, fontweight='bold', color='black',
                family='sans-serif')

# Axes styling
ax.set_xticks(x)
ax.set_xticklabels([d['gene'] for d in plot_data],
                    fontsize=12, fontweight='bold', family='sans-serif')
ax.set_ylabel('Expression (TPM)',
              fontsize=13, fontweight='bold', family='sans-serif')
ax.set_title('Co-upregulation of INSR and insulin signaling effectors in TED orbital fat',
             fontsize=13, fontweight='bold', family='sans-serif', pad=15)
ax.legend(loc='upper right', fontsize=11, frameon=True, edgecolor='black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.2)
ax.spines['bottom'].set_linewidth(1.2)
ax.tick_params(axis='both', labelsize=11, width=1.2)
ax.set_ylim(0, ymax * 1.28)

# P-values below x-axis labels
for i, d in enumerate(plot_data):
    col = C_RED if d['sig'] == 'YES' else '#888888'
    weight = 'bold' if d['sig'] == 'YES' else 'normal'
    ax.text(i, -ymax * 0.085, f"P = {d['p']:.3f}",
            ha='center', va='top', fontsize=10, color=col,
            fontweight=weight, family='sans-serif')

# Novelty callout for INSR
ax.annotate('First report of INSR\nupregulation in TED\norbital tissue',
            xy=(0 + w/2, plot_data[0]['ted']),
            xytext=(1.8, ymax * 0.72),
            fontsize=10.5, fontweight='bold', color=C_RED,
            family='sans-serif',
            bbox=dict(boxstyle='round,pad=0.4', facecolor=C_LIGHT_GOLD,
                      edgecolor=C_RED, linewidth=1.5),
            arrowprops=dict(arrowstyle='->', color=C_RED, lw=1.8))

# Bottom caption
fig.text(0.5, 0.015,
         '* P < 0.05 vs. control. Bulk RNA-seq of orbital adipose tissue, 4 inactive TED patients and 1 control, 2024.',
         ha='center', fontsize=9, style='italic', color='#555555')

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelC_journal.pdf',
            format='pdf', bbox_inches='tight', facecolor='white')
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelC_journal.png',
            format='png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelC_journal_preview.png',
            format='png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("Panel C saved.")
print("Exact values used (from Excel):")
for d in plot_data:
    print(f"  {d['gene']}: Control = {d['ctrl']:.3f} TPM, TED = {d['ted']:.3f} TPM, "
          f"P = {d['p']:.4f}, sig = {d['sig']}")
