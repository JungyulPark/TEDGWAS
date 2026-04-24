"""
Figure 3 Panel B — TSHR expression in TED orbital fat (standalone)
Journal-quality bar plot. Data sourced from Excel Transcriptome sheet.
"""
import matplotlib.pyplot as plt
import pandas as pd

# ===== Colors =====
C_RED = '#C00000'
C_GOLD = '#BF9000'
C_LIGHT_GOLD = '#FFF2CC'

# ===== Exact data (from Excel Transcriptome sheet) =====
# TSHR: Ctrl_TPM = 0.099755, TED_mean_TPM = 0.651395
# log2FC = 2.707071, P = 0.019465
xl = 'c:/ProjectTEDGWAS/Final_Numbers_Frozen_20260422_v2.xlsx'
trans = pd.read_excel(xl, sheet_name='Transcriptome')
tshr = trans[trans['Gene'] == 'TSHR'].iloc[0]

ctrl_val = tshr['Ctrl_TPM']     # 0.100
ted_val = tshr['TED_mean_TPM']  # 0.651
fold = ted_val / ctrl_val        # 6.5x
pval = tshr['P_value']           # 0.019

# ===== Figure =====
fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

bars = ax.bar(['Control\n(n = 1)', 'TED\n(n = 4, mean)'],
              [ctrl_val, ted_val],
              color=['#8C8C8C', C_RED],
              edgecolor='black', linewidth=1.5, width=0.55)

# Value labels on bars
for bar, val in zip(bars, [ctrl_val, ted_val]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
            f'{val:.2f}', ha='center', va='bottom',
            fontsize=14, fontweight='bold', family='sans-serif')

# Fold-change callout
ax.annotate('', xy=(1, ted_val + 0.18), xytext=(0, ctrl_val + 0.18),
            arrowprops=dict(arrowstyle='->', color='black', lw=2))
ax.text(0.5, ted_val + 0.28,
        f'{fold:.1f}× upregulation\n(P = {pval:.3f})',
        ha='center', fontsize=13, fontweight='bold',
        family='sans-serif',
        bbox=dict(boxstyle='round,pad=0.4', facecolor=C_LIGHT_GOLD,
                  edgecolor=C_GOLD, linewidth=2))

# Axes styling
ax.set_ylabel('TSHR expression (TPM)',
              fontsize=13, fontweight='bold', family='sans-serif')
ax.set_title('TSHR upregulation in TED orbital adipose tissue',
             fontsize=14, fontweight='bold', family='sans-serif',
             pad=15)
ax.set_ylim(0, max(ctrl_val, ted_val) * 2.0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.2)
ax.spines['bottom'].set_linewidth(1.2)
ax.tick_params(axis='both', labelsize=11, width=1.2)
ax.tick_params(axis='y', labelsize=11)

# Caption below axes
fig.text(0.5, 0.02,
         'Bulk RNA-seq of orbital adipose tissue from 4 inactive TED patients and 1 control, 2024.',
         ha='center', fontsize=9, style='italic', color='#555555')

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelB_journal.pdf',
            format='pdf', bbox_inches='tight', facecolor='white')
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelB_journal.png',
            format='png', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelB_journal_preview.png',
            format='png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("Panel B saved.")
print(f"Exact values used: Control = {ctrl_val:.3f} TPM, TED = {ted_val:.3f} TPM, "
      f"fold = {fold:.2f}x, P = {pval:.4f}")
