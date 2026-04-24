"""
Figure 3 Panel A — Triangulation Matrix (standalone)
Journal-quality matplotlib version for direct use or as Nanobananapro reference.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ===========================================================
# Colors
# ===========================================================
C_GOLD = '#BF9000'
C_LIGHT_GOLD = '#FFF2CC'
C_RED = '#C00000'
C_DARK_GREEN = '#38761D'
C_LIGHT_GREEN = '#93C47D'
C_YELLOW = '#FFD966'
C_CELL_GRAY = '#CCCCCC'

# ===========================================================
# Data (from Final_Numbers_Frozen Excel, Triangulation sheet)
# Each row: (gene, zone, mr_primary_p, mr_rep_p, mr_finn_p,
#           coloc_pp4, rna_log2fc, rna_p, score, note)
# ===========================================================
rows = [
    ('TSHR',   'TSHR-only',   6.79e-17, 3.88e-29, 6.81e-7, 0.985, 2.71,  0.019, 4, ''),
    ('HAS3',   'Shared',      None,     None,     None,    None,  1.15,  0.028, 2, ''),
    ('HAS2',   'Shared',      None,     None,     None,    None,  0.66,  0.048, 2, ''),
    ('INSR',   'Off-target',  None,     None,     None,    None,  0.54,  0.027, 2, '*'),
    ('IRS2',   'Off-target',  None,     None,     None,    None,  0.65,  0.036, 2, '*'),
    ('IGF1R',  'IGF1R-only',  0.053,    0.041,    0.216,   None,  0.80,  0.131, 1, ''),
    ('TNF',    'Shared',      0.429,    3.57e-5,  0.111,   None,  0.35,  0.611, 1, '‡'),
    ('IL6',    'Shared',      None,     None,     None,    None,  3.55,  0.129, 0, ''),
    ('IL1B',   'IGF1R-only',  None,     None,     None,    None,  2.96,  0.313, 0, ''),
    ('CXCL8',  'IGF1R-only',  None,     None,     None,    None,  1.74,  0.349, 0, ''),
    ('Others', '(7 genes)',   None,     None,     None,    None,  None,  None,  0, ''),
]

# ===========================================================
# Formatters with color classification
# ===========================================================
def p_color_text(p):
    if p is None: return C_CELL_GRAY, '—'
    if p < 0.001: return C_DARK_GREEN, f'{p:.2e}'
    if p < 0.05:  return C_LIGHT_GREEN, f'{p:.3f}'
    if p < 0.10:  return C_YELLOW, f'{p:.3f}'
    return C_CELL_GRAY, f'{p:.3f}'

def coloc_color_text(pp4):
    if pp4 is None: return C_CELL_GRAY, '—'
    if pp4 > 0.95:  return C_DARK_GREEN, f'{pp4:.3f}'
    if pp4 > 0.80:  return C_LIGHT_GREEN, f'{pp4:.3f}'
    return C_YELLOW, f'{pp4:.3f}'

def rna_color_text(log2fc, p):
    if log2fc is None or p is None: return C_CELL_GRAY, '—'
    if p < 0.05:  return C_LIGHT_GREEN, f'{log2fc:+.2f} ({p:.3f})'
    if p < 0.10:  return C_YELLOW, f'{log2fc:+.2f} ({p:.3f})'
    return C_CELL_GRAY, f'{log2fc:+.2f} ({p:.3f})'

def score_color_text(s):
    if s == 4:  return C_GOLD, f'{s} / 4'
    if s >= 2:  return C_LIGHT_GREEN, f'{s} / 4'
    if s == 1:  return C_YELLOW, f'{s} / 4'
    return C_CELL_GRAY, f'{s} / 4'

# ===========================================================
# Figure setup — landscape, wide for readable table
# ===========================================================
fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_axes([0.02, 0.08, 0.96, 0.88])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Title
ax.text(50, 96, 'Multi-modal Triangulation Matrix',
        ha='center', fontsize=16, fontweight='bold',
        family='sans-serif')
ax.text(50, 92,
        'Convergence of genetic (MR, colocalization) and transcriptomic evidence across candidate TED genes',
        ha='center', fontsize=10, style='italic', color='#555555',
        family='sans-serif')

# Column layout
col_widths = [9, 13, 12, 12, 12, 11, 15, 10]  # sum = 94
col_xs = [3]
for w in col_widths[:-1]:
    col_xs.append(col_xs[-1] + w)

col_headers = [
    'Gene', 'Zone',
    'MR Primary\n(Graves)',
    'MR Replication\n(Hyperthyroid)',
    'MR Sensitivity\n(FinnGen GO)',
    'Colocalization\nPP.H4',
    'RNA-seq\nlog2FC (P)',
    'Convergence\nScore'
]

row_h = 5.5
header_h = 7.5
y_header = 82
y_top = y_header - header_h

# Header row
for i, (x, w, h) in enumerate(zip(col_xs, col_widths, col_headers)):
    ax.add_patch(Rectangle((x, y_header - header_h), w, header_h,
                            facecolor='#2C2C2C', edgecolor='black',
                            linewidth=0.8, zorder=3))
    ax.text(x + w/2, y_header - header_h/2, h,
            ha='center', va='center', fontsize=8,
            fontweight='bold', color='white',
            family='sans-serif', zorder=4)

# Data rows
for ri, (gene, zone, mr_p, mr_rp, mr_fn, coloc, rna_fc, rna_p, score, note) in enumerate(rows):
    y = y_top - ri * row_h - row_h
    is_tshr = gene == 'TSHR'
    row_bg = C_LIGHT_GOLD if is_tshr else ('white' if ri % 2 == 0 else '#F5F5F5')

    # Special handling for "Others" row
    if gene == 'Others':
        span_w1 = col_widths[0] + col_widths[1]
        ax.add_patch(Rectangle((col_xs[0], y), span_w1, row_h,
                                facecolor=row_bg, edgecolor='#CCCCCC',
                                linewidth=0.5, zorder=2))
        ax.text(col_xs[0] + span_w1/2, y + row_h/2,
                'Other 7 intersection genes',
                ha='center', va='center', fontsize=8, style='italic',
                color='#222222', family='sans-serif', zorder=3)
        rest_x = col_xs[2]
        rest_w = sum(col_widths[2:])
        ax.add_patch(Rectangle((rest_x, y), rest_w, row_h,
                                facecolor=row_bg, edgecolor='#CCCCCC',
                                linewidth=0.5, zorder=2))
        ax.text(rest_x + rest_w/2, y + row_h/2,
                'IGF1, ADIPOQ, FABP4, HAS1, ARRB1, PPARG, CEBPA   —   all scores 0 / 4',
                ha='center', va='center', fontsize=7.5,
                style='italic', color='#666666',
                family='sans-serif', zorder=3)
        continue

    # Gene column
    ax.add_patch(Rectangle((col_xs[0], y), col_widths[0], row_h,
                            facecolor=row_bg, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[0] + col_widths[0]/2, y + row_h/2, gene + note,
            ha='center', va='center', fontsize=9.5,
            fontweight='bold' if score >= 2 or is_tshr else 'normal',
            color='black', family='sans-serif', zorder=3)

    # Zone column
    ax.add_patch(Rectangle((col_xs[1], y), col_widths[1], row_h,
                            facecolor=row_bg, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[1] + col_widths[1]/2, y + row_h/2, zone,
            ha='center', va='center', fontsize=8.5,
            style='italic', color='#333333',
            family='sans-serif', zorder=3)

    # MR Primary
    c, t = p_color_text(mr_p)
    ax.add_patch(Rectangle((col_xs[2], y), col_widths[2], row_h,
                            facecolor=c, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[2] + col_widths[2]/2, y + row_h/2, t,
            ha='center', va='center', fontsize=8.5,
            color='white' if c == C_DARK_GREEN else 'black',
            fontweight='bold' if c in [C_DARK_GREEN, C_LIGHT_GREEN] else 'normal',
            family='sans-serif', zorder=3)

    # MR Rep
    c, t = p_color_text(mr_rp)
    ax.add_patch(Rectangle((col_xs[3], y), col_widths[3], row_h,
                            facecolor=c, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[3] + col_widths[3]/2, y + row_h/2, t,
            ha='center', va='center', fontsize=8.5,
            color='white' if c == C_DARK_GREEN else 'black',
            fontweight='bold' if c in [C_DARK_GREEN, C_LIGHT_GREEN] else 'normal',
            family='sans-serif', zorder=3)

    # MR FinnGen
    c, t = p_color_text(mr_fn)
    ax.add_patch(Rectangle((col_xs[4], y), col_widths[4], row_h,
                            facecolor=c, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[4] + col_widths[4]/2, y + row_h/2, t,
            ha='center', va='center', fontsize=8.5,
            color='white' if c == C_DARK_GREEN else 'black',
            fontweight='bold' if c in [C_DARK_GREEN, C_LIGHT_GREEN] else 'normal',
            family='sans-serif', zorder=3)

    # Coloc
    c, t = coloc_color_text(coloc)
    ax.add_patch(Rectangle((col_xs[5], y), col_widths[5], row_h,
                            facecolor=c, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[5] + col_widths[5]/2, y + row_h/2, t,
            ha='center', va='center', fontsize=8.5,
            color='white' if c == C_DARK_GREEN else 'black',
            fontweight='bold' if c in [C_DARK_GREEN, C_LIGHT_GREEN] else 'normal',
            family='sans-serif', zorder=3)

    # RNA-seq
    c, t = rna_color_text(rna_fc, rna_p)
    ax.add_patch(Rectangle((col_xs[6], y), col_widths[6], row_h,
                            facecolor=c, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[6] + col_widths[6]/2, y + row_h/2, t,
            ha='center', va='center', fontsize=8,
            color='white' if c == C_DARK_GREEN else 'black',
            fontweight='bold' if c == C_LIGHT_GREEN else 'normal',
            family='sans-serif', zorder=3)

    # Score
    c, t = score_color_text(score)
    ax.add_patch(Rectangle((col_xs[7], y), col_widths[7], row_h,
                            facecolor=c, edgecolor='#CCCCCC',
                            linewidth=0.5, zorder=2))
    ax.text(col_xs[7] + col_widths[7]/2, y + row_h/2, t,
            ha='center', va='center', fontsize=9,
            fontweight='bold',
            color='white' if c == C_GOLD else 'black',
            family='sans-serif', zorder=3)

# Footnotes
y_foot = y_top - len(rows) * row_h - 3
ax.text(3, y_foot,
        '* Off-target gene (not in original 15-gene intersection)   '
        '‡ Significant heterogeneity (Cochran Q P = 0.028)',
        ha='left', va='top', fontsize=8, style='italic',
        color='#444444', family='sans-serif')
ax.text(3, y_foot - 4,
        'TSHR achieves 4/4 convergence across all evidence types. No other gene exceeds 2/4.',
        ha='left', va='top', fontsize=9, fontweight='bold',
        color=C_RED, family='sans-serif')

# Legend box (bottom right)
legend_x = 65
legend_y = y_foot - 1
ax.text(legend_x, legend_y, 'Color key:',
        ha='left', va='top', fontsize=8, fontweight='bold',
        family='sans-serif')
legend_items = [
    (C_DARK_GREEN, 'P < 0.001'),
    (C_LIGHT_GREEN, 'P < 0.05'),
    (C_YELLOW, '0.05 ≤ P < 0.10'),
    (C_CELL_GRAY, 'not significant / —'),
]
for i, (c, lbl) in enumerate(legend_items):
    box_x = legend_x + 10 + i * 7.5
    ax.add_patch(Rectangle((box_x, legend_y - 3), 2, 2,
                            facecolor=c, edgecolor='black', linewidth=0.5))
    ax.text(box_x + 2.3, legend_y - 2, lbl,
            ha='left', va='center', fontsize=7,
            family='sans-serif')

# Save
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelA_journal.pdf',
            format='pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelA_journal.png',
            format='png', dpi=600, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig('c:/ProjectTEDGWAS/Manuscript/figures_final/Fig3_PanelA_journal_preview.png',
            format='png', dpi=150, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("Panel A saved: PDF + PNG (600 DPI) + preview")
