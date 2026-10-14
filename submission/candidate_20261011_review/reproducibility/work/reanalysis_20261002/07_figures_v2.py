from pathlib import Path
import sys
T=Path(__file__).resolve().parent;W=T.parent
sys.path.insert(0,str(W/'pydeps'))
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
O=W.parent/'outputs/TED_Reanalysis_20261002';FIG=O/'Figures'
focusdf=pd.read_csv(O/'Tables/Table2_14_genes_MR_and_colocalization.csv')
focus=list(focusdf.gene_symbol.drop_duplicates());col=pd.read_csv(O/'Tables/Colocalization_840_settings.csv')
order=['BBJ_Graves','E4_GRAVES_STRICT','UKB_hyperthyroid','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
labels=['BBJ GD','FinnGen GD strict','UKB hyperthyroidism','FinnGen GO','FinnGen GO strict']
# Shareable scientific plots, generated directly from the exported numeric table.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','pdf.fonttype':42,'ps.fonttype':42})
matrix=focusdf.pivot(index='gene_symbol',columns='outcome',values='PP.H4').reindex(index=focus,columns=order)
fig,ax=plt.subplots(figsize=(11,8));im=ax.imshow(matrix.to_numpy(),vmin=0,vmax=1,cmap='Blues',aspect='auto')
ax.set_xticks(range(5),['BBJ\nGD','FinnGen\nGD strict (new)','UKB\nhyperthyroidism','FinnGen\nGO','FinnGen\nGO strict (new)']);ax.set_yticks(range(14),focus)
for y,g in enumerate(focus):
    for x,oc in enumerate(order):
        val=matrix.loc[g,oc];ax.text(x,y,f'{val:.3f}',ha='center',va='center',color='black',bbox=dict(facecolor='white',alpha=.82,edgecolor='none',pad=1.5))
        if val>=.8:ax.add_patch(Rectangle((x-.48,y-.48),.96,.96,fill=False,edgecolor='black',linewidth=1.8))
ax.set_title('Shared-variant support across five disease datasets',pad=18,weight='bold');fig.colorbar(im,ax=ax,label='PP.H4 (posterior probability)',fraction=.035,pad=.03)
fig.text(.12,.025,'Outline: PP.H4 >= 0.80 at p12 = 1e-5. FinnGen endpoints overlap; they are not independent replications.',fontsize=9)
fig.subplots_adjust(bottom=.12,top=.9,left=.15,right=.92)
for ext in ['png','pdf']:fig.savefig(FIG/f'Figure_R1_Colocalization.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(14.5,5.6),sharey=True)
for ax,g in zip(axs,['TSHR','IGF1R','CTLA4']):
    x=focusdf[focusdf.gene_symbol.eq(g)].set_index('outcome').reindex(order)
    for i,(oc,r) in enumerate(x.iterrows()):
        color='#126a86' if oc in ['E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT'] else '#555555'
        ax.errorbar(r.OR,i,xerr=[[r.OR-r.CI_lower],[r.CI_upper-r.OR]],fmt='o',color=color,capsize=3,markersize=5)
        ax.text(.03,i+.22,f'OR {r.OR:.2f} ({r.CI_lower:.2f}-{r.CI_upper:.2f}); P={r.pvalue:.2g}',transform=ax.get_yaxis_transform(),fontsize=8,color='black',bbox=dict(facecolor='white',edgecolor='none',alpha=.95,pad=1.2),zorder=5)
    ax.axvline(1,color='#777777',ls='--',lw=1);ax.set_xscale('log');ax.set_title(g,fontstyle='italic',weight='bold');ax.set_xlabel('Odds ratio (95% CI), log scale');ax.grid(axis='x',alpha=.18)
axs[0].set_yticks(range(5),['BBJ GD','FinnGen GD strict [new]','UKB hyperthyroidism','FinnGen GO','FinnGen GO strict [new]']);axs[0].set_ylim(4.6,-.45)
fig.suptitle('Blood expression MR: estimates by disease definition',weight='bold')
fig.text(.12,.015,'Approximate 1-SD expression units. Colours distinguish newly analysed endpoints. No cross-endpoint difference test or meta-analysis.',fontsize=9)
fig.tight_layout(rect=[0,.04,1,.93])
for ext in ['png','pdf']:fig.savefig(FIG/f'Figure_R2_Anchor_MR.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(figsize=(8,5))
for oc,label,color,marker in zip(order,labels,['#222222','#167b8c','#8757a3','#b86319','#286729'],['o','s','^','D','v']):
    x=col[col.gene.eq('TSHR')&col.outcome.eq(oc)&col.scenario.eq('original_reference')].sort_values('p12')
    ax.plot(x.p12,x['PP.H4'],label=label,marker=marker,color=color)
ax.axhline(.8,color='black',ls='--',lw=1);ax.set_xscale('log');ax.set_ylim(0,1.05);ax.set_xticks([1e-6,5e-6,1e-5],['1e-6','5e-6','1e-5']);ax.set_xlabel('Prior probability of a shared variant (p12)');ax.set_ylabel('PP.H4');ax.set_title('TSHR: sensitivity to the shared-variant prior',weight='bold');ax.legend(loc='lower right',fontsize=9);ax.grid(alpha=.15);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(FIG/f'Figure_R3_TSHR_prior.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
