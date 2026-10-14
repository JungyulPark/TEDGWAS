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
fig,ax=plt.subplots(figsize=(10,5))
for oc,label,color,marker in zip(order,labels,['#222222','#167b8c','#8757a3','#b86319','#286729'],['o','s','^','D','v']):
    x=col[col.gene.eq('TSHR')&col.outcome.eq(oc)&col.scenario.eq('original_reference')].sort_values('p12')
    ax.plot(x.p12,x['PP.H4'],label=label,marker=marker,color=color)
ax.axhline(.8,color='black',ls='--',lw=1);ax.set_xscale('log');ax.set_ylim(0,1.05);ax.set_xticks([1e-6,5e-6,1e-5],['1e-6','5e-6','1e-5']);ax.set_xlabel('Prior probability of a shared variant (p12)');ax.set_ylabel('PP.H4');ax.set_title('TSHR: sensitivity to the shared-variant prior',weight='bold');ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),fontsize=9);ax.grid(alpha=.15);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(FIG/f'Figure_R3_TSHR_prior.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
