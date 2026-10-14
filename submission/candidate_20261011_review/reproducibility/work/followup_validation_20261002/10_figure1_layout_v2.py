from coloc_core_v1 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import shutil
O=W.parent/'outputs/TED_Followup_Validation_20261002'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':230})
t=pd.read_csv(T/'TSHR_v2/TSHR_240_coloc_settings.csv');order=['BBJ','E4_GRAVES_STRICT','UKB','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
labels=['BBJ GD','FinnGen GD strict','UKB hyperthyroidism','FinnGen GO','FinnGen GO strict'];colors=['#246991','#b65322','#606b7d','#497b45','#735a93']
windows=['full_source_cis','GenePos_pm500kb','GenePos_pm250kb','five_outcome_intersection']
matrix=t[t.scenario.eq('original_reference')&t.p12.eq(1e-5)].pivot(index='outcome',columns='window',values='PP.H4').loc[order,windows]
fig,(ax,bx)=plt.subplots(1,2,figsize=(14,6.3),gridspec_kw={'width_ratios':[1.1,1]})
ax.imshow(matrix,cmap=LinearSegmentedColormap.from_list('pale',['#ffffff','#b5d2e2']),vmin=0,vmax=1,aspect='auto')
for i in range(5):
    for j in range(4):ax.text(j,i,f'{matrix.iloc[i,j]:.3f}',ha='center',va='center',fontsize=12,color='black')
ax.set_yticks(range(5),labels,fontsize=10);ax.set_xticks(range(4),['Full source\ncis region','Gene position\n±500 kb','Gene position\n±250 kb','Common SNPs\nacross 5 outcomes'],fontsize=9)
ax.set_xticks(np.arange(-.5,4,1),minor=True);ax.set_yticks(np.arange(-.5,5,1),minor=True);ax.grid(which='minor',color='white',lw=2);ax.tick_params(which='minor',bottom=False,left=False)
for s in ax.spines.values():s.set_visible(False)
ax.set_title('A   PP.H4 across predefined windows',loc='left',fontweight='bold',pad=15)
for oc,label,color in zip(order,labels,colors):
    q=t[t.scenario.eq('original_reference')&t.window.eq('full_source_cis')&t.outcome.eq(oc)].sort_values('p12')
    bx.plot(q.p12,q['PP.H4'],marker='o',color=color,label=label,lw=1.8)
bx.set_ylim(0,1.06);bx.axhline(.8,color='black',ls='--',lw=1);bx.grid(axis='y',alpha=.15);bx.set_xscale('log');bx.set_xticks([1e-6,5e-6,1e-5],['10⁻⁶','5 × 10⁻⁶','10⁻⁵']);bx.minorticks_off()
bx.set_ylabel('PP.H4: shared-variant posterior probability');bx.set_xlabel('Shared-variant prior probability per SNP (p12)');bx.set_title('B   Sensitivity to the shared-variant prior',loc='left',fontweight='bold',pad=15)
fig.legend(*bx.get_legend_handles_labels(),loc='lower center',ncol=3,bbox_to_anchor=(.55,.10),frameon=False,labelcolor='black')
fig.suptitle('Follow-up Figure 1. TSHR colocalization differs between disease datasets',fontweight='bold',fontsize=14,y=.98)
fig.subplots_adjust(left=.17,right=.97,top=.84,bottom=.32,wspace=.28)
fig.text(.17,.035,'Reference frequencies. Panel A: p12 = 10⁻⁵. Panel B: full source cis region.\nDashed line: PP.H4 = 0.80, not a P-value threshold. FinnGen endpoints overlap.',fontsize=10)
for ext in ['png','pdf']:fig.savefig(O/'Figures'/f'Figure_F1_TSHR_windows_and_priors.{ext}',bbox_inches='tight',facecolor='white')
plt.close(fig)
shutil.copy2(Path(__file__),O/'Verification/Scripts'/Path(__file__).name)
with (O/'Verification/VERSION_NOTES.md').open('a',encoding='utf8') as f:f.write('\nFigure F1 was relaid out by 10_figure1_layout_v2.py after visual inspection to separate nearby value labels; all data are unchanged. Run this after 09_build_report_v1.py to reproduce the final figure.\n')
print('Figure F1 layout revised; data unchanged.')
