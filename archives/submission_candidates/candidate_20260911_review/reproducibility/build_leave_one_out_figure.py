"""Figure S2: all 15 SNP omissions and five full-set comparators."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.lines import Line2D
D=Path(__file__).resolve().parents[1];P=D/'provenance';F=D/'figures'
data=D/'Supplementary_Data_6_Leave_one_out.csv';df=pd.read_csv(data)
assert len(df)==20 and sum(df.excluded_SNP!='None (all instruments)')==15
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'text.color':'black',
                     'axes.labelcolor':'black','xtick.color':'black','ytick.color':'black',
                     'pdf.fonttype':42,'mathtext.fontset':'dejavusans'})
fig=plt.figure(figsize=(8.3,10.0))
gs=fig.add_gridspec(1,4,left=.035,right=.98,top=.855,bottom=.105,width_ratios=[2.28,1.85,1.84,1.08],wspace=.07)
label,forest,effect,pvals=[fig.add_subplot(gs[0,i]) for i in range(4)]
groups=list(df.groupby(['gene','outcome'],sort=False))
labels={'BBJ_Graves':'BBJ','UKB_hyperthyroid':'UKB','FinnGen_GO':'FinnGen'}
cursor=0;positions=[];headers=[]
for (gene,outcome),g in groups:
    headers.append((cursor,gene+' / '+labels[outcome]));cursor+=.98
    for r in g.itertuples():positions.append((cursor,r));cursor+=1
    cursor+=.6
for ax in [label,forest,effect,pvals]:ax.set_ylim(cursor-.2,-.45)
for ax in [label,effect,pvals]:ax.set_xlim(0,1);ax.axis('off')
forest.set_xscale('log');forest.set_xlim(.04,25);forest.set_yticks([])
forest.axvline(1,color='black',ls=(0,(4,3)),lw=.85,zorder=1)
forest.set_xticks([.05,.2,1,5,20],['0.05','0.2','1','5','20']);forest.tick_params(axis='x',labelsize=9)
forest.set_xlabel('Odds ratio (log scale)',fontsize=9.5)
for side in ['top','left','right']:forest.spines[side].set_visible(False)
fig.text(.035,.975,'Figure S2',fontsize=14,fontweight='bold',va='top')
fig.text(.035,.943,'Sensitivity to excluding individual instruments',fontsize=12,fontweight='bold',va='top')
fig.text(.035,.913,'Original reference-frequency analysis · 15 omissions and 5 full-set estimates',fontsize=9.4,va='top')
for ax,x,title in [(label,0,'SNP excluded'),(effect,.5,'OR (95% CI)'),(pvals,.5,r'$P$ value')]:
    ax.text(x,1.025,title,transform=ax.transAxes,ha='left' if ax==label else 'center',fontweight='bold',fontsize=10.2)
sources=[];rowtexts=[]
for y,title in headers:
    rowtexts.append(label.text(0,y,title,fontweight='bold',fontsize=10.5,va='center'))
def pdisplay(p):
    if p>=.001:return f'{p:.3g}'
    b,e=f'{p:.2e}'.split('e')
    b=b.rstrip('0').rstrip('.')
    return b+'×10'+str(int(e)).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
def pplot(p):
    if p>=.001:return f'{p:.3g}'
    b,e=f'{p:.2e}'.split('e');expr=b.rstrip('0').rstrip('.')+r'{\times}10^{'+str(int(e))+'}'
    return r'$\mathbf{'+expr+'}$'
for y,r in positions:
    full=r.excluded_SNP=='None (all instruments)'
    if full:
        for ax in [label,forest,effect,pvals]:ax.axhspan(y-.40,y+.40,color='#F0F0F0',zorder=0)
    name='All SNPs (n='+str(r.n_iv)+')' if full else r.excluded_SNP
    rowtexts.append(label.text(.025,y,name,va='center',fontsize=10.2,fontweight='bold' if full else 'normal'))
    forest.errorbar(r.OR,y,xerr=[[r.OR-r.CI_lower],[r.CI_upper-r.OR]],fmt='D' if full else 'o',
                    markerfacecolor='black' if full else 'white',markeredgecolor='black',
                    color='black',capsize=2,lw=1.0,markersize=4.6,zorder=3)
    shown=f'{r.OR:.3f} ({r.CI_lower:.3f}–{r.CI_upper:.3f})'
    rowtexts.append(effect.text(.5,y,shown,ha='center',va='center',fontsize=9.5))
    rowtexts.append(pvals.text(.5,y,pplot(r.pvalue),ha='center',va='center',fontsize=10,fontweight='bold' if r.pvalue<.05 else 'normal'))
    sources.append({'gene':r.gene,'outcome':r.outcome,'excluded_SNP':r.excluded_SNP,'n_iv':r.n_iv,
                    'OR':r.OR,'CI_lower':r.CI_lower,'CI_upper':r.CI_upper,'pvalue':r.pvalue,
                    'displayed_effect':shown,'displayed_p':pdisplay(r.pvalue),'p_bold':bool(r.pvalue<.05)})
handles=[Line2D([],[],marker='D',color='black',lw=0,label='All SNPs',markersize=5),
         Line2D([],[],marker='o',color='black',mfc='white',lw=0,label='One SNP excluded',markersize=5)]
fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.025,.035),frameon=False,ncol=2,fontsize=9.5)
fig.text(.035,.025,r'Intervals are 95% CIs. Bold $P$ values indicate nominal $P$ < 0.05; all omission results are shown.',fontsize=8.5)
fig.text(.035,.010,'Single-SNP starting sets cannot undergo omission. One remaining SNP is analysed by Wald ratio.',fontsize=8.5)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
alltext=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
assert all(matplotlib.colors.to_hex(t.get_color())=='#000000' for t in alltext)
assert all(0<t.get_window_extent(renderer).x0<t.get_window_extent(renderer).x1<fig.bbox.width for t in rowtexts)
assert all(.04<r.CI_lower<=r.OR<=r.CI_upper<25 for _,r in positions)
fig.savefig(F/'FigureS2.png',dpi=300,facecolor='white')
fig.savefig(F/'FigureS2.pdf',facecolor='white')
source={'rows':sources,'omission_rows':15,'full_set_rows':5,'data_sha256':hashlib.sha256(data.read_bytes()).hexdigest(),
        'all_text_black':True,'all_CIs_within_axis':True,'row_labels_within_canvas':True,
        'pdf_sha256':hashlib.sha256((F/'FigureS2.pdf').read_bytes()).hexdigest(),
        'png_sha256':hashlib.sha256((F/'FigureS2.png').read_bytes()).hexdigest()}
(P/'leave_one_out_figure_sources.json').write_text(json.dumps(source,indent=2,ensure_ascii=False)+'\n',encoding='utf8',newline='\r\n')
print('Figure S2: 15 omissions + 5 baselines; all intervals shown in full; black text.')
