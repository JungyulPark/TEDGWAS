from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
O=Path(__file__).resolve().parents[1];P=O/'provenance';F=O/'figures';F.mkdir(exist_ok=True)
mr=pd.read_csv(P/'MR_primary_canonical.csv');inst=pd.read_csv(P/'instruments_verified.csv')
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif','DejaVu Serif'],'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'mathtext.fontset':'stix'})
def save(fig,name):
    fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight',facecolor='white');fig.savefig(F/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)
assert inst.gene_symbol.nunique()==2544
bbj=mr[mr.outcome=='BBJ_Graves'];hits=bbj[bbj.pvalue<.05/2544];assert len(bbj)==2234 and len(hits)==13
fig,ax=plt.subplots(figsize=(7.2,7.1));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
def box(x,y,w,h,text,size=11):
    ax.add_patch(Rectangle((x,y),w,h,fill=False,lw=.85,edgecolor='#555555'));ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=size)
def arrow(x1,y1,x2,y2):ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops={'arrowstyle':'->','lw':.9,'color':'#555555'})
box(.19,.865,.62,.09,'4,462 druggable genes')
box(.19,.72,.62,.09,'2,544 genes with selected\nblood-expression instruments')
box(.19,.575,.62,.09,'2,234 genes with estimable\nBBJ Graves disease associations')
box(.19,.43,.62,.09,'13 discovery genes\n$P$ < 0.05 / 2,544')
for y1,y2 in [(.865,.81),(.72,.665),(.575,.52)]:arrow(.5,y1,.5,y2)
for y,lab in [(.8375,'1,918 without a valid instrument'),(.6925,'310 not estimable in BBJ')]:
    arrow(.5,y,.66,y);ax.text(.67,y,lab,ha='left',va='center',fontsize=8.5)
for x,w,label in [(.015,.285,'2 known loci\nTSHR and CTLA4'),(.3575,.285,'5 MHC-region genes'),(.70,.285,'6 additional candidates')]:
    box(x,.265,w,.10,label,10);arrow(.5,.43,x+w/2,.365)
box(.575,.10,.41,.10,'0 of 6 met shared-variant support\nin both BBJ and FinnGen',9.5);arrow(.8425,.265,.80,.20)
box(.015,.10,.53,.10,'Three selected genes: TSHR · IGF1R · CTLA4,\ncompared in every outcome regardless of discovery $P$',8.8)
ax.text(.5,.025,'Discovery: BBJ GD   |   Comparisons: UKB hyperthyroidism and FinnGen ophthalmopathy',ha='center',fontsize=8.5)
save(fig,'Figure1')
fig=plt.figure(figsize=(9.4,6.7));gs=fig.add_gridspec(1,4,width_ratios=[2.5,2.7,2.1,1.1],wspace=.035)
la,fa,ea,pa=[fig.add_subplot(gs[0,i]) for i in range(4)]
for a in [la,fa,ea,pa]:a.set_ylim(.35,12.45)
for a in [la,ea,pa]:a.set_xlim(0,1);a.axis('off')
fa.set_xscale('log');fa.set_xlim(.035,3.0);fa.set_yticks([]);fa.axvline(1,color='#999999',ls='--',lw=.8);fa.spines['left'].set_visible(False);fa.set_xticks([.05,.1,.2,.5,1,2],[.05,.1,.2,.5,1,2]);fa.tick_params(axis='x',labelsize=9);fa.set_xlabel('Odds ratio (95% CI)\nLogarithmic scale',fontsize=10)
ea.text(.5,12.2,'OR (95% CI)',ha='center',fontweight='bold');pa.text(.5,12.2,'$\\bf{\\it{P}}$ value',ha='center',fontweight='bold')
colors=['#202020','#28799e','#a95038'];labels=['BBJ Graves disease','UKB hyperthyroidism','FinnGen ophthalmopathy'];ocs=['BBJ_Graves','UKB_hyperthyroid','FinnGen_GO']
records=[]
def pv(p):
    """The displayed value, as it reads. This is what provenance records."""
    if p>=.001:return f'{p:.3g}'
    b,e=f'{p:.2e}'.split('e');b=b.rstrip('0').rstrip('.')
    return b+'\u00d710'+str(int(e)).translate(str.maketrans('-0123456789','\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079'))
def pv_draw(p):
    """The same value as mathtext. Liberation Serif -- the metric-compatible Times
    substitute on Linux -- has no U+207B, so a Unicode superscript minus renders as a
    missing-glyph box there. STIX math ships with matplotlib, so this draws
    identically on every machine."""
    if p>=.001:return f'{p:.3g}'
    b,e=f'{p:.2e}'.split('e');b=b.rstrip('0').rstrip('.')
    return '$\\mathbf{'+b+'{\\times}10^{'+str(int(e))+'}}$'
for g,base in zip(['TSHR','IGF1R','CTLA4'],[11,7,3]):
    la.text(.02,base+.7,g,fontstyle='italic',fontweight='bold',fontsize=11)
    for i,(oc,label) in enumerate(zip(ocs,labels)):
        r=mr[(mr.gene_symbol==g)&(mr.outcome==oc)].iloc[0];y=base-i
        odds,lo,hi=np.exp([r.beta,r.beta-1.96*r.se,r.beta+1.96*r.se])
        la.text(.07,y,label,va='center',fontsize=9.5)
        fa.errorbar(odds,y,xerr=[[odds-lo],[hi-odds]],fmt=['o','s','^'][i],color=colors[i],capsize=3,lw=1.15,markersize=5)
        ea.text(.5,y,f'{odds:.2f} ({lo:.2f}–{hi:.2f})',ha='center',va='center',fontsize=10)
        pa.text(.5,y,pv_draw(r.pvalue),ha='center',va='center',fontweight='bold',fontsize=10,math_fontfamily='stix')
        records.append({'gene':g,'outcome':oc,'or':odds,'lower95':lo,'upper95':hi,'pvalue':r.pvalue,'displayed_p':pv(r.pvalue)})
fig.subplots_adjust(left=.015,right=.995,top=.99,bottom=.12);save(fig,'Figure2')
(P/'clinical_figure_sources.json').write_text(json.dumps({'Figure1':{'druggable_genes':4462,'instrumented_genes':2544,'BBJ_estimable_genes':2234,'discovery_hits':13,'known_loci':2,'MHC_genes':5,'additional_candidates':6,'additional_qualifiers':0},'Figure2':records,'FigureS1':'Preserved, previously inspected descriptive figure; no inferential P values.'},indent=2),encoding='utf-8')
print('Built clinical Figures 1 and 2; preserved Figure S1.')
