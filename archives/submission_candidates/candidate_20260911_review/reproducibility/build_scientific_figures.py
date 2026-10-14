"""Reproduce all three main figures from preserved aggregate results."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Rectangle,Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
O=Path(__file__).resolve().parents[1];P=O/'provenance';F=O/'figures'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--figures',nargs='+',choices=['Figure1','Figure2','Figure3'],default=['Figure1','Figure2','Figure3'])
selected=set(parser.parse_args().figures)
mr=pd.read_csv(P/'MR_primary_canonical.csv');inst=pd.read_csv(P/'instruments_verified.csv');co=pd.read_csv(P/'coloc_canonical_v2.csv')
OC=['BBJ_Graves','UKB_hyperthyroid','FinnGen_GO'];LABEL=['BBJ Graves disease','UKB hyperthyroidism','FinnGen ophthalmopathy']
GENES=['TSHR','IGF1R','CTLA4'];CAND=['TNFSF14','IFNGR1','MAPKAPK5','HSD3B7','VKORC1','PRSS36']
GC={'TSHR':'#0072B2','IGF1R':'#D55E00','CTLA4':'#6B5B95'};OCOLOR=list(GC.values())
INK='#000000';MUTED='#000000';GRID='#DEE5E9';H=['PP.H2','PP.H3','PP.H4'];HC=['#8DA8BE','#E6B966','#258E83'];OTHER='#E5E9ED'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,'text.color':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.edgecolor':MUTED,'mathtext.fontset':'stix'})
layout_checks={}
def save(fig,name):
    if name in selected:
        # Colour conveys data categories only; every title, label and legend is black.
        fig.canvas.draw()
        labels=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
        for t in labels:t.set_color('#000000')
        fig.canvas.draw()
        layout_checks[name]={'text_items':len(labels),'all_text_black':all(matplotlib.colors.to_hex(t.get_color())=='#000000' for t in labels)}
        fig.savefig(F/(name+'.png'),dpi=300,facecolor='white');fig.savefig(F/(name+'.pdf'),facecolor='white')
    plt.close(fig)
def panel(fig,x,y,letter,title):
    fig.text(x,y,letter,fontsize=14,fontweight='bold',va='top');fig.text(x+.028,y-.001,title,fontsize=10.5,fontweight='bold',va='top')
def rec(g,o,pr=1e-5):return co[(co.gene==g)&(co.outcome==o)&np.isclose(co.p12,pr,rtol=1e-9,atol=0)].iloc[0]
def pv(p):
    if p>=.001:return f'{p:.3g}'
    b,e=f'{p:.2e}'.split('e');return b.rstrip('0').rstrip('.')+'×10'+str(int(e)).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
def pdrawing(p,bold=False):
    if p>=.001:return f'{p:.3g}'
    b,e=f'{p:.2e}'.split('e');expr=b.rstrip('0').rstrip('.')+r'{\times}10^{'+str(int(e))+'}'
    return '$'+(r'\mathbf{'+expr+'}' if bold else expr)+'$'
source={'Figure1':{'druggable_genes':4462,'instrumented_genes':2544,'BBJ_estimable_genes':2234,'discovery_hits':13,'known_loci':2,'MHC_genes':5,'additional_candidates':6,'additional_qualifiers':0},'Figure2':[],'Figure2_posteriors':[],'Figure3_primary':[],'Figure3_priors':[],'FigureS1':'Preserved descriptive figure; no inferential P values.'}
# All gene-level associations, anchored to the strongest selected eQTL instrument.
bbj=mr[mr.outcome=='BBJ_Graves'].copy();hits=bbj[bbj.pvalue<.05/2544];assert len(bbj)==2234 and len(hits)==13
assert inst.gene_symbol.nunique()==2544 and inst.groupby('gene_symbol').chr.nunique().max()==1
anchors=inst.loc[inst.groupby('gene_symbol').zscore.apply(lambda s:s.abs().idxmax()),['gene_symbol','chr','pos_hg19']]
screen=bbj.merge(anchors,on='gene_symbol',validate='one_to_one').sort_values(['chr','pos_hg19']);assert len(screen)==len(bbj)
offset=0;ticks=[];chrom=[]
for c in sorted(screen.chr.unique()):
    ix=screen.chr==c;pos=screen.loc[ix,'pos_hg19']/1e6;screen.loc[ix,'x']=pos+offset;ticks.append(offset+pos.max()/2);chrom.append(str(int(c)));offset+=pos.max()+18
screen['y']=-np.log10(screen.pvalue);source['Figure1_screen']=screen[['gene_symbol','chr','pos_hg19','pvalue','x','y']].to_dict('records')
fig=plt.figure(figsize=(11.8,6.4));flow=fig.add_axes([.035,.09,.275,.81]);ax=fig.add_axes([.385,.18,.59,.69])
panel(fig,.03,.965,'A','Selection and discovery');panel(fig,.36,.965,'B','Druggable-gene MR associations in BBJ')
flow.set(xlim=(0,1),ylim=(0,1));flow.axis('off')
def box(y,n,label,fill='#EEF3F6'):
    flow.add_patch(FancyBboxPatch((.04,y),.90,.13,boxstyle='round,pad=0.012,rounding_size=0.015',fc=fill,ec='none'));flow.text(.085,y+.085,n,fontsize=19,fontweight='bold',va='center');flow.text(.085,y+.030,label,fontsize=9,va='center')
box(.80,'4,462','Druggable genes');box(.585,'2,544','Genes with selected instruments');box(.37,'2,234','Genes estimable in BBJ');box(.155,'13','Discovery genes','#E4F1EE')
for a,b,lab in [(.80,.715,'1,918 without a valid instrument'),(.585,.50,'310 not estimable in BBJ'),(.37,.285,r'$P$ < 0.05 / 2,544')]:
    flow.annotate('',(.49,b),(.49,a),arrowprops={'arrowstyle':'-|>','lw':1,'color':MUTED});flow.text(.51,(a+b)/2,lab,fontsize=8,va='center',ha='left',color=MUTED)
flow.text(.06,.06,'2 known loci  ·  5 MHC genes\n6 additional candidates',fontsize=9.2,linespacing=1.6)
for i,c in enumerate(sorted(screen.chr.unique())):
    r=screen[screen.chr==c];ax.scatter(r.x,r.y,s=8,color=['#CCD6DD','#91A5B2'][i%2],alpha=.85,linewidth=0,zorder=2)
dy={'HLA-A':10,'HLA-DQA2':16,'C4A':22,'TUBB':24,'PSMB8':-20,'IFNGR1':22,'HSD3B7':54,'VKORC1':32,'PRSS36':12,'MAPKAPK5':24,'TNFSF14':23,'CTLA4':12,'TSHR':14}
dx={'HLA-A':0,'HLA-DQA2':-24,'C4A':-25,'TUBB':-34,'PSMB8':-28,'IFNGR1':20,'HSD3B7':0,'VKORC1':-8,'PRSS36':18,'MAPKAPK5':-12,'TNFSF14':14,'CTLA4':0,'TSHR':0}
for g in hits.gene_symbol:
    r=screen[screen.gene_symbol==g].iloc[0];ax.scatter(r.x,r.y,s=25,fc=GC.get(g,'#344E5E'),ec='white',lw=.4,zorder=4);ax.annotate(g,(r.x,r.y),xytext=(dx[g],dy[g]),textcoords='offset points',ha='center',fontsize=7.8,fontstyle='italic',color=GC.get(g,INK),arrowprops={'arrowstyle':'-','lw':.55,'color':MUTED},zorder=5)
r=screen[screen.gene_symbol=='IGF1R'].iloc[0];ax.scatter(r.x,r.y,s=48,fc='#A33E00',ec='black',lw=.45,zorder=4);ax.annotate('IGF1R',(r.x,r.y),xytext=(0,15),textcoords='offset points',ha='center',va='bottom',fontsize=8.8,fontstyle='italic',fontweight='bold',color=INK,arrowprops={'arrowstyle':'-','color':INK,'lw':.55},zorder=6)
ax.axhline(-np.log10(.05/2544),color='#D55E00',ls=(0,(4,3)),lw=.85,zorder=1)
ax.set(xlim=(-25,offset),ylim=(-.4,25.5),ylabel=r'Association strength ($-\log_{10}P$)',xlabel='Chromosome of strongest selected expression instrument');ax.set_xticks(ticks,chrom);ax.tick_params(axis='x',labelsize=8);ax.set_yticks([0,5,10,15,20,25]);ax.grid(axis='y',color=GRID,lw=.6,zorder=0)
fig.text(.385,.06,'One point per gene  |  All 2,234 estimable genes shown  |  13 discovery genes labelled',fontsize=8.5,color=MUTED);save(fig,'Figure1')
# Separate grids and a reserved gutter distinguish A from B while preserving matched rows.
with plt.rc_context({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans','text.color':'#202020','axes.labelcolor':'#202020','xtick.color':'#555555','axes.edgecolor':'#777777'}):
    fig=plt.figure(figsize=(12.5,7.2))
    ag=fig.add_gridspec(1,4,left=.025,right=.650,top=.85,bottom=.20,width_ratios=[2.2,2.05,1.75,1.1],wspace=.055)
    bg=fig.add_gridspec(1,2,left=.730,right=.982,top=.85,bottom=.20,width_ratios=[2.65,.60],wspace=.08)
    la,fa,ea,pa=[fig.add_subplot(ag[0,i]) for i in range(4)]
    ba,ha=[fig.add_subplot(bg[0,i]) for i in range(2)]
    panel(fig,.025,.965,'A','Association estimates')
    panel(fig,.730,.965,'B','Regional genetic evidence')
    fig.add_artist(Line2D([.690,.690],[.045,.965],transform=fig.transFigure,color='#D1D1D1',lw=.8))
    for ax in [la,fa,ea,pa,ba,ha]:ax.set_ylim(-.9,10.6)
    for ax in [la,ea,pa,ha]:ax.set_xlim(0,1);ax.axis('off')
    fa.set_xscale('log');fa.set_xlim(.035,3.0);fa.set_yticks([]);fa.axvline(1,color='#777777',ls='--',lw=.8);fa.spines['left'].set_visible(False)
    fa.set_xticks([.05,.1,.2,.5,1,2],[.05,.1,.2,.5,1,2]);fa.tick_params(axis='x',labelsize=8);fa.set_xlabel('Odds ratio (log scale)',fontsize=8.5)
    ba.set_xlim(0,1);ba.set_yticks([]);ba.spines['left'].set_visible(False);ba.set_xticks([0,.5,1]);ba.set_xlabel('Posterior probability',fontsize=8.5)
    ea.text(.5,10.3,'OR (95% CI)',ha='center',fontweight='bold',fontsize=9)
    pa.text(.5,10.3,r'$P$ value',ha='center',fontweight='bold',fontsize=9)
    ha.text(.5,10.3,'H4',ha='center',fontweight='bold',fontsize=9)
    # Neutral rules separate gene groups; all A-panel marks and statistics are monochrome.
    for y in [6.45,2.95]:
        fy=la.transData.transform((0,y))[1]/fig.bbox.height
        for left,right in [(.025,.650),(.730,.982)]:
            fig.add_artist(Line2D([left,right],[fy,fy],transform=fig.transFigure,color='#DDDDDD',lw=.6))
    for g,base in zip(GENES,[9,5.5,2]):
        la.text(.025,base+.60,g,fontstyle='italic',fontweight='bold',color='#202020',fontsize=10.5)
        for i,(oc,label) in enumerate(zip(OC,LABEL)):
            r=mr[(mr.gene_symbol==g)&(mr.outcome==oc)].iloc[0];y=base-i;odds,lo,hi=np.exp([r.beta,r.beta-1.96*r.se,r.beta+1.96*r.se])
            la.text(.045,y,label,va='center',fontsize=8.0)
            fa.errorbar(odds,y,xerr=[[odds-lo],[hi-odds]],fmt=['o','s','^'][i],color='#202020',capsize=2.5,lw=1.15,markersize=5)
            ea.text(.5,y,f'{odds:.2f} ({lo:.2f}–{hi:.2f})',ha='center',va='center',fontsize=8.3)
            nominal=bool(r.pvalue<.05);discovery=bool(oc=='BBJ_Graves' and r.pvalue<.05/2544)
            marks=('*' if nominal else '')+('†' if discovery else '')
            pa.text(.5,y,pdrawing(r.pvalue,bold=nominal)+(' '+marks if marks else ''),ha='center',va='center',fontweight='bold' if nominal else 'normal',fontsize=8.8)
            source['Figure2'].append({'gene':g,'outcome':oc,'or':odds,'lower95':lo,'upper95':hi,'pvalue':r.pvalue,'displayed_p':pv(r.pvalue),'nominal_significant':nominal,'BBJ_discovery_significant':discovery,'p_bold':nominal,'significance_markers':marks})
            c=rec(g,oc);left=0
            for key,color in zip(H,HC):ba.barh(y,c[key],left=left,height=.46,color=color,edgecolor='white',lw=.35);left+=c[key]
            ba.barh(y,c['PP.H0']+c['PP.H1'],left=left,height=.46,color=OTHER,edgecolor='none')
            ha.text(.5,y,f"{c['PP.H4']:.3f}",ha='center',va='center',fontsize=8.7,fontweight='bold' if c['PP.H4']>=.8 else 'normal',color='#202020')
            source['Figure2_posteriors'].append({'gene':g,'outcome':oc,**{key:float(c[key]) for key in ['PP.H0','PP.H1',*H]}})
    fig.legend(handles=[Patch(fc=c,label=l) for c,l in zip(HC+[OTHER],['Expression only (H2)','Distinct variants (H3)','Shared variant (H4)','Other (H0/H1)'])],loc='lower left',bbox_to_anchor=(.730,.015),ncol=1,frameon=False,fontsize=7.6,borderaxespad=0,labelspacing=.35)
    fig.text(.03,.075,'* Nominal $P$ < 0.05 (bold); † BBJ discovery $P$ < 0.05/2,544\nThe asterisk alone does not indicate significance after multiple-testing correction.',fontsize=8.1,color=INK,linespacing=1.6)
    save(fig,'Figure2')
# Complete candidate matrix and the selected-gene prior sensitivity.
fig=plt.figure(figsize=(11.8,6.5));heat=fig.add_axes([.13,.16,.30,.66]);panel(fig,.03,.965,'A','Shared-variant support across genes');panel(fig,.48,.965,'B','Sensitivity to the shared-association prior')
allg=GENES+CAND;mat=np.array([[rec(g,o)['PP.H4'] for o in OC] for g in allg]);cmap=LinearSegmentedColormap.from_list('shared',['#F7FAFA','#D3E9E5','#8EC7BD'])
row_y=np.arange(len(allg),dtype=float);row_y[3:]+=.28
heat.set(xlim=(-.51,2.51),ylim=(row_y[-1]+.51,-.51))
heat.set_xticks(range(3),['BBJ GD','UKB\nhyperthyroidism','FinnGen\nophthalmopathy'],fontsize=8.5);heat.set_yticks(row_y,allg,fontstyle='italic',fontsize=9.5);heat.tick_params(length=0,pad=7)
for s in heat.spines.values():s.set_visible(False)
cell_labels=[];strong_cells=[]
for i,g in enumerate(allg):
    for k,o in enumerate(OC):
        v=mat[i,k];y=row_y[i]
        heat.add_patch(Rectangle((k-.5,y-.5),1,1,facecolor=cmap(v),edgecolor='white',lw=.5,zorder=1))
        label=heat.text(k,y,'<0.001' if v<.001 else f'{v:.3f}',ha='center',va='center',fontsize=9,color=INK,fontfamily='DejaVu Sans',clip_on=False,zorder=5)
        cell_labels.append((g,o,k,y,label))
        if v>=.8:strong_cells.append((k,y))
        source['Figure3_primary'].append({'gene':g,'outcome':o,'p12':1e-5,'PP.H4':float(v)})
for k,y in strong_cells:heat.add_patch(Rectangle((k-.5,y-.5),1,1,fill=False,ec='black',lw=1.35,zorder=3,clip_on=False))
fig.text(.032,.65,'Selected\ngenes',rotation=90,fontsize=8,color=INK,ha='center');fig.text(.032,.33,'Additional\ncandidates',rotation=90,fontsize=8,color=INK,ha='center')
axes=[]
for i,g in enumerate(GENES):
    ax=fig.add_axes([.50+i*.166,.27,.145,.50]);axes.append(ax)
    for k,o in enumerate(OC):
        vals=[float(rec(g,o,pr)['PP.H4']) for pr in [1e-5,5e-6,1e-6]];ax.plot(range(3),vals,color=OCOLOR[k],marker=['o','s','^'][k],ms=4,lw=1.3)
        source['Figure3_priors'].extend({'gene':g,'outcome':o,'p12':pr,'PP.H4':v} for pr,v in zip([1e-5,5e-6,1e-6],vals))
    ax.axhline(.8,ls='--',color=MUTED,lw=.8);ax.set(xlim=(-.16,2.16),ylim=(-.025,1.07));ax.set_title(g,fontstyle='italic',fontweight='bold',fontsize=11,pad=12,color=GC[g]);ax.set_xticks([0,1,2],['D','I','C']);ax.set_yticks([0,.2,.4,.6,.8,1]);ax.grid(axis='y',color=GRID,lw=.5)
    if i:ax.set_yticklabels([])
    else:ax.set_ylabel('Shared-variant probability',fontsize=9)
axes[0].annotate('0.661',(2,float(rec('TSHR',OC[0],1e-6)['PP.H4'])),xytext=(-5,-24),textcoords='offset points',ha='right',fontsize=8,color=OCOLOR[0],arrowprops={'arrowstyle':'-','color':OCOLOR[0],'lw':.6})
fig.legend(handles=[Line2D([0],[0],color=c,marker=m,ms=4,label=l) for c,m,l in zip(OCOLOR,['o','s','^'],LABEL)],loc='lower center',bbox_to_anchor=(.74,.075),ncol=1,frameon=False,fontsize=8)
fig.add_artist(Rectangle((.13,.083),.013,.014,transform=fig.transFigure,fc='none',ec='black',lw=1.35))
fig.text(.151,.09,'Black border: PP.H4 ≥ 0.80',fontsize=8.1,color=INK,va='center')
fig.text(.13,.059,'Strong shared-variant support (not a $P$ value)\nCombined criterion: BBJ AND FinnGen',fontsize=8.1,color=INK,linespacing=1.5,va='top')
fig.text(.50,.20,'D  Default     I  Intermediate     C  Conservative',fontsize=8.5,color=INK)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
for g,o,k,y,label in cell_labels:
    text_box=label.get_window_extent(renderer);bounds=heat.transData.transform([[k-.5,y-.5],[k+.5,y+.5]])
    lo=bounds.min(axis=0);hi=bounds.max(axis=0)
    assert text_box.x0>lo[0]+2 and text_box.x1<hi[0]-2 and text_box.y0>lo[1]+2 and text_box.y1<hi[1]-2,'Cell label padding: '+g+'/'+o
save(fig,'Figure3')
if 'Figure3' in selected:layout_checks['Figure3'].update({'cell_labels_with_padding':len(cell_labels),'strong_support_borders':len(strong_cells),'border_matches_cell_geometry':True})
source['input_sha256']={n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['MR_primary_canonical.csv','instruments_verified.csv','coloc_canonical_v2.csv']}
(P/'clinical_figure_sources.json').write_text(json.dumps(source,indent=2),encoding='utf-8');print('Built '+', '.join(sorted(selected))+' from preserved results; other figure files retained.')
layout_path=P/'figure_layout_checks.json'
prior=json.loads(layout_path.read_text(encoding='utf-8')) if layout_path.exists() else {}
prior.update(layout_checks);layout_path.write_text(json.dumps(prior,indent=2)+'\n',encoding='utf-8')
