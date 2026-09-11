import sys,shutil,importlib.util,json
from pathlib import Path
W=Path(__file__).resolve().parent;sys.path.insert(0,str(W/'pydeps'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=W.parent;F=O/'figures';F.mkdir(exist_ok=True)
P=O/'provenance'
co=pd.read_csv(P/'coloc_canonical_v2.csv');co=co[np.isclose(co.p12,1e-5)]
mr=pd.read_csv(P/'MR_primary_canonical.csv')
ti=pd.read_csv(P/'tissue_descriptive.csv')
ti=ti[['gene','ctrl_tpm','ted_tpm','ted_per_sample','log2fc_tpm']]
ti.to_csv(P/'tissue_descriptive.csv',index=False)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],'font.size':10,'mathtext.fontset':'stix','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
colors=['#202020','#28799e','#a95038'];ocs=['BBJ_Graves','UKB_hyperthyroid','FinnGen_GO'];labs=['BBJ Graves','UKB hyperthyroid','FinnGen GO']
genes=['TSHR','IGF1R','CTLA4'];candidates=['TNFSF14','IFNGR1','MAPKAPK5','HSD3B7','VKORC1','PRSS36']
def save(fig,name):
    fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight',facecolor='white')
    plt.close(fig)
# Figure 1 derives from the frozen journal-layout builder; only paths and wording change.
source=(W/'figure1_layout.py').read_text(encoding='utf-8')
source=source.replace('open(MASTER).read()',"open(MASTER, encoding='utf-8').read()")
source=source.replace('"Shared-variant support", "both Graves\\nphenotypes", "none"','"PP.H4 ≥ 0.80", "both Graves\\nphenotypes", "no outcome"')
source=source.replace('"discovery-only colocalization"','"discovery-only coloc"')
code=W/'figure1_layout.py';code.write_text(source,encoding='utf-8')
spec=importlib.util.spec_from_file_location('figure1',code);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.MASTER=str(O.parents[1]/'MANUSCRIPT_TED_TRAP_v5_MASTER.md');mod.COLOC=str(P/'coloc_canonical_v2.csv');mod.MR=str(P/'MR_primary_canonical.csv');mod.OUT=[str(F/'Figure1.png'),str(F/'Figure1.pdf')];mod.build()
plt.rcParams.update({'font.size':10})
# Full data above, enlarged hit region below; no undisclosed clipping of discovery estimates.
b=mr[mr.outcome==ocs[0]].copy();b['logp']=-np.log10(b.pvalue);thr=-np.log10(.05/2544)
hits=b[b.pvalue<.05/2544];assert len(hits)==13
mhc={'HLA-A','HLA-DQA2','C4A','TUBB','PSMB8'};cluster={'HSD3B7','VKORC1','PRSS36'}
def color(g):return '#ae3838' if g in {'TSHR','CTLA4'} else '#888888' if g in mhc else '#a999c1' if g in cluster else '#bc8a2e'
fig,(a,z)=plt.subplots(2,1,figsize=(7.3,8.8),gridspec_kw={'height_ratios':[1,1.55]})
for ax in (a,z):
    ax.scatter(b.beta,b.logp,s=9,c='#cccccc',linewidths=0,rasterized=True)
    ax.axhline(thr,c='#555555',ls='--',lw=.8);ax.axvline(0,c='#cccccc',lw=.7,zorder=0)
    for _,r in hits.iterrows():ax.scatter(r.beta,r.logp,s=30,c=color(r.gene_symbol),edgecolors='white',linewidths=.35,zorder=4)
    r=b[b.gene_symbol=='IGF1R'].iloc[0];ax.scatter(r.beta,r.logp,s=38,c='#28799e',marker='D',zorder=5)
    ax.set_ylabel(r'$-\log_{10}(P)$');ax.set_ylim(-.3,25)
a.set_xlim(-10,8);a.set_title('(A) All 2,234 estimable genes',loc='left',fontweight='bold',fontsize=11)
a.set_xlabel('MR effect on log-odds (β)');a.text(-9.6,thr+.55,'Bonferroni P < 0.05/2,544',fontsize=9)
z.set_xlim(-3.35,3.15);z.set_title('(B) Enlarged region containing all 13 discovery hits',loc='left',fontweight='bold',fontsize=11)
z.set_xlabel('MR effect on log-odds (β)')
pos={'HLA-A':(1.0,23.5),'HLA-DQA2':(1.1,18.2),'CTLA4':(-1.8,17.0),'TSHR':(-2.7,14.9),'C4A':(-1.6,11.3),'HSD3B7':(-2.35,8.6),'VKORC1':(-2.8,6.7),'TNFSF14':(-1.05,8.3),'MAPKAPK5':(-.7,3.4),'PSMB8':(.48,2.3),'IFNGR1':(1.45,4.15),'TUBB':(.9,9.1),'PRSS36':(2.1,7.5),'IGF1R':(.8,.8)}
for _,r in b[b.gene_symbol.isin(pos)].iterrows():
    z.annotate(r.gene_symbol,xy=(r.beta,r.logp),xytext=pos[r.gene_symbol],fontsize=9,fontstyle='italic',ha='center',arrowprops={'arrowstyle':'-','lw':.6,'color':'#777777'},zorder=6)
fig.tight_layout(h_pad=2.0);save(fig,'Figure2')
def h4plot(ax,gs):
    ticks=[];labels=[]
    for k,g in enumerate(gs):
        center=(len(gs)-1-k)*4;ticks.append(center);labels.append(g)
        for i,oc in enumerate(ocs):
            y=center+(.8-i*.8);r=co[(co.gene==g)&(co.outcome==oc)].iloc[0];v=r['PP.H4']
            ax.scatter(v,y,color=colors[i],marker=['o','s','^'][i],s=42,label=labs[i] if k==0 else None,zorder=3)
            dx=-.025 if v>.85 else .025;ax.text(v+dx,y,f'{v:.3f}',va='center',ha='right' if dx<0 else 'left',fontsize=9)
    ax.set_yticks(ticks,labels,fontstyle='italic');ax.axvline(.8,color='#888888',ls='--',lw=1);ax.set_xlim(-.045,1.045);ax.set_ylim(-1.5,(len(gs)-1)*4+1.5)
    ax.set_xlabel('Shared-variant posterior probability (PP.H4)');ax.grid(axis='x',alpha=.12);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
fig,(a,f)=plt.subplots(1,2,figsize=(9.5,5.2),gridspec_kw={'width_ratios':[1.1,1]})
h4plot(a,genes);a.set_title('(A) Colocalization',loc='left',fontweight='bold',fontsize=12)
for k,g in enumerate(genes):
    center=(2-k)*4
    for i,oc in enumerate(ocs):
        r=mr[(mr.gene_symbol==g)&(mr.outcome==oc)].iloc[0];y=center+(.8-i*.8)
        f.errorbar(r.beta,y,xerr=1.96*r.se,fmt=['o','s','^'][i],color=colors[i],markersize=5,capsize=3,lw=1.1)
        pp=f'{r.pvalue:.2g}';f.text(1.45,y,pp,va='center',fontsize=8.5)
f.axvline(0,c='#888888',ls='--',lw=1);f.set_yticks([8,4,0],genes,fontstyle='italic');f.set_ylim(-1.5,9.5);f.set_xlim(-3.65,2.2);f.text(1.45,9.25,'P value',fontsize=8.5)
f.set_xlabel('MR effect on log-odds (β, 95% CI)');f.set_title('(B) Mendelian randomization',loc='left',fontweight='bold',fontsize=12);f.spines['left'].set_visible(False);f.tick_params(axis='y',length=0)
handles,labels=a.get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,-.005));fig.tight_layout(rect=[0,.07,1,1],w_pad=2);save(fig,'Figure3')
fig,axes=plt.subplots(1,3,figsize=(7.2,3.4))
for ax,g in zip(axes,genes):
    r=ti[ti.gene==g].iloc[0];ted=np.array([float(v) for v in r.ted_per_sample.split(';')])
    ax.scatter([0],[r.ctrl_tpm],s=48,c='#202020',marker='s',zorder=3);ax.scatter(1+np.linspace(-.1,.1,4),ted,s=40,c='#a95038',zorder=3)
    ax.plot([.78,1.22],[r.ted_tpm]*2,c='#a95038',lw=1.5)
    ax.set_xlim(-.45,1.45);ax.set_ylim(bottom=0);ax.set_xticks([0,1],['Control\n(n = 1)','TED\n(n = 4)']);ax.set_title(g,fontstyle='italic',pad=24);ax.set_ylabel('Transcript abundance (TPM)')
    ax.text(.5,1.01,f'log₂ fold change = {r.log2fc_tpm:+.2f}',ha='center',transform=ax.transAxes,fontsize=8)
fig.tight_layout(w_pad=1.8);save(fig,'FigureS1')
fig,ax=plt.subplots(figsize=(7.0,7.4));h4plot(ax,candidates);ax.legend(frameon=False,loc='lower center',bbox_to_anchor=(.5,-.17),ncol=3,fontsize=9);fig.tight_layout(rect=[0,.04,1,1]);save(fig,'FigureS2')
print('Built 5 figures, each PNG at 300 dpi and vector PDF')
