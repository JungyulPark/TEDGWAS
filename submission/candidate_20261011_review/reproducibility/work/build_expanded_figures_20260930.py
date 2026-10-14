from pathlib import Path
import sys,json,hashlib
W=Path(__file__).resolve().parent;sys.path.insert(0,str(W/'pydeps'))
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
from matplotlib.ticker import NullLocator,NullFormatter
R=W/'expanded_analysis_20260930/results';O=W.parent/'outputs/TED_TwoStage_Analysis_20260930';F=O/'Figures';F.mkdir(exist_ok=True)
co=pd.read_csv(R/'Stage2_all_504_colocalization_settings.csv');co=co[co.scenario.eq('original_reference')];mr=pd.read_csv(R/'Stage1_13_discovery_genes.csv')
OC=['BBJ_Graves','UKB_hyperthyroid','FinnGen_GO'];labels=['BBJ GD','UKB hyperthyroidism','FinnGen ophthalmopathy'];GC=['#0072B2','#D55E00','#6B5B95'];genes=['TSHR','IGF1R','CTLA4'];mhc=['HLA-A','HLA-DQA2','C4A','TUBB','PSMB8'];allg=genes+['TNFSF14','IFNGR1','MAPKAPK5','HSD3B7','VKORC1','PRSS36']+mhc
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black'})
source={'matrix':[],'prior':[],'MHC':[]};checks={}
def rec(g,o,p=1e-5):return co[(co.gene==g)&(co.outcome==o)&co.p12.eq(p)].iloc[0]
def fmt(v):return '<0.001' if v<.001 else f'{v:.3f}'
def save(fig,name):
 for t in fig.findobj(Text):t.set_color('black')
 fig.canvas.draw();texts=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
 checks[name]={'all_text_black':all(matplotlib.colors.to_hex(t.get_color())=='#000000' for t in texts),'text_items':len(texts)}
 for ext in ['png','pdf']:fig.savefig(F/(name+'.'+ext),dpi=300,facecolor='white')
def heatmap(ax):
 cmap=LinearSegmentedColormap.from_list('shared',['#F7FAFA','#D3E9E5','#8EC7BD']);ys=np.arange(14,dtype=float);ys[3:]+=.25;ys[9:]+=.35
 ax.set(xlim=(-.51,2.51),ylim=(ys[-1]+.51,-.51));ax.set_xticks(range(3),['BBJ GD','UKB\nhyperthyroidism','FinnGen\nophthalmopathy']);ax.set_yticks(ys,allg,fontstyle='italic');ax.tick_params(length=0,pad=7)
 for spine in ax.spines.values():spine.set_visible(False)
 texts=[];strong=0
 for i,g in enumerate(allg):
  for k,o in enumerate(OC):
   v=rec(g,o)['PP.H4'];ax.add_patch(Rectangle((k-.5,ys[i]-.5),1,1,facecolor=cmap(v),edgecolor='white',lw=.5))
   if v>=.8:ax.add_patch(Rectangle((k-.5,ys[i]-.5),1,1,fill=False,edgecolor='black',lw=1.2));strong+=1
   t=ax.text(k,ys[i],fmt(v),ha='center',va='center',fontsize=9);texts.append((t,k,ys[i]));source['matrix'].append({'gene':g,'outcome':o,'H4':float(v)})
 ax.figure.canvas.draw();renderer=ax.figure.canvas.get_renderer()
 for t,k,y in texts:
  b=t.get_window_extent(renderer);xy=ax.transData.transform([[k-.5,y-.5],[k+.5,y+.5]]);lo=xy.min(axis=0);hi=xy.max(axis=0);assert b.x0>lo[0]+2 and b.x1<hi[0]-2 and b.y0>lo[1]+2 and b.y1<hi[1]-2
 assert strong==6
 return ax
fig=plt.figure(figsize=(12.4,8.6));heat=fig.add_axes([.12,.18,.31,.71]);heatmap(heat)
fig.text(.03,.965,'A',fontsize=16,fontweight='bold');fig.text(.07,.965,'Shared-variant support: all 14 genes',fontsize=12,fontweight='bold')
fig.text(.48,.965,'B',fontsize=16,fontweight='bold');fig.text(.52,.965,'Prior sensitivity: selected genes',fontsize=12,fontweight='bold')
for i,g in enumerate(genes):
 ax=fig.add_axes([.51+i*.16,.37,.138,.44])
 for k,o in enumerate(OC):
  vals=[float(rec(g,o,p)['PP.H4']) for p in [1e-5,5e-6,1e-6]];ax.plot([0,1,2],vals,color=GC[k],marker=['o','s','^'][k],lw=1.4,ms=4)
  source['prior'] += [dict(gene=g,outcome=o,p12=p,H4=v) for p,v in zip([1e-5,5e-6,1e-6],vals)]
 ax.axhline(.8,ls='--',color='black',lw=.8);ax.set(xlim=(-.18,2.18),ylim=(-.02,1.06));ax.set_xticks([0,1,2],['D','I','C']);ax.set_yticks(np.arange(0,1.01,.2));ax.set_title(g,fontstyle='italic',fontweight='bold',pad=12);ax.grid(axis='y',color='#dfe5e8',lw=.5)
 if i:ax.set_yticklabels([])
 else:ax.set_ylabel('Shared-variant posterior probability');ax.annotate('0.661',(2,rec('TSHR',OC[0],1e-6)['PP.H4']),xytext=(-8,-22),textcoords='offset points',ha='right',fontsize=9)
fig.legend(handles=[Line2D([],[],color=c,marker=m,label=l) for c,m,l in zip(GC,['o','s','^'],labels)],loc='lower left',bbox_to_anchor=(.50,.18),frameon=False,fontsize=9)
fig.text(.51,.31,'D: p12 = 10⁻⁵     I: 5×10⁻⁶     C: 10⁻⁶',fontsize=9)
fig.text(.12,.065,'Black border: PP.H4 ≥ 0.80 (not a P value)\nCombined criterion: BBJ AND FinnGen',fontsize=9,linespacing=1.7)
fig.text(.51,.065,'MHC genes added after review.\nSingle-causal-variant model; multiple signals unresolved.',fontsize=9,linespacing=1.7)
save(fig,'Figure3_Expanded_14Genes')
# A legible single-panel replacement for the old nine-gene presentation crop.
f2=plt.figure(figsize=(6.2,8.4));heatmap(f2.add_axes([.27,.17,.67,.75]));f2.text(.05,.96,'Colocalization: all 14 evaluated genes',fontsize=13,fontweight='bold');f2.text(.27,.055,'Black border: PP.H4 ≥ 0.80\nPrimary prior p12 = 10⁻⁵; MHC extension is post hoc.',fontsize=9,linespacing=1.6);save(f2,'Figure3A_Presentation_14Genes');plt.close(f2);plt.close(fig)
fig=plt.figure(figsize=(13.4,8.1));left=fig.add_axes([.19,.16,.22,.72]);nums=fig.add_axes([.435,.16,.285,.72]);bar=fig.add_axes([.755,.16,.17,.72]);nums.axis('off')
fig.text(.035,.955,'A',fontsize=16,fontweight='bold');fig.text(.073,.955,'MHC genes: MR associations',fontsize=12,fontweight='bold');fig.text(.73,.955,'B',fontsize=16,fontweight='bold');fig.text(.77,.955,'Colocalization',fontsize=12,fontweight='bold')
ys=[];ylabels=[]
for i,g in enumerate(mhc):
 for k,o in enumerate(OC):
  y=i*3.6+k;ys.append(y);ylabels.append(g+' / '+['BBJ','UKB','FinnGen'][k]);r=mr[(mr.gene_symbol==g)&(mr.outcome==o)].iloc[0];c=rec(g,o)
  left.errorbar(r['or'],y,xerr=[[r['or']-r.or_ci_lower],[r.or_ci_upper-r['or']]],fmt=['o','s','^'][k],color=GC[k],ms=4.5,capsize=2,lw=1.1)
  if r.pvalue<.001:
   mant,exp=f'{r.pvalue:.2e}'.split('e');pv=mant+'×10'+str(int(exp)).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
  else:pv=f'{r.pvalue:.3g}'
  nums.text(0,y,f'{r["or"]:.2f} ({r.or_ci_lower:.2f}–{r.or_ci_upper:.2f})',va='center',fontsize=9)
  nums.text(.66,y,pv,va='center',fontsize=9,fontweight='bold' if r.pvalue<.05 else 'normal')
  xpos=0
  for h,color in [('PP.H2','#8DA8BE'),('PP.H3','#E6B966'),('PP.H4','#258E83')]:bar.barh(y,c[h],left=xpos,height=.68,color=color,edgecolor='white',lw=.15);xpos+=c[h]
  bar.text(1.05,y,fmt(c['PP.H4']),va='center',ha='left',fontsize=9,clip_on=False)
  source['MHC'].append({'gene':g,'outcome':o,'or':float(r['or']),'lower':float(r.or_ci_lower),'upper':float(r.or_ci_upper),'pvalue':float(r.pvalue),'H2':float(c['PP.H2']),'H3':float(c['PP.H3']),'H4':float(c['PP.H4'])})
for ax in [left,nums,bar]:ax.set_ylim(ys[-1]+.65,-.9)
left.set_xscale('log');left.set_xlim(.2,8);left.set_xticks([.25,.5,1,2,4,8],[.25,.5,1,2,4,8]);left.set_yticks(ys,ylabels,fontsize=9);left.axvline(1,color='black',ls='--',lw=.8);left.set_xlabel('Odds ratio (log scale)');left.grid(axis='x',color='#e6e6e6',lw=.5)
left.xaxis.set_minor_locator(NullLocator());left.xaxis.set_minor_formatter(NullFormatter())
nums.set_xlim(0,1);nums.text(0,1.025,'OR (95% CI)',transform=nums.transAxes,fontweight='bold');nums.text(.66,1.025,'P value',transform=nums.transAxes,fontweight='bold')
bar.set_xlim(0,1);bar.set_yticks([]);bar.set_xticks([0,.5,1]);bar.set_xlabel('Posterior probability');bar.text(1.05,1.025,'H4',transform=bar.transAxes,fontweight='bold')
fig.legend(handles=[Patch(fc=c,label=l) for c,l in [('#8DA8BE','H2: expression only'),('#E6B966','H3: distinct variants'),('#258E83','H4: shared variant')]],loc='lower center',bbox_to_anchor=(.59,.057),frameon=False,ncol=3,fontsize=9)
fig.text(.035,.025,'Bold P < 0.05 is nominal. All five genes met the corrected BBJ discovery threshold. H3 dominance is model-dependent, not proof of no gene effect.',fontsize=9)
save(fig,'FigureS3_MHC_MR_and_Colocalization');plt.close(fig)
(F/'figure_source_values.json').write_text(json.dumps(source,indent=2),encoding='utf8');(F/'figure_layout_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf8')
print('Created three figure variants; all labels black, 42 heatmap values and 15 MHC rows from verified results.')
