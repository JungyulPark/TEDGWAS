from pathlib import Path
import json,sys,shutil,datetime,hashlib
T=Path(__file__).resolve().parent;W=T.parent;ROOT=W.parent
sys.path.insert(0,str(W/'pydeps'))
import numpy as np,pandas as pd
from scipy.stats import norm,chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

O=ROOT/'outputs/TED_Reanalysis_20261002';O.mkdir(exist_ok=True)
TAB=O/'Tables';TAB.mkdir(exist_ok=True);FIG=O/'Figures';FIG.mkdir(exist_ok=True)
assert (T/'coloc_v1/R_coloc_verification.csv').exists()
mr=pd.read_csv(T/'mr_v1/MR_all_genes_all_scenarios.csv')
old=pd.read_csv(W/'expanded_analysis_20260930/results/Stage1_all_2544_genes_3_outcomes.csv')
old=old.rename(columns={'or':'OR','or_ci_lower':'CI_lower','or_ci_upper':'CI_upper'})
new=mr[mr.outcome.ne('FinnGen_GO')&mr.scenario.eq('original_reference')].copy()
columns=['gene_symbol','outcome','method','n_iv','beta','se','pvalue','OR','CI_lower','CI_upper']
full=pd.concat([old[columns],new[columns]],ignore_index=True)
assert len(full)==2544*5 and not full.duplicated(['gene_symbol','outcome']).any()
full['estimable']=full.pvalue.notna();full['nominal_p_lt_0_05']=full.pvalue.lt(.05);full['P_lt_0_05_over_2544']=full.pvalue.lt(.05/2544)
hits=old[old.outcome.eq('BBJ_Graves')&old.pvalue.lt(.05/2544)].sort_values('pvalue').gene_symbol.tolist();focus=hits+['IGF1R']
assert len(hits)==13
col_old=pd.read_csv(W/'expanded_analysis_20260930/results/Stage2_all_504_colocalization_settings.csv')
col_new=pd.read_csv(T/'coloc_v1/Coloc_new_336_settings.csv')
col=pd.concat([col_old,col_new],ignore_index=True);assert len(col)==840
cp=col[col.scenario.eq('original_reference')&col.p12.eq(1e-5)].drop(columns=['scenario','p12']).rename(columns={'gene':'gene_symbol'})
focusdf=full[full.gene_symbol.isin(focus)].merge(cp,on=['gene_symbol','outcome'],validate='one_to_one')
low=col[col.scenario.eq('original_reference')&col.p12.eq(1e-6)][['gene','outcome','PP.H4']].rename(columns={'gene':'gene_symbol','PP.H4':'H4_low_prior'})
af=col[col.scenario.eq('available_eqtlgen')&col.p12.eq(1e-5)][['gene','outcome','PP.H4']].rename(columns={'gene':'gene_symbol','PP.H4':'H4_eqtlgen'})
focusdf=focusdf.merge(low,on=['gene_symbol','outcome']).merge(af,on=['gene_symbol','outcome'])
focusdf['BBJ_discovery']=focusdf.gene_symbol.isin(hits)
bb=full[full.outcome.eq('BBJ_Graves')][['gene_symbol','beta']].rename(columns={'beta':'BBJ_beta'})
focusdf=focusdf.merge(bb,on='gene_symbol');focusdf['same_direction_BBJ']=np.sign(focusdf.beta)==np.sign(focusdf.BBJ_beta)
focusdf['focused_GD_corroboration']=focusdf.outcome.eq('E4_GRAVES_STRICT')&focusdf.BBJ_discovery&focusdf.same_direction_BBJ&focusdf.pvalue.lt(.05/13)
order=['BBJ_Graves','E4_GRAVES_STRICT','UKB_hyperthyroid','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
labels=['BBJ GD','FinnGen GD strict','UKB hyperthyroidism','FinnGen GO','FinnGen GO strict']
focusdf['gene_order']=focusdf.gene_symbol.map({g:i for i,g in enumerate(focus)});focusdf['outcome_order']=focusdf.outcome.map(dict(zip(order,range(5))))
focusdf=focusdf.sort_values(['gene_order','outcome_order']).drop(columns=['gene_order','outcome_order'])
full.to_csv(TAB/'All_2544_genes_five_outcomes.csv',index=False)
focusdf.to_csv(TAB/'Table2_14_genes_MR_and_colocalization.csv',index=False)
col.to_csv(TAB/'Colocalization_840_settings.csv',index=False)
mr[mr.outcome.ne('FinnGen_GO')].to_csv(TAB/'New_MR_all_frequency_scenarios.csv',index=False)
shutil.copy2(T/'mr_v1/leave_one_out_14genes.csv',TAB/'Leave_one_out_new_outcomes.csv')
shutil.copy2(T/'mr_v1/harmonization_counts.csv',TAB/'MR_instrument_coverage.csv')
shutil.copy2(T/'coloc_v1/coverage.csv',TAB/'Colocalization_variant_coverage.csv')
coh=pd.DataFrame([
 ['BBJ_Graves','BBJ GD','GCST90018627',2809,172656,'Original discovery; East Asian'],
 ['E4_GRAVES_STRICT','FinnGen GD strict','R12 E4_GRAVES_STRICT',3962,496386,'New closer-phenotype comparison; Finnish'],
 ['UKB_hyperthyroid','UKB hyperthyroidism','GCST90038636',3731,480867,'Broader phenotype; approximate linear-to-log-OR conversion'],
 ['FinnGen_GO','FinnGen GO','R12 GRAVES_OPHT',858,499490,'Original ophthalmopathy vs population noncases'],
 ['E4_GRAVES_OPHT_STRICT','FinnGen GO strict','R12 E4_GRAVES_OPHT_STRICT',753,499595,'New phenotype sensitivity; overlapping FinnGen participants']],columns=['outcome','label','accession','cases','controls','role'])
coh.to_csv(TAB/'Table1_Datasets.csv',index=False)
counts=full.groupby('outcome').agg(eligible=('gene_symbol','size'),estimated=('estimable','sum'),nominal=('nominal_p_lt_0_05','sum'),screen_bonferroni=('P_lt_0_05_over_2544','sum')).reindex(order).reset_index()
counts['not_estimable']=counts.eligible-counts.estimated;counts.to_csv(TAB/'Screen_counts.csv',index=False)
confirm=focusdf[focusdf.outcome.eq('E4_GRAVES_STRICT')].copy();confirm.to_csv(TAB/'Table3_BBJ_discovery_followup_in_FinnGen_GD.csv',index=False)
newhits=new[new.pvalue.lt(.05/2544)].copy();newhits['previous_BBJ_discovery']=newhits.gene_symbol.isin(hits);newhits['coloc_evaluated']=newhits.gene_symbol.isin(focus)
newhits.sort_values(['outcome','pvalue']).to_csv(TAB/'New_endpoint_exploratory_screen_hits.csv',index=False)
paired=mr[mr.scenario.eq('paired_reference')].merge(mr[mr.scenario.eq('paired_eqtlgen')],on=['gene_symbol','outcome'],suffixes=('_ref','_new'))
reh=mr[mr.scenario.eq('original_reference')&mr.outcome.ne('FinnGen_GO')].merge(mr[mr.scenario.eq('reharmonized_eqtlgen')],on=['gene_symbol','outcome'],suffixes=('_ref','_new'))
freqstats=[]
for oc in order[1::3]: # overwritten explicitly below for clarity
    pass
for oc in ['E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']:
    x=reh[reh.outcome.eq(oc)];p=paired[paired.outcome.eq(oc)]
    c=col_new[col_new.outcome.eq(oc)&col_new.p12.eq(1e-5)]
    a=c[c.scenario.eq('paired_reference')].merge(c[c.scenario.eq('paired_eqtlgen')],on=['gene','outcome'],suffixes=('_ref','_new'))
    b=c[c.scenario.eq('original_reference')].merge(c[c.scenario.eq('available_eqtlgen')],on=['gene','outcome'],suffixes=('_ref','_new'))
    freqstats.append(dict(outcome=oc,MR_changed_instrument_count=int(x.n_iv_ref.ne(x.n_iv_new).sum()),MR_changed_nominal_status=int(x.pvalue_ref.lt(.05).ne(x.pvalue_new.lt(.05)).sum()),MR_changed_full_screen_status=int(x.pvalue_ref.lt(.05/2544).ne(x.pvalue_new.lt(.05/2544)).sum()),coloc_max_paired_delta_H4=float((a['PP.H4_new']-a['PP.H4_ref']).abs().max()),coloc_max_total_delta_H4=float((b['PP.H4_new']-b['PP.H4_ref']).abs().max()),coloc_threshold_crossings=b.loc[b['PP.H4_new'].ge(.8).ne(b['PP.H4_ref'].ge(.8)),'gene'].tolist()))
pd.DataFrame(freqstats).to_csv(TAB/'Frequency_sensitivity_summary.csv',index=False)

# Independent check of saved SNP-level MR results (Python), including the baseline identity control.
h=pd.read_csv(T/'mr_v1/harmonized_snps.csv');h=h[h.mr_keep]
checks=[]
for (g,oc,sc),x in h.groupby(['exposure','outcome','scenario']):
    bx=x['beta.exposure'].to_numpy();by=x['beta.outcome'].to_numpy();sy=x['se.outcome'].to_numpy();n=len(x)
    assert not x.SNP.duplicated().any()
    b=np.sum(bx*by/sy**2)/np.sum(bx**2/sy**2);q=np.sum((by-b*bx)**2/sy**2)
    se=np.sqrt((max(1,q/(n-1)) if n>1 else 1)/np.sum(bx**2/sy**2));p=2*norm.sf(abs(b/se))
    r=mr[mr.gene_symbol.eq(g)&mr.outcome.eq(oc)&mr.scenario.eq(sc)].iloc[0]
    checks.append(max(abs(b-r.beta),abs(se-r.se),abs(p-r.pvalue)));assert n==r.n_iv
assert max(checks)<1e-8

# Shareable scientific plots, generated directly from the exported numeric table.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','pdf.fonttype':42,'ps.fonttype':42})
matrix=focusdf.pivot(index='gene_symbol',columns='outcome',values='PP.H4').reindex(index=focus,columns=order)
fig,ax=plt.subplots(figsize=(10,7.5));im=ax.imshow(matrix.to_numpy(),vmin=0,vmax=1,cmap='Blues')
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
        ax.text(.03,i+.22,f'OR {r.OR:.2f} ({r.CI_lower:.2f}-{r.CI_upper:.2f}); P={r.pvalue:.2g}',transform=ax.get_yaxis_transform(),fontsize=8,color='black')
    ax.axvline(1,color='#777777',ls='--',lw=1);ax.set_xscale('log');ax.set_title(g,fontstyle='italic',weight='bold');ax.set_xlabel('Odds ratio (95% CI), log scale');ax.grid(axis='x',alpha=.18)
axs[0].set_yticks(range(5),['BBJ GD','FinnGen GD strict [new]','UKB hyperthyroidism','FinnGen GO','FinnGen GO strict [new]']);axs[0].invert_yaxis()
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
summary={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'counts':counts.to_dict('records'),'focus':json.loads(focusdf.to_json(orient='records',double_precision=15)),'corroborated_GD':confirm[confirm.focused_GD_corroboration].gene_symbol.tolist(),'new_GD_H4_ge_0_8':confirm[confirm['PP.H4'].ge(.8)].gene_symbol.tolist(),'frequency_sensitivity':freqstats,'new_exploratory_hits':json.loads(newhits.to_json(orient='records',double_precision=15)),'independent_python_MR_checks':len(checks),'independent_python_MR_max_difference':float(max(checks))}
(O/'analysis_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8')
print(json.dumps({k:v for k,v in summary.items() if k not in ['focus','new_exploratory_hits']},indent=2,ensure_ascii=False))
print(focusdf[focusdf.gene_symbol.isin(['TSHR','IGF1R','CTLA4'])][['gene_symbol','outcome','n_iv','OR','CI_lower','CI_upper','pvalue','PP.H4','H4_low_prior','H4_eqtlgen']].to_string(index=False))
