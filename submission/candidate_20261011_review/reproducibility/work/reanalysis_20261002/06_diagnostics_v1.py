from pathlib import Path
import sys,json,datetime
T=Path(__file__).resolve().parent;W=T.parent;O=W.parent/'outputs/TED_Reanalysis_20261002'
sys.path.insert(0,str(W/'pydeps'))
import numpy as np,pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
log={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'Post hoc diagnostic after new TSHR GD H4=0.264: inspect retained regional signal patterns and exact sensitivity changes. No model/threshold reselection, no SuSiE, no new causal claims.'}
(T/'diagnostic_plan.json').write_text(json.dumps(log,indent=2),encoding='utf8')
mr=pd.read_csv(T/'mr_v1/MR_all_genes_all_scenarios.csv')
base=mr[mr.outcome.ne('FinnGen_GO')&mr.scenario.eq('original_reference')]
sens=mr[mr.scenario.eq('reharmonized_eqtlgen')]
m=base.merge(sens,on=['gene_symbol','outcome'],suffixes=('_ref','_new'))
m['changed_screen_significance']=m.pvalue_ref.lt(.05/2544).ne(m.pvalue_new.lt(.05/2544))
m['changed_nominal_significance']=m.pvalue_ref.lt(.05).ne(m.pvalue_new.lt(.05))
m['changed_instrument_count']=m.n_iv_ref.ne(m.n_iv_new)
m[m.changed_screen_significance|m.changed_nominal_significance|m.changed_instrument_count].to_csv(O/'Tables/MR_frequency_changed_results.csv',index=False)
print('Full-screen frequency changes:');print(m[m.changed_screen_significance][['gene_symbol','outcome','n_iv_ref','n_iv_new','pvalue_ref','pvalue_new']].to_string(index=False))
print('GD 13-gene confirmation changes:');print(m[m.focused_GD_corroboration_ref.ne(m.focused_GD_corroboration_new)][['gene_symbol','outcome','pvalue_ref','pvalue_new']].to_string(index=False))
loo=pd.read_csv(T/'mr_v1/leave_one_out_14genes.csv')
print('Anchor leave-one-out:');print(loo[loo.gene_symbol.isin(['IGF1R','CTLA4'])][['gene_symbol','outcome','removed_snp','n_iv','OR','pvalue']].to_string(index=False))
eq=pd.read_csv(T/'inputs_v2/eqtl_14genes.tsv',sep='\t');eq=eq[eq.GeneSymbol.eq('TSHR')].copy()
pos=eq.set_index('SNP').SNPPos;leadpos=int(pos.loc['rs179252']);stats=[];patterns=[]
order=['BBJ_Graves','E4_GRAVES_STRICT','UKB_hyperthyroid','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
labels=['BBJ GD','FinnGen GD strict [new]','UKB hyperthyroidism','FinnGen GO','FinnGen GO strict [new]']
def logp(z):return -(np.log(2)+norm.logsf(np.abs(z)))/np.log(10)
for oc in order:
    p=T/'coloc_v1/r_validation'/f'TSHR__{oc}__original_reference.csv'
    if not p.exists():p=W/'expanded_analysis_20260930/inputs/r_validation'/f'TSHR__{oc}__original_reference.csv'
    d=pd.read_csv(p);d['pos_hg19']=d.snp.map(pos);assert d.pos_hg19.notna().all()
    d['genome_build']='GRCh37';d['minus_log10_p_disease_normal']=logp(d.beta1/np.sqrt(d.varbeta1));d['minus_log10_p_expression_normal']=logp(d.beta2/np.sqrt(d.varbeta2));d['outcome']=oc
    lead=d.loc[d.minus_log10_p_disease_normal.idxmax()];iv=d[d.snp.eq('rs179252')].iloc[0]
    for label,r in [('MR_instrument',iv),('strongest_disease_Z_in_analyzed_set',lead)]:
        stats.append({'outcome':oc,'role':label,'snp':r.snp,'pos_hg19':int(r.pos_hg19),'genome_build':'GRCh37','disease_abs_Z':float(abs(r.beta1/np.sqrt(r.varbeta1))),'disease_normal_P':float(2*norm.sf(abs(r.beta1/np.sqrt(r.varbeta1)))),'expression_abs_Z':float(abs(r.beta2/np.sqrt(r.varbeta2))),'expression_normal_P':float(2*norm.sf(abs(r.beta2/np.sqrt(r.varbeta2))))})
    patterns.append(d)
pd.DataFrame(stats).to_csv(O/'Tables/TSHR_regional_diagnostic_lead_variants.csv',index=False)
pd.concat(patterns).to_csv(O/'Tables/TSHR_regional_diagnostic_all_variants.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','pdf.fonttype':42})
fig,axs=plt.subplots(6,1,figsize=(11,10),sharex=True)
axs[0].scatter(eq.SNPPos/1e6,logp(eq.Zscore),s=5,color='#555555',alpha=.45);axs[0].set_title('TSHR regional patterns: gene expression and disease association',weight='bold',fontsize=13,pad=18);axs[0].text(.01,.85,'eQTLGen blood expression',transform=axs[0].transAxes,weight='bold')
for ax,d,label in zip(axs[1:],patterns,labels):
    ax.scatter(d.pos_hg19/1e6,d.minus_log10_p_disease_normal,s=5,color='#43859b',alpha=.5)
    lead=d.loc[d.minus_log10_p_disease_normal.idxmax()];iv=d[d.snp.eq('rs179252')].iloc[0]
    ax.scatter([iv.pos_hg19/1e6],[iv.minus_log10_p_disease_normal],s=32,color='black',marker='D',zorder=4)
    ax.scatter([lead.pos_hg19/1e6],[lead.minus_log10_p_disease_normal],s=45,edgecolors='#9e4215',facecolors='none',marker='s',zorder=4)
    ax.text(.01,.83,label,transform=ax.transAxes,weight='bold')
    ax.text(.99,.83,'Lead: '+lead.snp,transform=ax.transAxes,ha='right',fontsize=8)
for ax in axs:
    ax.axvline(leadpos/1e6,color='black',lw=.7,ls='--');ax.grid(axis='y',alpha=.15);ax.set_ylim(bottom=0);ax.spines[['top','right']].set_visible(False)
axs[-1].set_xlabel('Chromosome 14 position from the eQTL file (GRCh37), Mb')
fig.supylabel('-log10 P (normal approximation from beta / SE)',x=.025)
fig.text(.08,.015,'Diamond / dashed line: MR instrument rs179252. Square: strongest disease Z among analysed variants.\nDifferent lead SNPs do not establish distinct causal variants; LD and multiple signals remain unresolved.',fontsize=9)
fig.tight_layout(rect=[.04,.06,1,1])
for ext in ['png','pdf']:fig.savefig(O/f'Figures/Figure_R4_TSHR_regional_patterns.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
print(pd.DataFrame(stats).to_string(index=False))
