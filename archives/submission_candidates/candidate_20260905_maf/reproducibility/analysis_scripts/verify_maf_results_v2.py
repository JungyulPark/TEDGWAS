"""Independent primary MR formula checks and quantitative sensitivity summary."""
from pathlib import Path
import sys,json,hashlib
W=Path(__file__).resolve().parent;sys.path.insert(0,str(W/'pydeps'))
import pandas as pd
import numpy as np
from scipy.stats import norm
T=W/'maf_sensitivity_20260905';S=T/'sensitivity_v1';C=T/'coloc_sensitivity_v1'
primary=pd.read_csv(S/'primary_MR_all_scenarios.csv')
base=pd.read_csv(T/'baseline_v1/primary_MR_reproduced.csv')
h=pd.read_csv(S/'harmonized_eqtlgen_validation.csv')
h=h[h.mr_keep].copy()
rows=[]
for (g,oc),v in h.groupby(['id.exposure','id.outcome']):
    bx=v['beta.exposure'].to_numpy();by=v['beta.outcome'].to_numpy();sy=v['se.outcome'].to_numpy();n=len(v)
    if n==1:b=by[0]/bx[0];se=abs(sy[0]/bx[0])
    else:
        ss=np.sum(bx*bx/sy**2);b=np.sum(bx*by/sy**2)/ss
        q=np.sum((by-b*bx)**2/sy**2);se=np.sqrt(max(1,q/(n-1))/ss)
    p=2*norm.sf(abs(b/se))
    rows.append(dict(gene_symbol=g,outcome=oc,n_iv=n,beta=b,se=se,pvalue=p))
ind=pd.DataFrame(rows)
new=primary[primary.scenario.eq('reharmonized_eqtlgen')&primary.n_iv.gt(0)].copy()
cmp=ind.merge(new,on=['gene_symbol','outcome'],validate='one_to_one',suffixes=('_independent','_R'))
assert len(cmp)==len(new)==len(ind)
assert (cmp.n_iv_independent==cmp.n_iv_R).all()
deltas={c:float((cmp[f'{c}_independent']-cmp[f'{c}_R']).abs().max()) for c in ('beta','se','pvalue')}
assert max(deltas.values())<1e-8,deltas
cmp.to_csv(S/'independent_MR_validation_v2.csv',index=False)
ex=pd.read_csv(S/'exposure_frequency_comparison.csv')
old_h=pd.read_csv(T/'baseline_v1/harmonized_all.csv')
new_h=pd.read_csv(S/'harmonized_eqtlgen_validation.csv')
# A rsID can label multiple outcome allele pairs. Match records one-to-one;
# rsID-only joins create spurious changes by pairing retained and rejected alleles.
record_keys=['id.exposure','id.outcome','SNP','outcome_allele_pair','se.outcome','abs_beta_outcome']
for table in (old_h,new_h):
    table['outcome_allele_pair']=['/'.join(sorted([str(a),str(b)])) for a,b in zip(table['effect_allele.outcome'],table['other_allele.outcome'])]
    table['abs_beta_outcome']=table['beta.outcome'].abs()
    assert not table.duplicated(record_keys).any()
keepcols=record_keys+['mr_keep','palindromic','ambiguous','beta.outcome']
qc=old_h[keepcols].merge(new_h[keepcols],on=record_keys,how='outer',validate='one_to_one',suffixes=('_ref','_new'))
qc=qc.rename(columns={'id.exposure':'gene_symbol','id.outcome':'outcome','mr_keep_ref':'keep_ref','mr_keep_new':'keep_new','beta.outcome_ref':'beta_out_ref','beta.outcome_new':'beta_out_new'})
qc['keep_ref']=qc.keep_ref.fillna(False).astype(bool);qc['keep_new']=qc.keep_new.fillna(False).astype(bool)
qc['kept_changed']=qc.keep_ref!=qc.keep_new
qc['outcome_sign_changed']=qc.keep_ref&qc.keep_new&(np.sign(qc.beta_out_ref)!=np.sign(qc.beta_out_new))
qc[qc.kept_changed|qc.outcome_sign_changed].to_csv(S/'changed_variant_decisions_v2.csv',index=False)
changes=[]
for oc,v in qc.groupby('outcome'):
    changes.append(dict(outcome=oc,added=int((~v.keep_ref&v.keep_new).sum()),removed=int((v.keep_ref&~v.keep_new).sum()),retained_outcome_sign_changed=int(v.outcome_sign_changed.sum()),affected_gene_outcome_pairs=int(v.loc[v.kept_changed|v.outcome_sign_changed,'gene_symbol'].nunique())))
count_changes=pd.DataFrame(changes);count_changes.to_csv(S/'harmonization_changes_summary_v2.csv',index=False)
qcref=qc[qc.kept_changed].copy()
qcref['reason']=np.where(qcref.SNP.isin(ex.loc[ex.eaf_ref.isna(),'snp']),'reference_AF_previously_missing',np.where(qcref.palindromic_ref.fillna(False)|qcref.palindromic_new.fillna(False),'palindromic_AF_decision','other_requires_review'))
qcref.to_csv(S/'changed_variant_decisions_classified_v2.csv',index=False)
screen={}
for scenario,data in [('original_reference',base)]+[(s,primary[primary.scenario.eq(s)&primary.n_iv.gt(0)]) for s in primary.scenario.unique()]:
    sig=data[data.outcome.eq('BBJ_Graves')&data.pvalue.lt(.05/2544)]
    screen[scenario]={'estimable_genes_by_outcome':data.groupby('outcome').size().to_dict(),'BBJ_discovery_hits':sorted(sig.gene_symbol.tolist())}
baseline_hits=set(screen['original_reference']['BBJ_discovery_hits'])
for scenario,v in screen.items():
    v['added_discovery_hits']=sorted(set(v['BBJ_discovery_hits'])-baseline_hits)
    v['lost_discovery_hits']=sorted(baseline_hits-set(v['BBJ_discovery_hits']))
coloc=pd.read_csv(C/'coloc_all_scenarios.csv')
conclusions={}
for scenario in coloc.scenario.unique():
    x=coloc[coloc.scenario.eq(scenario)&coloc.p12.eq(1e-5)]
    strong=x[x['PP.H4']>=.8]
    both=set(strong.loc[strong.outcome.eq('BBJ_Graves'),'gene'])&set(strong.loc[strong.outcome.eq('FinnGen_GO'),'gene'])
    conclusions[scenario]={'BBJ_and_FinnGen_H4_ge_0_8':sorted(both)}
power=pd.read_csv(S/'power_all_scenarios.csv')
expected=pd.read_csv(W/'revision/submission/candidate_20260905/provenance/power_summary_verified.csv')
pm=power[power.scenario.eq('original_reference')].merge(expected,on='outcome',validate='one_to_one',suffixes=('_rerun','_canonical'))
assert (pm.n_genes_rerun==pm.n_genes_canonical).all()
powerdiff=max((pm[f'{k}_rerun']-pm[f'{k}_canonical']).abs().max() for k in ('or_q1','or_median','or_q3','frac_OR1_5','frac_OR2','frac_OR3'))
assert powerdiff<1e-8
summary=dict(status='PASS',eligible_genes=2544,instrument_rows=6135,instruments_with_AF=int(ex.eligible_af.sum()),max_abs_delta_Z=float(ex.delta_z.abs().max()),independent_MR_rows=len(cmp),independent_MR_max_abs_differences=deltas,baseline_power_max_abs_difference=float(powerdiff),harmonization_changes=changes,screen=screen,colocalization_conclusions=conclusions,unexpected_variant_decisions=int((qcref.reason=='other_requires_review').sum()),coloc_summary=json.loads((C/'coloc_summary.json').read_text()))
assert summary['max_abs_delta_Z']<1e-8
assert summary['unexpected_variant_decisions']==0,'Unexplained variant decisions need review'
(T/'sensitivity_verification_v2.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
