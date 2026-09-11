"""Coloc AF comparison: exact frozen baseline, paired SNPs, all usable AF SNPs.
Uses coloc 5.2.3 beta/varbeta equations, outcome dataset1, eQTL dataset2 (sdY=1).
Never imputes unavailable cohort frequencies. Retains the original p12 grid.
"""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'pydeps'))
import numpy as np
import pandas as pd
from scipy.special import logsumexp

root,afpath,canonpath=map(Path,sys.argv[1:4])
out=root/'coloc_sensitivity_v1'
if out.exists(): raise RuntimeError('Versioned result already exists; do not overwrite')
out.mkdir()
eq=pd.read_csv(root/'inputs/eqtlgen_9genes.tsv',sep='\t')
ref=pd.read_csv(root/'inputs/eur_freq_selected.tsv',sep='\t')
af=pd.read_csv(afpath)
assert not ref.SNP.duplicated().any() and not af.snp.duplicated().any()
assert af.genome_build.eq('GRCh37').all()
e=eq.merge(ref[['SNP','A1','A2','MAF']],on='SNP',how='left',validate='many_to_one')
e=e.merge(af,left_on='SNP',right_on='snp',how='left',validate='many_to_one')
e['eaf_ref']=np.where(e.AssessedAllele==e.A1,e.MAF,np.where(e.AssessedAllele==e.A2,1-e.MAF,np.nan))
ok=((e.AssessedAllele==e.allele_b)&(e.OtherAllele==e.allele_a))|((e.AssessedAllele==e.allele_a)&(e.OtherAllele==e.allele_b))
e['eaf_new']=np.where(ok,np.where(e.AssessedAllele==e.allele_b,e.af_b,1-e.af_b),np.nan)
for kind in ('ref','new'):
    f=e[f'eaf_{kind}'];e[f'valid_{kind}']=f.between(0,1,inclusive='neither')
    e[f'se_{kind}']=np.where(e[f'valid_{kind}'],1/np.sqrt(2*f*(1-f)*(e.NrSamples+e.Zscore**2)),np.nan)
e['genome_build']='GRCh37'
e.to_csv(out/'eqtl_AF_coverage.csv',index=False)
genes=dict(TSHR=14,IGF1R=15,CTLA4=2,TNFSF14=19,IFNGR1=6,MAPKAPK5=12,HSD3B7=16,VKORC1=16,PRSS36=16)
ocs=dict(BBJ_Graves='BBJ',UKB_hyperthyroid='UKB',FinnGen_GO='FinnGen')
rows=[];coverage=[]
def analyze(m,g,oc,scenario,kind):
    if len(m)<50: raise ValueError(f'Too few variants: {g} {oc} {scenario}: {len(m)}')
    z1=m.beta.to_numpy()/m.se.to_numpy();v1=m.se.to_numpy()**2
    z2=m.Zscore.to_numpy();v2=m[f'se_{kind}'].to_numpy()**2
    r1=.2**2/(.2**2+v1);r2=.15**2/(.15**2+v2)
    l1=.5*(np.log1p(-r1)+r1*z1*z1);l2=.5*(np.log1p(-r2)+r2*z2*z2)
    s1=logsumexp(l1);s2=logsumexp(l2);ss=logsumexp(l1+l2)
    h3=s1+s2+np.log1p(-np.exp(ss-s1-s2))
    for p12 in (1e-6,5e-6,1e-5):
        lh=np.array([0,np.log(1e-4)+s1,np.log(1e-4)+s2,2*np.log(1e-4)+h3,np.log(p12)+ss])
        pp=np.exp(lh-logsumexp(lh))
        assert abs(pp.sum()-1)<1e-12
        rows.append(dict(gene=g,outcome=oc,scenario=scenario,p12=p12,n_overlap=len(m),top_snp=m.SNP.iloc[int(np.argmax(l1+l2))],**{f'PP.H{i}':pp[i] for i in range(5)}))
for oc,name in ocs.items():
    d=pd.read_csv(root/f'inputs/{name}_selected.tsv',sep='\t')
    d=d[d.beta.notna()&d.se.gt(0)&d.snp.notna()].copy()
    if oc=='UKB_hyperthyroid':
        mu=3731/484598;d[['beta','se']]=d[['beta','se']]/(mu*(1-mu))
    for g,ch in genes.items():
        m=e[e.GeneSymbol.eq(g)].merge(d[d.chr.eq(ch)],left_on='SNP',right_on='snp',suffixes=('_e','_o'))
        same=(m.AssessedAllele==m.ea)&(m.OtherAllele==m.oa)
        swap=(m.AssessedAllele==m.oa)&(m.OtherAllele==m.ea)
        m=m[same|swap].drop_duplicates('SNP').sort_values('SNP')
        old=m[m.valid_ref];new=m[m.valid_new];common=m[m.valid_ref&m.valid_new]
        coverage.append(dict(gene=g,outcome=oc,n_ref=len(old),n_new=len(new),n_common=len(common),n_ref_missing_AF=len(old)-len(common),n_new_added=len(new)-len(common)))
        for data,scenario,kind in ((old,'original_reference','ref'),(common,'paired_reference','ref'),(common,'paired_eqtlgen','new'),(new,'available_eqtlgen','new')):
            analyze(data,g,oc,scenario,kind)
got=pd.DataFrame(rows)
got.to_csv(out/'coloc_all_scenarios.csv',index=False)
pd.DataFrame(coverage).to_csv(out/'coloc_variant_coverage.csv',index=False)
frozen=pd.read_csv(canonpath)
base=got[got.scenario.eq('original_reference')].merge(frozen,on=['gene','outcome','p12'],validate='one_to_one',suffixes=('_rerun','_canonical'))
base_delta=max((base[f'PP.H{i}_rerun']-base[f'PP.H{i}_canonical']).abs().max() for i in range(5))
assert len(base)==81 and base_delta<1e-9
assert (base.n_overlap_rerun==base.n_overlap_canonical).all() and (base.top_snp_rerun==base.top_snp_canonical).all()
summaries=[]
for sa,sb,label in (('paired_reference','paired_eqtlgen','frequency_only'),('original_reference','available_eqtlgen','frequency_and_coverage')):
    cmp=got[got.scenario.eq(sa)].merge(got[got.scenario.eq(sb)],on=['gene','outcome','p12'],validate='one_to_one',suffixes=('_ref','_new'))
    cmp['delta_H4']=cmp['PP.H4_new']-cmp['PP.H4_ref']
    cmp['crosses_0_8']=(cmp['PP.H4_new']>=.8)!=(cmp['PP.H4_ref']>=.8)
    cmp.to_csv(out/f'coloc_comparison_{label}.csv',index=False)
    summaries.append(dict(comparison=label,max_abs_delta_H4=float(cmp.delta_H4.abs().max()),n_threshold_crossings=int(cmp.crosses_0_8.sum()),n_primary_prior_crossings=int(cmp.loc[cmp.p12.eq(1e-5),'crosses_0_8'].sum())))
summary=dict(baseline_rows=81,max_baseline_posterior_difference=float(base_delta),comparisons=summaries,dataset1='outcome',dataset2='eQTL',eQTL_sdY=1,p1=1e-4,p2=1e-4,p12=[1e-6,5e-6,1e-5])
(out/'coloc_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
