"""Independent reproduction of the frozen coloc 5.2.3 beta/varbeta output.
Uses the published coloc equations with the original extracted input sets.
Never alters the original R script, its inputs, or its 81-row output.
"""
import sys, json, hashlib
from pathlib import Path
W=Path(__file__).resolve().parent
sys.path.insert(0,str(W/'pydeps'))
import numpy as np
import pandas as pd
from scipy.special import logsumexp

BASE=W.parent/'provenance'
IN=W/'inputs'
OUT=W/'rerun_verification';OUT.mkdir(exist_ok=True)
src=BASE/'coloc_canonical_v2.csv'
required=['eqtlgen_9genes.tsv','eur_freq_subset.tsv','BBJ_byrsid.tsv','UKB_byrsid.tsv','FinnGen_byrsid.tsv']
missing=[n for n in required if not (IN/n).is_file()]
if missing:
    sys.exit('Local-only input files are required. See inputs/README.md and inputs_manifest.json. Missing: '+', '.join(missing))
frozen=pd.read_csv(src)
eq=pd.read_csv(IN/'eqtlgen_9genes.tsv',sep='\t')
frq=pd.read_csv(IN/'eur_freq_subset.tsv',sep=r'\s+')
frq=frq[frq.SNP.isin(eq.SNP)]
assert not frq.SNP.duplicated().any()
e=eq.merge(frq[['SNP','A1','A2','MAF']],on='SNP',validate='many_to_one')
e['eaf']=np.where(e.AssessedAllele==e.A1,e.MAF,np.where(e.AssessedAllele==e.A2,1-e.MAF,np.nan))
e=e[e.eaf.between(0,1,inclusive='neither')].copy()
e['se_e']=1/np.sqrt(2*e.eaf*(1-e.eaf)*(e.NrSamples+e.Zscore**2))
genes=dict(TSHR=14,IGF1R=15,CTLA4=2,TNFSF14=19,IFNGR1=6,MAPKAPK5=12,HSD3B7=16,VKORC1=16,PRSS36=16)
ocs=dict(BBJ_Graves='BBJ',UKB_hyperthyroid='UKB',FinnGen_GO='FinnGen')
rows=[]
for oc,filename in ocs.items():
    d=pd.read_csv(IN/f'{filename}_byrsid.tsv',sep='\t')
    d=d[d.beta.notna() & d.se.gt(0) & d.rsid.notna()].copy()
    if oc=='UKB_hyperthyroid':
        mu=3731/484598;d[['beta','se']]=d[['beta','se']]/(mu*(1-mu))
    for g,ch in genes.items():
        m=e[e.GeneSymbol==g].merge(d[d.chr==ch],left_on='SNP',right_on='rsid',suffixes=('_e','_o'))
        same=(m.AssessedAllele==m.ea)&(m.OtherAllele==m.oa)
        swap=(m.AssessedAllele==m.oa)&(m.OtherAllele==m.ea)
        m=m[same|swap].drop_duplicates('SNP').sort_values('SNP')
        assert len(m)>=50
        z1=m.beta.to_numpy()/m.se.to_numpy();v1=m.se.to_numpy()**2
        z2=m.Zscore.to_numpy();v2=m.se_e.to_numpy()**2
        r1=.2**2/(.2**2+v1);r2=.15**2/(.15**2+v2)
        l1=.5*(np.log1p(-r1)+r1*z1*z1)
        l2=.5*(np.log1p(-r2)+r2*z2*z2)
        s1=logsumexp(l1);s2=logsumexp(l2);ss=logsumexp(l1+l2)
        h3=s1+s2+np.log1p(-np.exp(ss-s1-s2))
        for p12 in (1e-6,5e-6,1e-5):
            lh=np.array([0,np.log(1e-4)+s1,np.log(1e-4)+s2,2*np.log(1e-4)+h3,np.log(p12)+ss])
            pp=np.exp(lh-logsumexp(lh))
            rows.append(dict(gene=g,outcome=oc,p12=p12,n_overlap=len(m),top_snp=m.SNP.iloc[int(np.argmax(l1+l2))],**{f'PP.H{i}':pp[i] for i in range(5)}))
got=pd.DataFrame(rows)
cmp=got.merge(frozen,on=['gene','outcome','p12'],suffixes=('_reproduced','_frozen'),validate='one_to_one')
assert len(cmp)==81
delta=max(np.max(np.abs(cmp[f'PP.H{i}_reproduced']-cmp[f'PP.H{i}_frozen'])) for i in range(5))
assert delta<1e-9,delta
assert (cmp.n_overlap_reproduced==cmp.n_overlap_frozen).all()
assert (cmp.top_snp_reproduced==cmp.top_snp_frozen).all()
assert np.max(np.abs(got[[f'PP.H{i}' for i in range(5)]].sum(axis=1)-1))<1e-12
got.to_csv(OUT/'coloc_independent_reproduction.csv',index=False)
summary={'rows':81,'genes':9,'outcomes':3,'priors':3,'max_abs_posterior_difference':float(delta),'all_overlap_counts_identical':True,'all_top_snps_identical':True,'canonical_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'dataset1':'outcome GWAS','dataset2':'eQTLGen cis-eQTL','H2':'eQTL only','verified_against':'coloc 5.2.3 beta/varbeta equations','equation_source':'https://raw.githubusercontent.com/chr1swallace/coloc/v5.2.3/R/claudia.R'}
(OUT/'coloc_verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
