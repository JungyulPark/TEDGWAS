from pathlib import Path
import sys
T=Path(__file__).resolve().parent;W=T.parent;I=T/'inputs'
sys.path.insert(0,str(W/'pydeps'))
import numpy as np,pandas as pd
from scipy.special import logsumexp
OUTCOMES=['BBJ','UKB','FinnGen_GO','E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']
ALIASES={'BBJ':'BBJ_Graves','UKB':'UKB_hyperthyroid'}
SCENARIOS=[('original_reference','ref'),('paired_reference','ref'),('paired_eqtlgen','new'),('available_eqtlgen','new')]
def expression():
    e=pd.read_csv(I/'eqtl_11genes.tsv',sep='\t');ref=pd.read_csv(I/'eur_freq.tsv',sep='\t');af=pd.read_csv(I/'eqtl_af.tsv',sep='\t')
    assert not ref.SNP.duplicated().any() and not af.SNP.duplicated().any()
    e=e.merge(ref[['SNP','A1','A2','MAF']],on='SNP',how='left',validate='many_to_one').merge(af[['SNP','hg19_chr','hg19_pos','AlleleA','AlleleB','AlleleB_all']],on='SNP',how='left',validate='many_to_one')
    refok=((e.AssessedAllele==e.A1)&(e.OtherAllele==e.A2))|((e.AssessedAllele==e.A2)&(e.OtherAllele==e.A1))
    newok=(((e.AssessedAllele==e.AlleleA)&(e.OtherAllele==e.AlleleB))|((e.AssessedAllele==e.AlleleB)&(e.OtherAllele==e.AlleleA)))&e.SNPChr.eq(e.hg19_chr)&e.SNPPos.eq(e.hg19_pos)
    e['eaf_ref']=np.where(refok,np.where(e.AssessedAllele.eq(e.A1),e.MAF,1-e.MAF),np.nan)
    e['eaf_new']=np.where(newok,np.where(e.AssessedAllele.eq(e.AlleleB),e.AlleleB_all,1-e.AlleleB_all),np.nan)
    for k in ['ref','new']:
        f=e['eaf_'+k];e['valid_'+k]=f.between(0,1,inclusive='neither');e['se_'+k]=1/np.sqrt(2*f*(1-f)*(e.NrSamples+e.Zscore**2))
    return e
def matched(e,oc):
    d=pd.read_csv(I/f'{oc}.tsv',sep='\t');d=d[d.beta.notna()&np.isfinite(d.beta)&d.se.gt(0)&np.isfinite(d.se)].copy()
    m=e.merge(d,left_on=['SNP','SNPChr'],right_on=['snp','chr'],how='inner')
    same=m.AssessedAllele.eq(m.ea)&m.OtherAllele.eq(m.oa);swap=m.AssessedAllele.eq(m.oa)&m.OtherAllele.eq(m.ea)
    m=m[same|swap].copy()
    assert not m.duplicated(['Gene','SNP']).any(), 'Ambiguous duplicates after allele-pair filtering'
    return m.sort_values(['GeneSymbol','SNP'])
def subset(m,sc):
    return m[m.valid_ref] if sc=='original_reference' else (m[m.valid_new] if sc=='available_eqtlgen' else m[m.valid_ref&m.valid_new])
def abf(z,v,sd):
    r=sd**2/(sd**2+v)
    return .5*(np.log1p(-r)+r*z*z)
def calculate(d,kind):
    assert len(d)>=50
    z1=d.beta.to_numpy()/d.se.to_numpy();v1=d.se.to_numpy()**2;z2=d.Zscore.to_numpy();v2=d['se_'+kind].to_numpy()**2
    assert np.isfinite(np.r_[z1,v1,z2,v2]).all()
    l1=abf(z1,v1,.2);l2=abf(z2,v2,.15);s1=logsumexp(l1);s2=logsumexp(l2);ss=logsumexp(l1+l2)
    diff=ss-s1-s2;assert diff<=1e-10
    h3=s1+s2+np.log(-np.expm1(min(diff,0)))
    # Exactly replicate coloc 5.2.3 adjust_prior when requested p*n >= 1.
    n=len(d);p1=1/(n+1) if n*1e-4>=1 else 1e-4;p2=p1
    rows=[]
    for p12 in [1e-6,5e-6,1e-5]:
        p12eff=1/(n+1) if n*p12>=1 else p12
        lh=np.array([0,np.log(p1)+s1,np.log(p2)+s2,np.log(p1)+np.log(p2)+h3,np.log(p12eff)+ss]);pp=np.exp(lh-logsumexp(lh))
        assert np.isfinite(pp).all() and abs(pp.sum()-1)<1e-10
        rows.append(dict(p12=p12,p1_effective=p1,p2_effective=p2,p12_effective=p12eff,n_overlap=n,top_snp=d.SNP.iloc[int(np.argmax(l1+l2))],**{f'PP.H{i}':pp[i] for i in range(5)}))
    validation=pd.DataFrame({'snp':d.SNP.to_numpy(),'beta1':d.beta.to_numpy(),'varbeta1':v1,'beta2':z2*np.sqrt(v2),'varbeta2':v2})
    return rows,validation
