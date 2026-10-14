from pathlib import Path
import sys,json,datetime,traceback
T=Path(__file__).resolve().parent;W=T.parent
sys.path.insert(0,str(W/'pydeps'))
import numpy as np,pandas as pd
from scipy.special import logsumexp
O=T/'coloc_v1'
def log(s):
    line=datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s
    print(line,flush=True)
    with (O/'run.log').open('a',encoding='utf8') as f:f.write(line+'\n')
def abf(z,v,sd):
    r=sd**2/(sd**2+v)
    return .5*(np.log1p(-r)+r*z*z)
def main():
    assert (T/'mr_v1/verification.csv').exists(), 'MR must pass first'
    O.mkdir();V=O/'r_validation';V.mkdir()
    old=W/'maf_sensitivity_20260905/inputs';ext=W/'expanded_analysis_20260930/inputs'
    rows=[];coverage=[];region=[]
    groups=[('existing',pd.read_csv(old/'eqtlgen_9genes.tsv',sep='\t'),pd.read_csv(old/'eur_freq_selected.tsv',sep='\t'),pd.read_csv(old/'eqtlgen_af_selected_v1.csv')),
      ('mhc',pd.read_csv(ext/'eqtl_mhc.tsv',sep='\t'),pd.read_csv(ext/'eur_freq.tsv',sep='\t'),pd.read_csv(ext/'eqtl_af.tsv',sep='\t').rename(columns={'SNP':'snp','hg19_chr':'chr','hg19_pos':'pos','AlleleA':'allele_a','AlleleB':'allele_b','AlleleB_all':'af_b'}))]
    for group,eq,ref,af in groups:
        assert not ref.SNP.duplicated().any() and not af.snp.duplicated().any()
        e=eq.merge(ref[['SNP','A1','A2','MAF']],on='SNP',how='left',validate='many_to_one').merge(af[['snp','allele_a','allele_b','af_b']],left_on='SNP',right_on='snp',how='left',validate='many_to_one')
        e['eaf_ref']=np.where(e.AssessedAllele.eq(e.A1),e.MAF,np.where(e.AssessedAllele.eq(e.A2),1-e.MAF,np.nan))
        ok=((e.AssessedAllele==e.allele_b)&(e.OtherAllele==e.allele_a))|((e.AssessedAllele==e.allele_a)&(e.OtherAllele==e.allele_b))
        e['eaf_new']=np.where(ok,np.where(e.AssessedAllele==e.allele_b,e.af_b,1-e.af_b),np.nan)
        for kind in ['ref','new']:
            f=e['eaf_'+kind];e['valid_'+kind]=f.between(0,1,inclusive='neither');e['se_'+kind]=1/np.sqrt(2*f*(1-f)*(e.NrSamples+e.Zscore**2))
        for g,ge in e.groupby('GeneSymbol'):
            assert ge.Pvalue.max()>.9 and not ge.SNP.duplicated().any()
            region.append(dict(gene=g,n_eqtl=len(ge),p_min=ge.Pvalue.min(),p_max=ge.Pvalue.max(),chr=int(ge.SNPChr.iloc[0]),min_GRCh37=int(ge.SNPPos.min()),max_GRCh37=int(ge.SNPPos.max())))
        for oc in ['FinnGen_GO','E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']:
            path=(old/'FinnGen_selected.tsv' if group=='existing' else ext/'FinnGen.tsv') if oc=='FinnGen_GO' else T/'inputs_v2'/f'{oc}_selected.tsv'
            d=pd.read_csv(path,sep='\t');d=d[d.beta.notna()&d.se.gt(0)&d.snp.notna()].copy()
            for g,ge in e.groupby('GeneSymbol'):
                ch=int(ge.SNPChr.iloc[0]);m=ge.merge(d[d.chr.eq(ch)],left_on='SNP',right_on='snp',suffixes=('_e','_o'))
                same=(m.AssessedAllele==m.ea)&(m.OtherAllele==m.oa);swap=(m.AssessedAllele==m.oa)&(m.OtherAllele==m.ea)
                before=len(m);m=m[same|swap].copy();dups=int(m.SNP.duplicated().sum())
                assert dups==0, f'{g}/{oc} ambiguous duplicate AFTER allele matching'
                m=m.sort_values('SNP');a=m[m.valid_ref];b=m[m.valid_new];common=m[m.valid_ref&m.valid_new]
                coverage.append(dict(gene=g,outcome=oc,n_eqtl=len(ge),n_matched_before_alleles=before,n_allele_compatible=len(m),duplicates_after_alleles=dups,n_ref=len(a),n_new=len(b),n_common=len(common)))
                for data,scenario,kind in [(a,'original_reference','ref'),(common,'paired_reference','ref'),(common,'paired_eqtlgen','new'),(b,'available_eqtlgen','new')]:
                    assert len(data)>=50
                    z1=data.beta.to_numpy()/data.se.to_numpy();v1=data.se.to_numpy()**2;z2=data.Zscore.to_numpy();v2=data['se_'+kind].to_numpy()**2
                    assert np.isfinite(np.r_[z1,v1,z2,v2]).all()
                    l1=abf(z1,v1,.2);l2=abf(z2,v2,.15);s1=logsumexp(l1);s2=logsumexp(l2);ss=logsumexp(l1+l2)
                    diff=ss-s1-s2;assert diff<=1e-10
                    h3=s1+s2+np.log(-np.expm1(min(diff,0)))
                    pd.DataFrame({'snp':data.SNP.to_numpy(),'beta1':data.beta.to_numpy(),'varbeta1':v1,'beta2':z2*np.sqrt(v2),'varbeta2':v2}).to_csv(V/f'{g}__{oc}__{scenario}.csv',index=False)
                    for p12 in [1e-6,5e-6,1e-5]:
                        lh=np.array([0,np.log(1e-4)+s1,np.log(1e-4)+s2,2*np.log(1e-4)+h3,np.log(p12)+ss]);pp=np.exp(lh-logsumexp(lh))
                        assert abs(pp.sum()-1)<1e-10
                        rows.append(dict(gene=g,outcome=oc,scenario=scenario,p12=p12,n_overlap=len(data),top_snp=data.SNP.iloc[int(np.argmax(l1+l2))],**{f'PP.H{i}':pp[i] for i in range(5)}))
            log(f'{group} {oc} completed')
    got=pd.DataFrame(rows);assert len(got)==504
    frozen=pd.read_csv(W/'expanded_analysis_20260930/results/Stage2_all_504_colocalization_settings.csv')
    matched=frozen[frozen.outcome.eq('FinnGen_GO')].merge(got,on=['gene','outcome','scenario','p12'],suffixes=('_old','_new'),validate='one_to_one')
    delta=max((matched[f'PP.H{i}_old']-matched[f'PP.H{i}_new']).abs().max() for i in range(5))
    assert len(matched)==168 and delta<1e-9 and matched.n_overlap_old.eq(matched.n_overlap_new).all() and matched.top_snp_old.eq(matched.top_snp_new).all()
    got.to_csv(O/'Coloc_all_settings_with_identity.csv',index=False)
    got[got.outcome.ne('FinnGen_GO')].to_csv(O/'Coloc_new_336_settings.csv',index=False)
    pd.DataFrame(coverage).to_csv(O/'coverage.csv',index=False);pd.DataFrame(region).to_csv(O/'eqtl_region_check.csv',index=False)
    (O/'python_verification.json').write_text(json.dumps({'identity_settings':168,'identity_max_abs_delta':float(delta),'new_settings':336,'status':'Python PASS; R coloc validation pending'},indent=2),encoding='utf8')
    log(f'PASS: original FinnGen 168 posterior vectors reproduced; max delta {delta}')
if __name__=='__main__':
    try:main()
    except Exception:
        if O.exists():log('FAILED '+traceback.format_exc())
        raise
