from pathlib import Path
import json, hashlib, datetime, time, traceback
import pandas as pd
import numpy as np

T=Path(__file__).resolve().parent; W=T.parent
I=T/'inputs'; I.mkdir(exist_ok=True)
LOG=T/'analysis.log'
def log(s):
    text=datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s
    print(text,flush=True)
    with LOG.open('a',encoding='utf8') as f:f.write(text+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    assert not (I/'extraction_manifest.json').exists(),'Version already complete'
    plan=T/'ANALYSIS_PLAN.md'
    frozen={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'plan_sha256':sha(plan),'script_sha256':sha(Path(__file__))}
    (T/'plan_freeze.json').write_text(json.dumps(frozen,indent=2),encoding='utf8')
    d=W/'review_8dce23c/submission/candidate_20260911_review'
    iv=pd.read_csv(d/'Supplementary_Data_1_Instruments.csv')
    eq1=pd.read_csv(W/'maf_sensitivity_20260905/inputs/eqtlgen_9genes.tsv',sep='\t')
    eq2=pd.read_csv(W/'expanded_analysis_20260930/inputs/eqtl_mhc.tsv',sep='\t')
    eq=pd.concat([eq1,eq2],ignore_index=True)
    assert iv.shape[0]==6135 and iv.gene_symbol.nunique()==2544 and eq.GeneSymbol.nunique()==14
    assert not eq.duplicated(['GeneSymbol','SNP']).any()
    assert eq.groupby('GeneSymbol').Pvalue.max().gt(.9).all()
    keys=set(iv.snp)|set(eq.SNP)
    eq.to_csv(I/'eqtl_14genes.tsv',sep='\t',index=False)
    log(f'PLAN FROZEN {frozen}; extract target SNPs={len(keys)}; no outcome P filter')
    manifest=[]
    for oc in ['E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']:
        src=Path('C:/ProjectTEDGWAS')/f'finngen_R12_{oc}.gz'
        start=time.time(); n=0; selected=[]; missing=0; aliases=0
        bad={k:0 for k in ['beta','se','pvalue','eaf','position','alleles']}; nsig=0; minp=1.; maxp=0.
        for c in pd.read_csv(src,sep='\t',compression='gzip',chunksize=400000,dtype={'#chrom':str,'rsids':str,'ref':str,'alt':str}):
            n+=len(c)
            bad['beta']+=int((~np.isfinite(c.beta)).sum())
            bad['se']+=int((~np.isfinite(c.sebeta)|c.sebeta.le(0)).sum())
            bad['pvalue']+=int((~np.isfinite(c.pval)|~c.pval.between(0,1)).sum())
            bad['eaf']+=int((~np.isfinite(c.af_alt)|~c.af_alt.between(0,1)).sum())
            bad['position']+=int((~np.isfinite(c.pos)|c.pos.le(0)).sum())
            bad['alleles']+=int((c.ref.isna()|c.alt.isna()|c.ref.eq(c.alt)).sum())
            missing+=int(c.rsids.isna().sum());nsig+=int(c.pval.lt(5e-8).sum())
            minp=min(minp,float(c.pval.min()));maxp=max(maxp,float(c.pval.max()))
            hits=c[c.rsids.isin(keys)].copy(); hits['matched_rsid']=hits.rsids
            multi=c[c.rsids.fillna('').str.contains(',',regex=False)]
            extra=[]
            for ix,row in multi.iterrows():
                for snp in set(str(row.rsids).split(','))&keys:
                    r=row.copy();r['matched_rsid']=snp;extra.append(r);aliases+=1
            if extra:hits=pd.concat([hits,pd.DataFrame(extra)],ignore_index=True)
            selected.append(hits)
            if n%4000000==0:log(f'{oc}: scanned {n:,}; selected {sum(len(x) for x in selected):,}')
        raw=pd.concat(selected,ignore_index=True)
        raw.to_csv(I/f'{oc}_selected_raw.tsv',sep='\t',index=False)
        out=raw.rename(columns={'matched_rsid':'snp','#chrom':'chr','alt':'ea','ref':'oa','sebeta':'se','af_alt':'eaf','pval':'pvalue'})
        out=out[['snp','chr','pos','ea','oa','beta','se','eaf','pvalue']].copy();out['genome_build']='GRCh38'
        duplicate=out[out.duplicated('snp',keep=False)]
        duplicate.to_csv(I/f'{oc}_duplicate_matches.csv',index=False)
        # Multiple incompatible records for one rsID must be resolved explicitly, not chosen by significance.
        assert duplicate.empty, f'{oc}: ambiguous rsID records; investigate before MR'
        assert not any(bad.values()),f'{oc}: invalid numeric records {bad}'
        out.to_csv(I/f'{oc}_selected.tsv',sep='\t',index=False)
        q={'outcome':oc,'source':str(src),'source_bytes':src.stat().st_size,'source_sha256':sha(src),'rows':n,'matched_rows':len(out),'multi_rsid_alias_hits':aliases,'missing_rsid':missing,'invalid':bad,'p_lt_5e_8_variant_rows':nsig,'p_min':minp,'p_max':maxp,'elapsed_seconds':time.time()-start,'selected_sha256':sha(I/f'{oc}_selected.tsv')}
        manifest.append(q);log(json.dumps(q))
        (I/'extraction_progress.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
    (I/'extraction_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
    log('EXTRACTION PASS')
if __name__=='__main__':
    try:main()
    except Exception:
        log('FAILED extraction: '+traceback.format_exc());raise
