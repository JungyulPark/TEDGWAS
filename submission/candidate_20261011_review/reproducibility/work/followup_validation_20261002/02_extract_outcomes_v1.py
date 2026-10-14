from pathlib import Path
import sys,json,time,hashlib,datetime,traceback
import pandas as pd
T=Path(__file__).resolve().parent; I=T/'inputs'; W=T.parent
def log(s):
    line=datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s
    print(line,flush=True)
    with (T/'run.log').open('a',encoding='utf8') as f:f.write(line+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def main():
    assert (I/'eqtl_manifest.json').exists()
    assert not (I/'outcomes_manifest.json').exists()
    keys=set(pd.read_csv(I/'eqtl_11genes.tsv',sep='\t',usecols=['SNP']).SNP)
    manifests=[]
    root=Path('C:/ProjectTEDGWAS'); base=root/'TrackA_MR/v5_upgrade/04_druggable_mr/data'
    sources=[('eur_freq',base/'g1000_eur_freq.frq','SNP',r'\s+'),
      ('eqtl_af',W/'data/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz','SNP','\t')]
    for name,src,key,sep in sources:
        n=0; parts=[]
        for d in pd.read_csv(src,sep=sep,chunksize=500000):
            n+=len(d);parts.append(d[d[key].isin(keys)].copy())
        out=pd.concat(parts,ignore_index=True);assert not out[key].duplicated().any()
        out.to_csv(I/f'{name}.tsv',sep='\t',index=False)
        if name=='eqtl_af':
            den=out.allA_total+out.allAB_total+out.allB_total
            diff=((out.allAB_total+2*out.allB_total)/(2*den)-out.AlleleB_all).abs().max()
            assert diff<1e-9,diff
        manifests.append(dict(name=name,path=str(src),bytes=src.stat().st_size,sha256=sha(src),rows=n,selected=len(out)))
        log(f'{name}: {n:,} rows scanned; {len(out):,} selected')
    outcomes=[('BBJ',base/'outcomes/GCST90018627_harmonised.tsv.gz','BBJ'),
      ('UKB',base/'outcomes/GCST90038636_harmonised.tsv.gz','UKB'),
      ('FinnGen_GO',root/'finngen_R12_GRAVES_OPHT.gz','FG'),
      ('E4_GRAVES_STRICT',root/'finngen_R12_E4_GRAVES_STRICT.gz','FG'),
      ('E4_GRAVES_OPHT_STRICT',root/'finngen_R12_E4_GRAVES_OPHT_STRICT.gz','FG')]
    for name,src,kind in outcomes:
        n=0; parts=[];regions=[];start=time.time()
        for raw in pd.read_csv(src,sep='\t',chunksize=400000,low_memory=False):
            n+=len(raw)
            if kind=='FG':
                d=raw.rename(columns={'#chrom':'chr','pos':'pos38','rsids':'snp','alt':'ea','ref':'oa','sebeta':'se','pval':'p','af_alt':'eaf'})
            elif kind=='BBJ':
                d=raw.rename(columns={'chromosome':'chr','base_pair_location':'pos38','rsid':'snp','effect_allele':'ea','other_allele':'oa','standard_error':'se','p_value':'p','effect_allele_frequency':'eaf'})
            else:
                d=raw.rename(columns={'chromosome':'chr','base_pair_location':'pos38','hm_rsid':'snp','effect_allele':'ea','other_allele':'oa','standard_error':'se','p_value':'p','effect_allele_frequency':'eaf'})
            d=d[['chr','pos38','snp','ea','oa','beta','se','p','eaf']].copy()
            d['chr']=pd.to_numeric(d['chr'],errors='coerce')
            if kind=='UKB':
                mu=3731/484598;d['beta']=d.beta/(mu*(1-mu));d['se']=d.se/(mu*(1-mu))
            parts.append(d[d.snp.isin(keys)].copy())
            regions.append(d[d.chr.eq(14)&d.pos38.between(80000000,82100000)].copy())
        selected=pd.concat(parts,ignore_index=True);region=pd.concat(regions,ignore_index=True)
        selected.to_csv(I/f'{name}.tsv',sep='\t',index=False);region.to_csv(I/f'{name}_TSHR_raw_region.tsv',sep='\t',index=False)
        assert len(region)>1000 and len(selected)>10000
        record=dict(name=name,path=str(src),bytes=src.stat().st_size,sha256=sha(src),rows=n,selected=len(selected),raw_audit_region_rows=len(region),genome_build='GRCh38',rsid_rule='exact source rsID string; no alias expansion',scale='log odds' if kind!='UKB' else 'linear beta and SE each divided once by mu*(1-mu)',seconds=time.time()-start)
        manifests.append(record);log(json.dumps(record))
    (I/'outcomes_manifest.json').write_text(json.dumps(manifests,indent=2),encoding='utf8')
    log('Input extraction PASS')
if __name__=='__main__':
    try:main()
    except Exception:log('FAILED '+traceback.format_exc());raise
