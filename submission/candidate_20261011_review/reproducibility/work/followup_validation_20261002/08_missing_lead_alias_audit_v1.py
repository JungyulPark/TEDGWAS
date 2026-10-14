from pathlib import Path
import pandas as pd,json,re,time
T=Path(__file__).resolve().parent;O=T/'alias_audit_v1';I=T/'inputs'
def main():
    assert (T/'final_checks_v1/verification.json').exists();O.mkdir()
    cov=pd.read_csv(T/'final_checks_v1/marginal_signal_coverage.csv');manifest=json.loads((I/'outcomes_manifest.json').read_text())
    records=[];hits=[]
    for rec in manifest:
        oc=rec['name']
        if oc not in ['BBJ','UKB','FinnGen_GO','E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']:continue
        targets=set(cov[(cov.outcome==oc)&cov.source_lead_status.eq('absent_exact_rsid_in_outcome')].source_eqtl_lead)
        col='rsid' if oc=='BBJ' else ('hm_rsid' if oc=='UKB' else 'rsids')
        pattern='(?:^|[,;])(?:'+'|'.join(sorted(targets))+')(?:[,;]|$)'
        n=0;matched=0;start=time.time()
        for d in pd.read_csv(rec['path'],sep='\t',usecols=[col],dtype=str,chunksize=750000):
            n+=len(d);sel=d[d[col].str.contains(pattern,regex=True,na=False)]
            for v in sel[col]:
                for t in targets.intersection(re.split('[,;]',v)):
                    hits.append(dict(outcome=oc,target_rsid=t,source_rsid_field=v));matched+=1
        records.append(dict(outcome=oc,missing_source_leads=len(targets),rows_scanned=n,exact_or_comma_semicolon_alias_hits=matched,seconds=time.time()-start))
        print(records[-1],flush=True)
    pd.DataFrame(records).to_csv(O/'summary.csv',index=False)
    pd.DataFrame(hits,columns=['outcome','target_rsid','source_rsid_field']).to_csv(O/'hits.csv',index=False)
    (O/'status.json').write_text(json.dumps(dict(complete=True,alias_hits=len(hits),interpretation='No missing lead recovered by source rsID aliases' if not hits else 'Alias recovery requires revised matching before final interpretation'),indent=2),encoding='utf8')
if __name__=='__main__':main()
