from coloc_core_v1 import *
import json,datetime,traceback
O=T/'TSHR_v1'
def log(s):
    line=datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s;print(line,flush=True)
    with (O/'run.log').open('a',encoding='utf8') as f:f.write(line+'\n')
def main():
    assert (I/'outcomes_manifest.json').exists();O.mkdir();V=O/'r_validation';V.mkdir()
    e=expression();e=e[e.GeneSymbol.eq('TSHR')].copy();e.to_csv(O/'expression_frequency_audit.csv',index=False)
    matches={oc:matched(e,oc) for oc in OUTCOMES}
    sets={sc:set.intersection(*[set(subset(m,sc).SNP) for m in matches.values()]) for sc,_ in SCENARIOS}
    center=json.loads((T/'frozen_windows.json').read_text())['source_GenePos_GRCh37']
    results=[];counts=[];leads=[];rawframes=[]
    for oc,m in matches.items():
        raw=pd.read_csv(I/f'{oc}_TSHR_raw_region.tsv',sep='\t');raw['outcome']=oc
        raw['abs_z']=abs(raw.beta/raw.se)
        raw['in_eqtl_by_exact_rsid']=raw.snp.isin(e.SNP)
        raw['allele_compatible']=raw.snp.isin(m.SNP)
        raw['in_primary_coloc']=raw.snp.isin(subset(m,'original_reference').SNP)
        raw['in_eqtlgen_coloc']=raw.snp.isin(subset(m,'available_eqtlgen').SNP)
        raw['within_mapped_cis_bounds']=raw.pos38.between(m.pos38.min(),m.pos38.max())
        raw['exclusion_reason']=np.select([~raw.in_eqtl_by_exact_rsid,~raw.allele_compatible,~raw.in_primary_coloc],['no_exact_rsid_in_source_eQTL','incompatible_alleles_or_invalid_beta_se','reference_frequency_unavailable_or_incompatible'],default='included_primary')
        rawframes.append(raw)
        for scope,df in [('broad_audit_interval',raw),('within_mapped_cis_bounds',raw[raw.within_mapped_cis_bounds])]:
            valid=df[df.se.gt(0)&np.isfinite(df.abs_z)].copy()
            top=valid.nlargest(10,'abs_z').copy();top['scope']=scope;top['rank']=range(1,len(top)+1);leads.append(top)
            for reason,group in valid.groupby('exclusion_reason'):
                counts.append(dict(outcome=oc,scope=scope,reason=reason,rows=len(group),max_abs_z=group.abs_z.max(),min_p=group.p.min(),mapped_cis_min38=m.pos38.min(),mapped_cis_max38=m.pos38.max()))
        for window in ['full_source_cis','GenePos_pm500kb','GenePos_pm250kb','five_outcome_intersection']:
            wm=m if window in ['full_source_cis','five_outcome_intersection'] else m[m.SNPPos.between(center-(500000 if '500' in window else 250000),center+(500000 if '500' in window else 250000))]
            for sc,kind in SCENARIOS:
                d=subset(wm,sc)
                if window=='five_outcome_intersection':d=d[d.SNP.isin(sets[sc])]
                rows,validation=calculate(d,kind)
                for row in rows:results.append(dict(gene='TSHR',outcome=oc,window=window,scenario=sc,**row))
                validation.to_csv(V/f'TSHR__{oc}__{window}__{sc}.csv',index=False)
        log(f'{oc}: raw regional rows {len(raw)}; primary matched {len(subset(m,"original_reference"))}')
    got=pd.DataFrame(results);assert len(got)==240
    old=pd.read_csv(W.parent/'outputs/TED_Reanalysis_20261002/Tables/Colocalization_840_settings.csv');old=old[old.gene.eq('TSHR')].copy()
    baseline=got[got.window.eq('full_source_cis')].copy();baseline['outcome']=baseline.outcome.replace(ALIASES)
    compare=baseline.merge(old,on=['gene','outcome','scenario','p12'],suffixes=('_new','_old'),validate='one_to_one')
    delta=max((compare[f'PP.H{i}_new']-compare[f'PP.H{i}_old']).abs().max() for i in range(5))
    assert len(compare)==60 and delta<1e-8 and compare.n_overlap_new.eq(compare.n_overlap_old).all()
    got.to_csv(O/'TSHR_240_coloc_settings.csv',index=False);compare.to_csv(O/'baseline_identity.csv',index=False)
    pd.DataFrame(counts).to_csv(O/'regional_coverage_counts.csv',index=False);pd.concat(leads,ignore_index=True).to_csv(O/'raw_GWAS_top10_by_scope.csv',index=False)
    pd.concat(rawframes,ignore_index=True).to_csv(O/'raw_GWAS_region_audit.csv',index=False)
    (O/'python_verification.json').write_text(json.dumps(dict(baseline_vectors=60,baseline_max_delta=delta,new_window_vectors=180,intersection_sizes={k:len(v) for k,v in sets.items()},status='Python PASS; R pending'),indent=2),encoding='utf8')
    print(got[got.scenario.eq('original_reference')&got.p12.eq(1e-5)][['outcome','window','n_overlap','PP.H3','PP.H4']].to_string(index=False));log('Python TSHR PASS')
if __name__=='__main__':
    try:main()
    except Exception:
        if O.exists():log('FAILED '+traceback.format_exc())
        raise
