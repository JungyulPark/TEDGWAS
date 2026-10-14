from coloc_core_v1 import *
import json,datetime,traceback
O=T/'candidates_v1'
def log(s):
    line=datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s;print(line,flush=True)
    with (O/'run.log').open('a',encoding='utf8') as f:f.write(line+'\n')
def main():
    assert (T/'TSHR_v2/R_verification.csv').exists();O.mkdir();V=O/'r_validation';V.mkdir()
    e=expression();e=e[e.GeneSymbol.ne('TSHR')];results=[];coverage=[]
    for oc in OUTCOMES:
        matched_all=matched(e,oc)
        for gene,m in matched_all.groupby('GeneSymbol'):
            ge=e[e.GeneSymbol.eq(gene)]
            coverage.append(dict(gene=gene,outcome=oc,n_eqtl=len(ge),p_eqtl_max=ge.Pvalue.max(),min_eQTL_N=ge.NrSamples.min(),max_eQTL_N=ge.NrSamples.max(),n_allele_compatible=len(m),n_reference=int(m.valid_ref.sum()),n_eqtlgen=int(m.valid_new.sum()),n_common=int((m.valid_ref&m.valid_new).sum()),min_GRCh37=ge.SNPPos.min(),max_GRCh37=ge.SNPPos.max()))
            for sc,kind in SCENARIOS:
                d=subset(m,sc);rows,validation=calculate(d,kind)
                for row in rows:results.append(dict(gene=gene,outcome=oc,window='full_source_cis',scenario=sc,**row))
                validation.to_csv(V/f'{gene}__{oc}__full_source_cis__{sc}.csv',index=False)
        log(f'{oc}: ten candidate loci completed')
    got=pd.DataFrame(results);assert len(got)==600
    got.to_csv(O/'Candidate_600_coloc_settings.csv',index=False);pd.DataFrame(coverage).to_csv(O/'coverage.csv',index=False)
    iv=pd.read_csv(W/'maf_sensitivity_20260905/inputs/instruments_verified.csv');iv=iv[iv.gene_symbol.isin(e.GeneSymbol.unique())]
    shared=iv.groupby('snp').gene_symbol.agg(lambda s:';'.join(sorted(s.unique()))).reset_index(name='genes');shared=shared[shared.genes.str.contains(';')]
    shared.to_csv(O/'shared_selected_instruments.csv',index=False)
    main=got[got.scenario.eq('original_reference')&got.p12.eq(1e-5)].copy()
    mr=pd.read_csv(W.parent/'outputs/TED_Reanalysis_20261002/Tables/All_2544_genes_five_outcomes.csv')
    main['outcome']=main.outcome.replace(ALIASES)
    print('MR columns',mr.columns.tolist())
    main=main.merge(mr,left_on=['gene','outcome'],right_on=['gene_symbol','outcome'],how='left',validate='one_to_one')
    main.to_csv(O/'Candidate_primary_MR_coloc.csv',index=False)
    (O/'python_verification.json').write_text(json.dumps(dict(settings=600,n_genes=10,n_outcomes=5,settings_with_adjusted_p1=int(got.p1_effective.ne(1e-4).sum()),status='Python PASS; R pending'),indent=2),encoding='utf8')
    print(got[got.scenario.eq('original_reference')&got.p12.eq(1e-5)].pivot(index='gene',columns='outcome',values='PP.H4').to_string());log('Python candidates PASS')
if __name__=='__main__':
    try:main()
    except Exception:
        if O.exists():log('FAILED '+traceback.format_exc())
        raise
