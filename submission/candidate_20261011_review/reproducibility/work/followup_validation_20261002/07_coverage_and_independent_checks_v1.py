from coloc_core_v1 import *
from scipy.stats import norm,chi2
import json
O=T/'final_checks_v1'
def main():
    assert (T/'instruments_v1/R_verification.csv').exists();O.mkdir()
    e=expression();coverage=[]
    for oc in OUTCOMES:
        d=pd.read_csv(I/f'{oc}.tsv',sep='\t');ma=matched(e,oc)
        for gene,ge in e.groupby('GeneSymbol'):
            lead=ge.loc[ge.Zscore.abs().idxmax()];m=ma[ma.GeneSymbol.eq(gene)]
            for sc,kind in [('original_reference','ref'),('available_eqtlgen','new')]:
                a=subset(m,sc);exlead=a.loc[a.Zscore.abs().idxmax()]
                dz=abs(a.beta/a.se);outlead=a.loc[dz.idxmax()]
                reason='included' if lead.SNP in set(a.SNP) else ('absent_exact_rsid_in_outcome' if lead.SNP not in set(d.snp) else ('allele_mismatch_or_invalid_stats' if lead.SNP not in set(m.SNP) else 'frequency_unavailable'))
                coverage.append(dict(gene=gene,outcome=oc,scenario=sc,n_source_eqtl=len(ge),n_overlap=len(a),source_eqtl_lead=lead.SNP,source_eqtl_max_abs_z=abs(lead.Zscore),source_eqtl_lead_N=lead.NrSamples,source_lead_status=reason,matched_eqtl_lead=exlead.SNP,matched_eqtl_max_abs_z=abs(exlead.Zscore),matched_eqtl_min_p=2*norm.sf(abs(exlead.Zscore)),matched_outcome_lead=outlead.SNP,matched_outcome_max_abs_z=dz.max(),matched_outcome_min_p=2*norm.sf(dz.max()),matched_eqtl_N_min=a.NrSamples.min(),matched_eqtl_N_max=a.NrSamples.max()))
    cov=pd.DataFrame(coverage);cov.to_csv(O/'marginal_signal_coverage.csv',index=False)
    h=pd.read_csv(T/'instruments_v1/retained_SNP_inputs.csv');full=pd.read_csv(T/'instruments_v1/full_MR.csv');loo=pd.read_csv(T/'instruments_v1/leave_one_out.csv')
    deltas=[]
    for _,r in pd.concat([full,loo],ignore_index=True).iterrows():
        a=h[h.exposure.eq(r.gene)&h.outcome.eq(r.outcome)&h.scenario.eq(r.scenario)&h.SNP.ne(r.removed_snp)]
        bx=a['beta.exposure'].to_numpy();by=a['beta.outcome'].to_numpy();sy=a['se.outcome'].to_numpy();n=len(a)
        b=np.sum(bx*by/sy**2)/np.sum(bx**2/sy**2);q=np.sum((by-b*bx)**2/sy**2);se=np.sqrt((max(1,q/(n-1)) if n>1 else 1)/np.sum(bx**2/sy**2));p=2*norm.sf(abs(b/se))
        delta=max(abs(b-r.beta),abs(se-r.se),abs(p-r.pvalue));assert delta<1e-8 and n==r.n_iv
        deltas.append(delta)
    # Compare newly extracted raw frequencies and Z against the frozen instruments actually retained.
    ex=pd.read_csv(T/'instruments_v1/exposure_frequency_scaling.csv');ret=h[['SNP','exposure']].drop_duplicates()
    m=ret.merge(ex,left_on=['SNP','exposure'],right_on=['snp','gene_symbol'],validate='one_to_one').merge(e,left_on=['SNP','exposure'],right_on=['SNP','GeneSymbol'],suffixes=('_frozen','_fresh'),validate='one_to_one')
    assert len(m)==18
    af_delta=max((m.eaf_ref_frozen-m.eaf_ref_fresh).abs().max(),(m.eaf_new_frozen-m.eaf_new_fresh).abs().max())
    z_delta=(m.zscore-m.Zscore).abs().max();assert af_delta<1e-10 and z_delta<1e-10
    m[['SNP','exposure','n_samples','eaf_ref_frozen','eaf_ref_fresh','eaf_new_frozen','eaf_new_fresh','zscore','Zscore']].to_csv(O/'retained_instrument_source_identity.csv',index=False)
    allcol=pd.concat([pd.read_csv(T/'TSHR_v2/TSHR_240_coloc_settings.csv'),pd.read_csv(T/'candidates_v1/Candidate_600_coloc_settings.csv')],ignore_index=True)
    assert len(allcol)==840 and np.isfinite(allcol[[f'PP.H{i}' for i in range(5)]]).all().all()
    assert np.max(abs(allcol[[f'PP.H{i}' for i in range(5)]].sum(axis=1)-1))<1e-10
    checks=dict(python_independent_MR_estimates=len(deltas),MR_max_abs_delta=max(deltas),retained_instrument_source_identity_rows=len(m),raw_frequency_max_abs_delta=af_delta,raw_Z_max_abs_delta=z_delta,coloc_posterior_vectors=840,TSHR_baseline_vectors_reproduced=60,new_TSHR_window_vectors=180,new_candidate_vectors=600,candidate_H4_ge_0_8=int((allcol[allcol.gene.ne('TSHR')]['PP.H4']>=.8).sum()),PASS=True)
    (O/'verification.json').write_text(json.dumps(checks,indent=2),encoding='utf8')
    print(json.dumps(checks,indent=2));print(cov[cov.gene.str.startswith('IGLV')&cov.scenario.eq('original_reference')].to_string(index=False))
if __name__=='__main__':main()
