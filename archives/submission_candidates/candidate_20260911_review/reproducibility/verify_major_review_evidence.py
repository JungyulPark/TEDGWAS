"""Verify reviewer claims against local inputs; export aggregate evidence only.

Example: python verify_major_review_evidence.py --input-dir LOCAL_MAF_INPUTS
The reference LD check is reproducible with PLINK --bfile EUR_chr14 (or
EAS_chr14) --extract variants.txt --r2 --ld-window 99999 --ld-window-kb 2000
--ld-window-r2 0 --out EUR (or EAS). variants.txt contains rs179252, rs179247,
rs12101255, one per line. No significance filter or phenotype is used.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd
parser=argparse.ArgumentParser();parser.add_argument('--input-dir',type=Path,required=True);a=parser.parse_args()
D=Path(__file__).resolve().parents[1];P=D/'provenance'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
mr=pd.read_csv(P/'MR_primary_canonical.csv');co=pd.read_csv(P/'coloc_canonical_v2.csv');co=co[np.isclose(co.p12,1e-5,atol=1e-12)]
eq=pd.read_csv(a.input_dir/'eqtlgen_9genes.tsv',sep='\t');eq=eq[(eq.GeneSymbol=='TSHR')&(eq.SNP=='rs179252')].iloc[0]
freq=pd.read_csv(a.input_dir/'eur_freq_selected.tsv',sep='\t');freq=freq[freq.SNP=='rs179252'].iloc[0]
exposure_se=1/np.sqrt(2*freq.MAF*(1-freq.MAF)*(eq.NrSamples+eq.Zscore**2));exposure_beta=eq.Zscore*exposure_se
results=[]
for prefix,outcome in [('BBJ','BBJ_Graves'),('UKB','UKB_hyperthyroid'),('FinnGen','FinnGen_GO')]:
    src=pd.read_csv(a.input_dir/(prefix+'_selected.tsv'),sep='\t');r=src[src.snp=='rs179252'].iloc[0]
    assert {r.ea,r.oa}=={eq.AssessedAllele,eq.OtherAllele}
    sign=1 if r.ea==eq.AssessedAllele else -1
    mu=3731/484598;scale=mu*(1-mu) if prefix=='UKB' else 1
    b=float(sign*r.beta/scale);se=float(r.se/scale)
    primary=mr[(mr.gene_symbol=='TSHR')&(mr.outcome==outcome)].iloc[0]
    assert np.isclose(b/exposure_beta,primary.beta,rtol=1e-12)
    assert np.isclose(se/abs(exposure_beta),primary.se,rtol=1e-12)
    results.append({'outcome':outcome,'variant_log_or':b,'variant_or':float(np.exp(b)),'exposure_beta':float(exposure_beta),'MR_beta':float(primary.beta),'MR_SE':float(primary.se)})
f=co[co.outcome=='FinnGen_GO'];h2=f[[f'PP.H{i}' for i in range(5)]].idxmax(axis=1)=='PP.H2';ct=f[f.gene=='CTLA4'].iloc[0]
assert h2.sum()==7 and ct['PP.H4']>=.8
ig=co[co.gene=='IGF1R'];assert (ig.outcome_min_p<.05).all()
ld=pd.read_csv(P/'TSHR_reference_LD_review.csv');assert len(ld)==6
j={'status':'PASS','TSHR_scale_check':results,'TSHR_effect_allele':eq.AssessedAllele,
   'TSHR_reference_LD':ld.to_dict('records'),'FinnGen_H2_dominant_loci':int(h2.sum()),
   'CTLA4_FinnGen_H4':float(ct['PP.H4']),'IGF1R_regional_model_and_marginal_P':ig[['outcome','PP.H2','PP.H3','PP.H4','outcome_min_p']].to_dict('records'),
   'input_hashes':{n:sha(a.input_dir/n) for n in ['eqtlgen_9genes.tsv','eur_freq_selected.tsv','BBJ_selected.tsv','UKB_selected.tsv','FinnGen_selected.tsv']},
   'top_snp_definition':'argmax of outcome log-ABF + exposure log-ABF; conditional shared-variant ranking, not necessarily the GWAS lead SNP',
   'not_completed':['Multi-signal colocalization for UKB TSHR and BBJ CTLA4','Formal genetic directionality','Full Phase 1 contributing-cohort table inspection; participant overlap remains unquantified'],
   'interpretation':'The aggregate checks and reference-panel LD do not validate a multi-causal model or establish thymic expression direction.'}
(P/'major_review_evidence.json').write_text(json.dumps(j,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:j[k] for k in ['status','TSHR_scale_check','FinnGen_H2_dominant_loci','not_completed']},indent=2))
