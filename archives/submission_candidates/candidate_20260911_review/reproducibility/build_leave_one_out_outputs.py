"""Validate the SNP-exclusion analysis and export full-precision Data 6.

Run validate_leave_one_out.R on the preserved local harmonized input first,
then pass its aggregate output CSV using --r-validation.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np,pandas as pd
from scipy import stats
D=Path(__file__).resolve().parents[1];P=D/'provenance'
p=argparse.ArgumentParser();p.add_argument('--r-validation',type=Path,required=True);a=p.parse_args()
original=P/'posthoc_20260914/leave_one_out_estimates.csv'
loo=pd.read_csv(original);r=pd.read_csv(a.r_validation);primary=pd.read_csv(P/'MR_primary_canonical.csv')
assert len(loo)==20 and len(r)==24
keys=['gene','outcome','excluded_SNP'];cols=['n_iv','beta','se','pvalue','OR','CI_lower','CI_upper']
paired=loo.merge(r,on=keys,validate='one_to_one',suffixes=('_python','_R'))
assert len(paired)==20
errors={}
for col in cols:
    x=paired[col+'_python'].to_numpy();y=paired[col+'_R'].to_numpy()
    assert np.allclose(x,y,rtol=1e-10,atol=1e-30),col
    errors[col]=float(np.max(abs(x-y)))
base=r[r.excluded_SNP=='None (all instruments)'].merge(primary,left_on=['gene','outcome'],right_on=['gene_symbol','outcome'],validate='one_to_one',suffixes=('_R','_primary'))
assert len(base)==9
for col in ['beta','se','pvalue','n_iv']:assert np.allclose(base[col+'_R'],base[col+'_primary'],rtol=1e-10,atol=1e-30),col
assert np.allclose(loo.pvalue,2*stats.norm.sf(abs(loo.beta/loo.se)),rtol=1e-12,atol=1e-30)
for field,sign in [('CI_lower',-1),('CI_upper',1)]:assert np.allclose(loo[field],np.exp(loo.beta+sign*stats.norm.ppf(.975)*loo.se),rtol=1e-12)
omitted=loo[loo.excluded_SNP!='None (all instruments)']
assert len(omitted)==15 and sum(omitted.gene=='IGF1R')==11 and sum(omitted.gene=='CTLA4')==4
assert (omitted.n_iv==omitted.original_n_iv-1).all()
assert ((omitted.method=='Wald ratio')==(omitted.n_iv==1)).all()
assert (omitted.loc[omitted.gene=='IGF1R','beta']>0).all()
out=loo.copy();out['analysis_timing']='post_hoc';out['frequency_scenario']='original_reference'
out['test_distribution']='normal';out['ci_level']=.95
out.to_csv(D/'Supplementary_Data_6_Leave_one_out.csv',index=False,na_rep='NA',lineterminator='\r\n')
(P/'leave_one_out_R_validation.csv').write_bytes(a.r_validation.read_bytes())
j={'status':'PASS','analysis_timing':'Post hoc; integrated into manuscript, Figure S2 and Supplementary Data 6',
   'genes':['TSHR','IGF1R','CTLA4'],'baseline_estimates_reproduced':9,'eligible_comparisons':5,
   'omission_estimates':15,'IGF1R_omissions':11,'CTLA4_omissions':4,
   'R_version':'4.3.3','independent_implementation':'Base-R lm and Wald formulas; Python weighted cross-products and whitened least-squares',
   'max_absolute_R_Python_differences':errors,
   'source_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),
   'R_validation_sha256':hashlib.sha256(a.r_validation.read_bytes()).hexdigest(),
   'data6_sha256':hashlib.sha256((D/'Supplementary_Data_6_Leave_one_out.csv').read_bytes()).hexdigest(),
   'scope':'Original reference-frequency harmonized SNP sets. All eligible omissions reported. No age/sex subgroup, additional covariate adjustment, or new independent participants.',
   'limitations':['TSHR in all outcomes and CTLA4 in BBJ have only one SNP and cannot undergo SNP omission.',
                  'For two-SNP CTLA4 sets an omission leaves a single-SNP Wald ratio.',
                  'P-value changes reflect changes in precision as well as effect estimates; repeated tests are not independent replication.',
                  'This leave-one-out analysis was not repeated under substituted eQTLGen frequencies.']}
(P/'leave_one_out_verification.json').write_text(json.dumps(j,indent=2)+'\n',encoding='utf8',newline='\r\n')
print(json.dumps(j,indent=2))
