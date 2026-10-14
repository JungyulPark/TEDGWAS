from pathlib import Path
import sys,json,math,re,hashlib
from statistics import NormalDist
import pandas as pd,numpy as np
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Extension_20261003';T=O/'Tables';S=O/'Sources'
assert (S/'disease_source_manifest.json').exists()
assert not (T/'Tissue_MR_50_planned_combinations.csv').exists(),'Preserve completed estimates'
iv=pd.read_csv(T/'Selected_instruments_verified.csv');av=pd.read_csv(T/'Tissue_availability_10_combinations.csv')
outcomes=['E4_GRAVES_STRICT','BBJ_Graves','UKB_hyperthyroid','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
checks=[];h=[];result=[];source_p=[]
def check(name,ok):
 checks.append({'check':name,'pass':bool(ok)});assert ok,name
normal=NormalDist();zcrit=normal.inv_cdf(.975);mu=3731/484598
for oc in outcomes:
 d=pd.read_csv(S/'disease_extracts'/(oc+'.tsv'),sep='\t',dtype=str)
 for _,v in iv.iterrows():
  if oc=='BBJ_Graves':
   rows=d[(d.chromosome.astype(int)==v.chromosome)&(d.base_pair_location.astype(int)==v.pos38)]
   fields=['effect_allele','other_allele','beta','standard_error','p_value','effect_allele_frequency','rsid']
  elif oc=='UKB_hyperthyroid':
   rows=d[(d.hm_chrom.astype(int)==v.chromosome)&(d.hm_pos.astype(int)==v.pos38)]
   fields=['hm_effect_allele','hm_other_allele','hm_beta','standard_error','p_value','hm_effect_allele_frequency','hm_rsid']
  else:
   rows=d[(d['#chrom'].astype(int)==v.chromosome)&(d.pos.astype(int)==v.pos38)]
   fields=['alt','ref','beta','sebeta','pval','af_alt','rsids']
  check(oc+' '+v.snp+' one coordinate record',len(rows)<=1)
  if rows.empty:continue
  r=rows.iloc[0];ea,oa,b,sy,p,eaf,ids=[r[x] for x in fields];b,sy,p,eaf=map(float,[b,sy,p,eaf])
  check(oc+' '+v.snp+' identifier agrees',v.snp in ids.split(','))
  check(oc+' '+v.snp+' valid source stats',math.isfinite(b) and sy>0 and 0<p<=1 and 0<eaf<1)
  if oc.startswith('E4') or oc=='FinnGen_GO':check(oc+' '+v.snp+' GRCh38 REF/ALT identical',ea==v.effect_allele and oa==v.other_allele)
  if oc=='UKB_hyperthyroid':
   rsign=1 if r.effect_allele==ea and r.other_allele==oa else -1 if r.effect_allele==oa and r.other_allele==ea else None
   check('UKB raw vs harmonized sign '+v.snp,rsign is not None and math.isclose(b,rsign*float(r.beta),rel_tol=1e-12,abs_tol=1e-12))
  scale=1/(mu*(1-mu)) if oc=='UKB_hyperthyroid' else 1
  sign=1 if ea==v.effect_allele and oa==v.other_allele else -1 if ea==v.other_allele and oa==v.effect_allele else None
  check(oc+' '+v.snp+' nonpalindromic exact allele alignment',sign is not None and set([ea,oa]) not in [set('AT'),set('CG')])
  h.append({'gene':v.gene,'tissue':v.tissue,'outcome':oc,'snp':v.snp,'variant_id':v.variant_id,'genome_build':'GRCh38','pos38':v.pos38,'effect_allele':v.effect_allele,'other_allele':v.other_allele,'exposure_increasing_allele':v.effect_allele if v.slope>0 else v.other_allele,'beta_exposure':v.slope,'se_exposure':v.slope_se,'p_exposure':v.pval_nominal,'F':v.F,'n_exposure':v.n_exposure,'eaf_exposure':v.exposure_eaf,'source_effect_allele':ea,'source_other_allele':oa,'source_beta_outcome':b,'source_se_outcome':sy,'source_p_outcome':p,'source_eaf_outcome':eaf,'scale_factor':scale,'allele_sign':sign,'beta_outcome':b*scale*sign,'se_outcome':sy*scale,'eaf_outcome':eaf if sign==1 else 1-eaf})
harm=pd.DataFrame(h);harm.to_csv(T/'Harmonized_instruments_all_outcomes.csv',index=False)
for _,a in av.iterrows():
 for oc in outcomes:
  family='primary' if oc=='E4_GRAVES_STRICT' else 'secondary';m=10 if family=='primary' else 40;threshold=.05/m
  d=harm[(harm.gene==a.gene)&(harm.tissue==a.tissue)&(harm.outcome==oc)]
  r={'gene':a.gene,'tissue':a.tissue,'outcome':oc,'family':family,'planned_family_tests':m,'corrected_threshold':threshold,'n_exposure':int(a.n_exposure),'n_selected_instruments':int(a.n_selected_instruments),'n_retained_instruments':len(d),'genome_build':'GRCh38','expression_unit':'one unit of inverse-normal-transformed normalized expression','coloc_status':'not_run_full_regional_statistics_unavailable','Q':None,'Q_df':None,'Q_pvalue':None,'LOO_status':'not_applicable_single_instrument' if len(d)==1 else 'not_estimable_no_instrument','direction_vs_blood':'not_estimable'}
  if len(d)==0:
   r.update(status='not_estimable',reason=a.availability_status if a.n_selected_instruments==0 else 'no_variant_retained_in_disease_data',statistical_category='not_estimable',method=None,beta=None,se=None,OR=None,CI_lower=None,CI_upper=None,P=None,P_adjusted=None,F_min=None,F_median=None,snps=None)
  else:
   bx=d.beta_exposure.to_numpy();by=d.beta_outcome.to_numpy();sy=d.se_outcome.to_numpy();w=1/sy**2;den=float(np.sum(w*bx**2));b=float(np.sum(w*bx*by)/den);q=float(np.sum(w*(by-b*bx)**2));se=math.sqrt((1 if len(d)==1 else max(1,q/(len(d)-1)))/den);p=math.erfc(abs(b/se)/math.sqrt(2))
   check(a.gene+' '+a.tissue+' '+oc+' single IV expected',len(d)==1)
   r.update(status='estimated',reason=None,method='Wald ratio',beta=b,se=se,OR=math.exp(b),CI_lower=math.exp(b-zcrit*se),CI_upper=math.exp(b+zcrit*se),P=p,P_adjusted=min(1,m*p),F_min=float(d.F.min()),F_median=float(d.F.median()),snps=';'.join(d.snp))
   r['statistical_category']='family_corrected' if p<threshold else 'nominal_only' if p<.05 else 'not_nominally_significant'
   r['direction_vs_blood']='concordant' if (b<0 if a.gene=='TSHR' else b>0) else 'opposite' if b!=0 else 'zero'
  result.append(r)
res=pd.DataFrame(result);res.to_csv(T/'Tissue_MR_50_planned_combinations.csv',index=False)
res[res.family=='primary'].to_csv(T/'Tissue_MR_primary_10.csv',index=False)
res[res.status=='estimated'].to_csv(T/'Tissue_MR_all_estimable.csv',index=False)
check('all 50 planned cells',len(res)==50 and not res.duplicated(['gene','tissue','outcome']).any())
check('10 primary and 40 secondary tests retained',(res.family=='primary').sum()==10 and (res.family=='secondary').sum()==40)
check('all 4 eligible combinations x 5 outcomes estimated',(res.status=='estimated').sum()==20)
check('all 30 unestimable cells retained',(res.status=='not_estimable').sum()==30)
check('all estimates F>10',res.loc[res.status=='estimated','F_min'].gt(10).all())
(O/'MR_checks.json').write_text(json.dumps({'all_pass':all(c['pass'] for c in checks),'checks':checks,'no_threshold_relaxation':True,'n_estimated':20,'n_not_estimable':30,'n_total':50},indent=2))
print(res[res.status=='estimated'][['gene','tissue','outcome','snps','OR','CI_lower','CI_upper','P','P_adjusted','statistical_category','direction_vs_blood']].to_string(index=False))
