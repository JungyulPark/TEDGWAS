from pathlib import Path
import json,hashlib,gzip,time,math,csv
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
O=ROOT/'outputs/TED_Data_Review_20261002';O.mkdir(parents=True,exist_ok=True)
L=Path('C:/ProjectTEDGWAS');D=ROOT/'work/expanded_analysis_20260930/results'
I=ROOT/'work/review_8dce23c/submission/candidate_20260911_review/Supplementary_Data_1_Instruments.csv'
def log(s):
 print(s,flush=True)
 with (O/'analysis_log.txt').open('a',encoding='utf8') as f:f.write(time.strftime('%Y-%m-%d %H:%M:%S ')+s+'\n')
def save(name,obj):(O/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
inst=pd.read_csv(I);snps=set(inst.snp);qc=[];examples=[];top=[]
paths={
'BBJ':L/'TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz',
'UKB':L/'TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz',
'FinnGen':L/'finngen_R12_GRAVES_OPHT.gz'}
maps={
'BBJ':dict(snp='rsid',chrom='chromosome',pos='base_pair_location',ea='effect_allele',oa='other_allele',beta='beta',se='standard_error',p='p_value',eaf='effect_allele_frequency'),
'UKB':dict(snp='hm_rsid',chrom='chromosome',pos='base_pair_location',ea='effect_allele',oa='other_allele',beta='beta',se='standard_error',p='p_value',eaf='effect_allele_frequency'),
'FinnGen':dict(snp='rsids',chrom='#chrom',pos='pos',ea='alt',oa='ref',beta='beta',se='sebeta',p='pval',eaf='af_alt')}
for name,path in paths.items():
 log('Full-file QC starts: '+name)
 m=maps[name];start=time.monotonic();before=path.stat();headers=list(pd.read_csv(path,sep='\t',nrows=0).columns)
 counts=dict(rows=0,missing_rsid=0,missing_allele=0,invalid_beta=0,invalid_se=0,invalid_p=0,p_zero=0,invalid_eaf=0,p_lt_0_05=0,p_lt_5e_8=0,valid_core_rows=0,non_snv_rows=0)
 missing={c:0 for c in headers};chrs={};hp={};sample_min=math.inf;sample_max=0;found=set();tops=[];first=None;case_fracs=[];ukb_coord_diff=0
 for chunk in pd.read_csv(path,sep='\t',dtype=str,chunksize=300000,low_memory=False):
  counts['rows']+=len(chunk)
  for c in headers:missing[c]+=int(chunk[c].isna().sum())
  x=chunk.rename(columns={v:k for k,v in m.items()})
  b=pd.to_numeric(x.beta,errors='coerce');se=pd.to_numeric(x.se,errors='coerce');p=pd.to_numeric(x.p,errors='coerce');f=pd.to_numeric(x.eaf,errors='coerce')
  bp=np.isfinite(b);sp=np.isfinite(se)&se.gt(0);pp=np.isfinite(p)&p.between(0,1);fp=np.isfinite(f)&f.between(0,1)
  counts['missing_rsid']+=int(x.snp.isna().sum());counts['missing_allele']+=int((x.ea.isna()|x.oa.isna()).sum())
  counts['invalid_beta']+=int((~bp).sum());counts['invalid_se']+=int((~sp).sum());counts['invalid_p']+=int((~pp).sum());counts['invalid_eaf']+=int((~fp).sum())
  counts['p_zero']+=int(p.eq(0).sum());counts['p_lt_0_05']+=int((pp&p.lt(.05)).sum());counts['p_lt_5e_8']+=int((pp&p.lt(5e-8)).sum())
  counts['valid_core_rows']+=int((bp&sp&pp&x.ea.notna()&x.oa.notna()).sum());counts['non_snv_rows']+=int((x.ea.str.len().ne(1)|x.oa.str.len().ne(1)).sum())
  for c,n in x.chrom.value_counts().items():chrs[c]=chrs.get(c,0)+int(n)
  if p[pp].size:sample_min=min(sample_min,float(p[pp].min()));sample_max=max(sample_max,float(p[pp].max()))
  q=pd.DataFrame({'cohort':name,'rsid':x.snp,'chrom':x.chrom,'position':x.pos,'effect_allele':x.ea,'other_allele':x.oa,'beta':b,'se':se,'p':p,'eaf':f})
  tops.append(q.loc[p[pp].nsmallest(20).index]);tops=[pd.concat(tops).nsmallest(20,'p')]
  if first is None:first=chunk.head(3).fillna('').to_dict('records')
  keep=x.snp.isin(['rs179252','rs2654980','rs13030124','rs148849825'])
  if keep.any():examples.extend(q[keep].to_dict('records'))
  found.update(x.snp[x.snp.isin(snps)].tolist())
  if name=='UKB':ukb_coord_diff+=int((chunk['hm_pos'].notna()&chunk['base_pair_location'].notna()&chunk['hm_pos'].ne(chunk['base_pair_location'])).sum())
  if name=='FinnGen' and len(case_fracs)<20000:
   fc=pd.to_numeric(chunk.af_alt_cases,errors='coerce');fu=pd.to_numeric(chunk.af_alt_controls,errors='coerce');den=fc-fu;ok=den.abs().gt(.1)&f.notna()
   case_fracs.extend(((f[ok]-fu[ok])/den[ok]).tolist())
  if counts['rows']%3000000==0:log(name+f': {counts["rows"]:,} rows')
 after=path.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'source changed while reading'
 rec={'cohort':name,'source':str(path),'source_bytes':before.st_size,'sha256':digest(path),**counts,'p_min':sample_min,'p_max':sample_max,'columns':headers,'missing_by_column':missing,'chromosome_counts':chrs,'selected_unique_instrument_rsids_present':len(found),'selected_unique_instrument_rsids_total':len(snps),'inspection':'all rows; compressed stream fully read; duplicates/LD independence not assessed','elapsed_seconds':round(time.monotonic()-start,1),'first_rows':first}
 if name=='UKB':rec['hm_vs_original_position_differences']=ukb_coord_diff
 if case_fracs:rec['case_fraction_from_af_median']=float(np.median(case_fracs));rec['case_fraction_af_rows']=len(case_fracs)
 qc.append(rec);top.extend(tops[0].to_dict('records'));save('Raw_GWAS_QC.json',qc);log(name+f' completed {counts["rows"]:,} rows, {rec["elapsed_seconds"]} sec')
pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in qc]).to_csv(O/'Raw_GWAS_QC.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(examples).to_csv(O/'Example_variants.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(top).to_csv(O/'Top20_variant_rows_each_cohort.csv',index=False,encoding='utf-8-sig')
# Re-tabulate and verify the published gene-level analysis; this does not rerun MR from SNPs.
mr=pd.read_csv(D/'Stage1_all_2544_genes_3_outcomes.csv');coloc=pd.read_csv(D/'Stage2_14_genes_primary.csv')
assert len(mr)==2544*3 and not mr.duplicated(['gene_symbol','outcome']).any()
est=mr[mr.pvalue.notna()].copy();calc_p=np.array([math.erfc(abs(z)/math.sqrt(2)) for z in (est.beta/est.se)])
errs={'p_max_abs_delta':float(np.max(np.abs(calc_p-est.pvalue))),'log_or_max_abs_delta':float(np.max(np.abs(np.log(est['or'])-est.beta))),'ci_lower_log_max_abs_delta':float(np.max(np.abs(np.log(est.or_ci_lower)-(est.beta-1.96*est.se)))),'ci_upper_log_max_abs_delta':float(np.max(np.abs(np.log(est.or_ci_upper)-(est.beta+1.96*est.se))))}
assert errs['p_max_abs_delta']<1e-8 and errs['log_or_max_abs_delta']<1e-8
summary=[]
for oc,g in mr.groupby('outcome',sort=False):
 e=g[g.pvalue.notna()];summary.append({'outcome':oc,'eligible_genes':len(g),'estimable_genes':len(e),'not_estimable_genes':int(g.pvalue.isna().sum()),'p_lt_0_05':int(e.pvalue.lt(.05).sum()),'p_lt_0_05_over_2544_descriptive':int(e.pvalue.lt(.05/2544).sum()),'single_instrument_genes':int(e.n_iv.eq(1).sum()),'multiple_instrument_genes':int(e.n_iv.gt(1).sum()),'median_se':float(e.se.median()),'median_n_iv':float(e.n_iv.median()),'max_n_iv':int(e.n_iv.max())})
hits=mr[(mr.outcome=='BBJ_Graves')&mr.pvalue.lt(.05/2544)].sort_values('pvalue').gene_symbol.tolist();assert len(hits)==13
selected=mr[mr.gene_symbol.isin(hits+['IGF1R'])].merge(coloc,left_on=['gene_symbol','outcome'],right_on=['gene','outcome'],validate='one_to_one')
selected['stage1_rank']=selected.gene_symbol.map({g:i for i,g in enumerate(hits+['IGF1R'])});selected=selected.sort_values(['stage1_rank','outcome'])
pd.DataFrame(summary).to_csv(O/'MR_cohort_summary.csv',index=False,encoding='utf-8-sig');selected.to_csv(O/'Selected_14_gene_results.csv',index=False,encoding='utf-8-sig');mr.to_csv(O/'All_2544_genes_3_outcomes.csv',index=False,encoding='utf-8-sig')
piv=mr.pivot(index='gene_symbol',columns='outcome',values='beta');direction={}
for a,b in [('BBJ_Graves','UKB_hyperthyroid'),('BBJ_Graves','FinnGen_GO'),('UKB_hyperthyroid','FinnGen_GO')]:
 v=piv[[a,b]].dropna();h=v.loc[v.index.isin(hits)];direction[a+'__'+b]={'all_common_n':len(v),'all_same_direction':int(((v[a]*v[b])>0).sum()),'bbj_hits_common_n':len(h),'bbj_hits_same_direction':int(((h[a]*h[b])>0).sum())}
report={'date':'2026-10-02','scope':'full-file descriptive QC + re-tabulation/consistency checks of existing MR; no new SNP-level MR or colocalization','raw_qc':qc,'mr_summary':summary,'discovery_genes':hits,'mr_numeric_checks':errs,'direction_counts_descriptive_only':direction,'instrument_records':len(inst),'instrument_genes':inst.gene_symbol.nunique(),'instrument_unique_snps':len(snps),'instrument_n_min':int(inst.n_samples.min()),'instrument_n_median':float(inst.n_samples.median()),'instrument_n_max':int(inst.n_samples.max()),'source_hashes':{str(p):digest(p) for p in [I,D/'Stage1_all_2544_genes_3_outcomes.csv',D/'Stage2_14_genes_primary.csv']}}
save('analysis_summary.json',report)
log('Gene re-tabulation complete. '+json.dumps(summary));log('Consistency checks: '+json.dumps(errs))
