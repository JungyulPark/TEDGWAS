from pathlib import Path
import json,sys,hashlib
W=Path(__file__).resolve().parent;sys.path.insert(0,str(W/'pydeps'))
import pandas as pd,numpy as np
from scipy.special import logsumexp
T=W/'expanded_analysis_20260930';I=T/'inputs';O=T/'results';O.mkdir(exist_ok=True)
D=W/'review_8dce23c/submission/candidate_20260911_review';old=W/'maf_sensitivity_20260905/inputs'
inst=pd.read_csv(D/'Supplementary_Data_1_Instruments.csv');mr=pd.read_csv(D/'Supplementary_Data_2_MR.csv');mr=mr[mr.scenario.eq('original_reference')].copy()
assert not mr.duplicated(['gene_symbol','outcome']).any()
genes=sorted(inst.gene_symbol.unique());assert len(genes)==2544
mhc=['HLA-A','HLA-DQA2','C4A','TUBB','PSMB8'];ocs={'BBJ_Graves':'BBJ','UKB_hyperthyroid':'UKB','FinnGen_GO':'FinnGen'}
full=pd.MultiIndex.from_product([genes,list(ocs)],names=['gene_symbol','outcome']).to_frame(index=False).merge(mr,on=['gene_symbol','outcome'],how='left',validate='one_to_one')
full['estimable']=full.pvalue.notna();full['nominal_p_lt_0_05']=full.pvalue.lt(.05);full['BBJ_discovery']=full.outcome.eq('BBJ_Graves')&full.pvalue.lt(.05/2544)
full['status']=np.where(full.estimable,'Estimated','Not estimable; not a negative finding');full.to_csv(O/'Stage1_all_2544_genes_3_outcomes.csv',index=False)
hits=full[full.BBJ_discovery].sort_values('pvalue').gene_symbol.tolist();assert len(hits)==13
disc=full[full.gene_symbol.isin(hits)].copy();disc['gene_order']=disc.gene_symbol.map({g:i for i,g in enumerate(hits)});disc.sort_values(['gene_order','outcome']).drop(columns='gene_order').to_csv(O/'Stage1_13_discovery_genes.csv',index=False)
existing=pd.read_csv(old/'eqtlgen_9genes.tsv',sep='\t');extra=pd.read_csv(I/'eqtl_mhc.tsv',sep='\t')
assert set(extra.GeneSymbol)==set(mhc)
assert extra.SNPChr.eq(6).all() and extra.GeneChr.eq(6).all()
assert (extra.SNPPos-extra.GenePos).abs().max()<=1000000
assert not extra.duplicated(['GeneSymbol','SNP']).any()
qc=extra.groupby('GeneSymbol').agg(n_eqtl=('SNP','size'),min_p=('Pvalue','min'),max_p=('Pvalue','max'),gene_pos_hg19=('GenePos','first'),min_pos_hg19=('SNPPos','min'),max_pos_hg19=('SNPPos','max'),ensembl_ids=('Gene',lambda x:';'.join(sorted(set(x)))))
assert qc.max_p.gt(.9).all();qc.to_csv(O/'MHC_input_coverage.csv')
rows=[];coverage=[];allgenes=list(existing.GeneSymbol.unique())+mhc
vdir=I/'r_validation';vdir.mkdir(exist_ok=True)
def abf(z,v,sd):
 r=sd**2/(sd**2+v)
 return .5*(np.log1p(-r)+r*z*z)
for group,eq,ref,af in [('existing',existing,pd.read_csv(old/'eur_freq_selected.tsv',sep='\t'),pd.read_csv(old/'eqtlgen_af_selected_v1.csv')),('mhc',extra,pd.read_csv(I/'eur_freq.tsv',sep='\t'),pd.read_csv(I/'eqtl_af.tsv',sep='\t').rename(columns={'SNP':'snp','hg19_chr':'chr','hg19_pos':'pos','AlleleA':'allele_a','AlleleB':'allele_b','AlleleB_all':'af_b'}))]:
 assert not ref.SNP.duplicated().any() and not af.snp.duplicated().any()
 if group=='mhc':
  calc=(af.allAB_total+2*af.allB_total)/(2*(af.allA_total+af.allAB_total+af.allB_total))
  assert np.max(np.abs(calc-af.af_b))<1e-12
  coords=eq[['SNP','SNPChr','SNPPos']].drop_duplicates().merge(af,left_on='SNP',right_on='snp')
  assert coords.SNPChr.eq(coords.chr).all() and coords.SNPPos.eq(coords.pos).all()
 e=eq.merge(ref[['SNP','A1','A2','MAF']],on='SNP',how='left',validate='many_to_one').merge(af[['snp','allele_a','allele_b','af_b']],left_on='SNP',right_on='snp',how='left',validate='many_to_one')
 e['eaf_ref']=np.where(e.AssessedAllele==e.A1,e.MAF,np.where(e.AssessedAllele==e.A2,1-e.MAF,np.nan))
 ok=((e.AssessedAllele==e.allele_b)&(e.OtherAllele==e.allele_a))|((e.AssessedAllele==e.allele_a)&(e.OtherAllele==e.allele_b))
 e['eaf_new']=np.where(ok,np.where(e.AssessedAllele==e.allele_b,e.af_b,1-e.af_b),np.nan)
 for kind in ['ref','new']:
  f=e['eaf_'+kind];e['valid_'+kind]=f.between(0,1,inclusive='neither');e['se_'+kind]=1/np.sqrt(2*f*(1-f)*(e.NrSamples+e.Zscore**2))
 for oc,label in ocs.items():
  d=pd.read_csv(old/f'{label}_selected.tsv' if group=='existing' else I/f'{label}.tsv',sep='\t')
  d=d[d.beta.notna()&d.se.gt(0)&d.snp.notna()].copy()
  if label=='UKB':mu=3731/484598;d[['beta','se']]=d[['beta','se']]/(mu*(1-mu))
  for g in eq.GeneSymbol.unique():
   ge=e[e.GeneSymbol.eq(g)];ch=int(ge.SNPChr.iloc[0]);m=ge.merge(d[d.chr.eq(ch)],left_on='SNP',right_on='snp',suffixes=('_e','_o'))
   same=(m.AssessedAllele==m.ea)&(m.OtherAllele==m.oa);swap=(m.AssessedAllele==m.oa)&(m.OtherAllele==m.ea)
   before=len(m);m=m[same|swap].copy();dup=m.SNP.duplicated().sum();m=m.drop_duplicates('SNP').sort_values('SNP')
   a=m[m.valid_ref];b=m[m.valid_new];common=m[m.valid_ref&m.valid_new]
   coverage.append(dict(gene=g,outcome=oc,n_eqtl=len(ge),n_matched_before_alleles=before,n_allele_compatible=len(m),duplicates_removed=int(dup),n_ref=len(a),n_new=len(b),n_common=len(common)))
   for data,scenario,kind in [(a,'original_reference','ref'),(common,'paired_reference','ref'),(common,'paired_eqtlgen','new'),(b,'available_eqtlgen','new')]:
    assert len(data)>=50
    z1=data.beta.to_numpy()/data.se.to_numpy();v1=data.se.to_numpy()**2;z2=data.Zscore.to_numpy();v2=data['se_'+kind].to_numpy()**2
    assert np.isfinite(np.r_[z1,v1,z2,v2]).all()
    l1=abf(z1,v1,.2);l2=abf(z2,v2,.15);s1=logsumexp(l1);s2=logsumexp(l2);ss=logsumexp(l1+l2)
    diff=ss-s1-s2;assert diff<=1e-10
    h3=s1+s2+np.log(-np.expm1(min(diff,0)))
    pd.DataFrame({'snp':data.SNP.to_numpy(),'beta1':data.beta.to_numpy(),'varbeta1':v1,'beta2':z2*np.sqrt(v2),'varbeta2':v2}).to_csv(vdir/f'{g}__{oc}__{scenario}.csv',index=False)
    for p12 in [1e-6,5e-6,1e-5]:
     lh=np.array([0,np.log(1e-4)+s1,np.log(1e-4)+s2,2*np.log(1e-4)+h3,np.log(p12)+ss]);pp=np.exp(lh-logsumexp(lh));assert abs(pp.sum()-1)<1e-10
     rows.append(dict(gene=g,outcome=oc,scenario=scenario,p12=p12,n_overlap=len(data),top_snp=data.SNP.iloc[int(np.argmax(l1+l2))],**{f'PP.H{i}':pp[i] for i in range(5)}))
got=pd.DataFrame(rows);assert len(got)==14*3*4*3;got.to_csv(O/'Stage2_all_504_colocalization_settings.csv',index=False);pd.DataFrame(coverage).to_csv(O/'Stage2_variant_coverage.csv',index=False)
frozen=pd.read_csv(D/'Supplementary_Data_3_Colocalization.csv');match=frozen.merge(got,on=['gene','outcome','scenario','p12'],suffixes=('_old','_new'),validate='one_to_one');assert len(match)==324
delta=max((match[f'PP.H{i}_old']-match[f'PP.H{i}_new']).abs().max() for i in range(5));assert delta<1e-9
assert match.n_overlap_old.eq(match.n_overlap_new).all() and match.top_snp_old.eq(match.top_snp_new).all()
primary=got[got.scenario.eq('original_reference')&got.p12.eq(1e-5)].copy();primary.to_csv(O/'Stage2_14_genes_primary.csv',index=False)
summary=[]
for g in hits+['IGF1R']:
 x=got[got.gene.eq(g)];r={'gene':g,'BBJ_discovery':g in hits,'MHC':g in mhc}
 for oc,label in ocs.items():
  z=x[x.outcome.eq(oc)];base=z[z.scenario.eq('original_reference')];r[label+'_H4']=float(base[base.p12.eq(1e-5)]['PP.H4'].iloc[0]);r[label+'_H4_low_prior']=float(base[base.p12.eq(1e-6)]['PP.H4'].iloc[0]);r[label+'_H4_eqtlgen']=float(z[z.scenario.eq('available_eqtlgen')&z.p12.eq(1e-5)]['PP.H4'].iloc[0])
 r['combined_primary']=r['BBJ_H4']>=.8 and r['FinnGen_H4']>=.8;r['combined_low_prior']=r['BBJ_H4_low_prior']>=.8 and r['FinnGen_H4_low_prior']>=.8;summary.append(r)
pd.DataFrame(summary).to_csv(O/'Stage2_integrated_summary.csv',index=False)
report={'status':'Python computation PASS; independent R verification pending','genes':14,'settings':len(got),'existing_324_max_abs_delta':float(delta),'stage1_rows':len(full),'stage1_genes':len(genes),'discovery_genes':hits,'combined_primary':[r['gene'] for r in summary if r['combined_primary']],'combined_low_prior':[r['gene'] for r in summary if r['combined_low_prior']],'MHC_eqtl_region':'All reported cis variants within 1 Mb of the source GenePos (GRCh37); no significance filtering','posthoc':'Five MHC genes added after review; no prospective selection claim'}
(O/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2));print(pd.DataFrame(summary).to_string(index=False))
