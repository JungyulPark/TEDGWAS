from pathlib import Path
import sys,json,urllib.request,urllib.parse,hashlib,datetime,subprocess,math
import pandas as pd,numpy as np
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Extension_20261003';S=O/'Sources';T=O/'Tables';T.mkdir(exist_ok=True);L=O/'LD';L.mkdir(exist_ok=True)
TI=['Thyroid','Cells_Cultured_fibroblasts','Adipose_Subcutaneous','Adipose_Visceral_Omentum','Whole_Blood']
G={x['geneSymbol']:x for x in json.loads((S/'genes.json').read_text())['data']}
N={x['tissueSiteDetailId']:x['eqtlSampleSummary']['totalCount'] for x in json.loads((S/'tissues.json').read_text())['data']}
avail=[];rows=[];checks=[]
def check(name,ok):
 checks.append({'check':name,'pass':bool(ok)})
 assert ok,name
for tissue in TI:
 eg=pd.read_csv(S/(tissue+'.v8.egenes.target_genes.tsv'),sep='\t');sig=pd.read_csv(S/(tissue+'.v8.signif_variant_gene_pairs.target_genes.tsv'),sep='\t')
 for gene in ['TSHR','IGF1R']:
  meta=G[gene];eid=meta['gencodeId'];e=eg[eg.gene_id.eq(eid)];d=sig[sig.gene_id.eq(eid)].copy()
  check(tissue+' '+gene+' single eGene row',len(e)<=1)
  tested=len(e)==1;minimum=float(e.iloc[0].pval_nominal) if tested else None;threshold=float(e.iloc[0].pval_nominal_threshold) if tested else None
  if tested and minimum<5e-8:check(tissue+' '+gene+' significant release complete at fixed threshold',threshold>=5e-8 and float(e.iloc[0].qval)<=.05 and len(d)>0)
  d['F']=(d.slope/d.slope_se)**2;d=d[(d.pval_nominal<5e-8)&(d.F>10)&(d.maf>=.01)&(d.tss_distance.abs()<=1000000)].copy()
  for _,row in d.iterrows():
   c,pos,ref,alt,build=row.variant_id.split('_');check('build '+row.variant_id,build=='b38')
   check('source TSS '+row.variant_id,int(pos)-int(row.tss_distance)==meta['tss'])
   rows.append({**row.to_dict(),'gene':gene,'tissue':tissue,'n_exposure':N[tissue],'chromosome':int(c[3:]),'pos38':int(pos),'genome_build':'GRCh38','effect_allele':alt,'other_allele':ref,'expression_unit':'inverse-normal-transformed normalized expression','exposure_eaf':None})
  reason='candidate_instruments_available' if len(d) else 'not_tested_in_GTEx_tissue' if not tested else 'no_eQTL_P_below_5e-8' if minimum>=5e-8 else 'no_instrument_after_QC'
  avail.append({'gene':gene,'tissue':tissue,'gencode_id':eid,'n_exposure':N[tissue],'genome_build':'GRCh38','tested':tested,'regional_variants_tested':int(e.iloc[0].num_var) if tested else None,'minimum_eQTL_P':minimum,'gene_qvalue':float(e.iloc[0].qval) if tested else None,'release_nominal_threshold':threshold,'n_GWS_candidates':len(d),'availability_status':reason,'full_regional_summary_available':False,'coloc_status':'not_run_full_regional_statistics_unavailable'})
c=pd.DataFrame(rows);c.to_csv(T/'Candidate_instruments_before_LD.csv',index=False)
api=[];api_dir=S/'variant_mapping';api_dir.mkdir(exist_ok=True)
for chrom,sub in c.groupby('chromosome'):
 pos=sorted(set(sub.pos38));query=[('datasetId','gtex_v8'),('chromosome','chr'+str(chrom)),('itemsPerPage','1000')]+[('pos',str(x)) for x in pos]
 url='https://gtexportal.org/api/v2/dataset/variant?'+urllib.parse.urlencode(query)
 p=api_dir/f'chr{chrom}.json'
 if not p.exists():
  data=urllib.request.urlopen(url,timeout=45).read();p.write_bytes(data)
  p.with_suffix('.provenance.json').write_text(json.dumps({'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(data).hexdigest()},indent=2))
 j=json.loads(p.read_text());check('variant mapping pagination '+str(chrom),j['paging_info']['numberOfPages']<=1)
 api.extend(j['data'])
mapping={v['variantId']:v for v in api};check('all candidate variant IDs resolved',all(v in mapping for v in c.variant_id))
c['snp']=c.variant_id.map(lambda v:mapping[v]['snpId']);c['b37_variant_id']=c.variant_id.map(lambda v:mapping[v]['b37VariantId']);c['pos37']=c.b37_variant_id.str.split('_').str[1].astype(int)
check('no duplicate variant within gene tissue',not c.duplicated(['gene','tissue','variant_id']).any())
ids=set(c.snp);pref=Path('C:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/g1000_eur');ref={}
with pref.with_suffix('.bim').open() as f:
 for line in f:
  cells=line.split()
  if cells[1] in ids:
   check('unique reference SNP '+cells[1],cells[1] not in ref);ref[cells[1]]=cells
c['in_LD_reference']=c.snp.isin(ref)
for _,row in c[c.in_LD_reference].iterrows():
 b=ref[row.snp];check('reference coordinate '+row.snp,int(b[0])==row.chromosome and int(b[3])==row.pos37)
 check('reference allele pair '+row.snp,set(b[4:6])==set(row.b37_variant_id.split('_')[2:4]))
c.to_csv(T/'Candidate_instruments_mapped.csv',index=False)
(L/'extract_snps.txt').write_text('\n'.join(sorted(ids&ref.keys()))+'\n')
plink=Path('C:/ProjectTEDGWAS/tools/plink.exe');small=L/'EUR_candidates'
def run(args,name):
 r=subprocess.run([str(plink)]+args,capture_output=True,text=True,encoding='utf-8',errors='replace');(L/(name+'.console.txt')).write_text(r.stdout+r.stderr,encoding='utf-8');check(name+' PLINK',r.returncode==0)
if not small.with_suffix('.bed').exists():run(['--bfile',str(pref),'--extract',str(L/'extract_snps.txt'),'--keep-allele-order','--make-bed','--threads','2','--memory','8192','--out',str(small)],'extract')
n=sum(1 for _ in small.with_suffix('.fam').open());width=(n+3)//4;genos={};bim={}
with small.with_suffix('.bed').open('rb') as f:
 check('BED header',f.read(3)==b'\x6c\x1b\x01')
 for line in small.with_suffix('.bim').open():
  b=line.split();data=f.read(width);genos[b[1]]=np.array([{0:2.,1:np.nan,2:1.,3:0.}[(data[i//4]>>(2*(i%4)))&3] for i in range(n)]);bim[b[1]]=b
selected=[];ldpairs=[]
for (gene,tissue),sub in c.groupby(['gene','tissue'],sort=False):
 eligible=sub[sub.in_LD_reference].sort_values(['pval_nominal','snp']);stem=gene+'__'+tissue;inp=L/(stem+'.clump_input.tsv')
 eligible[['snp','pval_nominal']].rename(columns={'snp':'SNP','pval_nominal':'P'}).to_csv(inp,sep='\t',index=False)
 run(['--bfile',str(small),'--clump',str(inp),'--clump-p1','5e-8','--clump-p2','5e-8','--clump-r2','0.001','--clump-kb','10000','--threads','2','--memory','1024','--out',str(L/stem)],stem)
 out=pd.read_csv(L/(stem+'.clumped'),sep=r'\s+');plink_snps=list(out.SNP);manual=[]
 for snp in eligible.snp:
  blocked=False
  for kept in manual:
   x,y=genos[snp],genos[kept];ok=np.isfinite(x)&np.isfinite(y);rr=float(np.corrcoef(x[ok],y[ok])[0,1]);r2=rr*rr
   if abs(int(bim[snp][3])-int(bim[kept][3]))<=10000000 and r2>=.001:blocked=True;break
  if not blocked:manual.append(snp)
 check(stem+' independent greedy clump matches PLINK',set(manual)==set(plink_snps))
 chosen=eligible[eligible.snp.isin(plink_snps)].copy();chosen['n_reference']=n
 for _,row in chosen.iterrows():
  a=genos[row.snp];b=bim[row.snp];chosen.loc[row.name,'reference_effect_allele_frequency']=float(np.nanmean(a)/2 if row.effect_allele==b[4] else 1-np.nanmean(a)/2)
 for i,s1 in enumerate(plink_snps):
  for s2 in plink_snps[i+1:]:
   x,y=genos[s1],genos[s2];ok=np.isfinite(x)&np.isfinite(y);rr=float(np.corrcoef(x[ok],y[ok])[0,1]);ldpairs.append({'gene':gene,'tissue':tissue,'snp1':s1,'snp2':s2,'r2':rr*rr,'n_complete':int(ok.sum())});check(stem+' retained r2',rr*rr<.001)
 selected.append(chosen)
 print(gene,tissue,'candidates',len(sub),'reference',len(eligible),'selected',','.join(plink_snps),flush=True)
selected=pd.concat(selected,ignore_index=True);selected.to_csv(T/'Selected_instruments.csv',index=False)
pd.DataFrame(ldpairs,columns=['gene','tissue','snp1','snp2','r2','n_complete']).to_csv(T/'Retained_instrument_LD.csv',index=False)
for row in avail:
 row['n_LD_reference_candidates']=int(((c.gene==row['gene'])&(c.tissue==row['tissue'])&c.in_LD_reference).sum());row['n_selected_instruments']=int(((selected.gene==row['gene'])&(selected.tissue==row['tissue'])).sum())
pd.DataFrame(avail).to_csv(T/'Tissue_availability_10_combinations.csv',index=False)
(O/'Instrument_checks.json').write_text(json.dumps({'all_pass':all(x['pass'] for x in checks),'checks':checks,'n_reference':n},indent=2))
