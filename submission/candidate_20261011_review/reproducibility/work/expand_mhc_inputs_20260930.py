from pathlib import Path
import gzip,csv,json,time,hashlib
W=Path(__file__).resolve().parent;T=W/'expanded_analysis_20260930';I=T/'inputs';I.mkdir(parents=True,exist_ok=True)
L=Path('C:/ProjectTEDGWAS');V=L/'TrackA_MR/v5_upgrade';genes={b'HLA-A',b'HLA-DQA2',b'C4A',b'TUBB',b'PSMB8'}
def log(s):
 print(s,flush=True)
 with (T/'analysis.log').open('a',encoding='utf8') as f:f.write(time.strftime('%Y-%m-%d %H:%M:%S ')+s+'\n')
def extract(path,out,selector,header_fn):
 meta=out.with_suffix('.json')
 if out.exists() and meta.exists():log('Reuse '+out.name);return
 log('Read '+str(path));n=k=0;start=time.monotonic()
 op=gzip.open if path.suffix=='.gz' else open
 with op(path,'rb') as f,out.with_suffix('.partial').open('wb') as dest:
  header=f.readline().decode().strip().split();outheader,convert=header_fn(header);dest.write(('\t'.join(outheader)+'\n').encode())
  for line in f:
   n+=1;c=line.rstrip(b'\r\n').split(b'\t') if path.suffix=='.gz' else line.split()
   if selector(c,header):dest.write(b'\t'.join(convert(c))+b'\n');k+=1
   if n%5000000==0:log(f'{out.name}: {n:,} scanned; {k:,} selected')
 out.with_suffix('.partial').rename(out)
 m={'source':str(path),'source_bytes':path.stat().st_size,'source_mtime_ns':path.stat().st_mtime_ns,'rows_scanned':n,'rows_selected':k,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start,'gzip_full_read':path.suffix=='.gz'}
 meta.write_text(json.dumps(m,indent=2),encoding='utf8');log('Completed '+out.name+': '+str(k))
eqpath=L/'TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz'
extract(eqpath,I/'eqtl_mhc.tsv',lambda c,h:c[8] in genes,lambda h:(h,lambda c:c))
eq=list(csv.DictReader((I/'eqtl_mhc.tsv').open(),delimiter='\t'));keys={r['SNP'].encode() for r in eq};assert {r['GeneSymbol'] for r in eq}=={g.decode() for g in genes}
freq=V/'04_druggable_mr/data/g1000_eur_freq.frq'
extract(freq,I/'eur_freq.tsv',lambda c,h:c[1] in keys,lambda h:(h,lambda c:c))
maps={'snp':['hm_rsid','rsid','rsids','variant_id','snp','snpid'],'chr':['chromosome','#chrom','chrom','hm_chrom'],'pos':['base_pair_location','pos','hm_pos'],'ea':['effect_allele','alt','hm_effect_allele'],'oa':['other_allele','ref','hm_other_allele'],'beta':['beta','hm_beta','effect'],'se':['standard_error','sebeta','se'],'eaf':['effect_allele_frequency','af_alt','eaf','maf'],'pvalue':['p_value','pval','p','pvalue']}
for name,path in {'BBJ':V/'04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz','UKB':V/'04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz','FinnGen':L/'finngen_R12_GRAVES_OPHT.gz'}.items():
 with gzip.open(path,'rt') as f:h=next(f).strip().split('\t')
 lower=[x.lower() for x in h];ix={k:next(lower.index(v) for v in vals if v in lower) for k,vals in maps.items()}
 extract(path,I/(name+'.tsv'),lambda c,h:c[ix['snp']] in keys,lambda h:(list(maps)+['genome_build'],lambda c:[c[ix[k]] for k in maps]+[b'GRCh38']))
af=W/'data/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz'
extract(af,I/'eqtl_af.tsv',lambda c,h:c[0] in keys,lambda h:(h,lambda c:c))
log('All five-gene input extractions complete; source files unchanged.')
