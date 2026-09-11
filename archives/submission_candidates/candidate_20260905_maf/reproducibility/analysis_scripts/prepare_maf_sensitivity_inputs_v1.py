from pathlib import Path
import csv,gzip,hashlib,json,time,datetime,shutil
W=Path(__file__).resolve().parent;T=W/'maf_sensitivity_20260905';I=T/'inputs';I.mkdir(parents=True,exist_ok=True)
P=W/'revision/submission/candidate_20260905/provenance'
L=Path('C:/ProjectTEDGWAS');V=L/'TrackA_MR/v5_upgrade'
coloc=W.parent/'outputs/TED_TRAP_submission_candidate_20260905/reproducibility/inputs'
def log(s):
    line=datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s
    print(line,flush=True)
    with (T/'preparation.log').open('a',encoding='utf-8') as f:f.write(line+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()
inst=list(csv.DictReader((P/'instruments_verified.csv').open(encoding='utf-8')))
eq=list(csv.DictReader((coloc/'eqtlgen_9genes.tsv').open(encoding='utf-8'),delimiter='\t'))
keys={r['snp'].encode() for r in inst}|{r['SNP'].encode() for r in eq}
log(f'Start: {len(inst)} instrument rows; {len(keys)} unique requested rsIDs including coloc variants.')
shutil.copy2(P/'instruments_verified.csv',I/'instruments_verified.csv')
shutil.copy2(coloc/'eqtlgen_9genes.tsv',I/'eqtlgen_9genes.tsv')
sources={'BBJ':V/'04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz','UKB':V/'04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz','FinnGen':L/'finngen_R12_GRAVES_OPHT.gz'}
records=[]
maps={'snp':['hm_rsid','rsid','rsids','variant_id','snp','snpid'],'chr':['chromosome','#chrom','chrom','hm_chrom'],'pos':['base_pair_location','pos','hm_pos'],'ea':['effect_allele','alt','hm_effect_allele'],'oa':['other_allele','ref','hm_other_allele'],'beta':['beta','hm_beta','effect'],'se':['standard_error','sebeta','se'],'eaf':['effect_allele_frequency','af_alt','eaf','maf'],'pvalue':['p_value','pval','p','pvalue']}
for name,path in sources.items():
    out=I/(name+'_selected.tsv');meta=I/(name+'_selected.meta.json')
    if out.exists() and meta.exists():log(name+' existing complete extraction reused');records.append(json.loads(meta.read_text()));continue
    log(name+' reading '+str(path));start=last=time.monotonic();n=kept=0
    with gzip.open(path,'rb') as f,out.with_suffix('.partial').open('wb') as dest:
        headers=f.readline().decode().strip().split('\t');lower=[v.lower() for v in headers]
        mapping={k:next((lower.index(v) for v in names if v in lower),None) for k,names in maps.items()}
        assert all(x is not None for x in mapping.values()),(name,headers,mapping)
        dest.write(('\t'.join(mapping)+'\tgenome_build\n').encode())
        for line in f:
            n+=1;cols=line.rstrip(b'\r\n').split(b'\t')
            if cols[mapping['snp']] in keys:
                dest.write(b'\t'.join(cols[j] for j in mapping.values())+b'\tGRCh38\n');kept+=1
            if n%1000000==0 and time.monotonic()-last>20:log(f'{name}: scanned {n:,}; selected {kept:,}');last=time.monotonic()
    out.with_suffix('.partial').rename(out)
    record={'dataset':name,'source_path':str(path),'source_bytes':path.stat().st_size,'source_sha256':sha(path),'source_header':headers,'column_mapping':{k:headers[v] for k,v in mapping.items()},'selection':'exact rsID membership; preserves source row order and alleles; no coordinate filtering','rows_scanned':n,'rows_selected':kept,'output':out.name,'output_sha256':sha(out),'seconds':round(time.monotonic()-start,2)}
    meta.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');records.append(record);log(f'{name}: completed {n:,} -> {kept:,} rows')
freq=V/'04_druggable_mr/data/g1000_eur_freq.frq';out=I/'eur_freq_selected.tsv';meta=I/'eur_freq_selected.meta.json'
if out.exists() and meta.exists():records.append(json.loads(meta.read_text()));log('EUR frequency complete extraction reused')
else:
    log('EUR reference frequency extraction');n=kept=0
    with freq.open('rb') as f,out.with_suffix('.partial').open('wb') as dest:
        header=f.readline().split();ix=header.index(b'SNP');dest.write(b'\t'.join(header)+b'\tgenome_build\n')
        for line in f:
            n+=1;cols=line.split()
            if cols[ix] in keys:dest.write(b'\t'.join(cols)+b'\tGRCh37\n');kept+=1
    out.with_suffix('.partial').rename(out)
    record={'dataset':'1000G_EUR','source_path':str(freq),'source_bytes':freq.stat().st_size,'source_sha256':sha(freq),'rows_scanned':n,'rows_selected':kept,'output':out.name,'output_sha256':sha(out)}
    meta.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');records.append(record);log(f'EUR reference: {n:,} -> {kept:,} rows')
(T/'input_manifest.json').write_text(json.dumps({'created_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'local_only':True,'instruments_rows':len(inst),'instrument_sha256':sha(I/'instruments_verified.csv'),'coloc_eqtl_sha256':sha(I/'eqtlgen_9genes.tsv'),'sources':records},indent=2)+'\n',encoding='utf-8')
log('Local input preparation completed. No source file was modified.')
