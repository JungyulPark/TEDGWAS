from pathlib import Path
import urllib.request,urllib.parse,json,sys,time,hashlib,datetime,gzip,csv,io
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Extension_20261003';S=O/'Sources';D=W/'download';D.mkdir(exist_ok=True)
TI=['Thyroid','Cells_Cultured_fibroblasts','Adipose_Subcutaneous','Adipose_Visceral_Omentum','Whole_Blood']
genes={'ENSG00000165409.16':'TSHR','ENSG00000140443.13':'IGF1R'}
url='https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar'
cache={};log=[]
def fetch(start,n):
    key=f'{start}_{n}';dest=D/(key+'.part')
    if dest.exists():return dest.read_bytes()
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{start+n-1}'})
            with urllib.request.urlopen(req,timeout=90) as response:
                assert response.status==206,(response.status,'server must support ranges')
                b=response.read();cr=response.headers.get('Content-Range')
            assert len(b)==n,(start,n,len(b))
            dest.write_bytes(b)
            log.append({'offset':start,'bytes':n,'content_range':cr,'sha256':hashlib.sha256(b).hexdigest()})
            return b
        except Exception:
            if attempt==2:raise
            time.sleep(2)
offset=0;index=[];selected=[]
indexfile=S/'tar_index.json'
if indexfile.exists():index=json.loads(indexfile.read_text())
else:
 while True:
    block=fetch(offset,512)
    if block==b'\x00'*512:break
    name=block[:100].rstrip(b'\x00').decode();size=int(block[124:136].strip(b'\x00 ') or b'0',8);typ=block[156:157]
    stored=int(block[148:156].strip(b'\x00 '),8);assert stored==sum(block[:148])+8*32+sum(block[156:])
    if typ==b'L':
        longname=fetch(offset+512,size).rstrip(b'\x00').decode();offset+=512+((size+511)//512)*512;block=fetch(offset,512);size=int(block[124:136].strip(b'\x00 ') or b'0',8);name=longname;typ=block[156:157]
    index.append({'name':name,'bytes':size,'offset':offset+512,'type':typ.decode(errors='replace')})
    print('indexed',len(index),name,size,flush=True)
    offset+=512+((size+511)//512)*512
    if len(index)>150:raise RuntimeError('Unexpected archive structure')
 indexfile.write_text(json.dumps(index,indent=2))
for item in index:
 name=Path(item['name']).name
 if not any(name in [t+'.v8.egenes.txt.gz',t+'.v8.signif_variant_gene_pairs.txt.gz'] for t in TI):continue
 dest=D/name
 if not dest.exists():dest.write_bytes(fetch(item['offset'],item['bytes']))
 rows=[]
 with gzip.open(dest,'rt',encoding='utf-8') as f:
    reader=csv.DictReader(f,delimiter='\t');head=reader.fieldnames
    for row in reader:
        if row['gene_id'] in genes:rows.append(row)
 out=S/(name.replace('.txt.gz','')+'.target_genes.tsv')
 with out.open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=head,delimiter='\t');writer.writeheader();writer.writerows(rows)
 selected.append({**item,'local_file':str(dest),'url':url,'archive_generation':'see cloud_objects.json','retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'extract':out.name,'extract_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'rows':len(rows)})
 print('extracted',name,len(rows),flush=True)
(S/'gtex_extraction_manifest.json').write_text(json.dumps(selected,indent=2))
(S/'http_range_log.json').write_text(json.dumps(log,indent=2))
assert len(selected)==10,len(selected)
print('GTEx extraction COMPLETE',flush=True)
