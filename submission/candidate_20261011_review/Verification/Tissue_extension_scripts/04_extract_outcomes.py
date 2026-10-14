from pathlib import Path
import sys,gzip,re,json,hashlib,datetime,csv
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Extension_20261003';S=O/'Sources';T=O/'Tables';D=S/'disease_extracts';D.mkdir(exist_ok=True)
iv=pd.read_csv(T/'Selected_instruments.csv')
for idx,row in iv.iterrows():
 e=pd.read_csv(S/(row.tissue+'.v8.egenes.target_genes.tsv'),sep='\t');e=e[(e.gene_id==row.gene_id)&(e.variant_id==row.variant_id)]
 assert len(e)==1
 ef=e.iloc[0];assert float(ef.maf)==row.maf
 assert int(ef.ref_factor) in [1,-1]
 iv.loc[idx,'exposure_eaf']=row.maf if ef.ref_factor==1 else 1-row.maf
 iv.loc[idx,'exposure_eaf_source']='GTEx tissue eGene MAF and ref_factor'
 assert abs(row.ma_count/(2*row.n_exposure)-row.maf)<1e-6
 assert row.effect_allele==row.variant_id.split('_')[3]
 assert set([row.effect_allele,row.other_allele]) not in [set('AT'),set('CG')]
iv.to_csv(T/'Selected_instruments_verified.csv',index=False)
targets={(str(r.chromosome),str(r.pos38)):r.snp for _,r in iv.iterrows()};snps=set(iv.snp)
pattern=re.compile(rb'\t(?:'+b'|'.join(str(p).encode() for p in set(iv.pos38))+rb')\t')
src=Path('C:/ProjectTEDGWAS')
specs=[('BBJ_Graves',src/'TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz','chromosome','base_pair_location','rsid'),('UKB_hyperthyroid',src/'TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz','hm_chrom','hm_pos','hm_rsid'),('FinnGen_GO',src/'finngen_R12_GRAVES_OPHT.gz','#chrom','pos','rsids'),('E4_GRAVES_STRICT',src/'finngen_R12_E4_GRAVES_STRICT.gz','#chrom','pos','rsids'),('E4_GRAVES_OPHT_STRICT',src/'finngen_R12_E4_GRAVES_OPHT_STRICT.gz','#chrom','pos','rsids')]
prior=json.loads((R/'outputs/TED_Allele_Direction_20261003/SOURCE_MANIFEST.json').read_text())
prior={Path(p['path']).resolve():p for p in prior}
class HashReader:
 def __init__(self,f):self.f=f;self.h=hashlib.sha256()
 def read(self,n=-1):
  b=self.f.read(n);self.h.update(b);return b
 def tell(self):return self.f.tell()
def extract(spec):
 name,path,chrcol,poscol,idcol=spec;dest=D/(name+'.tsv');manifest=D/(name+'.manifest.json')
 if manifest.exists():
  m=json.loads(manifest.read_text());assert hashlib.sha256(dest.read_bytes()).hexdigest()==m['extract_sha256'];return m
 found=[];rows=0
 with path.open('rb') as raw,dest.open('wb') as out:
  hashed=HashReader(raw)
  with gzip.GzipFile(fileobj=hashed) as f:
   header=f.readline();out.write(header);fields=header.decode().strip().split('\t');ci=fields.index(chrcol);pi=fields.index(poscol);si=fields.index(idcol)
   for line in f:
    rows+=1
    if rows%5000000==0:print(name,'scanned',rows,flush=True)
    if not pattern.search(line):continue
    cells=line.decode().rstrip('\r\n').split('\t');key=(cells[ci].removeprefix('chr'),cells[pi])
    if key not in targets:continue
    out.write(line);found.append({'chromosome':key[0],'pos38':key[1],'expected_snp':targets[key],'source_ids':cells[si]})
  digest=hashed.h.hexdigest()
 assert digest==prior[path.resolve()]['sha256'],(name,'original source differs from allele-audit source')
 m={'outcome':name,'path':str(path),'source_sha256':digest,'rows_scanned':rows,'full_source_scan':True,'genome_build':'GRCh38','found':found,'extract':dest.name,'extract_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 manifest.write_text(json.dumps(m,indent=2));print(name,'complete',len(found),'records',flush=True);return m
with ThreadPoolExecutor(max_workers=3) as pool:ans=list(pool.map(extract,specs))
(S/'disease_source_manifest.json').write_text(json.dumps(ans,indent=2))
