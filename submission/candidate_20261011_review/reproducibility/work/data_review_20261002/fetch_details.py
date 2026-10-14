from pathlib import Path
import urllib.request,json,re,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2];M=ROOT/'outputs/TED_Data_Review_20261002/provenance'
urls={
'UKB_harmonized_yaml':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038636/harmonised/33959723-GCST90038636-EFO_0009189.h.tsv.gz-meta.yaml',
'BBJ_harmonized_yaml':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018627/harmonised/GCST90018627.h.tsv.gz-meta.yaml',
'UKB_paper_fulltext':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7610725/fullTextXML',
'FinnGen_R12_client':'https://r12.finngen.fi/static/js/main.5d056c54.js'}
for name,url in urls.items():
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=25) as r:b=r.read()
  (M/(name+'.txt')).write_bytes(b);print(name,len(b),flush=True)
  if name.endswith('fulltext'):
   root=ET.fromstring(b)
   text='\n'.join(''.join(p.itertext()) for p in root.findall('.//p'))
   (M/'UKB_paper_paragraphs.txt').write_text(text,encoding='utf8')
  if name.endswith('client'):
   s=b.decode();matches=[s[max(0,m.start()-60):m.end()+140] for m in re.finditer(r'/api/|/phenoinfo|/phenotypes|/pheno/',s)]
   (M/'FinnGen_api_routes.txt').write_text('\n'.join(matches),encoding='utf8')
 except Exception as e:print(name,str(e)[:200],flush=True)
