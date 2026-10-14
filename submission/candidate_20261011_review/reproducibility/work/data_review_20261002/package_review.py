from pathlib import Path
import json,hashlib,zipfile,shutil,csv,math,re
ROOT=Path(__file__).resolve().parents[2];B=Path(__file__).parent;O=ROOT/'outputs/TED_Data_Review_20261002'
data=json.loads((O/'review_data.json').read_text(encoding='utf8'));checks=0
def ck(ok,label):
 global checks
 assert ok,label
 checks+=1
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=list(csv.DictReader((O/'All_2544_genes_3_outcomes.csv').open(encoding='utf-8-sig')))
ck(len(source)==len(data['genes'])==7632,'all gene-outcome rows')
for a,b in zip(source,data['genes']):
 ck(a['gene_symbol']==b['gene_symbol'] and a['outcome']==b['outcome'],'row identity')
 for field in ['n_iv','beta','se','pvalue','or','or_ci_lower','or_ci_upper']:
  if not a[field]:ck(b[field] is None,'null stays null')
  else:ck(math.isclose(float(a[field]),b[field],rel_tol=1e-14,abs_tol=0),'full precision '+field)
html=(O/'START_HERE.html').read_text(encoding='utf8')
embedded=json.loads(re.search(r'<script id="data" type="application/json">(.*?)</script>',html,re.S).group(1))
ck(embedded==data,'HTML JSON matches payload')
q=json.loads((B/'qa/verification.json').read_text(encoding='utf8'));ck(q['status']=='PASS','browser checks')
checks+=len(q['checks'])
R=O/'reproducibility';R.mkdir(exist_ok=True)
for name in ['inspect_sources.py','fetch_metadata.py','fetch_details.py','verify_identity.py','profile_data.py','build_review.py','review_template.html','verify_dashboard.cjs','package_review.py']:shutil.copy2(B/name,R/name)
(R/'README.txt').write_text('These are copies of the scripts executed from work/data_review_20261002. They preserve the original workspace-relative paths and absolute raw input paths for traceability; they are not a portable pipeline. Large raw GWAS/eQTL files are not included. No new SNP-level MR or colocalization was run in this review. See analysis_log.txt and provenance/.\n',encoding='utf8')
report={'status':'PASS','numeric_and_browser_checks':checks,'full_gwas_rows':sum(x['rows'] for x in data['analysis']['raw_qc']),'gene_outcome_rows':len(source),'source_files_modified':False,'manuscript_modified':False,'new_snp_level_MR':False,'new_colocalization':False,'visual_review':'Cohort, QC and gene exploration screens reviewed; no overlap or clipping seen. Desktop/mobile overflow checks passed.','browser_checks':q['checks'],'remaining':'Detailed UKB source phenotype, ancestry and covariates; additional outcome full QC and SNP-level reanalysis; duplicate/INFO/LD checks beyond this descriptive QC.'}
(O/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
manifest={p.relative_to(O).as_posix():sha(p) for p in O.rglob('*') if p.is_file() and p.name!='PACKAGE_SHA256.json'}
(O/'PACKAGE_SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8')
zpath=O.with_suffix('.zip')
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for p in O.rglob('*'):
  if p.is_file():z.write(p,p.relative_to(O).as_posix())
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
SHORT=Path('${LOCAL_USER_HOME}/Documents/Codex/TED_Data_Review_20261002')
shutil.copytree(O,SHORT,dirs_exist_ok=True)
shutil.copy2(zpath,SHORT.with_suffix('.zip'))
for p in O.rglob('*'):
 if p.is_file():assert sha(p)==sha(SHORT/p.relative_to(O))
print(json.dumps(report,ensure_ascii=False));print('Delivered: '+str(SHORT))
