from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from xml.etree import ElementTree as E
import pandas as pd,json
W=Path(__file__).resolve().parent;R=W.parents[1];S=R/'outputs/TED_Allele_Integrated_20261003/Supplementary_Dataset.xlsx';A=R/'outputs/TED_Tissue_Extension_20261003'
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main';P='http://schemas.openxmlformats.org/package/2006/relationships';CT='http://schemas.openxmlformats.org/package/2006/content-types'
with ZipFile(S) as z:parts={n:z.read(n) for n in z.namelist()}
for file,ns in [('xl/workbook.xml',N),('xl/_rels/workbook.xml.rels',P),('[Content_Types].xml',CT)]:
 root=E.fromstring(parts[file]);E.register_namespace('',ns)
 if file=='xl/workbook.xml':
  items=root.find('{'+N+'}sheets')
  for x in list(items)[2:]:items.remove(x)
 elif file.endswith('.rels'):
  for x in list(root):
   if x.attrib['Type'].endswith('/worksheet') and x.attrib['Id'] not in ['rId1','rId2']:root.remove(x)
 else:
  for x in list(root):
   name=x.attrib.get('PartName','')
   if name.startswith('/xl/worksheets/') and name not in ['/xl/worksheets/sheet1.xml','/xl/worksheets/sheet2.xml']:root.remove(x)
 parts[file]=E.tostring(root)
with ZipFile(W/'edit_subset.xlsx','w',ZIP_DEFLATED) as z:
 for n,b in parts.items():
  if n.startswith('xl/worksheets/') and n not in ['xl/worksheets/sheet1.xml','xl/worksheets/sheet2.xml']:continue
  z.writestr(n,b)
blocks=[];row=5
for title,file in [('A Fixed tissue availability','Tissue_availability_10_combinations.csv'),('B All 50 planned MR comparisons','Tissue_MR_50_planned_combinations.csv'),('C Four selected instruments','Selected_instruments_verified.csv'),('D All aligned disease records','Harmonized_instruments_all_outcomes.csv'),('E All candidates before clumping','Candidate_instruments_mapped.csv')]:
 df=pd.read_csv(A/'Tables'/file)
 if title.startswith('B '):
  first=['gene','tissue','outcome','family','status','n_retained_instruments','OR','CI_lower','CI_upper','P','P_adjusted','statistical_category','direction_vs_blood','reason']
  df=df[first+[x for x in df.columns if x not in first]]
 blocks.append({'title':title,'file':file,'row':row,'headers':list(df.columns),'values':df.astype(object).where(pd.notna(df),None).values.tolist()});row+=len(df)+4
(W/'blocks.json').write_text(json.dumps(blocks,allow_nan=False),encoding='utf-8')
print([(x['title'],x['row'],len(x['values']),len(x['headers'])) for x in blocks])
