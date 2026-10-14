from pathlib import Path
import sys,csv,re,json,math,hashlib
from docx import Document
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Association_Manuscript_20261003';S=R/'outputs/TED_Manuscript_Reframed_20261002'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def rows(rel):
    with (R/rel).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
mr_path='outputs/TED_Reanalysis_20261002/Tables/All_2544_genes_five_outcomes.csv'
col_path='outputs/TED_Followup_Validation_20261002/Tables/Table_F1_24_genes_five_outcomes.csv'
freq_path='outputs/TED_Reanalysis_20261002/Tables/New_MR_all_frequency_scenarios.csv'
mr=rows(mr_path);col=rows(col_path);freq=rows(freq_path)
idx={(r['gene_symbol'],r['outcome']):r for r in mr}
ci={(r['gene'],r['outcome']):r for r in col}
fi={(r['gene_symbol'],r['outcome'],r['scenario']):r for r in freq}
t=(O/'Manuscript.md').read_text(encoding='utf-8');results=t.split('## Results\n',1)[1].split('## Discussion\n',1)[0]
checks=[]
def check(label,value):
    assert value,label
    checks.append(label)
check('12,720 primary rows for 2,544 genes and five datasets',len(mr)==12720)
outcomes=['BBJ_Graves','E4_GRAVES_STRICT','UKB_hyperthyroid','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
for outcome,n in zip(outcomes,[2234,2480,2505,2480,2480]):check('Estimable count '+outcome,sum(r['outcome']==outcome and r['estimable']=='True' for r in mr)==n)
threshold=.05/2544
hits=lambda o:{r['gene_symbol'] for r in mr if r['outcome']==o and r['estimable']=='True' and float(r['pvalue'])<threshold}
bbj=hits('BBJ_Graves');fg=hits('E4_GRAVES_STRICT')
check('13 BBJ discoveries',len(bbj)==13);check('13 FinnGen GD screen hits',len(fg)==13)
check('Four shared full-screen genes',bbj&fg=={'TSHR','CTLA4','HLA-DQA2','TUBB'})
focused={g for g in bbj if idx[g,'E4_GRAVES_STRICT']['estimable']=='True' and float(idx[g,'E4_GRAVES_STRICT']['pvalue'])<.05/13 and float(idx[g,'E4_GRAVES_STRICT']['beta'])*float(idx[g,'BBJ_Graves']['beta'])>0}
check('Seven direction-consistent focused comparisons',focused=={'TSHR','CTLA4','HLA-DQA2','TUBB','C4A','TNFSF14','IFNGR1'})
check('Nine additional FinnGen GD genes',fg-bbj=={'APOM','NEU1','IGLV2-18','IGLV3-12','TNF','HIST1H4E','LTA','LTB','C4B'})

# Validate the values newly given in running text, not just unchanged tables.
super_digits=str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def pstr(v):
    v=float(v)
    if v<.001:
        a,b=f'{v:.2e}'.split('e');return a+'×10'+str(int(b)).translate(super_digits)
    return f'{v:.3g}'
for gene in ['TSHR','IGF1R']:
    for outcome in outcomes:
        row=idx[gene,outcome]
        vals=[f"{float(row['OR']):.2f}",f"{float(row['CI_lower']):.2f}–{float(row['CI_upper']):.2f}",pstr(row['pvalue'])]
        for val in vals:check(f'{gene} {outcome} displayed {val}',val in results)
        check(f'{gene} {outcome} direction',float(row['OR'])<1 if gene=='TSHR' else float(row['OR'])>1)
check('IGF1R nominal only BBJ/UKB',{o for o in outcomes if float(idx['IGF1R',o]['pvalue'])<.05}=={'BBJ_Graves','UKB_hyperthyroid'})
check('IGF1R none full-screen',all(float(idx['IGF1R',o]['pvalue'])>threshold for o in outcomes))
check('TSHR strict GO nominal but not full-screen',threshold<float(idx['TSHR','E4_GRAVES_OPHT_STRICT']['pvalue'])<.05)
for gene in ['CTLA4','HLA-DQA2','TUBB']:
    row=idx[gene,'E4_GRAVES_STRICT']
    for val in [f"{float(row['OR']):.2f}",f"{float(row['CI_lower']):.2f}–{float(row['CI_upper']):.2f}",pstr(row['pvalue'])]:check(f'{gene} GD displayed {val}',val in results)
for gene in ['CTLA4','TUBB']:
    check(gene+' all five nominal',all(float(idx[gene,o]['pvalue'])<.05 for o in outcomes))
for gene in ['CTLA4','TUBB','C4A']:
    val=pstr(fi[gene,'E4_GRAVES_STRICT','original_reference']['Q_pvalue'])
    check(gene+' displayed heterogeneity '+val,val in results)
for gene in ['C4A','TNFSF14','IFNGR1','HLA-A']:
    val=pstr(idx[gene,'E4_GRAVES_STRICT']['pvalue']);check(gene+' displayed GD P '+val,val in results)
for gene in ['TSHR','CTLA4']:
    for o in outcomes:
        value=f"{float(ci[gene,o]['PP.H4']):.3f}";check(gene+' H4 '+o+' '+value,value in results)
for o in ['UKB_hyperthyroid','E4_GRAVES_STRICT']:
    value=f"{float(ci['TSHR',o]['PP.H3']):.3f}";check('TSHR H3 '+o,value in results)
for gene in ['HLA-DQA2','TUBB']:check(gene+' H4 below 0.001 in all five',all(float(ci[gene,o]['PP.H4'])<.001 for o in outcomes))
check('No selected FinnGen GD H4 reaches 0.80',all(float(r['PP.H4'])<.8 for r in col if r['outcome']=='E4_GRAVES_STRICT'))
for o in ['BBJ_Graves','FinnGen_GO','E4_GRAVES_OPHT_STRICT']:
    val=f"{float(ci['TSHR',o]['H4_low_prior']):.3f}";check('TSHR lower prior '+o,val in results)

used=[]
body,refs=t.split('## References\n',1)
for group in re.findall(r'\[([0-9, ]+)\]',body):
    for n in map(int,re.findall(r'\d+',group)):
        if n not in used:used.append(n)
check('All 32 references cited in first-appearance order',used==list(range(1,33)))
check('Reference list sequential',list(map(int,re.findall(r'(?m)^(\d+)\. ',refs)))==used)
check('Secondary designation disclosed as post hoc','reporting emphasis was adopted during post hoc reappraisal' in t)
check('External A-subunit hypothesis distinguished','They do not show that an A-subunit shed by human orbital fibroblasts, acting alone, causes TED independently of GD.' in t)
check('Eight discussion paragraphs',len(t.split('## Discussion\n',1)[1].split('### Limitations',1)[0].strip().split('\n\n'))==8)
check('Two limitation paragraphs',len(t.split('### Limitations\n',1)[1].split('### Conclusions',1)[0].strip().split('\n\n'))==2)
abstract=t.split('## Abstract\n',1)[1].split('**Keywords:**',1)[0]
check('Abstract at most 250 whitespace-delimited words',len(abstract.replace('*','').split())<=250)

# Exact DOCX content equivalence, including all references and symbols.
plain=lambda s:re.sub(r'\s+',' ',re.sub(r'(?m)^#{1,3} ','',s).replace('\\*','*').replace('*','')).strip()
for name in ['Manuscript','Cover_Letter_Draft']:
    doc=Document(O/(name+'.docx'));dt='\n'.join(p.text for p in doc.paragraphs)
    check(name+' DOCX agrees with Markdown',plain(dt)==plain((O/(name+'.md')).read_text(encoding='utf-8')))
unchanged={}
for name in ['Tables.docx','Tables.pdf','Tables.md','Supplementary_Material.docx','Supplementary_Material.pdf','Supplementary_Material.md','Supplementary_Dataset.xlsx','Figure_Legends.md','All_Figures.pdf']:
    check(name+' unchanged bytes',sha(O/name)==sha(S/name));unchanged[name]=sha(O/name)
for folder in ['Figures','Tables']:
    for p in (S/folder).iterdir():
        rel=p.relative_to(S);check(str(rel)+' unchanged',sha(O/rel)==sha(p));unchanged[str(rel).replace('\\','/')]=sha(p)
(O/'Verification/unchanged_artifacts.json').write_text(json.dumps(unchanged,indent=2),encoding='utf-8')
audit={'status':'PASS','scope':'Focused source-to-prose and unchanged-data audit; not a fresh GWAS/MR analysis or causal-validity certification','checks_count':len(checks),'checks':checks,'source_hashes':{p:sha(R/p) for p in [mr_path,col_path,freq_path]},'unchanged_supporting_artifacts':len(unchanged),'prior_workbook_validation':'Original_package_verification/final_numeric_audit.json; workbook bytes unchanged'}
(O/'Verification/final_editorial_numeric_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:audit[k] for k in ['status','checks_count','unchanged_supporting_artifacts']},indent=2))
