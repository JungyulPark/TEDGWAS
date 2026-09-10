from pathlib import Path
import re,json,hashlib,sys
import numpy as np,pandas as pd
from docx import Document
O=Path(__file__).resolve().parents[1];P=O/'provenance';master=O.parents[1]/'MANUSCRIPT_TED_TRAP_v5_MASTER.md';text=master.read_text(encoding='utf-8')
mr=pd.read_csv(P/'MR_primary_canonical.csv');raw=pd.read_csv(P/'MR_all_estimators_verified.csv');co=pd.read_csv(P/'coloc_canonical_v2.csv');inst=pd.read_csv(P/'instruments_verified.csv')
allmr=pd.read_csv(O/'Supplementary_Data_2_MR.csv');allco=pd.read_csv(O/'Supplementary_Data_3_Colocalization.csv');power=pd.read_csv(O/'Supplementary_Data_4_Power.csv')
display=json.loads((P/'clinical_display_sources.json').read_text(encoding='utf-8'));tables=display['tables'];checks=[];errors=[];numeric=0
def check(ok,label):
    checks.append(label)
    if not ok:errors.append(label)
def clean(s):return re.sub(r'\s+',' ',str(s).replace('*','').replace('\u00a0',' ')).strip()
tr=str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹−','-0123456789-')
def agree(shown,value,label):
    global numeric
    numeric+=1;s=clean(shown).replace(',','').translate(tr)
    if s=='NA':check(pd.isna(value),label+' unavailable');return
    if '×10' in s:
        mant,exp=s.split('×10');actual=float(mant)*10**int(exp);unit=10**int(exp)*10**(-len(mant.split('.')[1]) if '.' in mant else 0)
    else:
        actual=float(s);unit=10**(-len(s.split('.')[1]) if '.' in s else 0)
    check(np.isfinite(value) and abs(actual-value)<=unit/2+1e-12,label)
def effect(s,r,label):
    v=re.findall(r'\d+\.\d+',s);check(len(v)==3,label+' three OR/CI values')
    for a,b in zip(v,np.exp([r.beta,r.beta-1.96*r.se,r.beta+1.96*r.se])):agree(a,b,label)
ocs={'BBJ Graves disease':'BBJ_Graves','UKB hyperthyroidism':'UKB_hyperthyroid','FinnGen Graves ophthalmopathy':'FinnGen_GO'}
def get(g,o):return mr[(mr.gene_symbol==g)&(mr.outcome==o)].iloc[0]
for name,keys in [('MANUSCRIPT_Submission',['1','2','3']),('SUPPLEMENTARY_MATERIAL',['S1','S2','S3','S4'])]:
    doc=Document(O/(name+'.docx'));check(len(doc.tables)==len(keys),name+' table count')
    for tb,key in zip(doc.tables,keys):
        expected=[tables[key]['header']]+tables[key]['rows'];check(len(tb.rows)==len(expected),key+' row count')
        for row,ex in zip(tb.rows,expected):
            check(len(row.cells)==len(ex),key+' column count')
            for c,e in zip(row.cells,ex):check(clean(c.text)==clean(e),key+' DOCX cell '+clean(e))
            check(bool(row._tr.xpath('./w:trPr/w:cantSplit')),key+' row does not split')
        check(bool(tb.rows[0]._tr.xpath('./w:trPr/w:tblHeader')),key+' repeating header')
        for ex in expected:check('| '+' | '.join(map(str,ex))+' |' in text,key+' source row exists')
    alltext='\n'.join([p.text for p in doc.paragraphs]+[c.text for tb in doc.tables for row in tb.rows for c in row.cells])
    check('�' not in alltext,name+' no replacement characters')
    check(not re.search(r'\bP\s*=\s*0\.000(?:\s|[;,)]|$)',alltext),name+' no zero P value display')
for row in tables['2']['rows']:
    g=clean(row[0]);o=ocs[row[1]];r=get(g,o);label=f'Table2 {g}/{o}'
    agree(row[2],r.n_iv,label+' instruments');effect(row[3],r,label);agree(row[4],r.pvalue,label+' P')
    c=co[(co.gene==g)&(co.outcome==o)&np.isclose(co.p12,1e-5,atol=1e-12)].iloc[0];agree(row[5],c['PP.H4'],label+' H4')
for row in tables['3']['rows']:
    r=get(clean(row[0]),'BBJ_Graves');agree(row[1],r.n_iv,'Table3 instruments');effect(row[2],r,'Table3 effect');agree(row[3],r.pvalue,'Table3 P');check(r.pvalue<.05/2544,'Table3 discovery threshold')
check(len(tables['3']['rows'])==13,'All thirteen discovery hits retained')
for row in tables['S1']['rows']:
    g=clean(row[0]);o=ocs[row[1]];r=get(g,o);rr=raw[(raw.gene_symbol==g)&(raw.outcome==o)]
    agree(row[2],r.n_iv,'S1 instruments');agree(row[3],r.pvalue,'S1 primary P')
    for col,method in [(4,'Weighted median'),(5,'Weighted mode')]:
        v=rr[rr.method==method];agree(row[col],v.iloc[0].pvalue if len(v) else np.nan,'S1 '+method)
    agree(row[6],r.egger_intercept_p,'S1 Egger intercept');agree(row[7],r.cochran_q_p,'S1 Q')
for row in tables['S2']['rows']:
    g=clean(row[0]);o=ocs[row[1]]
    for col,p12 in [(2,1e-5),(3,5e-6),(4,1e-6)]:
        c=co[(co.gene==g)&(co.outcome==o)&np.isclose(co.p12,p12,atol=1e-12)].iloc[0];agree(row[col],c['PP.H4'],'S2 H4')
for row in tables['S3']['rows']:
    g=clean(row[0]);o=ocs[row[1]];r=get(g,o);new=allmr[(allmr.gene_symbol==g)&(allmr.outcome==o)&(allmr.scenario=='reharmonized_eqtlgen')].iloc[0]
    effect(row[2],r,'S3 reference');agree(row[3],r.pvalue,'S3 reference P');effect(row[4],new,'S3 cohort');agree(row[5],new.pvalue,'S3 cohort P')
for row in tables['S4']['rows']:
    o=ocs[row[0]];scenario='original_reference' if row[1]=='Reference' else 'reharmonized_eqtlgen';r=power[(power.outcome==o)&(power.scenario==scenario)].iloc[0]
    agree(row[2],r.n_genes,'S4 genes')
    for a,b in zip(re.findall(r'\d+\.\d+',row[3]),[r.or_median,r.or_q1,r.or_q3]):agree(a,b,'S4 detectable OR')
    agree(row[4],100*r.frac_OR1_5,'S4 OR1.5 percentage');agree(row[5],100*r.frac_OR2,'S4 OR2 percentage')
# Confirm the consolidated data retain the old primary results and full scenarios.
ref=allmr[allmr.scenario=='original_reference'].merge(mr,on=['gene_symbol','outcome'],suffixes=('_new','_old'),validate='one_to_one');check(len(ref)==7219,'All original MR estimates in consolidated file')
for col in ['beta','se','pvalue','n_iv']:check(np.allclose(ref[col+'_new'],ref[col+'_old'],atol=1e-12,rtol=1e-12),col+' primary data unchanged')
check(len(allmr)==30115,'Complete merged MR row count');check(set(allmr.scenario)=={'original_reference','paired_reference','paired_eqtlgen','reharmonized_eqtlgen'},'All MR scenarios retained')
check(len(allco)==324 and set(allco.scenario)=={'original_reference','paired_reference','paired_eqtlgen','available_eqtlgen'},'All coloc rows/scenarios retained')
check(len(power)==12,'All power scenarios retained')
check(inst.gene_symbol.nunique()==2544 and len(inst)==6135,'Instrument manifest unchanged')
refco=allco[allco.scenario=='original_reference'].merge(co,on=['gene','outcome','p12'],suffixes=('_new','_old'),validate='one_to_one');check(len(refco)==81,'All original coloc settings retained')
for k in range(5):check(np.allclose(refco[f'PP.H{k}_new'],refco[f'PP.H{k}_old'],atol=1e-12,rtol=1e-12),'Original H'+str(k)+' unchanged')
v=json.loads((P/'maf/sensitivity_verification_v2.json').read_text(encoding='utf-8'));check(v['status']=='PASS' and v['max_abs_delta_Z']<1e-8,'Completed sensitivity verification')
check(all(r['n_threshold_crossings']==0 for r in v['coloc_summary']['comparisons']),'Zero frequency threshold crossings retained')
for r in display['narrative_estimates']:
    m=get(r['gene'],r['outcome']);agree(r['P'],m.pvalue,'Narrative P');agree(r['OR'],np.exp(m.beta),'Narrative OR')
    check(r['P'] in text and r['OR'] in text and r['CI'] in text,'Narrative rendered values present')
fs=json.loads((P/'clinical_figure_sources.json').read_text(encoding='utf-8'))
for r in fs['Figure2']:
    m=get(r['gene'],r['outcome']);check(np.isclose(r['pvalue'],m.pvalue,rtol=1e-12,atol=0),'Figure2 source P');check(np.isclose(r['or'],np.exp(m.beta)),'Figure2 source OR');agree(r['displayed_p'],m.pvalue,'Figure2 display P')
check(not re.search(r'\b(?:Table|Tables) S(?:[5-9]|1[0-2])\b',text),'No old supplementary table references')
check('Figure 3.' not in text and 'Figure S2.' not in text,'No removed figure references')
check('{{' not in text and 'TODO' not in text,'No unresolved generation placeholders')
check('Vattikuti S, Vattikuti S' not in text,'Duplicated reference author corrected')
check('All authors read and approved the final manuscript.' not in text,'No claim that this new revision has already been approved')
check('not a P value' in text and 'not treatment effects' in text,'Probability and treatment-effect interpretation retained')
check('P* = 0.0056' in text,'TNFSF14 UKB nominal result retained')
counts=json.loads((P/'word_counts.json').read_text(encoding='utf-8'));check(counts['abstract_words']<=250,'Abstract within working 250-word ceiling');check(counts['main_words_including_headings']<=5000,'Main within working 5000-word ceiling')
checklist=Document(O/'STROBE_MR_CHECKLIST.docx');check(len(checklist.tables)==1,'Separate reporting checklist exists')
result={'status':'PASS' if not errors else 'FAIL','numeric_cells_compared':numeric,'checks':len(checks),'errors':errors,'scope':'Rounded values versus original analytical outputs, complete consolidated data, DOCX table transcription and essential limitations. Visual layout, author declarations and live journal requirements require separate review.'}
(P/'integrity_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2));sys.exit(bool(errors))
