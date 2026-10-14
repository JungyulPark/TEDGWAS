from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from xml.etree import ElementTree as E
import copy,json,math,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1]
S=R/'outputs/TED_Research_Contribution_20261003/Supplementary_Dataset.xlsx'
O=R/'outputs/TED_Allele_Integrated_20261003';V=O/'Verification'
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main';P='http://schemas.openxmlformats.org/package/2006/relationships';RR='http://schemas.openxmlformats.org/officeDocument/2006/relationships';CT='http://schemas.openxmlformats.org/package/2006/content-types'
q=lambda x:'{'+N+'}'+x
def xml(e):
    E.register_namespace('',e.tag.split('}')[0][1:]);E.register_namespace('r',RR)
    return E.tostring(e,encoding='utf-8',xml_declaration=True)
with ZipFile(S) as z:old={n:z.read(n) for n in z.namelist()}
with ZipFile(W/'edited_subset.xlsx') as z:new={n:z.read(n) for n in z.namelist()}
out=old.copy();os=E.fromstring(old['xl/styles.xml']);ns=E.fromstring(new['xl/styles.xml'])
maps={}
for name in ['numFmts','fonts','fills','borders','cellStyleXfs','cellXfs']:
    target=os.find(q(name));source=ns.find(q(name));mapping={}
    if target is None:
        target=E.Element(q(name));os.insert(0,target)
    if source is not None:
        for i,elem in enumerate(source):
            ne=copy.deepcopy(elem)
            if name=='numFmts':
                nextid=max([163]+[int(x.attrib['numFmtId']) for x in target])+1
                mapping[int(ne.attrib['numFmtId'])]=nextid;ne.set('numFmtId',str(nextid))
            else:mapping[i]=len(target)
            for attr,kind in [('fontId','fonts'),('fillId','fills'),('borderId','borders'),('numFmtId','numFmts'),('xfId','cellStyleXfs')]:
                if attr in ne.attrib and kind in maps:
                    val=int(ne.attrib[attr]);ne.set(attr,str(maps[kind].get(val,val)))
            target.append(ne)
    target.set('count',str(len(target)));maps[name]=mapping
out['xl/styles.xml']=xml(os)
shared=[]
if 'xl/sharedStrings.xml' in new:
    shared=list(E.fromstring(new['xl/sharedStrings.xml']))
def normalize_sheet(b):
    e=E.fromstring(b)
    for cell in e.iter(q('c')):
        if 's' in cell.attrib:cell.set('s',str(maps['cellXfs'][int(cell.attrib['s'])]))
        if cell.attrib.get('t')=='s':
            si=copy.deepcopy(shared[int(cell.find(q('v')).text)]);si.tag=q('is')
            for ch in list(cell):cell.remove(ch)
            cell.set('t','inlineStr');cell.append(si)
    for elem in e.iter():
        for att in ['style','s']:
            if elem.tag!=q('c') and att in elem.attrib:elem.set(att,str(maps['cellXfs'][int(elem.attrib[att])]))
    return e
edited={i:normalize_sheet(new[f'xl/worksheets/sheet{i}.xml']) for i in [1,2,3]}
for i,last in [(1,24),(2,278)]:
    e=E.fromstring(old[f'xl/worksheets/sheet{i}.xml']);data=e.find(q('sheetData'));nd=edited[i].find(q('sheetData'))
    for row in nd:
        if int(row.attrib['r'])>last:data.append(copy.deepcopy(row))
    if i==1:
        r=data.find(q('row'));cell=next(c for c in r if c.attrib['r']=='B1');idx=list(r).index(cell)
        replacement=next(c for c in nd.find(q('row')) if c.attrib['r']=='B1')
        r.remove(cell);r.insert(idx,copy.deepcopy(replacement))
    e.find(q('dimension')).set('ref',f'A1:{"B" if i==1 else "C"}{max(int(r.attrib["r"]) for r in data)}')
    out[f'xl/worksheets/sheet{i}.xml']=xml(e)
out['xl/worksheets/sheet20.xml']=xml(edited[3])
wb=E.fromstring(old['xl/workbook.xml']);rels=E.fromstring(old['xl/_rels/workbook.xml.rels']);ct=E.fromstring(old['[Content_Types].xml'])
rid='rId'+str(max(int(x.attrib['Id'][3:]) for x in rels)+1)
E.SubElement(wb.find(q('sheets')),q('sheet'),{'name':'TSHR_allele_audit','sheetId':'20','state':'visible','{'+RR+'}id':rid})
E.SubElement(rels,'{'+P+'}Relationship',{'Id':rid,'Type':RR+'/worksheet','Target':'/xl/worksheets/sheet20.xml'})
E.SubElement(ct,'{'+CT+'}Override',{'PartName':'/xl/worksheets/sheet20.xml','ContentType':'application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml'})
out['xl/workbook.xml']=xml(wb);out['xl/_rels/workbook.xml.rels']=xml(rels);out['[Content_Types].xml']=xml(ct)
with ZipFile(O/'Supplementary_Dataset.xlsx','w',ZIP_DEFLATED) as z:
    for n,b in out.items():z.writestr(n,b)

checks=[]
def check(name,ok):
    assert ok,name
    checks.append(name)
for i in range(3,20):check(f'Existing data worksheet {i} byte preserved',old[f'xl/worksheets/sheet{i}.xml']==out[f'xl/worksheets/sheet{i}.xml'])
def cells(e):
    d={}
    for c in e.iter(q('c')):
        if c.attrib.get('t')=='inlineStr':v=''.join(c.find(q('is')).itertext())
        elif c.attrib.get('t')=='str':v=c.find(q('v')).text if c.find(q('v')) is not None else ''
        elif c.attrib.get('t')=='b':v=c.find(q('v')).text=='1'
        elif c.find(q('v')) is not None:v=float(c.find(q('v')).text)
        else:v=None
        d[c.attrib['r']]=v
    return d
for i in [1,2]:
    a=cells(E.fromstring(old[f'xl/worksheets/sheet{i}.xml']));b=cells(E.fromstring(out[f'xl/worksheets/sheet{i}.xml']))
    for addr,val in a.items():
        if i==1 and addr=='B1':continue
        check(f'Original metadata {i}:{addr}',b[addr]==val)
cs=cells(edited[3]);blocks=json.loads((W/'blocks.json').read_text())
num=0
for block in blocks:
    for ri,row in enumerate(block['values']):
        for ci,val in enumerate(row):
            addr=f'{chr(65+ci)}{block["row"]+2+ri}';actual=cs.get(addr)
            if block['headers'][ci]=='method':val='PLINK ML, unphased'
            if block['headers'][ci]=='phase_basis':val='Inferred from unphased genotypes'
            if isinstance(val,(int,float)) and not isinstance(val,bool):
                check('Audit number '+addr,actual is not None and math.isclose(actual,val,rel_tol=1e-13,abs_tol=0));num+=1
            else:check('Audit text '+addr,actual==val)
check('Final workbook sheet count',len(wb.find(q('sheets')))==20)
report={'status':'PASS','checks':len(checks),'new_numeric_cells_verified':num,'unchanged_data_worksheets':17,'changed_metadata':['Readme B1 and row25','Column_dictionary new rows279:353'],'new_sheet':'TSHR_allele_audit','preservation_method':'Artifact Tool authoring of changed content; OPC merge preserves unrelated original worksheet bytes and existing style definitions','sha256':hashlib.sha256((O/'Supplementary_Dataset.xlsx').read_bytes()).hexdigest()}
(V/'workbook_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
