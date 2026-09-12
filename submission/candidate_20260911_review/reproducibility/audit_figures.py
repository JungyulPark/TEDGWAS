"""Check every plotted gene and posterior against preserved analytical sources."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
O=Path(__file__).resolve().parents[1];P=O/'provenance'
j=json.loads((P/'clinical_figure_sources.json').read_text(encoding='utf-8'))
mr=pd.read_csv(P/'MR_primary_canonical.csv');ins=pd.read_csv(P/'instruments_verified.csv');co=pd.read_csv(P/'coloc_canonical_v2.csv')
errors=[];checks=0
def check(ok,label):
    global checks
    checks+=1
    if not ok:errors.append(label)
for n,h in j['input_sha256'].items():check(hashlib.sha256((P/n).read_bytes()).hexdigest()==h,'Input hash '+n)
screen=pd.DataFrame(j['Figure1_screen']);bbj=mr[mr.outcome=='BBJ_Graves'].set_index('gene_symbol')
check(len(screen)==len(bbj)==2234 and set(screen.gene_symbol)==set(bbj.index),'Complete 2234-gene screen')
check(not screen.gene_symbol.duplicated().any(),'No duplicated plotted genes')
anchors=ins.loc[ins.groupby('gene_symbol').zscore.apply(lambda s:s.abs().idxmax())].set_index('gene_symbol')
for r in screen.itertuples():
    check(np.isclose(r.pvalue,bbj.loc[r.gene_symbol,'pvalue'],rtol=1e-12,atol=0),'Screen P '+r.gene_symbol)
    check(np.isclose(r.y,-np.log10(r.pvalue)),'Screen ordinate '+r.gene_symbol)
    a=anchors.loc[r.gene_symbol];check(r.chr==a.chr and r.pos_hg19==a.pos_hg19,'Strongest instrument position '+r.gene_symbol)
check(sum(screen.pvalue<.05/2544)==13,'Thirteen plotted discovery hits')
for key,count in [('Figure2_posteriors',9),('Figure3_primary',27),('Figure3_priors',27)]:
    check(len(j[key])==count,key+' complete comparisons')
    for r in j[key]:
        c=co[(co.gene==r['gene'])&(co.outcome==r['outcome'])&np.isclose(co.p12,r.get('p12',1e-5),rtol=1e-9,atol=0)].iloc[0]
        for name,v in r.items():
            if name.startswith('PP.H'):check(np.isclose(v,c[name],rtol=1e-12,atol=1e-15),key+' '+r['gene']+'/'+r['outcome']+'/'+name)
        if key=='Figure2_posteriors':check(abs(sum(r[h] for h in ['PP.H0','PP.H1','PP.H2','PP.H3','PP.H4'])-1)<1e-10,'Posterior bars sum to one')
result={'status':'FAIL' if errors else 'PASS','checks':checks,'screen_genes':len(screen),'posterior_comparisons':63,'errors':errors}
(P/'figure_numeric_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2));raise SystemExit(bool(errors))
