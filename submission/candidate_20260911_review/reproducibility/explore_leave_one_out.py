"""Post-review SNP influence analysis of the three selected genes.

Uses local harmonized summary statistics; exports aggregate estimates only.
This is post hoc sensitivity work, not a patient subgroup analysis, multivariable
adjustment, independent replication, or proof of absence of pleiotropy.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import stats

p = argparse.ArgumentParser()
p.add_argument('--harmonized', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
D = Path(__file__).resolve().parents[1]
a.out.mkdir(parents=True, exist_ok=False)
plan = {'status': 'Post hoc exploratory sensitivity; not incorporated in manuscript',
        'genes': ['TSHR', 'IGF1R', 'CTLA4'],
        'outcomes': ['BBJ_Graves', 'UKB_hyperthyroid', 'FinnGen_GO'],
        'selection': 'All three manuscript-selected genes and all three outcomes; no selection on new P values',
        'method': 'Drop each retained SNP in turn. Use multiplicative random-effects IVW with no under-dispersion for >=2 retained SNPs; Wald ratio if one remains. Single-SNP starting sets cannot be examined.',
        'frequency_scenario': 'original_reference',
        'interpretation': 'Report all estimates and intervals. P<0.05 crossings are descriptive stability checks, not additional discoveries.',
        'input_sha256': hashlib.sha256(a.harmonized.read_bytes()).hexdigest()}
(a.out/'analysis_plan.json').write_text(json.dumps(plan, indent=2)+'\n', encoding='utf8')
h = pd.read_csv(a.harmonized)
h = h[h.mr_keep.astype(str).str.upper().eq('TRUE') & h.exposure.isin(plan['genes'])].copy()
c = pd.read_csv(D/'provenance/MR_primary_canonical.csv')

def estimate(g):
    x = g['beta.exposure'].to_numpy()
    y = g['beta.outcome'].to_numpy()
    sy = g['se.outcome'].to_numpy()
    assert np.isfinite([x,y,sy]).all() and (sy>0).all() and (x!=0).all()
    if len(g)==1:
        b = y[0]/x[0]; se = sy[0]/abs(x[0]); q = np.nan
        method = 'Wald ratio'
    else:
        w = 1/sy**2
        b = np.dot(w*x,y)/np.dot(w*x,x)
        q = np.dot(w,(y-b*x)**2)
        se = np.sqrt(max(1,q/(len(g)-1))/np.dot(w*x,x))
        # Independent least-squares solution on whitened data.
        qr_b = np.linalg.lstsq((x/sy)[:,None],y/sy,rcond=None)[0][0]
        assert np.isclose(b,qr_b,rtol=1e-12,atol=1e-14)
        method = 'Inverse variance weighted'
    return {'method': method, 'n_iv':len(g), 'beta':float(b), 'se':float(se),
            'pvalue':float(2*stats.norm.sf(abs(b/se))), 'OR':float(np.exp(b)),
            'CI_lower':float(np.exp(b-stats.norm.ppf(.975)*se)),
            'CI_upper':float(np.exp(b+stats.norm.ppf(.975)*se)), 'Q':float(q)}

rows=[]; summaries=[]; not_estimable=[]
for gene in plan['genes']:
    for outcome in plan['outcomes']:
        g=h[(h.exposure==gene)&(h.outcome==outcome)]
        assert not g.SNP.duplicated().any() and len(g)>0
        base=estimate(g)
        original=c[(c.gene_symbol==gene)&(c.outcome==outcome)].iloc[0]
        for k in ['beta','se','pvalue']:
            assert np.isclose(base[k],original[k],rtol=1e-8,atol=1e-30), (gene,outcome,k)
        assert base['n_iv']==original.n_iv
        if len(g)<2:
            not_estimable.append({'gene':gene,'outcome':outcome,'reason':'Only one original instrument'})
            continue
        prefix={'gene':gene,'outcome':outcome,'original_n_iv':len(g)}
        rows.append(prefix|{'excluded_SNP':'None (all instruments)'}|base)
        omitted=[]
        for snp in g.SNP:
            r=estimate(g[g.SNP!=snp]);omitted.append(r)
            rows.append(prefix|{'excluded_SNP':snp}|r)
        summaries.append(prefix|{'all_OR':base['OR'],'all_P':base['pvalue'],
                         'LOO_OR_min':min(r['OR'] for r in omitted),
                         'LOO_OR_max':max(r['OR'] for r in omitted),
                         'LOO_P_min':min(r['pvalue'] for r in omitted),
                         'LOO_P_max':max(r['pvalue'] for r in omitted),
                         'direction_reversals':int(sum(np.sign(r['beta'])!=np.sign(base['beta']) for r in omitted)),
                         'nominal_P_crossings':int(sum((r['pvalue']<.05)!=(base['pvalue']<.05) for r in omitted))})
pd.DataFrame(rows).to_csv(a.out/'leave_one_out_estimates.csv',index=False,na_rep='NA')
pd.DataFrame(summaries).to_csv(a.out/'leave_one_out_summary.csv',index=False)
result={'status':'PASS','baseline_comparisons':9,'eligible_comparisons':len(summaries),
        'omission_estimates':sum(r['excluded_SNP']!='None (all instruments)' for r in rows),
        'not_estimable':not_estimable,'summary':summaries,
        'not_tested':['Age or sex subgroups','Multivariable clinical adjustment','Multi-signal colocalization','Formal directionality']}
(a.out/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
print(json.dumps(result,indent=2))
