import csv,json
from pathlib import Path
p=Path(__file__).resolve().parent/'maf_sensitivity_20260905/identity_control_v2'
base={(r['gene_symbol'],r['outcome']):r for r in csv.DictReader((p/'baseline_v1/primary_MR_reproduced.csv').open())}
new=list(csv.DictReader((p/'sensitivity_v1/primary_MR_all_scenarios.csv').open()))
summary={}
for scenario in sorted({r['scenario'] for r in new}):
    rows=[r for r in new if r['scenario']==scenario and int(r['n_iv'])>0]
    assert len(rows)==len(base)==7219
    delta=0.
    for r in rows:
        b=base[(r['gene_symbol'],r['outcome'])]
        assert r['n_iv']==b['n_iv'] and r['method']==b['method']
        delta=max(delta,*[abs(float(r[f])-float(b[f])) for f in ('beta','se','pvalue')])
    assert delta<1e-8
    summary[scenario]={'rows':len(rows),'max_abs_difference':delta}
(p/'identity_verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
