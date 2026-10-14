#!/usr/bin/env python3
"""Second implementation of coloc.abf posteriors in 50-digit decimal arithmetic
(direct Bayes-factor sums, no log-space tricks), recomputed from the SNP-level
files written by 03_analysis.py, as an internal numerical check.
Usage: python -I 04_check_decimal.py WORKDIR"""
import sys
from decimal import Decimal as D, getcontext
from pathlib import Path
import pandas as pd
getcontext().prec = 50
work = Path(sys.argv[1])
res = pd.read_csv(work / 'results' / 'coloc_variantA.csv')
worst = 0.0
W1, W2 = D('0.2') ** 2, D('0.15') ** 2
for f in sorted((work / 'results' / 'snp_level').glob('*__*.csv')):
    if f.name.startswith('MR_'):
        continue
    gene, ep = f.stem.split('__')
    d = pd.read_csv(f)
    S1 = S2 = S12 = D(0)
    for b, s, z2, v2 in zip(d.beta_assessed, d.sebeta, d.Zscore, d.V2):
        b, s, z2, v2 = D(repr(b)), D(repr(s)), D(repr(z2)), D(repr(v2))
        z1 = b / s; v1 = s * s
        r1 = W1 / (W1 + v1); r2 = W2 / (W2 + v2)
        a1 = (1 - r1).sqrt() * (r1 * z1 * z1 / 2).exp()
        a2 = (1 - r2).sqrt() * (r2 * z2 * z2 / 2).exp()
        S1 += a1; S2 += a2; S12 += a1 * a2
    for p12s in ['1e-5', '5e-6', '1e-6']:
        p12 = D(p12s); p1 = p2 = D('1e-4')
        h = [D(1), p1 * S1, p2 * S2, p1 * p2 * (S1 * S2 - S12), p12 * S12]
        tot = sum(h)
        pp = [float(x / tot) for x in h]
        row = res[(res.gene == gene) & (res.endpoint == ep) & ((res.p12 - float(p12)).abs() < 1e-12)].iloc[0]
        dmax = max(abs(pp[i] - row[f'PP.H{i}']) for i in range(5))
        worst = max(worst, dmax)
        print(gene, ep, p12s, ' '.join(f'{x:.6g}' for x in pp), f'maxdiff={dmax:.2e}')
print('WORST', worst)
