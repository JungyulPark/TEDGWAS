#!/usr/bin/env python3
"""Second, arbitrary-precision implementation of coloc.abf posteriors (no logsumexp),
recomputed from the SNP-level files written by 03_analysis.py, as an internal check.
Usage: python -I 04_check_mpmath.py WORKDIR"""
import sys
from pathlib import Path
import mpmath as mp
import pandas as pd
mp.mp.dps = 60
work = Path(sys.argv[1])
res = pd.read_csv(work / 'results' / 'coloc_variantA.csv')
worst = 0.0
for f in sorted((work / 'results' / 'snp_level').glob('*__*.csv')):
    if f.name.startswith('MR_'):
        continue
    gene, ep = f.stem.split('__')
    d = pd.read_csv(f)
    S1 = S2 = S12 = mp.mpf(0)
    for b, s, z2, v2 in zip(d.beta_assessed, d.sebeta, d.Zscore, d.V2):
        z1 = mp.mpf(b) / mp.mpf(s); v1 = mp.mpf(s) ** 2
        r1 = mp.mpf('0.04') / (mp.mpf('0.04') + v1)
        r2 = mp.mpf('0.0225') / (mp.mpf('0.0225') + mp.mpf(v2))
        a1 = mp.sqrt(1 - r1) * mp.exp(r1 * z1 ** 2 / 2)
        a2 = mp.sqrt(1 - r2) * mp.exp(r2 * mp.mpf(z2) ** 2 / 2)
        S1 += a1; S2 += a2; S12 += a1 * a2
    for p12 in [mp.mpf('1e-5'), mp.mpf('5e-6'), mp.mpf('1e-6')]:
        p1 = p2 = mp.mpf('1e-4')
        h = [mp.mpf(1), p1 * S1, p2 * S2, p1 * p2 * (S1 * S2 - S12), p12 * S12]
        tot = sum(h)
        pp = [float(x / tot) for x in h]
        row = res[(res.gene == gene) & (res.endpoint == ep) & ((res.p12 - float(p12)).abs() < 1e-12)].iloc[0]
        dmax = max(abs(pp[i] - row[f'PP.H{i}']) for i in range(5))
        worst = max(worst, dmax)
        print(gene, ep, float(p12), ' '.join(f'{x:.6g}' for x in pp), f'maxdiff={dmax:.2e}')
print('WORST', worst)
