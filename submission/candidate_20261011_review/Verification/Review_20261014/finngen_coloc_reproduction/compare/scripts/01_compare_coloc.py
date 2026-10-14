"""Compare reproduction A and B coloc outputs with each other, September identity control and Codex.
Written from scratch for the comparison; reads only CSVs. Writes only under compare_AB/results."""
import gzip, json
from pathlib import Path
import numpy as np, pandas as pd

S = Path('<scratch>')
A = S / 'coloc_repro_A/results'; B = S / 'coloc_repro_B/results'; OUT = S / 'compare_AB/results'
R = Path('/home/user/TEDGWAS')
CODEX = R / 'submission/candidate_20261011_review/reproducibility/work/reanalysis_20261002/coloc_v1/Coloc_all_settings_with_identity.csv'
SEPT = R / 'archives/submission_candidates/candidate_20260911_review/provenance/coloc_canonical_v2.csv'
GENES = ['TSHR', 'IGF1R', 'CTLA4']; EPS = ['FinnGen_GO', 'E4_GRAVES_STRICT', 'E4_GRAVES_OPHT_STRICT']
H = [f'PP.H{i}' for i in range(5)]
LP1 = np.log(1e-4)

def key(df):
    df = df.copy(); df['p12k'] = df.p12.map(lambda x: f'{float(x):.0e}'); return df

def lse(x):
    m = np.max(x); return m + np.log(np.sum(np.exp(x - m)))

def labf(z, v, w):
    r = w * w / (w * w + v); return 0.5 * (np.log1p(-r) + r * z * z)

def posteriors(l1, l2, p12):
    s1, s2, s12 = lse(l1), lse(l2), lse(l1 + l2)
    d = s12 - (s1 + s2)
    h3 = s1 + s2 + np.log(-np.expm1(d))
    lh = np.array([0.0, LP1 + s1, LP1 + s2, 2 * LP1 + h3, np.log(p12) + s12])
    return np.exp(lh - lse(lh))

# ---------- load ----------
codex = key(pd.read_csv(CODEX)); codex = codex[(codex.scenario == 'original_reference') & codex.gene.isin(GENES)].rename(columns={'outcome': 'endpoint', 'n_overlap': 'n'})
sept = key(pd.read_csv(SEPT)); sept = sept[(sept.outcome == 'FinnGen_GO') & sept.gene.isin(GENES)].rename(columns={'outcome': 'endpoint', 'n_overlap': 'n'})
a = key(pd.read_csv(A / 'coloc_variantA.csv')).rename(columns={'n_snps': 'n'}); a['variant'] = 'A_composite'
asens = key(pd.read_csv(A / 'coloc_frequency_sensitivity_variantA.csv')).rename(columns={'n_snps': 'n'}); asens['variant'] = 'A_' + asens.sensitivity
b = key(pd.read_csv(B / 'coloc_variantB_all.csv')).rename(columns={'n_snps': 'n'}); b['variant'] = b.freq_scenario
cols = ['variant', 'gene', 'endpoint', 'p12k', 'n', 'top_snp'] + H
allrep = pd.concat([a[cols], asens[cols], b[cols]], ignore_index=True)
K = ['gene', 'endpoint', 'p12k']
out = {}

# ---------- 1. identity: Codex GO vs September ----------
m = codex[codex.endpoint == 'FinnGen_GO'].merge(sept, on=K, suffixes=('_c', '_s'))
out['codex_vs_sept_GO'] = dict(n_rows=len(m), max_abs_dPP=float(max((m[h + '_c'] - m[h + '_s']).abs().max() for h in H)),
                               n_equal=bool((m.n_c == m.n_s).all()), top_equal=bool((m.top_snp_c == m.top_snp_s).all()))

# ---------- 2. A vs B under identical frequency definitions ----------
pairs = [('A_composite', 'B_composite_404nfe_99fin'), ('A_AF_nfe_only', 'B_nfe_primary'), ('A_gnomAD_PASS_only', 'B_composite_404nfe_99fin')]
ab = {}
for x, y in pairs:
    mx = allrep[allrep.variant == x].merge(allrep[allrep.variant == y], on=K, suffixes=('_x', '_y'))
    ab[f'{x}__vs__{y}'] = dict(n_rows=len(mx), max_abs_dPP=float(max((mx[h + '_x'] - mx[h + '_y']).abs().max() for h in H)),
                               n_diff=mx.loc[mx.n_x != mx.n_y, K + ['n_x', 'n_y']].to_dict('records'),
                               top_equal=bool((mx.top_snp_x == mx.top_snp_y).all()))
out['A_vs_B_same_frequency'] = ab
print('n per A sensitivity:', asens.groupby(['variant', 'gene']).n.first().to_dict())

# ---------- 3. deltas vs Codex (all endpoints) & vs September (GO) ----------
rows = []
for v in ['A_composite', 'B_nfe_primary', 'B_composite_404nfe_99fin', 'A_AF_nfe_only', 'B_nfe_PASS_only', 'A_gnomAD_PASS_only', 'B_nfe_naive_ref_1minusAF']:
    mv = allrep[allrep.variant == v].merge(codex, on=K, suffixes=('', '_c'))
    for _, r in mv.iterrows():
        d = dict(variant=v, gene=r.gene, endpoint=r.endpoint, p12=r.p12k, n=r.n, n_codex=r.n_c, top=r.top_snp, top_codex=r.top_snp_c)
        for h in H: d['d' + h[3:]] = r[h] - r[h + '_c']
        d['H3'] = r['PP.H3']; d['H4'] = r['PP.H4']; d['H2'] = r['PP.H2']
        d['H3_codex'] = r['PP.H3_c']; d['H4_codex'] = r['PP.H4_c']; d['H2_codex'] = r['PP.H2_c']
        # log-Bayes-factor decomposition (relative to H0)
        d['dL1'] = np.log(r['PP.H1'] / r['PP.H0']) - np.log(r['PP.H1_c'] / r['PP.H0_c'])
        d['dL2'] = np.log(r['PP.H2'] / r['PP.H0']) - np.log(r['PP.H2_c'] / r['PP.H0_c'])
        d['dL12'] = np.log(r['PP.H4'] / r['PP.H0']) - np.log(r['PP.H4_c'] / r['PP.H0_c'])
        d['dlog_H4_over_H3'] = np.log(r['PP.H4'] / r['PP.H3']) - np.log(r['PP.H4_c'] / r['PP.H3_c'])
        d['dlogit_H4'] = np.log(r['PP.H4'] / (1 - r['PP.H4'])) - np.log(r['PP.H4_c'] / (1 - r['PP.H4_c']))
        d['round3_H4_match'] = round(r['PP.H4'], 3) == round(r['PP.H4_c'], 3)
        d['round3_H3_match'] = round(r['PP.H3'], 3) == round(r['PP.H3_c'], 3)
        d['round3_H2_match'] = round(r['PP.H2'], 3) == round(r['PP.H2_c'], 3)
        d['side_080_same'] = (r['PP.H4'] >= 0.8) == (r['PP.H4_c'] >= 0.8)
        rows.append(d)
dd = pd.DataFrame(rows)
dd.to_csv(OUT / 'coloc_deltas_vs_codex.csv', index=False)
summ = []
for (v, ep), g in dd.groupby(['variant', dd.endpoint.eq('FinnGen_GO').map({True: 'identity_GO', False: 'new_GD_sGO'})]):
    summ.append(dict(variant=v, set=ep, max_abs_dH4=g.dH4.abs().max(), max_abs_dH3=g.dH3.abs().max(), max_abs_dH2=g.dH2.abs().max(),
                     max_abs_dPP=g[['dH0', 'dH1', 'dH2', 'dH3', 'dH4']].abs().max().max(), max_abs_dlogit_H4=g.dlogit_H4.abs().max(),
                     max_abs_dlog_H4_H3=g.dlog_H4_over_H3.abs().max(), all_top_equal=bool((g.top == g.top_codex).all()),
                     all_side080_same=bool(g.side_080_same.all()), n_round3_H4_mismatch=int((~g.round3_H4_match).sum()),
                     n_round3_H3_mismatch=int((~g.round3_H3_match).sum()), n_round3_H2_mismatch=int((~g.round3_H2_match).sum())))
summ = pd.DataFrame(summ); summ.to_csv(OUT / 'coloc_delta_summary.csv', index=False)
print(summ.to_string())
# endpoint-invariance of eQTL-only term dL2 (must be identical across endpoints for a pure eQTL-side difference)
inv = dd.groupby(['variant', 'gene']).dL2.agg(['min', 'max']); inv['range'] = inv['max'] - inv['min']
print('\n dL2 (eQTL-only log BF) range across endpoints/priors:\n', inv.to_string())
out['dL2_endpoint_invariance_max_range'] = float(inv['range'].max())

# ---------- 4. Recompute posteriors from each agent's per-SNP inputs with my own code ----------
rec = []; snpcmp = []
for g in GENES:
    for ep in EPS:
        fa = pd.read_csv(A / 'snp_level' / f'{g}__{ep}.csv')
        fb = pd.read_csv(B / f'snp_inputs__{g}__{ep}.csv.gz')
        z1a = fa.beta_assessed / fa.sebeta; v1a = fa.sebeta ** 2
        l1a = labf(z1a.values, v1a.values, 0.2)
        v2a = 1 / (2 * fa.f_assessed * (1 - fa.f_assessed) * (fa.NrSamples + fa.Zscore ** 2))
        l2a = labf(fa.Zscore.values, v2a.values, 0.15)
        l1b = labf(fb.z_disease.values, fb.V_disease.values, 0.2); l2b = labf(fb.z_eqtl.values, fb.V_eqtl.values, 0.15)
        mm = pd.DataFrame({'SNP': fa.SNP, 'z1a': z1a, 'v1a': v1a, 'Za': fa.Zscore, 'fa': fa.f_assessed}).merge(
            fb.rename(columns={'z_disease': 'z1b', 'V_disease': 'v1b', 'z_eqtl': 'Zb', 'f_assessed': 'fb'}), on='SNP', how='outer', indicator=True)
        both = mm[mm._merge == 'both']
        snpcmp.append(dict(gene=g, endpoint=ep, nA=len(fa), nB=len(fb), n_both=len(both), onlyA=mm.loc[mm._merge == 'left_only', 'SNP'].tolist(),
                           onlyB=mm.loc[mm._merge == 'right_only', 'SNP'].tolist(), max_rel_dz1=float(((both.z1a - both.z1b) / both.z1b.abs().clip(lower=1e-12)).abs().max()),
                           max_abs_dz1=float((both.z1a - both.z1b).abs().max()), max_rel_dv1=float(((both.v1a - both.v1b) / both.v1b).abs().max()),
                           max_abs_dZ=float((both.Za - both.Zb).abs().max()), sign_flips_z1=int((np.sign(both.z1a) != np.sign(both.z1b)).sum()),
                           median_abs_df=float((both.fa - both.fb).abs().median()), max_abs_df=float((both.fa - both.fb).abs().max())))
        for p12 in [1e-5, 5e-6, 1e-6]:
            pa = posteriors(l1a, l2a, p12); pb = posteriors(l1b, l2b, p12)
            ra = a[(a.gene == g) & (a.endpoint == ep) & (a.p12k == f'{p12:.0e}')].iloc[0]
            rb = b[(b.variant == 'B_nfe_primary') & (b.gene == g) & (b.endpoint == ep) & (b.p12k == f'{p12:.0e}')].iloc[0]
            rec.append(dict(gene=g, endpoint=ep, p12=p12, nA=len(l1a), nB=len(l1b),
                            maxdiff_A_recomputed_vs_reported=float(np.max(np.abs(pa - ra[H].values.astype(float)))),
                            maxdiff_B_recomputed_vs_reported=float(np.max(np.abs(pb - rb[H].values.astype(float)))),
                            topA=fa.SNP.iloc[int(np.argmax(l1a + l2a))], topB=fb.SNP.iloc[int(np.argmax(l1b + l2b))]))
rec = pd.DataFrame(rec); rec.to_csv(OUT / 'recompute_from_snp_inputs.csv', index=False)
snpcmp = pd.DataFrame(snpcmp); snpcmp.to_csv(OUT / 'A_vs_B_snp_inputs.csv', index=False)
print('\nRecompute max |diff| A:', rec.maxdiff_A_recomputed_vs_reported.max(), ' B:', rec.maxdiff_B_recomputed_vs_reported.max())
print(snpcmp.drop(columns=['onlyA', 'onlyB']).to_string()); print('onlyA/onlyB:', snpcmp[['gene', 'endpoint', 'onlyA', 'onlyB']].to_dict('records'))
out['recompute_max_diff_A'] = float(rec.maxdiff_A_recomputed_vs_reported.max()); out['recompute_max_diff_B'] = float(rec.maxdiff_B_recomputed_vs_reported.max())
(OUT / 'compare_summary_part1.json').write_text(json.dumps(out, indent=1, default=str))
print(json.dumps(out, indent=1, default=str))
