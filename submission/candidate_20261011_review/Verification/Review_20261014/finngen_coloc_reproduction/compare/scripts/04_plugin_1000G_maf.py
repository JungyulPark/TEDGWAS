"""Partial frequency substitution: put the recovered September 1000G-EUR MAFs (instrument SNPs only) into the gnomAD-composite
per-SNP coloc inputs and measure how much of the Codex-vs-reproduction residual disappears. For TSHR also restrict to the
retained 1000G EUR LD-reference SNP set (identified in 02)."""
import json
from pathlib import Path
import numpy as np, pandas as pd

S = Path('<scratch>'); A = S / 'coloc_repro_A/results'; OUT = S / 'compare_AB/results'
R = Path('/home/user/TEDGWAS')
CODEX = R / 'submission/candidate_20261011_review/reproducibility/work/reanalysis_20261002/coloc_v1/Coloc_all_settings_with_identity.csv'
LDREF = R / 'TrackA_MR/v5_upgrade/02_tshr_finemap/results/TaskB_04c_allele_harmonisation_20260904.csv'
IV = R / 'archives/submission_candidates/candidate_20260911_review/provenance/instruments_verified.csv'
GENES = ['TSHR', 'IGF1R', 'CTLA4']; EPS = ['FinnGen_GO', 'E4_GRAVES_STRICT', 'E4_GRAVES_OPHT_STRICT']; P12 = [1e-5, 5e-6, 1e-6]
H = [f'PP.H{i}' for i in range(5)]; LP1 = np.log(1e-4)
MAF1KG = {'rs179252': 0.494, 'rs148849825': 0.01988, 'rs13030124': 0.4503, 'rs2654980': 0.1839, 'rs59467480': 0.2724, 'rs7181339': 0.04771}
def lse(x):
    m = np.max(x); return float(m + np.log(np.sum(np.exp(x - m))))
def labf(z, v, w):
    r = w * w / (w * w + v); return 0.5 * (np.log1p(-r) + r * z * z)
def pp(l1, l2, p12):
    s1, s2, s12 = lse(l1), lse(l2), lse(l1 + l2); h3 = s1 + s2 + np.log(-np.expm1(s12 - s1 - s2))
    lh = np.array([0, LP1 + s1, LP1 + s2, 2 * LP1 + h3, np.log(p12) + s12]); return np.exp(lh - lse(lh))
codex = pd.read_csv(CODEX); codex = codex[(codex.scenario == 'original_reference') & codex.gene.isin(GENES)]
iv = pd.read_csv(IV).set_index('snp')
ld = pd.read_csv(LDREF); eur = set(ld.loc[ld.ld_reference == 'EUR', 'snp'])
rows = []
for g in GENES:
    for ep in EPS:
        f = pd.read_csv(A / 'snp_level' / f'{g}__{ep}.csv')
        l1 = labf((f.beta_assessed / f.sebeta).values, (f.sebeta ** 2).values, 0.2)
        fa = f.f_assessed.values.copy(); fsub = fa.copy(); nsub = 0
        for s, m in MAF1KG.items():
            i = np.where(f.SNP.values == s)[0]
            if len(i) == 1 and iv.loc[s, 'gene_symbol'] == g:
                # orient minor-allele frequency to the eQTL AssessedAllele using the gnomAD side to choose minor/major
                cur = fa[i[0]]; fsub[i[0]] = m if cur < 0.5 else 1 - m
                # sanity: for rs179252 (MAF 0.494) orientation is ambiguous; use EAF(G)=0.506 from Table A4, G is AssessedAllele?
                if s == 'rs179252': fsub[i[0]] = 0.506 if f.AssessedAllele.values[i[0]] == 'G' else 0.494
                nsub += 1
        Z, N = f.Zscore.values, f.NrSamples.values
        sets = {'all': np.ones(len(f), bool)}
        if g == 'TSHR': sets['EUR_ldref'] = f.SNP.isin(eur).values
        for sn, msk in sets.items():
            for fn, fv in [('gnomAD_composite', fa), ('composite+1000G_instrument_MAF', fsub)]:
                l2 = labf(Z, 1 / (2 * fv * (1 - fv) * (N + Z * Z)), 0.15)
                for p12 in P12:
                    c = codex[(codex.gene == g) & (codex.outcome == ep) & np.isclose(codex.p12, p12)].iloc[0]
                    v = pp(l1[msk], l2[msk], p12)
                    rows.append(dict(gene=g, endpoint=ep, p12=p12, snpset=sn, n=int(msk.sum()), freq=fn, n_substituted=nsub,
                                     H2=v[2], H3=v[3], H4=v[4], dH2=v[2] - c['PP.H2'], dH3=v[3] - c['PP.H3'], dH4=v[4] - c['PP.H4'],
                                     max_abs_dPP=float(np.max(np.abs(v - c[H].values.astype(float)))),
                                     dL2=lse(l2[msk]) - (np.log(c['PP.H2'] / c['PP.H0']) - LP1)))
d = pd.DataFrame(rows); d.to_csv(OUT / 'plugin_1000G_instrument_maf.csv', index=False)
pd.set_option('display.width', 250)
print(d.to_string(float_format=lambda z: f'{z:.4g}'))
print(d.groupby(['gene', 'snpset', 'freq']).agg(max_abs_dPP=('max_abs_dPP', 'max'), max_abs_dH4=('dH4', lambda x: x.abs().max()), dL2=('dL2', 'first')).to_string())
