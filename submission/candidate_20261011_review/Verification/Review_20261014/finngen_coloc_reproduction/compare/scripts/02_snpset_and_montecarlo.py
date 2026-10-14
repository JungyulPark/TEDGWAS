"""(a) TSHR: identify Codex/September 1000G-EUR SNP set via the retained 1000G EUR LD-reference SNP list and test the
frequency-free disease-only Bayes factor sum S1 against Codex. (b) Monte Carlo calibration: what posterior deviation does
1000G-like sampling (503 EUR = 1006 chromosomes) around the gnomAD composite frequency produce, jointly across endpoints?"""
import json
from pathlib import Path
import numpy as np, pandas as pd

S = Path('<scratch>')
A = S / 'coloc_repro_A/results'; OUT = S / 'compare_AB/results'
R = Path('/home/user/TEDGWAS')
CODEX = R / 'submission/candidate_20261011_review/reproducibility/work/reanalysis_20261002/coloc_v1/Coloc_all_settings_with_identity.csv'
LDREF = R / 'TrackA_MR/v5_upgrade/02_tshr_finemap/results/TaskB_04c_allele_harmonisation_20260904.csv'
GENES = ['TSHR', 'IGF1R', 'CTLA4']; EPS = ['FinnGen_GO', 'E4_GRAVES_STRICT', 'E4_GRAVES_OPHT_STRICT']; P12 = [1e-5, 5e-6, 1e-6]
H = [f'PP.H{i}' for i in range(5)]; LP1 = np.log(1e-4)

def lse(x, axis=None):
    m = np.max(x, axis=axis, keepdims=True); return (m + np.log(np.sum(np.exp(x - m), axis=axis, keepdims=True))).squeeze(axis) if axis is not None else float(m.squeeze() + np.log(np.sum(np.exp(x - m))))
def labf(z, v, w):
    r = w * w / (w * w + v); return 0.5 * (np.log1p(-r) + r * z * z)
def pp_from_sums(s1, s2, s12, p12):
    h3 = s1 + s2 + np.log(-np.expm1(s12 - s1 - s2))
    lh = np.stack([np.zeros_like(s1), LP1 + s1, LP1 + s2, 2 * LP1 + h3, np.log(p12) + s12], axis=-1)
    return np.exp(lh - lse(lh, axis=-1)[..., None])

codex = pd.read_csv(CODEX); codex = codex[(codex.scenario == 'original_reference') & codex.gene.isin(GENES)]
def codex_row(g, ep, p12):
    return codex[(codex.gene == g) & (codex.outcome == ep) & np.isclose(codex.p12, p12)].iloc[0]

def load(g):
    d = {}
    for ep in EPS:
        f = pd.read_csv(A / 'snp_level' / f'{g}__{ep}.csv').sort_values('SNP').reset_index(drop=True)
        d[ep] = f
    snps = d[EPS[0]].SNP
    assert all((d[e].SNP.values == snps.values).all() for e in EPS)
    base = d[EPS[0]]
    l1 = {e: labf((d[e].beta_assessed / d[e].sebeta).values, (d[e].sebeta ** 2).values, 0.2) for e in EPS}
    return base, l1

res = {}
# ---------------- (a) TSHR SNP set ----------------
base, l1 = load('TSHR')
ld = pd.read_csv(LDREF); eur = set(ld.loc[ld.ld_reference == 'EUR', 'snp'])
inset = base.SNP.isin(eur).values
res['TSHR_ldref'] = dict(n_A=len(base), n_in_EUR_ldref=int(inset.sum()), codex_n=int(codex_row('TSHR', 'FinnGen_GO', 1e-5).n_overlap))
for ep in EPS:
    c = codex_row('TSHR', ep, 1e-5)
    L1c = np.log(c['PP.H1'] / c['PP.H0']) - LP1
    res['TSHR_ldref'][ep] = dict(L1_codex=L1c, L1_A_all=lse(l1[ep]), L1_A_ldrefset=lse(l1[ep][inset]),
                                 dL1_all=lse(l1[ep]) - L1c, dL1_ldrefset=lse(l1[ep][inset]) - L1c)
# posteriors restricted to the LD-reference set but with gnomAD composite frequencies -> isolates frequency-value effect
f = base.f_assessed.values; Z = base.Zscore.values; N = base.NrSamples.values
l2 = labf(Z, 1 / (2 * f * (1 - f) * (N + Z * Z)), 0.15)
rows = []
for ep in EPS:
    for p12 in P12:
        c = codex_row('TSHR', ep, p12)
        for nm, msk in [('all_gnomAD_SNPs', np.ones_like(inset)), ('EUR_ldref_subset', inset)]:
            pp = pp_from_sums(np.array(lse(l1[ep][msk])), np.array(lse(l2[msk])), np.array(lse(l1[ep][msk] + l2[msk])), p12)
            rows.append(dict(gene='TSHR', endpoint=ep, p12=p12, snpset=nm, n=int(msk.sum()), H2=pp[2], H3=pp[3], H4=pp[4],
                             dH3_vs_codex=pp[3] - c['PP.H3'], dH4_vs_codex=pp[4] - c['PP.H4'],
                             dL2_vs_codex=lse(l2[msk]) - (np.log(c['PP.H2'] / c['PP.H0']) - LP1)))
tshr_set = pd.DataFrame(rows); tshr_set.to_csv(OUT / 'TSHR_snpset_decomposition.csv', index=False)
print(json.dumps(res, indent=1, default=float)); print(tshr_set.to_string())

# ---------------- (b) Monte Carlo ----------------
rng = np.random.default_rng(20261014)
D = 4000
mc_rows = []; mc_summary = []
for g in GENES:
    base, l1 = load(g)
    f = base.f_assessed.values; Z = base.Zscore.values; N = base.NrSamples.values
    l2_0 = labf(Z, 1 / (2 * f * (1 - f) * (N + Z * Z)), 0.15)
    pp0 = {(ep, p): pp_from_sums(np.array(lse(l1[ep])), np.array(lse(l2_0)), np.array(lse(l1[ep] + l2_0)), p) for ep in EPS for p in P12}
    k = rng.binomial(1006, np.clip(f, 0, 1), size=(D, len(f))); fs = k / 1006.0
    keep = (fs > 0) & (fs < 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        l2s = labf(Z[None, :], 1 / (2 * fs * (1 - fs) * (N + Z * Z)[None, :]), 0.15)
    l2s = np.where(keep, l2s, -np.inf)
    s2 = lse(l2s, axis=1)
    ndrop = (~keep).sum(1)
    sim = {}
    for ep in EPS:
        s1 = lse(np.where(keep, l1[ep][None, :], -np.inf), axis=1)
        s12 = lse(l1[ep][None, :] + l2s, axis=1)
        for p in P12:
            sim[(ep, p)] = pp_from_sums(s1, s2, s12, p)  # D x 5
    for ep in EPS:
        for p in P12:
            c = codex_row(g, ep, p)
            dsim = sim[(ep, p)] - pp0[(ep, p)][None, :]
            obs = np.array([c[h] for h in H]) - pp0[(ep, p)]
            for i in (2, 3, 4):
                lo, hi = np.quantile(dsim[:, i], [0.005, 0.995]); lo95, hi95 = np.quantile(dsim[:, i], [0.025, 0.975])
                mc_rows.append(dict(gene=g, endpoint=ep, p12=p, PP=f'H{i}', A_composite=pp0[(ep, p)][i], codex=c[H[i]], obs_delta_codex_minus_A=obs[i],
                                    sim_q025=lo95, sim_q975=hi95, sim_q005=lo, sim_q995=hi, sim_sd=dsim[:, i].std(),
                                    obs_percentile=float((dsim[:, i] < obs[i]).mean()), within_95=bool(lo95 <= obs[i] <= hi95), within_99=bool(lo <= obs[i] <= hi)))
    # joint consistency: regress new-endpoint dH4 on GO dH4 (same prior) across draws; predict from observed GO delta
    for ep in EPS[1:]:
        for p in P12:
            x = sim[('FinnGen_GO', p)][:, 4] - pp0[('FinnGen_GO', p)][4]; y = sim[(ep, p)][:, 4] - pp0[(ep, p)][4]
            xo = codex_row(g, 'FinnGen_GO', p)['PP.H4'] - pp0[('FinnGen_GO', p)][4]; yo = codex_row(g, ep, p)['PP.H4'] - pp0[(ep, p)][4]
            bcoef = np.polyfit(x, y, 1); resid = y - np.polyval(bcoef, x); pred = np.polyval(bcoef, xo)
            mc_summary.append(dict(gene=g, endpoint=ep, p12=p, corr_GO_vs_ep=float(np.corrcoef(x, y)[0, 1]), observed_GO_delta=xo, observed_ep_delta=yo,
                                   predicted_ep_delta_given_GO=float(pred), resid_sd=float(resid.std()), z_resid=float((yo - pred) / resid.std())))
    print(g, 'simulated dropped SNPs (mean, 5-95%):', ndrop.mean(), np.quantile(ndrop, [0.05, 0.95]))
    res[f'{g}_mc_dropped'] = dict(mean=float(ndrop.mean()), q05=float(np.quantile(ndrop, 0.05)), q95=float(np.quantile(ndrop, 0.95)), n_base=len(f))
mc = pd.DataFrame(mc_rows); mc.to_csv(OUT / 'montecarlo_1000G_sampling_vs_observed.csv', index=False)
mcs = pd.DataFrame(mc_summary); mcs.to_csv(OUT / 'montecarlo_joint_GO_conditioning.csv', index=False)
pd.set_option('display.width', 250)
print(mc.to_string(float_format=lambda z: f'{z:.4g}'))
print(mcs.to_string(float_format=lambda z: f'{z:.4g}'))
(OUT / 'compare_summary_part2.json').write_text(json.dumps(res, indent=1, default=float))
