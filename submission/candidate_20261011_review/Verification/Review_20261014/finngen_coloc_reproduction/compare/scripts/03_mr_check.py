"""MR comparison. Reads raw FinnGen regional extracts saved by A and B (independent fetches), harmonises instruments
myself, recovers the September 1000G-EUR exposure scale from September's own FinnGen_GO IVW + leave-one-out outputs
(closed-form, no optimisation), checks implied MAFs are k/1006, and tests that scale out-of-sample on the new endpoints."""
import gzip, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm

S = Path('<scratch>'); OUT = S / 'compare_AB/results'
R = Path('/home/user/TEDGWAS')
SEPTP = R / 'archives/submission_candidates/candidate_20260911_review/provenance'
CODEX_MR = R / 'submission/candidate_20261011_review/reproducibility/work/reanalysis_20261002/mr_v1/MR_all_genes_all_scenarios.csv'
A4 = R / 'submission/candidate_20261011_review/Verification/Previous_allele_audit/Tables/Table_A4_MR_reproduction.csv'
A3 = R / 'submission/candidate_20261011_review/Verification/Previous_allele_audit/Tables/Table_A3_both_variants_five_disease_GWAS.csv'
GENES = ['TSHR', 'IGF1R', 'CTLA4']; EPS = ['FinnGen_GO', 'E4_GRAVES_STRICT', 'E4_GRAVES_OPHT_STRICT']
FILE_EP = {'FinnGen_GO': 'GRAVES_OPHT', 'E4_GRAVES_STRICT': 'E4_GRAVES_STRICT', 'E4_GRAVES_OPHT_STRICT': 'E4_GRAVES_OPHT_STRICT'}
WIN = {'TSHR': 'TSHR_chr14_78551662_83550983', 'IGF1R': 'IGF1R_chr15_96307070_101306900', 'CTLA4': 'CTLA4_chr2_201371096_206369831'}

iv = pd.read_csv(SEPTP / 'instruments_verified.csv'); iv = iv[iv.gene_symbol.isin(GENES)].copy()
print(iv[['gene_symbol', 'snp', 'effect_allele', 'other_allele', 'zscore', 'n_samples']].to_string())

def finn_rows(src, g, ep):
    if src == 'A':
        p = S / f'coloc_repro_A/raw_finngen/finngen_R12_{FILE_EP[ep]}_{WIN[g]}.tsv'; d = pd.read_csv(p, sep='\t')
    else:
        p = S / f'coloc_repro_B/data/finngen/{FILE_EP[ep]}__{g}.tsv.gz'; d = pd.read_csv(p, sep='\t').rename(columns={'#chrom': 'chrom'})
    return d

# ---- 1. harmonised outcome rows from both independent fetches ----
rows = []
for g in GENES:
    for ep in EPS:
        dA, dB = finn_rows('A', g, ep), finn_rows('B', g, ep)
        for _, r in iv[iv.gene_symbol == g].iterrows():
            rec = dict(gene=g, endpoint=ep, snp=r.snp, ea=r.effect_allele, oa=r.other_allele, z_exp=r.zscore, n_exp=r.n_samples)
            for src, d in [('A', dA), ('B', dB)]:
                exact = d[d.rsids.astype(str) == r.snp]
                alias = d[d.rsids.astype(str).str.contains(r.snp + '(,|$)', regex=True) & (d.rsids.astype(str) != r.snp)]
                rec[f'n_rows_{src}'] = len(exact); rec[f'n_alias_{src}'] = len(alias)
                if len(exact) == 1:
                    x = exact.iloc[0]
                    s = 1 if (x.alt == r.effect_allele and x.ref == r.other_allele) else (-1 if (x.ref == r.effect_allele and x.alt == r.other_allele) else 0)
                    rec[f'by_{src}'] = s * x.beta if s else np.nan; rec[f'sy_{src}'] = x.sebeta; rec[f'pos38_{src}'] = x.pos; rec[f'sign_{src}'] = s
            rows.append(rec)
oc = pd.DataFrame(rows)
oc['A_equals_B'] = (oc.by_A == oc.by_B) & (oc.sy_A == oc.sy_B) | (oc.n_rows_A.eq(0) & oc.n_rows_B.eq(0))
oc.to_csv(OUT / 'mr_outcome_rows_A_vs_B.csv', index=False)
print(oc.to_string())
# cross-check rs179252 against the allele-audit Table A3 (Codex package)
a3 = pd.read_csv(A3); a3 = a3[a3.snp == 'rs179252'].set_index('outcome')
for ep in EPS:
    r = oc[(oc.snp == 'rs179252') & (oc.endpoint == ep)].iloc[0]
    print('rs179252', ep, 'mine', r.by_A, r.sy_A, ' TableA3', a3.loc[ep, 'aligned_G_beta_logOR'], a3.loc[ep, 'aligned_G_se'])

# ---- 2. recover September exposure scale from September FinnGen_GO outputs (closed form) ----
mrp = pd.read_csv(SEPTP / 'MR_primary_canonical.csv'); mrp = mrp[(mrp.outcome == 'FinnGen_GO') & mrp.gene_symbol.isin(GENES)].set_index('gene_symbol')
loo = pd.read_csv(SEPTP / 'posthoc_20260914/leave_one_out_estimates.csv'); loo = loo[(loo.outcome == 'FinnGen_GO') & loo.gene.isin(GENES)]
go = oc[(oc.endpoint == 'FinnGen_GO') & oc.by_A.notna()].set_index('snp')
bexp = {}
# TSHR: Wald ratio
bexp['rs179252'] = go.loc['rs179252', 'by_A'] / mrp.loc['TSHR', 'beta']
# CTLA4: leave-one-out leaves single-IV Wald ratios
for _, r in loo[(loo.gene == 'CTLA4') & (loo.n_iv == 1)].iterrows():
    kept = [s for s in ['rs13030124', 'rs148849825'] if s != r.excluded_SNP][0]
    bexp[kept] = go.loc[kept, 'by_A'] / r.beta
# IGF1R: W_full - W_{-j} = bx_j^2/sy_j^2 ; SE_fixed = SE_RE / max(1, sqrt(Q/(k-1)))
def wsum(se, q, k):
    sef = se / max(1.0, np.sqrt(q / (k - 1))); return 1 / sef ** 2
full = loo[(loo.gene == 'IGF1R') & (loo.n_iv == 3)].iloc[0]; Wf = wsum(full.se, full.Q, 3)
for _, r in loo[(loo.gene == 'IGF1R') & (loo.n_iv == 2)].iterrows():
    s = r.excluded_SNP; Wm = wsum(r.se, r.Q, 2)
    z = iv.set_index('snp').loc[s, 'zscore']
    bexp[s] = np.sign(z) * go.loc[s, 'sy_A'] * np.sqrt(Wf - Wm)
ivs = iv.set_index('snp')
rec = []
for s, b in bexp.items():
    z, n = ivs.loc[s, 'zscore'], ivs.loc[s, 'n_samples']
    q = z * z / (2 * b * b * (n + z * z)); maf = (1 - np.sqrt(1 - 4 * q)) / 2
    rec.append(dict(snp=s, gene=ivs.loc[s, 'gene_symbol'], bexp_sept=b, implied_maf=maf, maf_x1006=maf * 1006, k_round=round(maf * 1006), k_resid=maf * 1006 - round(maf * 1006)))
scale = pd.DataFrame(rec); print(scale.to_string())
a4 = pd.read_csv(A4); print('Table A4 TSHR beta_exposure:', a4.beta_exposure.iloc[0], ' recovered:', bexp['rs179252'], ' diff', bexp['rs179252'] - a4.beta_exposure.iloc[0])

# ---- 3. MR fits ----
def fit(bx, by, sy):
    bx, by, sy = map(np.asarray, (bx, by, sy)); k = len(bx)
    W = np.sum(bx ** 2 / sy ** 2); b = np.sum(bx * by / sy ** 2) / W
    q = np.sum((by - b * bx) ** 2 / sy ** 2) if k > 1 else np.nan
    se = np.sqrt((max(1, q / (k - 1)) if k > 1 else 1) / W)
    return dict(k=k, beta=b, se=se, Q=q, p=2 * norm.sf(abs(b / se)), OR=np.exp(b), lo196=np.exp(b - 1.96 * se), hi196=np.exp(b + 1.96 * se),
                lo_q=np.exp(b - norm.ppf(.975) * se), hi_q=np.exp(b + norm.ppf(.975) * se))

# gnomAD frequency of effect allele from B's gnomAD extracts (AF_nfe, AF_fin by allele pair)
def gnomad_f(g, snp):
    d = pd.read_csv(S / f'coloc_repro_B/data/gnomad/{g}.tsv.gz', sep='\t')
    pos = int(ivs.loc[snp, 'pos_hg19']); ea, oa = ivs.loc[snp, 'effect_allele'], ivs.loc[snp, 'other_allele']
    x = d[(d.pos == pos) & (((d.ref == oa) & (d.alt == ea)) | ((d.ref == ea) & (d.alt == oa)))]
    if len(x) != 1: return np.nan, np.nan
    x = x.iloc[0]; nfe = x.AF_nfe if x.alt == ea else 1 - x.AF_nfe; fin = x.AF_fin if x.alt == ea else 1 - x.AF_fin
    return nfe, (404 * nfe + 99 * fin) / 503

codex = pd.read_csv(CODEX_MR); codex = codex[(codex.scenario == 'original_reference') & codex.gene_symbol.isin(GENES)].set_index(['gene_symbol', 'outcome'])
sept_scale = dict(bexp)
res = []
for g in GENES:
    for ep in EPS:
        d = oc[(oc.gene == g) & (oc.endpoint == ep) & oc.by_A.notna()]
        scen = {}
        scen['sept_1000G_scale'] = [sept_scale[s] for s in d.snp]
        fs = [gnomad_f(g, s) for s in d.snp]
        for nm, idx in [('gnomAD_NFE', 0), ('gnomAD_composite', 1)]:
            scen[nm] = [ivs.loc[s, 'zscore'] / np.sqrt(2 * f[idx] * (1 - f[idx]) * (ivs.loc[s, 'n_samples'] + ivs.loc[s, 'zscore'] ** 2)) for s, f in zip(d.snp, fs)]
        c = codex.loc[(g, ep)]
        for nm, bx in scen.items():
            r = fit(bx, d.by_A.values, d.sy_A.values)
            res.append(dict(gene=g, endpoint=ep, scale=nm, ivs=';'.join(d.snp), **r, codex_beta=c.beta, codex_se=c.se, codex_p=c.pvalue, codex_OR=c.OR,
                            codex_lo=c.CI_lower, codex_hi=c.CI_upper, codex_Q=c.Q, d_beta=r['beta'] - c.beta, d_se=r['se'] - c.se, rel_d_p=r['p'] / c.pvalue - 1,
                            d_lo_q=r['lo_q'] - c.CI_lower, d_hi_q=r['hi_q'] - c.CI_upper))
res = pd.DataFrame(res); res.to_csv(OUT / 'mr_recomputed_vs_codex.csv', index=False)
pd.set_option('display.width', 250)
print(res[['gene', 'endpoint', 'scale', 'k', 'OR', 'lo_q', 'hi_q', 'p', 'Q', 'codex_OR', 'codex_lo', 'codex_hi', 'codex_p', 'codex_Q', 'd_beta', 'd_se', 'rel_d_p', 'd_lo_q', 'd_hi_q']].to_string(float_format=lambda z: f'{z:.6g}'))
s = res[res.scale == 'sept_1000G_scale']
print('\nSept-scale max |d_beta|', s.d_beta.abs().max(), ' max |d_se|', s.d_se.abs().max(), ' max |rel d_p|', s.rel_d_p.abs().max(), ' max |dCI|', max(s.d_lo_q.abs().max(), s.d_hi_q.abs().max()))
# gnomAD frequencies vs recovered 1000G for each instrument
for snp in bexp:
    g = ivs.loc[snp, 'gene_symbol']; nfe, comp = gnomad_f(g, snp); m = scale.set_index('snp').loc[snp, 'implied_maf']
    print(snp, g, 'EA', ivs.loc[snp, 'effect_allele'], 'gnomAD NFE f_EA', round(nfe, 5), 'composite', round(comp, 5), ' 1000G MAF', round(m, 5), '(k=%d)' % round(m * 1006))
json.dump({k: float(v) for k, v in bexp.items()}, open(OUT / 'sept_exposure_scale_recovered.json', 'w'), indent=1)
