#!/usr/bin/env python3
"""Independent reproduction (variant A) of coloc.abf posteriors and MR estimates.

Genes: TSHR, IGF1R, CTLA4.  Endpoints (FinnGen R12):
  FinnGen_GO            = GRAVES_OPHT            (858 / 499,490)  identity control
  E4_GRAVES_STRICT      = FinnGen GD             (3,962 / 496,386)
  E4_GRAVES_OPHT_STRICT = strict GO              (753 / 499,595)

Expression: eQTLGen cis-eQTL full regions (GRCh37).  Expression SE reconstructed
from Z, N and a European reference frequency.  Variant A frequency:
  f_ALT = (404 * AF_nfe + 99 * AF_fin) / 503   from gnomAD v2.1.1 genomes,
matched to eQTLGen by GRCh37 position AND allele pair, then expressed for the
eQTLGen AssessedAllele.  (The September analysis used 1000 Genomes Phase 3 EUR,
n = 503 = 404 non-Finnish + 99 FIN; not reachable here.)

Written from scratch (numpy/scipy/pandas); no code from the previous pipeline is
imported or executed.

Usage: python -I 03_analysis.py WORKDIR
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm
from scipy.optimize import least_squares

GENES = ['TSHR', 'IGF1R', 'CTLA4']
ENDPOINTS = {  # label -> FinnGen R12 phenocode
    'FinnGen_GO': 'GRAVES_OPHT',
    'E4_GRAVES_STRICT': 'E4_GRAVES_STRICT',
    'E4_GRAVES_OPHT_STRICT': 'E4_GRAVES_OPHT_STRICT',
}
P12S = [1e-5, 5e-6, 1e-6]
P1 = P2 = 1e-4
W1, W2 = 0.2, 0.15
EQTL_DIR = Path('/home/user/TEDGWAS/TrackA_MR/data')
SEPT = Path('/home/user/TEDGWAS/archives/submission_candidates/candidate_20260911_review/provenance')
CODEX = Path('/home/user/TEDGWAS/submission/candidate_20261011_review/reproducibility/work/reanalysis_20261002')

# Codex Table 2 (Tables.md) values as printed, for the FinnGen rows.
TABLE2 = {
    ('TSHR', 'E4_GRAVES_STRICT'): '0.08 (0.05-0.12); P=3.19e-32',
    ('TSHR', 'FinnGen_GO'): '0.10 (0.04-0.24); P=2.82e-7',
    ('TSHR', 'E4_GRAVES_OPHT_STRICT'): '0.15 (0.06-0.39); P=9.02e-5',
    ('IGF1R', 'E4_GRAVES_STRICT'): '1.07 (0.86-1.33); P=0.539',
    ('IGF1R', 'FinnGen_GO'): '1.41 (0.85-2.33); P=0.182',
    ('IGF1R', 'E4_GRAVES_OPHT_STRICT'): '1.52 (0.92-2.49); P=0.0996',
    ('CTLA4', 'E4_GRAVES_STRICT'): '0.17 (0.08-0.36); P=2.15e-6',
    ('CTLA4', 'FinnGen_GO'): '0.17 (0.04-0.68); P=0.0118',
    ('CTLA4', 'E4_GRAVES_OPHT_STRICT'): '0.17 (0.05-0.64); P=0.00888',
}


# ----------------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------------
def labf(z, v, w):
    """Wakefield log approximate Bayes factor."""
    r = w ** 2 / (w ** 2 + v)
    return 0.5 * (np.log1p(-r) + r * z ** 2)


def logdiff(a, b):
    """log(exp(a) - exp(b)), a >= b."""
    return a + np.log1p(-np.exp(b - a))


def coloc_abf(l1, l2, p12):
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    lh = np.array([
        0.0,
        np.log(P1) + s1,
        np.log(P2) + s2,
        np.log(P1) + np.log(P2) + logdiff(s1 + s2, s12),
        np.log(p12) + s12,
    ])
    return np.exp(lh - logsumexp(lh))


def wald_or_ivw(bx, by, sy):
    bx, by, sy = map(np.asarray, (bx, by, sy))
    k = len(bx)
    if k == 1:
        b = by[0] / bx[0]
        se = sy[0] / abs(bx[0])
        q = np.nan
        method = 'Wald ratio'
    else:
        w = bx ** 2 / sy ** 2
        b = np.sum(bx * by / sy ** 2) / np.sum(w)
        se_f = np.sqrt(1 / np.sum(w))
        q = float(np.sum((by - b * bx) ** 2 / sy ** 2))
        se = se_f * max(1.0, np.sqrt(q / (k - 1)))
        method = 'IVW (multiplicative random effects)'
    p = 2 * norm.sf(abs(b / se))
    return dict(method=method, n_iv=k, beta=b, se=se, OR=np.exp(b), CI_low=np.exp(b - 1.96 * se),
                CI_high=np.exp(b + 1.96 * se), p=p, Q=q)


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------
def load_gnomad_freq(work, gene, eq):
    g = pd.read_csv(work / 'raw_gnomad' / gene / f'gnomad_r2.1.1_genomes_{gene}_eqtl_positions.tsv', sep='\t',
                    dtype={'chrom': str, 'ref': str, 'alt': str, 'id': str, 'filter': str},
                    na_values=['.'])
    for c in ['AF_nfe', 'AF_fin', 'AN_nfe', 'AN_fin', 'AF']:
        g[c] = pd.to_numeric(g[c], errors='coerce')
    g = g[(g.ref.str.len() == 1) & (g.alt.str.len() == 1)].copy()  # eQTLGen is SNV-only
    m = eq[['SNP', 'SNPPos', 'AssessedAllele', 'OtherAllele']].merge(g, left_on='SNPPos', right_on='pos', how='left')
    same = (m.AssessedAllele == m.alt) & (m.OtherAllele == m.ref)
    swap = (m.AssessedAllele == m.ref) & (m.OtherAllele == m.alt)
    m = m[same | swap].copy()
    m['assessed_is_alt'] = (m.AssessedAllele == m.alt)
    # QC: duplicates after allele matching (should be none for split biallelic records)
    dup = int(m.SNP.duplicated().sum())
    assert dup == 0, f'{gene}: {dup} duplicate gnomAD matches after allele matching'
    usable_an = (m.AN_nfe > 0) & (m.AN_fin > 0) & m.AF_nfe.notna() & m.AF_fin.notna()
    falt = (404 * m.AF_nfe + 99 * m.AF_fin) / 503
    m['f_alt_composite'] = np.where(usable_an, falt, np.nan)
    m['f_assessed'] = np.where(m.assessed_is_alt, m.f_alt_composite, 1 - m.f_alt_composite)
    # sensitivity frequencies
    m['f_assessed_nfe'] = np.where(m.AN_nfe > 0, np.where(m.assessed_is_alt, m.AF_nfe, 1 - m.AF_nfe), np.nan)
    m['gnomad_pass'] = m['filter'].eq('PASS')
    m['gnomad_rsid_agrees'] = [snp in str(i).split(';') for snp, i in zip(m.SNP, m['id'])]
    qc = dict(n_eqtl=len(eq), n_pos_allele_matched=len(m), n_freq_defined=int(usable_an.sum()),
              n_gnomad_nonpass=int((~m.gnomad_pass).sum()), n_rsid_disagree=int((~m.gnomad_rsid_agrees).sum()))
    return m[['SNP', 'f_assessed', 'f_assessed_nfe', 'gnomad_pass', 'filter', 'AF_nfe', 'AF_fin', 'AN_nfe',
              'AN_fin', 'id']].rename(columns={'id': 'gnomad_id'}), qc


def load_finngen(work, gene, phenocode):
    f = list((work / 'raw_finngen').glob(f'finngen_R12_{phenocode}_{gene}_chr*.tsv'))
    assert len(f) == 1, f
    d = pd.read_csv(f[0], sep='\t', dtype={'chrom': str, 'ref': str, 'alt': str, 'rsids': str})
    n_rows = len(d)
    d = d[d.rsids.notna()]
    alias = d.rsids.str.contains(',', regex=False)
    d_alias = d[alias]
    d = d[~alias].copy()
    d = d[np.isfinite(d.beta) & np.isfinite(d.sebeta) & (d.sebeta > 0)]
    return d, dict(window_rows=n_rows, alias_rows=int(alias.sum())), d_alias


def harmonise(eq, fg):
    m = eq.merge(fg, left_on='SNP', right_on='rsids', how='inner')
    n_before = len(m)
    n_rsid_matched = m.SNP.nunique()
    same = (m.AssessedAllele == m.alt) & (m.OtherAllele == m.ref)
    swap = (m.AssessedAllele == m.ref) & (m.OtherAllele == m.alt)
    m = m[same | swap].copy()
    m['beta_assessed'] = np.where(m.AssessedAllele == m.alt, m.beta, -m.beta)
    dup = int(m.SNP.duplicated().sum())
    return m, dict(n_rows_rsid_matched=n_before, n_unique_rsid_matched=n_rsid_matched,
                   n_allele_compatible=len(m), duplicates_after_alleles=dup)


# ----------------------------------------------------------------------------
def main():
    work = Path(sys.argv[1])
    out = work / 'results'
    (out / 'snp_level').mkdir(parents=True, exist_ok=True)

    sept_c = pd.read_csv(SEPT / 'coloc_canonical_v2.csv')
    sept_c = sept_c[(sept_c.outcome == 'FinnGen_GO') & sept_c.gene.isin(GENES)]
    codex_c = pd.read_csv(CODEX / 'coloc_v1' / 'Coloc_all_settings_with_identity.csv')
    codex_c = codex_c[(codex_c.scenario == 'original_reference') & codex_c.gene.isin(GENES)]
    codex_cov = pd.read_csv(CODEX / 'coloc_v1' / 'coverage.csv')
    sept_mr = pd.read_csv(SEPT / 'MR_primary_canonical.csv')
    sept_mr = sept_mr[(sept_mr.outcome == 'FinnGen_GO') & sept_mr.gene_symbol.isin(GENES)]
    codex_mr = pd.read_csv(CODEX / 'mr_v1' / 'MR_14genes_all_scenarios.csv')
    codex_mr = codex_mr[(codex_mr.scenario == 'original_reference') & codex_mr.gene_symbol.isin(GENES)]
    ivs = pd.read_csv(SEPT / 'instruments_verified.csv')
    ivs = ivs[ivs.gene_symbol.isin(GENES)].copy()

    coloc_rows, cov_rows, sens_rows, mr_rows, freq_qc = [], [], [], [], {}
    gnomad_freqs = {}
    for gene in GENES:
        eq = pd.read_csv(EQTL_DIR / f'eqtl_{gene}_all.csv')
        assert not eq.SNP.duplicated().any()
        fr, qc = load_gnomad_freq(work, gene, eq)
        freq_qc[gene] = qc
        gnomad_freqs[gene] = fr
        eqf = eq.merge(fr, on='SNP', how='left')
        eqf['f_valid'] = (eqf.f_assessed > 0) & (eqf.f_assessed < 1)
        f = eqf.f_assessed
        eqf['V2'] = 1.0 / (2 * f * (1 - f) * (eqf.NrSamples + eqf.Zscore ** 2))
        for label, pheno in ENDPOINTS.items():
            fg, fqc, _ = load_finngen(work, gene, pheno)
            h, hqc = harmonise(eqf, fg)
            assert hqc['duplicates_after_alleles'] == 0, (gene, label, hqc)
            d = h[h.f_valid].sort_values('SNP').copy()
            z1 = (d.beta_assessed / d.sebeta).to_numpy()
            v1 = (d.sebeta ** 2).to_numpy()
            z2 = d.Zscore.to_numpy()
            v2 = d.V2.to_numpy()
            assert np.isfinite(np.r_[z1, v1, z2, v2]).all()
            l1, l2 = labf(z1, v1, W1), labf(z2, v2, W2)
            top = d.SNP.iloc[int(np.argmax(l1 + l2))]
            d.assign(lABF1=l1, lABF2=l2, beta2=z2 * np.sqrt(v2)).to_csv(
                out / 'snp_level' / f'{gene}__{label}.csv', index=False,
                columns=['SNP', 'SNPPos', 'chrom', 'pos', 'AssessedAllele', 'OtherAllele', 'ref', 'alt', 'beta',
                         'sebeta', 'beta_assessed', 'Zscore', 'NrSamples', 'f_assessed', 'V2', 'beta2', 'lABF1',
                         'lABF2'])
            cc = codex_cov[(codex_cov.gene == gene) & (codex_cov.outcome == label)].iloc[0]
            cov_rows.append(dict(gene=gene, endpoint=label, finngen_phenocode=pheno, n_eqtl=len(eq), **fqc, **hqc,
                                 n_with_valid_freq=len(d), n_snps_used=len(d),
                                 n_allele_compatible_no_gnomad_freq=int((~h.f_valid).sum()),
                                 codex_n_matched_before_alleles=int(cc.n_matched_before_alleles),
                                 codex_n_allele_compatible=int(cc.n_allele_compatible),
                                 codex_n_ref=int(cc.n_ref)))
            for p12 in P12S:
                pp = coloc_abf(l1, l2, p12)
                coloc_rows.append(dict(gene=gene, endpoint=label, finngen_phenocode=pheno, p12=p12, n_snps=len(d),
                                       **{f'PP.H{i}': pp[i] for i in range(5)}, top_snp=top))
            # sensitivities on the frequency definition (not primary)
            for sname, mask, fcol in [('gnomAD_PASS_only', d.gnomad_pass, 'f_assessed'),
                                      ('AF_nfe_only', None, 'f_assessed_nfe')]:
                dd = h.copy()
                ff = dd[fcol]
                ok = (ff > 0) & (ff < 1)
                if sname == 'gnomAD_PASS_only':
                    ok &= dd.gnomad_pass.fillna(False).astype(bool)
                dd = dd[ok]
                vv = 1.0 / (2 * dd[fcol] * (1 - dd[fcol]) * (dd.NrSamples + dd.Zscore ** 2))
                a1 = labf((dd.beta_assessed / dd.sebeta).to_numpy(), (dd.sebeta ** 2).to_numpy(), W1)
                a2 = labf(dd.Zscore.to_numpy(), vv.to_numpy(), W2)
                for p12 in P12S:
                    pp = coloc_abf(a1, a2, p12)
                    sens_rows.append(dict(gene=gene, endpoint=label, sensitivity=sname, p12=p12, n_snps=len(dd),
                                          **{f'PP.H{i}': pp[i] for i in range(5)},
                                          top_snp=dd.SNP.iloc[int(np.argmax(a1 + a2))]))

    coloc = pd.DataFrame(coloc_rows)
    cov = pd.DataFrame(cov_rows)
    sens = pd.DataFrame(sens_rows)

    # ---- comparisons: September identity (FinnGen_GO) and Codex (all three endpoints)
    key = ['gene', 'endpoint', 'p12']
    s = sept_c.rename(columns={'outcome': 'endpoint'})[['gene', 'endpoint', 'p12', 'n_overlap', 'top_snp'] +
                                                         [f'PP.H{i}' for i in range(5)]]
    s = s.rename(columns={c: f'sept_{c}' for c in s.columns if c not in key})
    c = codex_c.rename(columns={'outcome': 'endpoint'})[['gene', 'endpoint', 'p12', 'n_overlap', 'top_snp'] +
                                                         [f'PP.H{i}' for i in range(5)]]
    c = c.rename(columns={c_: f'codex_{c_}' for c_ in c.columns if c_ not in key})
    coloc['p12_key'] = coloc.p12.map(lambda x: f'{x:.0e}')
    s['p12_key'] = s.p12.map(lambda x: f'{x:.0e}')
    c['p12_key'] = c.p12.map(lambda x: f'{x:.0e}')
    coloc = coloc.merge(s.drop(columns='p12'), on=['gene', 'endpoint', 'p12_key'], how='left') \
                 .merge(c.drop(columns='p12'), on=['gene', 'endpoint', 'p12_key'], how='left')
    for i in range(5):
        coloc[f'delta_vs_sept_PP.H{i}'] = coloc[f'PP.H{i}'] - coloc[f'sept_PP.H{i}']
        coloc[f'delta_vs_codex_PP.H{i}'] = coloc[f'PP.H{i}'] - coloc[f'codex_PP.H{i}']
    coloc = coloc.drop(columns='p12_key')
    coloc.to_csv(out / 'coloc_variantA.csv', index=False)
    cov.to_csv(out / 'coverage_variantA.csv', index=False)
    sens.to_csv(out / 'coloc_frequency_sensitivity_variantA.csv', index=False)

    # ---- MR
    for gene in GENES:
        g_iv = ivs[ivs.gene_symbol == gene].copy()
        fr = gnomad_freqs[gene].set_index('SNP')
        g_iv['f_effect_gnomad'] = np.nan
        for ix, r in g_iv.iterrows():
            # instrument effect allele frequency from the same gnomAD composite, oriented to effect_allele
            eqrow = pd.read_csv(EQTL_DIR / f'eqtl_{gene}_all.csv').set_index('SNP').loc[r.snp]
            fa = fr.loc[r.snp, 'f_assessed'] if r.snp in fr.index else np.nan
            if eqrow.AssessedAllele == r.effect_allele and eqrow.OtherAllele == r.other_allele:
                fe = fa
            elif eqrow.AssessedAllele == r.other_allele and eqrow.OtherAllele == r.effect_allele:
                fe = 1 - fa
            else:
                fe = np.nan
            g_iv.loc[ix, 'f_effect_gnomad'] = fe
        h = 2 * g_iv.f_effect_gnomad * (1 - g_iv.f_effect_gnomad)
        g_iv['se_exp'] = 1 / np.sqrt(h * (g_iv.n_samples + g_iv.zscore ** 2))
        g_iv['beta_exp'] = g_iv.zscore * g_iv.se_exp
        for label, pheno in ENDPOINTS.items():
            fg, _, fg_alias = load_finngen(work, gene, pheno)
            rows = []
            for _, r in g_iv.iterrows():
                cand = fg[fg.rsids == r.snp]
                hit = None
                for _, o in cand.iterrows():
                    if o.alt == r.effect_allele and o.ref == r.other_allele:
                        hit = (o.beta, o.sebeta, o.alt, o.ref, o.pos)
                    elif o.ref == r.effect_allele and o.alt == r.other_allele:
                        hit = (-o.beta, o.sebeta, o.alt, o.ref, o.pos)
                in_alias = bool(fg_alias.rsids.str.split(',').apply(lambda x: r.snp in x).any())
                rows.append(dict(gene=gene, endpoint=label, snp=r.snp, effect_allele=r.effect_allele,
                                 other_allele=r.other_allele, zscore=r.zscore, n_samples=r.n_samples,
                                 f_effect_gnomad=r.f_effect_gnomad, beta_exp=r.beta_exp, se_exp=r.se_exp,
                                 n_finngen_rows_same_rsid=len(cand), present_only_as_alias=in_alias and len(cand) == 0,
                                 retained=hit is not None,
                                 beta_out=hit[0] if hit else np.nan, se_out=hit[1] if hit else np.nan,
                                 finngen_alt=hit[2] if hit else None, finngen_ref=hit[3] if hit else None,
                                 finngen_pos_GRCh38=hit[4] if hit else np.nan))
            ivt = pd.DataFrame(rows)
            ivt.to_csv(out / 'snp_level' / f'MR_instruments__{gene}__{label}.csv', index=False)
            k = ivt[ivt.retained]
            res = wald_or_ivw(k.beta_exp, k.beta_out, k.se_out)
            mr_rows.append(dict(gene=gene, endpoint=label, finngen_phenocode=pheno, ivs=';'.join(k.snp),
                                excluded_ivs=';'.join(ivt[~ivt.retained].snp), **res))
            # Same-SNP check using the September exposure scale where it can be recovered exactly.
    mr = pd.DataFrame(mr_rows)

    # ---- recover September exposure betas from the September FinnGen_GO MR (pipeline check isolating frequency)
    loo = pd.read_csv(SEPT / 'leave_one_out_R_validation.csv')
    sept_bexp = {}
    for gene in GENES:
        ins = pd.read_csv(out / 'snp_level' / f'MR_instruments__{gene}__FinnGen_GO.csv')
        ins = ins[ins.retained].set_index('snp')
        sm = sept_mr[sept_mr.gene_symbol == gene].iloc[0]
        lg = loo[(loo.gene == gene) & (loo.outcome == 'FinnGen_GO')]
        snps = list(ins.index)
        sign = np.sign(ins.zscore.to_numpy())

        def model(log_abs_bx):
            bx = sign * np.exp(log_abs_bx)
            resid = []
            targets = [(None, sm.beta, sm.se)] + [(r.excluded_SNP, r.beta, r.se) for r in lg.itertuples()
                                                   if r.excluded_SNP.startswith('rs')]
            for excl, tb, tse in targets:
                keep = [i for i, sn in enumerate(snps) if sn != excl]
                rr = wald_or_ivw(bx[keep], ins.beta_out.to_numpy()[keep], ins.se_out.to_numpy()[keep])
                resid += [(rr['beta'] - tb) / abs(tb), (rr['se'] - tse) / tse]
            return np.array(resid)

        x0 = np.log(np.abs(ins.beta_exp.to_numpy()))
        sol = least_squares(model, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        bx = sign * np.exp(sol.x)
        sept_bexp[gene] = dict(zip(snps, bx))
        sept_bexp[gene + '_max_rel_resid'] = float(np.max(np.abs(sol.fun)))
    sept_scale_rows = []
    for gene in GENES:
        for label in ENDPOINTS:
            ins = pd.read_csv(out / 'snp_level' / f'MR_instruments__{gene}__{label}.csv')
            ins = ins[ins.retained]
            bx = np.array([sept_bexp[gene][sn] for sn in ins.snp])
            rr = wald_or_ivw(bx, ins.beta_out, ins.se_out)
            sept_scale_rows.append(dict(gene=gene, endpoint=label, sept_scale_beta=rr['beta'], sept_scale_se=rr['se'],
                                        sept_scale_OR=rr['OR'], sept_scale_CI_low=rr['CI_low'],
                                        sept_scale_CI_high=rr['CI_high'], sept_scale_p=rr['p'],
                                        sept_scale_bexp=';'.join(f'{sn}:{sept_bexp[gene][sn]:.6f}' for sn in ins.snp),
                                        sept_scale_fit_max_rel_resid=sept_bexp[gene + '_max_rel_resid']))
    mr = mr.merge(pd.DataFrame(sept_scale_rows), on=['gene', 'endpoint'])

    sm = sept_mr.rename(columns={'gene_symbol': 'gene', 'outcome': 'endpoint'})[
        ['gene', 'endpoint', 'n_iv', 'beta', 'se', 'or', 'or_ci_lower', 'or_ci_upper', 'pvalue']]
    sm.columns = ['gene', 'endpoint'] + ['sept_' + x for x in ['n_iv', 'beta', 'se', 'OR', 'CI_low', 'CI_high', 'p']]
    cm = codex_mr.rename(columns={'gene_symbol': 'gene', 'outcome': 'endpoint'})[
        ['gene', 'endpoint', 'n_iv', 'beta', 'se', 'OR', 'CI_lower', 'CI_upper', 'pvalue', 'Q']]
    cm.columns = ['gene', 'endpoint'] + ['codex_' + x for x in ['n_iv', 'beta', 'se', 'OR', 'CI_low', 'CI_high', 'p',
                                                                  'Q']]
    mr = mr.merge(sm, on=['gene', 'endpoint'], how='left').merge(cm, on=['gene', 'endpoint'], how='left')
    mr['codex_table2_reported'] = [TABLE2[(g, e)] for g, e in zip(mr.gene, mr.endpoint)]
    mr['delta_beta_vs_codex'] = mr.beta - mr.codex_beta
    mr['delta_OR_vs_codex'] = mr.OR - mr.codex_OR
    mr['ratio_p_vs_codex'] = mr.p / mr.codex_p
    mr['delta_sept_scale_beta_vs_codex'] = mr.sept_scale_beta - mr.codex_beta
    mr['delta_sept_scale_se_vs_codex'] = mr.sept_scale_se - mr.codex_se
    mr['delta_sept_scale_p_rel_vs_codex'] = mr.sept_scale_p / mr.codex_p - 1
    mr.to_csv(out / 'mr_variantA.csv', index=False)

    # ---- final combined CSV
    fin_c = coloc.assign(analysis='coloc')
    fin_m = mr.assign(analysis='MR')
    final = pd.concat([fin_c, fin_m], ignore_index=True, sort=False)
    final.to_csv(out / 'final_results_variantA.csv', index=False)

    summary = dict(
        freq_qc=freq_qc,
        identity_max_abs_delta_PP_vs_september=float(np.nanmax(np.abs(
            coloc.loc[coloc.endpoint == 'FinnGen_GO', [f'delta_vs_sept_PP.H{i}' for i in range(5)]].to_numpy()))),
        identity_max_abs_delta_H4_vs_september=float(np.nanmax(np.abs(
            coloc.loc[coloc.endpoint == 'FinnGen_GO', 'delta_vs_sept_PP.H4']))),
        max_abs_delta_H4_vs_codex_new_endpoints=float(np.nanmax(np.abs(
            coloc.loc[coloc.endpoint != 'FinnGen_GO', 'delta_vs_codex_PP.H4']))),
        max_abs_delta_any_PP_vs_codex_new_endpoints=float(np.nanmax(np.abs(
            coloc.loc[coloc.endpoint != 'FinnGen_GO', [f'delta_vs_codex_PP.H{i}' for i in range(5)]].to_numpy()))),
        max_abs_delta_H4_vs_codex_all=float(np.nanmax(np.abs(coloc['delta_vs_codex_PP.H4']))),
        sept_exposure_fit_residuals={g: sept_bexp[g + '_max_rel_resid'] for g in GENES},
    )
    (out / 'summary_variantA.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
