#!/usr/bin/env python3
"""Variant B, step 2: coloc.abf posteriors + cis-MR for TSHR/IGF1R/CTLA4 in three
FinnGen R12 endpoints, from the regional extracts written by 01_fetch_regions.py.

Independent implementation (written from scratch; log-sum-exp implemented here,
no scipy.special).  Frequency for the expression SE:
  primary      AF_nfe (gnomAD v2.1.1 genomes, non-Finnish Europeans)
  sensitivity  composite (404*AF_nfe + 99*AF_fin)/503  (1000G EUR composition proxy)
  sensitivity  AF_nfe restricted to gnomAD FILTER == PASS
Run:  python3 -I scripts/02_coloc_mr.py
"""
import gzip, json, math, os, sys, datetime
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
EQTL = "/home/user/TEDGWAS/TrackA_MR/data/eqtl_{g}_all.csv"
SEPT = "/home/user/TEDGWAS/archives/submission_candidates/candidate_20260911_review/provenance"
CODEX = "/home/user/TEDGWAS/submission/candidate_20261011_review/reproducibility/work/reanalysis_20261002"
GENES = ["TSHR", "IGF1R", "CTLA4"]
# file endpoint name -> label used by Codex / September tables
ENDPOINTS = {"GRAVES_OPHT": "FinnGen_GO", "E4_GRAVES_STRICT": "E4_GRAVES_STRICT",
             "E4_GRAVES_OPHT_STRICT": "E4_GRAVES_OPHT_STRICT"}
CASES = {"GRAVES_OPHT": (858, 499490), "E4_GRAVES_STRICT": (3962, 496386),
         "E4_GRAVES_OPHT_STRICT": (753, 499595)}
P1 = P2 = 1e-4
P12S = [1e-5, 5e-6, 1e-6]
W1, W2 = 0.2, 0.15
Z975 = 1.96
COMP = {"nfe": 404.0, "fin": 99.0}


# ----------------------------------------------------------------- numerics
def lse(x):
    """log(sum(exp(x))) computed directly (own implementation)."""
    x = np.asarray(x, dtype=float)
    m = float(np.max(x))
    if not math.isfinite(m):
        return m
    return m + math.log(float(np.sum(np.exp(x - m))))


def ldiff(a, b):
    """log(exp(a) - exp(b)) for a >= b."""
    if b > a:
        return float("-inf")
    return a + math.log1p(-math.exp(b - a)) if b > -math.inf else a


def labf(z, v, w):
    r = w * w / (w * w + v)
    return 0.5 * (np.log(1.0 - r) + r * z * z)


def coloc(z1, v1, z2, v2, p12):
    l1 = labf(z1, v1, W1); l2 = labf(z2, v2, W2); l12 = l1 + l2
    s1, s2, s12 = lse(l1), lse(l2), lse(l12)
    lh = np.array([0.0,
                   math.log(P1) + s1,
                   math.log(P2) + s2,
                   math.log(P1) + math.log(P2) + ldiff(s1 + s2, s12),
                   math.log(p12) + s12])
    pp = np.exp(lh - lse(lh))
    top = int(np.argmax(l12))
    return pp, top, float(math.exp(l12[top] - s12))


def norm_sf2(z):
    """two-sided normal P using erfc (no scipy)."""
    return math.erfc(abs(z) / math.sqrt(2.0))


# ----------------------------------------------------------------- inputs
def load_eqtl(g):
    e = pd.read_csv(EQTL.format(g=g), dtype={"SNP": str, "AssessedAllele": str, "OtherAllele": str})
    assert e.SNP.is_unique and e.GeneSymbol.eq(g).all() and e.SNPChr.nunique() == 1
    return e


def load_gnomad(g):
    d = pd.read_csv(os.path.join(BASE, "data", "gnomad", f"{g}.tsv.gz"), sep="\t",
                    dtype={"ref": str, "alt": str, "filter": str, "id": str, "chrom": str})
    for c in ["AF_nfe", "AF_fin", "AN_nfe", "AN_fin", "AC_nfe", "AC_fin", "AF", "AN"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def assessed_freq(e, gn):
    """Frequency of the eQTL AssessedAllele in gnomAD NFE and FIN, matched by GRCh37
    position + allele pair (either orientation).  Records are one ALT per line; when the
    assessed allele is the REF of a multi-allelic site, its frequency is 1 - sum(ALT AF)
    over all SNV/indel ALTs at that position sharing that REF."""
    gn = gn.copy()
    # per-position sum of ALT AF (same REF) for REF-allele frequency
    sums = gn.groupby(["pos", "ref"]).agg(sum_nfe=("AF_nfe", lambda s: s.sum(min_count=1)),
                                         sum_fin=("AF_fin", lambda s: s.sum(min_count=1)),
                                         n_alt_records=("alt", "size")).reset_index()
    gn = gn.merge(sums, on=["pos", "ref"], how="left")
    m = e[["SNP", "SNPPos", "AssessedAllele", "OtherAllele"]].merge(gn, left_on="SNPPos", right_on="pos", how="left")
    fwd = (m.AssessedAllele == m.alt) & (m.OtherAllele == m.ref)   # assessed = ALT
    rev = (m.AssessedAllele == m.ref) & (m.OtherAllele == m.alt)   # assessed = REF
    m = m[fwd | rev].copy()
    m["assessed_is_alt"] = ((m.AssessedAllele == m.alt) & (m.OtherAllele == m.ref))
    m["f_nfe"] = np.where(m.assessed_is_alt, m.AF_nfe, 1.0 - m.sum_nfe)
    m["f_fin"] = np.where(m.assessed_is_alt, m.AF_fin, 1.0 - m.sum_fin)
    m["f_nfe_naive"] = np.where(m.assessed_is_alt, m.AF_nfe, 1.0 - m.AF_nfe)
    # unusable AN -> missing
    m.loc[~(m.AN_nfe > 0), ["f_nfe", "f_nfe_naive"]] = np.nan
    m.loc[~(m.AN_fin > 0), "f_fin"] = np.nan
    dup = int(m.SNP.duplicated().sum())
    assert dup == 0, f"{dup} eQTL SNPs matched >1 gnomAD record after allele matching"
    # strand-flip diagnostic: eQTL SNPs at a gnomAD position that only match after complement
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    un = e[~e.SNP.isin(m.SNP)][["SNP", "SNPPos", "AssessedAllele", "OtherAllele"]].merge(
        gn[["pos", "ref", "alt"]], left_on="SNPPos", right_on="pos", how="inner")
    flip_only = un[(un.AssessedAllele.map(comp) == un.alt) & (un.OtherAllele.map(comp) == un.ref) |
                   (un.AssessedAllele.map(comp) == un.ref) & (un.OtherAllele.map(comp) == un.alt)].SNP.nunique()
    diag = dict(eqtl_snps=len(e), gnomad_pos_present=int(e.SNPPos.isin(gn.pos).sum()),
                allele_matched=len(m), assessed_is_ref=int((~m.assessed_is_alt).sum()),
                multiallelic_ref_assessed=int(((~m.assessed_is_alt) & (m.n_alt_records > 1)).sum()),
                strand_complement_only=int(flip_only), non_pass=int((m["filter"] != "PASS").sum()))
    out = m[["SNP", "filter", "f_nfe", "f_fin", "f_nfe_naive", "AN_nfe", "AN_fin", "n_alt_records"]].copy()
    out["f_comp"] = (COMP["nfe"] * out.f_nfe + COMP["fin"] * out.f_fin) / (COMP["nfe"] + COMP["fin"])
    return out, diag


def load_finngen(ep, g, chrom):
    d = pd.read_csv(os.path.join(BASE, "data", "finngen", f"{ep}__{g}.tsv.gz"), sep="\t",
                    dtype={"#chrom": str, "ref": str, "alt": str, "rsids": str, "nearest_genes": str})
    d = d.rename(columns={"#chrom": "chrom"})
    assert d.chrom.eq(str(chrom)).all()
    return d


def align(e, fg):
    """Join FinnGen rows to eQTL SNPs by exact rsID (single-ID rows only), check allele
    pairs, and express the disease beta per eQTL AssessedAllele."""
    fg = fg[fg.rsids.notna()]
    alias = fg[fg.rsids.str.contains(",", regex=False)]
    alias_hits = sum(1 for s in alias.rsids if set(s.split(",")) & set(e.SNP))
    single = fg[~fg.rsids.str.contains(",", regex=False)]
    m = e.merge(single, left_on="SNP", right_on="rsids", how="inner")
    n_before = len(m); n_rsid = m.SNP.nunique()
    same = (m.AssessedAllele == m.alt) & (m.OtherAllele == m.ref)
    swap = (m.AssessedAllele == m.ref) & (m.OtherAllele == m.alt)
    m = m[same | swap].copy()
    m["beta_d"] = np.where(m.AssessedAllele == m.alt, m.beta, -m.beta)
    m["flipped"] = m.AssessedAllele == m.ref
    dups = int(m.SNP.duplicated().sum())
    if dups:
        # never chosen by P value; report and keep none of the ambiguous SNPs
        bad = m[m.SNP.duplicated(keep=False)].SNP.unique()
        m = m[~m.SNP.isin(bad)].copy()
    m = m[np.isfinite(m.beta_d) & (m.sebeta > 0)].copy()
    cov = dict(n_eqtl=len(e), alias_rows_hitting_eqtl=alias_hits, n_matched_before_alleles=n_before,
               n_unique_rsid_matched=n_rsid, n_allele_compatible=len(m), n_flipped=int(m.flipped.sum()),
               duplicates_after_alleles=dups,
               n_palindromic=int(((m.AssessedAllele + m.OtherAllele).isin(["AT", "TA", "CG", "GC"])).sum()))
    return m, cov


# ----------------------------------------------------------------- MR
def ivw(bx, by, sy):
    bx, by, sy = map(np.asarray, (bx, by, sy))
    k = len(bx)
    if k == 1:
        b = by[0] / bx[0]; se = sy[0] / abs(bx[0])
        return dict(method="Wald ratio", n_iv=1, beta=b, se=se, Q=np.nan, Q_df=np.nan)
    w = 1.0 / sy ** 2
    b = float(np.sum(w * bx * by) / np.sum(w * bx * bx))
    se_f = math.sqrt(1.0 / np.sum(w * bx * bx))
    Q = float(np.sum(w * (by - b * bx) ** 2))
    se = se_f * max(1.0, math.sqrt(Q / (k - 1)))
    return dict(method="IVW (multiplicative random effects)", n_iv=k, beta=b, se=se, Q=Q, Q_df=k - 1)


def finish(r):
    r["p"] = norm_sf2(r["beta"] / r["se"])
    r["or"] = math.exp(r["beta"]); r["ci_low"] = math.exp(r["beta"] - Z975 * r["se"])
    r["ci_high"] = math.exp(r["beta"] + Z975 * r["se"])
    return r


def main():
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    sept = pd.read_csv(os.path.join(SEPT, "coloc_canonical_v2.csv"))
    sept = sept[sept.outcome.eq("FinnGen_GO") & sept.gene.isin(GENES)]
    codex = pd.read_csv(os.path.join(CODEX, "coloc_v1", "Coloc_all_settings_with_identity.csv"))
    codex = codex[codex.gene.isin(GENES) & codex.scenario.eq("original_reference")]
    codex_cov = pd.read_csv(os.path.join(CODEX, "coloc_v1", "coverage.csv"))
    codex_mr = pd.read_csv(os.path.join(CODEX, "mr_v1", "MR_14genes_all_scenarios.csv"))
    codex_mr = codex_mr[codex_mr.gene_symbol.isin(GENES) & codex_mr.scenario.eq("original_reference")]
    sept_mr = pd.read_csv(os.path.join(SEPT, "MR_primary_canonical.csv"))
    sept_mr = sept_mr[sept_mr.gene_symbol.isin(GENES) & sept_mr.outcome.eq("FinnGen_GO")]

    rows, covs, freq_diag, snp_tables = [], [], {}, {}
    fg_cache = {}
    for g in GENES:
        e = load_eqtl(g); chrom = int(e.SNPChr.iloc[0])
        fr, diag = assessed_freq(e, load_gnomad(g)); freq_diag[g] = diag
        for ep, lab in ENDPOINTS.items():
            fg = load_finngen(ep, g, chrom); fg_cache[(ep, g)] = fg
            m, cov = align(e, fg)
            m = m.merge(fr, on="SNP", how="left")
            m = m.sort_values("SNP").reset_index(drop=True)
            scen = {
                "B_nfe_primary": m.f_nfe,
                "B_composite_404nfe_99fin": m.f_comp,
                "B_nfe_PASS_only": m.f_nfe.where(m["filter"].eq("PASS")),
                "B_nfe_naive_ref_1minusAF": m.f_nfe_naive,
            }
            cov.update(gene=g, endpoint=lab, finngen_file=ep,
                       n_gnomad_allele_matched=int(m.f_nfe.notna().sum()) + int((m.f_nfe.isna() & m.SNP.isin(fr.SNP)).sum()),
                       n_no_gnomad_record=int((~m.SNP.isin(fr.SNP)).sum()))
            for sname, f in scen.items():
                ok = (f > 0) & (f < 1)
                d = m[ok].copy(); fv = f[ok].to_numpy()
                z1 = (d.beta_d / d.sebeta).to_numpy(); v1 = (d.sebeta ** 2).to_numpy()
                z2 = d.Zscore.to_numpy(); v2 = 1.0 / (2 * fv * (1 - fv) * (d.NrSamples.to_numpy() + z2 ** 2))
                assert np.isfinite(np.r_[z1, v1, z2, v2]).all()
                cov[f"n_used_{sname}"] = int(ok.sum())
                if sname == "B_nfe_primary":
                    snp_tables[(g, lab)] = pd.DataFrame(dict(SNP=d.SNP, z_disease=z1, V_disease=v1, z_eqtl=z2,
                                                             V_eqtl=v2, f_assessed=fv))
                for p12 in P12S:
                    pp, top, top_pp = coloc(z1, v1, z2, v2, p12)
                    rows.append(dict(gene=g, endpoint=lab, finngen_file=ep, freq_scenario=sname, p12=p12,
                                     n_snps=int(ok.sum()), top_snp=d.SNP.iloc[top], top_snp_pp_h4=top_pp,
                                     **{f"PP.H{i}": float(pp[i]) for i in range(5)}))
            covs.append(cov)
    res = pd.DataFrame(rows); cov = pd.DataFrame(covs)

    # ---------------- comparisons
    cx = codex.rename(columns={"outcome": "endpoint"})
    res = res.merge(cx[["gene", "endpoint", "p12", "n_overlap", "top_snp"] + [f"PP.H{i}" for i in range(5)]]
                    .rename(columns={"n_overlap": "codex_n_overlap", "top_snp": "codex_top_snp",
                                     **{f"PP.H{i}": f"codex_PP.H{i}" for i in range(5)}}),
                    on=["gene", "endpoint", "p12"], how="left")
    sp = sept.rename(columns={"outcome": "endpoint"})
    res = res.merge(sp[["gene", "endpoint", "p12", "n_overlap", "top_snp"] + [f"PP.H{i}" for i in range(5)]]
                    .rename(columns={"n_overlap": "sept_n_overlap", "top_snp": "sept_top_snp",
                                     **{f"PP.H{i}": f"sept_PP.H{i}" for i in range(5)}}),
                    on=["gene", "endpoint", "p12"], how="left")
    for i in range(5):
        res[f"delta_vs_codex_PP.H{i}"] = res[f"PP.H{i}"] - res[f"codex_PP.H{i}"]
        res[f"delta_vs_sept_PP.H{i}"] = res[f"PP.H{i}"] - res[f"sept_PP.H{i}"]
    res.to_csv(os.path.join(RES, "coloc_variantB_all.csv"), index=False)
    cov = cov.merge(codex_cov[codex_cov.gene.isin(GENES)].rename(columns={"outcome": "endpoint"})
                    .add_prefix("codex_").rename(columns={"codex_gene": "gene", "codex_endpoint": "endpoint"}),
                    on=["gene", "endpoint"], how="left")
    cov.to_csv(os.path.join(RES, "coverage_variantB.csv"), index=False)
    for (g, lab), t in snp_tables.items():
        t.to_csv(os.path.join(RES, f"snp_inputs__{g}__{lab}.csv.gz"), index=False)

    # ---------------- MR
    inst = pd.read_csv(os.path.join(SEPT, "instruments_verified.csv"))
    inst = inst[inst.gene_symbol.isin(GENES)].copy()
    mr_rows, mr_snp = [], []
    for g in GENES:
        e = load_eqtl(g)
        fr, _ = assessed_freq(e, load_gnomad(g))
        ig = inst[inst.gene_symbol.eq(g)].merge(e[["SNP", "AssessedAllele", "OtherAllele", "Zscore", "NrSamples"]],
                                                left_on="snp", right_on="SNP", how="left")
        assert ig.SNP.notna().all() and (ig.effect_allele == ig.AssessedAllele).all() and \
            np.allclose(ig.zscore, ig.Zscore) and (ig.n_samples == ig.NrSamples).all()
        ig = ig.merge(fr, on="SNP", how="left")
        for ep, lab in ENDPOINTS.items():
            fg = fg_cache[(ep, g)]
            fg = fg[fg.rsids.notna() & ~fg.rsids.str.contains(",", regex=False)]
            for fname, fcol in [("B_nfe_primary", "f_nfe"), ("B_composite_404nfe_99fin", "f_comp")]:
                bx, by, sy, used = [], [], [], []
                for _, s in ig.iterrows():
                    f = s[fcol]
                    sex = 1.0 / math.sqrt(2 * f * (1 - f) * (s.n_samples + s.zscore ** 2))
                    bexp = s.zscore * sex
                    o = fg[fg.rsids.eq(s.snp)]
                    o = o[((o.alt == s.effect_allele) & (o.ref == s.other_allele)) |
                          ((o.ref == s.effect_allele) & (o.alt == s.other_allele))]
                    status = "absent_in_finngen" if len(o) == 0 else ("ambiguous" if len(o) > 1 else "used")
                    rec = dict(gene=g, endpoint=lab, freq=fname, snp=s.snp, effect_allele=s.effect_allele,
                               other_allele=s.other_allele, f_effect=f, beta_exp=bexp, se_exp=sex, status=status)
                    if status == "used":
                        o = o.iloc[0]
                        bo = o.beta if o.alt == s.effect_allele else -o.beta
                        rec.update(finngen_ref=o.ref, finngen_alt=o.alt, finngen_pos38=int(o.pos),
                                   beta_out=bo, se_out=o.sebeta, p_out=o.pval, af_alt=o.af_alt)
                        bx.append(bexp); by.append(bo); sy.append(o.sebeta); used.append(s.snp)
                    mr_snp.append(rec)
                r = finish(ivw(bx, by, sy))
                r.update(gene=g, endpoint=lab, freq=fname, snps=";".join(used))
                mr_rows.append(r)
    mr = pd.DataFrame(mr_rows); mrs = pd.DataFrame(mr_snp)

    # ---- September-implied exposure scale (calibration, to isolate the frequency effect)
    # TSHR: beta 0.10546 stated in the task; also implied by Sept FinnGen_GO Wald SE.
    # CTLA4: single-SNP Wald estimates in Sept FinnGen_GO leave-one-out table.
    # IGF1R: solve the 3 exposure betas from Sept FinnGen_GO full + leave-one-out IVW betas.
    loo = pd.read_csv(os.path.join(SEPT, "posthoc_20260914", "leave_one_out_estimates.csv")).drop_duplicates()
    go = mrs[(mrs.endpoint == "FinnGen_GO") & (mrs.freq == "B_nfe_primary") & (mrs.status == "used")].set_index("snp")
    implied = {}
    t = sept_mr[sept_mr.gene_symbol.eq("TSHR")].iloc[0]
    implied["rs179252"] = dict(beta_exp=go.loc["rs179252", "se_out"] / t.se,
                               beta_exp_from_point=go.loc["rs179252", "beta_out"] / t.beta, stated=0.10546)
    for _, rr in loo[(loo.gene == "CTLA4") & (loo.outcome == "FinnGen_GO") & (loo.n_iv == 1)].iterrows():
        kept = [s for s in ["rs13030124", "rs148849825"] if s != rr.excluded_SNP][0]
        implied[kept] = dict(beta_exp=math.copysign(go.loc[kept, "se_out"] / rr.se, go.loc[kept, "beta_out"] / rr.beta),
                             beta_exp_from_point=go.loc[kept, "beta_out"] / rr.beta)
    ig_snps = ["rs2654980", "rs7181339", "rs59467480"]
    sub = loo[(loo.gene == "IGF1R") & (loo.outcome == "FinnGen_GO")]
    targets = [(ig_snps, float(sub[sub.excluded_SNP.str.startswith("None")].beta.iloc[0]))]
    for s in ig_snps:
        targets.append(([x for x in ig_snps if x != s], float(sub[sub.excluded_SNP == s].beta.iloc[0])))
    by = {s: go.loc[s, "beta_out"] for s in ig_snps}; sy = {s: go.loc[s, "se_out"] for s in ig_snps}
    x = np.array([go.loc[s, "beta_exp"] for s in ig_snps], float)  # start at gnomAD-NFE scale

    def resid(x):
        bx = dict(zip(ig_snps, x)); out = []
        for snps, tb in targets:
            num = sum(bx[s] * by[s] / sy[s] ** 2 for s in snps); den = sum(bx[s] ** 2 / sy[s] ** 2 for s in snps)
            out.append(num / den - tb)
        return np.array(out)
    for _ in range(200):  # Gauss-Newton with numerical Jacobian (no scipy)
        r0 = resid(x); J = np.zeros((len(r0), 3))
        for j in range(3):
            h = 1e-7 * max(1.0, abs(x[j])); xp = x.copy(); xp[j] += h; J[:, j] = (resid(xp) - r0) / h
        step, *_ = np.linalg.lstsq(J, -r0, rcond=None); x = x + step
        if np.max(np.abs(step)) < 1e-13:
            break
    for s, v in zip(ig_snps, x):
        implied[s] = dict(beta_exp=float(v), solve_max_abs_residual=float(np.max(np.abs(resid(x)))))
    inst_all = inst.set_index("snp")
    for s, d in implied.items():
        z, n = inst_all.loc[s, "zscore"], inst_all.loc[s, "n_samples"]
        q = z * z / (2 * d["beta_exp"] ** 2 * (n + z * z))  # = f(1-f)
        d["implied_maf_1000G"] = 0.5 - math.sqrt(max(0.0, 0.25 - q))
        fr_row = mrs[(mrs.snp == s) & (mrs.freq == "B_nfe_primary")].iloc[0]
        fc_row = mrs[(mrs.snp == s) & (mrs.freq == "B_composite_404nfe_99fin")].iloc[0]
        d["gnomad_nfe_f_effect"] = fr_row.f_effect; d["composite_f_effect"] = fc_row.f_effect
        d["beta_exp_nfe"] = fr_row.beta_exp; d["beta_exp_composite"] = fc_row.beta_exp
    # calibrated MR in all endpoints
    for ep, lab in ENDPOINTS.items():
        for g in GENES:
            u = mrs[(mrs.gene == g) & (mrs.endpoint == lab) & (mrs.freq == "B_nfe_primary") & (mrs.status == "used")]
            u = u[u.snp.isin(implied)]
            r = finish(ivw([implied[s]["beta_exp"] for s in u.snp], u.beta_out, u.se_out))
            r.update(gene=g, endpoint=lab, freq="Sept_implied_exposure_scale", snps=";".join(u.snp))
            mr = pd.concat([mr, pd.DataFrame([r])], ignore_index=True)

    cmr = codex_mr.rename(columns={"gene_symbol": "gene", "outcome": "endpoint"})
    mr = mr.merge(cmr[["gene", "endpoint", "n_iv", "beta", "se", "pvalue", "OR", "CI_lower", "CI_upper", "Q"]]
                  .add_prefix("codex_").rename(columns={"codex_gene": "gene", "codex_endpoint": "endpoint"}),
                  on=["gene", "endpoint"], how="left")
    smr = sept_mr.rename(columns={"gene_symbol": "gene", "outcome": "endpoint"})
    mr = mr.merge(smr[["gene", "endpoint", "n_iv", "beta", "se", "pvalue", "or"]]
                  .add_prefix("sept_").rename(columns={"sept_gene": "gene", "sept_endpoint": "endpoint"}),
                  on=["gene", "endpoint"], how="left")
    mr["delta_beta_vs_codex"] = mr.beta - mr.codex_beta
    mr["delta_or_vs_codex"] = mr["or"] - mr.codex_OR
    mr["ratio_p_vs_codex"] = mr.p / mr.codex_pvalue
    mr.to_csv(os.path.join(RES, "mr_variantB.csv"), index=False)
    mrs.to_csv(os.path.join(RES, "mr_snp_level_variantB.csv"), index=False)
    with open(os.path.join(RES, "diagnostics_variantB.json"), "w") as fh:
        json.dump(dict(created_utc=now, numpy=np.__version__, pandas=pd.__version__, python=sys.version,
                       gnomad_frequency_matching=freq_diag, september_implied_exposure=implied), fh, indent=1,
                  default=float)
    print("done", len(res), "coloc rows;", len(mr), "MR rows")


if __name__ == "__main__":
    main()
