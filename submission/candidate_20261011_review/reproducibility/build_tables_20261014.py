#!/usr/bin/env python3
"""Generate the table Markdown sources from the workbook (2026-10-14, after the
three-reviewer check of the displays).

Writes Tables.md (Tables 1-3) and SupplementaryTables.md (Tables S1-S4), numbered
in order of first citation in the manuscript:
  Table 1  data sources            Table S1 additional FinnGen GD discoveries (was 3B)
  Table 2  TSHR / IGF1R / CTLA4    Table S2 frequency-sensitive candidates
  Table 3  13 BBJ discoveries      Table S3 TSHR regional and prior sensitivity (was S1)
                                   Table S4 tissue extension (was S3)
Every number comes from Supplementary_Dataset.xlsx or the retained aggregate MR
export; nothing is typed by hand. Exponents use pandoc superscript (^−14^) so Word
sets them as true superscripts.

    python3 submission/candidate_20261011_review/reproducibility/build_tables_20261014.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
X = PKG / "Supplementary_Dataset.xlsx"
Q14 = PKG / "reproducibility/work/reanalysis_20261002/mr_v1/MR_14genes_all_scenarios.csv"
TIS = PKG / "Verification/Tissue_extension/Tables"
ORDER = ["BBJ_Graves", "E4_GRAVES_STRICT", "UKB_hyperthyroid", "FinnGen_GO", "E4_GRAVES_OPHT_STRICT"]
LABEL = {"BBJ_Graves": "BBJ GD", "E4_GRAVES_STRICT": "FinnGen GD", "UKB_hyperthyroid": "UKB hyperthyroidism",
         "FinnGen_GO": "FinnGen GO", "E4_GRAVES_OPHT_STRICT": "FinnGen GO strict"}
SCREEN = 0.05 / 2544
MHC = {"HLA-A", "HLA-DQA2", "C4A", "TUBB", "PSMB8", "APOM", "NEU1", "TNF", "HIST1H4E", "LTA", "LTB", "C4B"}
CLUSTER16 = {"HSD3B7", "VKORC1", "PRSS36"}


def sci(x, digits=3):
    m, e = f"{x:.{digits - 1}e}".split("e")
    return f"{m}×10^−{abs(int(e))}^" if int(e) < 0 else f"{m}×10^{int(e)}^"


def fp(x):
    """P value: 3 significant figures, trailing zeros kept."""
    if pd.isna(x):
        return "NA"
    if x < 0.001:
        return sci(x)
    d = 2 - int(np.floor(np.log10(x)))
    return f"{x:.{d}f}"


def pp(x, bold=False):
    s = "<0.001" if x < 0.001 else f"{x:.3f}"
    return f"**{s}**" if bold and x >= 0.8 else s


def orci(o, lo, hi, d=2):
    return f"{o:.{d}f} ({lo:.{d}f}–{hi:.{d}f})"


def table(cols, rows, w):
    """Pipe table; pandoc sets relative column widths from the dash counts in w."""
    assert len(w) == len(cols)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("-" * n for n in w) + "|"]
    out += ["| " + " | ".join(str(v) for v in r) + " |" for r in rows]
    return "\n".join(out)


def load():
    d = {s: pd.read_excel(X, sheet_name=s) for s in
         ["MR_primary", "Colocalization", "MR_estimators_original", "Instruments", "Region_coverage", "TSHR_windows",
          "Frequency_diagnostics"]}
    d["q14"] = pd.read_csv(Q14)
    a = pd.read_excel(X, sheet_name="TSHR_allele_audit", header=None).values.tolist()
    hi = next(i for i, r in enumerate(a) if "aligned_G_OR" in [str(v) for v in r])
    hdr = [str(v) for v in a[hi]]
    rows = []
    for r in a[hi + 1:]:
        if pd.isna(r[0]) or "outcome" in [str(v) for v in r]:
            break
        rows.append(dict(zip(hdr, r)))
    d["allele"] = pd.DataFrame(rows)
    return d


def mrrow(d, g, o):
    return d["MR_primary"][(d["MR_primary"].gene_symbol == g) & (d["MR_primary"].outcome == o)].iloc[0]


def coloc(d, g, o, p12=1e-5):
    c = d["Colocalization"]
    return c[(c.gene == g) & (c.outcome == o) & (c.scenario == "original_reference") & np.isclose(c.p12, p12)].iloc[0]


def qp(d, g, o):
    if o in ("E4_GRAVES_STRICT", "E4_GRAVES_OPHT_STRICT"):
        r = d["q14"][(d["q14"].gene_symbol == g) & (d["q14"].outcome == o) & (d["q14"].scenario == "original_reference")]
        v = r.Q_pvalue.iloc[0] if len(r) else np.nan
    else:
        e = d["MR_estimators_original"]
        r = e[(e.gene_symbol == g) & (e.outcome == o) & e.method.isin(["Inverse variance weighted", "Wald ratio"])]
        v = r.cochran_q_p.iloc[0] if len(r) else np.nan
    return "NA" if pd.isna(v) else fp(v)


def chrom(d, g):
    c = int(d["Instruments"][d["Instruments"].gene_symbol == g].chr.iloc[0])
    if g in MHC:
        return f"{c} (MHC)"
    if g in CLUSTER16:
        return f"{c} (16p11.2 cluster)"
    return str(c)


ABBR_DATA = "BBJ, Biobank Japan; FinnGen GD, E4_GRAVES_STRICT; GD, Graves disease; GO, Graves ophthalmopathy; UKB, UK Biobank"
CONTROLS = ("FinnGen ophthalmopathy controls are all other participants, including people with GD but no recorded eye "
            "disease, so GO estimates do not estimate TED risk among patients with GD.")


def table1():
    rows = [
        ["eQTLGen blood cis-eQTL", "—", "—", "Exposure; up to 31,684 participants (median 28,092 per variant); predominantly European"],
        ["BBJ GD (GCST90018627)", "2,809", "172,656", "Discovery; East Asian"],
        ["FinnGen GD (E4_GRAVES_STRICT)", "3,962", "496,386", "Post hoc same-disease comparison; Finnish"],
        ["UKB hyperthyroidism (GCST90038636)", "3,731", "480,867",
         "Broader-phenotype comparison; self-reported hyperthyroidism or thyrotoxicosis; ancestry not established as exclusively European"],
        ["FinnGen GO (GRAVES_OPHT)", "858", "499,490", "Original ophthalmopathy (TED-enriched) comparison; Finnish"],
        ["FinnGen GO strict (E4_GRAVES_OPHT_STRICT)", "753", "499,595", "Post hoc ophthalmopathy definition sensitivity; Finnish"],
        ["GTEx v8 tissue eQTL", "—", "—", "Post hoc tissue extension (thyroid, cultured fibroblasts, subcutaneous and visceral adipose, whole blood; n = 469–670; Table S4)"],
        ["1000 Genomes Phase 3 European reference (n = 503)", "—", "—", "LD clumping of instruments and reference allele frequencies"],
    ]
    note = ("The three FinnGen R12 endpoints share participants and are not independent replications. " + CONTROLS +
            " UKB hyperthyroidism is broader than GD. 'Strict' is the FinnGen endpoint label, not proof of perfect "
            "specificity. The descriptive orbital RNA-sequencing series (four TED patients, one control) is separate from "
            "these datasets. eQTL, expression quantitative trait locus; LD, linkage disequilibrium; TED, thyroid eye disease; "
            + ABBR_DATA + ".")
    return "## Table 1. Data sources and their roles\n\n" + table(["Dataset / endpoint", "Cases", "Controls", "Role"], rows, [28, 8, 9, 48]) + "\n\n" + note


def table2(d):
    rows = []
    for g in ["TSHR", "IGF1R", "CTLA4"]:
        for o in ORDER:
            m, c, c6 = mrrow(d, g, o), coloc(d, g, o), coloc(d, g, o, 1e-6)
            mark = ("\\*" if m.pvalue < 0.05 else "") + ("†" if m.pvalue < SCREEN else "")
            rows.append([f"*{g}*", LABEL[o], int(m.n_iv), orci(m.OR, m.CI_lower, m.CI_upper), fp(m.pvalue) + mark,
                         qp(d, g, o), pp(c["PP.H2"]), pp(c["PP.H3"]), pp(c["PP.H4"], True), pp(c6["PP.H4"], True)])
    f = d["allele"]
    pa = f[(f.snp == "rs179252")].aligned_G_OR
    note = ("MR estimates are ORs per approximately one SD higher genetically predicted blood expression: Wald ratio for one "
            "IV, multiplicative random-effects IVW for several. *TSHR* used one instrument, rs179252, in every dataset; its G "
            f"allele raises blood expression by 0.105 SD per allele (F = 176), and its per-allele disease ORs were "
            f"{pa.min():.2f}–{pa.max():.2f} (Table S3). The per-SD ORs therefore rescale a moderate per-allele association, and "
            "each single-IV *P* value equals the variant's disease-association z-test *P* value. *CTLA4* multi-IV estimates "
            "are concentrated in rs13030124 and heterogeneous (Q-test *P*); in BBJ and UKB, two of four single-instrument "
            "omissions removed nominal *IGF1R* significance (Figure S2). \\* Nominal *P* < 0.05; † *P* < 0.05/2,544, the "
            "screen threshold, applied descriptively outside BBJ. PP.H2–H4 are colocalization posterior probabilities "
            "(H2, expression signal only; H3, distinct variants; H4, one shared variant) with p1 = p2 = 10^−4^ and "
            "p12 = 10^−5^ unless stated; bold, ≥0.80. They are not *P* values. BBJ colocalization pairs European eQTL "
            "with East Asian disease statistics. " + CONTROLS + " CI, confidence interval; IV, instrumental variable; IVW, "
            "inverse-variance weighted; MR, Mendelian randomization; NA, not applicable (one IV); OR, odds ratio; PP, "
            "posterior probability; SD, standard deviation; " + ABBR_DATA + ".")
    cols = ["Gene", "Dataset", "IVs", "OR (95% CI)", "*P* value", "Q-test *P*", "PP.H2", "PP.H3", "PP.H4",
            "PP.H4, p12 = 10^−6^"]
    return ("## Table 2. MR and colocalization results for *TSHR*, *IGF1R* and *CTLA4* across five disease datasets\n\n"
            + table(cols, rows, [7, 17, 4, 15, 13, 9, 7, 7, 7, 11]) + "\n\n" + note)


def table3(d, hits):
    rows = []
    for g in hits:
        b, f = mrrow(d, g, "BBJ_Graves"), mrrow(d, g, "E4_GRAVES_STRICT")
        dig = 3 if g == "MAPKAPK5" else 2
        est = bool(f.estimable) and not pd.isna(f.OR)
        focused = "NA" if not est else ("Yes" if (f.pvalue < 0.05 / 13 and np.sign(np.log(f.OR)) == np.sign(np.log(b.OR))) else "No")
        rows.append([f"*{g}*", chrom(d, g), f"{int(b.n_iv)}/{int(f.n_iv) if est else 0}", orci(b.OR, b.CI_lower, b.CI_upper, dig),
                     fp(b.pvalue), orci(f.OR, f.CI_lower, f.CI_upper, dig) if est else "NA", fp(f.pvalue) if est else "NA",
                     focused, pp(coloc(d, g, "E4_GRAVES_STRICT")["PP.H4"], True)])
    yes = sum(r[7] == "Yes" for r in rows)
    note = (f"Genes passed the BBJ screen threshold (*P* < 0.05/2,544) and are listed by BBJ *P* value. Focused support in "
            f"FinnGen GD requires the same direction and *P* < 0.05/13; {yes} of 13 qualify. IVs are counted separately for "
            "each dataset because instrument availability differs. The five MHC genes and the three chromosome 16p11.2 genes "
            "are not independent signals; the latter are compatible with one regional association. *MAPKAPK5* is shown to "
            "three decimals to preserve its small effect. NA, not estimable (no retained instrument), not a negative result; "
            "*VKORC1* colocalization does not require an MR instrument. PP.H4, shared-variant posterior probability under "
            "p12 = 10^−5^; bold, ≥0.80. Chr, chromosome of the strongest instrument; MHC, major histocompatibility complex; "
            "CI, confidence interval; IV, instrumental variable; OR, odds ratio; " + ABBR_DATA + ".")
    cols = ["BBJ discovery gene", "Chr", "IVs, BBJ/FinnGen GD", "BBJ OR (95% CI)", "BBJ *P*", "FinnGen GD OR (95% CI)",
            "FinnGen GD *P*", "Focused support", "FinnGen GD PP.H4"]
    return "## Table 3. Thirteen BBJ discoveries and their FinnGen GD comparison\n\n" + table(cols, rows, [10, 12, 11, 15, 10, 16, 11, 8, 9]) + "\n\n" + note


def tableS1(d, hits):
    mr = d["MR_primary"]
    gd = mr[(mr.outcome == "E4_GRAVES_STRICT") & (mr.pvalue < SCREEN)].sort_values("pvalue")
    extra = [g for g in gd.gene_symbol if g not in hits]
    assert extra == ["APOM", "NEU1", "IGLV2-18", "IGLV3-12", "TNF", "HIST1H4E", "LTA", "LTB", "C4B"], extra
    rows = []
    for g in extra:
        f = mrrow(d, g, "E4_GRAVES_STRICT")
        rows.append([f"*{g}*", chrom(d, g), int(f.n_iv), orci(f.OR, f.CI_lower, f.CI_upper), fp(f.pvalue),
                     pp(coloc(d, g, "BBJ_Graves")["PP.H4"], True), pp(coloc(d, g, "E4_GRAVES_STRICT")["PP.H4"], True),
                     pp(coloc(d, g, "E4_GRAVES_OPHT_STRICT")["PP.H4"], True)])
    note = ("Genes outside the BBJ discovery set that passed *P* < 0.05/2,544 in the post hoc full FinnGen GD screen (1000 "
            "Genomes European reference frequencies). These are not independent of the screen and seven of nine are in the "
            "MHC. *APOM*/*NEU1* share rs3101018 and *IGLV2-18*/*IGLV3-12* share rs5759410. The lead *IGLV* expression "
            "variant is absent in BBJ, weakening interpretation of its low PP.H4. None reached PP.H4 ≥0.80 (p12 = 10^−5^) "
            "in any dataset. Chr, chromosome; IV, instrumental variable; MHC, major histocompatibility complex; OR, odds "
            "ratio; PP.H4, shared-variant posterior probability; " + ABBR_DATA + ".")
    cols = ["Gene", "Chr", "FinnGen GD IVs", "FinnGen GD OR (95% CI)", "FinnGen GD *P*", "BBJ PP.H4", "FinnGen GD PP.H4",
            "FinnGen GO strict PP.H4"]
    return "## Table S1. Additional FinnGen GD screen discoveries\n\n" + table(cols, rows, [10, 8, 8, 15, 12, 9, 10, 11]) + "\n\n" + note


def tableS2(d):
    fdg = d["Frequency_diagnostics"]
    lab = {"original_reference": "1000G EUR", "reharmonized_eqtlgen": "eQTLGen cohort"}
    rows = []
    for g in ["C4B", "HLA-G"]:
        for o in ["E4_GRAVES_STRICT", "E4_GRAVES_OPHT_STRICT"]:
            for s in ["original_reference", "reharmonized_eqtlgen"]:
                r = fdg[(fdg.gene == g) & (fdg.outcome == o) & (fdg.scenario == s) & (fdg.removed_snp == "none")].iloc[0]
                rows.append([f"*{g}*", LABEL[o], lab[s], int(r.n_iv), orci(r.OR, r.CI_lower, r.CI_upper), fp(r.pvalue),
                             fp(r.Q_pvalue)])
    note = ("Exposure effects were reconstructed with 1000 Genomes European reference allele frequencies (primary) or "
            "eQTLGen cohort frequencies; variant sets and outcome effects were held fixed, so changes reflect rescaled "
            "expression effects and IVW weights. *HLA-G* was added because its eQTLGen-frequency *P* crossed the screen "
            "threshold (*P* < 0.05/2,544). Q tests heterogeneity, not directional pleiotropy; its *P* value is not a "
            "gene-association *P* value. All 72 omission estimates and instrument weights are in the workbook. CI, confidence "
            "interval; IV, instrumental variable; IVW, inverse-variance weighted; OR, odds ratio; " + ABBR_DATA + ".")
    cols = ["Gene", "Dataset", "Frequency source", "IVs", "OR (95% CI)", "*P* value", "Q-test *P*"]
    return "## Table S2. Frequency-sensitive exploratory candidates\n\n" + table(cols, rows, [8, 15, 13, 5, 15, 12, 12]) + "\n\n" + note


def tableS3(d):
    rc, w, al = d["Region_coverage"], d["TSHR_windows"], d["allele"]
    rcname = {"BBJ_Graves": "BBJ", "UKB_hyperthyroid": "UKB"}
    aname = {"BBJ_Graves": "BBJ", "UKB_hyperthyroid": "UKB"}
    rows = []
    for o in ORDER:
        r = rc[(rc.gene == "TSHR") & (rc.outcome == rcname.get(o, o)) & (rc.scenario == "original_reference")].iloc[0]
        a = al[(al.snp == "rs179252") & (al.outcome == aname.get(o, o))].iloc[0]
        c5, c6 = coloc(d, "TSHR", o), coloc(d, "TSHR", o, 1e-6)
        ww = w[(w.gene == "TSHR") & (w.outcome == rcname.get(o, o)) & (w.scenario == "original_reference") & np.isclose(w.p12, 1e-5)]["PP.H4"]
        assert len(ww) == 4, (o, len(ww))
        rows.append([LABEL[o], f"{int(r.n_overlap):,}", f"{r.matched_outcome_lead} ({fp(r.matched_outcome_min_p)})",
                     fp(a.recomputed_normal_p), orci(a.aligned_G_OR, a.aligned_G_CI_low, a.aligned_G_CI_high),
                     c5.top_snp, f"{pp(c5['PP.H4'], True)} / {pp(c6['PP.H4'], True)}", f"{ww.min():.3f}–{ww.max():.3f}"])
    note = ("Disease lead, the strongest disease-associated variant among the matched cis variants, with its *P* value; "
            "rs179252 *P*, z-test *P* of the aligned G-allele disease association (the UKB source BOLT-LMM *P* was "
            "1.3×10^−27^). H4 top variant, the variant with the highest posterior under the shared-variant hypothesis; it "
            "is not necessarily the disease lead. PP.H4 under p12 = 10^−5^ and 10^−6^ (p1 = p2 = 10^−4^); bold, ≥0.80. "
            "Window range covers the full source cis-region, ±500 kb, ±250 kb and the five-outcome intersection at "
            "p12 = 10^−5^; it is not a confidence interval. All 240 window settings are in TSHR_windows. BBJ pairs "
            "European eQTL with East Asian disease statistics. CI, confidence interval; OR, odds ratio; PP.H4, "
            "shared-variant posterior probability; " + ABBR_DATA + ".")
    cols = ["Dataset", "Matched variants", "Disease lead (*P*)", "rs179252 *P*", "rs179252-G OR (95% CI)",
            "H4 top variant", "PP.H4, 10^−5^ / 10^−6^", "PP.H4 range across windows"]
    return "## Table S3. *TSHR* regional colocalization and prior sensitivity\n\n" + table(cols, rows, [14, 9, 20, 10, 17, 11, 12, 12]) + "\n\n" + note


def tableS4():
    avail = pd.read_csv(TIS / "Tissue_availability_10_combinations.csv")
    est = pd.read_csv(TIS / "Tissue_MR_all_estimable.csv")
    tname = {"Thyroid": "Thyroid", "Cells_Cultured_fibroblasts": "Cultured fibroblasts", "Adipose_Subcutaneous":
             "Subcutaneous adipose", "Adipose_Visceral_Omentum": "Visceral adipose", "Whole_Blood": "Whole blood"}
    rows_a = []
    for _, r in avail.iterrows():
        e = est[(est.gene == r.gene) & (est.tissue == r.tissue)]
        if len(e):
            inst = f"{e.snps.iloc[0]} (F = {e.F_min.iloc[0]:.1f})"
            status = "Estimable"
        elif r.availability_status == "not_tested_in_GTEx_tissue":
            inst, status = "—", "Not tested by GTEx in this tissue"
        else:
            assert r.availability_status == "no_eQTL_P_below_5e-8", r.availability_status
            inst, status = "—", f"No instrument (minimum eQTL *P* = {fp(r.minimum_eQTL_P)})"
        rows_a.append([f"*{r.gene}* / {tname[r.tissue]}", int(r.n_exposure), inst, status])
    rows_b = []
    tord = {k: i for i, k in enumerate(tname)}
    est = est.assign(_g=est.gene.map({"TSHR": 0, "IGF1R": 1}), _t=est.tissue.map(tord), _o=est.outcome.map({o: i for i, o in enumerate(ORDER)}))
    for _, r in est.sort_values(["_g", "_t", "_o"]).iterrows():
        thr = 0.005 if r.family == "primary" else 0.00125
        rows_b.append([f"*{r.gene}* / {tname[r.tissue]}", LABEL[r.outcome], "Primary" if r.family == "primary" else "Secondary",
                       orci(r.OR, r.CI_lower, r.CI_upper), fp(r.P), "Yes" if r.P < thr else "No",
                       "Same" if r.direction_vs_blood == "concordant" else "Opposite"])
    note = ("Post hoc fixed-panel extension using GTEx v8; tissue colocalization was not performed. A, availability of the "
            "ten gene–tissue combinations; F statistics are for the single retained instrument. 'Not tested' means the "
            "gene is not among the genes GTEx tested in that tissue. GTEx whole blood (n = 670) yielded no instrument at "
            "the fixed threshold, whereas the main analysis used eQTLGen (up to 31,684 participants). B, all 20 estimates. "
            "The primary family is the ten FinnGen GD comparisons (*P* < 0.005); the secondary family is the 40 comparisons "
            "with the other datasets (*P* < 0.00125); unestimable combinations stay in both denominators. ORs are per unit "
            "of inverse-normal-transformed tissue expression and are not comparable in size with the blood ORs. With one "
            "instrument, each *P* value equals that variant's disease-association *P* value. The adipose *TSHR* instruments "
            "(rs12101255, rs74064796; 50 bp apart in *TSHR* intron 1) are themselves blood *TSHR* eQTLs near rs179252 and "
            "re-express the regional association rather than providing independent evidence. Cultured fibroblasts and adipose "
            "tissue are not orbital tissue. Direction, compared with the blood estimate in the same dataset. CI, confidence "
            "interval; GTEx, Genotype-Tissue Expression; OR, odds ratio; " + ABBR_DATA + ".")
    return ("## Table S4. *TSHR* and *IGF1R* tissue expression extension\n\n### A. Availability\n\n"
            + table(["Gene / tissue", "GTEx N", "Instrument (F)", "Status"], rows_a, [25, 8, 20, 40])
            + "\n\n### B. Estimates\n\n"
            + table(["Gene / tissue", "Dataset", "Family", "OR (95% CI)", "*P* value", "Below family threshold", "Direction vs blood"], rows_b,
                    [25, 16, 9, 15, 11, 10, 10])
            + "\n\n" + note)


def main():
    d = load()
    mr = d["MR_primary"]
    hits = mr[(mr.outcome == "BBJ_Graves") & mr.P_lt_0_05_over_2544].sort_values("pvalue").gene_symbol.tolist()
    assert len(hits) == 13
    main_md = "# Main tables\n\n" + "\n\n".join([table1(), table2(d), table3(d, hits)]) + "\n"
    supp_md = "\n\n".join([tableS1(d, hits), tableS2(d), tableS3(d), tableS4()]) + "\n"
    (PKG / "Tables.md").write_bytes(main_md.replace("\n", "\r\n").encode("utf-8"))
    (PKG / "SupplementaryTables.md").write_bytes(supp_md.replace("\n", "\r\n").encode("utf-8"))
    print("Tables.md and SupplementaryTables.md written")


if __name__ == "__main__":
    main()
