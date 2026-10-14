#!/usr/bin/env python3
"""Apply the fixes confirmed by the independent verification of the 2026-10-14 revision.

Findings (scratch verification run wf_7c26b909): orbital Abstract sentence lost
"Mean" (false for one IGF1R sample); prior sensitivity flagged for BBJ only although
strict GO falls further; "two largest datasets" true only by cases; thyroid IGF1R
direction generalised; wrong workbook sheet cited; ambiguous antecedents and
undefined "GO" in the main text; reproduction paragraph counts unexplained; roman
gene names/P and flat subscripts in the supplement; cover letter overstated the
STROBE crosswalk and restated unconfirmed declarations. Each replacement must
match exactly once. Applied once; later edits are recorded in git history.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "submission" / "candidate_20261011_review"

MASTER = [
    ("Jungyul Park¹, Min-Seon Kim², Kyung-Hwa Shin³⁎, Suk-Woo Yang¹⁎",
     "Jungyul Park¹, Min-Seon Kim², Kyung-Hwa Shin³\\*, Suk-Woo Yang¹\\*"),
    # Abstract (word-neutral overall)
    ("**Objective:** Graves disease (GD), an autoimmune cause of hyperthyroidism, can be accompanied by thyroid eye disease (TED).",
     "**Objective:** Thyroid eye disease (TED) can accompany Graves disease (GD), an autoimmune cause of hyperthyroidism."),
    ("UK Biobank (UKB) hyperthyroidism and three overlapping FinnGen endpoints: GD, Graves ophthalmopathy (GO) and strict GO.",
     "UK Biobank (UKB) hyperthyroidism and overlapping FinnGen GD, Graves ophthalmopathy (GO) and strict GO endpoints."),
    ("Tissue eQTL and orbital RNA-sequencing (four TED patients, one control) analyses were supplementary.",
     "Tissue eQTL and orbital RNA-sequencing (four TED, one control) were supplementary."),
    ("Shared-variant probability exceeded 0.80 in BBJ (prior-sensitive) and both GO endpoints, but UKB (0.226) and FinnGen GD (0.264) favoured distinct variants.",
     "Shared-variant probability exceeded 0.80 in FinnGen GO and, prior-sensitively, BBJ and strict GO, but UKB (0.226) and FinnGen GD (0.264) favored distinct variants."),
    ("the most precise, FinnGen GD, was 1.07 (0.86–1.33).",
     "the most precise estimate, FinnGen GD, was OR 1.07 (95% CI 0.86–1.33)."),
    ("Orbital *TSHR* and *IGF1R* transcript levels exceeded the single control.",
     "Mean orbital *TSHR* and *IGF1R* transcript levels exceeded the control."),
    ("**Conclusions:** Association evidence was stronger and more consistent for *TSHR* than *IGF1R*, but absent colocalization in the two largest datasets leaves expression-mediated causality unproven. Population-controlled data cannot isolate TED susceptibility within GD; weaker *IGF1R* evidence does not bear on IGF-1R–targeted treatment.",
     "**Conclusions:** Association evidence was stronger for *TSHR* than *IGF1R*, but absent colocalization in UKB and FinnGen GD leaves expression-mediated causality unproven. Population-controlled data cannot isolate TED susceptibility within GD; weaker *IGF1R* evidence does not argue against IGF-1R–targeted treatment."),
    # Methods
    ("We subsequently assessed both genes in a fixed panel of five GTEx v8 tissues",
     "In a post hoc extension, we assessed both genes in a fixed panel of five GTEx v8 tissues"),
    ("An in-house orbital RNA-sequencing series of four TED patients and one control provided descriptive context only.",
     "An in-house orbital RNA-sequencing series of four patients with inactive TED and one control provided descriptive context only (Supplementary Methods)."),
    # Results
    ("(*P* = 0.66, 0.56 and 0.43 in BBJ, UKB and FinnGen GO)",
     "(*P* = 0.66, 0.56 and 0.43 in BBJ, UKB and FinnGen ophthalmopathy)"),
    ("Their cross-dataset associations should therefore be distinguished from agreement among instruments.",
     "For the first three, cross-dataset associations should therefore be distinguished from agreement among instruments."),
    ("removed nominal *IGF1R* significance: rs2654980 (*P* = 0.0751 and 0.468) and, respectively, rs59467480 (*P* = 0.0503) and rs117212126 (*P* = 0.0529);",
     "removed nominal *IGF1R* significance: rs2654980 in both (*P* = 0.0751 and 0.468), rs59467480 in BBJ (*P* = 0.0503) and rs117212126 in UKB (*P* = 0.0529);"),
    ("(Supplementary Dataset, MR_new_freq)",
     "(Supplementary Dataset, MR_original_freq and MR_new_freq)"),
    ("the corresponding FinnGen GD proportions were 24.8% and 52.2% (Supplementary Dataset, Detection_original; Supplementary Methods).",
     "the corresponding FinnGen GD proportions were 24.8% and 52.2% (BBJ: Supplementary Dataset, Detection_original; FinnGen GD: Supplementary Methods)."),
    ("It was lower in UKB (0.226) and FinnGen GD (0.264), where PP.H3 was 0.774 and 0.736, respectively. The latter patterns favoured distinct variants",
     "It was lower in UKB (0.226) and FinnGen GD (0.264), where PP.H3 was 0.774 and 0.736, respectively. The latter patterns favored distinct variants"),
    ("and both are themselves blood *TSHR* eQTLs (eQTLGen *P* ≈ 10⁻³⁴). With one instrument, each tissue *P* value equals that variant's disease-association *P* value,",
     "and both are themselves blood *TSHR* eQTLs (eQTLGen *P* = 4.1×10⁻³⁴ and 6.5×10⁻³⁴). With one instrument, each tissue *P* value equals that variant's disease-association z-test *P* value,"),
    # Discussion
    ("but favoured distinct variants in UKB and FinnGen GD, the two datasets with the most cases.",
     "but favored distinct variants in UKB and FinnGen GD, the two datasets with the most cases; BBJ and strict ophthalmopathy support was prior-sensitive."),
    ("with that caveat, its direction-consistent associations across all five outcomes",
     "with that caveat, *CTLA4*'s direction-consistent associations across all five outcomes"),
    ("but the thymic step was not tested here, and it does not demonstrate that thymic expression",
     "but the thymic step was not tested here, and this consistency does not demonstrate that thymic expression"),
    ("The marked departure of the MR ORs from 1 reflects scaling rather than a large per-variant effect: rs179252 shifts",
     "The marked departure of the MR ORs from 1 reflects scaling: rs179252 shifts"),
    ("0.78 (FinnGen GO), 0.77 (FinnGen GD) and 0.82 (strict GO). Dividing these small effects by the small expression effect",
     "0.78 (FinnGen ophthalmopathy), 0.77 (FinnGen GD) and 0.82 (strict ophthalmopathy). Dividing these moderate per-allele effects by the small expression effect"),
    ("The loss of nominal significance after removing rs2654980 further limits confidence",
     "Loss of nominal significance after omitting single variants further limits confidence"),
    ("and the thyroid estimate ran opposite to the blood direction.",
     "and directions were inconsistent across datasets."),
    ("teprotumumab efficacy rests on randomized trials in established TED [3], and experimental TSHR–IGF-1R crosstalk [5] concerns receptor signalling",
     "teprotumumab efficacy rests on a randomized trial in active TED [3], and experimental TSHR–IGF-1R crosstalk [5] concerns receptor signaling"),
    ("whereas our primary settings used 10⁻⁵ and at least 0.80.",
     "whereas our primary settings used p12 = 10⁻⁵ and PP.H4 ≥ 0.80."),
    ("Colocalization supported a shared variant for *TSHR* in BBJ and the ophthalmopathy endpoints but not in UKB or FinnGen GD,",
     "Colocalization supported a shared variant for *TSHR* in BBJ and the ophthalmopathy endpoints, partly prior-dependently, but not in UKB or FinnGen GD,"),
]

SUPP = [
    ("## Detailed methods", "## Supplementary methods"),
    ("The 24 colocalization genes comprise 13 BBJ discoveries, IGF1R, nine additional primary FinnGen GD discoveries, and frequency-sensitive HLA-G.",
     "The 24 colocalization genes comprise 13 BBJ discoveries, *IGF1R*, nine additional primary FinnGen GD discoveries, and frequency-sensitive *HLA-G*."),
    ("Significant cis-eQTL instruments (P < 5×10⁻⁸ within", "Significant cis-eQTL instruments (*P* < 5×10⁻⁸ within"),
    ("Single-instrument TSHR and BBJ CTLA4 were ineligible.", "Single-instrument *TSHR* and BBJ *CTLA4* were ineligible."),
    ("and 72 C4B/HLA-G omission estimates", "and 72 *C4B*/*HLA-G* omission estimates"),
    ("C4B/HLA-G instrument weights use βX²/SEY², normalized", "*C4B*/*HLA-G* instrument weights use β~X~²/SE~Y~², normalized"),
    ("A single-SNP P value can be much smaller than the full IVW P value", "A single-SNP *P* value can be much smaller than the full IVW *P* value"),
    ("TSHR clumping sensitivity used eQTL thresholds", "*TSHR* clumping sensitivity used eQTL thresholds"),
    ("Original GWAS P values and recalculated MR P values were retained separately.",
     "Original GWAS *P* values and recalculated MR *P* values were retained separately."),
    ("C4B and HLA-G were audited with identical retained variants", "*C4B* and *HLA-G* were audited with identical retained variants"),
    ("The post hoc tissue extension retained a frozen panel of TSHR and IGF1R in five",
     "The post hoc tissue extension retained a frozen panel of *TSHR* and *IGF1R* in five"),
    ("Variants were restricted to ±1 Mb of the source GENCODE v26 TSS, P < 5×10⁻⁸,",
     "Variants were restricted to ±1 Mb of the source GENCODE v26 TSS, *P* < 5×10⁻⁸,"),
    ("One of 20 thyroid IGF1R candidates was absent", "One of 20 thyroid *IGF1R* candidates was absent"),
    ("Wald estimates used βY/βX and SEY/|βX| with normal two-sided P values", "Wald estimates used β~Y~/β~X~ and SE~Y~/|β~X~| with normal two-sided *P* values"),
    ("Primary adjusted P = min(1,10P) for FinnGen GD; secondary adjusted P = min(1,40P)",
     "Primary adjusted *P* = min(1, 10*P*) for FinnGen GD; secondary adjusted *P* = min(1, 40*P*)"),
    ("and no differential-expression P value was used.", "and no differential-expression *P* value was used."),
    ("Four TED patients and one control contributed orbital RNA-sequencing data. Technical replicates were combined by biological sample.",
     "Four patients with inactive TED and one control contributed orbital RNA-sequencing data. Orbital adipose tissue (approximately 2 ml, inferolateral compartment) was obtained during orbital decompression for persistent proptosis, diplopia or exposure keratopathy from patients with a documented active phase (Clinical Activity Score ≥3/7) who had completed immunosuppressive treatment, had remained inactive (Clinical Activity Score <3/7) for at least 6 months and had been euthyroid for at least 3 months. The control tissue was the nasal (medial) fat pad of one person without thyroid or autoimmune disease who underwent cosmetic upper blepharoplasty. Samples were snap-frozen; RNA was extracted with the RNeasy Lipid Tissue Mini kit (QIAGEN), and paired-end libraries (200–400 bp inserts) were sequenced on an Illumina NovaSeq platform, with two technical-replicate libraries per donor. Reads were checked with FastQC 0.11.7, trimmed with Trimmomatic 0.38, aligned to GRCh38 with HISAT2 2.1.0 and quantified as TPM with StringTie 2.1.3b. Technical replicates were combined by biological sample."),
    # reproduction paragraph (counts explained; 27 includes the control)
    ("For *TSHR*, *IGF1R* and *CTLA4*, colocalization and MR in FinnGen GD and strict ophthalmopathy were independently reproduced from re-downloaded FinnGen R12 regional statistics using two separately written implementations. Variant matching was identical (6,708/6,675, 7,894/7,850 and 5,023/4,993 rows matched/allele-compatible), and the top variant was the same in all 27 gene–endpoint–prior settings. Because the 1000 Genomes European panel was not accessible for this check, gnomAD v2.1.1 European frequencies were substituted; on the original ophthalmopathy endpoint this changed PP.H4 by at most 4.1×10⁻⁴. Reported new-endpoint PP.H4 values were reproduced within 1.2×10⁻³, with no setting crossing 0.80; *IGF1R* PP.H2 and PP.H3 differed by up to 0.003, the same shift observed for the original endpoint. MR estimates were reproduced to numerical precision on the original exposure scale.",
     "For *TSHR*, *IGF1R* and *CTLA4*, colocalization and MR in FinnGen GD and strict ophthalmopathy were independently reproduced from re-downloaded FinnGen R12 regional statistics using two separately written implementations, with the original ophthalmopathy endpoint as a control. The allele-compatible variant sets were the same (6,675, 7,850 and 4,993 variants), and the top variant was the same in all 27 gene–endpoint–prior settings, including the control. The reproduction used gnomAD v2.1.1 European frequencies, which were available for more of these variants than the 1000 Genomes European frequencies used in the reported analysis (6,674, 7,850 and 4,992 analysed, versus 6,579, 7,740 and 4,924). On the control endpoint, PP.H4 changed by at most 4.1×10⁻⁴; for *TSHR* this difference came from the additional variants. Reported new-endpoint PP.H4 values were reproduced within 1.2×10⁻³, and no setting changed side of 0.80. The *TSHR* strict-ophthalmopathy value at p12 = 10⁻⁶ was 0.449 rather than 0.451 (0.451 on the reported variant set), and *IGF1R* PP.H2 and PP.H3 differed by up to 0.003, as on the control endpoint. MR estimates were reproduced to numerical precision on the original exposure scale."),
]

TABLES = [
    ("| *CTLA4* | UKB HT | 2 | 0.21 (0.08–0.57) | **0.00220** |", "| *CTLA4* | UKB HT | 2 | 0.21 (0.08–0.57) | **0.0022** |"),
    ("IVs, instrumental variables (retained genetic variants). Wald ratio",
     "GD, Graves disease; GO, Graves ophthalmopathy; HT, hyperthyroidism/thyrotoxicosis. IVs, instrumental variables (retained genetic variants). Wald ratio"),
]

S12 = [
    ("| BBJ GD | 4336 |", "| BBJ GD | 4,336 |"),
    ("| FinnGen GD | 6579 |", "| FinnGen GD | 6,579 |"),
    ("| UKB HT | 5667 |", "| UKB HT | 5,667 |"),
    ("| FinnGen GO | 6579 |", "| FinnGen GO | 6,579 |"),
    ("| FinnGen GO strict | 6579 |", "| FinnGen GO strict | 6,579 |"),
]

LEGENDS = [
    ("H0/H1 denotes neither trait or disease only; H2 expression only;",
     "H0/H1 denotes no association or disease association only; H2 expression only;"),
]

COVER = [
    ("That success is sometimes read as evidence that *IGF1R* also carries inherited susceptibility to Graves disease (GD). Therapeutic efficacy and disease susceptibility are different claims, and we asked",
     "That success raises the question of whether *IGF1R* also carries inherited susceptibility to Graves disease (GD). Therapeutic efficacy and disease susceptibility are different claims, and we asked"),
    ("Colocalization supported a shared causal variant in Biobank Japan and both ophthalmopathy endpoints but favoured distinct variants",
     "Colocalization supported a shared causal variant in Biobank Japan and both ophthalmopathy endpoints (in Biobank Japan and the strict endpoint only under the primary prior) but favored distinct variants"),
    ("a single results workbook and a completed STROBE-MR checklist.",
     "a single results workbook and a STROBE-MR reporting crosswalk."),
    ("This manuscript has not been published and is not under consideration elsewhere. The authors report no specific funding and no competing interests. The orbital component received Pusan National University Hospital IRB approval 2104-018-102, with written informed consent.",
     "This manuscript has not been published and is not under consideration elsewhere. Funding, competing interests, ethics approval (Pusan National University Hospital IRB 2104-018-102) and written informed consent are declared in the manuscript."),
]


def apply(path: Path, edits) -> None:
    raw = path.read_bytes().decode("utf-8")
    for old, new in edits:
        n = raw.count(old)
        if n != 1:
            raise SystemExit(f"{path.name}: expected one match, found {n}: {old[:80]!r}")
        raw = raw.replace(old, new)
    path.write_bytes(raw.encode("utf-8"))
    print(f"{path.name}: {len(edits)} edits")


def main() -> int:
    apply(ROOT / "MANUSCRIPT_TED_TRAP_v5_MASTER.md", MASTER)
    apply(PKG / "Supplementary_Material.md", SUPP)
    apply(PKG / "Tables.md", TABLES)
    apply(PKG / "SupplementaryTables_S1_S2.md", S12)
    apply(PKG / "Figure_Legends.md", LEGENDS)
    apply(PKG / "Cover_Letter_Draft.md", COVER)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
