#!/usr/bin/env python3
"""Correct the package's reviewer-facing Markdown sources (2026-10-14 review).

Removes internal version-history and author to-do notes from the supplement and
legends, restores subscripts lost in rendering, extends the detection benchmark to
FinnGen GD (values recomputed from Supplementary_Dataset.xlsx MR_primary with the
same formula as Detection_original), restores PP.H2 to Table 2, italicises gene
names and harmonises P-value notation in Table S3. Each replacement must match
exactly once. Applied once (2026-10-14); later edits (Table 2 PP.H2 column,
table italics, P-value precision, cover letter, reproduction paragraph) are
recorded in git history.
"""
from pathlib import Path

PKG = Path(__file__).resolve().parents[1] / "submission" / "candidate_20261011_review"

SUPP = [
    ("The expanded plan was frozen before those new analyses, but was not a prospective registration of the complete study.",
     "The expanded analysis plan was recorded before the new endpoint associations were computed, but it was not a prospective registration of the complete study."),
    ("For a single instrument, the Wald coefficient is βY/βX and its first-order SE is SEY/|βX|.",
     "For a single instrument, the Wald coefficient is β~Y~/β~X~ and its first-order SE is SE~Y~/|β~X~|."),
    ("Previous TSHR clumping sensitivity used eQTL thresholds",
     "TSHR clumping sensitivity used eQTL thresholds"),
    ("Previous chromosome 16p11.2 analyses found high East Asian LD",
     "Chromosome 16p11.2 analyses found high East Asian LD"),
    (" These earlier clumping and conditional analyses were preserved; the separate allele-direction audit below newly examined rs179252 and rs179247.",
     ""),
    ("Reconstructing βX and SEX with fixed Z and N preserves their ratio up to floating-point precision.",
     "Reconstructing β~X~ and SE~X~ with fixed Z and N preserves their ratio up to floating-point precision."),
    ("The earlier detection benchmarks use (z₁₋α/₂ + z₀.₈₀) × SE as the approximate 80% minimum detectable absolute log-OR.",
     "The detection benchmarks use (z₁₋α/₂ + z₀.₈₀) × SE as the approximate 80% minimum detectable absolute log-OR."),
    ("The Detection_original worksheet retains this original three-outcome scope; no unperformed new-endpoint power analysis is implied.",
     "The Detection_original worksheet retains these three-outcome results. Applying the same formula to the MR_primary standard errors for FinnGen GD at alpha 0.05/2,544 gave 24.8% and 52.2% of 2,480 estimable genes with 80% power for ORs of 1.5 and 2.0. At alpha 0.05, the detectable *IGF1R* ORs were 1.72 in BBJ, 1.39 in UKB and 1.36 in FinnGen GD."),
    ("; independent R calculations verified all 20 estimates.",
     "; independent calculations reproduced all 20 estimates."),
    ("Figure S1 retains the previously reviewed sample-level transcript abundance and log2[(mean TED TPM + 0.01)/(control TPM + 0.01)].",
     "Figure S1 shows sample-level transcript abundance and log2[(mean TED TPM + 0.01)/(control TPM + 0.01)]."),
    (" Collection and laboratory metadata should be confirmed by the tissue investigators before submission; this reanalysis did not independently reconstruct tissue collection or sequencing procedures.",
     ""),
    ("Readme and Column_dictionary identify each sheet, its scope and its actual columns; no outcome_min_p field is described unless present in the source. Source hashes and numerical checks accompany the review package.",
     "Readme and Column_dictionary identify each sheet, its scope and its columns."),
]

LEGENDS = [
    ("Previously reviewed sample-level transcript abundance for TSHR, IGF1R and CTLA4 in four TED patients and one control.",
     "Sample-level transcript abundance for *TSHR*, *IGF1R* and *CTLA4* in four TED patients and one control."),
    (" The figure is unchanged from the earlier review package.", ""),
]

S3 = [
    ("| 0.3591 |", "| 0.359 |"),
    ("| 1.28e-34 | 1.28e-33 |", "| 1.28×10⁻³⁴ | 1.28×10⁻³³ |"),
    ("| 4.09e-19 | 4.09e-18 |", "| 4.09×10⁻¹⁹ | 4.09×10⁻¹⁸ |"),
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
    apply(PKG / "Supplementary_Material.md", SUPP)
    apply(PKG / "Figure_Legends.md", LEGENDS)
    apply(PKG / "TableS3_Tissue_extension.md", S3)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
