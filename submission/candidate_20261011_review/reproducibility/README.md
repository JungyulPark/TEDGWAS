# Reproducibility and current checks

Run from the repository root: `python scripts/audit_paper1_integrity.py`.
This read-only handoff check validates tracked file hashes, DOCX/PDF/XLSX presence,
the 50-cell tissue plan, and recomputes the 20 tissue Wald ratios from the retained
harmonized aggregate instrument table. It does not reproduce unprovided regional
source data, certify causal assumptions, or resolve the pending scientific review.

`work/` contains the original executed September/October analysis and build scripts.
They retain original workspace layouts and are provenance, not a turnkey installer.
Restore them under a matching `work/` directory, supply authorized local datasets,
and adapt runtime/output paths before executing. Some historical scripts reference
earlier output directories. Do not run them blindly over the reviewed candidate.
The existing repository's archived September reproducibility package retains the
earlier blood MR, frequency, colocalization and leave-one-out calculation scripts.

Tissue scripts are also in `../Verification/Tissue_extension_scripts/`; the allele
audit is in `../Verification/Previous_allele_audit/Reproduce/`. Dataset providers,
source hashes, aggregate results and protocol freezes are retained. Original source
extracts, reference genotypes and execution logs listed in GIT_EXPORT.json are
local-only. Machine-specific home paths in provenance/code are placeholders.

Archived manifest/audit reports describe the original local delivery and may list
intentionally omitted files. The current Git checkout manifest is GIT_MANIFEST.json.
