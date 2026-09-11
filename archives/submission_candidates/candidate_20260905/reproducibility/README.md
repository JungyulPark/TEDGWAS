# GitHub copy — input availability

The five source-data extracts are excluded from this GitHub copy under the repository data rules. Numerical audit and document/figure builds use committed aggregate results. Independent colocalization reproduction additionally needs authorized local inputs; see `inputs/README.md` and `inputs_manifest.json`. The sole master Markdown is `../../../MANUSCRIPT_TED_TRAP_v5_MASTER.md` relative to this directory.

# Reproduction and optional frequency sensitivity

The submitted data files are the source for the manuscript and rebuilt figures. Run the following from this directory in an environment containing the packages listed in requirements.txt:

```powershell
python .\audit_submission.py
python .\reproduce_colocalization.py
```

The first command checks displayed manuscript values against full-precision source values and checks the Word tables. Its pass does not certify author declarations, journal acceptance, or the unperformed allele-frequency sensitivity. The second independently reproduces 81 colocalization rows from the authorized local public-data extracts and compares posteriors, overlap counts and top SNPs with the frozen v2 CSV. It writes a new `rerun_verification` folder.

`build_figures.py` recreates Figures 1–3 and S1–S2. `build_documents.py` rebuilds the three Word files from the master Markdown and the included cover-letter text. Rebuilding overwrites these candidate outputs; run it in a copy if retaining the current package. A Word rebuild must be rendered and visually checked again before submission. Pandoc and the listed Python packages are required; `pypandoc_binary` can supply Pandoc.

The historical source used R 4.3.3, coloc 5.2.3 and TwoSampleMR 0.7.4. Independent colocalization verification used the beta/varbeta equations in coloc 5.2.3 (`R/claudia.R`), with dataset 1 = outcome and dataset 2 = eQTL. The exposure prior SD is 0.15 with sdY = 1; the case-control outcome prior SD is 0.20. UKB beta and SE are both divided by mu(1−mu), where mu = 3731/484598. No outcome allele-frequency substitution is used in these Bayes factors.

The locally preserved EUR frequency extract is restricted to SNPs in the nine-gene eQTL extract. This filtering reduces storage and leaves all verified coloc inputs unchanged. These data are public summary-statistic extracts, not individual participant records. Source identifiers and hashes are in `../provenance/source_manifest.json`.

The MAF analysis remains unperformed. The bundled `Prepare_MAF_inputs.ps1` only downloads the official eQTLGen frequency file and checks installed R packages. From PowerShell in this directory:

```powershell
& .\Prepare_MAF_inputs.ps1
```

If local script policy prevents execution, use the official [eQTLGen cis-eQTL page](https://www.eqtlgen.org/cis-eqtls.html) to download its allele-frequency file manually. HTTPS validation failed in the preparation environment, and the existing local R user-library folder could not be read. Neither security checks nor filesystem permissions were bypassed.

Before a future sensitivity analysis is called complete, it must align the frequency to the assessed allele, compare the same variants first, verify max absolute change in reconstructed Z < 1e-8, reproduce the baseline estimates, rerun harmonization, report all three IGF1R IVW estimates, compare all 81 coloc rows and 0.80 threshold crossings, and recalculate the full-screen Table S8 metrics across the verified 2,544-gene denominator. Missing variants and changed harmonization exclusions must be reported separately from frequency-only changes. Keep the manuscript's unperformed-analysis limitation until those checks have actually run. The older repository sensitivity script was not validated for this package and should not be treated as a completed analysis.

After any rebuild or edit, repeat the relevant checks and regenerate the package manifest. The supplied manifest binds only this delivered version. If the manuscript title or author details change, update the separate cover-letter Markdown as well.
