> **15 September 2026 compact supplementary review:** Leave-one-out is now reported in Methods/Results/Discussion, Figure S2, Data 6 and STROBE item 13. Fifteen omissions comprise 11 IGF1R and 4 CTLA4 results; all nine full-set estimates and every omission were independently checked in base R. The European CTLA4 dependence on rs13030124 is disclosed in Table 2. No new multi-signal colocalization was undertaken: its absence and the single-causal-variant assumption remain explicit limitations. Formal directionality and participant overlap remain unverified. Read the current REVIEW_REPORT_KO.md for review status and audit counts.

# Reproduction and verification

`audit_submission.py` compares the manuscript and Word tables with preserved original analytical outputs, checks full data consolidation and tests that essential limitations remain. `build_documents.py` reads the repository-root Markdown master and the local cover-letter/checklist Markdown files. `build_figures.py` regenerates the three main scientific figures (via `build_scientific_figures.py`) from the verified primary results. Figure S1 is preserved from the previously inspected descriptive figure; point-level tissue CSV records are not distributed.

Install the packages in requirements.txt and use Pandoc to build Word files. After a build, render and inspect every page before delivery. Builders do not make the author declarations true or confirm journal acceptance criteria.

The executed frequency-sensitivity scripts are preserved in analysis_scripts. Their historical workspace layout is described in `archives/submission_candidates/candidate_20260905_maf/reproducibility/README.md`, at the analytical baseline commit recorded in provenance/clinical_revision.json. The locally authorized raw inputs are listed with their sizes and SHA-256 hashes in `inputs_manifest.json`; see `inputs/README.md`. No statistical rerun was performed merely to simplify presentation on 11 September. Consolidated files preserve primary and cohort-frequency results, and the original verification reports and input hashes remain in provenance/maf.

Raw GWAS/eQTL extracts, cohort frequency rows, genotype counts and harmonized variant records are local-only. The AF download previously required an explicitly authorized one-file certificate exception because the server certificate had expired. The received hash is an integrity record, not a publisher authenticity checksum. This package does not include an automatic TLS-bypass downloader. The new narrative, table numbering and author-review documents supersede older package instructions.

`audit_release.py` verifies current manuscript and packaged file hashes, numerical audit summary, word counts and the identity/page counts of all ten visually reviewed DOCX files. Run it after the content and figure audits. `audit_figures.py` compares every plotted gene and posterior with the canonical aggregate sources; it does not inspect appearance. The release audit additionally checks the reviewed PNG/PDF identities and source-manifest identity. Results/Discussion expansion on 12 September did not alter analytical estimates. Visual review remains a human inspection step; this script cannot perform it. `build_documents.py` also emits separate editable Tables 1–3 and S1–S3.

To update only Figure 2, run `python submission/candidate_20260911_review/reproducibility/build_figures.py --figures Figure2` from the repository root. The generator retains the other figure files and recomputes the common source record. Refresh the Figure 2 review identity and release manifest after visual inspection.

The figure audit uses PyMuPDF to verify black exported text and the actual P-column markers/font weight. `figure_layout_checks.json` records label colours and the Figure 3 cell-padding checks performed by the builder; `figure_pdf_text_checks.json` links the PDF text checks to each exported file hash. Rebuild and rerun the audits after any figure change.

`build_complete_statistics.py` exports Data 5 with test-specific CIs and verifies every original P value. `verify_major_review_evidence.py --input-dir LOCAL_MAF_INPUTS` verifies the new scale/H2/LD report; see its docstring for the PLINK LD command. Local input rows and reference panels must not be committed. No fresh multi-signal or Steiger analysis is represented as completed.

## Reported post hoc leave-one-out analysis (14 September 2026)

`explore_leave_one_out.py` preserves the original Python calculation. `validate_leave_one_out.R` independently fits weighted regressions using base R and single-SNP Wald ratios. Run it with the locally authorized `harmonized_all.csv` and a new output path:

```text
Rscript --vanilla validate_leave_one_out.R LOCAL_HARMONIZED_ALL.csv NEW_R_VALIDATION.csv
python build_leave_one_out_outputs.py --r-validation NEW_R_VALIDATION.csv
python build_leave_one_out_figure.py
python build_documents.py MANUSCRIPT_Submission SUPPLEMENTARY_MATERIAL STROBE_MR_CHECKLIST
```

The output builder requires 9 reproduced baseline estimates and all 15 eligible omissions, checks Python/R agreement, and writes Data 6 plus aggregate verification records. Figure S2 uses every Data 6 row: 15 omissions and five full-set comparators. It has its own source/hash record, `leave_one_out_figure_sources.json`; the original figures are preserved. The content and figure audits now check these rows as well. Re-render and inspect every changed Word page and the new plot, then update the release manifest and run the release audit. Do not record visual PASS automatically from an uninspected build.

These are original-reference-frequency results only. No new multi-signal colocalization was undertaken in this revision, as agreed in the final review; its absence remains a manuscript limitation. Directionality and participant overlap are still unverified. The older `provenance/posthoc_20260914` records describe the preliminary stage and are retained as historical records; the current reported state is in `leave_one_out_verification.json` and the manuscript. Audit JSON uses explicit CRLF to preserve stable bytes across platforms, retaining the remote 10b14ef fix.
