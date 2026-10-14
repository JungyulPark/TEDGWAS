# Reproducing this audit

These are archived copies of the executed scripts. Their original working location is `<workspace>/work/allele_direction_20261003/`; copy them there before running. The workspace contains the prior canonical analysis inputs in `work/` and `outputs/`. Original large data are read from `C:/ProjectTEDGWAS`; those files are not bundled. Their SHA-256 hashes are in SOURCE_MANIFEST.json.

1. Preserve the delivered package before a rerun, since audit outputs use the same versioned directory.
2. Do not rerun freeze_protocol.py over the frozen plan. Its intentional existence guard prevents replacing the pre-analysis record.
3. Run 01_audit.py to extract the original raw records again, compute LD, and verify MR. The `--reuse-extracts` option is only for recovery using the exact already-generated raw extracts.
4. Run 02_report.py to build the interpretive report, final checks, and ZIP.

The original protocol-freeze and outcome files are provenance records. A future tissue analysis or amended protocol must use a new version and retain these files. Dependencies: Python with numpy and pandas; the configured PLINK 1.9 executable. The local run used the bundled Codex Python; no new packages were installed.
