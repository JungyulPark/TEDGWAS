# Execution record

Scripts 01–06 retrieve the fixed sources, select instruments, align disease records, calculate the 50-cell grid and independently verify the 20 estimates in base R. 07 creates the review forest plot. Later scripts edit, render, validate and package the manuscript.

These are the executed workspace scripts with their original relative workspace assumptions, not a portable one-command installer. To reproduce, place them under work/tissue_extension_20261003 beneath the workspace root; restore source datasets at the recorded paths or explicitly adapt paths. Sources, protocol freezes, request metadata, hashes and full selected records are in the adjacent Tissue_extension folder. Python pandas/numpy, base R 4.3.3 and PLINK 1.9 were used; document/spreadsheet integration also used python-docx, pypdf, bundled LibreOffice and Artifact Tool. Original large GWAS files, cached GTEx compressed downloads and individual reference genotypes are not redistributed here.

The first extraction attempt with PLINK's 1 GB memory allocation failed; the preserved log identifies this resource limit. An 8 GB extraction then succeeded. It did not alter the selection criteria. Subsequent clumping and independent genotype-based selection agreed. Initial document assembly assertions were corrected before final output; final source, numeric and rendering checks are recorded in final_checks.json. R_sessionInfo.txt records the statistical verification runtime.

The workbook integration reuses the OPC merge routine from the previous allele extension after Artifact Tool authored all changed cells. The earlier helper is preserved as previous_merge_workbook.py. To run 13_prepare_merge.py unchanged, restore that helper to work/allele_integration_20261003/merge_workbook.py.

Completed analytical CSVs are immutable. To reproduce into a fresh run, change output directory names rather than overwriting this reviewed run.
