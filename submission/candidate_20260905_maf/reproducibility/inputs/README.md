# Local-only reproduction inputs

The repository's data rules allow scripts and aggregate analysis results, but prohibit redistribution of source datasets. The five extracted eQTL/GWAS/reference input files used for independent verification are therefore retained locally and are not committed. Their filenames, sizes and SHA-256 hashes are recorded in `../inputs_manifest.json`.

For reproduction, place authorized local copies here. Obtain source data through the original eQTLGen, GWAS Catalog (BBJ GCST90018627; UKB GCST90038636), FinnGen R12 and 1000 Genomes channels under their applicable terms. The verification report records the prior successful local calculation; a fresh checkout alone does not contain those inputs. `audit_submission.py` can run using the committed aggregate results without these files.
