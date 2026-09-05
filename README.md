# TEDGWAS — TED–TRAP

Druggable-gene-wide Mendelian randomization and colocalization comparing genetically proxied TSHR and IGF1R expression across Graves disease, hyperthyroidism and a TED-enriched outcome.

**Current manuscript: [5 September 2026 submission candidate](submission/candidate_20260905/README.md).** Prepared for Endocrine Connections; final author information, declarations and journal checks remain pending. The GitHub repository is public. Data-management rules are in `CLAUDE.md`.

## Current evidence

- TSHR: PP.H4 0.951 in BBJ and 0.986 in FinnGen; UKB instead favours distinct variants.
- IGF1R: nominal MR associations in BBJ/UKB, with unresolved shared-variant evidence. UKB splits support between H2 (0.400) and H4 (0.404). This does not establish an exclusive effector role or exclude a genetic contribution.
- CTLA4: BBJ PP.H3 0.799 / PP.H4 0.201; colocalization is supported in the two European outcomes. It does not pass the combined BBJ-plus-TED criterion.
- 2,544 eligible genes, with estimable primary MR results for 2,234 BBJ, 2,505 UKB and 2,480 FinnGen genes. No novel candidate passed the combined colocalization filter.
- Orbital expression observations from four TED samples and one control are descriptive. eQTLGen-specific frequency sensitivity is unperformed and remains a limitation.

## Files and verification

| Path | Purpose |
|---|---|
| `MANUSCRIPT_TED_TRAP_v5_MASTER.md` | Sole editable manuscript master |
| `submission/candidate_20260905/` | Current Word documents, figures, supplementary data and revision report |
| `submission/candidate_20260905/provenance/` | Source hashes, numerical checks, independent-coloc and visual-review records |
| `submission/candidate_20260905/reproducibility/` | Audit, figure/document builders and reproduction instructions |
| `TrackA_MR/v5_upgrade/` | Historical analysis outputs used as evidence |
| `internal/` | Internal records, not journal upload files |

The current numerical audit compares 463 displayed values and passes 580 checks. All 41 rendered Word pages were visually reviewed. Independent local recalculation reproduced all 81 coloc rows to floating-point precision. These are provenance and layout checks, not proof of clinical causality or a rerun of the full MR pipeline.

Run the current audit from the repository root after installing the packages in the candidate's `reproducibility/requirements.txt`:

```powershell
python submission/candidate_20260905/reproducibility/audit_submission.py
```

See the candidate reproduction README before rebuilding. Source-data extracts are local-only; the repository contains scripts and aggregate results. IRB raw RNA-seq, patent materials and restricted source datasets must not be committed. Earlier submission binaries, scientific summaries and build/audit instructions elsewhere are historical and may disagree with this reviewed candidate.
