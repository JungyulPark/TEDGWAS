> **13 September 2026 review status:** Tables 1–3/S1–S4 and Data 1–5 are complete for presentation; 33 references, 1,250 content checks / 350 displayed values, 6,864 figure checks and 44 reviewed Word pages. New TSHR reference-panel LD was calculated. Multi-signal colocalization, formal directionality, cohort-overlap verification and author/journal checks remain unresolved. This supersedes earlier statements that only author/journal checks remain. Read the current REVIEW_REPORT_KO.md.

# TEDGWAS — TED–TRAP

Druggable-gene-wide Mendelian randomization and colocalization comparing genetically proxied TSHR and IGF1R expression across Graves disease, hyperthyroidism and a TED-enriched outcome.

**Current manuscript: [13 September 2026 final review candidate](submission/candidate_20260911_review/README.md).** Prepared for review toward Endocrine Connections submission; multi-signal colocalization, directionality, cohort-overlap verification and final author/journal checks remain pending. The GitHub repository is public. Data-management rules are in `CLAUDE.md`.

The clinical revision contains 3,326 main-text words, 231 abstract words, 33 references, three main tables, three main figures, four supplementary tables and one supplementary figure. ORs, 95% CIs and P values lead the presentation. All 44 Word pages were visually reviewed, including seven separately supplied tables. The content audit passes 1,250 checks covering 350 displayed numerical values; the release audit also verifies current file identities and review records.

## Current evidence

- TSHR: PP.H4 0.951 in BBJ and 0.986 in FinnGen; UKB favours distinct variants under the single-causal-variant model. Multiple signals remain an untested alternative.
- IGF1R: nominal MR associations in BBJ/UKB, with unresolved shared-variant evidence. UKB splits support between H2 (0.400) and H4 (0.404). This does not establish an exclusive effector role or exclude a genetic contribution.
- CTLA4: BBJ PP.H3 0.799 / PP.H4 0.201; colocalization is supported in the two European outcomes. It does not pass the combined BBJ-plus-TED criterion.
- 2,544 eligible genes, with estimable primary MR results for 2,234 BBJ, 2,505 UKB and 2,480 FinnGen genes. No novel candidate passed the combined colocalization filter.
- Orbital expression observations from four TED samples and one control are descriptive. eQTLGen-specific frequency sensitivity was completed: discovery hits and colocalization classifications were retained, while some MR estimates and detection thresholds changed.

## Files and verification

| Path | Purpose |
|---|---|
| `MANUSCRIPT_TED_TRAP_v5_MASTER.md` | Sole editable manuscript master |
| `submission/candidate_20260911_review/` | Current Word documents, figures, supplementary data and revision report |
| `submission/candidate_20260911_review/provenance/` | Source hashes, numerical checks, independent-coloc and visual-review records |
| `submission/candidate_20260911_review/reproducibility/` | Audit, figure/document builders and reproduction instructions |
| `TrackA_MR/v5_upgrade/` | Historical analysis outputs used as evidence |
| `internal/` | Internal records, not journal upload files |

The current figure audit passes 6,864 checks across all 2,234 plotted genes and 63 posterior comparisons. Results now contains 982 words; Discussion contains 1,235 words in seven interpretation paragraphs, two limitations paragraphs and one conclusion. All 7,219 baseline primary MR estimates and 81 coloc rows were reproduced. Cohort-frequency sensitivity reran primary MR and harmonization, compared all 81 coloc settings and recalculated whole-screen power; independent primary-MR formula checks passed. These checks do not establish clinical causality or certify author declarations.

Run the current audit from the repository root after installing the packages in the candidate's `reproducibility/requirements.txt`:

```powershell
python scripts/audit_paper1_integrity.py
python submission/candidate_20260911_review/reproducibility/audit_figures.py
python submission/candidate_20260911_review/reproducibility/audit_release.py
```

See the candidate reproduction README before rebuilding. Source-data extracts are local-only; the repository contains scripts and aggregate results. IRB raw RNA-seq, patent materials and restricted source datasets must not be committed. Earlier submission binaries, scientific summaries and build/audit instructions elsewhere are historical and may disagree with this reviewed candidate.
