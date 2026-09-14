> **15 September 2026 compact supplementary review:** Leave-one-out is now reported in Methods/Results/Discussion, Figure S2, Data 6 and STROBE item 13. Fifteen omissions comprise 11 IGF1R and 4 CTLA4 results; all nine full-set estimates and every omission were independently checked in base R. The European CTLA4 dependence on rs13030124 is disclosed in Table 2. No new multi-signal colocalization was undertaken: its absence and the single-causal-variant assumption remain explicit limitations. Formal directionality and participant overlap remain unverified. Read the current REVIEW_REPORT_KO.md for review status and audit counts.

# TEDGWAS — TED–TRAP

Druggable-gene-wide Mendelian randomization and colocalization comparing genetically proxied TSHR and IGF1R expression across Graves disease, hyperthyroidism and a TED-enriched outcome.

**Current manuscript: [14 September 2026 leave-one-out integrated candidate](submission/candidate_20260911_review/README.md).** Prepared for final author review toward Endocrine Connections submission. The multi-signal model is explicitly outside this revision; directionality and participant overlap remain unverified, with final author/journal checks pending. The GitHub repository is public. Data-management rules are in `CLAUDE.md`.

The clinical revision contains 3,546 main-text words, 231 abstract words, 33 references, three main tables, three main figures, three supplementary tables and two supplementary figures. ORs, 95% CIs and P values lead the presentation. The package includes six separately supplied tables. Current audit and visual-review totals are recorded in the candidate review report; the release audit also verifies current file identities and review records.

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

The current figure audit passes 7,032 checks across all 2,234 plotted genes and 63 posterior comparisons, plus all 20 leave-one-out plot rows. Results now contains 1,100 words; Discussion contains 1,292 words in seven interpretation paragraphs, two limitations paragraphs and one conclusion. All 7,219 baseline primary MR estimates and 81 coloc rows were reproduced. Cohort-frequency sensitivity reran primary MR and harmonization, compared all 81 coloc settings and recalculated whole-screen power; independent primary-MR formula checks passed. These checks do not establish clinical causality or certify author declarations.

Run the current audit from the repository root after installing the packages in the candidate's `reproducibility/requirements.txt`:

```powershell
python scripts/audit_paper1_integrity.py
python submission/candidate_20260911_review/reproducibility/audit_figures.py
python submission/candidate_20260911_review/reproducibility/audit_release.py
```

See the candidate reproduction README before rebuilding. Source-data extracts are local-only; the repository contains scripts and aggregate results. IRB raw RNA-seq, patent materials and restricted source datasets must not be committed. Earlier submission binaries, scientific summaries and build/audit instructions elsewhere are historical and may disagree with this reviewed candidate.
