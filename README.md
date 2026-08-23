# TEDGWAS — TED-TRAP

Druggable-gene-wide Mendelian randomization, colocalization, and orbital
transcriptomics distinguishing **TSHR-anchored genetic susceptibility** from
**IGF1R pharmacologic effector biology** in Graves disease (GD) and thyroid eye
disease (TED).

**Status:** v5 — prepared for submission (Journal of Endocrinological Investigation).
**Private repository.** See `CLAUDE.md` for the full data/scientific rule set.

## Key result
- **TSHR** — expression-colocalized susceptibility anchor (rs179252; coloc PP.H4
  0.951 BBJ / 0.986 FinnGen; consistent protective MR across 3 outcomes).
- **IGF1R** — nominal risk-direction MR without expression colocalization
  (H2-dominant) → interpreted as a pharmacologic effector axis.
- **CTLA4** — positive-control autoimmune locus.
- Screen of 2,545 MR-testable druggable genes yielded **no robust novel target**
  after cross-outcome colocalization (an informative, conservative result).

## Where things are
| Path | Contents |
|---|---|
| `submission/` | **Final JEI package** — manuscript .docx, cover letter, 5 figures, checklist |
| `MANUSCRIPT_TED_TRAP_v5_MASTER.md` | Markdown source of truth (md5 `a672b8e2…`) |
| `FIGURE_VERIFICATION.md` | Figure-vs-master data verification record |
| `TrackA_MR/v5_upgrade/` | Core v5 analysis results + canonical figures |
| `scripts/` | Numbered R analysis + public-data download scripts (no data) |
| `internal/` | Notes **not for submission** (external-GEO decision, venue cascade, security) |

## Data sources (not stored here)
| Dataset | Role | Access |
|---|---|---|
| eQTLGen blood cis-eQTL | Exposure instruments | Public (download via script) |
| Biobank Japan Graves (GCST90018627) | Primary discovery | GWAS Catalog |
| UK Biobank hyperthyroidism (GCST90038636) | Replication | GWAS Catalog |
| FinnGen R12 Graves ophthalmopathy | TED-enriched sensitivity | FinnGen |
| In-house orbital RNA-seq (4 TED + 1 ctrl) | Exploratory tissue | **IRB-restricted, on request** |

> ⚠️ Raw RNA-seq, patent material, and license-bound data are **never** committed
> (`.gitignore`). The in-house `data.txt` was purged from git history — see
> `internal/GITHUB_SECURITY_CLEANUP.md`.

## Reproduce
1. Run download scripts in `scripts/` to fetch public summary statistics locally.
2. Run MR / coloc / fine-mapping / tissue scripts (R 4.3.3; TwoSampleMR, coloc,
   susieR, DESeq2; PLINK v1.9).
3. Build figures from the `build_figure*.R` scripts; verify inputs against the
   master (`FIGURE_VERIFICATION.md`).

IRB: Pusan National University Hospital 2104-018-102.
