#!/bin/bash
# =============================================================================
# Task B.3 — GCTA-COJO Conditional Analysis at TSHR Locus
# =============================================================================
# Purpose : Identify independent cis-eQTL signals at TSHR locus via conditional
#           stepwise selection (GCTA-COJO).
# Output  : 02_tshr_finemap/results/TaskB_03_tshr_cojo_EUR.jma.cojo
#           02_tshr_finemap/results/TaskB_03_tshr_cojo_EAS.jma.cojo
# Spec    : WEEK1_TASK_SPEC_v1.1.md, Section 4, Step B.3
#
# REFERENCE:
#   Yang J et al. 2012 Nat Genet 44:369-375 (GCTA-COJO method paper)
#   Documentation: https://yanglab.westlake.edu.cn/software/gcta/#COJO
#
# IMPORTANT: GCTA-COJO requires:
#   - .ma file with columns: SNP A1 A2 freq b se p N
#   - PLINK bfile for LD reference (matched ancestry to outcome)
#   - Sufficient sample size in LD reference (>500 recommended)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
WORK_DIR="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
INPUT_TSV="${WORK_DIR}/02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv"
COJO_MA_FILE="${WORK_DIR}/02_tshr_finemap/data/TSHR_eqtlgen_locus.ma"
OUT_DIR="${WORK_DIR}/02_tshr_finemap/results"

# GCTA binary (Gemini to confirm path; on Windows likely gcta64.exe)
GCTA_BIN="c:/ProjectTEDGWAS/tools/gcta-1.95.1-windows-x86_64/bin/gcta64.exe"

# LD reference bfiles
BFILE_EUR="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EUR_phase3/EUR_chr14"
BFILE_EAS="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_EAS_phase3/EAS_chr14"

# COJO parameters (per spec v1.1, Yang 2012 defaults)
COJO_P=5e-6        # P threshold for stepwise selection (relaxed from 5e-8 to allow secondary signals)
COJO_WIND=10000    # window in kb
COJO_COLLINEAR=0.9 # collinearity threshold

# -----------------------------------------------------------------------------
# STEP 1 — Convert TSV to COJO .ma format
# -----------------------------------------------------------------------------
# .ma format requires beta and SE. eQTLGen provides only Z + N.
# We delegate this conversion to an R helper script that derives beta/se using
# 1000G EUR MAF as proxy for EAF (consistent with TwoSampleMR convention).

echo "Step 1: Converting eQTLGen locus to COJO .ma format..."
c:/R/R-4.3.3/bin/Rscript.exe "${WORK_DIR}/scripts/taskB_03_helper_make_ma.R" "$INPUT_TSV" "$COJO_MA_FILE"

if [ ! -f "$COJO_MA_FILE" ]; then
  echo "ERROR: .ma file not generated at $COJO_MA_FILE"
  exit 1
fi

echo "  Generated: $COJO_MA_FILE ($(wc -l < "$COJO_MA_FILE") lines)"

# -----------------------------------------------------------------------------
# STEP 2 — Run COJO with EUR LD reference
# -----------------------------------------------------------------------------
echo ""
echo "Step 2: Running GCTA-COJO with EUR LD reference..."

"$GCTA_BIN" \
  --bfile "$BFILE_EUR" \
  --cojo-file "$COJO_MA_FILE" \
  --cojo-slct \
  --cojo-p "$COJO_P" \
  --cojo-wind "$COJO_WIND" \
  --cojo-collinear "$COJO_COLLINEAR" \
  --out "${OUT_DIR}/TaskB_03_tshr_cojo_EUR_raw"

echo "  Done. Results: ${OUT_DIR}/TaskB_03_tshr_cojo_EUR_raw.jma.cojo"

# -----------------------------------------------------------------------------
# STEP 3 — Run COJO with EAS LD reference
# -----------------------------------------------------------------------------
echo ""
echo "Step 3: Running GCTA-COJO with EAS LD reference..."

"$GCTA_BIN" \
  --bfile "$BFILE_EAS" \
  --cojo-file "$COJO_MA_FILE" \
  --cojo-slct \
  --cojo-p "$COJO_P" \
  --cojo-wind "$COJO_WIND" \
  --cojo-collinear "$COJO_COLLINEAR" \
  --out "${OUT_DIR}/TaskB_03_tshr_cojo_EAS_raw"

echo "  Done. Results: ${OUT_DIR}/TaskB_03_tshr_cojo_EAS_raw.jma.cojo"

# -----------------------------------------------------------------------------
# STEP 4 — Parse + combine into spec-formatted output
# -----------------------------------------------------------------------------
echo ""
echo "Step 4: Parsing COJO results into spec format..."
c:/R/R-4.3.3/bin/Rscript.exe "${WORK_DIR}/scripts/taskB_03_helper_parse_cojo.R" "$OUT_DIR"

echo ""
echo "Task B.3 complete."
echo "Output: ${OUT_DIR}/TaskB_03_tshr_cojo_independent_signals_v1.tsv"
