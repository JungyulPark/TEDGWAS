#!/bin/bash
# =============================================================================
# Week 1 Master Run Script
# =============================================================================
# Purpose : Orchestrate Week 1 Task A (pQTL) and Task B (TSHR fine-mapping)
#           in correct sequence with parallel execution where possible.
# Usage   : bash run_week1.sh
# =============================================================================

set -e

WORK_DIR="c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade"
SCRIPTS="${WORK_DIR}/scripts"
LOG="${WORK_DIR}/00_meta/analytical_log.md"

START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo "============================================"
echo "TED-TRAP v5 Week 1 — Master Run"
echo "Started: $START_TIME"
echo "============================================"
echo "" >> "$LOG"
echo "## Week 1 master run started: $START_TIME" >> "$LOG"

# -----------------------------------------------------------------------------
# Pre-flight checks
# -----------------------------------------------------------------------------
echo ""
echo "Pre-flight checks:"

REQUIRED_INPUTS=(
  "c:/ProjectTEDGWAS/TrackA_MR/data/eqtl_TSHR_all.csv"
)
for f in "${REQUIRED_INPUTS[@]}"; do
  if [ -f "$f" ]; then
    echo "  [OK] $f"
  else
    echo "  [MISSING] $f"
    echo "ABORT: required input missing"
    exit 1
  fi
done

REQUIRED_TOOLS=(
  "c:/ProjectTEDGWAS/tools/gcta-1.95.1-windows-x86_64/bin/gcta64.exe"
  "c:/ProjectTEDGWAS/tools/plink.exe"
)
for t in "${REQUIRED_TOOLS[@]}"; do
  if [ -f "$t" ]; then
    echo "  [OK] $t"
  else
    echo "  [MISSING] $t — Gemini to install and update paths in scripts"
  fi
done

# -----------------------------------------------------------------------------
# Folder structure
# -----------------------------------------------------------------------------
echo ""
echo "Ensuring folder structure..."
mkdir -p "${WORK_DIR}/00_meta"
mkdir -p "${WORK_DIR}/01_pqtl_check/data" "${WORK_DIR}/01_pqtl_check/results" "${WORK_DIR}/01_pqtl_check/docs"
mkdir -p "${WORK_DIR}/02_tshr_finemap/data" "${WORK_DIR}/02_tshr_finemap/results" "${WORK_DIR}/02_tshr_finemap/docs"
mkdir -p "${WORK_DIR}/03_decision"
mkdir -p "${WORK_DIR}/scripts"
echo "  Done."

# -----------------------------------------------------------------------------
# Task A and Task B can run in parallel (friend mod 1)
# Day 1-2: Task A pQTL check + Task B.1 locus extraction
# Day 3:   Task B.2 LD clumping
# Day 3-4: Task B.3 GCTA-COJO
# Day 4:   Task B.4 SuSiE
# Day 5:   Task B.5 integration + Task C decision memo
# -----------------------------------------------------------------------------

# Task A (independent of Task B)
echo ""
echo "============================================"
echo "Task A: pQTL availability check"
echo "============================================"
c:/R/R-4.3.3/bin/Rscript.exe "${SCRIPTS}/taskA_01_pqtl_check.R"

# Task B sequence
echo ""
echo "============================================"
echo "Task B.1: TSHR locus extraction from eQTLGen"
echo "============================================"
c:/R/R-4.3.3/bin/Rscript.exe "${SCRIPTS}/taskB_01_locus_extract.R"

echo ""
echo "============================================"
echo "Task B.2: LD clumping sensitivity"
echo "============================================"
c:/R/R-4.3.3/bin/Rscript.exe "${SCRIPTS}/taskB_02_ld_clumping.R"

echo ""
echo "============================================"
echo "Task B.3: GCTA-COJO conditional analysis (EUR + EAS)"
echo "============================================"
bash "${SCRIPTS}/taskB_03_gcta_cojo.sh"

echo ""
echo "============================================"
echo "Task B.4: SuSiE fine-mapping (EUR + EAS)"
echo "============================================"
c:/R/R-4.3.3/bin/Rscript.exe "${SCRIPTS}/taskB_04_susie_finemap.R"

echo ""
echo "============================================"
echo "Task B.5: Integration of COJO + SuSiE"
echo "============================================"
c:/R/R-4.3.3/bin/Rscript.exe "${SCRIPTS}/taskB_05_integration.R"

# -----------------------------------------------------------------------------
# Final manifest check
# -----------------------------------------------------------------------------
END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo ""
echo "============================================"
echo "All Week 1 analytical tasks complete"
echo "Started: $START_TIME"
echo "Ended:   $END_TIME"
echo "============================================"
echo "" >> "$LOG"
echo "## Week 1 master run ended: $END_TIME" >> "$LOG"

# Check that all 11 required output files exist
REQUIRED_OUTPUTS=(
  "${WORK_DIR}/01_pqtl_check/results/TaskA_01_pqtl_availability_v1.csv"
  "${WORK_DIR}/02_tshr_finemap/results/TaskB_00_ld_reference_metadata_v1.csv"
  "${WORK_DIR}/02_tshr_finemap/results/TaskB_01_tshr_locus_extraction_v1.tsv"
  "${WORK_DIR}/02_tshr_finemap/results/TaskB_02_ld_clumping_sensitivity_v1.csv"
  "${WORK_DIR}/02_tshr_finemap/results/TaskB_03_tshr_cojo_independent_signals_v1.tsv"
  "${WORK_DIR}/02_tshr_finemap/results/TaskB_04_tshr_susie_credible_sets_v1.tsv"
  "${WORK_DIR}/02_tshr_finemap/results/TaskB_05_tshr_finemap_summary_v1.csv"
)

echo ""
echo "Final output manifest:"
all_ok=1
for f in "${REQUIRED_OUTPUTS[@]}"; do
  if [ -f "$f" ]; then
    echo "  [OK] $f"
  else
    echo "  [MISSING] $f"
    all_ok=0
  fi
done

if [ $all_ok -eq 1 ]; then
  echo ""
  echo "All Week 1 deliverables present."
  echo "Next: run Task C — write decision_table, pre_analysis_plan, summary_memo."
  echo "Then upload deliverables to Claude for verification."
else
  echo ""
  echo "WARNING: Some deliverables missing. Investigate before proceeding to Task C."
fi
