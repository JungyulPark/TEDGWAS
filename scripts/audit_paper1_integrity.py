"""Compatibility entry point for the current reviewed candidate.

Historical phrase checks remain available at commit 33f8e25. The current audit
retains its author-decision and prespecification checks and additionally checks
source values, Word tables, consolidated results and essential limitations.
"""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).resolve().parents[1]/'submission/candidate_20260911_review/reproducibility/audit_submission.py'),run_name='__main__')
