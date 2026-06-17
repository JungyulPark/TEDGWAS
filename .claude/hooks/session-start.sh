#!/bin/bash
# SessionStart hook: surface the Paper 1 submission-ready integrity audit.
# Read-only and non-blocking — never fails the session.
set -uo pipefail

# Web sessions only (no-op locally)
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "${CLAUDE_PROJECT_DIR:-.}" || exit 0

echo "── TED-TRAP Paper 1 submission-ready audit (SessionStart) ──"
if [ -f scripts/audit_paper1_integrity.py ]; then
  python3 scripts/audit_paper1_integrity.py \
    || echo "⚠️  Integrity audit reported FAIL — review before submitting."
else
  echo "(audit script not found; skipping)"
fi
exit 0
