#!/usr/bin/env python3
"""Regenerate submission/candidate_20261011_review/GIT_MANIFEST.json deliberately.

Run only after an authorised change has been rebuilt, inspected and audited.
The manifest lists the root manuscript master and every file in the live
candidate (excluding the manifest itself) with byte size and SHA-256.

Usage: python3 scripts/29_refresh_manifest.py YYYY-MM-DD
"""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "submission" / "candidate_20261011_review"


def entry(p: Path) -> dict:
    data = p.read_bytes()
    return {"path": p.relative_to(ROOT).as_posix(), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main() -> int:
    if len(sys.argv) != 2:
        sys.exit("usage: 29_refresh_manifest.py YYYY-MM-DD")
    files = [ROOT / "MANUSCRIPT_TED_TRAP_v5_MASTER.md"]
    files += sorted(p for p in CANDIDATE.rglob("*") if p.is_file() and p.name != "GIT_MANIFEST.json")
    manifest = {
        "date": sys.argv[1],
        "scope": "Current candidate files plus the sole root manuscript master; excluding this manifest.",
        "files": [entry(p) for p in files],
    }
    out = CANDIDATE / "GIT_MANIFEST.json"
    out.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"{out.relative_to(ROOT)}: {len(manifest['files'])} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
