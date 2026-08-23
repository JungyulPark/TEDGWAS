#!/usr/bin/env python3
"""Word count of the Paper 1 main text, section by section.

Main text = Introduction + Methods + Results + Discussion (the sections journals
count toward a word limit). Abstract, Declarations, References, Figure Legends
and Tables are counted separately and reported but excluded from the total.

Markdown emphasis/heading markers are stripped before counting so that
`*TSHR*` counts as one word and `### Heading` lines are dropped, matching how a
Word-processor count behaves on the rendered .docx.

Usage: python3 scripts/26_wordcount_main_text.py [path/to/master.md]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

MAIN = ["Introduction", "Methods", "Results", "Discussion"]
OTHER = ["Abstract", "Declarations", "References", "Figure Legends", "Tables"]


def split_sections(text: str) -> dict[str, list[str]]:
    """Split on level-2 headings (## X); subheadings stay inside their section."""
    sections: dict[str, list[str]] = {}
    current = "_front"
    sections[current] = []
    for line in text.splitlines():
        m = re.match(r"^##\s+(?!#)(.+?)\s*$", line)
        if m:
            current = m.group(1).strip()
            sections.setdefault(current, [])
            continue
        sections.setdefault(current, []).append(line)
    return sections


def count_words(lines: list[str]) -> int:
    body = []
    for line in lines:
        s = line.strip()
        if not s or s == "---":
            continue
        if s.startswith("#"):          # sub-headings do not count
            continue
        if s.startswith("|"):          # markdown table rows
            continue
        body.append(s)
    text = " ".join(body)
    text = re.sub(r"[*_`]", "", text)  # emphasis markers
    text = re.sub(r"\s+", " ", text)
    return len([w for w in text.split(" ") if w])


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("MANUSCRIPT_TED_TRAP_v5_MASTER.md")
    sections = split_sections(path.read_text(encoding="utf-8"))

    print(f"== Word count: {path} ==")
    total = 0
    for name in MAIN:
        if name not in sections:
            print(f"  MISSING section: {name}")
            continue
        n = count_words(sections[name])
        total += n
        print(f"  {name:<14} {n:>6}")
    print(f"  {'MAIN TEXT':<14} {total:>6}   (Endocrine Connections limit 5,000)")

    print("  --- not counted toward the limit ---")
    for name in OTHER:
        if name in sections:
            print(f"  {name:<14} {count_words(sections[name]):>6}")

    if total > 5000:
        print(f"\nOVER LIMIT by {total - 5000} words")
        return 1
    print(f"\nUnder limit with {5000 - total} words of headroom")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
