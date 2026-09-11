#!/usr/bin/env python3
"""Word count of the Paper 1 main text, section by section.

Main text = Introduction + Methods + Results + Discussion (the sections journals
count toward a word limit). Abstract, Declarations, References, Figure Legends
and Tables are counted separately and reported but excluded from the total.

Markdown emphasis markers are stripped so `*TSHR*` counts as one word. The main
total drops sub-heading lines; the "+ headings" line adds them back, which is what
a Word-processor count reports on the rendered .docx. The `**Keywords:**` line is
metadata and is excluded from the abstract, which has its own 250-word limit.

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


def count_words(lines: list[str], headings: bool = False) -> int:
    """Words in a section. With headings=True the sub-headings count too, which
    is what a Word-processor count does on the rendered .docx."""
    body = []
    for line in lines:
        s = line.strip()
        if not s or s == "---":
            continue
        if s.startswith("**Keywords:**"):  # metadata, not part of the abstract
            continue
        if s.startswith("#"):              # sub-headings
            if not headings:
                continue
            s = s.lstrip("#").strip()
        if s.startswith("|"):              # markdown table rows
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
    total = head_total = 0
    for name in MAIN:
        if name not in sections:
            print(f"  MISSING section: {name}")
            continue
        n = count_words(sections[name])
        total += n
        head_total += count_words(sections[name], headings=True)
        print(f"  {name:<14} {n:>6}")
    print(f"  {'MAIN TEXT':<14} {total:>6}   (Endocrine Connections limit 5,000)")
    print(f"  {'  + headings':<14} {head_total:>6}   (what a Word count reports)")

    print("  --- not counted toward the limit ---")
    for name in OTHER:
        if name in sections:
            n = count_words(sections[name])
            flag = "   OVER 250" if name == "Abstract" and n > 250 else ""
            print(f"  {name:<14} {n:>6}{flag}")

    if total > 5000:
        print(f"\nOVER LIMIT by {total - 5000} words")
        return 1
    print(f"\nUnder limit with {5000 - total} words of headroom")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
