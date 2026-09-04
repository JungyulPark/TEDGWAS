#!/usr/bin/env python3
"""Build the submission package from the locked master.

Journals want the main text and the supplementary material as separate files, so
this splits the single master at "## Supplementary Material" and renders each half
through pandoc, alongside the cover letter. It never edits the master — the master
stays the one source of truth and everything here is derived from it.

Produces, in submission/:
    MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx   title page through Tables 1-3
    SUPPLEMENTARY_MATERIAL.docx              Tables S1-S9 and supplementary legends
    COVER_LETTER_EndocrineConnections.docx   from its markdown source
    figures/Figure{1,2,3,S1,S2}.png          copied from the canonical figure directory

Usage:  python3 scripts/29_build_submission_package.py
Exit 0 if every file was written and passed its content check.
"""
from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER = f"{ROOT}/MANUSCRIPT_TED_TRAP_v5_MASTER.md"
SUB = f"{ROOT}/submission"
FIGSRC = f"{ROOT}/TrackA_MR/v5_upgrade/07_manuscript/figures"
BUILD = f"{ROOT}/submission/_build"

SPLIT = "## Supplementary Material"
FIGURES = {"Figure1.png": "Figure1_study_design.png",
           "Figure2.png": None,                       # built outside this repo
           "Figure3.png": "Figure3.png",
           "FigureS1.png": "FigureS1.png",
           "FigureS2.png": "FigureS2.png"}

fails: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(f"  {'ok ' if ok else 'FAIL'}  {msg}")
    if not ok:
        fails.append(msg)


def docx_text(path: str) -> str:
    x = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    return re.sub(r"<[^>]+>", "", x)


def docx_tables(path: str) -> int:
    return zipfile.ZipFile(path).read("word/document.xml").decode("utf-8").count("<w:tbl>")


def pandoc(src: str, out: str) -> None:
    subprocess.run(["pandoc", src, "-o", out, "--standalone"], check=True)


def main() -> int:
    text = open(MASTER, encoding="utf-8").read()
    md5 = hashlib.md5(text.encode()).hexdigest()
    print(f"== Submission package ==\nmaster md5: {md5}\n")

    assert SPLIT in text, "master has no Supplementary Material section"
    main_md, supp_md = text.split(SPLIT, 1)
    os.makedirs(BUILD, exist_ok=True)
    os.makedirs(f"{SUB}/figures", exist_ok=True)

    # ---------------- main text ----------------
    main_path = f"{BUILD}/main.md"
    open(main_path, "w", encoding="utf-8").write(main_md.rstrip() + "\n")
    main_docx = f"{SUB}/MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx"
    pandoc(main_path, main_docx)

    mt = docx_text(main_docx)
    check("Running head:" in mt, "main: opens with the running head")
    for s in ("## Abstract".lstrip("# "), "Introduction", "Methods", "Results",
              "Discussion", "Declarations", "References", "Figure Legends", "Tables"):
        check(s in mt, f"main: section present — {s}")
    check(docx_tables(main_docx) == 3, f"main: Tables 1-3 present — got {docx_tables(main_docx)}")
    # the main text legitimately cross-references the supplementary tables; what must
    # not be there is the supplementary section itself
    check("Supplementary Material" not in mt, "main: supplementary section not included")
    check("Detectable effect range of the druggable-gene-wide screen" not in mt,
          "main: supplementary table captions not included")

    # ---------------- supplementary ----------------
    supp_path = f"{BUILD}/supplementary.md"
    header = ("# Supplementary Material\n\n"
              "**Genetic susceptibility and therapeutic target biology diverge at "
              "*TSHR* and *IGF1R* in Graves disease and thyroid eye disease**\n\n"
              "Jungyul Park, Min-Seon Kim, Kyung-Hwa Shin, Suk-Woo Yang\n\n"
              "---\n")
    open(supp_path, "w", encoding="utf-8").write(header + supp_md.lstrip() + "\n")
    supp_docx = f"{SUB}/SUPPLEMENTARY_MATERIAL.docx"
    pandoc(supp_path, supp_docx)

    st = docx_text(supp_docx)
    for i in range(1, 10):
        check(f"Table S{i}." in st, f"supplementary: Table S{i} present")
    for f in ("Figure S1.", "Figure S2."):
        check(f in st, f"supplementary: {f} legend present")
    check("## Introduction" not in st and "Bahn RS" not in st,
          "supplementary: main text and reference list not duplicated")

    # ---------------- cover letter ----------------
    cl_md = f"{SUB}/COVER_LETTER_EndocrineConnections.md"
    cl_docx = f"{SUB}/COVER_LETTER_EndocrineConnections.docx"
    if os.path.exists(cl_md):
        pandoc(cl_md, cl_docx)
        ct = docx_text(cl_docx)
        check("Yae-Eun Kang" not in ct, "cover letter: no superseded author name")
        check(len(ct) > 1500, "cover letter: non-trivial length")
    else:
        check(False, "cover letter markdown missing")

    # ---------------- figures ----------------
    for dest, src in FIGURES.items():
        d = f"{SUB}/figures/{dest}"
        if src:
            s = f"{FIGSRC}/{src}"
            if os.path.exists(s):
                shutil.copy2(s, d)
        check(os.path.exists(d) and os.path.getsize(d) > 20_000, f"figure present — {dest}")

    print("\n  -- package --")
    for f in sorted(os.listdir(SUB)):
        p = f"{SUB}/{f}"
        if os.path.isfile(p):
            print(f"    {os.path.getsize(p):>9,} B  {f}")
    for f in sorted(os.listdir(f"{SUB}/figures")):
        print(f"    {os.path.getsize(f'{SUB}/figures/{f}'):>9,} B  figures/{f}")

    print("\nRESULT:", "ALL PASS" if not fails else f"{len(fails)} FAIL")
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
