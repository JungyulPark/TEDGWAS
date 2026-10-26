# Resident presentation review (25–26 October 2026)

A resident built a 17-slide conference deck from the paper. It was revised to match the final manuscript
(`MANUSCRIPT_TED_TRAP_v5_MASTER.md`, commit f5efd35) and the package tables, figures and workbook. The result is
21 slides: 16 main and 5 supplementary, with Korean speaker notes (v4). The current deck is v6 (26 October):
24 slides with a findings-first Discussion (four slides: TSHR and protection from GD; pinpointing the TSHR signal;
a shared TSHR signal in GD and ophthalmopathy; IGF1R and IGF-1R blockade), a "What we found / What it suggests"
Conclusion and a supplementary Limitations slide. The change report below describes v2 to v4.

## Contents
- [`PRESENTATION_CHANGE_REPORT_20261026.md`](PRESENTATION_CHANGE_REPORT_20261026.md) (renders on GitHub) and
  `PRESENTATION_CHANGE_REPORT_20261026.pdf`: change report in Korean, original deck vs final deck. It covers the
  slide mapping, the 11 required corrections, number checks, tables and figures, format, notes, verification,
  open items and side-by-side thumbnails. The thumbnails mask the third-party mechanism figure on slide 2.
- `review_checklist_20261025.md`: corrections requested for the original deck (first review).
- `deck_review_findings_20261025.md`: independent review of deck v3 (6 reviewers plus adversarial verification;
  70 confirmed, 5 rejected). All confirmed findings were applied in v4, the final deck.
- `report_review_findings_20261026.md`: independent check of the change report itself.
- `thumbs/`: slide images used by the report (original and final; slide 2 mechanism figure masked).
- `figures/`: slide versions of Figures 1B, 2, 3A and S2. They are redrawn from `Supplementary_Dataset.xlsx` with the
  logic of `reproducibility/build_figures_20261014.py`, at slide size, in Carlito (metric-compatible with Calibri).
- `scripts/`:
  - `make_slide_figures.py` draws the figures.
  - `build_deck.py` builds the final deck from the original deck.
  - `make_report.py` and `print_pdf.js` build the report.
  - They were run from the session scratchpad, so their paths are absolute.

## Not committed
Neither deck (`.pptx`) is in the repository. Both carry the hospital's Fasoo document-tracking metadata (a user code and
access times), and this repository is public. Slide 2 also has a mechanism figure whose source has not been
identified. The final deck was delivered to the author directly. Before any external sharing, the institutional
policy on the tracking metadata must be checked.

## Open items for the presenter
- Source of the slide 2 figure.
- Acknowledgement wording and orbital data-use permission.
- Manuscript status on the title slide.
- Presenter affiliation in the notes.
- Figure S1 font: it is the manuscript figure, unchanged.
