# Corrections requested for the resident's original deck (to be checked in the revised deck)

Original deck text: ../out/text.md (17 slides). Revised deck text: deck_text.md (21 slides).

## Must fix
1. Authors: Kyung-Hwa Shin removed as author; Suk-Woo Yang sole corresponding author; Jungyul Park and Min-Seon Kim
   co-first (†); spelling Min-Seon (not Min-Sun); "Manuscript under review" is wrong (not yet submitted); acknowledge Shin if
   orbital data are shown.
2. Five disease datasets, not three: add FinnGen GD (E4_GRAVES_STRICT, 3,962/496,386) and strict ophthalmopathy
   (E4_GRAVES_OPHT_STRICT, 753/499,595); FinnGen endpoints overlap (not independent replications); UKB ancestry not established
   as exclusively European; FinnGen Finnish; GO controls include GD without recorded eye disease.
3. TSHR colocalization: FinnGen GD PP.H4 0.264 and UKB 0.226 (H3 dominant); BBJ 0.951->0.661 and strict 0.891->0.451 at
   p12 = 1e-6; "strong colocalization in BBJ and FinnGen" must say FinnGen ophthalmopathy (GO), not FinnGen GD.
4. IGF1R: add FinnGen GD 1.07 (0.86-1.33), P = 0.539, PP.H4 0.010 and strict 1.52 (0.92-2.49), P = 0.0996, 0.038; FinnGen GD
   is the most precise and near null (detectable OR ~1.36); PP.H2 dominant in four of five (0.571-0.731); two of four
   single-instrument omissions removed nominal significance in BBJ and UKB.
5. CTLA4: BBJ MR P = 5.45e-15 passed the screen threshold; only PP.H4 (0.201) was below 0.80; UKB/GO P nominal and dependent on
   rs13030124 (omitting it P = 0.323-0.972).
6. Figures: replace old Figure 1 and Figure 2 (three datasets, numbers in figure, old stacking) with current versions; no callout
   covering the plot.
7. Treatment wording: remove the unsupported claim that TSHR inhibitors are effective in TED; teprotumumab "reduced proptosis in
   randomized trials in active TED"; mechanism figure needs a source.
8. Causal language: association wording throughout (slides and notes); exposure is blood mRNA, not protein.
9. OR 0.12 per SD: explain per-allele ORs 0.77-0.82 (moderate).
10. Prior studies: Ji et al. (FinnGen R10, strict ophthalmopathy) and Zhang 2026 (FinnGen R12, same TSHR pattern) must appear.
11. Conclusions: findings only, include non-colocalization; limitations: FinnGen overlap, European eQTL with East Asian BBJ,
    single TSHR instrument.

## Minor
- Gene symbols and P italic; protein IGF-1R roman; "Graves disease" spelling consistent.
- Slide numbers as automatic fields, not typed.
- No shapes covering content (old slides 5, 9, 11, 14).
- Colocalization slide: "causal" not "casual"; priors p1 = p2 = 1e-4, p12 = 1e-5; single-causal-variant assumption.
- MHC slide / 13 discoveries: add FinnGen GD comparison; MHC and 16p11.2 genes not independent.
- Notes typos (IGFR1, FinGen, etc.); internal file names removed from notes.
