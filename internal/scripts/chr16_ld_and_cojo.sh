#!/bin/bash
# =============================================================================
# Task 5 — chr16p11.2 cluster: pairwise LD and conditional analysis
# =============================================================================
# Table S5 states pairwise r2 from the 1000 Genomes EUROPEAN panel
# (0.90 / 0.19 / 0.20) and then concludes the three hits are a single LD signal.
# The r2 values are correct, but r2 = 0.19-0.20 does not support that
# conclusion -- it argues against it.
#
# The candidates were discovered in BBJ (East Asian), so the ancestry-matched
# panel is EAS. This script computes both panels and then settles the question
# with a GCTA-COJO stepwise conditional analysis on the BBJ summary statistics.
#
# Expected output:
#   EUR  r2 = 0.8957 / 0.1897 / 0.2044   (matches Table S5)
#   EAS  r2 = 0.8537 / 0.7612 / 0.8630   (one block)
#   COJO 1 independent signal (rs8050588); conditioning on it leaves the trio at
#        P = 0.928 / 0.839 / 0.806
#
# Requires internal/scripts/prepare_coloc_v2_inputs.sh to have been run, for
# BBJ_byrsid.tsv. Set IN_DIR to that output directory.
#
# Usage : bash internal/scripts/chr16_ld_and_cojo.sh [IN_DIR] [OUT_DIR]
# =============================================================================

set -e
cd "$(dirname "$0")/../.."      # -> c:/ProjectTEDGWAS

IN="${1:-$PWD/internal/coloc_v2_inputs}"
OUT="${2:-$PWD/internal/chr16_task5}"
mkdir -p "$OUT"

PLINK="./tools/plink.exe"
GCTA="./tools/gcta-1.95.1-windows-x86_64/bin/gcta64.exe"
N_BBJ=175465                    # GWAS Catalog GCST90018627

TRIO="$OUT/chr16_trio.txt"
printf "rs4889606\nrs34649473\nrs78924645\n" > "$TRIO"

# ---- 1. pairwise LD in both panels -----------------------------------------
for anc in eur eas; do
  BF="TrackA_MR/v5_upgrade/data/1000G_${anc^^}_phase3/g1000_${anc}"
  if [ ! -f "$BF.bed" ]; then echo "[MISSING] $BF"; continue; fi
  echo "########## 1000G ${anc^^} ##########"
  "$PLINK" --bfile "$BF" --extract "$TRIO" --r2 inter-chr --ld-window-r2 0 \
           --out "$OUT/trio_$anc" --silent
  awk 'NR==1{next}{ printf "  %-12s %-12s r2 = %.4f\n", $3, $6, $7 }' "$OUT/trio_$anc.ld"
  echo "  --- positions (hg19, from .bim) ---"
  awk '$2=="rs4889606"||$2=="rs34649473"||$2=="rs78924645"{
         printf "    %-12s chr%s:%s  %s/%s\n", $2, $1, $4, $5, $6 }' "$BF.bim"
done

# ---- 2. COJO conditional analysis on BBJ, EAS reference ---------------------
echo
echo "########## GCTA-COJO (BBJ chr16, EAS reference) ##########"
# .ma format: SNP A1 A2 freq b se p N
awk -F'\t' -v n="$N_BBJ" '
  NR==1{ print "SNP\tA1\tA2\tfreq\tb\tse\tp\tN"; next }
  $2==16 && $6!="NA" && $7!="NA" && $8!="NA" && $1!="" {
    print $1"\t"$4"\t"$5"\t"$8"\t"$6"\t"$7"\t"$9"\t"n }' \
  "$IN/BBJ_byrsid.tsv" > "$OUT/bbj_chr16.ma"
echo "  .ma rows: $(( $(wc -l < "$OUT/bbj_chr16.ma") - 1 ))"

"$GCTA" --bfile TrackA_MR/v5_upgrade/data/1000G_EAS_phase3/g1000_eas \
        --chr 16 --cojo-file "$OUT/bbj_chr16.ma" --cojo-slct --cojo-p 5e-6 \
        --out "$OUT/bbj_chr16_cojo_eas" 2>&1 | grep -E "selected|independent signals"

echo
echo "  --- independently selected signal(s) ---"
awk 'NR==1{next}{ printf "    %-12s chr%s:%-10s b=%+.4f p=%.2e\n", $2, $1, $3, $11, $13 }' \
    "$OUT/bbj_chr16_cojo_eas.jma.cojo"

echo
echo "  --- trio: marginal vs conditional ---"
printf "    %-12s %-11s %-11s %-11s %-11s\n" SNP b_marginal p_marginal b_cond p_cond
awk -F'\t' 'NR==1{next}
  $2=="rs4889606"||$2=="rs34649473"||$2=="rs78924645"{
    printf "    %-12s %+11.4f %-11.2e %+11.4f %-11.2e\n", $2, $6, $8, $11, $13 }' \
    "$OUT/bbj_chr16_cojo_eas.cma.cojo"

echo
echo "Outputs in: $OUT"
