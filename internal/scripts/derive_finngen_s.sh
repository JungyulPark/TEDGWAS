#!/bin/bash
# =============================================================================
# Task 2 — Derive the case fraction s directly from FinnGen summary-stat files
# =============================================================================
# Purpose : Establish s = cases/(cases+controls) WITHOUT relying on external
#           documentation, using only the allele-frequency columns present in
#           the FinnGen .gz files themselves.
#
# Method  : af_alt is the sample-size-weighted mean of the case and control
#           allele frequencies:
#
#               af_alt = s * af_alt_cases + (1 - s) * af_alt_controls
#           =>  s      = (af_alt - af_alt_controls)
#                        / (af_alt_cases - af_alt_controls)
#
#           Rows are filtered to those where the case/control frequency
#           difference is large enough that 6-significant-digit rounding in the
#           reported frequencies contributes negligible error.
#
# Usage   : bash internal/scripts/derive_finngen_s.sh
# Output  : stdout only (median / Q1 / Q3 of per-SNP s for each file)
#
# Result recorded in internal/FINNGEN_endpoint_provenance.md section 2.
# =============================================================================

set -e
cd "$(dirname "$0")/../.."   # -> c:/ProjectTEDGWAS

FILES=(
  finngen_R12_GRAVES_OPHT.gz
  finngen_R12_E4_GRAVES_OPHT_STRICT.gz
  finngen_R12_E4_GRAVES_STRICT.gz
)

# FinnGen R12 total cohort size (finngen.fi release announcement)
N_R12=500348

# Columns: 1 #chrom 2 pos 3 ref 4 alt 5 rsids 6 nearest_genes 7 pval 8 mlogp
#          9 beta 10 sebeta 11 af_alt 12 af_alt_cases 13 af_alt_controls
AF=11; AF_CASE=12; AF_CTRL=13

ROWS=400000        # rows scanned per file
MIN_DIFF=0.0005    # minimum |af_cases - af_controls| to keep a row
MAX_AF=0.05        # restrict to lower-frequency variants (best precision)

for f in "${FILES[@]}"; do
  if [ ! -f "$f" ]; then
    echo "[MISSING] $f — skipped"
    continue
  fi
  echo "########## $f ##########"

  tmp=$(mktemp)
  gzip -dc "$f" \
    | head -"$ROWS" \
    | awk -F'\t' -v a="$AF" -v c="$AF_CASE" -v o="$AF_CTRL" \
                 -v mind="$MIN_DIFF" -v maxaf="$MAX_AF" '
        NR == 1 { next }                      # header
        {
          af = $a + 0; ac = $c + 0; ao = $o + 0;
          d  = ac - ao;
          if (d > mind && af > 0 && af < maxaf) {
            s = (af - ao) / d;
            if (s > 0 && s < 0.2) print s     # sanity bound
          }
        }' \
    | sort -g > "$tmp"

  n=$(wc -l < "$tmp")
  echo "  informative SNPs : $n"
  awk -v n="$n" -v ntot="$N_R12" '
    NR == int(n*0.25) { q1  = $1 }
    NR == int(n*0.50) { med = $1 }
    NR == int(n*0.75) { q3  = $1 }
    { sum += $1 }
    END {
      printf "  median s         = %.8f\n", med
      printf "  Q1 / Q3          = %.8f / %.8f\n", q1, q3
      printf "  mean s           = %.8f\n", sum/NR
      printf "  implied cases at N=%d : %.2f  -> %d\n", ntot, med*ntot, int(med*ntot + 0.5)
    }' "$tmp"
  rm -f "$tmp"
  echo
done

echo "Cross-check against FinnGen R12 PheWeb (official GWAS analysis set):"
echo "  GRAVES_OPHT            858 cases / 499490 controls"
echo "  E4_GRAVES_OPHT_STRICT  753 cases / 499595 controls"
echo "  E4_GRAVES_STRICT      3962 cases / 496386 controls"
