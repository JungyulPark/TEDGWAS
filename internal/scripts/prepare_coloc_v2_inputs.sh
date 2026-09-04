#!/bin/bash
# =============================================================================
# Task 4 — Prepare pre-filtered inputs for taskE_01b_coloc_v2.R
# =============================================================================
# The eQTLGen release is 3.9 GB gzipped and the three outcome files are
# 0.4-0.8 GB each. Loading them whole into R for nine genes is wasteful, so we
# extract once here.
#
# IMPORTANT: outcome rows are selected by rsID membership, NOT by a position
# window. The outcome files are on GRCh38 while eQTLGen and the 1000G panels are
# on GRCh37 (rs179252 is 14:80,969,641 in the outcomes and 14:81,435,985 in
# eQTLGen). taskE_01_coloc.R applied hg19 windows to hg38 positions, an offset
# of ~466 kb on chr14. Matching on rsID removes the problem entirely, and the
# coloc analysis set is the eQTLGen-outcome intersection in any case.
#
# Usage : bash internal/scripts/prepare_coloc_v2_inputs.sh [OUT_DIR]
# =============================================================================

set -e
cd "$(dirname "$0")/../.."      # -> c:/ProjectTEDGWAS

OUT="${1:-$PWD/internal/coloc_v2_inputs}"
mkdir -p "$OUT"

EQTLGEN="TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
BBJ="TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz"
UKB="TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz"
FIN="finngen_R12_GRAVES_OPHT.gz"
FRQ="TrackA_MR/v5_upgrade/04_druggable_mr/data/g1000_eur_freq.frq"

GENES="TSHR IGF1R CTLA4 TNFSF14 IFNGR1 MAPKAPK5 HSD3B7 PRSS36 VKORC1"

# ---- 1. eQTLGen, nine genes ------------------------------------------------
# cols: 1 Pvalue 2 SNP 3 SNPChr 4 SNPPos 5 AssessedAllele 6 OtherAllele
#       7 Zscore 8 Gene 9 GeneSymbol 10 GeneChr 11 GenePos 12 NrCohorts
#       13 NrSamples 14 FDR 15 BonferroniP
echo "[1/6] eQTLGen -> nine genes"
gzip -dc "$EQTLGEN" | awk -F'\t' -v G="$GENES" '
  BEGIN{ n=split(G,a," "); for(i=1;i<=n;i++) K[a[i]]=1
         print "Pvalue\tSNP\tSNPChr\tSNPPos\tAssessedAllele\tOtherAllele\tZscore\tGeneSymbol\tNrSamples" }
  NR==1{next}
  ($9 in K){ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$13 }' \
  > "$OUT/eqtlgen_9genes.tsv"
echo "      rows: $(( $(wc -l < "$OUT/eqtlgen_9genes.tsv") - 1 ))"

# ---- 2. rsID key set --------------------------------------------------------
echo "[2/6] rsID key set"
cut -f2 "$OUT/eqtlgen_9genes.tsv" | tail -n +2 | sort -u > "$OUT/eqtl_rsids.txt"
echo "      unique rsIDs: $(wc -l < "$OUT/eqtl_rsids.txt")"

emit_header() { printf "rsid\tchr\tpos\tea\toa\tbeta\tse\teaf\tpval\n"; }

# ---- 3. BBJ -----------------------------------------------------------------
# cols: 1 chromosome 2 base_pair_location 3 effect_allele 4 other_allele
#       5 beta 6 standard_error 7 effect_allele_frequency 8 p_value 12 rsid
echo "[3/6] BBJ (GCST90018627)"
gzip -dc "$BBJ" | awk -F'\t' -v L="$OUT/eqtl_rsids.txt" '
  BEGIN{ while((getline r < L)>0) K[r]=1
         print "rsid\tchr\tpos\tea\toa\tbeta\tse\teaf\tpval" }
  NR==1{next}
  ($12 in K){ print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 }' \
  > "$OUT/BBJ_byrsid.tsv"
echo "      rows: $(( $(wc -l < "$OUT/BBJ_byrsid.tsv") - 1 ))"

# ---- 4. UKB -----------------------------------------------------------------
# cols: 2 hm_rsid 14 chromosome 15 base_pair_location 16 effect_allele
#       17 other_allele 18 effect_allele_frequency 19 beta 20 standard_error
#       21 p_value      (the non-hm_ block; beta is on the BOLT-LMM linear
#                        scale and is rescaled to log-odds inside the R script)
echo "[4/6] UKB (GCST90038636)"
gzip -dc "$UKB" | awk -F'\t' -v L="$OUT/eqtl_rsids.txt" '
  BEGIN{ while((getline r < L)>0) K[r]=1
         print "rsid\tchr\tpos\tea\toa\tbeta\tse\teaf\tpval" }
  NR==1{next}
  ($2 in K){ print $2"\t"$14"\t"$15"\t"$16"\t"$17"\t"$19"\t"$20"\t"$18"\t"$21 }' \
  > "$OUT/UKB_byrsid.tsv"
echo "      rows: $(( $(wc -l < "$OUT/UKB_byrsid.tsv") - 1 ))"

# ---- 5. FinnGen -------------------------------------------------------------
# cols: 1 #chrom 2 pos 3 ref 4 alt 5 rsids 7 pval 9 beta 10 sebeta 11 af_alt
# effect allele is alt.
echo "[5/6] FinnGen R12 GRAVES_OPHT"
gzip -dc "$FIN" | awk -F'\t' -v L="$OUT/eqtl_rsids.txt" '
  BEGIN{ while((getline r < L)>0) K[r]=1
         print "rsid\tchr\tpos\tea\toa\tbeta\tse\teaf\tpval" }
  NR==1{next}
  ($5 in K){ print $5"\t"$1"\t"$2"\t"$4"\t"$3"\t"$9"\t"$10"\t"$11"\t"$7 }' \
  > "$OUT/FinnGen_byrsid.tsv"
echo "      rows: $(( $(wc -l < "$OUT/FinnGen_byrsid.tsv") - 1 ))"

# ---- 6. EUR allele frequencies (for the eQTLGen beta/se reconstruction) -----
# eQTLGen is European, so the European panel is the correct reference here.
# Only the chromosomes carrying the nine genes are kept.
echo "[6/6] 1000G EUR frequencies (chr 2,6,12,14,15,16,19)"
awk 'NR==1{ print "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS"; next }
     { c=$1+0
       if(c==2||c==6||c==12||c==14||c==15||c==16||c==19)
         print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' "$FRQ" \
  > "$OUT/eur_freq_subset.tsv"
echo "      rows: $(( $(wc -l < "$OUT/eur_freq_subset.tsv") - 1 ))"

echo
echo "Inputs ready in: $OUT"
echo "Run:  COLOC_V2_IN=\"$OUT\" Rscript TrackA_MR/v5_upgrade/scripts/taskE_01b_coloc_v2.R"
