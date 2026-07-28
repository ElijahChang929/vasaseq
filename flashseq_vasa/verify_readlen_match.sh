#!/bin/bash
###############################################################################
# verify_readlen_match.sh -- did the vasalen arm actually reproduce VASA's
# read-length distribution, and did it move the quantity the containment rule
# turns on?
#
# The vasalen arm draws a target length per read from VASA's own distribution and
# truncates to it, so its OUTPUT is the pointwise minimum of (input, target) --
# it can only be shorter than VASA, never longer. That is the conservative
# direction for this control, but "conservative" is not a number, so both arms'
# actual distributions are measured here against the target rather than asserted.
#
# Two quantities, both measured the same way on both sides:
#
#   1. STAR-input read length, over the WHOLE prepped fastq of every library.
#      Directly comparable with vasa_starinput_len_lut.txt.hist.tsv, which is the
#      pooled distribution of VASA's 12 real cells measured from its own step-3
#      fastqs.
#
#   2. BED interval span (End - Start) at stride 64. This is what step 5's
#      containment test actually compares against a feature, and it is NOT read
#      length -- a spliced alignment's interval spans its intron, which is why
#      the VASA distribution has p50 = 130 but a max in the millions. The
#      comparable published figure is 43.24% of VASA rows at >= 140 nt (pooled
#      real cells, same stride, same quantity).
#
# Usage: verify_readlen_match.sh <scratch_root> <libs...>
###############################################################################
set -eo pipefail

ROOT=${1:?scratch root, e.g. /nemo/lab/turnerj/scratch/zhangg/vasaseq/flashseq_vasa}
shift
LIBS="$*"
[ -n "$LIBS" ] || { echo "give at least one library" >&2; exit 1; }

for arm in native vasalen; do
    S="$ROOT/$arm/cells"
    [ -d "$S" ] || { echo "no $S -- skipping $arm"; continue; }

    # --- read length, whole files, 10 libraries at a time --------------------
    for lib in $LIBS; do
        echo "$lib"
    done | xargs -P 10 -I{} bash -c \
        "zcat '$S/{}_cbc_noumi_R1.fq.gz' | awk 'NR%4==2{h[length(\$0)]++} END{for(l in h) print l\"\t\"h[l]}' > './len_${arm}_{}.tsv'"
    awk -F'\t' '{h[$1]+=$2} END{for(l in h) print l"\t"h[l]}' ./len_${arm}_*.tsv \
        | sort -n > "./readlen_pooled_${arm}.tsv"
    rm -f ./len_${arm}_*.tsv

    # --- BED span, stride 64, same quantity/stride as the VASA measurement ---
    for lib in $LIBS; do
        echo "$lib"
    done | xargs -P 10 -I{} bash -c \
        "zcat '$S/{}_cbc_noumi_E99_Aligned.out.singlemappers_genes.bed.gz' | awk 'NR%64==1{h[\$3-\$2]++} END{for(l in h) print l\"\t\"h[l]}' > './sp_${arm}_{}.tsv'"
    awk -F'\t' '{h[$1]+=$2} END{for(l in h) print l"\t"h[l]}' ./sp_${arm}_*.tsv \
        | sort -n > "./bedspan_pooled_${arm}.tsv"
    rm -f ./sp_${arm}_*.tsv

    echo "=== $arm ==="
    awk -F'\t' '{n+=$2; s+=$1*$2} END{printf "  readlen: n=%d mean=%.3f\n",n,s/n}' "./readlen_pooled_${arm}.tsv"
    awk -F'\t' '{n+=$2; if($1>=140) g+=$2} END{printf "  bedspan: n=%d  >=140nt=%d (%.2f%%)\n",n,g,100*g/n}' "./bedspan_pooled_${arm}.tsv"
done
ls -la ./readlen_pooled_*.tsv ./bedspan_pooled_*.tsv
