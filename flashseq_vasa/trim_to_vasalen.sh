#!/bin/bash
###############################################################################
# trim_to_vasalen.sh -- hard-trim an adapter-trimmed FLASH-seq R1 fastq so that
# its read-length distribution reproduces VASA's STAR-input distribution.
#
# This is the READ-LENGTH CONTROL for the species axis of the comparison, and
# the only thing that makes any claim about tRNA/snoRNA/snRNA/miRNA testable.
# The confound it controls for is structural, not biological:
#
#   step 5 tags a read jS:IN only if the read is CONTAINED in the feature; step 6
#   keeps a non-spliceable biotype ONLY when jS == IN; and in the v2 BED 98.5% of
#   tRNA features, 99.2% of miRNA, 96.5% of snoRNA and 84.6% of snRNA are shorter
#   than one 151 nt read. So at native length those species are suppressed by
#   arithmetic before biology gets a vote.
#
# HOW THE MATCH IS DONE, exactly
# ------------------------------
# One length is drawn PER READ from a 10,000-entry lookup table built by
# measure_vasa_readlen.py from the pooled STAR-input read lengths of VASA's 12
# real cells, and the read is truncated to it. Drawing per read rather than
# trimming to one fixed number is the whole point: VASA's ability to satisfy
# jS:IN on a sub-151 nt feature lives in the SHORT TAIL of its distribution, and
# trimming everything to VASA's median (~110-130 nt) would still exceed nearly
# every short feature and would recover nothing.
#
# Truncation keeps the FIRST L bases, i.e. it cuts from the 3' end. The read's 5'
# start is therefore unchanged, which is the end STAR anchors an unspliced
# alignment on.
#
# A read already shorter than its drawn L is kept whole. That is not a fudge --
# it is what "match the distribution" means when the source is already short --
# but it does mean the output distribution is the pointwise MINIMUM of the input
# and the target, so it can only be shorter than the target, never longer. The
# script counts these and prints the rate; at 151 nt input against a target whose
# p99 is well under 151 the rate is small, and the stats block is where to check
# that rather than assume it.
#
# REPRODUCIBLE, NOT RANDOM. srand(seed) with a fixed seed over a fixed input in
# fixed order gives byte-identical output on every run. The seed is a parameter
# so the draw is auditable and so a second draw can be taken deliberately.
#
# Usage: trim_to_vasalen.sh <in_R1.fq.gz> <lut.txt> <out_R1.fq.gz> <seed> [stats.tsv]
###############################################################################
set -euo pipefail

if [ $# -lt 4 ]; then
    sed -n '2,45p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 1
fi

IN=$1
LUT=$2
OUT=$3
SEED=$4
STATS=${5:-}

[ -s "$IN" ]  || { echo "FATAL: no input fastq $IN" >&2; exit 1; }
[ -s "$LUT" ] || { echo "FATAL: no LUT $LUT -- run measure_vasa_readlen.py first" >&2; exit 1; }

# gzip -1: these are large intermediates on scratch, read once by STAR. The
# decompression cost dominates either way and -1 halves the write time.
#
# NB the awk reads the LUT in its own BEGIN block rather than taking it on stdin,
# so the fastq stream stays clean.
zcat "$IN" \
  | awk -v lut="$LUT" -v seed="$SEED" -v statsf="$STATS" '
    BEGIN {
        n=0
        while ((getline line < lut) > 0) { if (line != "") { L[++n]=line+0 } }
        close(lut)
        if (n == 0) { print "FATAL: LUT is empty" > "/dev/stderr"; exit 1 }
        srand(seed)
        nrec=0; ntrunc=0; nshort=0; bin=0; bout=0
    }
    # Records arrive in fours. The draw is taken on the sequence line and reused
    # on the quality line, so seq and qual are always cut to the same length.
    NR%4==1 { print; next }
    NR%4==2 {
        seq=$0; li=length(seq)
        # int(rand()*n)+1 is a uniform draw over 1..n; rand() never returns 1.0.
        want=L[int(rand()*n)+1]
        keep = (li > want) ? want : li
        if (li > want) { ntrunc++ } else if (li < want) { nshort++ }
        nrec++; bin+=li; bout+=keep
        cut=keep
        print substr(seq,1,keep)
        next
    }
    NR%4==3 { print; next }
    NR%4==0 { print substr($0,1,cut); next }
    END {
        if (statsf != "") {
            printf "reads\t%d\nbases_in\t%d\nbases_out\t%d\nmean_len_in\t%.3f\nmean_len_out\t%.3f\ntruncated\t%d\nalready_shorter_than_draw\t%d\nlut_entries\t%d\nseed\t%s\n", \
                   nrec, bin, bout, (nrec>0?bin/nrec:0), (nrec>0?bout/nrec:0), ntrunc, nshort, n, seed > statsf
            close(statsf)
        }
        printf "trim_to_vasalen: %d reads, mean %.2f -> %.2f nt, %d truncated (%.2f%%), %d already shorter than their draw (%.2f%%)\n", \
               nrec, (nrec>0?bin/nrec:0), (nrec>0?bout/nrec:0), \
               ntrunc, (nrec>0?100.0*ntrunc/nrec:0), \
               nshort, (nrec>0?100.0*nshort/nrec:0) > "/dev/stderr"
    }' \
  | gzip -1 > "$OUT"

# The pipeline is under pipefail, so a truncated zcat or a failed awk fails here
# rather than silently producing a short fastq. Verify the record count anyway:
# a fastq whose line count is not a multiple of 4 would be silently accepted by
# STAR's reader on some versions.
nl=$(zcat "$OUT" | wc -l)
if [ $(( nl % 4 )) -ne 0 ]; then
    echo "FATAL: $OUT has $nl lines, not a multiple of 4" >&2
    exit 1
fi
echo "trim_to_vasalen: wrote $OUT ($(( nl / 4 )) records)"
