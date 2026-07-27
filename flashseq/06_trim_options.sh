#!/bin/bash
###############################################################################
# 06_trim_options.sh -- would better trimming recover anything, and which fix?
#
# README.md reports two artefacts that survive nf-core's TrimGalore step and
# reach STAR: poly-G (the NovaSeq X two-colour dark-cycle artefact, up to 8% of
# a library) and adapter read-through whose read BEGINS mid-mosaic-end
# ('CTTATACACATCT...') and therefore does not contain 'CTGTCTCTTATA', the
# pattern cutadapt anchors on. It proposed two fixes: `--nextseq-trim 20` for
# the first, and a second adapter pattern without the `CTGTCTCT` prefix for the
# second.
#
# This script tests both against the current settings on one library, because
# "add a second adapter" is the kind of change that looks obviously right and
# is not. Three schemes, identical except for the flags under test:
#
#   cur    the exact call nf-core's TrimGalore made (see 05_rrna_bwa.sh's
#          header for where that command line comes from)
#   polyg  cur + --nextseq-trim=20
#   both   polyg + -a/-A CTTATACACATCT
#
# and four measurements on the results: pairs surviving, bases written, reads
# still carrying mosaic-end sequence at their 3' end, and -- the one that
# settles it -- how far the two artefact populations overlap.
#
# WHAT IT FOUND (ZHA8833A10, the worst poly-G library and a 30 pg rung,
# 90,404 read pairs sampled at stride 256, 2026-07-27):
#
#   scheme  pairs kept   bases written   R1 with mosaic end left at 3'
#   cur     90,130       22,532,184      0.10%
#   polyg   87,689       21,755,261      0.11%
#   both    87,680       22,134,623      37.18%
#
#   * --nextseq-trim=20 removes 2,441 pairs (2.7%) that currently reach STAR.
#     That matches the 2.80% of R2 reads which are >=80% G almost exactly, so
#     it is removing the poly-G population and not much else.
#
#   * THE TWO ARTEFACTS ARE THE SAME READS. Only 1.14% of R1 begins
#     mid-mosaic-end (1,032 reads) -- far less than FastQC's table suggests,
#     because FastQC reads only the first 50 bp -- and 1,019 of those 1,032
#     have a poly-G R2 mate. They are one no-insert artefact seen from two
#     ends. --nextseq-trim=20 already removes them; the extra adapter adds 9
#     pairs.
#
#   * THE SECOND ADAPTER IS ACTIVELY HARMFUL, so README's proposal 2 is
#     withdrawn. 'CTTATACACATCT' is 13 nt against 'CTGTCTCTTATA''s 12, so on
#     ordinary read-through -- 37% of reads, which the current settings trim
#     correctly -- it wins cutadapt's best-match contest and cuts 7 nt LATER,
#     leaving mosaic-end sequence on 37.18% of R1 reads instead of 0.10%. That
#     is also why `both` writes MORE bases than `polyg`, which otherwise looks
#     like it recovered something.
#
# So the recommendation is `--nextseq-trim 20` (cutadapt) / `--2colour 20`
# (TrimGalore) ALONE. What it buys is a cleaner denominator -- 2.7% fewer junk
# reads reaching STAR -- not recovered data: these pairs have no insert to
# recover. Whether that justifies a re-run is a judgement call, and quantifying
# the mapping-rate gain needs the STAR index this run did not save.
#
#   bash code/flashseq/06_trim_options.sh [library]      # default ZHA8833A10
#
# Runs in ~3 min on a login node: one decompress pass at stride 256, then three
# cutadapt calls over 90k pairs.
###############################################################################
set -eu

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$HERE/config.sh"

LIB="${1:-ZHA8833A10}"
STRIDE="${FS_TRIMOPT_STRIDE:-256}"   # coarser than FS_STRIDE: this compares
                                     # schemes against each other, and 90k
                                     # pairs puts the s.e. on a 3% rate at 0.06%
WORK="${FS_TRIMOPT_DIR:-$FS_OUT/trim_options}/$LIB"

MOSAIC="CTGTCTCTTATA"       # Nextera transposase mosaic end, what TrimGalore uses
MIDMOSAIC="CTTATACACATCT"   # the same read-through starting 7 nt in

say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

r1=$(ls "$FS_FASTQ/${LIB}_S"*"_R1_001.fastq.gz")
r2=$(ls "$FS_FASTQ/${LIB}_S"*"_R2_001.fastq.gz")
rm -rf "$WORK"; mkdir -p "$WORK"

say "sampling $LIB at stride $STRIDE"
for m in 1 2; do
    src=$([ "$m" = 1 ] && echo "$r1" || echo "$r2")
    ( set -o pipefail
      zcat "$src" | awk -v s="$STRIDE" 'int((NR-1)/4)%s==0' | gzip -1 > "$WORK/sub_R$m.fq.gz" )
done
N=$(( $(zcat "$WORK/sub_R1.fq.gz" | wc -l) / 4 ))
say "$N read pairs"

run() {
    local name=$1; shift
    "$FS_CUTADAPT" -j "${FS_NCORES:-4}" -e 0.1 -q 20 -O 1 -m 20:20 --pair-filter=any "$@" \
        -o "$WORK/${name}_R1.fq.gz" -p "$WORK/${name}_R2.fq.gz" \
        "$WORK/sub_R1.fq.gz" "$WORK/sub_R2.fq.gz" > "$WORK/${name}.log" 2>&1
}
say "trimming: cur / polyg / both"
run cur   -a "$MOSAIC" -A "$MOSAIC"
run polyg -a "$MOSAIC" -A "$MOSAIC" --nextseq-trim=20
run both  -a "$MOSAIC" -A "$MOSAIC" -a "$MIDMOSAIC" -A "$MIDMOSAIC" --nextseq-trim=20

echo
printf '%-7s %12s %15s %12s\n' scheme pairs_kept bases_written mosaic_left_R1
printf '%-7s %12s %15s %12s\n' ------- ------------ --------------- ------------
for n in cur polyg both; do
    kept=$(awk -F'[ (]+' '/^Pairs written/{print $5}' "$WORK/$n.log")
    bp=$(awk -F'[ (]+' '/^Total written/{print $4}' "$WORK/$n.log")
    # A trimmed read that still ENDS in mosaic-end sequence was cut too late.
    left=$(zcat "$WORK/${n}_R1.fq.gz" \
        | awk 'NR%4==2{n++; if ($0 ~ /CTGTC$|CTGTCT$|CTGTCTC$|CTGTCTCT$/) k++}
               END{printf "%.2f%%", 100*k/n}')
    printf '%-7s %12s %15s %12s\n' "$n" "$kept" "$bp" "$left"
done

echo
say "are the two artefacts the same reads?"
paste <(zcat "$WORK/sub_R1.fq.gz" | awk 'NR%4==2') \
      <(zcat "$WORK/sub_R2.fq.gz" | awk 'NR%4==2') \
    | awk -v mm="$MIDMOSAIC" -v n="$N" '{
        s1 = (index($1, mm) == 1)
        g = $2; gsub(/[^G]/, "", g); p2 = (length(g) >= 0.8 * length($2))
        if (s1) a++; if (p2) b++; if (s1 && p2) c++
      } END {
        printf "  R1 begins mid-mosaic-end : %6d (%.2f%%)\n", a, 100*a/n
        printf "  R2 is >=80%% G            : %6d (%.2f%%)\n", b, 100*b/n
        printf "  both, i.e. same pair     : %6d (%.1f%% of the first)\n", c, 100*c/a
      }'
echo
say "outputs kept in $WORK"
