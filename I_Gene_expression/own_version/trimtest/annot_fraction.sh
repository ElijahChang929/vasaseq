#!/bin/bash
###############################################################################
# annot_fraction.sh [variant ...] -- where do the alignments actually land?
#
# aligned_composition.py catches one artefact class: the poly-A tail aligning
# to a genomic A-tract. Round 7 exposed another that it cannot see -- a read
# like
#
#     TCACATTCGA AAAAAAAA TGTCAC TTTGTA
#     |--insert--|-polyA--|--rc(CBC)--rc(UMI)--|
#
# is 10 nt of real sequence and 19 nt of junk, and STAR will happily align the
# 29 nt as a unit. That read is not A-rich, so it passes the composition check,
# and it inflates the uniquely-mapped count with an alignment that is mostly
# barcode.
#
# The test that does see it: junk has no reason to land inside a gene. This
# counts, for uniquely mapped reads only, the fraction overlapping any feature
# in the annotation BED that step 5 uses. A trimming setting whose extra reads
# are real should hold or raise that fraction; one whose extra reads are junk
# will dilute it.
#
#   ./annot_fraction.sh v1 v17 v22
###############################################################################
set -uo pipefail
cd "${TRIMTEST_DIR:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest}" || exit 1

EBROOT=/camp/apps/eb/software
BED=${TRIMTEST_BED:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/refsorted.bed}
source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0 2>/dev/null
SAM=${EBROOT}/SAMtools/1.11-GCC-10.2.0/bin/samtools
BT=${EBROOT}/BEDTools/2.30.0-GCC-11.2.0/bin/bedtools

CELL=${CELL:-011}
printf "%-6s %10s %12s %7s %12s %7s\n" var uniq inAnnot % exonPC %
for v in "$@"; do
  bam="bam/${v}_${CELL}_Aligned.out.bam"
  [ -s "$bam" ] || { printf "%-6s %s\n" "$v" "-- no bam --"; continue; }
  tmp=$(mktemp -d); trap 'rm -rf "$tmp"' RETURN
  $SAM view -h "$bam" | awk '/^@/ || /NH:i:1\t/ || /NH:i:1$/' \
    | $SAM view -b - 2>/dev/null | $BT bamtobed -i - 2>/dev/null | sort -k1,1 -k2,2n > "$tmp/r.bed"
  n=$(wc -l < "$tmp/r.bed")
  a=$($BT intersect -u -sorted -a "$tmp/r.bed" -b "$BED" 2>/dev/null | wc -l)
  e=$($BT intersect -u -sorted -a "$tmp/r.bed" \
        -b <(awk '$5 ~ /ProteinCoding_exon/' "$BED") 2>/dev/null | wc -l)
        # NB $5, not $4 -- this BED is chr/start/end/STRAND/name/..., so the
        # gene name with its biotype is the fifth column, not the fourth.
  awk -v v="$v" -v n="$n" -v a="$a" -v e="$e" 'BEGIN{
    printf "%-6s %10d %12d %6.1f%% %12d %6.1f%%\n", v, n, a, n?100*a/n:0, e, n?100*e/n:0 }'
  rm -rf "$tmp"; trap - RETURN
done
