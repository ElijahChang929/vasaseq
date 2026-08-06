#!/bin/bash
###############################################################################
# residual_rrna.sh <CELLDIR> <ANNOT_BED> <OUTFILE>
#
# How much rRNA survived the in-silico depletion?
#
# Step 3 removes reads that align to the rRNA FASTA. Whatever it misses goes on
# to STAR and lands wherever it lands -- so the leak is measurable at step 4, by
# asking how many of the reads that PASSED depletion align to an rRNA feature in
# the genome annotation. That is a different question from step 3's own report,
# which can only describe the reads it caught.
#
#   rRNA        biotype rRNA                   nuclear rRNA that leaked through
#   Mt_rRNA     biotype MtRrna / Mt_rRNA       mitochondrial, depleted on purpose
#   other       everything else                the reads we wanted
#
# The two annotation builds spell the mitochondrial biotype differently --
# `MtRrna` in the mouse BED, `Mt_rRNA` in the mixed one -- so both are matched.
# Getting this wrong is silent: the mouse pattern simply returned 0 features and
# the column would have read 0.000% rather than failing.
#
# Uniquely mapped reads only (NH:i:1). Multimappers are exactly what a repeated
# rDNA array produces, so counting them here would inflate the rRNA share by an
# amount that depends on the aligner's reporting cap rather than on the biology.
# The number this prints is therefore a floor, and it is the comparable one.
#
# Output: cell <TAB> uniq_mapped <TAB> on_rRNA <TAB> on_Mt_rRNA
#
#   sbatch -c 16 --mem=32G -t 1:00:00 --wrap "./residual_rrna.sh <CELLDIR> <BED> out.tsv"
###############################################################################
set -eu

CELLDIR=${1:?need CELLDIR}
BED=${2:?need annotation BED}
OUT=${3:-residual_rrna.tsv}

# samtools needs its module for libbz2; calling the binary by path alone fails
# with "error while loading shared libraries".
source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/camp/apps/eb/modules/all
module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0 2>/dev/null

EBROOT=/camp/apps/eb/software
SAM=${EBROOT}/SAMtools/1.11-GCC-10.2.0/bin/samtools
BT=${EBROOT}/BEDTools/2.30.0-GCC-11.2.0/bin/bedtools
export SAM BT

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# Split the annotation once, not per cell: 384 cells x a 1.2M-row BED is the
# difference between minutes and hours. Sorted for `intersect -sorted`.
awk -F'\t' '$5 ~ /_(MtRrna|Mt_rRNA)_/' "$BED" | sort -k1,1 -k2,2n > "$TMP/mt.bed"
awk -F'\t' '$5 ~ /_(rRNA|Rrna)_/ && $5 !~ /_(MtRrna|Mt_rRNA)_/' "$BED" \
  | sort -k1,1 -k2,2n > "$TMP/rrna.bed"
[ -s "$TMP/rrna.bed" ] || { echo "no rRNA features matched in $BED -- check the biotype spelling" >&2; exit 1; }
export TMP
echo "rRNA features: $(wc -l < "$TMP/rrna.bed")   Mt_rRNA: $(wc -l < "$TMP/mt.bed")" >&2

one() {
    bam=$1
    b=$(basename "$bam"); cell=$(echo "$b" | sed 's/.*_\([0-9]\{3\}\)_.*/\1/')
    d=$(mktemp -d "$TMP/c.XXXXXX")
    # NH:i:1 only -- see the header note on multimappers
    "$SAM" view -h "$bam" | awk '/^@/ || /NH:i:1\t/ || /NH:i:1$/' \
      | "$SAM" view -b - 2>/dev/null | "$BT" bamtobed -i - 2>/dev/null \
      | sort -k1,1 -k2,2n > "$d/r.bed"
    n=$(wc -l < "$d/r.bed")
    r=$("$BT" intersect -u -sorted -a "$d/r.bed" -b "$TMP/rrna.bed" 2>/dev/null | wc -l)
    m=$("$BT" intersect -u -sorted -a "$d/r.bed" -b "$TMP/mt.bed"   2>/dev/null | wc -l)
    printf '%s\t%s\t%s\t%s\n' "$cell" "$n" "$r" "$m"
    rm -rf "$d"
}
export -f one

printf 'cell\tuniq_mapped\ton_rRNA\ton_Mt_rRNA\n' > "$OUT"
ls "$CELLDIR"/*Aligned.out.bam | xargs -P 16 -I{} bash -c 'one "$@"' _ {} \
  | sort -k1,1 >> "$OUT"

awk -F'\t' 'NR>1{u+=$2; r+=$3; m+=$4}
  END{printf "  uniquely mapped %d\n  on rRNA        %d (%.3f%%)\n  on Mt_rRNA     %d (%.3f%%)\n",
             u, r, 100*r/u, m, 100*m/u}' "$OUT"
echo "wrote $OUT"
