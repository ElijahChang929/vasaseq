#!/bin/bash
###############################################################################
# read_length_dist.sh [--prefilter] [CELLDIR] [OUTFILE]
#
# Length of every read once trimming has nothing left to remove.
#
# TWO MODES, AND THE DIFFERENCE IS THE WHOLE POINT
# ------------------------------------------------
#   default      tallies the pass-2 OUTPUT, i.e. the reads that survived -m.
#                The length filter has already removed everything below the
#                floor, so the left tail of this distribution does not exist --
#                which is exactly the part worth seeing.
#
#   --prefilter  re-runs pass 2 with the production adapters but NO -m, so every
#                read is measured at its final trimmed length, including the
#                ones the floor discards. This is the honest distribution: the
#                discarded reads are the bars below TRIM_MINLEN, not a gap.
#
# --prefilter does not touch the run directory. cutadapt writes to a pipe and
# only lengths are kept; nothing is stored and the pipeline's own output is
# never overwritten.
#
# Reads pass 1's output (*_cbc_trimmed.fq.gz) in --prefilter mode and pass 2's
# (*_cbc_trimmed_homoATCG.fq.gz) otherwise.
#
# Output: cell <TAB> length <TAB> reads
#
#   sbatch -c 16 --mem=16G -t 3:00:00 --wrap "./read_length_dist.sh --prefilter <CELLDIR> out.tsv"
###############################################################################
set -eu           # NB not pipefail: a `zcat | awk` that exits early would abort

PRE=no
MODE=vasa
while :; do
  case "${1:-}" in
    --prefilter) PRE=yes; shift ;;
    --mode)      MODE=$2; shift 2 ;;
    *)           break ;;
  esac
done
case "$MODE" in vasa|legacy) ;; *) echo "--mode must be vasa|legacy" >&2; exit 2 ;; esac
CELLDIR=${1:-cells}
OUT=${2:-read_length_dist.tsv}

CUTADAPT=${TRIM_CUTADAPT:-/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/cutadapt}
ADAPTER3=${TRIM_ADAPTER3:-GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT}

# The two pass-2s share nothing but -m, so --prefilter has to be told which one
# produced the data it is re-deriving. Using the wrong set would measure a
# trimming that never ran.
#   vasa    own_version/trim.sh -- measured read-through adapter, 20-mers with
#           min_overlap, a 5' poly-T, --poly-a, -n 10
#   legacy  a_Mapping/trim.sh (the published pipeline, and what the VASA-plate
#           library was run with) -- four 6-mer homopolymers, default -n 1
if [ "$MODE" = legacy ]; then
    CUT_ARGS="--trim-n -a polyG1=GG{5} -a polyC1=CC{5} -a polyT1=TT{5} -a polyA1=AA{5}"
else
    CUT_ARGS="--trim-n -n 10 --poly-a -a rt=${ADAPTER3};min_overlap=8 -a polyA=A{20};min_overlap=10 -a polyG=G{20};min_overlap=10 -g polyT5=T{20};min_overlap=10"
fi
export CUTADAPT CUT_ARGS

# post-filter: just count what is already on disk
one_post() {
    f=$1
    b=$(basename "$f" _cbc_trimmed_homoATCG.fq.gz); cell=${b##*_}
    zcat "$f" | awk -v c="$cell" 'NR%4==2{n[length($0)]++}
                                  END{for (L in n) printf "%s\t%s\t%s\n", c, L, n[L]}'
}

# pre-filter: trim again with no -m, measure, throw the reads away
one_pre() {
    f=$1
    b=$(basename "$f" _cbc_trimmed.fq.gz); cell=${b##*_}
    # -m omitted entirely: cutadapt's default minimum is 0, so nothing is
    # dropped and even fully consumed reads are written at length 0.
    # CUT_ARGS is deliberately unquoted -- it is a pre-split argument list.
    env -u PYTHONPATH "$CUTADAPT" $CUT_ARGS -o - "$f" 2>/dev/null \
      | awk -v c="$cell" 'NR%4==2{n[length($0)]++}
                          END{for (L in n) printf "%s\t%s\t%s\n", c, L, n[L]}'
}
export -f one_post one_pre

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
printf 'cell\tlength\treads\n' > "$OUT"
if [ "$PRE" = yes ]; then
    ls "$CELLDIR"/*_cbc_trimmed.fq.gz \
      | xargs -P 16 -I{} bash -c 'one_pre "$@"' _ {}
else
    ls "$CELLDIR"/*_cbc_trimmed_homoATCG.fq.gz \
      | xargs -P 16 -I{} bash -c 'one_post "$@"' _ {}
fi | sort -k1,1 -k2,2n >> "$OUT"

echo "wrote $OUT  (prefilter=$PRE)"
awk -F'\t' 'NR>1{t+=$3} END{printf "  %d reads over %d rows\n", t, NR-1}' "$OUT"
