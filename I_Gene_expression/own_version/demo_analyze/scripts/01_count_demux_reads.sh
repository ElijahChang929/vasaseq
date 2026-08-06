#!/bin/bash
# Exact read count per cell barcode, straight from the step-1 demultiplexed FASTQs.
# This is the ground-truth "reads assigned to this barcode" number that the
# barplot is built on -- measured here rather than taken from any report.
#
# Output: demux_read_counts.tsv  (cell <TAB> reads)
set -euo pipefail
CELLS=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out/cells
OUT="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)/tables/own130/demux_read_counts.tsv"

count_one() {
    f="$1"
    b=$(basename "$f" _cbc.fastq.gz)          # ZHA9292A1_001
    cell="${b##*_}"                            # 001
    n=$(zcat "$f" | wc -l)
    printf '%s\t%s\n' "$cell" "$((n / 4))"
}
export -f count_one

printf 'cell\treads\n' > "$OUT"
ls "$CELLS"/*_cbc.fastq.gz | xargs -P 8 -I{} bash -c 'count_one "$@"' _ {} | sort -k1,1 >> "$OUT"
echo "wrote $OUT"
