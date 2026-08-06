#!/bin/bash
# Exact read count per cell barcode for the published VASA-plate library
# (SRR14783059 / GSM5369495, vasaplate-HEK293T-mESC), straight from the step-1
# demultiplexed FASTQs. Same measurement as ../count_demux_reads.sh, 384 cells
# instead of 16.
#
# Step 1 is shared across every run of this library (stages 1-4 are symlinked
# between run directories), so vasaplate_out/ is the right source regardless of
# which downstream run you are looking at.
#
# Output: demux_read_counts.tsv  (well <TAB> reads)
set -euo pipefail
CELLS=/nemo/lab/turnerj/working/guangxin/vasaseq/data/ref/fastq_vasaplate/vasaplate_out
OUT="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)/tables/plate/demux_read_counts.tsv"

count_one() {
    f="$1"
    b=$(basename "$f" _cbc.fastq.gz)          # SRR14783059_001
    well="${b##*_}"                            # 001
    n=$(zcat "$f" | wc -l)
    printf '%s\t%s\n' "$well" "$((n / 4))"
}
export -f count_one

printf 'well\treads\n' > "$OUT"
ls "$CELLS"/*_cbc.fastq.gz | xargs -P 16 -I{} bash -c 'count_one "$@"' _ {} | sort -k1,1 >> "$OUT"
echo "wrote $OUT  ($(( $(wc -l < "$OUT") - 1 )) wells)"
