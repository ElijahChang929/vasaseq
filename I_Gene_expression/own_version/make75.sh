#!/bin/bash
# Build a 75 nt copy of step 1's output, so steps 2-7 can be re-run on reads the
# same length as the published plate library (75 nt biological).
#
# Truncates the BIOLOGICAL read: step 1 has already removed the 21 nt technical
# prefix (SKIP5=21), so these files are 130 nt and -l 75 leaves 75 nt of insert
# -- directly comparable to the plate. Cutting the raw 151 nt R2 to 75 instead
# would leave 54 nt after the skip and would not be comparable.
#
# Step 1 is NOT re-run: barcode extraction reads R1 and is unaffected by R2's
# length, so its output is only passed through the truncation. Read names, and
# therefore the CB/RX tags every later stage parses, are untouched.
set -eu          # NB not pipefail: `zcat | head` below would SIGPIPE and abort
SRC=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out/cells
DST=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out75/cells
CA=/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/cutadapt
mkdir -p "$DST" "${DST%/cells}/logs"

# .cells is step 1's manifest and is what pipeline.sh's cell_list() reads -- not
# a glob over the FASTQs. Copying only the FASTQs makes every later step report
# "no cell list yet -- run step1 first" and quietly process 0 cells, which looks
# like success in the queue. Copy it with them.
cp "${SRC}/.cells" "${DST}/.cells"

one() {
  f=$1; DST=$2; CA=$3
  b=$(basename "$f")
  [ -s "$DST/$b" ] && return 0
  env -u PYTHONPATH "$CA" -l 75 -o "$DST/$b.part" "$f" >/dev/null 2>&1
  mv "$DST/$b.part" "$DST/$b"          # atomic: a killed job leaves no half file
}
export -f one
ls "$SRC"/*_cbc.fastq.gz | xargs -P 16 -I{} bash -c 'one "$@"' _ {} "$DST" "$CA"

echo "cells: $(ls $DST/*_cbc.fastq.gz | wc -l)"
echo -n "read length: "
zcat "$DST"/ZHA9292A1_011_cbc.fastq.gz 2>/dev/null | awk 'NR%4==2{print length($0); if(NR>4000) exit}' | sort -n | uniq -c | tail -1
a=$(zcat "$SRC"/ZHA9292A1_011_cbc.fastq.gz | wc -l)
b=$(zcat "$DST"/ZHA9292A1_011_cbc.fastq.gz | wc -l)
echo "reads 011: $((a/4)) -> $((b/4))"
