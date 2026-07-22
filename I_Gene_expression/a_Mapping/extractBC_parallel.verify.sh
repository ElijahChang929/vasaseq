#!/bin/bash
################################################################################
# extractBC_parallel.verify.sh
#
# Proves extractBC_parallel.sh produces output EQUIVALENT to serial extractBC.sh
# on a slice of the real library: same reads in every per-cell file (compared as
# sorted sets, since parallel merge changes within-cell order) and same log totals.
#
# Usage: extractBC_parallel.verify.sh [nreads] [ncores] [reads_per_shard]
################################################################################
set -euo pipefail

P2S=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/a_Mapping
SRC=/nemo/lab/turnerj/working/guangxin/vasaseq/data/ref/fastq
NREADS=${1:-300000}
NCORES=${2:-8}
RPS=${3:-40000}          # small shards so the slice really is split many ways

T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
cd "$T"
echo "workdir: $T   nreads=$NREADS ncores=$NCORES reads_per_shard=$RPS"

echo "== making a ${NREADS}-read slice of SRR14783059 =="
head_lines=$((NREADS * 4))
# head closes the pipe early -> zcat gets SIGPIPE; don't let pipefail treat that as failure
set +o pipefail
zcat "$SRC/SRR14783059_1.fastq.gz" | head -n "$head_lines" | gzip > SRR14783059_R1.fastq.gz
zcat "$SRC/SRR14783059_2.fastq.gz" | head -n "$head_lines" | gzip > SRR14783059_R2.fastq.gz
set -o pipefail

echo "== serial extractBC.sh =="
mkdir -p out_serial
"$P2S/extractBC.sh" SRR14783059 vasaplate "$P2S" out_serial >/dev/null 2>&1

echo "== parallel extractBC_parallel.sh =="
"$P2S/extractBC_parallel.sh" SRR14783059 vasaplate "$P2S" out_par "$NCORES" "$RPS"

echo "== compare per-cell content as sorted sets =="
mism=0; cells=0; nonempty=0
for f in out_serial/SRR14783059_*_cbc.fastq.gz; do
    cell=$(basename "$f")
    cells=$((cells+1))
    a=$(zcat "out_serial/$cell" | paste - - - - | sort | md5sum | cut -d' ' -f1)
    if [ -f "out_par/$cell" ]; then
        b=$(zcat "out_par/$cell" | paste - - - - | sort | md5sum | cut -d' ' -f1)
    else
        b="MISSING"
    fi
    n=$(zcat "out_serial/$cell" | wc -l); n=$((n/4))
    [ "$n" -gt 0 ] && nonempty=$((nonempty+1))
    if [ "$a" != "$b" ]; then
        echo "  MISMATCH  $cell  (serial reads=$n)"
        mism=$((mism+1))
    fi
done

echo "== compare log totals =="
grep -E "total sequenced|proper barcodes" out_serial/SRR14783059.log | sed 's/^/  serial:   /'
grep -E "total sequenced|proper barcodes" out_par/SRR14783059.log    | sed 's/^/  parallel: /'

st_tot=$(awk -F',' '/total sequenced/{gsub(/ /,"",$2);print $2}' out_serial/SRR14783059.log)
pl_tot=$(awk -F',' '/total sequenced/{gsub(/ /,"",$2);print $2}' out_par/SRR14783059.log)
st_kep=$(awk -F',' '/proper barcodes/{gsub(/ /,"",$2);print $2}' out_serial/SRR14783059.log)
pl_kep=$(awk -F',' '/proper barcodes/{gsub(/ /,"",$2);print $2}' out_par/SRR14783059.log)

echo "== RESULT =="
echo "  cells=$cells (nonempty=$nonempty)  content mismatches=$mism"
echo "  totals: serial tot=$st_tot kept=$st_kep | parallel tot=$pl_tot kept=$pl_kep"
if [ "$mism" -eq 0 ] && [ "$st_tot" = "$pl_tot" ] && [ "$st_kep" = "$pl_kep" ]; then
    echo "  PASS: parallel output is equivalent to serial."
    exit 0
else
    echo "  FAIL: divergence detected."
    exit 1
fi
