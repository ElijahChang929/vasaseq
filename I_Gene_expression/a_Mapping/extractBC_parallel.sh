#!/bin/bash
################################################################################
# extractBC_parallel.sh
#
# Data-parallel drop-in for extractBC.sh (vasaplate protocol).
#
# extractBC.sh runs ONE single-core concatenator.py over the whole library
# (hours for a full plate). This wrapper does the identical work in parallel:
#
#   1. split the raw R1/R2 into shards of N reads each (both mates cut at the
#      SAME read boundaries so they stay paired), gzipping each shard;
#   2. run one concatenator.py per shard, up to $ncores at a time (each
#      concatenator is itself single-core — parallelism is across shards);
#   3. merge the per-cell outputs across shards (gzip members concatenate, so a
#      plain `cat a.gz b.gz` is a valid .gz) and aggregate the log.
#
# OUTPUT EQUIVALENCE: every read lands in the same per-cell file as the serial
# version (the per-read barcode logic is untouched — it IS concatenator.py).
# Only the ORDER of reads within a cell file differs, because shards are merged
# in shard order. Downstream UMI counting is order-independent, so this is
# equivalent, though not byte-identical. Verified by extractBC_parallel.verify.sh.
#
# USAGE (same first four args as extractBC.sh, plus two optional tuning args):
#   extractBC_parallel.sh <inputroot> <protocol> <path2scripts> <outdir> \
#                         [ncores] [reads_per_shard]
#   - inputroot      : fastq prefix; globs ${inputroot}*_R1*.fastq.gz in CWD
#   - protocol       : only 'vasaplate' is implemented here
#   - path2scripts   : dir holding concatenator.py + bc_celseq2.tsv
#   - outdir         : where the final <inputroot>_<cellID>_cbc.fastq.gz go
#   - ncores         : concatenators to run concurrently (default 8)
#   - reads_per_shard: reads per shard (default 4000000 -> ~56 shards for a
#                      224M-read plate; more shards = finer load balancing)
#
# Requires python3 (with numpy+pandas) already on PATH — the caller activates
# the conda env, exactly as for the serial extractBC.sh.
################################################################################
set -euo pipefail

if [ $# -lt 4 ]; then
    echo "usage: extractBC_parallel.sh <inputroot> <protocol> <path2scripts> <outdir> [ncores] [reads_per_shard]"
    exit 1
fi

inroot=$1
protocol=$2
p2s=$3
outdir=$4
ncores=${5:-8}
reads_per_shard=${6:-4000000}

if [ "$protocol" != "vasaplate" ]; then
    echo "extractBC_parallel.sh: only the 'vasaplate' protocol is implemented (got '$protocol')."
    echo "Use extractBC.sh for other protocols, or extend this wrapper."
    exit 1
fi

lines_per_shard=$((reads_per_shard * 4))   # keep whole 4-line fastq records intact

r1=$(ls ${inroot}*_R1*.fastq.gz 2>/dev/null | head -1 || true)
r2=$(ls ${inroot}*_R2*.fastq.gz 2>/dev/null | head -1 || true)
if [ -z "$r1" ] || [ -z "$r2" ]; then
    echo "extractBC_parallel.sh: could not find ${inroot}*_R1*/_R2*.fastq.gz in $PWD"
    exit 1
fi

mkdir -p "$outdir"
# absolute paths: the per-shard workers cd into the shards dir, so every path
# handed to them (outdir, scripts dir) must be absolute to survive the cd.
outdir=$(readlink -f "$outdir")
p2s=$(readlink -f "$p2s")
work="${outdir}/.ptmp_${inroot}_$$"
mkdir -p "$work/shards" "$work/out"

# ---- exported for the split --filter subshell and the xargs worker ----
export work p2s outdir

echo "[extractBC_parallel] $(date) splitting $r1 / $r2 into shards of ${reads_per_shard} reads ..."
# Split each mate at identical read boundaries, gzip each shard (-1 = fast; these
# are throwaway intermediates). R1 and R2 use the SAME numeric prefix 'sh' with
# different suffixes so shard k is the pair sh<k>_R1.fastq.gz / sh<k>_R2.fastq.gz.
( zcat "$r1" | split -l "$lines_per_shard" -d -a 5 --additional-suffix=_R1.fastq \
      --filter='gzip -1 > "$work/shards/$FILE.gz"' - sh ) &
sp1=$!
( zcat "$r2" | split -l "$lines_per_shard" -d -a 5 --additional-suffix=_R2.fastq \
      --filter='gzip -1 > "$work/shards/$FILE.gz"' - sh ) &
sp2=$!
wait $sp1
wait $sp2

# shard prefixes present (e.g. sh00000, sh00001, ...)
prefixes=$(cd "$work/shards" && ls sh*_R1.fastq.gz | sed 's/_R1.fastq.gz$//' | sort -u)
nsh=$(echo "$prefixes" | wc -l)
echo "[extractBC_parallel] $(date) $nsh shards; running concatenator.py ${ncores}-way ..."

# one concatenator per shard; glob resolves in the shards dir; each writes its
# own out/<prefix>/ so no two workers touch the same file.
run_one() {
    local pfx=$1
    mkdir -p "$work/out/$pfx"
    ( cd "$work/shards" && \
      python3 "$p2s/concatenator.py" --fqf "$pfx" --cbcfile "$p2s/bc_celseq2.tsv" \
              --cbchd 0 --lenumi 6 --umifirst --demux --outdir "$work/out/$pfx" \
      ) > "$work/out/$pfx.stdout" 2>&1
}
export -f run_one
echo "$prefixes" | xargs -P "$ncores" -I{} bash -c 'run_one "$@"' _ {}

echo "[extractBC_parallel] $(date) merging per-cell files across shards ..."
# Collect the per-shard output dirs ONCE (sorted, so the merge is deterministic).
# NB: never glob all shards*cells (~43k files for a full plate) on one command
# line -- that overflows ARG_MAX. We walk shard dirs instead, and pipe each
# cell's shard files through xargs so `cat` is batched regardless of shard count.
mapfile -t shdirs < <(find "$work"/out -mindepth 1 -maxdepth 1 -type d | sort)
# every shard opens all whitelist cells, so one shard dir lists the full cell set
cellids=$(ls "${shdirs[0]}"/*_cbc.fastq.gz | sed -E 's/.*_([0-9]{3})_cbc\.fastq\.gz$/\1/' | sort -u)
for cell in $cellids; do
    # one glob per shard dir -> one path each; xargs cat keeps it ARG_MAX-safe.
    # concatenated gzip members form one valid .gz
    for d in "${shdirs[@]}"; do
        printf '%s\n' "$d"/*_"${cell}"_cbc.fastq.gz
    done | xargs cat > "${outdir%/}/${inroot}_${cell}_cbc.fastq.gz"
done

echo "[extractBC_parallel] $(date) aggregating log ..."
tot=0; kept=0
while IFS= read -r lg; do
    t=$(awk -F',' '/total sequenced reads/{gsub(/ /,"",$2); print $2}' "$lg")
    k=$(awk -F',' '/reads with proper barcodes/{gsub(/ /,"",$2); print $2}' "$lg")
    tot=$((tot + ${t:-0}))
    kept=$((kept + ${k:-0}))
done < <(find "$work"/out -name '*.log' | sort)
frac=$(awk -v k="$kept" -v t="$tot" 'BEGIN{ if (t>0) printf "%.16f", k/t; else print "0" }')
{
    echo "=> to generate cbc file <="
    echo "fastq file:, ${inroot}, "
    echo "full barcode in:, R1, "
    echo "biological read in:, R2, "
    echo "cell specific barcode length:, 8, "
    echo "umi length:, 6, "
    echo "umi goes first:, True, "
    echo "total sequenced reads:, ${tot}, "
    echo "reads with proper barcodes:, ${kept}, ${frac}, "
    echo "# produced by extractBC_parallel.sh (${nsh} shards, ${ncores}-way)"
} > "${outdir%/}/${inroot}.log"

rm -rf "$work"
echo "[extractBC_parallel] $(date) done. total=${tot} kept=${kept} (${frac}) -> ${outdir}"
