#!/bin/bash
# Driver for shortspecies_from_bed.sh: one worker per library (the helper loops
# both arms internally), then concatenate and pool. Kept as a file rather than
# inlined in a job script because the xargs {} placeholders and awk braces do not
# survive being embedded in a formatted string.
set -eo pipefail
ROOT=${1:?scratch root}
HELPER=${2:?path to shortspecies_from_bed.sh}
OUT=${3:?output tsv}
shift 3
LIBS="$*"

rm -rf parts; mkdir -p parts
echo "$LIBS" | tr ' ' '\n' | xargs -P 10 -I@ bash -c \
    "'$HELPER' '$ROOT' 'parts/ss_@.tsv' @ > 'parts/ss_@.log' 2>&1 || echo 'FAILED @'"

first=$(ls parts/ss_*.tsv | head -1)
head -1 "$first" > "$OUT"
for f in parts/ss_*.tsv; do tail -n +2 "$f"; done >> "$OUT"
echo "rows written: $(( $(wc -l < "$OUT") - 1 ))"

echo "=== pooled per arm x biotype ==="
awk -F'\t' 'NR>1{t[$1"\t"$3]+=$4; i[$1"\t"$3]+=$5}
  END{printf "%-9s %-10s %14s %14s %9s\n","arm","biotype","rows_total","rows_jsIN","pctIN";
      for(k in t){split(k,p,"\t"); printf "%-9s %-10s %14d %14d %8.2f%%\n",p[1],p[2],t[k],i[k],100*i[k]/t[k]}}' "$OUT" | sort
