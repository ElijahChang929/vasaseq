#!/bin/bash
# 01_build_intervals.sh -- make the probe-target interval file for BOTH references.
#
# The probe-addressable intervals were computed once (see README section "how the
# reference was made") against the MOUSE-ONLY rRNA reference. The published plate
# was mapped against the HUMAN+MOUSE mixed reference, whose mouse contigs carry a
# "mouse_" prefix AND drop the gene symbol, e.g.
#     mouse-only : ENSMUSG00000064337_mt-Rnr1_Mt_rRNA(+)
#     mixed      : mouse_ENSMUSG00000064337_Mt_rRNA(+)
# so the join key is the ENSMUSG id, not the string. This script translates the
# interval file to mixed-reference names and refuses to proceed on any unmatched
# contig.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

MO_FA=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa
MIX_FA=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mixed/unique_rRNA_human_mouse.v2.fa
IN=$HERE/probe_target_intervals.mouse.tsv
OUT=$HERE/probe_target_intervals.mixed.tsv

for f in "$MO_FA" "$MIX_FA" "$IN"; do [ -s "$f" ] || { echo "MISSING: $f" >&2; exit 1; }; done

# mixed-reference mouse contig names, keyed by ENSMUSG id (47S keyed by itself)
grep '^>' "$MIX_FA" | sed 's/^>//' | awk '{print $1}' | grep '^mouse' \
  | awk '{ name=$0
           if (match(name, /ENSMUSG[0-9]+/)) key=substr(name, RSTART, RLENGTH)
           else key=name                       # e.g. mouse_rDNA_47S_BK000964.3_1-13403
           print key "\t" name }' > /tmp/mixmap.$$

awk -F'\t' -v MAP=/tmp/mixmap.$$ '
  BEGIN{ while ((getline line < MAP) > 0) { split(line, a, "\t"); m[a[1]]=a[2] } close(MAP) }
  NR==1 { print; next }
  { name=$1
    if (match(name, /ENSMUSG[0-9]+/)) key=substr(name, RSTART, RLENGTH); else key=name
    if (!(key in m)) { print "UNMATCHED CONTIG: " name > "/dev/stderr"; bad++; next }
    $1=m[key]; OFS="\t"; print
  }
  END{ if (bad) { print "FAILED: " bad " contigs could not be translated" > "/dev/stderr"; exit 1 } }
' "$IN" > "$OUT"
rm -f /tmp/mixmap.$$

echo "wrote $OUT"
echo "  mouse-only intervals : $(( $(wc -l < "$IN") - 1 ))"
echo "  mixed-ref intervals  : $(( $(wc -l < "$OUT") - 1 ))"
echo "  contigs (mouse-only) : $(tail -n +2 "$IN"  | cut -f1 | sort -u | wc -l)"
echo "  contigs (mixed)      : $(tail -n +2 "$OUT" | cut -f1 | sort -u | wc -l)"
echo "  bp addressable       : $(tail -n +2 "$OUT" | awk -F'\t' '{s+=$3-$2+1} END{print s}')"
