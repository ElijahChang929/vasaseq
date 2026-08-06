#!/bin/bash
###############################################################################
# 03_mapped_length_dist.sh -- read length split by what STAR did with the read,
# over all four datasets.
#
# Same tally as demo_analyze/scripts/10_mapped_length_dist.sh, which is
# validated; the only change is that the job list comes from scripts/datasets.sh
# instead of being written out three times inside the script.
#
# NOTHING IS RE-RUN. Every STAR BAM in all four runs was written with
# `--outSAMunmapped Within`, so the unmapped reads and their sequences are on
# disk next to the mapped ones. Read-only tally.
#
# HOW A READ IS CLASSIFIED
#   `-F 0x900` drops secondary (0x100) and supplementary (0x800), leaving exactly
#   ONE record per read -- that is what makes these totals reconcile with STAR's
#   own Log.final.txt. Then:
#     flag 0x4 + uT:A:3 -> toomany     flag 0x4 -> unmapped
#     NH:i:1            -> unique      NH:i:>1  -> multi
#   uT is STAR's unmapped-reason tag and the only way to separate "too many loci"
#   from genuinely unmapped: STAR writes both as flag-0x4 records.
#
# Length is `length(SEQ)`. STAR soft-clips and never hard-clips, so SEQ is the
# full read for mapped and unmapped records alike.
#
# ALL FOUR ARE NOW ON THE SAME PIPELINE INPUT. FLASH-seq reaches STAR through
# 00_fs_ribo.sh, i.e. with its rRNA removed by the same script and the same
# reference as the three VASA runs. Before that stage existed this comparison
# would have put a pipeline difference (rRNA in one STAR input, not the other)
# on the same axis as the protocol difference -- and rRNA multimaps, which is
# precisely what the step-4 figures measure.
#
# SELF-CHECK, AND IT IS NOT DECORATION
#   Per dataset the tally is compared, category by category, against the sum of
#   the same units' *_E99_Log.final.txt. Any disagreement is a hard error -- a
#   wrong flag filter or a missed tag would otherwise produce a perfectly
#   plausible-looking curve.
#   When parsing Log.final.txt anchor on "Number of ...": the "% of ..." lines
#   carry the same words and a naive sum silently adds ~29 per unit.
#
# Output: tables/cross/mapped_length_dist.tsv    dataset unit category length reads
#         tables/cross/step4_mapping.tsv         the four totals per dataset
#
#   sbatch -c 16 --mem=8G -t 120 --wrap="scripts/03_mapped_length_dist.sh"
###############################################################################
set -euo pipefail
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
cd "$ROOT"
source scripts/datasets.sh

NPAR=${NPAR:-8}
THREADS=${THREADS:-2}
OUT=$OUTROOT/tables/cross/mapped_length_dist.tsv
REC=$OUTROOT/tables/cross/step4_mapping.tsv
mkdir -p "$OUTROOT/tables/cross"

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/camp/apps/eb/modules/all
module load SAMtools/1.11-GCC-10.2.0 2>/dev/null
SAM=/camp/apps/eb/software/SAMtools/1.11-GCC-10.2.0/bin/samtools
[ -x "$SAM" ] || { echo "no samtools at $SAM" >&2; exit 1; }

TMP=$(mktemp -d "${TMPDIR:-/tmp}/mlen4.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

tally() {
    ds=$1; unit=$2; bam=$3; out=$4
    "$SAM" view -@ "$THREADS" -F 0x900 "$bam" | awk -v ds="$ds" -v cell="$unit" '
      { nh = 0; ut = "-"
        for (i = 12; i <= NF; i++) {
            if      (substr($i,1,5) == "NH:i:") nh = substr($i,6) + 0
            else if (substr($i,1,5) == "uT:A:") ut = substr($i,6) }
        if (int($2/4) % 2) cat = (ut == "3") ? "toomany"  : "unmapped"
        else               cat = (nh == 1)   ? "unique"   : "multi"
        n[cat SUBSEP length($10)]++ }
      END { for (k in n) { split(k, a, SUBSEP)
              printf "%s\t%s\t%s\t%s\t%d\n", ds, cell, a[1], a[2], n[k] } }' > "$out"
}
export -f tally
export SAM THREADS

# --- job list ---------------------------------------------------------------
# Tab-separated and read with `xargs -d '\n'` + a shell that splits on tabs:
# dataset LABELS contain spaces, so the whitespace-splitting `xargs -n 4` that
# demo_analyze used only worked because it passed space-free KEYS and resolved
# the label inside the worker. Carrying the label directly is simpler and is
# safe as long as the delimiter is not whitespace.
: > "$TMP/jobs"
i=0
for key in $DS_KEYS; do
    label=$(ds_label "$key")
    while IFS=$'\t' read -r unit map_stem ribo_stem; do
        bam="${map_stem}_E99_Aligned.out.bam"
        [ -s "$bam" ] || { echo "missing $bam" >&2; exit 1; }
        i=$((i+1))
        printf '%s\t%s\t%s\t%s\n' "$label" "$unit" "$bam" "$TMP/t.$i" >> "$TMP/jobs"
    done < <(ds_units "$key")
done

echo "scanning $(wc -l < "$TMP/jobs") BAMs, $NPAR at a time..."
# `set -eo pipefail` inside the child, not just out here: without it a samtools
# that fails to start leaves awk to exit 0 over an empty stream and the shard
# file is silently empty.
xargs -d '\n' -P "$NPAR" -I{} bash -c 'set -eo pipefail; IFS=$'"'"'\t'"'"' read -r a b c d <<< "{}"; tally "$a" "$b" "$c" "$d"' < "$TMP/jobs"

{ printf 'dataset\tunit\tcategory\tlength\treads\n'
  cat "$TMP"/t.* | sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4n
} > "$OUT"
echo "wrote $OUT ($(( $(wc -l < "$OUT") - 1 )) rows)"

# --- reconcile against STAR's own logs --------------------------------------
logsum() {   # -> "input unique multi toomany unmapped" over one dataset's units
    awk -F'\t' -v d="$1" '$1==d {print $3}' "$TMP/jobs" \
      | sed 's/_Aligned\.out\.bam$/_Log.final.txt/' \
      | while read -r f; do cat "$f"; done \
      | awk -F'\t' '
          /Number of input reads/                          { i += $2 }
          /Uniquely mapped reads number/                   { u += $2 }
          /Number of reads mapped to multiple loci/        { m += $2 }
          /Number of reads mapped to too many loci/        { t += $2 }
          /Number of reads unmapped: too many mismatches/  { x += $2 }
          /Number of reads unmapped: too short/            { x += $2 }
          /Number of reads unmapped: other/                { x += $2 }
          END { printf "%d %d %d %d %d\n", i, u, m, t, x }'
}
bamsum() {
    awk -F'\t' -v d="$1" '$1==d { s[$3] += $5; tot += $5 }
      END { printf "%d %d %d %d %d\n", tot, s["unique"], s["multi"], s["toomany"], s["unmapped"] }' "$OUT"
}

printf 'dataset\tunits\tinput\tunique\tmulti\ttoomany\tunmapped\tpct_multi\tpct_unmapped\n' > "$REC"
fail=0
for key in $DS_KEYS; do
    ds=$(ds_label "$key")
    n=$(awk -F'\t' -v d="$ds" '$1==d' "$TMP/jobs" | wc -l)
    read -r bi bu bm bt bx <<< "$(bamsum "$ds")"
    read -r li lu lm lt lx <<< "$(logsum "$ds")"
    if [ "$bi $bu $bm $bt $bx" != "$li $lu $lm $lt $lx" ]; then
        echo "RECONCILE FAILED  $ds" >&2
        echo "  bam: $bi $bu $bm $bt $bx" >&2
        echo "  log: $li $lu $lm $lt $lx" >&2
        fail=1
    fi
    printf '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n' "$ds" "$n" "$bi" "$bu" "$bm" "$bt" "$bx" \
        "$(echo "100*($bm+$bt)/$bi" | bc -l)" "$(echo "100*$bx/$bi" | bc -l)" >> "$REC"
done
[ "$fail" = 0 ] || { echo "tally does not match STAR's logs -- not trusting $OUT" >&2; exit 1; }
echo "reconciled against Log.final.txt on every dataset."
column -t -s$'\t' "$REC"
echo "wrote $REC"
