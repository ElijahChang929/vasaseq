#!/bin/bash
###############################################################################
# 01_insilico_depletion.sh -- how much of each library step 3 called ribosomal.
#
# Reads nothing but the per-unit `.ribo-map.log` that step 3 already wrote, so
# this is seconds on the login node.
#
# THIS REPLACES A HAND-MADE TABLE. demo_analyze/tables/cross/insilico_depletion.tsv
# has no generator script anywhere in the repo -- the same defect its own README
# records for step4_mapping.tsv. This reproduces it to 16 reads in 25.6M (see
# below) and can be re-run.
#
# RIBOSOMAL = reads - unmapped, NOT the sum of the log's aln/mem categories.
# The log prints "Number of mapped reads:" as a HEADER over five aln/mem
# subcategories, and those sum to one less than reads-unmapped per unit -- a
# category it does not print. Only `unmapped` is definitionally right: it is the
# read count of the .nonRibo fastq, and STAR's "Number of input reads"
# reproduces it exactly (VASA 130 nt: 72,632,400 both ways). The hand-made table
# used the category sum, which is why it reads 25,642,356 where this reads
# 25,642,372.
#
# THE STRAND FLAG IS NOT THE SAME ON THE TWO PROTOCOLS, AND IT MOVES THIS NUMBER
# ------------------------------------------------------------------------------
# The three VASA runs ran `riboread-selection.py` with stranded=y, FLASH-seq with
# stranded=n (00_fs_ribo.sh's header says why: FLASH-seq is genuinely unstranded
# at 49.1-50.5% forward, VASA is not at 76.1%). `y` discards ~24% of VASA's
# ribosomal reads and would discard ~50% of FLASH-seq's, so the FLASH-seq row is
# measured under the flag that is correct for it, not under VASA's.
# **Quote this table with that stated.** The `.nsorted.all-ribo.bam` are kept on
# both sides, so the other flag stays recoverable without re-running bwa.
#
# Output: tables/cross/insilico_depletion.tsv          per dataset
#         tables/cross/insilico_depletion_per_unit.tsv per unit
###############################################################################
set -euo pipefail
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
cd "$ROOT"
source scripts/datasets.sh

mkdir -p "$OUTROOT/tables/cross"
PER=$OUTROOT/tables/cross/insilico_depletion_per_unit.tsv
AGG=$OUTROOT/tables/cross/insilico_depletion.tsv

printf 'dataset\tunit\treads_in\tribosomal\tto_star\tpct\n' > "$PER"
missing=0
for key in $DS_KEYS; do
    label=$(ds_label "$key")
    while IFS=$'\t' read -r unit map_stem ribo_stem; do
        log="${ribo_stem}.ribo-map.log"
        if [ ! -s "$log" ]; then
            echo "  missing $log" >&2; missing=$((missing+1)); continue
        fi
        awk -F'\t' -v ds="$label" -v u="$unit" '
            /Number of reads:/          { gsub(/[^0-9]/,"",$NF); n = $NF }
            /Number of unmapped reads:/ { gsub(/[^0-9]/,"",$NF); m = $NF }
            END { if (n == "" || m == "") { print "unparsable log" > "/dev/stderr"; exit 1 }
                  printf "%s\t%s\t%d\t%d\t%d\t%.4f\n", ds, u, n, n-m, m, 100*(n-m)/n }' "$log" >> "$PER"
    done < <(ds_units "$key")
done
[ "$missing" -eq 0 ] || { echo "$missing units have no ribo-map.log -- table would be partial" >&2; exit 1; }
echo "wrote $PER ($(( $(wc -l < "$PER") - 1 )) units)"

# --- pooled per dataset ------------------------------------------------------
# Pooled, not a mean of per-unit rates: the units differ 20-60x in depth here and
# a unit mean would let the smallest blank weigh as much as the deepest cell.
{ printf 'dataset\tunits\treads_in\tribosomal\tto_star\tpct\n'
  awk -F'\t' 'NR>1 { n[$1]++; a[$1]+=$3; b[$1]+=$4; c[$1]+=$5; if (!($1 in seen)) { seen[$1]=1; ord[++k]=$1 } }
    END { for (i=1;i<=k;i++) { d=ord[i]
            printf "%s\t%d\t%d\t%d\t%d\t%.2f\n", d, n[d], a[d], b[d], c[d], 100*b[d]/a[d] } }' "$PER"
} > "$AGG"
echo "wrote $AGG"

# --- cross-check: to_star must equal STAR's own input count ------------------
# The two numbers are written by different programs two stages apart. If they
# disagree the .nonRibo fastq that STAR read is not the one step 3 wrote.
fail=0
for key in $DS_KEYS; do
    label=$(ds_label "$key")
    s=0; l=0
    while IFS=$'\t' read -r unit map_stem ribo_stem; do
        f="${map_stem}_E99_Log.final.txt"
        [ -s "$f" ] || continue
        l=$((l + $(awk -F'\t' '/Number of input reads/{gsub(/[^0-9]/,"",$2); print $2; exit}' "$f")))
        s=$((s + $(awk -F'\t' -v d="$label" -v u="$unit" '$1==d && $2==u {print $5; exit}' "$PER")))
    done < <(ds_units "$key")
    if [ "$l" -gt 0 ] && [ "$s" != "$l" ]; then
        echo "  MISMATCH $label: step3 says $s reads to STAR, STAR logs say $l" >&2; fail=1
    elif [ "$l" -gt 0 ]; then
        printf '  %-24s %d reads to STAR, agrees with STAR own logs\n' "$label" "$s"
    else
        printf '  %-24s no STAR logs yet -- cross-check skipped\n' "$label"
    fi
done
[ "$fail" = 0 ] || { echo "step 3 and step 4 disagree -- not trusting $AGG" >&2; exit 1; }

echo
column -t -s$'\t' "$AGG"
