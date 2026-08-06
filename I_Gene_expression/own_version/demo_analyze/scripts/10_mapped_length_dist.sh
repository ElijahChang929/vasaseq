#!/bin/bash
###############################################################################
# 10_mapped_length_dist.sh -- read length split by what STAR did with the read.
#
# THE QUESTION
# ------------
# figures/04_mapping/step4_mapping.png says the own library multimaps far more
# than the published plate: 28.6% (130 nt) and 33.0% (75 nt) against 13.3%.
# The obvious suspect is length -- this library's inserts are short and short
# reads place ambiguously. This script measures that instead of assuming it.
#
# NOTHING IS RE-RUN. Every STAR BAM in all three runs was written with
# `--outSAMunmapped Within` (pipeline.sh, a_Mapping/map_star.sh), so the
# unmapped reads and their sequences are already on disk next to the mapped
# ones. This is a read-only tally.
#
# HOW A READ IS CLASSIFIED
# ------------------------
# `-F 0x900` drops secondary (0x100) and supplementary (0x800) records, leaving
# exactly ONE record per read -- that is what makes these totals reconcile with
# STAR's own Log.final.txt. Then:
#
#   flag 0x4 set + uT:A:3  -> toomany    (more loci than --outFilterMultimapNmax)
#   flag 0x4 set           -> unmapped   (uT:A:0/1/2 -- no seed, too short, too many MM)
#   NH:i:1                 -> unique
#   NH:i:>1                -> multi
#
# uT is STAR's unmapped-reason tag, and it is the only way to separate "too many
# loci" from genuinely unmapped: STAR writes both as flag-0x4 records.
#
# Length is `length(SEQ)`. STAR soft-clips and never hard-clips, so SEQ is the
# full read for mapped and unmapped records alike -- verified: SEQ is "*" on
# none of them, including secondaries.
#
# Output: tables/cross/mapped_length_dist.tsv
#             dataset <TAB> cell <TAB> category <TAB> length <TAB> reads
#         tables/cross/step4_mapping_from_bam.tsv
#             the four totals per dataset, reconciled against Log.final.txt
#
# The four categories are kept separate in the TSV. The figures fold `toomany`
# into `multi` -- a read at >20 loci is a multimapper -- but folding here would
# throw the split away and it is 0.4-0.9% worth keeping.
#
# SELF-CHECK, AND IT IS NOT DECORATION
# ------------------------------------
# Per dataset the script sums its own tally and compares it, category by
# category, against the sum of the same cells' *_E99_Log.final.txt. Any
# disagreement is a hard error -- a wrong flag filter or a missed tag would
# otherwise produce a perfectly plausible-looking curve. Measured on
# ZHA9292A1_002 before this script existed: 1196866 / 644160 / 8910 / 352534,
# exact against that cell's log.
#
# When parsing Log.final.txt, anchor on "Number of ...": the "% of ..." lines
# carry the same words and a naive sum silently adds ~29 per cell.
#
# COST (measured on the real run, job 51245883)
#   205 BAMs / 22 GB (9.2 own130 + 7.8 own75 + 5.0 plate) at -P 8, -@ 2:
#   2 min 33 s wall, 86 MB peak RSS. The memory ask below is generous on purpose
#   -- nothing here is held in RAM, only the per-length tally of one BAM.
#
#   sbatch -c 16 --mem=8G -t 60 --wrap="scripts/10_mapped_length_dist.sh"
###############################################################################
set -euo pipefail
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
cd "$ROOT"

VASA=/nemo/lab/turnerj/working/guangxin/vasaseq
OWN130=$VASA/data/PM26037/out/cells
OWN75=$VASA/data/PM26037/out75/cells
PLATE=$VASA/data/ref/fastq_vasaplate/vasaplate_out_v3
PERCELL=$VASA/res/vasaplate/per_cell.tsv

# Which species call types the plate's mouse wells.
#
# `ours_v3` is our own anchor run and is what scripts/04_build_tables_plate.R
# and everything in tables/plate/ uses -- 173 wells. NOTE that the two other
# cross/ products disagree: scripts/05_probe_qc.sh and the hand-made
# tables/cross/step4_mapping.tsv both use `published` (the paper's own calls,
# 172 wells). The two sets differ by exactly one well and 0.09% of the reads,
# so nothing here turns on it, but it is a one-word switch and the README says
# which was used.
PSOURCE=${PSOURCE:-ours_v3}
RULE=call_fig1d          # the paper's Fig. 1d UFI-fraction doublet rule

NPAR=${NPAR:-8}          # BAMs scanned at once
THREADS=${THREADS:-2}    # samtools decompression threads per BAM

OUT=tables/cross/mapped_length_dist.tsv
REC=tables/cross/step4_mapping_from_bam.tsv

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/camp/apps/eb/modules/all
module load SAMtools/1.11-GCC-10.2.0 2>/dev/null
SAM=/camp/apps/eb/software/SAMtools/1.11-GCC-10.2.0/bin/samtools
[ -x "$SAM" ] || { echo "no samtools at $SAM" >&2; exit 1; }

TMP=$(mktemp -d "${TMPDIR:-/tmp}/mlen.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# The job list is whitespace-split by `xargs -n 4`, so the dataset travels as a
# space-free key and only becomes its label here. Changing a label is a one-line
# edit that cannot desynchronise from the keys used to group the scan.
label_of() {
    case "$1" in
        own130) echo "own library, 130 nt"    ;;
        own75)  echo "own library, 75 nt"     ;;
        plate)  echo "published, mouse wells" ;;
        *)      echo "unknown dataset key: $1" >&2; exit 1 ;;
    esac
}

# --- one BAM ---------------------------------------------------------------

tally() {
    key=$1; cell=$2; bam=$3; out=$4
    "$SAM" view -@ "$THREADS" -F 0x900 "$bam" | awk -v ds="$(label_of "$key")" -v cell="$cell" '
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
export -f tally label_of
export SAM THREADS

# --- job list: key cell bam outfile ----------------------------------------
# Paths carry no spaces (checked), so the list is space-separated.

: > "$TMP/jobs"
i=0
add() {   # add <key> <bam>
    b=$(basename "$2"); c=${b%%_cbc*}; c=${c##*_}
    i=$((i+1))
    printf '%s %s %s %s\n' "$1" "$c" "$2" "$TMP/t.$i" >> "$TMP/jobs"
}

for bam in "$OWN130"/*_E99_Aligned.out.bam; do add own130 "$bam"; done
for bam in "$OWN75"/*_E99_Aligned.out.bam;  do add own75  "$bam"; done

# plate: mouse wells only, by the chosen call source
mouse=$(awk -F'\t' -v src="$PSOURCE" -v rule="$RULE" '
    NR==1 { for (i=1;i<=NF;i++) { if ($i=="source") s=i; if ($i=="well") w=i; if ($i==rule) c=i } }
    NR>1 && $s==src && $c=="mouse" { printf "%03d\n", $w }' "$PERCELL")
[ -n "$mouse" ] || { echo "no mouse wells in $PERCELL for source=$PSOURCE" >&2; exit 1; }
for w in $mouse; do
    bam=$PLATE/SRR14783059_${w}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam
    [ -s "$bam" ] || { echo "missing $bam" >&2; exit 1; }
    add plate "$bam"
done

echo "scanning $(wc -l < "$TMP/jobs") BAMs, $NPAR at a time (plate source=$PSOURCE)..."
# `set -eo pipefail` inside the child, not just out here: without it a samtools
# that fails to start (the classic is a missing module, i.e. no libbz2) leaves
# awk to exit 0 over an empty stream, and the shard file is silently empty.
xargs -P "$NPAR" -n 4 bash -c 'set -eo pipefail; tally "$0" "$1" "$2" "$3"' < "$TMP/jobs"

{ printf 'dataset\tcell\tcategory\tlength\treads\n'
  cat "$TMP"/t.* | sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4n
} > "$OUT"
echo "wrote $OUT ($(( $(wc -l < "$OUT") - 1 )) rows)"

# --- reconcile against STAR's own logs -------------------------------------

# The logs are streamed through one `cat` rather than handed to awk as an
# argument list: 384-cell runs are what made ARG_MAX a real bug here once.
logsum() {   # -> "input unique multi toomany unmapped" over one dataset's cells
    awk -v k="$1" '$1==k {print $3}' "$TMP/jobs" \
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
bamsum() {   # the same five numbers out of the tally we just wrote
    awk -F'\t' -v d="$1" '$1==d { s[$3] += $5; tot += $5 }
      END { printf "%d %d %d %d %d\n", tot, s["unique"], s["multi"], s["toomany"], s["unmapped"] }' "$OUT"
}

printf 'dataset\tcells\tinput\tunique\tmulti\ttoomany\tunmapped\tpct_multi\tpct_unmapped\n' > "$REC"
fail=0
for key in own130 own75 plate; do
    ds=$(label_of "$key")
    n=$(awk -v k="$key" '$1==k' "$TMP/jobs" | wc -l)
    read -r bi bu bm bt bx <<< "$(bamsum "$ds")"
    read -r li lu lm lt lx <<< "$(logsum "$key")"
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
