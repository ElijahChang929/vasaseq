#!/bin/bash
###############################################################################
# datasets.sh -- the four datasets, and where each one's per-unit files are.
#
# SOURCED, not run. Every scanner in this folder gets its job list from here, so
# adding or moving a dataset is one edit and cannot leave two scanners disagreeing
# about which wells or which run directory they scanned. demo_analyze's scanners
# each carried their own copy of that list and they had already drifted -- 173
# wells in one, 172 in another.
#
# THE TWO STEMS, AND WHY IT IS TWO
# --------------------------------
# Every per-unit file in this pipeline is <stem><fixed suffix>, but the stem
# changes at step 3: VASA's step-4 outputs are named after the .nonRibo fastq and
# its step-3 outputs are named after the trimmed one. So a unit carries
#
#   map_stem   ${map_stem}_E99_Aligned.out.bam
#              ${map_stem}_E99_Log.final.txt
#              ${map_stem}_E99_Aligned.out.singlemappers_genes.bed.gz
#              ${map_stem}_E99_Aligned.out.nsorted.multimappers_genes.bed.gz
#   ribo_stem  ${ribo_stem}.ribo-map.log
#              ${ribo_stem}.Ribo.bam
#              ${ribo_stem}.nsorted.all-ribo.bam
#
# On the three VASA runs map_stem = ribo_stem + ".nonRibo"; on FLASH-seq they are
# the same string, because 00_fs_ribo.sh names its outputs off the library stem.
# Nothing else about the four differs at the file level.
#
# ds_units <key> emits, one unit per line:   unit <TAB> map_stem <TAB> ribo_stem
###############################################################################

VASA=/nemo/lab/turnerj/working/guangxin/vasaseq
OWN130_CELLS=$VASA/data/PM26037/out/cells
OWN75_CELLS=$VASA/data/PM26037/out75/cells
# Overridable so the whole folder can be re-run against the E116-remapped plate
# (scripts/00b_plate_e116.sh) without editing a line:
#   PLATE_CELLS=.../plate_e116/cells FOURWAY_OUT=.../fourway_e116 ./run.sh
# The default stays the E99 run, so a bare invocation still reproduces the
# original comparison rather than silently changing reference under someone.
PLATE_CELLS=${PLATE_CELLS:-$VASA/data/ref/fastq_vasaplate/vasaplate_out_v3}
FS_CELLS=$VASA/data/flashseq_vasa/run/nonribo/cells
PERCELL=$VASA/res/vasaplate/per_cell.tsv

# The plate's mouse wells. `ours_v3` is our own anchor run, which is what
# demo_analyze's tables/plate/ and its step-4 scan used -- 173 wells. Keep it
# here too so the two folders' plate rows are the same reads.
PSOURCE=${PSOURCE:-ours_v3}
PRULE=${PRULE:-call_fig1d}

# Where products go. Code root and product root are SEPARATE: every scanner
# does `cd "$ROOT"; source scripts/datasets.sh`, so ROOT has to stay the code
# checkout or the source line breaks. FOURWAY_OUT moves only tables/ and
# figures/, which is what a re-run against a different reference needs.
OUTROOT=${FOURWAY_OUT:-${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}
mkdir -p "$OUTROOT/tables/cross"

DS_KEYS="own130 own75 plate fs"

ds_label() {
    case "$1" in
        own130) echo "VASA own, 130 nt"       ;;
        own75)  echo "VASA own, 75 nt"        ;;
        plate)  echo "VASA published, mouse"  ;;
        fs)     echo "FLASH-seq"              ;;
        *)      echo "unknown dataset key: $1" >&2; return 1 ;;
    esac
}

# The VASA runs are globbed off the BAMs that exist; the plate is restricted to
# its mouse wells by an external call table, and FLASH-seq to the ten delivered
# libraries. A glob that silently returns nothing is the failure mode here, so
# every branch checks it found something.
ds_units() {
    local key=$1 n=0 b u
    case "$key" in
        own130|own75|plate)
            local dir suffix=_cbc_trimmed_homoATCG.nonRibo
            case "$key" in
                own130) dir=$OWN130_CELLS ;;
                own75)  dir=$OWN75_CELLS  ;;
                plate)  dir=$PLATE_CELLS  ;;
            esac
            local want=''
            if [ "$key" = plate ]; then
                # Space-separated, NOT newline-separated: the membership test
                # below is a `case` glob on " $want ", which never matches when
                # the separators are newlines -- and it fails by silently
                # selecting zero wells, not by erroring.
                want=$(awk -F'\t' -v s="$PSOURCE" -v r="$PRULE" '
                    NR==1 { for (i=1;i<=NF;i++) { if ($i=="source") si=i; if ($i=="well") wi=i; if ($i==r) ci=i } }
                    NR>1 && $si==s && $ci=="mouse" { printf "%03d ", $wi }' "$PERCELL")
                [ -n "$want" ] || { echo "no mouse wells in $PERCELL for source=$PSOURCE" >&2; return 1; }
            fi
            for b in "$dir"/*"${suffix}"_E99_Aligned.out.bam; do
                [ -s "$b" ] || continue
                u=$(basename "$b"); u=${u%%_cbc*}; u=${u##*_}
                if [ -n "$want" ]; then case " $want " in *" $u "*) ;; *) continue ;; esac; fi
                local m=${b%_E99_Aligned.out.bam}
                printf '%s\t%s\t%s\n' "$u" "$m" "${m%.nonRibo}"
                n=$((n+1))
            done
            ;;
        fs)
            for b in "$FS_CELLS"/*_cbc_noumi_E99_Aligned.out.bam; do
                [ -s "$b" ] || continue
                u=$(basename "$b"); u=${u%_cbc_noumi_E99_Aligned.out.bam}
                local m=${b%_E99_Aligned.out.bam}
                printf '%s\t%s\t%s\n' "$u" "$m" "$m"
                n=$((n+1))
            done
            ;;
        *) echo "unknown dataset key: $key" >&2; return 1 ;;
    esac
    [ "$n" -gt 0 ] || { echo "ds_units $key: found no units" >&2; return 1; }
}
