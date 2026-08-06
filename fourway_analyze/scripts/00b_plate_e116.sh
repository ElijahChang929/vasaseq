#!/bin/bash
###############################################################################
# 00b_plate_e116.sh -- re-map the published plate's mouse wells onto the SAME
# reference the other three datasets use, then submit stages 4-7.
#
# WHY
# ---
# Three of the four datasets are quantified on Ensembl 116 / GRCm39; the
# published plate was on Ensembl 99 / GRCm38 (human+mouse). That single
# difference sat on every biotype axis in this folder and was NOT small:
# E116 carries 32,889 mouse lncRNA genes against E99's 9,959 (3.30x), and the
# measured lncRNA fraction tracked it -- plate 2.66% on E99 against our
# 8.3-8.4% on E116, with the plate landing inside the published filter window
# ([0.01,0.03], filterParams.py) and ours three times outside it.
#
# So the plate is re-mapped rather than the comparison being caveated.
#
# MOUSE-ONLY, AND WHAT THAT COSTS -- MEASURED, NOT ASSUMED
# --------------------------------------------------------
# This library is a HEK293T/mESC mixing control, so a mouse-only index gives
# human reads nowhere correct to go. The question is whether they vanish or
# mismap onto mouse. Measured on 2026-08-06 by mapping two human-called wells
# (97.3% and 97.8% human by UFI) against this exact index:
#
#   well 225 (human)  2,824,850 reads  ->  16.88% uniquely mapped
#   well 121 (human)  2,198,198 reads  ->  17.24% uniquely mapped
#   well 002 (mouse)     46,559 reads  ->  76.28% uniquely mapped
#   well 004 (mouse)    491,999 reads  ->  81.08% uniquely mapped
#
# Human wells still hold ~2.5% mouse, which at 78% contributes ~1.9 pp, so
# roughly 15% of genuinely human reads mismap onto mouse. Folding that into
# each mouse well's own off-species share gives the spurious fraction of its
# uniquely mapped pool: median 1.08%, 90th pct 1.79%, max 6.56%.
#
# ~1% against biotype differences of 5-15 pp. Kept, and reported, rather than
# filtered: the alternative (a mixed GRCh38+GRCm39 index) would put the plate
# on a different index from the other three and trade a measured 1% for an
# unmeasured asymmetry in multimapping competition.
#
# All 173 mouse-called wells are carried, matching scripts/datasets.sh's plate
# definition exactly so the two folders' plate rows stay the same wells. The
# ten wells above 10% off-species can be dropped at table level afterwards --
# baking that filter into a 10-hour re-run would be the wrong place for it.
#
# STAGES 1-3 ARE REUSED, NOT RE-RUN
# ----------------------------------
# rRNA depletion does not depend on the annotation, and the plate's 384
# nonRibo FASTQs are on disk (55 GB; the 173 mouse wells are 2.4 GB). Only
# stages 4-7 are re-run, against own_version/config.sh's own defaults --
# which already point at star_index_151_r116 and the E116 v2 BED, so the plate
# ends up on literally the same index and BED as own130/own75/fs.
#
# One deliberate deviation to record: that index has sjdbOverhang 150 and the
# published run used 73. 150 is long for a 75 nt read and costs a little
# junction sensitivity -- but own75 already runs at 75 nt on this same index,
# so the deviation cancels between the two datasets being compared.
#
# COST (from the previous full-plate run, 384 wells; this is 173)
#   stage 4 gmap   16 chunks x ~5 min, one genome load
#   stage 5 assign minutes
#   stage 6 pickle 8 h 54 m / 57 GB at 384 wells
#   stage 7 tables 4 h 38 m / 86 GB at 384 wells
#
# Usage:  scripts/00b_plate_e116.sh [setup|submit|status]
###############################################################################
set -euo pipefail

VASA=/nemo/lab/turnerj/working/guangxin/vasaseq
OWN=$VASA/code/I_Gene_expression/own_version
SRC=$VASA/data/ref/fastq_vasaplate/vasaplate_out_v3
OUT=${PLATE_E116_OUT:-$VASA/data/ref/fastq_vasaplate/plate_e116}
LIB=SRR14783059
PERCELL=$VASA/res/vasaplate/per_cell.tsv

CELLS=$OUT/cells
export SAMPLE=$LIB OUTDIR=$OUT

setup() {
    mkdir -p "$CELLS" "$OUT/logs"
    # The mouse wells, by the same call table and rule scripts/datasets.sh uses
    # (source ours_v3, rule call_fig1d) -- derived here rather than copied, so
    # the two cannot drift apart silently.
    awk -F'\t' 'NR>1 && $2=="ours_v3" && $10=="mouse"{print $1}' "$PERCELL" | sort -u > "$OUT/.wells"
    local n; n=$(wc -l < "$OUT/.wells")
    [ "$n" -gt 0 ] || { echo "no mouse wells parsed from $PERCELL" >&2; exit 1; }
    echo "mouse wells: $n"

    : > "$CELLS/.cells"
    local miss=0 w f
    while read -r w; do
        f=$SRC/${LIB}_${w}_cbc_trimmed_homoATCG.nonRibo.fastq.gz
        if [ -s "$f" ]; then
            ln -sf "$f" "$CELLS/${LIB}_${w}_cbc_trimmed_homoATCG.nonRibo.fastq.gz"
            echo "${LIB}_${w}" >> "$CELLS/.cells"
        else
            echo "  MISSING stage-3 output for well $w" >&2; miss=$((miss+1))
        fi
    done < "$OUT/.wells"
    [ "$miss" -eq 0 ] || { echo "$miss wells have no nonRibo FASTQ -- refusing to submit" >&2; exit 1; }
    echo "linked $(wc -l < "$CELLS/.cells") cells into $CELLS"
    echo "total input: $(du -Lsh "$CELLS" | cut -f1)"
}

submit() {
    [ -s "$CELLS/.cells" ] || { echo "run setup first" >&2; exit 1; }
    local E="SAMPLE=$LIB OUTDIR=$OUT"
    local j4 j5 j6 j7
    j4=$(sbatch --parsable -J pe116_4 --chdir="$OWN" -c 16 --mem=80G -t 6:00:00 \
         -o "$OUT/logs/step4.%j.out" --wrap "$E STAR_THREADS=16 ./pipeline.sh step4")
    j5=$(sbatch --parsable -J pe116_5 --chdir="$OWN" -c 16 --mem=32G -t 4:00:00 \
         --dependency=afterok:$j4 -o "$OUT/logs/step5.%j.out" --wrap "$E NCORES=16 ./pipeline.sh step5")
    # stage 6 was 8h54/57GB on 384 wells; 173 wells with headroom for the peak
    # being driven by distinct UMI x gene combinations rather than well count.
    j6=$(sbatch --parsable -J pe116_6 --chdir="$OWN" -c 8 --mem=160G -t 16:00:00 \
         --dependency=afterok:$j5 -o "$OUT/logs/step6.%j.out" --wrap "$E ./pipeline.sh step6")
    j7=$(sbatch --parsable -J pe116_7 --chdir="$OWN" -c 4 --mem=160G -t 10:00:00 \
         --dependency=afterok:$j6 -o "$OUT/logs/step7.%j.out" --wrap "$E ./pipeline.sh step7")
    echo "submitted: step4=$j4 step5=$j5 step6=$j6 step7=$j7"
    echo "$j4 $j5 $j6 $j7" > "$OUT/.jobs"
    echo
    echo "The four-way figures need only stages 4-5, so they can be redrawn"
    echo "once $j5 finishes -- do not wait for $j7."
}

status() {
    printf "  %-22s %s\n" "cells linked"  "$(wc -l < "$CELLS/.cells" 2>/dev/null || echo 0)"
    printf "  %-22s %s\n" "step4 BAM"     "$(ls "$CELLS"/*_E99_Aligned.out.bam 2>/dev/null | wc -l)"
    printf "  %-22s %s\n" "step5 single"  "$(ls "$CELLS"/*singlemappers_genes.bed.gz 2>/dev/null | wc -l)"
    printf "  %-22s %s\n" "step5 multi"   "$(ls "$CELLS"/*multimappers_genes.bed.gz 2>/dev/null | wc -l)"
    printf "  %-22s %s\n" "step6 pickle"  "$([ -s "$OUT/$LIB.pickle.gz" ] && echo yes || echo no)"
    printf "  %-22s %s\n" "step7 tables"  "$(ls "$OUT/$LIB"*Counts.tsv 2>/dev/null | wc -l)"
    [ -s "$OUT/.jobs" ] && sacct -j "$(tr ' ' ',' < "$OUT/.jobs")" \
        --format=JobID%12,JobName%10,State%14,Elapsed,MaxRSS -X 2>/dev/null
}

case "${1:-status}" in
    setup)  setup ;;
    submit) submit ;;
    status) status ;;
    *) echo "usage: $0 [setup|submit|status]" >&2; exit 2 ;;
esac
