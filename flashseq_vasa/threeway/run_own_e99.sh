#!/bin/bash
###############################################################################
# run_own_e99.sh -- re-quantify the own VASA plate (ZHA9292A1) under the
#                   PUBLISHED plate's reference: Ensembl 99, GRCh38+GRCm38.
#
# WHY
# ---
# The published plate is GRCm38/Ensembl 99 (human+mouse); the own plate is
# GRCm39/Ensembl 116 (mouse only). Every published-vs-own gap therefore mixes
# protocol biology with annotation release. The shared-gene-universe control in
# res/threeway/threeway_release_control.tsv is only partial: it restricts to
# genes present in both releases but still quantifies each dataset under its OWN
# annotation, gene models and genome build. This script removes the release term
# outright by running the own plate through the published reference.
#
# HOW -- and what is deliberately NOT re-done
# -------------------------------------------
# Stages 1-3 (barcode extraction, trimming, rRNA depletion) are reference-free
# with ONE exception, and are reused as-is from the existing run:
#   - stage 3 depletes against unique_rRNA_mouse.v2.fa (mouse 47S). The mixed
#     branch uses unique_rRNA_human_mouse.v3.fa. Both carry the same mouse
#     BK000964.3 47S unit, so for a mouse-only library the depleted output is
#     equivalent; the human sequences can only remove reads that are not there.
#     This is asserted, not assumed: see ASSERTION 3 below.
# Stages 4-7 are re-run:
#   4 STAR   -> the new Ensembl 99 index at sjdbOverhang 150
#   5 assign -> the Ensembl 99 mixed BED
#   6 pickle -> unchanged upstream code
#   7 tables -> unchanged upstream code
#
# ISOLATION (this is the part that must not go wrong)
# --------------------------------------------------
# The existing E116 outputs are a VALIDATED result and are read-only here. Three
# separate mechanisms keep them that way:
#   1. OUT points at a NEW tree, out_E99/, never at out/.
#   2. out_E99/cells holds SYMLINKS to the step-3 fastqs, so stages 1-3 cannot
#      be re-run by accident and the 30 GB of fastq is not duplicated.
#   3. OUT is asserted not to overlap SRC before anything is written.
# NB the existing files are named *_E99_* -- that label is a MISNOMER: their
# STAR Log.txt records genomeDir=GRCm39/star_index_151_r116, i.e. Ensembl 116.
# This script never writes into that directory, so the misnomer is left alone
# rather than renamed underneath a result that is already in use.
#
# Fork discipline (Rule 1): a_Mapping/ is called, never modified. Step 5 uses
# the SAME upstream scripts the published plate used, not own_version's fork --
# see ASSERTION 4.
#
# Usage: run_own_e99.sh <prepare|check|step4|step5|step6|step7>
###############################################################################
set -euo pipefail

W=/nemo/lab/turnerj/working/guangxin
REF=${W}/reference/vasaseq
SRC=${W}/vasaseq/data/PM26037/out                 # existing E116 run, READ ONLY
OUT=${W}/vasaseq/data/PM26037/out_E99             # this run
CELLDIR=${OUT}/cells
SAMPLE=ZHA9292A1
SCRATCH=/nemo/lab/turnerj/scratch/zhangg/vasaseq

STAR_INDEX=${REF}/mixed/star_index_151
REF_BED=${REF}/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed
STRANDED=y                                        # Rule 6: VASA is stranded
VASA_SCRIPTS=${W}/vasaseq/code/I_Gene_expression/a_Mapping
CONDA_ENV=${W}/envs/vasa
NCORES=16
STAR_THREADS=16

EBROOT=/camp/apps/eb/software
P2STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin
P2SAMTOOLS=${EBROOT}/SAMtools/1.11-GCC-10.2.0/bin
P2BEDTOOLS=${EBROOT}/BEDTools/2.30.0-GCC-11.2.0/bin
ML_INIT="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=${EBROOT%/software}/modules/all"
ML_STAR="${ML_INIT}; module load STAR/2.7.7a-GCC-10.2.0 SAMtools/1.11-GCC-10.2.0"
ML_BED="${ML_INIT}; module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0"
CONDA_ACTIVATE="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate ${CONDA_ENV}"

# Suffix for stage-4+ products. Deliberately NOT "_E99_": that string is already
# on the E116 run's filenames and reusing it invites exactly the confusion this
# script exists to resolve.
TAG=_e99mixed_

say() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
rm_stale() { [ $# -gt 0 ] && rm -f "$@"; return 0; }

cell_list() { cat "${CELLDIR}/.cells"; }

# ---------------------------------------------------------------------------
# prepare: build out_E99/cells as symlinks to the step-3 products
# ---------------------------------------------------------------------------
prepare() {
    say "prepare: ${OUT}"
    [ -d "$SRC" ] || die "source run missing: $SRC"
    # ASSERTION 1: never write into the existing validated run.
    case "$OUT" in
        "$SRC"|"$SRC"/*) die "OUT would overlap SRC -- refusing" ;;
    esac
    mkdir -p "$CELLDIR"
    cp -f "${SRC}/cells/.cells" "${CELLDIR}/.cells"
    local n=0 cell f
    while read -r cell; do
        f=${cell}_cbc_trimmed_homoATCG.nonRibo.fastq.gz
        [ -s "${SRC}/cells/${f}" ] || die "missing step-3 output: ${f}"
        ln -sf "${SRC}/cells/${f}" "${CELLDIR}/${f}"
        n=$((n+1))
    done < "${CELLDIR}/.cells"
    say "prepare: linked ${n} step-3 fastqs (no copy, no re-trim)"
}

# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------
step_check() {
    local bad=0
    say "check: reference"
    local f
    for f in "${STAR_INDEX}/SA" "${STAR_INDEX}/Genome" "${STAR_INDEX}/SAindex" "$REF_BED"; do
        if [ -s "$f" ]; then echo "  ok   $f"; else echo "  MISS $f"; bad=1; fi
    done
    # ASSERTION 2: the index must be the Ensembl 99 MIXED one, at an overhang
    # large enough for 130 nt reads. Both, not either.
    local ovh
    ovh=$(awk '$1=="sjdbOverhang"{print $2; exit}' "${STAR_INDEX}/genomeParameters.txt" 2>/dev/null || echo "")
    echo "  index sjdbOverhang: ${ovh:-unreadable}"
    if [ -n "$ovh" ] && [ "$ovh" -ge 129 ]; then
        echo "  ok   overhang >= 129 (151 nt raw - 21 nt skip - 1)"
    else
        echo "  PROBLEM: need sjdbOverhang >= 129"; bad=1
    fi
    local nm nh
    nm=$(grep -c '^GRCm38_' "${STAR_INDEX}/chrName.txt" 2>/dev/null || true)
    nh=$(grep -c '^GRCh38_' "${STAR_INDEX}/chrName.txt" 2>/dev/null || true)
    echo "  index contigs: GRCm38_ ${nm}, GRCh38_ ${nh}"
    if [ "${nm:-0}" -gt 0 ] && [ "${nh:-0}" -gt 0 ]; then
        echo "  ok   index is the species-prefixed mixed build the BED expects"
    else
        echo "  PROBLEM: index is not the mixed one"; bad=1
    fi

    # ASSERTION 3: stage-3 equivalence. The mouse 47S unit that does the actual
    # depleting must be byte-identical between the mouse-only rRNA reference the
    # existing depletion used and the mixed one the published plate used --
    # otherwise reusing stage 3 silently changes the denominator.
    say "check: stage-3 rRNA reference equivalence (mouse 47S unit)"
    local mo=${REF}/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa
    local mx=${REF}/mixed/unique_rRNA_human_mouse.v3.fa
    # The 47S record is named `>mouse_rDNA_47S_BK000964.3_1-13403` in BOTH files
    # (verified 2026-07-30: line 1 of the mouse-only ref, line 967 of the mixed
    # one). Matching on `^>BK000964` instead extracts NOTHING from either file
    # and then compares two empty strings, whose md5 is
    # d41d8cd98f00b204e9800998ecf8427e -- a check that always passes and proves
    # nothing. Hence the emptiness guard below: an extraction that yields no
    # sequence is a FAILED check, not a passed one.
    if [ -s "$mo" ] && [ -s "$mx" ]; then
        local a b na nb
        a=$(awk '/^>mouse_rDNA_47S_BK000964/{p=1;next} /^>/{p=0} p' "$mo" | tr -d '\n')
        b=$(awk '/^>mouse_rDNA_47S_BK000964/{p=1;next} /^>/{p=0} p' "$mx" | tr -d '\n')
        na=${#a}; nb=${#b}
        echo "  47S extracted: mouse-only ${na} nt, mixed ${nb} nt (expect 13403)"
        if [ "$na" -eq 0 ] || [ "$nb" -eq 0 ]; then
            echo "  PROBLEM: extracted no 47S sequence -- the header pattern is wrong,"
            echo "           so this check would be vacuous. Refusing to pass it."
            bad=1
        elif [ "$na" -ne 13403 ] || [ "$nb" -ne 13403 ]; then
            echo "  PROBLEM: 47S length is not 13403 nt on both sides"; bad=1
        else
            a=$(printf '%s' "$a" | md5sum | cut -d' ' -f1)
            b=$(printf '%s' "$b" | md5sum | cut -d' ' -f1)
            echo "  47S md5 mouse-only = ${a}"
            echo "  47S md5 mixed      = ${b}"
            if [ "$a" = "$b" ]; then
                echo "  ok   the 47S unit doing the depletion is identical -> stage 3 reusable"
            else
                echo "  PROBLEM: 47S units differ -- stage 3 is NOT reusable"; bad=1
            fi
        fi
    else
        echo "  MISS one of the rRNA references"; bad=1
    fi

    # ASSERTION 4: step-5 code identity with the published run. The published
    # plate used a_Mapping/ (data/ref/fastq_vasaplate/rerun_cout_pick.sh:25,
    # p2s=.../a_Mapping). own_version/ forked those two scripts for the
    # NH:i:10-19 selection bug. Using the fork here would make the own-E99 arm
    # differ from the published arm by BOTH annotation and step-5 code, which
    # would defeat the control this script exists to provide.
    say "check: step-5 code provenance (must match the published run)"
    local s
    for s in deal_with_singlemappers.sh deal_with_multimappers.sh; do
        if [ -s "${VASA_SCRIPTS}/${s}" ]; then
            echo "  ok   upstream ${s}"
        else
            echo "  MISS ${VASA_SCRIPTS}/${s}"; bad=1
        fi
    done

    say "check: inputs"
    local n=0 miss=0 cell
    while read -r cell; do
        n=$((n+1))
        [ -s "${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo.fastq.gz" ] || miss=$((miss+1))
    done < "${CELLDIR}/.cells"
    echo "  cells ${n}, missing step-3 input ${miss}"
    [ "$miss" -eq 0 ] || bad=1

    if [ "$bad" -eq 0 ]; then say "CHECK PASSED"; else die "fix the above first"; fi
}

# ---------------------------------------------------------------------------
# step4 -- STAR against the Ensembl 99 mixed index
# ---------------------------------------------------------------------------
release_genome() {
    say "step4: releasing genome from shared memory"
    "${P2STAR}/STAR" --genomeDir "$STAR_INDEX" --genomeLoad Remove \
        --outFileNamePrefix "${STAR_TMP}/rm_" > /dev/null 2>&1 || true
    rm -rf "$STAR_TMP"
}
step4_map() {
    local n; n=$(wc -l < "${CELLDIR}/.cells")
    say "step4: STAR mapping ${n} cells against ${STAR_INDEX}"
    eval "$ML_STAR"
    mkdir -p "$SCRATCH"
    STAR_TMP=$(mktemp -d "${SCRATCH}/star_e99.XXXXXX")
    trap release_genome EXIT
    say "step4a: loading genome into shared memory (once)"
    "${P2STAR}/STAR" --genomeDir "$STAR_INDEX" --genomeLoad LoadAndExit \
        --outFileNamePrefix "${STAR_TMP}/load_" > /dev/null || die "STAR could not load the index"
    local i=0 cell fq pfx
    while read -r cell; do
        i=$((i+1))
        fq="${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo.fastq.gz"
        pfx="${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo${TAG}"
        if [ ! -s "$fq" ]; then echo "  SKIP (no input): $cell"; continue; fi
        printf "  mapping %d/%d  %s\n" "$i" "$n" "$cell"
        # STAR parameters identical to both existing runs; only genomeDir differs.
        "${P2STAR}/STAR" --runThreadN "$STAR_THREADS" --genomeDir "$STAR_INDEX" \
            --genomeLoad LoadAndKeep \
            --readFilesIn "$fq" --readFilesCommand zcat \
            --outFilterMultimapNmax 20 --outSAMunmapped Within \
            --outSAMtype BAM Unsorted --outSAMattributes All \
            --outFileNamePrefix "$pfx" > /dev/null 2>&1 \
            || echo "  FAILED star: $cell"
        rm -rf "${pfx}_STARtmp" "${pfx}Log.progress.out"
        mv -f "${pfx}Log.final.out" "${pfx}Log.final.txt" 2>/dev/null || true
        mv -f "${pfx}Log.out"       "${pfx}Log.txt"       2>/dev/null || true
    done < "${CELLDIR}/.cells"
    local got
    got=$(find "$CELLDIR" -name "*${TAG}Aligned.out.bam" | wc -l)
    say "step4 done: ${got} BAMs"
    [ "$got" -eq "$n" ] || die "expected ${n} BAMs, got ${got}"
}

# ---------------------------------------------------------------------------
# step5 -- assign against the Ensembl 99 BED, with the UPSTREAM scripts
# ---------------------------------------------------------------------------
do_assign() {
    local cell=$1
    local stem="${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo${TAG}Aligned.out"
    local bam="${stem}.bam"
    if [ ! -s "$bam" ]; then echo "  SKIP (no bam): $cell"; return 0; fi
    # Both deal_with_*mappers.sh end with a bare `gzip`; clear their targets
    # first, or a stale file is silently folded into the count tables by step 6.
    rm_stale "${stem}.singlemappers_genes.bed.gz" "${stem}.singlemappers_genes.bed" \
             "${stem}.nsorted.multimappers_genes.bed.gz" "${stem}.nsorted.multimappers_genes.bed"
    "${VASA_SCRIPTS}/deal_with_singlemappers.sh" "$bam" "$REF_BED" "$STRANDED" \
        "$P2SAMTOOLS" "$P2BEDTOOLS" > /dev/null 2>&1 || echo "  FAILED single: $cell"
    "${VASA_SCRIPTS}/deal_with_multimappers.sh"  "$bam" "$REF_BED" "$STRANDED" \
        "$P2SAMTOOLS" "$P2BEDTOOLS" > /dev/null 2>&1 || echo "  FAILED multi: $cell"
}
step5_assign() {
    say "step5: assigning against ${REF_BED}"
    eval "$ML_BED"
    # rm_stale must be exported too: workers run in a fresh `bash -c`, which
    # inherits neither functions nor unexported variables.
    export -f do_assign say die rm_stale
    export CELLDIR REF_BED STRANDED P2SAMTOOLS P2BEDTOOLS VASA_SCRIPTS TAG
    xargs -a "${CELLDIR}/.cells" -P "$NCORES" -I{} bash -c 'do_assign "$@"' _ {}
    local bad5=0 cell stem out
    while read -r cell; do
        stem="${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo${TAG}Aligned.out"
        for out in "${stem}.singlemappers_genes.bed.gz" "${stem}.nsorted.multimappers_genes.bed.gz"; do
            if [ ! -s "$out" ]; then
                echo "  BAD step5: $cell -- missing $(basename "$out")"; bad5=$((bad5+1))
            elif [ "$out" -ot "${stem}.bam" ]; then
                echo "  BAD step5: $cell -- $(basename "$out") is OLDER than its BAM"; bad5=$((bad5+1))
            fi
        done
    done < "${CELLDIR}/.cells"
    [ "$bad5" -eq 0 ] || die "step5: ${bad5} missing/stale outputs -- do NOT run step6"
    say "step5 done: $(find "$CELLDIR" -name "*singlemappers_genes.bed.gz" | wc -l) single, $(find "$CELLDIR" -name "*multimappers_genes.bed.gz" | wc -l) multi (all newer than their BAMs)"
}

# ---------------------------------------------------------------------------
# step6 / step7 -- unchanged upstream code
# ---------------------------------------------------------------------------
step6_pickle() {
    say "step6: building the count pickle (memory-hungry stage)"
    eval "$CONDA_ACTIVATE"
    rm_stale "${OUT}/${SAMPLE}.pickle.gz" "${OUT}/${SAMPLE}.pickle"
    # Run from OUT with a RELATIVE folder name: the cell id is derived from the
    # path, so an absolute path would end up inside every column name.
    ( cd "$OUT" && "${VASA_SCRIPTS}/countTables_2pickle_cellsSpliced.py" \
        "cells" "$SAMPLE" vasa r ) || die "step6 failed"
    say "step6 done: ${OUT}/${SAMPLE}.pickle.gz"
}
step7_tables() {
    say "step7: writing count tables"
    eval "$CONDA_ACTIVATE"
    ( cd "$OUT" && "${VASA_SCRIPTS}/countTables_fromPickle.py" \
        "${SAMPLE}.pickle.gz" "$SAMPLE" vasa y ) || die "step7 failed"
    say "step7 done:"
    find "$OUT" -maxdepth 1 -name "${SAMPLE}*.tsv" -printf "    %f\n" | sort
}

case "${1:-help}" in
    prepare) prepare ;;
    check)   prepare; step_check ;;
    step4)   step4_map ;;
    step5)   step5_assign ;;
    step6)   step6_pickle ;;
    step7)   step7_tables ;;
    *) echo "usage: $0 <prepare|check|step4|step5|step6|step7>"; exit 1 ;;
esac
