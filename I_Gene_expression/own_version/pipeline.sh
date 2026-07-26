#!/bin/bash
###############################################################################
# pipeline.sh -- minimal, step-by-step VASA-plate mapping pipeline
#
# HOW TO USE
#   1. Edit config.sh
#   2. ./pipeline.sh check          <- verifies every path and tool. Do this first.
#   3. ./pipeline.sh step1          <- run one step at a time, look at the output
#      ./pipeline.sh step2          ...
#   or ./pipeline.sh all            <- run steps 1-7 back to back
#
# DEBUGGING HELPERS
#   MAXCELLS=4 ./pipeline.sh step2  <- only process the first 4 cells (fast!)
#   ./pipeline.sh status            <- how many files exist at each stage
#   ./pipeline.sh                   <- this help
#
#   Every step writes a log to $LOGDIR/<step>.log as well as the screen.
#   Every step is re-runnable: it overwrites its own outputs.
#
# THE SEVEN STEPS (each one's output is the next one's input)
#   step1 extract  fastq            -> cells/<sample>_<cell>_cbc.fastq.gz
#   step2 trim     _cbc             -> _cbc_trimmed_homoATCG.fq.gz
#   step3 ribo     trimmed          -> .nonRibo.fastq.gz      (rRNA removed)
#   step4 map      nonRibo          -> _E99_Aligned.out.bam   (STAR)
#   step5 assign   bam              -> *_genes.bed.gz         (reads -> genes)
#   step6 pickle   all cells' beds  -> <sample>.pickle.gz
#   step7 tables   pickle           -> <sample>_*.tsv         (final counts)
###############################################################################

set -u -o pipefail

# Locate this script's directory so config.sh (which lives next to it) can be
# sourced no matter the working directory. Normally BASH_SOURCE gives it. But
# `sbatch pipeline.sh` copies the script into a spool dir (/tmp/slurmd/...) and
# runs the copy, so BASH_SOURCE points there and config.sh is absent; in that
# case fall back to SLURM_SUBMIT_DIR, the directory sbatch was invoked from.
# (die() is not defined until below, so error out inline here.)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ ! -f "${HERE}/config.sh" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.sh" ]; then
    HERE="${SLURM_SUBMIT_DIR}"
fi
[ -f "${HERE}/config.sh" ] || { echo "ERROR: cannot find config.sh next to pipeline.sh (looked in ${HERE}). If you ran 'sbatch pipeline.sh', submit from the own_version/ directory so SLURM_SUBMIT_DIR points at it." >&2; exit 1; }
source "${HERE}/config.sh"

MAXCELLS="${MAXCELLS:-0}"     # 0 = all cells; >0 = only the first N (debugging)

# --- tiny helpers ------------------------------------------------------------
say()  { echo "[$(date '+%H:%M:%S')] $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }
rule() { echo "-------------------------------------------------------------"; }

# The list of cells we are working on, derived from step1's output.
# Written to a manifest so every later step agrees on the order and the set.
cell_list() {
    [ -s "${CELLDIR}/.cells" ] || die "no cell list yet -- run step1 first"
    if [ "$MAXCELLS" -gt 0 ]; then head -n "$MAXCELLS" "${CELLDIR}/.cells"
    else cat "${CELLDIR}/.cells"; fi
}

# Run a shell function once per cell, NCORES at a time.
# Usage: parallel_over_cells <function-name>
parallel_over_cells() {
    local fn=$1
    export -f "$fn" say die
    export CELLDIR OUTDIR LOGDIR SAMPLE VASA_SCRIPTS STRANDED REF_BED RRNA_FASTA
    export P2TRIMGALORE P2CUTADAPT P2BWA P2SAMTOOLS P2BEDTOOLS
    export TRIM_SH TRIM_MODE TRIM_ADAPTER3 TRIM_MINLEN TRIM_POLYA TRIM_POLYG TRIM_CUTADAPT
    export TRIM_ANCHOR_BC TRIM_LEN_UMI LEN_UMI TRIM_POLYT5 TRIM_ANCHOR_ADLEN
    cell_list | xargs -P "$NCORES" -I{} bash -c "$fn {}" _
}

###############################################################################
# check -- verify everything before you waste hours finding out otherwise
###############################################################################
step_check() {
    local bad=0
    rule; say "checking configuration"; rule

    echo "sample     : ${SAMPLE}"
    echo "fastq dir  : ${FASTQ_DIR}"
    echo "outdir     : ${OUTDIR}"
    echo "read setup : skip5=${SKIP5}  umi=${LEN_UMI}  cbc=${LEN_CBC}  umifirst=${UMI_FIRST}"
    echo

    # -- input fastqs --
    local r1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    local r2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
    for f in "$r1" "$r2"; do
        if [ -s "$f" ]; then echo "  OK   $f"
        else
            echo "  MISS $f"; bad=1
            # the most common cause: files named _1/_2 instead of _R1/_R2
            local alt="${f/_R1.fastq.gz/_1.fastq.gz}"; alt="${alt/_R2.fastq.gz/_2.fastq.gz}"
            if [ -s "$alt" ]; then
                echo "       but found  $alt"
                echo "       fix with:  ln -s $alt $f"
            fi
        fi
    done

    # -- references --
    for f in "$BC_WHITELIST" "$RRNA_FASTA" "$REF_BED"; do
        [ -s "$f" ] && echo "  OK   $f" || { echo "  MISS $f"; bad=1; }
    done
    # Test for SA, not just the directory: build_mouse_reference.sh creates the
    # directory up front, so a still-building or failed index would otherwise
    # pass check and only blow up hours later in step4.
    if [ -s "${STAR_INDEX}/SA" ] && [ -s "${STAR_INDEX}/Genome" ]; then
        echo "  OK   $STAR_INDEX"
    else
        echo "  MISS ${STAR_INDEX}/SA  (index absent or still building -- run build_mouse_reference.sh)"; bad=1
    fi
    # bwa index files must sit beside the rRNA fasta
    for ext in amb ann bwt pac sa; do
        [ -s "${RRNA_FASTA}.${ext}" ] || { echo "  MISS ${RRNA_FASTA}.${ext} (run: bwa index ${RRNA_FASTA})"; bad=1; }
    done

    # -- our own + borrowed scripts --
    for own in "$CONCATENATOR" "$TRIM_SH"; do
        [ -x "$own" ] && echo "  OK   $own" || { echo "  MISS/NOT-EXEC $own"; bad=1; }
    done
    # trim.sh is ours now, so it is not in this list
    for s in ribo-bwamem.sh riboread-selection.py deal_with_singlemappers.sh \
             deal_with_multimappers.sh countTables_2pickle_cellsSpliced.py countTables_fromPickle.py; do
        [ -x "${VASA_SCRIPTS}/$s" ] && echo "  OK   ${VASA_SCRIPTS}/$s" \
            || { echo "  MISS/NOT-EXEC ${VASA_SCRIPTS}/$s"; bad=1; }
    done

    # -- tools --
    for t in "${P2STAR}/STAR" "${P2BWA}/bwa" "${P2SAMTOOLS}/samtools" \
             "${P2BEDTOOLS}/bedtools" "${P2TRIMGALORE}/trim_galore" "${P2CUTADAPT}/cutadapt"; do
        [ -x "$t" ] && echo "  OK   $t" || { echo "  MISS $t"; bad=1; }
    done
    [ -x "${CONDA_ENV}/bin/python3" ] && echo "  OK   ${CONDA_ENV}/bin/python3" \
        || { echo "  MISS ${CONDA_ENV}/bin/python3"; bad=1; }
    # step 2 pass 2 needs cutadapt >= 4.4 for per-adapter `;min_overlap=`; the
    # module tree only has 1.18, hence the copy in the conda env.
    if [ -x "$TRIM_CUTADAPT" ]; then
        echo "  OK   $TRIM_CUTADAPT ($("$TRIM_CUTADAPT" --version 2>/dev/null))"
    else
        echo "  MISS $TRIM_CUTADAPT  (pip install cutadapt into ${CONDA_ENV})"; bad=1
    fi

    # -- read geometry sanity, using the actual first read --
    if [ -s "$r1" ] && [ -s "$r2" ]; then
        echo
        local l1 l2 need s1 s2 sbc sbio lbio
        # '2{p;q}' not '2p': sed must QUIT at line 2, otherwise it drains the
        # whole stream and zcat decompresses all 15 GB just to measure one read.
        l1=$(zcat "$r1" | sed -n '2{p;q}' | tr -d '\n' | wc -c)
        l2=$(zcat "$r2" | sed -n '2{p;q}' | tr -d '\n' | wc -c)
        # effective per-mate skip (per-mate override wins over SKIP5)
        s1="${SKIP5R1:-$SKIP5}"; s2="${SKIP5R2:-$SKIP5}"
        need=$(( LEN_UMI + LEN_CBC ))
        echo "  first read: R1=${l1} nt, R2=${l2} nt"
        echo "  5' skip   : R1=${s1} nt, R2=${s2} nt"
        # which mate holds the barcode, and which the cDNA
        if [ "$BC_READ" = "R1" ]; then sbc=$s1; lbio=$l2; sbio=$s2
        else                           sbc=$s2; lbio=$l1; sbio=$s1; fi
        local lbc; [ "$BC_READ" = "R1" ] && lbc=$l1 || lbc=$l2
        if [ "$lbc" -lt $(( sbc + need )) ]; then
            echo "  PROBLEM: barcode read ${BC_READ} is ${lbc} nt but needs >= $(( sbc + need ))"
            echo "           (skip ${sbc} + umi ${LEN_UMI} + cbc ${LEN_CBC})"
            echo "           either the skip is wrong, or the barcode is not behind the prefix."
            echo "           if only the cDNA mate carries the prefix, try: SKIP5R1=0 SKIP5R2=${SKIP5}"
            bad=1
        else
            echo "  geometry OK: after skipping ${sbc}, the barcode block still fits in ${BC_READ}"
        fi
        # Verify the index against the geometry rather than just printing the
        # number. The index is now shared with the nf-core FLASH-seq run, so its
        # overhang (150) is deliberately LARGER than this library needs (129) --
        # that is fine and must not be reported as a mismatch. Too small is the
        # only failure: STAR still maps and still writes a BAM, just worse,
        # silently. That is the one thing this check exists to catch.
        local need_ovh=$(( lbio - sbio - 1 ))
        local gp="${STAR_INDEX}/genomeParameters.txt"
        local have_ovh=""
        [ -r "$gp" ] && have_ovh=$(awk '$1=="sjdbOverhang"{print $2; exit}' "$gp")
        case "$have_ovh" in
            ''|*[!0-9]*)
                echo "  -> build/verify your STAR index with sjdbOverhang = ${need_ovh}"
                echo "     (could not read sjdbOverhang from ${gp})"
                ;;
            *)
                if [ "$have_ovh" -ge "$need_ovh" ]; then
                    echo "  STAR index sjdbOverhang = ${have_ovh} (needs >= ${need_ovh}) -- OK"
                else
                    echo "  PROBLEM: STAR index sjdbOverhang = ${have_ovh}, needs >= ${need_ovh}"
                    echo "           ${STAR_INDEX}"
                    echo "           rebuild it, or point STAR_INDEX at one built for >= ${need_ovh}."
                    bad=1
                fi
                ;;
        esac
    fi

    rule
    [ $bad -eq 0 ] && say "ALL CHECKS PASSED -- you can run ./pipeline.sh step1" \
                   || die "fix the MISS/PROBLEM lines above before running anything"
}

###############################################################################
# step1 -- extract barcodes and split into one fastq per cell
#
# Moves the cell barcode + UMI off the barcode read and onto the read NAME of
# the biological read, then writes one fastq per whitelist cell. Our forked
# concatenator.py also strips the ${SKIP5} nt prefix from both mates first.
#
# concatenator.py is single-core, so for speed we cut the input into shards,
# run one concatenator per shard in parallel, and glue the per-cell outputs
# back together. Set NCORES=1 in config.sh to take the simple one-process path
# instead -- much easier to follow when debugging.
###############################################################################

# sbatch -c 32 --mem=64G --time=8:00:00 pipeline.sh step1

step1_extract() {
    mkdir -p "$CELLDIR" "$LOGDIR"
    local r1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    local r2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
    [ -s "$r1" ] && [ -s "$r2" ] || die "input fastqs not found -- run ./pipeline.sh check"

    # arguments shared by both paths
    local args=(--cbcfile "$BC_WHITELIST" --cbchd "$BC_HAMMING"
                --lenumi "$LEN_UMI" --lencbc "$LEN_CBC"
                --bcread "$BC_READ" --bioread "$BIO_READ"
                --skip5 "$SKIP5" --demux)
    [ "$UMI_FIRST" = "yes" ] && args+=(--umifirst)
    # per-mate prefix overrides, only passed when actually set
    [ -n "${SKIP5R1:-}" ] && args+=(--skip5r1 "$SKIP5R1")
    [ -n "${SKIP5R2:-}" ] && args+=(--skip5r2 "$SKIP5R2")

    eval "$CONDA_ACTIVATE"

    if [ "$NCORES" -le 1 ]; then
        # ---------- simple path: one process ----------
        say "step1: single-process demultiplex (NCORES=1)"
        local work="${OUTDIR}/.step1_simple"; rm -rf "$work"; mkdir -p "$work"
        ln -sf "$r1" "${work}/${SAMPLE}_R1.fastq.gz"
        ln -sf "$r2" "${work}/${SAMPLE}_R2.fastq.gz"
        ( cd "$work" && "$CONCATENATOR" --fqf "$SAMPLE" "${args[@]}" --outdir "$CELLDIR" ) \
            || die "concatenator.py failed"
        mv -f "${CELLDIR}/${SAMPLE}.log" "${LOGDIR}/step1_extract.summary" 2>/dev/null
        rm -rf "$work"
    else
        # ---------- parallel path: shard, run, merge ----------
        local work="${OUTDIR}/.step1_work"; rm -rf "$work"; mkdir -p "$work/shards" "$work/out"
        local lines=$(( SHARD_READS * 4 ))

        say "step1a: splitting into shards of ${SHARD_READS} reads"
        # Split R1 and R2 the same way, so shard N of R1 pairs with shard N of R2.
        ( zcat "$r1" | split -l "$lines" -d -a 5 --additional-suffix=_R1.fastq \
              --filter="gzip -1 > ${work}/shards/\$FILE.gz" - sh ) &
        local p1=$!
        ( zcat "$r2" | split -l "$lines" -d -a 5 --additional-suffix=_R2.fastq \
              --filter="gzip -1 > ${work}/shards/\$FILE.gz" - sh ) &
        local p2=$!
        wait $p1 || die "splitting R1 failed"
        wait $p2 || die "splitting R2 failed"

        local shards; shards=$(cd "${work}/shards" && ls sh*_R1.fastq.gz | sed 's/_R1.fastq.gz$//' | sort)
        say "step1b: running $(echo "$shards" | wc -l) shards, ${NCORES} at a time"

        # one concatenator per shard, each writing into its own directory so no
        # two workers ever touch the same output file
        run_shard() {
            local pfx=$1
            mkdir -p "${WORKDIR}/out/${pfx}"
            ( cd "${WORKDIR}/shards" && \
              "$CONCATENATOR" --fqf "$pfx" "${CATARGS[@]}" --outdir "${WORKDIR}/out/${pfx}" \
            ) > "${WORKDIR}/out/${pfx}.log" 2>&1
        }
        export -f run_shard
        export WORKDIR="$work" CONCATENATOR
        export CATARGS_STR="${args[*]}"
        # rebuild the array inside the subshell (arrays don't survive export)
        echo "$shards" | xargs -P "$NCORES" -I{} bash -c \
            'read -r -a CATARGS <<< "$CATARGS_STR"; run_shard "$1"' _ {}

        say "step1c: merging per-cell files across shards"
        # Never glob every shard x cell at once -- that overflows the command
        # line. Walk shard dirs and let xargs batch the cat calls.
        local shdirs=(); mapfile -t shdirs < <(find "${work}/out" -mindepth 1 -maxdepth 1 -type d | sort)
        [ ${#shdirs[@]} -gt 0 ] || die "no shard output produced -- check ${work}/out/*.log"
        local cells; cells=$(ls "${shdirs[0]}"/*_cbc.fastq.gz | sed -E 's/.*_([0-9]+)_cbc\.fastq\.gz$/\1/' | sort -u)
        for cell in $cells; do
            for d in "${shdirs[@]}"; do printf '%s\n' "$d"/*_"${cell}"_cbc.fastq.gz; done \
                | xargs cat > "${CELLDIR}/${SAMPLE}_${cell}_cbc.fastq.gz"
        done

        # add up the per-shard logs into one summary
        awk -F',' '/total sequenced reads/{gsub(/ /,"",$2); t+=$2}
                   /reads with proper barcodes/{gsub(/ /,"",$2); k+=$2}
                   END{printf "total sequenced reads: %d\nreads with proper barcodes: %d (%.4f)\n", t, k, (t>0?k/t:0)}' \
            $(find "${work}/out" -name '*.log' | sort) > "${LOGDIR}/step1_extract.summary"
        rm -rf "$work"
    fi

    # the manifest every later step reads
    ls "${CELLDIR}/${SAMPLE}"_*_cbc.fastq.gz 2>/dev/null \
        | sed 's#.*/##; s#_cbc\.fastq\.gz$##' | sort > "${CELLDIR}/.cells"
    local n; n=$(wc -l < "${CELLDIR}/.cells")
    [ "$n" -gt 0 ] || die "step1 produced no cells -- check ${LOGDIR}/step1_extract.summary"
    say "step1 done: ${n} cells in ${CELLDIR}"
    cat "${LOGDIR}/step1_extract.summary"
}

###############################################################################
# step2 -- trim adapters/quality, then the 3' read-through
# TrimGalore first (adapters + Phred<20 from the 3' end), then cutadapt for the
# read-through construct this library actually carries:
#     [insert][poly-A][12 nt = revcomp(CBC+UMI)][revcomp(R1 5' prefix)+Nextera]
# This uses own_version/trim.sh, NOT ../a_Mapping/trim.sh -- see that file's
# header and README.md "Step 2" for why, and TRIM_MODE=legacy in config.sh to
# get the upstream behaviour back.
###############################################################################
do_trim() {
    local cell=$1
    # the TRIM_* settings reach trim.sh through parallel_over_cells' export list
    "$TRIM_SH" "${CELLDIR}/${cell}_cbc.fastq.gz" "$CELLDIR" \
        "$P2TRIMGALORE" "$P2CUTADAPT" > "${LOGDIR}/step2_${cell}.log" 2>&1 \
        || echo "  FAILED trim: $cell (see ${LOGDIR}/step2_${cell}.log)"
}
step2_trim() {
    say "step2: trimming $(cell_list | wc -l) cells, ${NCORES} at a time (TRIM_MODE=${TRIM_MODE})"
    eval "$ML_TRIM"
    parallel_over_cells do_trim
    say "step2 done: $(ls ${CELLDIR}/*_trimmed_homoATCG.fq.gz 2>/dev/null | wc -l) trimmed files"
}

###############################################################################
# step3 -- in-silico rRNA depletion
# Maps against the rRNA reference with BOTH `bwa aln` and `bwa mem` and keeps
# only reads that neither aligned. rRNA fragments span too wide a length range
# for either aligner alone to catch them all.
###############################################################################
do_ribo() {
    local cell=$1
    "${VASA_SCRIPTS}/ribo-bwamem.sh" "$RRNA_FASTA" \
        "${CELLDIR}/${cell}_cbc_trimmed_homoATCG.fq.gz" \
        "${CELLDIR}/${cell}_cbc_trimmed_homoATCG" \
        "$P2BWA" "$P2SAMTOOLS" "$STRANDED" "$VASA_SCRIPTS" > /dev/null 2>&1 \
        || echo "  FAILED ribo: $cell"
}
step3_ribo() {
    say "step3: rRNA depletion on $(cell_list | wc -l) cells"
    eval "$ML_RIBO"
    parallel_over_cells do_ribo
    say "step3 done: $(ls ${CELLDIR}/*.nonRibo.fastq.gz 2>/dev/null | wc -l) nonRibo files"
}

###############################################################################
# step4 -- map to the genome with STAR
#
# The index is tens of GB and loading it takes far longer than mapping one
# small cell. So we load it into shared memory ONCE, map every cell against
# the loaded copy, then free it. The trap guarantees the memory is released
# even if you Ctrl-C or a cell fails -- otherwise it would sit there until
# the node reboots.
###############################################################################
step4_map() {
    say "step4: STAR mapping $(cell_list | wc -l) cells"
    eval "$ML_STAR"
    local tmp; tmp=$(mktemp -d "${OUTDIR}/.star.XXXXXX")

    release_genome() {
        say "step4: releasing genome from shared memory"
        "${P2STAR}/STAR" --genomeDir "$STAR_INDEX" --genomeLoad Remove \
            --outFileNamePrefix "${tmp}/rm_" > /dev/null 2>&1
        rm -rf "$tmp"
    }
    trap release_genome EXIT

    say "step4a: loading genome into shared memory (slow, happens once)"
    "${P2STAR}/STAR" --genomeDir "$STAR_INDEX" --genomeLoad LoadAndExit \
        --outFileNamePrefix "${tmp}/load_" > /dev/null || die "STAR could not load the index"

    local i=0 n; n=$(cell_list | wc -l)
    while read -r cell; do
        i=$((i+1))
        local fq="${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo.fastq.gz"
        local pfx="${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo_E99_"
        if [ ! -s "$fq" ]; then echo "  SKIP (no input): $cell"; continue; fi
        printf "\r  mapping %d/%d  %s   " "$i" "$n" "$cell"
        "${P2STAR}/STAR" --runThreadN "$STAR_THREADS" --genomeDir "$STAR_INDEX" \
            --genomeLoad LoadAndKeep \
            --readFilesIn "$fq" --readFilesCommand zcat \
            --outFilterMultimapNmax 20 --outSAMunmapped Within \
            --outSAMtype BAM Unsorted --outSAMattributes All \
            --outFileNamePrefix "$pfx" > /dev/null 2>&1 \
            || echo "  FAILED star: $cell"
        rm -rf "${pfx}_STARtmp" "${pfx}Log.progress.out"
        mv -f "${pfx}Log.final.out" "${pfx}Log.final.txt" 2>/dev/null
        mv -f "${pfx}Log.out"       "${pfx}Log.txt"       2>/dev/null
    done < <(cell_list)
    echo
    say "step4 done: $(ls ${CELLDIR}/*_E99_Aligned.out.bam 2>/dev/null | wc -l) BAMs"
    # trap fires here and frees the shared memory
}

###############################################################################
# step5 -- assign reads to genes
# Two passes over the same BAM: uniquely-mapping reads (NH:i:1) and
# multi-mapping reads (NH:i:2-9). Both intersect the annotation BED and write
# *_genes.bed.gz. Multimappers are RESCUED, not discarded.
###############################################################################
do_assign() {
    local cell=$1
    local bam="${CELLDIR}/${cell}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam"
    [ -s "$bam" ] || { echo "  SKIP (no bam): $cell"; return; }
    "${VASA_SCRIPTS}/deal_with_singlemappers.sh" "$bam" "$REF_BED" "$STRANDED" \
        "$P2SAMTOOLS" "$P2BEDTOOLS" > /dev/null 2>&1 || echo "  FAILED single: $cell"
    "${VASA_SCRIPTS}/deal_with_multimappers.sh"  "$bam" "$REF_BED" "$STRANDED" \
        "$P2SAMTOOLS" "$P2BEDTOOLS" > /dev/null 2>&1 || echo "  FAILED multi: $cell"
}
step5_assign() {
    say "step5: assigning reads to genes for $(cell_list | wc -l) cells"
    eval "$ML_BED"
    parallel_over_cells do_assign
    say "step5 done: $(ls ${CELLDIR}/*singlemappers_genes.bed.gz 2>/dev/null | wc -l) single, $(ls ${CELLDIR}/*multimappers_genes.bed.gz 2>/dev/null | wc -l) multi"
}

###############################################################################
# step6 -- collapse every cell into one UMI-aware structure
# Runs ONCE over the whole folder (not per cell), so steps 2-5 must be finished
# for every cell first.
###############################################################################
step6_pickle() {
    say "step6: building the count pickle (this is the memory-hungry step)"
    eval "$CONDA_ACTIVATE"
    # run from OUTDIR with a relative folder name: the cell id is derived from
    # the path, so an absolute path would end up inside every column name
    ( cd "$OUTDIR" && "${VASA_SCRIPTS}/countTables_2pickle_cellsSpliced.py" \
        "cells" "$SAMPLE" vasa "$CELLID_FROM" ) || die "step6 failed"
    say "step6 done: ${OUTDIR}/${SAMPLE}.pickle.gz"
}

###############################################################################
# step7 -- final count tables
# Produces spliced / unspliced / total / tRNA tables as Read, UMI (UFI) and
# Transcript counts.
###############################################################################
step7_tables() {
    say "step7: writing count tables"
    eval "$CONDA_ACTIVATE"
    ( cd "$OUTDIR" && "${VASA_SCRIPTS}/countTables_fromPickle.py" \
        "${SAMPLE}.pickle.gz" "$SAMPLE" vasa y ) || die "step7 failed"
    say "step7 done. Tables:"
    ls -1 "${OUTDIR}/${SAMPLE}"*.tsv 2>/dev/null | sed 's/^/    /'
}

###############################################################################
# status / all / dispatcher
###############################################################################
step_status() {
    rule; say "pipeline status for ${SAMPLE}"; rule
    local n=0
    [ -s "${CELLDIR}/.cells" ] && n=$(wc -l < "${CELLDIR}/.cells")
    printf "  %-34s %s\n" "cells from step1"        "$n"
    printf "  %-34s %s\n" "step2 trimmed"           "$(ls ${CELLDIR}/*_trimmed_homoATCG.fq.gz 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "step3 nonRibo"           "$(ls ${CELLDIR}/*.nonRibo.fastq.gz 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "step4 BAM"               "$(ls ${CELLDIR}/*_E99_Aligned.out.bam 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "step5 singlemappers bed" "$(ls ${CELLDIR}/*singlemappers_genes.bed.gz 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "step5 multimappers bed"  "$(ls ${CELLDIR}/*multimappers_genes.bed.gz 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "step6 pickle"            "$([ -s ${OUTDIR}/${SAMPLE}.pickle.gz ] && echo yes || echo no)"
    printf "  %-34s %s\n" "step7 tables"            "$(ls ${OUTDIR}/${SAMPLE}*.tsv 2>/dev/null | wc -l)"
    rule
}

step_all() {
    step1_extract; step2_trim; step3_ribo; step4_map; step5_assign; step6_pickle; step7_tables
    say "ALL STEPS COMPLETE"
}

usage() {
    sed -n '2,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 0
}

mkdir -p "$LOGDIR" 2>/dev/null || true
cmd="${1:-help}"
case "$cmd" in
    check)  step_check ;;
    step1|extract) step1_extract 2>&1 | tee "${LOGDIR}/step1.log" ;;
    step2|trim)    step2_trim    2>&1 | tee "${LOGDIR}/step2.log" ;;
    step3|ribo)    step3_ribo    2>&1 | tee "${LOGDIR}/step3.log" ;;
    step4|map)     step4_map     2>&1 | tee "${LOGDIR}/step4.log" ;;
    step5|assign)  step5_assign  2>&1 | tee "${LOGDIR}/step5.log" ;;
    step6|pickle)  step6_pickle  2>&1 | tee "${LOGDIR}/step6.log" ;;
    step7|tables)  step7_tables  2>&1 | tee "${LOGDIR}/step7.log" ;;
    all)    step_all 2>&1 | tee "${LOGDIR}/all.log" ;;
    status) step_status ;;
    help|-h|--help) usage ;;
    *) echo "unknown command: $cmd"; echo; usage ;;
esac
