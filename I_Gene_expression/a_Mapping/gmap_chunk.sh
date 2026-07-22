#!/bin/bash
################################################################################
# gmap_chunk.sh
#
# PURPOSE
#   Batched replacement for map_star.sh (STAGE 4, gmap) when the pipeline is
#   run per-cell. It maps a CONTIGUOUS CHUNK of cells with a single STAR
#   genome load instead of one load per cell.
#
# WHY
#   VASA-plate demultiplexes to one fastq per cell (concatenator.py --demux),
#   so the driver's stage loop runs 384 times. The MIXED STAR index is 53 GB
#   (SA alone is 48 GB) while the median cell is only ~18 MB gz (~400k reads),
#   i.e. well under a minute of actual alignment. Loading the index once per
#   cell would therefore spend ~95% of the stage's wall time and ~20 TB of
#   Lustre reads on repeated index loads. Chunking 24 cells per job cuts that
#   to 16 loads.
#
# EQUIVALENCE TO map_star.sh
#   The STAR command line is byte-for-byte the one in map_star.sh with the
#   single addition of `--genomeLoad LoadAndKeep`. Alignment is per-read
#   independent, so the resulting BAMs are identical to the per-cell path.
#   Output names are unchanged --
#       ${folder}/${cell}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam
#   -- which matters because STAGE 6 (countTables_2pickle_cellsSpliced.py)
#   globs *.singlemappers_genes.bed.gz and recovers the cell id from the
#   filename prefix before '_cbc'.
#
# SHARED MEMORY
#   --genomeLoad LoadAndKeep puts the index in a SysV shared-memory segment
#   keyed on the genome directory, so two array tasks that happen to land on
#   the SAME node share one copy -- a bonus, except that whichever finishes
#   first would otherwise `Remove` the segment out from under the other. A
#   flock-guarded reference count in node-local /tmp makes the load/remove
#   pair safe: the first task on a node loads, the last one removes, and the
#   EXIT trap guarantees the decrement even if STAR fails or the job is
#   killed. Compute nodes have effectively unlimited shmmax/shmall (verified
#   on cn087), so a 56 GB segment is fine.
################################################################################

if [ $# -ne 7 ]
then
    echo "Please, give:"
    echo "1) manifest file (one cell basename per line)"
    echo "2) chunk id (1-based; use \$SLURM_ARRAY_TASK_ID)"
    echo "3) cells per chunk"
    echo "4) folder holding the per-cell fastqs"
    echo "5) STAR genome index dir"
    echo "6) path to STAR"
    echo "7) path to samtools"
    exit 1
fi

manifest=$1
chunk=$2
percell=$3
folder=$4
genome=$5
p2star=$6
p2samtools=$7

start=$(( (chunk - 1) * percell + 1 ))
end=$(( chunk * percell ))
mapfile -t cells < <(sed -n "${start},${end}p" "$manifest")

if [ ${#cells[@]} -eq 0 ]
then
    echo "[gmap_chunk] chunk ${chunk}: no cells in range ${start}-${end}, nothing to do"
    exit 0
fi

echo "[gmap_chunk] $(date) chunk ${chunk}: cells ${start}-$(( start + ${#cells[@]} - 1 )) (${#cells[@]} cells) on $(hostname)"

### ---- refcounted shared-memory genome load -------------------------------
# Key the lock on the genome path so unrelated indices never collide. /tmp is
# node-local, which is exactly the scope we need (one segment per node).
key=$(printf '%s' "$genome" | md5sum | cut -c1-12)
lockf=/tmp/vasa_star_shm.$(id -u).${key}.lock
cntf=/tmp/vasa_star_shm.$(id -u).${key}.cnt
tmpd=$(mktemp -d "${TMPDIR:-/tmp}/gmap_chunk.XXXXXX")

released=0
release() {
    [ "$released" -eq 1 ] && return
    released=1
    exec 9>>"$lockf" || return
    flock 9
    n=$(cat "$cntf" 2>/dev/null || echo 1)
    n=$(( n - 1 )); [ "$n" -lt 0 ] && n=0
    if [ "$n" -eq 0 ]
    then
        echo "[gmap_chunk] $(date) last user on $(hostname): removing shared-memory genome"
        "${p2star}"/STAR --genomeDir "$genome" --genomeLoad Remove \
            --outFileNamePrefix "${tmpd}/rm_" > /dev/null 2>&1
        rm -f "$cntf"
    else
        echo "$n" > "$cntf"
    fi
    flock -u 9
    exec 9>&-
    rm -rf "$tmpd"
}
trap release EXIT

exec 9>>"$lockf" || { echo "[gmap_chunk] cannot open lock $lockf"; exit 1; }
flock 9
n=$(cat "$cntf" 2>/dev/null || echo 0)
if [ "$n" -eq 0 ]
then
    echo "[gmap_chunk] $(date) loading genome into shared memory on $(hostname) ..."
    "${p2star}"/STAR --genomeDir "$genome" --genomeLoad LoadAndExit \
        --outFileNamePrefix "${tmpd}/load_"
    if [ $? -ne 0 ]
    then
        echo "[gmap_chunk] STAR LoadAndExit FAILED"
        flock -u 9; exit 1
    fi
    echo "[gmap_chunk] $(date) genome loaded"
fi
echo $(( n + 1 )) > "$cntf"
flock -u 9

### ---- map each cell in the chunk ------------------------------------------
# Command line identical to map_star.sh plus --genomeLoad LoadAndKeep.
rc=0
for cell in "${cells[@]}"
do
    inputfq=${folder}/${cell}_cbc_trimmed_homoATCG.nonRibo.fastq.gz
    outprefix=${folder}/${cell}_cbc_trimmed_homoATCG.nonRibo_E99_

    if [ ! -s "$inputfq" ]
    then
        echo "[gmap_chunk] SKIP ${cell}: missing or empty ${inputfq}"
        rc=1
        continue
    fi

    echo "[gmap_chunk] $(date) mapping ${cell}"
    "${p2star}"/STAR --runThreadN 8 --genomeDir "${genome}" \
        --genomeLoad LoadAndKeep \
        --readFilesIn "${inputfq}" --readFilesCommand zcat \
        --outFilterMultimapNmax 20 --outSAMunmapped Within \
        --outSAMtype BAM Unsorted --outSAMattributes All \
        --outFileNamePrefix "${outprefix}"
    if [ $? -ne 0 ]
    then
        echo "[gmap_chunk] STAR FAILED for ${cell}"
        rc=1
        continue
    fi

    # same tidy-up as map_star.sh
    rm -rf "${outprefix}_STARtmp"
    rm -f "${outprefix}Log.progress.out"
    mv "${outprefix}Log.out" "${outprefix}Log.txt"
    mv "${outprefix}Log.final.out" "${outprefix}Log.final.txt"
done

echo "[gmap_chunk] $(date) chunk ${chunk} done (rc=${rc})"
exit $rc
