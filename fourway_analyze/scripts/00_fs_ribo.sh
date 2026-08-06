#!/bin/bash
###############################################################################
# 00_fs_ribo.sh -- the stage FLASH-seq was missing, so the four-way is a
# four-way.
#
# WHY THIS EXISTS
# ---------------
# `flashseq_vasa/pipeline_fs.sh` goes prep -> map -> assign. There is no rRNA
# depletion in it, so FLASH-seq's STAR input contains its rRNA and VASA's does
# not -- VASA's step 3 removes 26.1% of its reads before STAR ever sees them.
#
# That is harmless for a count table (rRNA has no gene to be assigned to) and
# fatal for everything this folder measures. rRNA is transcribed off a repeat
# array, so rRNA reads multimap by construction; leaving them in one library's
# STAR input and not the other's puts a protocol difference and a pipeline
# difference on the same axis. The demo_analyze step-4/step-5 conclusions --
# multimapping rate, biotype composition -- are exactly what that would corrupt.
#
# So this runs VASA's own step 3 on FLASH-seq's trimmed reads, with the same
# script, the same rRNA reference and the same two aligners, and the map/assign
# stages then run over the depleted fastq.
#
#   flashseq_vasa/pipeline_fs.sh prep            (native arm, unchanged)
#     -> run/native/cells/<LIB>_cbc_noumi_R1.fq.gz
#   THIS SCRIPT
#     -> run/nonribo/cells/<LIB>_cbc_noumi.Ribo.bam
#                          <LIB>_cbc_noumi.nsorted.all-ribo.bam
#                          <LIB>_cbc_noumi.ribo-map.log
#                          <LIB>_cbc_noumi.nonRibo.fastq.gz
#     -> run/nonribo/cells/<LIB>_cbc_noumi_R1.fq.gz  = SYMLINK to the above
#   flashseq_vasa/pipeline_fs.sh map assign  with FSV_OUTDIR=run/nonribo
#
# THE SYMLINK IS DELIBERATE AND IS THE HONEST FORM
# ------------------------------------------------
# `prep_fq()` in pipeline_fs.sh hard-codes `<stem>_R1.fq.gz`, and pipeline_fs.sh
# is not edited (it is the FLASH-seq side's own fork and this folder does not
# own it). A copy under that name would be a silent substitution -- the file
# would claim to be prep's output and nothing on disk would say otherwise. A
# symlink says exactly what it is in `ls -l`, costs no disk, and cannot drift
# from its target.
#
# STRANDED = n, AND THAT IS NOT THE SAME FLAG VASA USED
# -----------------------------------------------------
# FLASH-seq is genuinely unstranded -- 49.1-50.5% of its ribosomal reads are on
# the forward strand across all ten libraries (flashseq_vasa/README.md), so `n`
# is the correct flag and it is the one pipeline_fs.sh already uses at assign.
# VASA ran `y`, and VASA is measurably NOT perfectly stranded either (76.1%
# forward), so `y` discards ~24% of its ribosomal reads where the same flag
# would discard ~50% of FLASH-seq's.
#
# **The two in-silico depletion percentages are therefore not on one footing,
# and the four-way README says so wherever the number appears.** The flag is not
# split the other way here because the fastq that feeds STAR can only be made
# once, and `n` is the one that is right for this protocol. `.nsorted.all-ribo.bam`
# is kept so the `y` count stays recoverable from the same alignment without
# re-running bwa.
#
# COST (measured, job 51307376: 9 of 10 libraries -- A1 failed on stale files)
#   1h38m wall at NPAR=4 on -c 32 --mem=48G, peak 15.6 GB.
#   bwa aln is 350-480 s per library and is the floor; the rest is samtools.
#   MEMORY IS NOT THE CONSTRAINT -- 15.6 GB is the whole cgroup (this cluster
#   runs jobacct_gather/cgroup, so MaxRSS is the job total, not one process),
#   i.e. 3x headroom at 48G. Do not raise --mem hoping to go faster.
#
# Resumable: a library whose output already passes verify_lib is skipped, so a
# rerun after an interruption costs only what was unfinished. FORCE=1 redoes all.
#
#   sbatch -c 32 --mem=48G -t 8:00:00 --wrap="NPAR=4 scripts/00_fs_ribo.sh"
###############################################################################
set -euo pipefail
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)

VASA=/nemo/lab/turnerj/working/guangxin/vasaseq
REF=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116
RRNA=$REF/unique_rRNA_mouse.v2.fa
SCRIPTS=$VASA/code/I_Gene_expression/a_Mapping        # ribo-bwamem.sh + riboread-selection.py
RUN=$VASA/data/flashseq_vasa/run
SRC=$RUN/native/cells
DST=$RUN/nonribo/cells
LOG=$RUN/nonribo/logs
STRANDED=n
NPAR=${NPAR:-4}     # ribo-bwamem.sh is internally threaded (bwa mem -t 8, two
                    # samtools --threads 8), so 4 concurrent libraries already
                    # oversubscribes 16 cores; the stage is gzip-bound anyway.

EBROOT=/camp/apps/eb/software
P2BWA=$EBROOT/BWA/0.7.17-GCC-10.3.0/bin
P2SAMTOOLS=$EBROOT/SAMtools/1.11-GCC-10.2.0/bin

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=$EBROOT/../modules/all
module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0 2>/dev/null
# riboread-selection.py needs pysam, which is only in the conda env; the binaries
# above are called by absolute path so the env leading on PATH is harmless.
source $EBROOT/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa

for f in "$RRNA" "$RRNA.bwt" "$SCRIPTS/ribo-bwamem.sh" "$SCRIPTS/riboread-selection.py" \
         "$P2BWA/bwa" "$P2SAMTOOLS/samtools"; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done
mkdir -p "$DST" "$LOG"

libs=$(ls "$SRC"/*_cbc_noumi_R1.fq.gz 2>/dev/null | xargs -n1 basename | sed 's/_cbc_noumi_R1\.fq\.gz$//')
[ -n "$libs" ] || { echo "no prep fastqs in $SRC -- run pipeline_fs.sh prep first" >&2; exit 1; }
echo "ribo-depleting $(echo "$libs" | wc -l) libraries, $NPAR at a time, stranded=$STRANDED"

one() {
    lib=$1
    out="$DST/${lib}_cbc_noumi"
    # EVERY intermediate ribo-bwamem.sh creates, not just the final ones. Three
    # separate tools in this stage refuse to overwrite and fail instead:
    #   gzip           riboread-selection.py ends with a bare `gzip`, which will
    #                  not replace an existing .gz and exits without the caller
    #                  noticing -- that cost a whole VASA run on 2026-07-27.
    #   samtools merge "File '...all-ribo.bam' exists. Please apply '-f'. Abort."
    #   samtools sort  "failed to create temporary file ...tmp.0000.bam: File
    #                  exists" -- the shards, which are named deterministically.
    # Missing the last two cost library ZHA8833A1 in job 51307376: leftovers from
    # an interrupted run made every later attempt fail on files it had itself
    # written. A stale intermediate must never be able to survive a retry.
    rm -f "${out}.nonRibo.fastq.gz" "${out}.nonRibo.fastq" "${out}.Ribo.bam" \
          "${out}.aln-ribo.bam" "${out}.mem-ribo.bam" "${out}.all-ribo.bam" \
          "${out}.nsorted.all-ribo.bam" "${out}.nsorted.all-ribo.bam".tmp.*.bam \
          "${out}_R1.fq.gz"
    "$SCRIPTS/ribo-bwamem.sh" "$RRNA" "$SRC/${lib}_cbc_noumi_R1.fq.gz" "$out" \
        "$P2BWA" "$P2SAMTOOLS" "$STRANDED" "$SCRIPTS" \
        > "$LOG/ribo_${lib}.log" 2>&1 || { echo "  FAILED ribo: $lib"; return 1; }
    ln -sfn "${out}.nonRibo.fastq.gz" "${out}_R1.fq.gz"
}
# `one` runs as a background job of THIS shell, not via `bash -c`, so it needs
# no `export -f` and the variables above are already in scope.

# verify_lib <lib> -- is this library's output complete, and newer than its
# input? Prints what is wrong, returns non-zero if anything is.
#
# ONE definition of "done", used twice: to skip already-finished libraries on a
# rerun, and as the final gate. Two separate implementations of the same
# predicate is precisely the thing that drifts silently.
verify_lib() {
    local lib=$1 bad=0 f
    local in="$SRC/${lib}_cbc_noumi_R1.fq.gz"
    local out="$DST/${lib}_cbc_noumi.nonRibo.fastq.gz"
    for f in "$out" "$DST/${lib}_cbc_noumi.Ribo.bam" "$DST/${lib}_cbc_noumi.ribo-map.log"; do
        [ -s "$f" ] || { echo "  BAD $lib: missing $(basename "$f")"; bad=$((bad+1)); }
    done
    if [ -s "$out" ] && [ "$out" -ot "$in" ]; then
        echo "  BAD $lib: nonRibo older than its input"; bad=$((bad+1))
    fi
    if [ -s "$DST/${lib}_cbc_noumi.nonRibo.fastq" ]; then
        echo "  BAD $lib: uncompressed .nonRibo.fastq left behind"; bad=$((bad+1))
    fi
    [ -e "$DST/${lib}_cbc_noumi_R1.fq.gz" ] || { echo "  BAD $lib: symlink does not resolve"; bad=$((bad+1)); }
    [ "$bad" -eq 0 ]
}

# --- what still needs doing --------------------------------------------------
# Resumable, so an interrupted run costs only the libraries it had not finished.
# FORCE=1 redoes everything. "Complete" here is verify_lib, i.e. evidence on
# disk -- not a marker file, which can outlive the output it claims.
todo=""; skipped=0
for lib in $libs; do
    if [ "${FORCE:-0}" = 0 ] && verify_lib "$lib" >/dev/null 2>&1; then
        skipped=$((skipped+1))
    else
        todo="$todo $lib"
    fi
done
[ "$skipped" -eq 0 ] || echo "  $skipped already complete, skipping (FORCE=1 to redo)"

# --- fan out -----------------------------------------------------------------
# NOT `xargs -P`. GNU xargs ABANDONS every remaining input the moment one child
# is killed by a signal. Reproduced 2026-08-06: 10 inputs at -P2, one child
# SIGTERMed -> 3 ran, 7 were silently never launched, xargs exits 125.
#
# That is what happened to job 51253372 on 2026-08-05. Something outside this
# script killed one of the four in-flight libraries; xargs dropped the other
# six; the `|| true` that used to sit on this line turned exit 125 into 0, and
# SLURM recorded the job COMPLETED having done 4 of 10 libraries. Had the
# verify gate below not existed, four-tenths of FLASH-seq would have gone into
# the four-way comparison looking like all of it.
#
# A plain job-control fan-out launches every library regardless of what happens
# to its siblings, and nothing here decides success from an exit code -- the
# gate is what is on disk.
for lib in $todo; do
    while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n 2>/dev/null || true; done
    one "$lib" &
done
wait 2>/dev/null || true

# --- verify, because counting files is not proof -----------------------------
bad=0; done_n=0
for lib in $libs; do
    if verify_lib "$lib"; then done_n=$((done_n+1)); else bad=$((bad+1)); fi
done
echo "  $done_n/$(echo "$libs" | wc -l) libraries complete"
[ "$bad" -eq 0 ] || { echo "ribo FAILED on $bad libraries -- do NOT run map" >&2; exit 1; }

# Ribosomal is `reads - unmapped`, NOT the sum of the log's aln/mem categories.
# The two differ by exactly 1 read per unit (a category the log does not print),
# and only `unmapped` is definitionally right: it is the read count of the
# nonRibo fastq, and STAR's "Number of input reads" reproduces it to the read on
# the VASA side (72,632,400).
printf '\n%-12s %14s %14s %8s\n' library reads_in ribosomal pct
tot_in=0; tot_r=0
for lib in $libs; do
    l="$DST/${lib}_cbc_noumi.ribo-map.log"
    read -r n u <<< "$(awk '/Number of reads:/{gsub(/[^0-9]/,"",$NF); n=$NF}
                            /Number of unmapped reads:/{gsub(/[^0-9]/,"",$NF); u=$NF}
                            END{print n, u}' "$l")"
    printf '%-12s %14d %14d %7.2f%%\n' "$lib" "$n" "$((n-u))" "$(echo "100*($n-$u)/$n" | bc -l)"
    tot_in=$((tot_in+n)); tot_r=$((tot_r+n-u))
done
printf '%-12s %14d %14d %7.2f%%\n' ALL "$tot_in" "$tot_r" "$(echo "100*$tot_r/$tot_in" | bc -l)"
echo
echo "ribo done. map/assign next, with FSV_OUTDIR=$RUN/nonribo"
