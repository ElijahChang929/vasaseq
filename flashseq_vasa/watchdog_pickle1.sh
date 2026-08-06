#!/bin/bash
###############################################################################
# watchdog_pickle1.sh -- resubmit FLASH-seq stage-6 jobs that run out of
# walltime, and rebuild the stage-6/7 job that depends on all ten.
#
# WHY THIS EXISTS AS A SLURM JOB AND NOT AS A PERSON WATCHING
# ------------------------------------------------------------
# The ten pickle1 jobs were submitted at -t 8:00:00 on evidence that turned out
# to be thin: the first two to finish took 1h50 and 2h18 on 0.44/0.45 GB of
# BED, so cost looked linear in BED size. It is not -- those two are the 1.5 ng
# and 3 ng libraries, and the later ones ran 6h02 and 6h23 on similar input.
# The remaining six sit at 6h32 with 1h27 of walltime left, which is inside the
# observed spread rather than safely past it.
#
# The failure mode is what makes this worth automating. `fsp67` waits on
# `afterok` for all ten. A TIMEOUT is not an ok exit, so fsp67 does not fail --
# it turns into DependencyNeverSatisfied and sits in the queue forever, quietly,
# with no error anywhere. Overnight that costs the whole night.
#
# SAFETY -- WHAT IT WILL NOT DO
# ------------------------------
# * It resubmits a given library AT MOST ONCE (a marker file per library), so a
#   job that fails for a real reason cannot loop.
# * It only ever acts on TIMEOUT. FAILED and OUT_OF_MEMORY are reported and
#   left alone, because those need a diagnosis, not a longer clock.
# * It decides "done" from the product on disk -- <lib>dict.pickle non-empty --
#   never from an exit code. This pipeline has already produced a COMPLETED job
#   that did 40% of its work (the xargs -P bug, see fourway_analyze/scripts/
#   00_fs_ribo.sh) and an OUT_OF_MEMORY job whose output was complete.
# * fsp67 is rebuilt only when every one of the ten has its pickle on disk.
#
# It exits as soon as there is nothing left to wait for, so it does not hold a
# node overnight for nothing.
#
# Usage:  sbatch -c 1 --mem=2G -t 24:00:00 -o <log> --wrap watchdog_pickle1.sh
###############################################################################
set -uo pipefail

FSV=/nemo/lab/turnerj/working/guangxin/vasaseq/code/flashseq_vasa
RUN=/nemo/lab/turnerj/working/guangxin/vasaseq/data/flashseq_vasa/run
OUTD=$RUN/nonribo
STATE=$RUN/watchdog; mkdir -p "$STATE"
LIBS="ZHA8833A1 ZHA8833A2 ZHA8833A3 ZHA8833A4 ZHA8833A5 ZHA8833A6 ZHA8833A7 ZHA8833A8 ZHA8833A9 ZHA8833A10"
E="FSV_OUTDIR=$OUTD FSV_ARM=native FSV_SAMPLE=FS_nonribo"
INTERVAL=${INTERVAL:-300}

say() { echo "[$(date '+%F %T')] $*"; }

# Product, not exit code.
done_lib() { [ -s "$OUTD/pickle/$1/${1}dict.pickle" ]; }

# Newest job for this library, whatever its state.
last_job() { sacct -n -X --name="p1_$1" --format=JobID,State -S 2026-08-06 2>/dev/null | tail -1; }

say "watching $(echo $LIBS | wc -w) libraries; poll ${INTERVAL}s"

while :; do
    pending=0 resubmitted=0
    for lib in $LIBS; do
        done_lib "$lib" && continue
        pending=$((pending+1))
        read -r jid state <<<"$(last_job "$lib")"
        [ -n "${state:-}" ] || { say "$lib: no job found yet"; continue; }
        case "$state" in
            RUNNING|PENDING) ;;
            TIMEOUT)
                if [ -e "$STATE/$lib.retried" ]; then
                    say "$lib: TIMEOUT again after a retry -- NOT resubmitting, needs a look"
                else
                    new=$(sbatch --parsable -J "p1_$lib" --chdir="$FSV" \
                        -o "$RUN/pickle1_${lib}.%j.out" -c 4 --mem=96G -t 20:00:00 \
                        --wrap "$E FSV_LIBS=$lib ./pipeline_fs.sh pickle1")
                    date > "$STATE/$lib.retried"
                    say "$lib: TIMEOUT on $jid -> resubmitted as $new (-t 20:00:00, --mem=96G)"
                    resubmitted=$((resubmitted+1))
                fi ;;
            *)  say "$lib: $state on $jid -- left alone (only TIMEOUT is auto-retried)" ;;
        esac
    done

    if [ "$pending" -eq 0 ]; then
        say "all 10 pickles present on disk"
        # fsp67's afterok cannot be satisfied once any of the ten has been
        # replaced, so it is rebuilt rather than waited on.
        if squeue -h -n fsp67 -u "$USER" -o %T 2>/dev/null | grep -q .; then
            old=$(squeue -h -n fsp67 -u "$USER" -o %i); scancel $old 2>/dev/null
            say "cancelled stale fsp67 ($old)"
        fi
        if [ ! -e "$STATE/fsp67.submitted" ]; then
            new=$(sbatch --parsable -J fsp67 --chdir="$FSV" \
                -o "$RUN/tables.%j.out" -c 8 --mem=160G -t 12:00:00 \
                --wrap "$E ./pipeline_fs.sh pickle_merge tables recon")
            date > "$STATE/fsp67.submitted"
            say "submitted fsp67 as $new (no dependency -- all inputs verified on disk)"
        fi
        say "nothing left to watch; exiting"
        exit 0
    fi

    say "$pending libraries outstanding, $resubmitted resubmitted this pass"
    sleep "$INTERVAL"
done
