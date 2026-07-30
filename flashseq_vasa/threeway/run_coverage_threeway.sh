#!/bin/bash
# run_coverage_threeway.sh -- driver for the three-way transcript-body coverage
# comparison. One job, sequential stages, because the stages depend on each
# other and the whole thing is CPU-bound single-threaded streaming.
#
# WHY ONE JOB AND NOT AN ARRAY: the gene set has to exist before any BAM is
# profiled, and every profile must use the SAME gene set in the SAME order or
# the bins do not mean the same thing. Sequencing them inside one job makes that
# structural rather than a thing to remember.
#
# Trap 1: `ls <no match> | wc -l` exits 2 and kills the job under pipefail.
# Everything here counts with `find`.
# Trap SIGPIPE: no `... | head -N` on a large stream anywhere in this script.
set -euo pipefail

W=/nemo/lab/turnerj/working/guangxin/vasaseq
REF=/nemo/lab/turnerj/working/guangxin/reference
SCR=/nemo/lab/turnerj/scratch/zhangg/vasaseq
PY=/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/python
SRC=$W/code/flashseq_vasa/threeway/mk_coverage_threeway.py

RES=$W/res/threeway
COV=$SCR/threeway/cov
mkdir -p "$RES" "$COV"

GTF99=$REF/vasaseq/mixed/build/combined.gtf
GTF116=$REF/vasaseq/mouse_GRCm39_E116/build/mouse.gtf

PUB=$W/data/ref/fastq_vasaplate
OWN=$W/data/PM26037/out
FSN=$SCR/flashseq_vasa/native
FSV=$SCR/flashseq_vasa/vasalen

echo "=================================================================="
echo "host      : $(hostname)   job ${SLURM_JOB_ID:-none}"
echo "started   : $(date -Is)"
$PY -c "import sys,numpy,pandas,pysam;print('python',sys.version.split()[0],
'numpy',numpy.__version__,'pandas',pandas.__version__,'pysam',pysam.__version__)"
echo "script md5: $(md5sum "$SRC" | cut -d' ' -f1)"
echo "=================================================================="

# ---------------------------------------------------------------- stage 1
echo; echo "### STAGE 1  published-plate cell selection (species purity)"
# The UFI table family matters and is not interchangeable. The established
# call table (res/threeway/threeway_published_cellcalls.tsv, written by
# threeway_paperform.call_published_cells) was computed from the uniaggGenes
# family, where multi-gene combination entries are collapsed. Reading the raw
# `total` family instead shifts the Methods rule's gene fractions and moves
# barcodes 090/208/296 mixed->mouse and 245 discarded->human, giving 144 mouse
# instead of 141. select() asserts agreement with the established table, so a
# wrong path here fails the job in stage 1 rather than silently profiling a
# different cell set.
$PY "$SRC" select \
    "$PUB/vasaplate_out_v3_uniaggGenes_total.UFICounts.tsv" \
    "$RES/coverage_threeway_pubcells.tsv"

# ---------------------------------------------------------------- stage 2
# ReadCounts (Rule 4: reads are the only unit all three protocols share), and
# the uniaggGenes family on every side -- the same family stage 1 called cells
# from. See mk_coverage_threeway.geneset()'s note.
echo; echo "### STAGE 2  gene set + transcript models for BOTH releases"
$PY "$SRC" geneset \
    "$OWN/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv" \
    "$FSN/FSall10_native_uniaggGenes_total.ReadCounts.tsv" \
    "$FSV/FSall10_vasalen_uniaggGenes_total.ReadCounts.tsv" \
    "$PUB/vasaplate_out_v3_uniaggGenes_total.ReadCounts.tsv" \
    "$RES/coverage_threeway_pubcells.tsv" \
    "$GTF99" "$GTF116" \
    "$RES/coverage_threeway_geneset"

M99=$RES/coverage_threeway_geneset.E99.npz
M116=$RES/coverage_threeway_geneset.E116.npz
test -s "$M99" && test -s "$M116"

# ---------------------------------------------------------------- stage 3
# One process per BAM. The units share nothing but the read-only .npz model
# file, so this is embarrassingly parallel; the precheck measured the serial
# cost at ~100 core-hours, which is not acceptable as wall time.
#
# `wait` alone returns 0 whatever the children did, so every pid is captured and
# waited on individually -- otherwise a failed unit would leave a silently
# incomplete merge.
NPROC=${SLURM_CPUS_PER_TASK:-8}
echo; echo "### STAGE 3  per-unit profiles, $NPROC at a time"

WORK=$COV/logs
mkdir -p "$WORK"
declare -a PIDS=() NAMES=()

launch () {   # launch <models> <bam> <label> <group>
  local M=$1 B=$2 LAB=$3 GRP=$4
  test -s "$B" || { echo "MISSING $B"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 5; done
  $PY "$SRC" profile "$M" "$B" "$LAB" "$GRP" "$COV/$LAB" \
      > "$WORK/$LAB.log" 2>&1 &
  PIDS+=($!); NAMES+=("$LAB")
  echo "  launched $LAB (pid ${PIDS[-1]})"
}

# published plate: the mouse-pure cells stage 1 selected, against E99 models
PUBCELLS=$($PY - "$RES/coverage_threeway_pubcells.tsv" <<'PY'
import sys, csv
with open(sys.argv[1]) as fh:
    print(' '.join(r['cell'] for r in csv.DictReader(fh, delimiter='\t')
                   if r['selected'] == 'True'))
PY
)
echo "published mouse-pure cells: $PUBCELLS"
for CELL in $PUBCELLS; do
  launch "$M99" \
    "$PUB/vasaplate_out_v3/SRR14783059_${CELL}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam" \
    "VASApub_$CELL" VASA_published
done

# own plate: the 6 deepest real cells -- the same six the upstream script used,
# so the own-plate number here is comparable to the already-reported one
for CELL in 007 009 010 011 012 013; do
  launch "$M116" \
    "$OWN/cells/ZHA9292A1_${CELL}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam" \
    "VASAown_$CELL" VASA_own
done

# FLASH-seq: the same four libraries upstream used, both arms.
# A8 is excluded by the user's qc_verdict (18.3% human CALB1) and is not among
# them; A7 (caveat, 3.6% CALB1) is likewise absent. Neither is silently dropped
# -- the selection matches upstream's so the regression check is meaningful.
for LIB in A1 A5 A9 A10; do
  for ARM in native vasalen; do
    launch "$M116" \
      "$SCR/flashseq_vasa/$ARM/cells/ZHA8833${LIB}_cbc_noumi_E99_Aligned.out.bam" \
      "FS${ARM}_ZHA8833$LIB" "FLASHseq_${ARM}"
  done
done

echo "waiting on ${#PIDS[@]} units"
FAILED=0
for i in "${!PIDS[@]}"; do
  if wait "${PIDS[$i]}"; then
    echo "  OK   ${NAMES[$i]}"
  else
    echo "  FAIL ${NAMES[$i]} -- last 20 lines:"
    tail -20 "$WORK/${NAMES[$i]}.log" | sed 's/^/       /'
    FAILED=$((FAILED+1))
  fi
done
test "$FAILED" -eq 0 || { echo "$FAILED unit(s) failed; refusing to merge a partial set"; exit 1; }

NCOV=$(find "$COV" -name '*.cov.tsv' | wc -l)
echo "per-unit profiles written: $NCOV of ${#PIDS[@]}"
test "$NCOV" -eq "${#PIDS[@]}" || { echo "count mismatch"; exit 1; }

echo "--- per-unit summary lines"
for N in "${NAMES[@]}"; do grep -h 'placed of' "$WORK/$N.log" || true; done

# ---------------------------------------------------------------- stage 4
echo; echo "### STAGE 4  merge, statistics, regression check"
$PY "$SRC" merge "$COV" "$RES" "$RES/coverage_threeway_geneset"

echo; echo "### STAGE 5  annotation-release bound"
$PY "$SRC" crossrel "$COV" "$RES"

echo; echo "finished : $(date -Is)"
find "$RES" -name 'coverage_threeway*' -printf '%10s  %f\n' | sort -k2
