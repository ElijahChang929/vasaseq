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
$PY "$SRC" select \
    "$PUB/vasaplate_out_v3_total.UFICounts.tsv" \
    "$RES/coverage_threeway_pubcells.tsv"

# ---------------------------------------------------------------- stage 2
echo; echo "### STAGE 2  gene set + transcript models for BOTH releases"
$PY "$SRC" geneset \
    "$OWN/ZHA9292A1_total.ReadCounts.tsv" \
    "$FSN/FSall10_native_total.ReadCounts.tsv" \
    "$FSV/FSall10_vasalen_total.ReadCounts.tsv" \
    "$PUB/vasaplate_out_v3_total.ReadCounts.tsv" \
    "$RES/coverage_threeway_pubcells.tsv" \
    "$GTF99" "$GTF116" \
    "$RES/coverage_threeway_geneset"

M99=$RES/coverage_threeway_geneset.E99.npz
M116=$RES/coverage_threeway_geneset.E116.npz
test -s "$M99" && test -s "$M116"

# ---------------------------------------------------------------- stage 3
# published plate: the cells stage 1 selected, against the E99 models
echo; echo "### STAGE 3a  published VASA plate (E99/GRCm38 models)"
PUBCELLS=$($PY - "$RES/coverage_threeway_pubcells.tsv" <<'PY'
import sys, csv
with open(sys.argv[1]) as fh:
    print(' '.join(r['cell'] for r in csv.DictReader(fh, delimiter='\t')
                   if r['selected'] == 'True'))
PY
)
echo "cells: $PUBCELLS"
for CELL in $PUBCELLS; do
  B=$PUB/vasaplate_out_v3/SRR14783059_${CELL}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam
  test -s "$B" || { echo "MISSING $B"; exit 1; }
  $PY "$SRC" profile "$M99" "$B" "VASApub_$CELL" VASA_published \
      "$COV/VASApub_$CELL"
done

# own plate: the 6 deepest real cells, matching the upstream script's selection
echo; echo "### STAGE 3b  own VASA plate (E116/GRCm39 models)"
for CELL in 007 009 010 011 012 013; do
  B=$OWN/cells/ZHA9292A1_${CELL}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam
  test -s "$B" || { echo "MISSING $B"; exit 1; }
  $PY "$SRC" profile "$M116" "$B" "VASAown_$CELL" VASA_own \
      "$COV/VASAown_$CELL"
done

# FLASH-seq: same four libraries the upstream script used, both arms
echo; echo "### STAGE 3c  FLASH-seq, both arms (E116/GRCm39 models)"
for LIB in A1 A5 A9 A10; do
  for ARM in native vasalen; do
    D=$SCR/flashseq_vasa/$ARM
    B=$D/cells/ZHA8833${LIB}_cbc_noumi_E99_Aligned.out.bam
    test -s "$B" || { echo "MISSING $B"; exit 1; }
    $PY "$SRC" profile "$M116" "$B" "FS${ARM}_ZHA8833$LIB" \
        "FLASHseq_${ARM}" "$COV/FS${ARM}_ZHA8833$LIB"
  done
done

NCOV=$(find "$COV" -name '*.cov.tsv' | wc -l)
echo; echo "per-unit profiles written: $NCOV"

# ---------------------------------------------------------------- stage 4
echo; echo "### STAGE 4  merge, statistics, regression check"
$PY "$SRC" merge "$COV" "$RES" "$RES/coverage_threeway_geneset"

echo; echo "### STAGE 5  annotation-release bound"
$PY "$SRC" crossrel "$COV" "$RES"

echo; echo "finished : $(date -Is)"
find "$RES" -name 'coverage_threeway*' -printf '%10s  %f\n' | sort -k2
