#!/usr/bin/env bash
#SBATCH --partition=ncpu
#SBATCH --account=u_turnerj
#SBATCH --job-name=vp_after7
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=04:00:00
set -uo pipefail
W=/nemo/lab/turnerj/working/guangxin/vasaseq
VC=$W/code/I_Gene_expression/vasaplate_check
V=$W/data/ref/fastq_vasaplate/vasaplate_out_v3
R=$W/res/vasaplate

echo "=== stage 6/7 outcome, from sacct ==="
sacct -j 51029417,51029418 --format=JobID%14,JobName%18,State%14,Elapsed%10,MaxRSS%10,ExitCode%8 \
  | grep -vE "\.batch|\.extern|^-"

echo
echo "=== completeness: mapStats must be 21 lines (trap 2) ==="
for f in $V/*mapStats.log; do [ -f "$f" ] && echo "  $(basename $f): $(wc -l < $f) lines"; done
echo "  tsv tables: $(find $V -maxdepth 1 -name '*.tsv' | wc -l)"

echo
echo "=== PREDICTION 1: step-3 depletion should RISE 5.06% -> ~5.34% ==="
# sum the ribo logs across all 384 cells, same accounting as run 4
tot=0; rib=0
for f in $V/*.ribo.log $V/ribo-*.out; do :; done
python3 - "$V" <<'PYEOF'
import glob, os, re, sys
V = sys.argv[1]
# ribo-bwamem.sh writes "<kept> <ribo> <total>"-style lines into the per-cell log;
# find whichever log carries the counts rather than assuming the name
cands = {}
for pat in ('*ribo*.log', 'ribo-*.out', '*.nonRibo.log'):
    fs = glob.glob(os.path.join(V, pat))
    if fs:
        cands[pat] = len(fs)
print("candidate log globs:", cands)
PYEOF

echo
echo "=== run the published-table comparison ==="
PY=/camp/apps/eb/software/Anaconda3/2024.10-1/bin/python
env -u PYTHONPATH $PY "$VC/01_compare.py" v3 2>&1 | tail -40 ||   echo "01_compare.py failed -- vp_common.RUNS probably needs a 'v3' entry"

echo
echo "=== comparison outputs ==="
ls -la $R/*.tsv 2>/dev/null | awk '{print "  ", $5, $9}'
