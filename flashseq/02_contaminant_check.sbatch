#!/bin/bash
#SBATCH --job-name=fs_contam
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=%j_fs_contam.log
#
# Streams the mouse genome (2.7 GB plain) and the human genome (3.1 GB gzipped)
# once each, then re-counts every classified sequence across all 10 libraries.
# Single-threaded and memory-flat by construction: the search holds one 4 MB
# block at a time, so 8 GB is generous.
#
# The 18 min this took at first was with stage 2 reading only the first 400k
# reads of each file. It now samples every FS_STRIDE-th read across the whole
# file, so all 20 fastqs (~36 GB compressed) are decompressed end to end.
# Measured after that change, job 50855331: 22m26s elapsed, MaxRSS 155 MB --
# less of a penalty than expected, because stage 1's two genome passes were
# always the larger half. 6 h is headroom, not a prediction.
#
#   sbatch code/flashseq/02_contaminant_check.sbatch

set -euo pipefail

FS_CODE="${FS_CODE:-/nemo/lab/turnerj/working/guangxin/vasaseq/code/flashseq}"
source "$FS_CODE/config.sh"
fs_check_python

fs_python "$FS_CODE/02_contaminant_check.py"
