#!/bin/bash
#SBATCH --job-name=fs_contam
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=%j_fs_contam.log
#
# Streams the mouse genome (2.7 GB plain) and the human genome (3.1 GB gzipped)
# once each, then re-counts every classified sequence across all 10 libraries.
# Single-threaded and memory-flat by construction: the search holds one 4 MB
# block at a time, so 8 GB is generous. 4 h covers gzip-decompression variance
# on a busy Lustre.
#
#   sbatch code/flashseq/02_contaminant_check.sbatch

set -euo pipefail

FS_CODE="${FS_CODE:-/nemo/lab/turnerj/working/guangxin/vasaseq/code/flashseq}"
source "$FS_CODE/config.sh"
fs_check_python

fs_python "$FS_CODE/02_contaminant_check.py"
