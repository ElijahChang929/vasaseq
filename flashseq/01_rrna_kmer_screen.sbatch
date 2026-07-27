#!/bin/bash
#SBATCH --job-name=fs_kmer
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=%j_fs_kmer.log
#
# Was a login-node script when it screened the first 400k reads of each library.
# It now samples every FS_STRIDE-th read across the WHOLE file -- head sampling
# is measurably biased, see config.sh -- so each of the ten 1.5-2.6 GB fastqs is
# decompressed end to end. Python's gzip is several times slower than zcat, so
# budget minutes per library rather than seconds.
#
# Single-threaded and memory-flat: the k-mer sets are built once (a few hundred
# thousand 31-mers from a 71 kB reference) and the reads stream past one at a
# time.
#
# Measured, job 50855262: 10m09s elapsed, 9m48s CPU, MaxRSS 157 MB. The 4 h ask
# is headroom for a busy Lustre, not a prediction.
#
#   sbatch code/flashseq/01_rrna_kmer_screen.sbatch

set -euo pipefail

FS_CODE="${FS_CODE:-/nemo/lab/turnerj/working/guangxin/vasaseq/code/flashseq}"
source "$FS_CODE/config.sh"
fs_check_python

fs_python "$FS_CODE/01_rrna_kmer_screen.py"
