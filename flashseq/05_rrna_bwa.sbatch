#!/bin/bash
#SBATCH --job-name=fs_ribo
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%j_fs_ribo.log
#
# Ten libraries x two arms (raw, trimmed) x (bwa aln + bwa mem) over an
# every-64th-read sample of each library, ~456k reads.
#
# Alignment is NOT the cost -- the rRNA index is 71 kB, and on 456k 151 bp reads
# `bwa aln` takes 4.8 s and `bwa mem -t 8` 8.7 s (measured 2026-07-27). The cost
# is decompressing each 1.5-2.6 GB fastq once to reach the whole file, at ~57 s
# per file, x2 mates x10 libraries. 8 cores is what `bwa mem -t 8` and
# samtools --threads 8 inside ribo-bwamem.sh actually ask for.
#
# Measured, job 50855065: 46m11s elapsed, 1h15m CPU, MaxRSS 778 MB -- about
# 4.5 min per library, of which ~2 min is the two decompress passes. TotalCPU
# over Elapsed is 1.6, so most of the 8 cores sit idle waiting on gzip; they are
# requested because ribo-bwamem.sh asks for them, not because they are used.
# The 4 h ask is headroom, not a prediction.
#
#   sbatch code/flashseq/05_rrna_bwa.sbatch

set -euo pipefail

FS_CODE="${FS_CODE:-/nemo/lab/turnerj/working/guangxin/vasaseq/code/flashseq}"
export FS_NCORES="${SLURM_CPUS_PER_TASK:-8}"

bash "$FS_CODE/05_rrna_bwa.sh"
