#!/bin/bash
###############################################################################
# bench_trim3.sh -- round 3.
#
# Round 2's v5-v8 all lost ~8k uniquely mapped reads per 300k versus v4 (no
# second pass at all). The culprit is the wildcard block in the "vasa3" adapter:
# A{8}NNNNNNNNNNNN... . cutadapt allows a 3' adapter to match partially at the
# read end, the 12 N's match ANY base, so the last ~20 bases of essentially any
# read that happens to end in a couple of A's are eaten for free. Wildcards and
# 3'-partial matching do not mix.
#
# Round 3 drops the wildcards and tests only literal, high-min-overlap adapters:
#   v9   --poly-a only
#   v10  the real read-through adapter rc(RA5), min_overlap 8
#   v11  v10 + --poly-a
#   v12  v11 + poly-G (2-colour artefact) at min_overlap 10, -m 25
# v4 (TrimGalore only) stays the control to beat.
#
#   sbatch bench_trim3.sh
###############################################################################
#SBATCH --job-name=trimbench3
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim3.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 007 001 016"
EBROOT=/camp/apps/eb/software
STAR_INDEX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
NT=16
VARIANTS="v9 v10 v11 v12"

ML_INIT="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all"
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
CONDA="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa"

RA5RC=GATCGTCGGACTGTAGAACTCTGAAC

say(){ echo "[$(date +%H:%M:%S)] $*"; }

run_v9() { cutadapt -j 4 -m 20 --trim-n --poly-a -o "$2" "$1"; }
run_v10(){ cutadapt -j 4 -m 20 --trim-n -a "RA5rc=${RA5RC};min_overlap=8" -o "$2" "$1"; }
run_v11(){ cutadapt -j 4 -m 20 --trim-n -a "RA5rc=${RA5RC};min_overlap=8" --poly-a -o "$2" "$1"; }
run_v12(){ cutadapt -j 4 -m 25 --trim-n -a "RA5rc=${RA5RC};min_overlap=8" \
             -a "polyG=G{20};min_overlap=10" --poly-a -o "$2" "$1"; }

( eval "$CONDA"
  for v in $VARIANTS; do
    mkdir -p "$v"
    for c in $CELLS; do
      [ -s "${v}/${c}.fq.gz" ] && continue
      say "$v $c"
      run_$v "tg/${c}_trimmed.fq.gz" "${v}/${c}.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
    done
  done )

( eval "$ML_INIT"; module load STAR/2.7.7a-GCC-10.2.0
  tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
  release(){ $STAR --genomeDir "$STAR_INDEX" --genomeLoad Remove \
               --outFileNamePrefix "${tmp}/rm_" >/dev/null 2>&1; rm -rf "$tmp"; }
  trap release EXIT
  say "loading STAR genome"
  $STAR --genomeDir "$STAR_INDEX" --genomeLoad LoadAndExit \
        --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1
  for v in $VARIANTS; do
    for c in $CELLS; do
      pfx="${v}/${c}_"
      [ -s "${pfx}Log.final.out" ] && continue
      say "STAR $v $c"
      $STAR --runThreadN $NT --genomeDir "$STAR_INDEX" --genomeLoad LoadAndKeep \
            --readFilesIn "${v}/${c}.fq.gz" --readFilesCommand zcat \
            --outFilterMultimapNmax 20 --outSAMtype None \
            --outFileNamePrefix "$pfx" >/dev/null 2>&1
      rm -rf "${pfx}_STARtmp" "${pfx}Log.progress.out" "${pfx}Log.out"
    done
  done )

say "bench3 done"
