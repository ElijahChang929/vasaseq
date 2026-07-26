#!/bin/bash
###############################################################################
# bench_trim2.sh -- round 2. Round 1 showed the current XX{5} homopolymer pass
# LOSES uniquely mapped reads versus doing nothing at all (v4). So instead of
# guessing at a whole new recipe, build up from v4 one modifier at a time and
# see which ones pay for themselves.
#
#   v4  TrimGalore only, -m 20                       (round-1 control, reused)
#   v5  v4 + trim the real 3' construct
#           polyA + 12 nt (rc CBC+UMI) + rc(RA5 prefix)
#   v6  v5 + --poly-a           (clean the residual poly-A tail)
#   v7  v6 + --nextseq-trim=20  (2-colour "no signal" poly-G, NovaSeq X)
#   v8  v6 + -a polyG=G{20}     (poly-G as a plain adapter instead)
#
#   sbatch bench_trim2.sh
###############################################################################
#SBATCH --job-name=trimbench2
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim2.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 007 001 016"
EBROOT=/camp/apps/eb/software
STAR_INDEX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
NT=16
VARIANTS="v5 v6 v7 v8"

ML_INIT="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all"
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
CONDA="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa"

RA5RC=GATCGTCGGACTGTAGAACTCTGAAC
VASA3="A{8}NNNNNNNNNNNN${RA5RC:0:16}"     # polyA + 12 random + start of rc(RA5)

say(){ echo "[$(date +%H:%M:%S)] $*"; }

run_v5(){ cutadapt -j 4 -m 20 --trim-n -n 2 \
            -a "vasa3=${VASA3}" -a "RA5rc=${RA5RC}" -o "$2" "$1"; }
run_v6(){ cutadapt -j 4 -m 20 --trim-n -n 2 \
            -a "vasa3=${VASA3}" -a "RA5rc=${RA5RC}" --poly-a -o "$2" "$1"; }
run_v7(){ cutadapt -j 4 -m 20 --trim-n -n 2 --nextseq-trim=20 \
            -a "vasa3=${VASA3}" -a "RA5rc=${RA5RC}" --poly-a -o "$2" "$1"; }
run_v8(){ cutadapt -j 4 -m 20 --trim-n -n 2 \
            -a "vasa3=${VASA3}" -a "RA5rc=${RA5RC}" -a "polyG=G{20}" --poly-a -o "$2" "$1"; }

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

say "bench2 done"
