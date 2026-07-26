#!/bin/bash
###############################################################################
# bench_trim4.sh -- round 4, with the read-through adapter MEASURED from the
# reads instead of inferred.
#
# Round 3's rc(RA5) guess ended ...GAACTCTGAAC. Anchoring 68,250 reads on the
# 16 nt that is unambiguously revcomp(R1's 5' prefix) and taking a per-position
# consensus gives ...GAACTCCTGTCTCTTATACACATCT... -- the last 5 bases of the
# guess were wrong, and what actually follows is the Nextera mosaic end. That
# still scored best in round 3 because cutadapt could match the correct 21 nt
# prefix partially; with the right sequence it should do better.
#
#   ADAPTER = revcomp(R1 5' prefix) + Nextera mosaic end
#
#   v13  ADAPTER, min_overlap 8                     (round-3 v10, fixed)
#   v14  v13 + poly-A as a strict A{20} adapter      (needs a real 20 nt run)
#   v15  v13 + --poly-a                             (cutadapt's greedy poly-A)
#   v16  v14 + poly-G, -m 25                        (candidate for production)
#
#   sbatch bench_trim4.sh
###############################################################################
#SBATCH --job-name=trimbench4
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim4.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 007 001 016"
EBROOT=/camp/apps/eb/software
STAR_INDEX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
NT=16
VARIANTS="v13 v14 v15 v16"

ML_INIT="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all"
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
CONDA="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa"

# measured consensus, see bench_trim4 header
ADAPTER=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT

say(){ echo "[$(date +%H:%M:%S)] $*"; }

run_v13(){ cutadapt -j 4 -m 20 --trim-n -a "rt=${ADAPTER};min_overlap=8" -o "$2" "$1"; }
run_v14(){ cutadapt -j 4 -m 20 --trim-n -n 2 -a "rt=${ADAPTER};min_overlap=8" \
             -a "polyA=A{20};min_overlap=10" -o "$2" "$1"; }
run_v15(){ cutadapt -j 4 -m 20 --trim-n -a "rt=${ADAPTER};min_overlap=8" --poly-a -o "$2" "$1"; }
run_v16(){ cutadapt -j 4 -m 25 --trim-n -n 2 -a "rt=${ADAPTER};min_overlap=8" \
             -a "polyA=A{20};min_overlap=10" -a "polyG=G{20};min_overlap=10" -o "$2" "$1"; }

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

say "bench4 done"
