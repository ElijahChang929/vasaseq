#!/bin/bash
###############################################################################
# bench_trim.sh -- compare homopolymer-trimming settings for step2 on ZHA9292A1
#
# Runs the SAME TrimGalore pass for every variant, then only varies the second
# (cutadapt homopolymer) pass, then maps each result with STAR. The score that
# decides the argument is "uniquely mapped reads per 300k input reads" -- not
# the mapping *rate*, which you can trivially inflate by throwing reads away.
#
#   sbatch bench_trim.sh
###############################################################################
#SBATCH --job-name=trimbench
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 007 001 016"
EBROOT=/camp/apps/eb/software
STAR_INDEX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
NT=16

ML_INIT="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all"
CA118=${EBROOT}/cutadapt/1.18-foss-2018b-Python-3.6.6/bin/cutadapt
TG=${EBROOT}/Trim_Galore/0.6.2-foss-2018b-Python-3.6.6/trim_galore
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
CONDA="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa"

# revcomp of the 21 nt RA5 prefix that sits at the 5' end of R1; on R2 it is
# what you read through into, 12 nt (revcomp CBC+UMI) after the poly-A ends.
RA5RC=GATCGTCGGACTGTAGAACTCTGAAC

say(){ echo "[$(date +%H:%M:%S)] $*"; }

###############################################################################
# 0. shared TrimGalore pass (adapters + Q<20), identical for every variant
###############################################################################
mkdir -p tg
if [ ! -s "tg/016_trimmed.fq.gz" ]; then
  ( eval "$ML_INIT"; module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
    for c in $CELLS; do
      say "trim_galore $c"
      $TG --path_to_cutadapt "$CA118" --cores 4 -o tg "in/${c}.fastq.gz" >/dev/null 2>&1
    done )
fi
ls -l tg/

###############################################################################
# 1. homopolymer variants
###############################################################################
run_v0() {   # what step2 does today: 6-mer homopolymer, cutadapt 1.18, -m 15
  $CA118 -m 15 --trim-n \
    -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" \
    -o "$2" "$1"
}
run_v1() {   # identical settings, cutadapt 5.1 -- isolates the version effect
  cutadapt -j 4 -m 15 --trim-n \
    -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" \
    -o "$2" "$1"
}
run_v2() {   # same idea, just a longer run required before cutting
  cutadapt -j 4 -m 20 --trim-n \
    -a "polyG=G{20}" -a "polyC=C{20}" -a "polyT=T{20}" -a "polyA=A{20}" \
    -o "$2" "$1"
}
run_v3() {   # model the actual 3' construct: polyA + 12 nt (rc CBC+UMI) + rc(RA5)
  cutadapt -j 4 -m 20 --trim-n --nextseq-trim=20 -n 3 \
    -a "vasa3=A{8}NNNNNNNNNNNN${RA5RC:0:16}" \
    -a "RA5rc=${RA5RC}" \
    -a "polyA=A{20}" -a "polyG=G{20}" \
    --poly-a \
    -o "$2" "$1"
}
run_v4() {   # control: no homopolymer pass at all
  cutadapt -j 4 -m 20 --trim-n -o "$2" "$1"
}

( eval "$CONDA"
  for v in v0 v1 v2 v3 v4; do
    mkdir -p "$v"
    for c in $CELLS; do
      [ -s "${v}/${c}.fq.gz" ] && continue
      say "$v $c"
      run_$v "tg/${c}_trimmed.fq.gz" "${v}/${c}.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
    done
  done )

###############################################################################
# 2. map every variant with STAR (genome held in shared memory once)
###############################################################################
( eval "$ML_INIT"; module load STAR/2.7.7a-GCC-10.2.0
  tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
  release(){ $STAR --genomeDir "$STAR_INDEX" --genomeLoad Remove \
               --outFileNamePrefix "${tmp}/rm_" >/dev/null 2>&1; rm -rf "$tmp"; }
  trap release EXIT
  say "loading STAR genome"
  $STAR --genomeDir "$STAR_INDEX" --genomeLoad LoadAndExit \
        --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1

  for v in v0 v1 v2 v3 v4; do
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

say "bench done -- summarise with ./summarise_trim.sh"
