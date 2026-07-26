#!/bin/bash
###############################################################################
# bench_trim6.sh -- round 6. What is still left in the reads after v17?
#
# Measured on v17/011.fq.gz (164,798 reads):
#   adapter still present            0.26%
#   contains a poly-A run >= 20 nt   0.00%   <- the A{20} adapter clears these
#   ends in A>=6                     1.88%
#   revcomp cell barcode in last 25 nt  13.58%
#
# So the adapter and every long poly-A are gone, but the 12 nt
# revcomp(CBC+UMI) remnant survives in 13.6% of reads. It only disappears today
# as a side effect: A{20} is a 3' adapter, so removing the poly-A takes
# everything after it with it. When the poly-A run is SHORTER than 20 nt there
# is no match, and both the short run and the 12 nt stay.
#
# The barcode and UMI are already on the read name from step 1, so this is pure
# junk -- STAR soft-clips it, but it still counts against
# outFilterMatchNminOverLread.
#
#   v18  v17 + A{10} followed by 12 wildcards, min_overlap = FULL LENGTH.
#        Round 2 showed wildcards are dangerous, but only because cutadapt
#        allows a 3' adapter to match PARTIALLY at the read end, where 12 N's
#        match anything. Requiring the full 22 nt forbids partial matches and
#        removes the hazard.
#   v19  v17 + --poly-a, to mop up short trailing runs
#   v20  v18 + --poly-a
#
#   sbatch bench_trim6.sh
###############################################################################
#SBATCH --job-name=trimbench6
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim6.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 007 001 016"
EBROOT=/camp/apps/eb/software
IDX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
VARIANTS="v18 v19 v20"

AD=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT
RT="rt=${AD};min_overlap=8"
PA="polyA=A{20};min_overlap=10"
PG="polyG=G{20};min_overlap=10"
# 10 A + 12 N, full-length match only (22) -- no partial 3'-end matching
PABC="polyAbc=A{10}NNNNNNNNNNNN;min_overlap=22"

say(){ echo "[$(date +%H:%M:%S)] $*"; }

source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
mkdir -p v18 v19 v20 bam

run_v18(){ cutadapt -j 4 -m 20 --trim-n -n 3 -a "$RT" -a "$PABC" -a "$PA" -a "$PG" -o "$2" "$1"; }
run_v19(){ cutadapt -j 4 -m 20 --trim-n -n 2 -a "$RT" -a "$PA" -a "$PG" --poly-a -o "$2" "$1"; }
run_v20(){ cutadapt -j 4 -m 20 --trim-n -n 3 -a "$RT" -a "$PABC" -a "$PA" -a "$PG" --poly-a -o "$2" "$1"; }

for v in $VARIANTS; do
  for c in $CELLS; do
    [ -s "${v}/${c}.fq.gz" ] && continue
    say "$v $c"
    run_$v "tg/${c}_trimmed.fq.gz" "${v}/${c}.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
  done
done

source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
module load STAR/2.7.7a-GCC-10.2.0
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
trap '$STAR --genomeDir $IDX --genomeLoad Remove --outFileNamePrefix ${tmp}/rm_ >/dev/null 2>&1; rm -rf $tmp' EXIT
$STAR --genomeDir "$IDX" --genomeLoad LoadAndExit --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1

for v in $VARIANTS; do
  for c in $CELLS; do
    [ -s "${v}/${c}_Log.final.out" ] && continue
    say "STAR $v $c"
    $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
      --readFilesIn "${v}/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
      --outSAMtype None --outFileNamePrefix "${v}/${c}_" >/dev/null 2>&1
    rm -rf "${v}/${c}__STARtmp" "${v}/${c}_Log.progress.out" "${v}/${c}_Log.out"
  done
  for c in 011 016; do
    [ -s "bam/${v}_${c}_Aligned.out.bam" ] && continue
    $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
      --readFilesIn "${v}/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
      --outSAMtype BAM Unsorted --outSAMattributes All \
      --outFileNamePrefix "bam/${v}_${c}_" >/dev/null 2>&1
    rm -rf "bam/${v}_${c}__STARtmp" "bam/${v}_${c}_Log.progress.out" "bam/${v}_${c}_Log.out"
  done
done
echo "bench6 done"
