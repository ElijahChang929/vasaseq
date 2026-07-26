#!/bin/bash
###############################################################################
# bench_trim5.sh -- round 5. Is poly-G trimming free?
#
# Round 4 picked v14 (read-through adapter + strict A{20} poly-A). v16 added
# poly-G but also raised -m 20 -> 25, so the two were confounded. v17 is v14
# plus poly-G and nothing else.
#
# Answer: free, and pointless on this run -- 84,878 -> 84,860 unique on cell
# 011, a difference of 18 reads in 300,000. Kept in the production settings
# anyway, because poly-G is the two-colour "no signal" artefact and this is a
# NovaSeq X run; the next flow cell may actually need it.
#
# v17 IS the adopted production setting. own_version/trim.sh with
# TRIM_MODE=vasa reproduces v17/{011,016}.fq.gz byte for byte (md5-checked).
#
#   sbatch bench_trim5.sh
###############################################################################
#SBATCH --job-name=trimbench5
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim5.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 007 001 016"
EBROOT=/camp/apps/eb/software
IDX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
AD=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT

source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
mkdir -p v17 bam
for c in $CELLS; do
  [ -s "v17/${c}.fq.gz" ] && continue
  cutadapt -j 4 -m 20 --trim-n -n 2 -a "rt=${AD};min_overlap=8" \
    -a "polyA=A{20};min_overlap=10" -a "polyG=G{20};min_overlap=10" \
    -o "v17/${c}.fq.gz" "tg/${c}_trimmed.fq.gz" > "v17/${c}.cutadapt.log" 2>&1
done

source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
module load STAR/2.7.7a-GCC-10.2.0
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
trap '$STAR --genomeDir $IDX --genomeLoad Remove --outFileNamePrefix ${tmp}/rm_ >/dev/null 2>&1; rm -rf $tmp' EXIT
$STAR --genomeDir "$IDX" --genomeLoad LoadAndExit --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1

for c in $CELLS; do
  [ -s "v17/${c}_Log.final.out" ] && continue
  $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
    --readFilesIn "v17/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
    --outSAMtype None --outFileNamePrefix "v17/${c}_" >/dev/null 2>&1
  rm -rf "v17/${c}__STARtmp" "v17/${c}_Log.progress.out" "v17/${c}_Log.out"
done

# BAMs, so aligned_composition.py can look at what actually aligned
for c in 011 016; do
  [ -s "bam/v17_${c}_Aligned.out.bam" ] && continue
  $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
    --readFilesIn "v17/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
    --outSAMtype BAM Unsorted --outSAMattributes All \
    --outFileNamePrefix "bam/v17_${c}_" >/dev/null 2>&1
  rm -rf "bam/v17_${c}__STARtmp" "bam/v17_${c}_Log.progress.out" "bam/v17_${c}_Log.out"
done
echo "bench5 done"
