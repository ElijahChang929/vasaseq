#!/bin/bash
###############################################################################
# bench_bam.sh -- is v13's extra yield real, or poly-A sticking to genomic
# A-tracts?
#
# v13 leaves the poly-A tail on the read (it only removes the read-through
# adapter), so a read that is 8 nt of insert followed by 100 nt of A could in
# principle align to a genomic poly-A stretch and be counted. This re-maps the
# three candidates keeping the BAM, so aligned_composition.py can look at the
# ALIGNED portion of every uniquely mapped read and say what fraction is A-rich.
#
#   sbatch bench_bam.sh
###############################################################################
#SBATCH --job-name=trimbam
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_bam.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

EBROOT=/camp/apps/eb/software
STAR_INDEX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
module load STAR/2.7.7a-GCC-10.2.0

mkdir -p bam
tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
release(){ $STAR --genomeDir "$STAR_INDEX" --genomeLoad Remove \
             --outFileNamePrefix "${tmp}/rm_" >/dev/null 2>&1; rm -rf "$tmp"; }
trap release EXIT
$STAR --genomeDir "$STAR_INDEX" --genomeLoad LoadAndExit \
      --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1

for v in v1 v13 v14 v15 v16; do
  for c in 011 016; do
    pfx="bam/${v}_${c}_"
    [ -s "${pfx}Aligned.out.bam" ] && continue
    echo "[$(date +%H:%M:%S)] BAM $v $c"
    $STAR --runThreadN 16 --genomeDir "$STAR_INDEX" --genomeLoad LoadAndKeep \
          --readFilesIn "${v}/${c}.fq.gz" --readFilesCommand zcat \
          --outFilterMultimapNmax 20 --outSAMtype BAM Unsorted \
          --outSAMattributes All --outFileNamePrefix "$pfx" >/dev/null 2>&1
    rm -rf "${pfx}_STARtmp" "${pfx}Log.progress.out" "${pfx}Log.out"
  done
done
echo "bench_bam done"
