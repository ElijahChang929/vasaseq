#!/bin/bash
###############################################################################
# bench_trim8.sh -- round 8. Does the reverse orientation exist, and is it
# worth trimming?
#
# The forward 3' construct is
#     [insert][poly-A][rc(CBC)][rc(UMI)][adapter]
# so a read in the OPPOSITE orientation reads its reverse complement:
#     [rc(adapter)][UMI][CBC][poly-T][rc(insert)]
# which is a 5' problem, not a 3' one. cutadapt's -g removes a 5' adapter and
# everything BEFORE it, so the same anchor trick applies mirrored -- except the
# UMI is not constant within a cell, so the wildcards sit in the middle again
# and min_overlap must again be the full length.
#
# Measured before building it: UMI+CBC appears in the first 60 nt of 0.40% of
# reads in cell 011 and 3.28% in the blank 016. So this is worth about 1,200
# reads in 300,000 for a real cell. v24 tests whether they are recoverable.
#
#   v24  v23 + 5' anchor + 5' poly-T
#
#   sbatch bench_trim8.sh
###############################################################################
#SBATCH --job-name=trimbench8
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim8.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 016"
EBROOT=/camp/apps/eb/software
IDX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116

AD=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT
RT="rt=${AD};min_overlap=8"
PA="polyA=A{20};min_overlap=10"
PG="polyG=G{20};min_overlap=10"
declare -A BC=( [011]=GTGACA [016]=TGTCGA )
revcomp(){ echo "$1" | tr ACGT TGCA | rev; }

source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
mkdir -p v24 bam

for c in $CELLS; do
  [ -s "v24/${c}.fq.gz" ] && continue
  fwd="bcumi=$(revcomp "${BC[$c]}")NNNNNN${AD:0:16};min_overlap=28"
  rev="revbcumi=$(revcomp "${AD:0:16}")NNNNNN${BC[$c]};min_overlap=28"
  echo "[$(date +%H:%M:%S)] v24 $c"
  echo "   -a $fwd"
  echo "   -g $rev"
  cutadapt -j 4 -m 20 --trim-n -n 3 \
    -a "$fwd" -a "$RT" -a "$PA" -a "$PG" --poly-a \
    -g "$rev" -g "polyT5=T{20};min_overlap=10" \
    -o "v24/${c}.fq.gz" "tg/${c}_trimmed.fq.gz" > "v24/${c}.cutadapt.log" 2>&1
done

source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
module load STAR/2.7.7a-GCC-10.2.0
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
trap '$STAR --genomeDir $IDX --genomeLoad Remove --outFileNamePrefix ${tmp}/rm_ >/dev/null 2>&1; rm -rf $tmp' EXIT
$STAR --genomeDir "$IDX" --genomeLoad LoadAndExit --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1
for c in $CELLS; do
  [ -s "bam/v24_${c}_Aligned.out.bam" ] && continue
  $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
    --readFilesIn "v24/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
    --outSAMtype BAM Unsorted --outSAMattributes All \
    --outFileNamePrefix "bam/v24_${c}_" >/dev/null 2>&1
  cp "bam/v24_${c}_Log.final.out" "v24/${c}_Log.final.out"
  rm -rf "bam/v24_${c}__STARtmp" "bam/v24_${c}_Log.progress.out" "bam/v24_${c}_Log.out"
done
echo "bench8 done"
