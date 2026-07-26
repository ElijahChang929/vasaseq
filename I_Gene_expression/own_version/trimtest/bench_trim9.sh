#!/bin/bash
###############################################################################
# bench_trim9.sh -- round 9. How much of the adapter should the barcode anchor
# carry?
#
# The anchor is  revcomp(CBC) + N x LEN_UMI + <first L nt of the adapter>  at
# min_overlap = its own full length. L was originally 16, which was arbitrary:
# it was the 16-mer used as a search anchor when the read-through consensus was
# derived, reused without justification. The principled value is 21, because
#
#     TRIM_ADAPTER3[:21] == revcomp(the 21 nt prefix stripped from R1)
#
# exactly, and everything past 21 is the Nextera mosaic end.
#
# But longer is not automatically better, and the reason is min_overlap. The
# whole pattern must fit INSIDE the read for a full-length match, so a bigger L
# means a read has to have sequenced more adapter before the anchor will fire.
# Specificity is not the binding constraint -- even L=8 leaves 14 specific
# bases, and 4^-14 over 300k x 130 positions is a fraction of one expected
# chance hit -- so the sweep is really about sensitivity.
#
#   L in 8 12 16 21 26      (anchor length = L + 12)
#
#   sbatch bench_trim9.sh
###############################################################################
#SBATCH --job-name=trimbench9
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim9.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 016"
LENS="8 12 16 21 26"
EBROOT=/camp/apps/eb/software
IDX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116

AD=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT
RT="rt=${AD};min_overlap=8"
PA="polyA=A{20};min_overlap=10"
PG="polyG=G{20};min_overlap=10"
PT="polyT5=T{20};min_overlap=10"
declare -A BC=( [011]=GTGACA [016]=TGTCGA )
revcomp(){ echo "$1" | tr ACGT TGCA | rev; }

source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
mkdir -p bam

for L in $LENS; do
  v="L${L}"; mkdir -p "$v"
  for c in $CELLS; do
    [ -s "${v}/${c}.fq.gz" ] && continue
    anchor="$(revcomp "${BC[$c]}")NNNNNN${AD:0:$L}"
    echo "[$(date +%H:%M:%S)] $v $c  anchor=${anchor} (len ${#anchor})"
    cutadapt -j 4 -m 20 --trim-n -n 3 \
      -a "bcumi=${anchor};min_overlap=${#anchor}" \
      -a "$RT" -a "$PA" -a "$PG" --poly-a -g "$PT" \
      -o "${v}/${c}.fq.gz" "tg/${c}_trimmed.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
  done
done

source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
module load STAR/2.7.7a-GCC-10.2.0
STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
trap '$STAR --genomeDir $IDX --genomeLoad Remove --outFileNamePrefix ${tmp}/rm_ >/dev/null 2>&1; rm -rf $tmp' EXIT
$STAR --genomeDir "$IDX" --genomeLoad LoadAndExit --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1
for L in $LENS; do
  v="L${L}"
  for c in $CELLS; do
    [ -s "bam/${v}_${c}_Aligned.out.bam" ] && continue
    echo "[$(date +%H:%M:%S)] STAR $v $c"
    $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
      --readFilesIn "${v}/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
      --outSAMtype BAM Unsorted --outSAMattributes All \
      --outFileNamePrefix "bam/${v}_${c}_" >/dev/null 2>&1
    cp "bam/${v}_${c}_Log.final.out" "${v}/${c}_Log.final.out"
    rm -rf "bam/${v}_${c}__STARtmp" "bam/${v}_${c}_Log.progress.out" "bam/${v}_${c}_Log.out"
  done
done
echo "bench9 done"
