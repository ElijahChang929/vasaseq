#!/bin/bash
###############################################################################
# bench_trim7.sh -- round 7. Anchor the trim on the cell barcode itself.
#
# Rounds 1-6 kept trying to recognise the poly-A tail by its own shape, and
# every attempt either missed short tails or truncated real mRNA at an internal
# A-run. The tail does not have to be recognised by shape: we KNOW what follows
# it. concatenator.py put the cell barcode and UMI on the read name at step 1,
# and the 3' construct is
#
#     [insert][poly-A][revcomp(CBC)][revcomp(UMI)][revcomp(R1 prefix)+Nextera]
#
# Within one cell's fastq the CBC is CONSTANT, so per cell the whole tail is a
# fixed 28 nt pattern with only the 6 nt UMI unknown:
#
#     bcumi = revcomp(CBC) NNNNNN <first 16 nt of the adapter>
#              6 specific   6 any      16 specific      = 22 specific of 28
#
# Requiring min_overlap = 28 (the full length) forbids partial 3'-end matching,
# which is what made the round-2 wildcard adapter promiscuous. 22 specific
# bases make a chance hit impossible.
#
# The point is what this unlocks: once the anchor is removed, the poly-A is at
# the very 3' END of the read, which is the only place --poly-a looks. Round 4
# had to reject --poly-a because after plain adapter removal the 3' end was the
# 12 nt barcode remnant and --poly-a could not reliably see past it. With the
# anchor gone it can, and it removes tails of ANY length -- including the runs
# shorter than 20 nt that A{20} structurally cannot catch.
#
#   v21  bcumi + rt + poly-G + --poly-a          (no A{20} -- --poly-a replaces it)
#   v22  bcumi + v17 as-is                       (anchor only, keep A{20})
#   v23  bcumi + v17 + --poly-a                  (both)
# v17 is the incumbent to beat.
#
# Measured anchor availability on cell 011, 300k reads: 23.2% carry the exact
# 12 nt rc(CBC)+rc(UMI), a further 8.6% carry rc(CBC) alone, 0.4% are in the
# reverse orientation (3.3% in the blank 016 -- not worth handling), and 67.8%
# carry no anchor at all, mostly because the insert is long enough that the
# read never reaches the tail.
#
#   sbatch bench_trim7.sh
###############################################################################
#SBATCH --job-name=trimbench7
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim7.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
cd "$ROOT" || exit 1

CELLS="011 007 001 016"
EBROOT=/camp/apps/eb/software
IDX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
VARIANTS="v21 v22 v23"

AD=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT
RT="rt=${AD};min_overlap=8"
PA="polyA=A{20};min_overlap=10"
PG="polyG=G{20};min_overlap=10"

# cell id -> barcode, from ../bc_PM26037_6nt.tsv
declare -A BC=( [001]=ACTCGA [007]=CATGAG [011]=GTGACA [016]=TGTCGA )
revcomp(){ echo "$1" | tr ACGT TGCA | rev; }

say(){ echo "[$(date +%H:%M:%S)] $*"; }

source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
mkdir -p v21 v22 v23 bam

for v in $VARIANTS; do
  for c in $CELLS; do
    [ -s "${v}/${c}.fq.gz" ] && continue
    bcumi="bcumi=$(revcomp "${BC[$c]}")NNNNNN${AD:0:16};min_overlap=28"
    say "$v $c   ($bcumi)"
    case $v in
      v21) args=(-n 3 -a "$bcumi" -a "$RT" -a "$PG" --poly-a) ;;
      v22) args=(-n 3 -a "$bcumi" -a "$RT" -a "$PA" -a "$PG") ;;
      v23) args=(-n 3 -a "$bcumi" -a "$RT" -a "$PA" -a "$PG" --poly-a) ;;
    esac
    cutadapt -j 4 -m 20 --trim-n "${args[@]}" \
      -o "${v}/${c}.fq.gz" "tg/${c}_trimmed.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
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
echo "bench7 done"
