#!/bin/bash
#SBATCH -J trimbench11
#SBATCH -p ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim11.%j.out
###############################################################################
# bench_trim11.sh -- how short a read is still worth keeping?
#
# TRIM_MINLEN is 20. This library's median insert before the poly-A is 10 nt, so
# the floor is not a detail here: it is the single setting that decides how much
# of the library survives at all.
#
# READ THIS BEFORE INTERPRETING ANY EARLIER RESULT: until trim.sh was fixed
# (2026-07-26), TRIM_MINLEN reached only cutadapt. trim_galore hardcodes
# $length_cutoff = 20 when --length is not passed, and trim.sh never passed it,
# so pass 1 had ALREADY deleted every read below 20 before pass 2's -m saw
# anything. Lowering TRIM_MINLEN did nothing. This sweep is only meaningful with
# that fix in place -- it varies --length and -m together.
#
# Why short reads are suspect, and why it has to be measured rather than
# reasoned: a mammalian genome is ~2.7 Gb, so a random k-mer is expected to
# occur 2.7e9 / 4^k times -- 0.63 times at k=16, 0.04 at k=18, 0.0025 at k=20.
# So 20 is roughly where uniqueness starts, and that is the naive calculation
# that ignores repeats, which make it worse. STAR does not save you either:
# outFilterMatchNminOverLread is RELATIVE (0.66), so a 16 nt read needs only 11
# matching bases to be reported.
#
# The score is protein-coding exonic reads, via annot_fraction.sh, for the
# reason established in round 7: raw uniquely-mapped counts are gameable by
# junk, and short reads are exactly the junk that games them.
#
#   sbatch bench_trim11.sh
###############################################################################
set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
OWN=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/own_version
cd "$ROOT" || exit 1

CELLS="011 016"
LENS="12 15 18 20 25 30"
EBROOT=/camp/apps/eb/software
IDX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
TG=${EBROOT}/Trim_Galore/0.6.2-foss-2018b-Python-3.6.6/trim_galore
CA118=${EBROOT}/cutadapt/1.18-foss-2018b-Python-3.6.6/bin/cutadapt

AD=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT
RT="rt=${AD};min_overlap=8"
PA="polyA=A{20};min_overlap=10"
PG="polyG=G{20};min_overlap=10"
PT="polyT5=T{20};min_overlap=10"

say(){ echo "[$(date +%H:%M:%S)] $*"; }
mkdir -p bam

# pass 0 is independent of the length floor -- reuse bench_trim10's output
for c in $CELLS; do
    [ -s "pass0/${c}.fastq.gz" ] || { echo "run bench_trim10.sh first (need pass0/${c}.fastq.gz)"; exit 1; }
done

for L in $LENS; do
  v="M${L}"; mkdir -p "$v"
  # pass 1 with the floor actually applied
  ( source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
    module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
    for c in $CELLS; do
      [ -s "${v}/${c}_trimmed.fq.gz" ] && continue
      say "$v trim_galore --length $L $c"
      $TG --path_to_cutadapt "$CA118" --cores 4 --length "$L" -o "$v" "pass0/${c}.fastq.gz" >/dev/null 2>&1
      mv "${v}/${c}_trimmed.fq.gz" "${v}/${c}_tg.fq.gz" 2>/dev/null || true
    done )
  # pass 2 with the same floor
  ( source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
    conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
    for c in $CELLS; do
      [ -s "${v}/${c}.fq.gz" ] && continue
      say "$v cutadapt -m $L $c"
      cutadapt -j 4 -m "$L" --trim-n -n 3 --poly-a \
        -a "$RT" -a "$PA" -a "$PG" -g "$PT" \
        -o "${v}/${c}.fq.gz" "${v}/${c}_tg.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
    done )
done

( source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
  module load STAR/2.7.7a-GCC-10.2.0
  STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
  tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
  trap '$STAR --genomeDir '"$IDX"' --genomeLoad Remove --outFileNamePrefix ${tmp}/rm_ >/dev/null 2>&1; rm -rf $tmp' EXIT
  $STAR --genomeDir "$IDX" --genomeLoad LoadAndExit --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1
  for L in $LENS; do
    v="M${L}"
    for c in $CELLS; do
      [ -s "bam/${v}_${c}_Aligned.out.bam" ] && continue
      say "STAR $v $c"
      $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
        --readFilesIn "${v}/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
        --outSAMtype BAM Unsorted --outSAMattributes All \
        --outFileNamePrefix "bam/${v}_${c}_" >/dev/null 2>&1
      cp "bam/${v}_${c}_Log.final.out" "${v}/${c}_Log.final.out" 2>/dev/null
      rm -rf "bam/${v}_${c}__STARtmp" "bam/${v}_${c}_Log.progress.out" "bam/${v}_${c}_Log.out"
    done
  done )
echo "bench11 done -- score with: ${OWN}/trimtest/annot_fraction.sh M12 M15 M18 M20 M25 M30"
