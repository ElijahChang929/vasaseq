#!/bin/bash
###############################################################################
# bench_trim10.sh -- round 10. The simple version of the barcode anchor.
#
# Rounds 7-9 did the anchor inside cutadapt, which cannot take a per-read
# pattern. That forced the UMI to be 6 wildcards, wildcards forced
# min_overlap = the pattern's full length, and that in turn forced ~21 nt of
# adapter into the pattern -- so it only fired on reads that had sequenced that
# much adapter. All of that machinery exists to work around one tool limitation.
#
# ../trim_bc_anchor.py does it directly: the UMI is on the read name, so the
# 12 nt after the poly-A are a literal string for THIS read. Find it, drop
# everything from there to the 3' end, walk back over the poly-A. No adapter
# needed, so it also fires on reads whose tail is truncated -- and it can match
# a partial anchor at the read end, which the cutadapt version structurally
# cannot.
#
#   L21  incumbent: anchor inside cutadapt (rounds 7-9)
#   v26  trim_bc_anchor.py -> TrimGalore -> cutadapt WITHOUT the bcumi pattern
#   v27  as v26 but keeping bcumi too, to see whether it still adds anything
#
# Note v26/v27 re-run TrimGalore, because the pre-pass changes its input. That
# is the point of the ordering -- pass 0, then the normal chain.
#
#   sbatch bench_trim10.sh
###############################################################################
#SBATCH --job-name=trimbench10
#SBATCH --partition=ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim10.%j.out

set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
OWN=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/own_version
cd "$ROOT" || exit 1

CELLS="011 016"
EBROOT=/camp/apps/eb/software
IDX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
TG=${EBROOT}/Trim_Galore/0.6.2-foss-2018b-Python-3.6.6/trim_galore
CA118=${EBROOT}/cutadapt/1.18-foss-2018b-Python-3.6.6/bin/cutadapt

AD=GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT
RT="rt=${AD};min_overlap=8"
PA="polyA=A{20};min_overlap=10"
PG="polyG=G{20};min_overlap=10"
PT="polyT5=T{20};min_overlap=10"
declare -A BC=( [011]=GTGACA [016]=TGTCGA )
revcomp(){ echo "$1" | tr ACGT TGCA | rev; }

say(){ echo "[$(date +%H:%M:%S)] $*"; }

mkdir -p pass0 tg0 v26 v27 bam

# --- pass 0: the python anchor, on the RAW demultiplexed reads ---------------
( source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
  conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
  for c in $CELLS; do
    [ -s "pass0/${c}.fastq.gz" ] && continue
    say "pass0 $c"
    /usr/bin/time -f "  %e s  %M KB" python "${OWN}/trim_bc_anchor.py" \
      "in/${c}.fastq.gz" "pass0/${c}.fastq.gz" --log "pass0/${c}.log" 2>&1 | tee -a "pass0/${c}.stderr"
  done )

# --- TrimGalore on the pre-passed reads --------------------------------------
( source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
  module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
  for c in $CELLS; do
    [ -s "tg0/${c}_trimmed.fq.gz" ] && continue
    say "trim_galore(pass0) $c"
    $TG --path_to_cutadapt "$CA118" --cores 4 -o tg0 "pass0/${c}.fastq.gz" >/dev/null 2>&1
  done )

# --- the regular cutadapt pass ----------------------------------------------
( source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
  conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
  for c in $CELLS; do
    bcumi="bcumi=$(revcomp "${BC[$c]}")NNNNNN${AD:0:21}"
    [ -s "v26/${c}.fq.gz" ] || { say "v26 $c"
      cutadapt -j 4 -m 20 --trim-n -n 3 -a "$RT" -a "$PA" -a "$PG" --poly-a -g "$PT" \
        -o "v26/${c}.fq.gz" "tg0/${c}_trimmed.fq.gz" > "v26/${c}.cutadapt.log" 2>&1; }
    [ -s "v27/${c}.fq.gz" ] || { say "v27 $c"
      cutadapt -j 4 -m 20 --trim-n -n 3 -a "${bcumi};min_overlap=33" \
        -a "$RT" -a "$PA" -a "$PG" --poly-a -g "$PT" \
        -o "v27/${c}.fq.gz" "tg0/${c}_trimmed.fq.gz" > "v27/${c}.cutadapt.log" 2>&1; }
  done )

# --- map --------------------------------------------------------------------
( source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
  module load STAR/2.7.7a-GCC-10.2.0
  STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
  tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
  trap '$STAR --genomeDir '"$IDX"' --genomeLoad Remove --outFileNamePrefix ${tmp}/rm_ >/dev/null 2>&1; rm -rf $tmp' EXIT
  $STAR --genomeDir "$IDX" --genomeLoad LoadAndExit --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1
  for v in v26 v27; do
    for c in $CELLS; do
      [ -s "bam/${v}_${c}_Aligned.out.bam" ] && continue
      say "STAR $v $c"
      $STAR --runThreadN 16 --genomeDir "$IDX" --genomeLoad LoadAndKeep \
        --readFilesIn "${v}/${c}.fq.gz" --readFilesCommand zcat --outFilterMultimapNmax 20 \
        --outSAMtype BAM Unsorted --outSAMattributes All \
        --outFileNamePrefix "bam/${v}_${c}_" >/dev/null 2>&1
      cp "bam/${v}_${c}_Log.final.out" "${v}/${c}_Log.final.out"
      rm -rf "bam/${v}_${c}__STARtmp" "bam/${v}_${c}_Log.progress.out" "bam/${v}_${c}_Log.out"
    done
  done )
echo "bench10 done"
