#!/bin/bash
#SBATCH -J trimbench12
#SBATCH -p ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim12.%j.out
###############################################################################
# bench_trim12.sh -- is --nextseq-trim the right way to handle poly-G here?
#
# THE QUESTION
# ------------
# pass 2 currently removes poly-G with `-a polyG=G{20};min_overlap=10`, i.e. as
# a 3' adapter. cutadapt documents a dedicated option for this instead:
#
#   --nextseq-trim 3'CUTOFF   "works like regular quality trimming ... except
#                              that the qualities of G bases are ignored"
#
# That matters because this is a two-colour instrument. The delivery flowcell is
# LH00442 = NovaSeq X Plus, where "basecalls without any signal are called as
# high-quality G bases" (trim_galore --2colour help text). Measured on 200,000
# real reads of cell 011:
#
#   0.217% of reads end in a run of >=10 G
#   75% of those runs have mean Phred >= 25
#
# High Phred is the whole problem: pass 1's `-q 20` cannot remove them because
# they do not look low-quality, and `-a G{20};min_overlap=10` only fires when the
# run is >=10 nt AND is a clean 3' suffix.
#
# WHY IT IS NOT OBVIOUSLY AN IMPROVEMENT, HENCE THIS BENCHMARK
# ------------------------------------------------------------
# Switching pass 2 to --nextseq-trim=20 drops ~1,970 more reads per 200,000
# (1.7%) than the current setting. Those reads are NOT mostly poly-G: of the
# 2,080 reads kept by the current setting and dropped by --nextseq-trim, 4% are
# >=50% G and none are >=80% G (mean length 72.6 nt, a mix of poly-C, adapter
# read-through, and reads that look real). So the question "does this remove
# junk or real data?" cannot be answered from base composition, exactly as in
# round 7. It is answered the same way that one was: by where the surviving
# reads land, via annot_fraction.sh.
#
# PLACEMENT IS ALSO A QUESTION
# ----------------------------
# --nextseq-trim is quality trimming, and quality trimming belongs BEFORE
# adapter removal. In pass 2 it runs after pass 1 has already trimmed. So G3
# tests the other placement: trim_galore --2colour 20 in pass 1. Note that flag
# is mutually exclusive with -q, so G3 REPLACES pass 1's -q 20 rather than
# adding to it. cutadapt 1.18 (what trim_galore drives) was verified to accept
# --nextseq-trim before writing this.
#
# THE VARIANTS
# ------------
#   G0  current production, unchanged                      <- baseline
#   G1  G0 + pass 2 --nextseq-trim=20
#   G2  G1 with -a polyG removed  (once --nextseq-trim is on, polyG changes
#                                  only 105 reads per 200k = 0.09%)
#   G3  pass 1 trim_galore --2colour 20 instead of -q 20, pass 2 unchanged
#
# Everything else is production as of 2026-08-03, in particular -n 10 (verified
# again here: -n 3 leaves 5,376 >=90%-T reads in blank cell 016, -n 10 leaves 0,
# -n 20 adds nothing, and -n 10 is faster than -n 3).
#
# READ THE SCORE, NOT THE COUNT: uniquely-mapped is gameable by junk. The
# decision metric is protein-coding exonic reads and the in-annotation fraction.
# A setting whose extra reads are real holds or raises those; one whose extra
# reads are junk dilutes them.
#
#   sbatch bench_trim12.sh
#   ./annot_fraction.sh G0 G1 G2 G3          # then, for both cells
#   CELL=016 ./annot_fraction.sh G0 G1 G2 G3
###############################################################################
set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
OWN=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/own_version
cd "$ROOT" || exit 1

CELLS="011 016"
VARIANTS="G0 G1 G2 G3"
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

for c in $CELLS; do
    [ -s "pass0/${c}.fastq.gz" ] || { echo "run bench_trim10.sh first (need pass0/${c}.fastq.gz)"; exit 1; }
done

# --- pass 1 -----------------------------------------------------------------
# Two flavours only: stock (G0/G1/G2 share it) and --2colour (G3).
( source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
  module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
  for flavour in stock 2colour; do
    d="tg_${flavour}"; mkdir -p "$d"
    for c in $CELLS; do
      [ -s "${d}/${c}_tg.fq.gz" ] && continue
      say "pass1 ${flavour} ${c}"
      if [ "$flavour" = "2colour" ]; then
        # --2colour is mutually exclusive with -q, so this REPLACES -q 20
        $TG --path_to_cutadapt "$CA118" --cores 4 --length 20 --2colour 20 \
            -o "$d" "pass0/${c}.fastq.gz" >/dev/null 2>&1
      else
        $TG --path_to_cutadapt "$CA118" --cores 4 --length 20 \
            -o "$d" "pass0/${c}.fastq.gz" >/dev/null 2>&1
      fi
      mv "${d}/${c}_trimmed.fq.gz" "${d}/${c}_tg.fq.gz" 2>/dev/null || true
    done
  done )

# --- pass 2 -----------------------------------------------------------------
( source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
  conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
  for v in $VARIANTS; do
    mkdir -p "$v"
    case $v in
      G0) src=tg_stock;   extra=(-a "$PG") ;;
      G1) src=tg_stock;   extra=(-a "$PG" --nextseq-trim=20) ;;
      G2) src=tg_stock;   extra=(--nextseq-trim=20) ;;
      G3) src=tg_2colour; extra=(-a "$PG") ;;
    esac
    for c in $CELLS; do
      [ -s "${v}/${c}.fq.gz" ] && continue
      say "pass2 $v $c"
      cutadapt -j 4 -m 20 --trim-n -n 10 --poly-a \
        -a "$RT" -a "$PA" "${extra[@]}" -g "$PT" \
        --json "${v}/${c}.cutadapt.json" \
        -o "${v}/${c}.fq.gz" "${src}/${c}_tg.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
    done
  done )

# --- map --------------------------------------------------------------------
( source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
  module load STAR/2.7.7a-GCC-10.2.0
  STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin/STAR
  tmp=$(mktemp -d "${ROOT}/.star.XXXXXX")
  trap '$STAR --genomeDir '"$IDX"' --genomeLoad Remove --outFileNamePrefix ${tmp}/rm_ >/dev/null 2>&1; rm -rf $tmp' EXIT
  $STAR --genomeDir "$IDX" --genomeLoad LoadAndExit --outFileNamePrefix "${tmp}/load_" >/dev/null || exit 1
  for v in $VARIANTS; do
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

echo
echo "bench12 done -- score with:"
echo "  ${OWN}/trimtest/annot_fraction.sh G0 G1 G2 G3"
echo "  CELL=016 ${OWN}/trimtest/annot_fraction.sh G0 G1 G2 G3"
