#!/bin/bash
#SBATCH -J trimbench13
#SBATCH -p ncpu
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/bench_trim13.%j.out
###############################################################################
# bench_trim13.sh -- does pass 0 need its own poly-A trimmer at all?
#
# THE QUESTION
# ------------
# trim_bc_anchor.py (pass 0) cuts at the barcode anchor and then calls its own
# strip_polya() to walk back over the poly-A. Pass 2 then runs cutadapt with
# --poly-a, which does the same job. The two overlap, and pass 0's copy is the
# weaker of the two.
#
# WHY IT IS WEAKER -- ONE MISSING LINE
# ------------------------------------
# cutadapt's documented algorithm has two parts:
#
#   (a) score each suffix +1 per A, -2 per non-A, take the max
#   (b) "exclude all suffixes from consideration that have more than 20% non-A"
#
# strip_polya() implements (a) and NOT (b). The scoring alone only enforces
# <33% non-A (a suffix of a A's and b non-A's scores a-2b, positive iff
# b/(a+b) < 1/3), so between 20% and 33% non-A the two disagree, and pass 0
# eats sequence cutadapt would keep.
#
# MEASURED, on 400,000 reads of cell 011 (118,163 anchor hits):
#
#   identical            118,057   99.910%
#   disagree                 106    0.090%   <- ALL of them pass 0 over-trimming
#   mean over-trim        24.0 nt            <- worst cases 89->3, 126->49, 125->59
#
# It is not dead code: strip_polya fires on ~97% of anchor hits (24.5% of all
# reads library-wide), so this is a real, if rare, defect on a hot path.
#
# END-TO-END, pass 0 -> TrimGalore -> pass 2, same 400,000 reads:
#
#   strip_polya ON     180,679 reads   18,208,958 bases
#   strip_polya OFF    181,175 reads   18,232,101 bases   (+496 reads)
#   548 reads recovered, 0% of them >=80% A -- real sequence, median 20 nt
#   of the 180,627 shared reads, 94.98% byte-identical; of the 9,067 that
#   differ, 8,325 (91.8%) are LONGER with it off
#
# WHY THAT IS NOT YET AN ANSWER
# -----------------------------
# The recovered reads have median length 20 nt, i.e. exactly at the -m 20 floor,
# and round 11 established that short reads are precisely the junk that inflates
# a raw count. "More reads" and "longer reads" are not the metric. annot_fraction
# is. Round 12 just made the same point in the other direction: --nextseq-trim
# was the documented, plausible change and still lost exonic reads.
#
# THE VARIANTS
# ------------
#   S0  current: pass 0 with strip_polya                  <- baseline
#   S1  pass 0 with strip_polya REMOVED (pass 2 does it)
#   S2  pass 0 with strip_polya FIXED (20% non-A guard added), kept
#
# S2 separates the two hypotheses: if S2 ~= S1 the function is simply
# redundant and should go; if S2 > S1 the walk-back is doing something pass 2
# cannot (e.g. shortening reads before TrimGalore's --length 20 sees them) and
# should be kept, fixed.
#
#   sbatch bench_trim13.sh
#   ./annot_fraction.sh S0 S1 S2
#   CELL=016 ./annot_fraction.sh S0 S1 S2
###############################################################################
set -uo pipefail
ROOT=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest
OWN=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/own_version
cd "$ROOT" || exit 1

CELLS="011 016"
VARIANTS="S0 S1 S2"
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
mkdir -p bam s13

# --- build the two pass-0 variants by patching a copy ------------------------
# Copies, never edits in place: trim_bc_anchor.py stays production.
cp "${OWN}/trim_bc_anchor.py" s13/anchor_S0.py

# S1: skip the call entirely
sed 's|^                    t = strip_polya(seq)$|                    t = seq  # S1: strip_polya removed|' \
    "${OWN}/trim_bc_anchor.py" > s13/anchor_S1.py

# S2: add cutadapt's 20% non-A guard to the accept condition
python3 - <<'PY'
src = open("/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/own_version/trim_bc_anchor.py").read()
old = """    best = score = 0
    cut = len(seq)
    for i in range(len(seq) - 1, -1, -1):
        score += 1 if seq[i] == "A" else -2
        if score > best:
            best, cut = score, i"""
new = """    best = score = 0
    cut = len(seq)
    na = nn = 0
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == "A":
            na += 1; score += 1
        else:
            nn += 1; score -= 2
        # S2: cutadapt's second rule -- a suffix with >20% non-A is not a tail
        if score > best and nn <= 0.2 * (na + nn):
            best, cut = score, i"""
assert old in src, "strip_polya body not found -- did the function change?"
open("s13/anchor_S2.py", "w").write(src.replace(old, new))
PY
chmod +x s13/anchor_S*.py
for v in S0 S1 S2; do [ -s "s13/anchor_${v}.py" ] || { echo "failed to build s13/anchor_${v}.py"; exit 1; }; done

# --- pass 0, one per variant -------------------------------------------------
( source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
  conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
  for v in $VARIANTS; do
    mkdir -p "s13/p0_${v}"
    for c in $CELLS; do
      [ -s "s13/p0_${v}/${c}.fastq.gz" ] && continue
      say "pass0 $v $c"
      python "s13/anchor_${v}.py" "in/${c}.fastq.gz" \
        "s13/p0_${v}/${c}.fastq.gz" --log "s13/p0_${v}/${c}.log" 2>/dev/null
    done
  done )

# --- pass 1 + pass 2, identical for every variant ----------------------------
( source /usr/share/lmod/lmod/init/bash; export MODULEPATH=/camp/apps/eb/modules/all
  module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
  for v in $VARIANTS; do
    mkdir -p "$v"
    for c in $CELLS; do
      [ -s "${v}/${c}_tg.fq.gz" ] && continue
      say "pass1 $v $c"
      $TG --path_to_cutadapt "$CA118" --cores 4 --length 20 \
          -o "$v" "s13/p0_${v}/${c}.fastq.gz" >/dev/null 2>&1
      mv "${v}/${c}_trimmed.fq.gz" "${v}/${c}_tg.fq.gz" 2>/dev/null || true
    done
  done )

( source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
  conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa
  for v in $VARIANTS; do
    for c in $CELLS; do
      [ -s "${v}/${c}.fq.gz" ] && continue
      say "pass2 $v $c"
      cutadapt -j 4 -m 20 --trim-n -n 10 --poly-a \
        -a "$RT" -a "$PA" -a "$PG" -g "$PT" \
        --json "${v}/${c}.cutadapt.json" \
        -o "${v}/${c}.fq.gz" "${v}/${c}_tg.fq.gz" > "${v}/${c}.cutadapt.log" 2>&1
    done
  done )

# --- map ---------------------------------------------------------------------
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
echo "bench13 done -- score with:"
echo "  ${OWN}/trimtest/annot_fraction.sh S0 S1 S2"
echo "  CELL=016 ${OWN}/trimtest/annot_fraction.sh S0 S1 S2"
