#!/bin/bash
# Rebuild every four-way table and figure, in order.
#
#   ./run.sh          steps 5-6  (tables from the scans, then every figure)
#   ./run.sh units    step 5 only
#   ./run.sh plots    step 6 only    <- the usual edit-rerun loop
#
# Steps 0-4 read FASTQ/BAM/BED and are sbatch one-offs, so they are NOT in here.
# Step 0 also has to run BEFORE FLASH-seq is mapped at all -- it is the stage
# pipeline_fs.sh was missing. The whole FLASH-seq chain, once:
#
#   # 1. trim, in flashseq_vasa/, writing to the lab share (scratch gets purged)
#   R=/nemo/lab/turnerj/working/guangxin/vasaseq/data/flashseq_vasa/run
#   sbatch --chdir=<flashseq_vasa> -c 32 --mem=32G -t 240 \
#     --wrap="FSV_SCRATCH=$R FSV_ARM=native FSV_NCORES=8 ./pipeline_fs.sh prep"
#
#   # 2. THE MISSING STAGE -- rRNA depletion, VASA's own step 3
#   sbatch -c 32 --mem=48G -t 8:00:00 --wrap="NPAR=4 scripts/00_fs_ribo.sh"
#
#   # 3. map + assign off the depleted reads. FSV_SAMPLE is set explicitly:
#   #    its default FS_$FSV_ARM would overwrite the existing FS_native_* tables.
#   E="FSV_OUTDIR=$R/nonribo FSV_ARM=native FSV_SAMPLE=FS_nonribo"
#   sbatch --chdir=<flashseq_vasa> -c 16 --mem=64G -t 8:00:00 --wrap="$E ./pipeline_fs.sh map"
#   sbatch --chdir=<flashseq_vasa> -c 10 --mem=32G -t 8:00:00 --wrap="$E FSV_NCORES=5 ./pipeline_fs.sh assign"
#
#   # 4. steps 6-7, the count tables
#   #    pickle1 is per library because step 6's cost is linear in BED size
#   for L in ZHA8833A1 ... ; do sbatch ... --wrap="$E FSV_LIBS=$L ./pipeline_fs.sh pickle1"; done
#   sbatch ... --wrap="$E ./pipeline_fs.sh pickle_merge tables recon"
#
# Then the four-way scans, which read all four datasets at once:
#
#   scripts/01_insilico_depletion.sh                         # seconds, login node
#   sbatch -c 8  --mem=8G  -t 4:00:00 --wrap="scripts/02_probe_qc.sh"
#   sbatch -c 16 --mem=8G  -t 120    --wrap="scripts/03_mapped_length_dist.sh"
#   sbatch -c 16 --mem=64G -t 240    --wrap="scripts/04_step5_biotype.sh"
#   sbatch -c 4  --mem=8G  -t 240    --wrap="scripts/07_genebody_coverage.py"
#     measured 29m10s, peak 190 MB -- it streams, so memory holds the gene models
#     and 100 counters per gene and does not grow with read count. The 32G this
#     was first submitted with was a guess; 8G is still ~40x the measurement.
#   sbatch -c 4  --mem=16G -t 240    --wrap="scripts/09_saturation.py"
#
# THE SAME FOLDER, AGAINST THE E116-REMAPPED PLATE
# -------------------------------------------------
# scripts/00b_plate_e116.sh puts the published plate's 173 mouse wells on the
# same reference as the other three. Nothing here is edited to use it -- two
# env vars move the inputs and the products:
#
#   export PLATE_CELLS=<...>/plate_e116/cells
#   export FOURWAY_OUT=<...>/fourway_analyze_e116
#   export MOUSECONTIG=''      # REQUIRED -- see below
#   scripts/01_insilico_depletion.sh; sbatch ... scripts/0{2,3,4}_*; ...
#   FOURWAY_OUT=$FOURWAY_OUT PLATE_CELLS=$PLATE_CELLS ./run.sh
#
# Scripts 01 and 02 score rRNA depletion, which is stage 3 and identical under
# either annotation; 00b symlinks the stage-3 artefacts into the new cells dir
# so they run unchanged rather than being special-cased.
#
# MOUSECONTIG='' is not optional. Script 04 restricts plate rows to `^GRCm38_`
# because the E99 plate is on a mixed reference; the E116 plate is mouse-only
# and its contigs are plain `1`, `2`, ..., so the default filter matches nothing
# and discards all 173 wells. 04 now refuses to finish when a whole dataset ends
# with zero BED reads, but the variable still has to be set.
#
# Layout: scripts/ code (numbered by run order), tables/ TSVs,
# figures/<step>/ PNG+PDF. Everything below is seconds on the login node.
set -euo pipefail
R=/nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/Rscript
cd "$(dirname "$(readlink -f "$0")")"
S=scripts

step() { printf '\n=== [%s] %s ===\n' "$1" "$2"; }
units() { step 5 "build_units.R -- one row per unit, four datasets"; "$R" $S/05_build_units.R; }
plots() { step 6  "plot_fourway.R    -- every cross-dataset figure"; "$R" $S/06_plot_fourway.R
          step 8  "plot_genebody.R   -- gene body coverage 5'->3'";  "$R" $S/08_plot_genebody.R
          step 10 "plot_saturation.R -- genes vs depth";             "$R" $S/10_plot_saturation.R; }

case "${1:-all}" in
  units) units ;;
  plots) plots ;;
  all)   units; plots ;;
  *) echo "usage: $0 [all|units|plots]" >&2; exit 2 ;;
esac
printf '\ndone.\n'
