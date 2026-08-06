#!/bin/bash
# Rebuild every table and figure, in order.
#
#   ./run.sh          steps 3-13
#   ./run.sh tables   steps 3-4 only
#   ./run.sh plots    steps 6-13 only   <- the usual edit-rerun loop
#
# Steps 1, 2, 5, 10 and 12 read FASTQ/BAM/BED and are sbatch one-offs, so they
# are NOT in here. Run them once, by hand:
#
#   sbatch -c 8  --mem=8G  -t 60      --wrap="scripts/01_count_demux_reads.sh"
#   sbatch -c 16 --mem=8G  -t 120     --wrap="scripts/02_count_demux_reads_plate.sh"
#   sbatch -c 4  --mem=8G  -t 0:40:00 --wrap="scripts/05_probe_qc.sh"
#   sbatch -c 16 --mem=8G  -t 60      --wrap="scripts/10_mapped_length_dist.sh"
#   sbatch -c 16 --mem=16G -t 90      --wrap="scripts/12_step5_biotype.sh"
#
# Layout: scripts/ code (numbered by run order), tables/<dataset>/ TSVs,
# figures/<step>/ PNG+PDF. Everything below is seconds on the login node.
set -euo pipefail
R=/nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/Rscript
cd "$(dirname "$(readlink -f "$0")")"
S=scripts
OUT75=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out75

step() { printf '\n=== [%s] %s ===\n' "$1" "$2"; }

tables() {
  step 3 "build_tables.R -- own library, 130 nt"
  "$R" $S/03_build_tables.R
  step 3 "build_tables.R -- own library, 75 nt"
  if [ -d "$OUT75/cells" ]; then
    [ -s tables/own75/demux_read_counts.tsv ] || cp tables/own130/demux_read_counts.tsv tables/own75/
    DEMO_RUN="$OUT75" DEMO_OUT="$PWD/tables/own75" "$R" $S/03_build_tables.R
  else
    echo "  skipped -- no $OUT75/cells"
  fi
  step 4 "build_tables_plate.R -- published VASA-plate"
  "$R" $S/04_build_tables_plate.R
}
plots() {
  step 6 "plot_all.R      -- reads / trim / classes / length / cross"
  "$R" $S/06_plot_all.R
  step 7 "plot_insilico.R -- step 3, in-silico rRNA depletion"
  "$R" $S/07_plot_insilico.R
  step 8 "plot_probe_qc.R -- probe-scoped residual"
  "$R" $S/08_plot_probe_qc.R
  step 9 "plot_step4.R    -- step 4, STAR mapping"
  "$R" $S/09_plot_step4.R
  step 11 "plot_step4_length.R -- step 4, read length by mapping outcome"
  "$R" $S/11_plot_step4_length.R
  step 13 "plot_step5.R    -- step 5, what the assigned reads are"
  "$R" $S/13_plot_step5.R
}
case "${1:-all}" in
  tables) tables ;;
  plots)  plots ;;
  all)    tables; plots ;;
  *) echo "usage: $0 [all|tables|plots]" >&2; exit 2 ;;
esac
printf '\ndone.\n'
