#!/bin/bash
###############################################################################
# probe_qc.sh -- how much probe-targeted rRNA the RNase H reaction left behind.
#
# ONE NUMBER PER DATASET: probe-target residual reads / reads entering step 3.
#
# Not "% rRNA". Most of the residual was never a probe target -- the 5'ETS /
# ITS1 / ITS2 / 3'ETS spacers carry no complementary probe, so a perfect
# reaction leaves every one of them. Scoring only the probe-bindable windows is
# what makes this a measure of the wet-lab reaction rather than of the library.
# Intervals: flashseq_vasa/probe_target_intervals.tsv (Adiconis 2013 human
# 50-mers projected onto mouse at >=90% identity per 50 nt window).
#
# SCOPED TO THE 47S CONTIG, DELIBERATELY
# --------------------------------------
# The two runs used different rRNA references, and `mouse_rDNA_47S_BK000964.3
# _1-13403` is the only contig they name identically -- the other 156 interval
# contigs (Rn18s-rs5, the 5S family, mito) have no counterpart name in the
# plate's mixed reference. Restricting both sides to the 47S makes the
# comparison exact instead of approximate. By the by-target breakdown in
# probe_scoped_qc_README.md that keeps ~98% of probe-target residual
# (28S 78.1% + 18S 16.8% + 5.8S 3.1%); 5S and mito are the ~2% left out.
#
# The published side is MOUSE WELLS ONLY (Fig. 1d rule, typed by the deposited
# table) -- the probes are human-designed, and scoring human cells against
# mouse-projected intervals would answer a different question.
#
#   scripts/probe_qc.sh   ->  tables/cross/probe_qc.tsv
###############################################################################
set -euo pipefail
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
cd "$ROOT"

IV=$ROOT/probe_reference/probe_target_intervals.mouse.tsv
C47=mouse_rDNA_47S_BK000964.3_1-13403
OWN=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out
PLATE=/nemo/lab/turnerj/working/guangxin/vasaseq/data/ref/fastq_vasaplate/vasaplate_out_v3
PERCELL=/nemo/lab/turnerj/working/guangxin/vasaseq/res/vasaplate/per_cell.tsv

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/camp/apps/eb/modules/all
module load SAMtools/1.11-GCC-10.2.0 2>/dev/null
SAM=/camp/apps/eb/software/SAMtools/1.11-GCC-10.2.0/bin/samtools

# probe-target reads on the 47S of one BAM
score() {
    # whole-BAM scan, contig filtered in awk: these Ribo.bam are not indexed
    # (nothing downstream needs them to be), so a region query errors out.
    "$SAM" view "$1" 2>/dev/null | awk -v IV="$IV" -v C="$C47" 'BEGIN{
        while ((getline l < IV) > 0){ split(l,f,"\t")
          if (f[1]!=C) continue; k=++n; S[k]=f[2]; E[k]=f[3] } close(IV) }
      { if ($3 != C) next
        if (and($2,256)||and($2,2048)) next
        span=0; num=""
        for(i=1;i<=length($6);i++){ c=substr($6,i,1)
          if(c~/[0-9]/) num=num c
          else { if(c~/[MDN=X]/) span+=num+0; num="" } }
        if(span<1) next
        e=$4+span-1
        for(k=1;k<=n;k++) if(!(e<S[k] || $4>E[k])) { hit++; break } }
      END{ print hit+0 }'
}
# reads entering step 3, from that unit's own ribo-map log
reads_in() { awk '/Number of reads:/{gsub(/[^0-9]/,"",$NF); print $NF; exit}' "$1"; }

printf 'dataset\tunits\tprobe_residual\treads_in\tpct\n' > tables/cross/probe_qc.tsv

# --- own library ------------------------------------------------------------
h=0; r=0; u=0
for bam in "$OWN"/cells/*.Ribo.bam; do
    log=${bam%.Ribo.bam}.ribo-map.log
    [ -s "$log" ] || continue
    h=$((h + $(score "$bam"))); r=$((r + $(reads_in "$log"))); u=$((u+1))
done
printf 'own library\t%d\t%d\t%d\t%.2f\n' "$u" "$h" "$r" "$(echo "100*$h/$r" | bc -l)" >> tables/cross/probe_qc.tsv

# --- published plate, mouse wells only --------------------------------------
mouse=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="source")s=i;if($i=="well")w=i;if($i=="call_fig1d")c=i}}
                    NR>1 && $s=="published" && $c=="mouse"{printf "%03d\n",$w}' "$PERCELL")
h=0; r=0; u=0
for w in $mouse; do
    bam=$(ls "$PLATE"/*_"${w}"_cbc_trimmed_homoATCG.Ribo.bam 2>/dev/null | head -1)
    log=${bam%.Ribo.bam}.ribo-map.log
    [ -n "$bam" ] && [ -s "$log" ] || continue
    h=$((h + $(score "$bam"))); r=$((r + $(reads_in "$log"))); u=$((u+1))
done
printf 'published, mouse wells\t%d\t%d\t%d\t%.2f\n' "$u" "$h" "$r" "$(echo "100*$h/$r" | bc -l)" >> tables/cross/probe_qc.tsv

column -t -s$'\t' tables/cross/probe_qc.tsv
