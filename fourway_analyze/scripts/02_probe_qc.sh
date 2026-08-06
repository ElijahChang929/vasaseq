#!/bin/bash
###############################################################################
# 02_probe_qc.sh -- how much probe-targeted rRNA the wet-lab reaction left
# behind, on all four datasets.
#
# ONE NUMBER PER DATASET: probe-target residual reads / reads entering step 3.
#
# NOT "% rRNA". Most residual was never a probe target -- the 5'ETS / ITS1 /
# ITS2 / 3'ETS spacers carry no complementary probe, so a perfect reaction
# leaves every one of them. Scoring only the probe-bindable windows is what
# makes this a measure of the reaction rather than of the library.
#
# FLASH-SEQ HAD NO PROBES, AND THAT IS THE POINT
# ----------------------------------------------
# The RNase H probe depletion is part of the VASA protocol. FLASH-seq does not
# do it, so its bar is not "a worse reaction" -- it is the no-reaction control
# this figure never had. Read it as the baseline the VASA bars are working
# down from, and note the FLASH-seq libraries are polyA-primed, which suppresses
# rRNA by a different mechanism entirely and before any probe is involved.
#
# SCOPED TO THE 47S CONTIG, DELIBERATELY
# --------------------------------------
# `mouse_rDNA_47S_BK000964.3_1-13403` is the only contig every rRNA reference
# here names identically -- the three mouse-referenced datasets carry
# unique_rRNA_mouse.v2.fa, the published plate carries the mixed
# unique_rRNA_human_mouse.v3.fa, and the other 156 interval contigs (Rn18s-rs5,
# the 5S family, mito) have no counterpart name across both. Restricting to the
# 47S makes the comparison exact instead of approximate, and by the by-target
# breakdown it keeps ~98% of probe-target residual (28S 78.1% + 18S 16.8% +
# 5.8S 3.1%).
#
# Intervals are demo_analyze's, unchanged: Adiconis 2013 human 50-mers projected
# onto mouse at >=90% identity per 50 nt window.
#
# Output: tables/cross/probe_qc.tsv
#         tables/cross/probe_qc_per_unit.tsv
#
#   sbatch -c 8 --mem=8G -t 4:00:00 --wrap="scripts/02_probe_qc.sh"
###############################################################################
set -euo pipefail
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
cd "$ROOT"
source scripts/datasets.sh

IV=$VASA/code/I_Gene_expression/own_version/demo_analyze/probe_reference/probe_target_intervals.mouse.tsv
C47=mouse_rDNA_47S_BK000964.3_1-13403
NPAR=${NPAR:-8}
[ -s "$IV" ] || { echo "no interval file at $IV" >&2; exit 1; }

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/camp/apps/eb/modules/all
module load SAMtools/1.11-GCC-10.2.0 2>/dev/null
SAM=/camp/apps/eb/software/SAMtools/1.11-GCC-10.2.0/bin/samtools
[ -x "$SAM" ] || { echo "no samtools at $SAM" >&2; exit 1; }

TMP=$(mktemp -d "${TMPDIR:-/tmp}/pqc4.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# Probe-target reads on the 47S of one Ribo.bam. Whole-BAM scan with the contig
# filtered in awk: these BAMs are not indexed (nothing downstream needs them to
# be), so a region query errors out.
score() {
    ds=$1; unit=$2; bam=$3; log=$4; out=$5
    hit=$("$SAM" view "$bam" 2>/dev/null | awk -v IV="$IV" -v C="$C47" 'BEGIN{
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
      END{ print hit+0 }')
    n=$(awk '/Number of reads:/{gsub(/[^0-9]/,"",$NF); print $NF; exit}' "$log")
    printf '%s\t%s\t%d\t%d\n' "$ds" "$unit" "$hit" "$n" > "$out"
}
export -f score
export SAM IV C47

: > "$TMP/jobs"
i=0
for key in $DS_KEYS; do
    label=$(ds_label "$key")
    while IFS=$'\t' read -r unit map_stem ribo_stem; do
        bam="${ribo_stem}.Ribo.bam"; log="${ribo_stem}.ribo-map.log"
        [ -s "$bam" ] || { echo "missing $bam" >&2; exit 1; }
        [ -s "$log" ] || { echo "missing $log" >&2; exit 1; }
        i=$((i+1))
        printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$unit" "$bam" "$log" "$TMP/p.$i" >> "$TMP/jobs"
    done < <(ds_units "$key")
done

echo "scoring $(wc -l < "$TMP/jobs") Ribo.bam, $NPAR at a time..."
xargs -d '\n' -P "$NPAR" -I{} bash -c 'set -eo pipefail
    IFS=$'"'"'\t'"'"' read -r a b c d e <<< "{}"; score "$a" "$b" "$c" "$d" "$e"' < "$TMP/jobs"

{ printf 'dataset\tunit\tprobe_residual\treads_in\tpct\n'
  cat "$TMP"/p.* | awk -F'\t' '{printf "%s\t%s\t%d\t%d\t%.4f\n", $1,$2,$3,$4, ($4?100*$3/$4:0)}' \
    | sort -t$'\t' -k1,1 -k2,2
} > $OUTROOT/tables/cross/probe_qc_per_unit.tsv
echo "wrote $OUTROOT/tables/cross/probe_qc_per_unit.tsv"

# Pooled, not a mean of per-unit rates: units differ 20-60x in depth.
{ printf 'dataset\tunits\tprobe_residual\treads_in\tpct\n'
  awk -F'\t' 'NR>1 { n[$1]++; h[$1]+=$3; r[$1]+=$4; if(!($1 in s)){s[$1]=1; o[++k]=$1} }
    END { for(i=1;i<=k;i++){ d=o[i]
            printf "%s\t%d\t%d\t%d\t%.2f\n", d, n[d], h[d], r[d], 100*h[d]/r[d] } }' \
    $OUTROOT/tables/cross/probe_qc_per_unit.tsv
} > $OUTROOT/tables/cross/probe_qc.tsv
echo "wrote $OUTROOT/tables/cross/probe_qc.tsv"
echo
column -t -s$'\t' $OUTROOT/tables/cross/probe_qc.tsv
