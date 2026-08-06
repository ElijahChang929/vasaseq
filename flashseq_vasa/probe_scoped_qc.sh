#!/bin/bash
# Probe-scoped rRNA depletion QC.
#
# Rationale: in-silico depletion removes ALL rRNA species (unchanged, deliberately).
# But judging the WET-LAB RNase H reaction against that total is unfair, because most
# of what survives (spacers, 47S precursor) was never a probe target. This script
# scores residual rRNA ONLY inside the probe-addressable footprint, so the number
# reflects the experiment rather than the biology of rRNA processing.
#
# Probe set: Adiconis 2013 Supplementary Table 6 (195 x 50-mer), as used by VASA-seq.
# Those probes are HUMAN-designed; probe_target_intervals.tsv holds the mouse
# footprint where a 50 nt window retains >=90% identity (i.e. can still hybridise).
set -euo pipefail

module load SAMtools/1.11-GCC-10.2.0
SAMTOOLS=$(which samtools)
echo "samtools: $SAMTOOLS"

CELLS=/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out/cells
IV=./probe_target_intervals.tsv
OUT=./probe_scoped_qc.tsv          # job workdir on SCRATCH; nothing written to working/
printf 'cell\tclass\tn\n' > "$OUT"

for bam in "$CELLS"/*.Ribo.bam; do
  cell=$(basename "$bam" _cbc_trimmed_homoATCG.Ribo.bam)
  "$SAMTOOLS" view "$bam" | awk -v cell="$cell" -v IV="$IV" 'BEGIN{
      OFS="\t"
      while ((getline line < IV) > 0) {
        n=split(line, f, "\t")
        if (f[1]=="contig") continue          # header
        k=++cnt[f[1]]
        S[f[1] SUBSEP k]=f[2]; E[f[1] SUBSEP k]=f[3]; T[f[1] SUBSEP k]=f[4]
      }
      close(IV)
    }
    {
      flag=$2; rname=$3; pos=$4; cig=$6
      if (and(flag,4) || and(flag,256) || and(flag,2048)) next
      span=0; num=""
      for(i=1;i<=length(cig);i++){ ch=substr(cig,i,1)
        if(ch ~ /[0-9]/){num=num ch}
        else{ if(ch=="M"||ch=="D"||ch=="N"||ch=="="||ch=="X") span+=num+0; num="" } }
      if(span<1) next
      rend=pos+span-1
      tot++
      # best-overlapping probe interval on this contig
      bestov=0; bestT=""
      m = (rname in cnt) ? cnt[rname] : 0
      for(k=1;k<=m;k++){
        s=S[rname,k]; e=E[rname,k]
        if (rend < s || pos > e) continue
        lo=(pos>s?pos:s); hi=(rend<e?rend:e); ov=hi-lo+1
        if (ov>bestov){ bestov=ov; bestT=T[rname,k] }
      }
      if (bestov>0) { inprobe++; byT[bestT]++
                      if (bestov==span) fully++ }
      else offprobe++
    }
    END{
      print cell, "all_rRNA_reads", tot
      print cell, "probe_target_reads", inprobe+0
      print cell, "probe_target_fully_contained", fully+0
      print cell, "off_probe_reads", offprobe+0
      for(t in byT) print cell, t, byT[t]
    }' >> "$OUT"
done

# per-cell library size, for denominators
printf 'cell\ttotal_trimmed\tnonRibo\n' > ./libsize.tsv
for lg in "$CELLS"/*.ribo-map.log; do
  cell=$(basename "$lg" _cbc_trimmed_homoATCG.ribo-map.log)
  tot=$(awk -F': ' '/^Number of reads/{print $2}' "$lg")
  unm=$(awk -F': ' '/^Number of unmapped reads/{print $2}' "$lg")
  printf '%s\t%s\t%s\n' "$cell" "$tot" "$unm" >> ./libsize.tsv
done

wc -l "$OUT" ./libsize.tsv
