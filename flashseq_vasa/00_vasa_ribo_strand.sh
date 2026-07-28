#!/bin/bash
# vasa_ribo_strand.sh -- make the VASA side of the stranded=y/n question a
# MEASUREMENT, symmetric to what 05_rrna_bwa.sh already did for FLASH-seq.
#
# 05_rrna_bwa.sh ran FLASH-seq with stranded=n and then re-selected the SAME
# merged BAM with stranded=y, so its factor-of-two claim is measured. The VASA
# run only ever ran stranded=y, so the mirror-image number (what VASA's rRNA %
# would be under stranded=n) has never been computed. Both all-ribo BAMs
# survive for all 16 barcodes, so it is measurable rather than arguable.
#
# We do NOT re-run riboread-selection.py (that would rewrite step-3 outputs the
# downstream tables were built from). We count, over the same name-sorted BAM
# it read, exactly the two predicates its two branches use:
#
#   stranded == 'n' :  mapreads = [x for x in reads if not x.is_unmapped]
#   stranded == 'y' :  mapreads = [x for x in reads if (not x.is_unmapped)
#                                                   and (not x.is_reverse)]
#   a read is ribosomal iff len(mapreads) >= 1
#
# and we replicate its group bookkeeping exactly, including the fact that its
# loop counts the FINAL qname group in nreads but never classifies it (the
# `else` branch only fires on a qname CHANGE). So closed = groups - 1, and
# any_fwd over closed groups must equal sum(nmapped) in the existing
# .ribo-map.log. That equality is the script's own self-check: if it fails,
# this counting is wrong and nothing downstream should be believed.
#
# Portable arithmetic (int(f/4)%2) rather than gawk's and(), so it runs under
# whatever awk the node has.
set -uo pipefail

W=/nemo/lab/turnerj/working/guangxin/vasaseq
CELLS="$W/data/PM26037/out/cells"
SAMTOOLS=/camp/apps/eb/software/SAMtools/1.11-GCC-10.2.0/bin/samtools

[ -x "$SAMTOOLS" ] || { echo "FATAL: no samtools at $SAMTOOLS" >&2; exit 1; }

# find, not ls: `ls <no match>` exits 2 and would kill the job
mapfile -t BAMS < <(find "$CELLS" -name '*.nsorted.all-ribo.bam' | sort)
echo "n_bams=${#BAMS[@]}"
[ "${#BAMS[@]}" -eq 16 ] || echo "WARNING: expected 16 all-ribo BAMs, found ${#BAMS[@]}" >&2

count_one () {
    local bam="$1" tag="$2"
    "$SAMTOOLS" view "$bam" | awk -v cell="$tag" '
    BEGIN{ FS="\t"; init=0; groups=0; closed=0; anym=0; anyf=0; nrec=0 }
    {
        nrec++
        q=$1; f=$2+0
        if (!init) { cur=q; init=1; groups=1; m=0; fw=0 }
        else if (q != cur) {
            closed++
            if (m)  anym++
            if (fw) anyf++
            cur=q; groups++; m=0; fw=0
        }
        if (int(f/4)%2 == 0) {          # not unmapped
            m=1
            if (int(f/16)%2 == 0) fw=1  # and not reverse
        }
    }
    END{ printf "%s\t%d\t%d\t%d\t%d\t%d\n", cell, nrec, groups, closed, anym, anyf }
    '
}

printf 'cell\tn_records\tn_groups\tn_closed\tribo_stranded_n\tribo_stranded_y\n' > ribo_strand_counts.tsv
pids=(); tmps=()
for bam in "${BAMS[@]}"; do
    b=$(basename "$bam")
    cell=$(echo "$b" | sed -E 's/^ZHA9292A1_([0-9]+)_.*/\1/')
    t=$(mktemp ./rs.XXXXXX); tmps+=("$t")
    count_one "$bam" "$cell" > "$t" &
    pids+=($!)
done
# `wait` alone returns 0 whatever the children did; wait on each pid
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
[ "$rc" -eq 0 ] || { echo "FATAL: at least one counting child failed" >&2; exit 1; }
cat "${tmps[@]}" | sort >> ribo_strand_counts.tsv
rm -f "${tmps[@]}"

# ---- self-check against the logs riboread-selection.py itself wrote --------
# sum of the per-RG mapped counters == our any_fwd over closed groups
printf 'cell\tlog_nreads\tlog_nunmapped\tlog_nmapped_sum\n' > ribo_log_counts.tsv
for lg in $(find "$CELLS" -name '*.ribo-map.log' | sort); do
    cell=$(basename "$lg" | sed -E 's/^ZHA9292A1_([0-9]+)_.*/\1/')
    awk -v cell="$cell" '
        /^Number of reads:/            { nr=$NF }
        /^Number of unmapped reads:/   { nu=$NF }
        /^\t/                          { s+=$NF }
        END{ printf "%s\t%d\t%d\t%d\n", cell, nr, nu, s }
    ' "$lg"
done | sort >> ribo_log_counts.tsv

echo "--- ribo_strand_counts.tsv ---"; cat ribo_strand_counts.tsv
echo "--- ribo_log_counts.tsv ---";    cat ribo_log_counts.tsv

# ---- annotation: how much does nf-core's GTF filter remove? ----------------
GTF_FULL=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/annotation/release-116/gtf/Mus_musculus.GRCm39.116.gtf
GTF_FILT="$W/data/flashseq/results/genome/Mus_musculus.GRCm39.116.filtered.gtf"
{
  printf 'gtf\tn_lines\tn_gene_rows\tn_transcript_rows\tn_exon_rows\n'
  for g in "$GTF_FULL" "$GTF_FILT"; do
      awk -v f="$g" 'BEGIN{FS="\t"} !/^#/{n++; if($3=="gene")gg++; else if($3=="transcript")tt++; else if($3=="exon")ee++}
          END{printf "%s\t%d\t%d\t%d\t%d\n", f, n, gg, tt, ee}' "$g"
  done
} > gtf_compare.tsv
echo "--- gtf_compare.tsv ---"; cat gtf_compare.tsv

# ---- which genomeDir did the nf-core STAR align actually use? --------------
{
  echo "## nf-core STAR align genomeDir / versionGenome / sjdbOverhang"
  for f in $(find "$W/data/flashseq/results/star_rsem/log" -name '*.Log.out' | sort | head -2); do
      echo "# $f"
      grep -E '^\s*(genomeDir|versionGenome|sjdbOverhang|sjdbGTFfile)' "$f" | head -8
      grep -m1 'STAR version' "$f"
  done
} > nfcore_star_provenance.txt
cat nfcore_star_provenance.txt

echo "DONE_PART_A"
