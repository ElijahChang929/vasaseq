#!/bin/bash
###############################################################################
# 04_diagnose_residual_rrna.sh -- why our count tables carry ~600x the published
# share of rRNA-biotype UFIs, all of it on one locus.
#
# THE OBSERVATION
#   Every biotype agrees with GSM5369495 within a factor of two except rRNA:
#       ENSMUSG00000106106_CT010467.1_rRNA    ours 95,823 UFIs   published 72
#   Every other rRNA-biotype gene agrees. So it is one locus, not a global
#   depletion difference.
#
# THE CAUSE -- two defects that only bite in combination
#
#   A. riboread-selection.py writes reverse-strand reads REVERSE-COMPLEMENTED.
#      It emits `r0.seq` where `r0 = reads[0]` is the first BAM record of the
#      group. pysam's .seq returns the sequence AS STORED IN THE BAM, which for
#      a reverse alignment is the reverse complement of the original read. The
#      same applies to .qual, which is also reversed. Reads kept by the
#      stranded reprieve (all hits reverse) therefore enter STAR flipped, and
#      for a stranded protocol that inverts sense and antisense.
#
#      Measured on cell 005: 2,222 of 45,880 reads with >=1 rRNA hit (4.8%)
#      are spared this way, and every one is written flipped.
#
#   B. Our rRNA FASTA contains an Ensembl entry stored ANTISENSE to the real
#      transcript. Component A is a gene_biotype dump extracted with
#      `bedtools getfasta -s`, so each entry is the sense sequence OF THE
#      ENSEMBL GENE. ENSMUSG00000106106 is annotated on the strand opposite the
#      rRNA transcript: aligned to the 47S unit its entry comes back flag=16
#      (reverse) at position 4007, i.e. it is the 18S region, stored backwards.
#
#      So a genuine 18S read aligns REVERSE to that entry. If it does not also
#      hit the 47S unit forward, all its hits are reverse -> spared by
#      stranded=y -> written reverse-complemented by defect A -> STAR now maps
#      it SENSE to the minus-strand gene -> counted as rRNA.
#
#      This entry is UNIQUE in that respect: of the 915 Ensembl-derived entries,
#      897 have no homology to the NCBI units (5S and dispersed fragments), 17
#      are stored sense, and exactly ONE is stored antisense -- this one. See
#      EVIDENCE 5. That is why a single locus blows up and no other does.
#
# WHAT WE CANNOT SAY, AND EARLIER DRAFTS OF THIS FILE SAID ANYWAY
#   Do not claim the authors "used only NCBI sequences and no Ensembl dump".
#   Their rRNA FASTA was never deposited, the Methods name no accessions, and
#   a_Mapping/README.md's list is prefixed "e.g." -- it is an example, not a
#   specification. The published table also argues against that story: on every
#   rRNA-biotype gene EXCEPT this one their residual is equal to or higher than
#   ours (8 vs 2, 13 and 11 on human 5.8S where we have 0), i.e. our depletion of
#   the dispersed rRNA genes is at least as thorough as theirs. What the evidence
#   does support is narrower and sufficient: whatever their reference contained,
#   it did not carry this locus in an orientation that produces this artefact.
#
# WHAT THIS SCRIPT DOES
#   Reproduces all five pieces of evidence from files already on disk. Read
#   only; it writes nothing outside $TMP.
#
# Usage:  ./04_diagnose_residual_rrna.sh [cell]      (default 005)
###############################################################################
set -euo pipefail

CELL=${1:-005}
FQDIR=${FQDIR:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/ref/fastq_vasaplate}
RUN=${RUN:-vasaplate_out_rrnav2}
REFDIR=${REFDIR:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mixed}
FASTA=${FASTA:-${REFDIR}/unique_rRNA_human_mouse.v2.fa}
GENE=${GENE:-ENSMUSG00000106106}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/camp/apps/eb/modules/all
module purge >/dev/null 2>&1 || true
module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0

P=${FQDIR}/${RUN}/SRR14783059_${CELL}_cbc_trimmed_homoATCG
BAM=${P}.nsorted.all-ribo.bam
for f in "$BAM" "${P}.fq.gz" "${P}.nonRibo.fastq.gz" "$FASTA"; do
    [ -s "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
done

echo "=== EVIDENCE 1: the reference entry for ${GENE} is stored ANTISENSE ==="
awk -v g="$GENE" '$0 ~ "^>mouse_"g {k=1;print;next} /^>/{k=0} k' "$FASTA" > "$TMP/gene.fa"
awk '/^>mouse_rDNA_47S/{k=1;print;next} /^>/{k=0} k' "$FASTA" > "$TMP/47s.fa"
bwa index "$TMP/47s.fa" 2>/dev/null
bwa mem "$TMP/47s.fa" "$TMP/gene.fa" 2>/dev/null | samtools view - \
  | awk '{printf "  flag=%s %s   pos=%s (47S: 5ETS 1-4007, 18S 4008-5877)\n",
          $2, ($2==16?"= REVERSE, i.e. antisense to the transcript":"= forward"), $4}'

echo
echo "=== EVIDENCE 2: reads on this locus hit the rRNA reference REVERSE-only ==="
zcat "${P}.nonRibo_E99_Aligned.out.singlemappers_genes.bed.gz" \
  | grep "$GENE" | cut -f4 | sort -u > "$TMP/names.txt"
echo "  reads assigned to ${GENE} in cell ${CELL}: $(wc -l < "$TMP/names.txt")"
samtools view "$BAM" | awk 'NR==FNR{w[$0]=1;next} ($1 in w){
      rg="?"; for(i=12;i<=NF;i++) if($i ~ /^RG:Z:/) rg=($i ~ /aln-ribo/)?"aln":"mem"
      st=($2==4?"UNMAPPED":($2==16?"REVERSE":"FORWARD"))
      print "  "rg"\t"st"\t"$3}' "$TMP/names.txt" - | sort | uniq -c | sort -rn

echo
echo "=== EVIDENCE 3: such a read is written REVERSE-COMPLEMENTED into nonRibo ==="
# NB: no early `exit` in these awks. Closing the pipe early sends SIGPIPE to
# zcat/samtools, the pipeline returns 141, and `set -e` kills the script -- which
# is exactly what happened the first time this was written.
n=$(samtools view "$BAM" | awk 'NR==FNR{w[$0]=1;next} ($1 in w) && $2==16 && !seen++ {print $1}' \
      "$TMP/names.txt" -)
echo "  read: ${n}"
grab() { zcat "$1" | awk -v k="@$2" '{h=$0;getline s;getline p;getline q; if(h==k && !g++) print s}'; }
orig=$(grab "${P}.fq.gz" "$n")
kept=$(grab "${P}.nonRibo.fastq.gz" "$n")
rc=$(printf '%s' "$orig" | rev | tr 'ACGTacgt' 'TGCAtgca')
echo "  trimmed fastq (truth) : ${orig}"
echo "  nonRibo fastq (step 3): ${kept}"
if [ "$kept" = "$rc" ]; then
    echo "  >>> CONFIRMED: step 3 wrote the reverse complement of the original read"
elif [ "$kept" = "$orig" ]; then
    echo "  >>> unchanged (pick a read whose records are reverse, not unmapped)"
else
    echo "  >>> neither -- investigate"
fi

echo
echo "=== EVIDENCE 4: how many reads this affects in cell ${CELL} ==="
samtools view "$BAM" | awk '
    {m=($2!=4); rev=($2==16)
     if($1!=p){ if(p!=""){ if(nm>0) withhit++; if(nm>0 && nf==0) spared++ } p=$1; nm=0; nf=0 }
     if(m){ nm++; if(!rev) nf++ } }
    END{ if(nm>0) withhit++; if(nm>0 && nf==0) spared++
         printf "  reads with >=1 rRNA hit          : %d\n", withhit
         printf "  all hits reverse -> spared+flipped: %d (%.1f%%)\n", spared, 100*spared/withhit }'

echo
echo "=== EVIDENCE 5: how many Ensembl entries are stored antisense? ==="
# The census that makes this a one-locus defect rather than a systemic one.
awk '/^>(mouse_rDNA_47S|human_45S)/{k=1;print;next} /^>/{k=0} k' "$FASTA" > "$TMP/units.fa"
awk '/^>(mouse_ENS|human_ENS)/{k=1;print;next} /^>/{k=0} k'        "$FASTA" > "$TMP/ens.fa"
bwa index "$TMP/units.fa" 2>/dev/null
echo "  Ensembl-derived entries: $(grep -c '^>' "$TMP/ens.fa")"
bwa mem -t 4 "$TMP/units.fa" "$TMP/ens.fa" 2>/dev/null | samtools view - \
  | awk '$2==0{s++} $2==16{a++} $2==4{n++}
         END{printf "  sense (correct)          : %d\n  ANTISENSE (the defect)   : %d\n  no homology to the units : %d\n", s+0, a+0, n+0}'
bwa mem -t 4 "$TMP/units.fa" "$TMP/ens.fa" 2>/dev/null | samtools view - \
  | awk '$2==16{print "  antisense entry: "$1" -> "$3" @ "$4}'

echo
echo "FIXES, in order of cost:"
echo "  B (cheap, one-line): the defect is ONE entry of 915. Drop or re-orient"
echo "    the antisense entry in build_rrna_reference_mixed.sh (or orient every"
echo "    Ensembl entry to the NCBI units instead of trusting the annotation"
echo "    strand). Rebuild + re-run steps 3-7."
echo "  A (correct, but changes every run): in riboread-selection.py emit the"
echo "    original orientation for reverse records -- revcomp r0.seq and reverse"
echo "    r0.qual before writing. This is an upstream published script, so under"
echo "    this repo's conventions it is a genuine-bug fix made in place and noted."
