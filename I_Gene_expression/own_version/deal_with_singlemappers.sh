#!/bin/bash
###############################################################################
# OWN_VERSION FORK -- differs from a_Mapping/deal_with_singlemappers.sh in
# exactly TWO places. Everything else is byte-identical to upstream.
#
# (1) READ SELECTION -- upstream matched the raw SAM text with
#         $0 ~ /NH:i:1\tHI:i:1\t/          (this script)
#         $0 ~ /NH:i:[2-9]/                (deal_with_multimappers.sh)
#     [2-9] is ONE character, so it only ever inspects the first digit after
#     "NH:i:". NH:i:10 .. NH:i:19 therefore match NEITHER pattern -- "NH:i:1"
#     is not followed by a TAB, and "1" is not in [2-9] -- and those alignments
#     were silently dropped by both scripts. (NH:i:20 matched only by accident,
#     "NH:i:2" being a prefix of it.) STAR here runs --outFilterMultimapNmax 20,
#     so 10-19 is a perfectly legal range: measured on cell 011, 273,092 of
#     3,000,000 alignments -- about 21,400 reads, ~5% of all multimapping reads.
#     Now the NH value is parsed as a NUMBER: nh==1 here, nh>=2 for multi.
#     Unmapped reads carry NH:i:0 (--outSAMunmapped Within) and are excluded by
#     both the old and the new test.
#
# (2) `if (readstrand="+")` -- a single `=` is an ASSIGNMENT in awk, not a
#     comparison. It is always true, so the `else` branch was dead code, and it
#     overwrote readstrand with "+". Consequences: jS:5/jS:3 came out backwards
#     on minus-strand genes; the (correctly written) `==` on the next line then
#     computed the Cov column from the wrong end of the gene; and the printed
#     Strand column became "+". Now `==`.
#     NB none of those three is read by countTables_2pickle_cellsSpliced.py --
#     it only ever tests jSs=='IN', and never touches Strand or Cov -- so this
#     one does NOT change the count tables. It is fixed because it is wrong,
#     and because anyone using Cov or Strand for QC would be misled.
#
# a_Mapping/ is left untouched: it is published code, and upstream's own
# single-pass workflow is not affected by either issue in a way that changes
# their results. See README.md, "Step 5".
###############################################################################


if [ $# -ne 5 ]
then
    echo "Please, give:"
    echo "1) input bam file"
    echo "2) bed file for introns, exons and tRNA"
    echo "3) stranded protocol (n/y)"
    echo "4) path to samtools"
    echo "5) path to bedtools"
    exit
fi

inbam=$1
refBED=$2
stranded=$3
p2samtools=$4
p2bedtools=$5
#refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed

${p2samtools}/samtools view -h $inbam | awk 'BEGIN{OFS="\t"} {
    if ($1 ~ /^@/) {print $0}
    else {
        nh=-1; nm=""
        for (i=1; i<=NF; i++) {
            if ($i ~ /^NH:i:[0-9]+$/) {nh = substr($i, 6) + 0}
            else if ($i ~ /nM:i:[0-9]/) {
                col=i; nm=substr($col, 6, length($col))
            }
        };
        if (nh == 1) {$1=$1";CG:"$6";nM:"nm; print $0}
    }
}' | ${p2samtools}/samtools view -Sb > ${inbam%.bam}.singlemappers.bam

${p2bedtools}/bedtools bamtobed -i ${inbam%.bam}.singlemappers.bam | ${p2bedtools}/bedtools sort > ${inbam%.bam}.singlemappers.bed

if [ $stranded == "y" ]
then
    ${p2bedtools}/bedtools intersect -a ${inbam%.bam}.singlemappers.bed -b $refBED -wa -wb | awk 'BEGIN {OFS="\t"; w="T"} {
        chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;
        sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        if (readstrand==refstrand) {
            if ((readstart >= refstart) && (readend <= refend)) {
                readname=readname";jS:IN"; w="T"
            } else if ((readstart < refstart) && (readend > refend)) {
                readname=readname";jS:OUT"; w="F";
            } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
                if (readstrand=="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
                } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
                if (readstrand=="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
            } else {print $0 > "checkme.txt"}
            if (readstrand=="+") {x=1-(geneend-readend)/genelen} else {x=1-(readstart-genestart)/genelen}
            sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
            if (w=="T") {print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
        }
    }' > ${inbam%.bam}.singlemappers_genes.bed

elif [ $stranded == 'n' ]
then
    ${p2bedtools}/bedtools intersect -a ${inbam%.bam}.singlemappers.bed -b $refBED -wa -wb | awk 'BEGIN {OFS="\t"; w="T"} {
        chr=$1; readstart=$2; readend=$3; readname=$4; readstrand=$6; refstart=$8; refend=$9; refstrand=$10; refname=$11; genelen=$12; genestart=$13; geneend=$14;
        sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        if ((readstart >= refstart) && (readend <= refend)) {
            readname=readname";jS:IN"; w="T"
        } else if ((readstart < refstart) && (readend > refend)) {
            readname=readname";jS:OUT"; w="F";
        } else if ( ((readstart < refstart)&&(readend <= refend)) || ((readstart <= refstart)&&(readend < refend)) ) {
            if (readstrand=="+") {readname=readname";jS:5"} else {readname=readname";jS:3"}; w="T"
            } else if ( ((readstart > refstart) && (readend >= refend)) || ((readstart >= refstart) && (readend > refend)) ) {
            if (readstrand=="+") {readname=readname";jS:3"} else {readname=readname";jS:5"}; w="T"
        } else {print $0 > "checkme.txt"}
        if (readstrand=="+") {x=1-(geneend-readend)/genelen} else {x=1-(readstart-genestart)/genelen}
        sx=match(readname, /;CG:/); rn=substr(readname, 0, sx-1); rq=substr(readname,sx+1,length(readname))
        if (w=="T") {print chr, readstart, readend, rn, readstrand, refname, rq, refend-refstart, x}
    }' > ${inbam%.bam}.singlemappers_genes.bed
fi

gzip ${inbam%.bam}.singlemappers_genes.bed

#sort -k4 ${inbam%.bam}.singlemappers_genes.bed > ${inbam%.bam}.nsorted.singlemappers_genes.bed # unnecessary

rm ${inbam%.bam}.singlemappers.bam ${inbam%.bam}.singlemappers.bed
