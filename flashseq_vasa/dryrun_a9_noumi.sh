#!/bin/bash
###############################################################################
# dryrun_a9_noumi.sh -- smallest genuine end-to-end slice that proves VASA's
# protocol='smartseq_noUMI' path runs on FLASH-seq data.
#
# Library ZHA8833A9 (30 pg rung, qc_verdict=ok). A deterministic stride
# subsample through steps 4-7 of the VASA pipeline, TWICE: once paired-end and
# once R1-only, both landing in ONE cells/ folder so step 6 makes them two
# columns of one table and step 7 quantifies them side by side. Nothing in
# a_Mapping/ is modified; every upstream script is called as-is.
#
# WHAT IS DELIBERATE HERE, AND WHY
# --------------------------------
# 1. STRIDE, NOT HEAD. The head of a fastq is one end of the flowcell and is
#    not a sample of the library (measured on ZHA8833A1: adapter rate 58.5% in
#    the first 20k reads vs 55.1% over the whole library). Every Nth read pair
#    is taken instead, R1/R2 in lockstep so the pairing survives -- and that is
#    checked, not assumed.
#
# 2. STRANDED=n. FLASH-seq is unstranded; the forward strand carries 49.1-50.5%
#    of ribosomal reads (measured, code/flashseq/README.md). Passing y would
#    halve every biotype figure. The VASA side runs with y.
#
# 3. PE AND SE, BOTH. In PE a fragment is two alignment records sharing one
#    QNAME; both mates carry the same NH so both survive the singlemapper
#    filter, and step 6 groups the BED by the Name column (= QNAME), so the two
#    mates should collapse into ONE assignment. That is an assumption about
#    upstream code written for single-end reads, so the SE arm is run as the
#    control and the two columns are compared directly.
#
# 4. NO SHARED MEMORY. Upstream map_star.sh passes no --genomeLoad, so STAR
#    defaults to NoSharedMemory; that is kept. It costs a second index load but
#    cannot leave 23 GB stranded in a node's shm if the job dies.
#
# 5. CELLID_FROM=f, and the '_cbc' in the filename. Step 6 mode 'f' does
#    cellfile[:cellfile.index('_cbc')], which raises ValueError without a
#    '_cbc'. FLASH-seq has no cell barcode, so the token is put in the filename
#    purely to satisfy that contract. Mode 'r' would need an SM tag injected
#    into every read name. See NOUMI_PATH.md.
#
# 6. filt_unigenes = n. At ncols=2 the step-7 filter is max(5, round(0.01*2))
#    = 5 columns out of 2, which no gene can pass.
#
# Usage: dryrun_a9_noumi.sh [outdir] [stride]
###############################################################################
set -euo pipefail

OUTDIR="${1:-$PWD/dryrun_a9}"
STRIDE="${2:-64}"

W=/nemo/lab/turnerj/working/guangxin/vasaseq
VASA_SCRIPTS=$W/code/I_Gene_expression/a_Mapping          # upstream, untouched
VASA_OWN=$W/code/I_Gene_expression/own_version            # the NH>=10 fork
FQDIR=/nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/RN26038/20260325_LH00442_0237_B23GT7GLT3/fastq
LIB=ZHA8833A9
SAMPLE=A9dry
R1=$FQDIR/${LIB}_S108_L007_R1_001.fastq.gz
R2=$FQDIR/${LIB}_S108_L007_R2_001.fastq.gz

STAR_INDEX=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116
REF_BED=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed
STRANDED=n

EBROOT=/camp/apps/eb/software
P2STAR=$EBROOT/STAR/2.7.7a-GCC-10.2.0/bin
P2SAMTOOLS=$EBROOT/SAMtools/1.11-GCC-10.2.0/bin
P2BEDTOOLS=$EBROOT/BEDTools/2.30.0-GCC-11.2.0/bin
CONDA_ENV=/nemo/lab/turnerj/working/guangxin/envs/vasa

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=$EBROOT/../modules/all
module load STAR/2.7.7a-GCC-10.2.0 SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0
source $EBROOT/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate $CONDA_ENV

say() { echo "[$(date +%H:%M:%S)] $*"; }

CELLDIR=$OUTDIR/cells
mkdir -p "$CELLDIR" "$OUTDIR/logs"
export TMPDIR="$OUTDIR/tmp"; mkdir -p "$TMPDIR"
cd "$OUTDIR"

SUB1=$CELLDIR/sub_R1.fastq.gz
SUB2=$CELLDIR/sub_R2.fastq.gz

###############################################################################
say "step A: stride-$STRIDE subsample (every ${STRIDE}th read PAIR, whole file)"
###############################################################################
if [ ! -s "$SUB1" ] || [ ! -s "$SUB2" ]; then
    for spec in "R1|$R1|$SUB1" "R2|$R2|$SUB2"; do
        IFS='|' read -r tag src dst <<< "$spec"
        say "  streaming $tag"
        zcat "$src" | awk -v s="$STRIDE" 'NR%4==1{keep=(int((NR-1)/4)%s==0)} keep' \
            | gzip -c > "$dst"
    done
fi
n1=$(( $(zcat "$SUB1" | wc -l) / 4 ))
n2=$(( $(zcat "$SUB2" | wc -l) / 4 ))
say "  subsampled R1=$n1 R2=$n2 read records"
[ "$n1" -eq "$n2" ] || { echo "FATAL: R1/R2 out of lockstep"; exit 1; }

# Prove the pairing survived: mate names must match record for record.
mism=$(paste <(zcat "$SUB1" | awk 'NR%4==1{print $1}') \
             <(zcat "$SUB2" | awk 'NR%4==1{print $1}') \
       | awk '$1!=$2{c++} END{print c+0}')
say "  mate-name mismatches: $mism"
[ "$mism" -eq 0 ] || { echo "FATAL: subsample broke pairing"; exit 1; }

###############################################################################
# steps B-C, once per arm. The '_cbc' token is load-bearing (note 5).
###############################################################################
run_arm() {
    local arm=$1; shift
    local stem="${LIB}${arm}_cbc_noumi_E99_"
    local bam="$CELLDIR/${stem}Aligned.out.bam"

    say "step B[$arm]: STAR, upstream map_star.sh flag set"
    if [ ! -s "$bam" ]; then
        "$P2STAR"/STAR --runThreadN 8 --genomeDir "$STAR_INDEX" \
            --readFilesIn "$@" --readFilesCommand zcat \
            --outFilterMultimapNmax 20 --outSAMunmapped Within \
            --outSAMtype BAM Unsorted --outSAMattributes All \
            --outFileNamePrefix "$CELLDIR/$stem" > /dev/null
        rm -rf "$CELLDIR/${stem}_STARtmp" "$CELLDIR/${stem}Log.progress.out"
        mv -f "$CELLDIR/${stem}Log.final.out" "$CELLDIR/${stem}Log.final.txt"
        mv -f "$CELLDIR/${stem}Log.out" "$CELLDIR/${stem}Log.txt"
    fi
    grep -E 'input reads|Uniquely mapped reads|multiple loci|too short' \
        "$CELLDIR/${stem}Log.final.txt" | sed 's/^/      /'

    say "step C[$arm]: assignment, own_version singlemappers fork + upstream multi, stranded=$STRANDED"
    # own_version/deal_with_singlemappers.sh parses NH as a NUMBER, so
    # NH:i:10..19 are not silently dropped, and fixes the `=`-for-`==` awk bug.
    # The multi side is upstream's, exactly as pipeline.sh calls it.
    local sgl="$CELLDIR/${stem}Aligned.out.singlemappers_genes.bed.gz"
    local mlt="$CELLDIR/${stem}Aligned.out.nsorted.multimappers_genes.bed.gz"
    [ -s "$sgl" ] || "$VASA_OWN"/deal_with_singlemappers.sh "$bam" "$REF_BED" "$STRANDED" "$P2SAMTOOLS" "$P2BEDTOOLS"
    [ -s "$mlt" ] || "$VASA_SCRIPTS"/deal_with_multimappers.sh "$bam" "$REF_BED" "$STRANDED" "$P2SAMTOOLS" "$P2BEDTOOLS"

    local nsgl nmlt nq m1 m2
    nsgl=$(zcat "$sgl" | wc -l); nmlt=$(zcat "$mlt" | wc -l)
    nq=$(zcat "$sgl" | cut -f4 | sort -u | wc -l)
    m1=$("$P2SAMTOOLS" view -c -f 64 "$bam"); m2=$("$P2SAMTOOLS" view -c -f 128 "$bam")
    say "      single BED rows=$nsgl  multi BED rows=$nmlt  distinct QNAMEs(single)=$nq"
    say "      BAM records flagged mate1=$m1 mate2=$m2"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$arm" "$nsgl" "$nmlt" "$nq" "$m1" "$m2" \
        >> "$OUTDIR/logs/arm_stats.tsv"
}

printf 'arm\tsingle_bed_rows\tmulti_bed_rows\tdistinct_qnames_single\tbam_mate1\tbam_mate2\n' \
    > "$OUTDIR/logs/arm_stats.tsv"
run_arm pe "$SUB1" "$SUB2"
run_arm se "$SUB1"

###############################################################################
say "step D: step 6 -- pickle, protocol=smartseq_noUMI, cellid from filename"
###############################################################################
# Both arms are in cells/, so step 6 makes them two columns of one structure.
rm -f "$OUTDIR/${SAMPLE}.pickle" "$OUTDIR/${SAMPLE}.pickle.gz" "$OUTDIR/${SAMPLE}dict.pickle"
( cd "$OUTDIR" && "$VASA_SCRIPTS"/countTables_2pickle_cellsSpliced.py \
    "cells" "$SAMPLE" smartseq_noUMI f )
ls -la "$OUTDIR/${SAMPLE}.pickle.gz" "$OUTDIR/${SAMPLE}dict.pickle"

###############################################################################
say "step E: step 7 -- count tables, protocol=smartseq_noUMI, filt_unigenes=n"
###############################################################################
( cd "$OUTDIR" && "$VASA_SCRIPTS"/countTables_fromPickle.py \
    "${SAMPLE}.pickle.gz" "$SAMPLE" smartseq_noUMI n )
say "  mapStats.log ($(wc -l < "$OUTDIR/${SAMPLE}_mapStats.log") lines; a complete VASA log is 21):"
sed 's/^/      /' "$OUTDIR/${SAMPLE}_mapStats.log"
say "  tables written: $(find "$OUTDIR" -maxdepth 1 -name "${SAMPLE}*.tsv" | wc -l)"

###############################################################################
say "step F: sanity report over the tables step 7 just wrote"
###############################################################################
python "$W/code/flashseq_vasa/noumi_dryrun_report.py" \
    "$OUTDIR" "$SAMPLE" "$CELLDIR" "$OUTDIR/noumi_dryrun_A9.txt"

say "DRY RUN COMPLETE (stride $STRIDE, arms pe+se)"
