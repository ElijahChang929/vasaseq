#!/bin/bash
################################################################################
# submit_vasaplate_map_array.sh
#
# PURPOSE
#   Submits STAGES 2-7 of the VASA-plate mapping pipeline (trim -> ribo ->
#   gmap -> b2bs/b2bm -> cout -> pick) as SLURM JOB ARRAYS, for the case where
#   STAGE 1 has already been run and the per-cell *_cbc.fastq.gz files exist.
#   It is the fqext=n half of submit_vasaplate_map.sh, restructured; STAGE 1
#   is not submitted here at all.
#
# WHY THIS EXISTS ALONGSIDE submit_vasaplate_map.sh
#   VASA-plate demultiplexes to one fastq per cell, so the upstream driver's
#   per-lane loop actually iterates 384 times and issues 384 x 5 = 1920
#   independent sbatch jobs with afterany dependencies. Two problems with that
#   on this cluster:
#
#     1. It submits gmap with --mem=40G, but the MIXED STAR index is 53 GB.
#        Every one of the 384 STAR jobs would OOM. (This is a real bug in the
#        upstream resource request, not a tuning choice.)
#     2. 384 separate 53 GB index loads is ~20 TB of Lustre reads to do a few
#        hours of actual alignment.
#
#   This script fixes both: gmap runs batched via gmap_chunk.sh (one genome
#   load per chunk of cells, --mem=80G), and every stage is a job array with a
#   %N concurrency cap so the run is polite to the shared queue and can be
#   stopped or resumed with a single scancel.
#
#   submit_vasaplate_map.sh is left untouched and remains the reference for
#   the upstream stage semantics and for STAGE 1.
#
# OUTPUT EQUIVALENCE
#   Every stage runs the SAME helper script with the SAME arguments as
#   upstream, so all intermediate and final filenames are unchanged. STAGE 6
#   (countTables_2pickle_cellsSpliced.py) globs the folder for
#   *.singlemappers_genes.bed.gz and takes the cell id from the filename
#   prefix before '_cbc', so nothing downstream can tell the difference.
#
# USAGE
#   Run from the directory that holds the raw fastqs / the output folder,
#   exactly as with submit_vasaplate_map.sh:
#
#     cd data/ref/fastq_vasaplate
#     ${p2s}/submit_vasaplate_map_array.sh SRR14783059 MIXED 74 \
#            vasaplate_HEK293T-mESC vasaplate_out f
#
#   Set DRYRUN=1 to print the sbatch commands without submitting.
################################################################################

### input paths -- kept in sync with submit_vasaplate_map.sh, which remains the
### reference for why each of these is what it is (module-vs-conda split, the
### per-stage module loads, the Trim_Galore/foss-2018b libstdc++ conflict).
p2s=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/a_Mapping
email=zhangg@crick.ac.uk

### array knobs
cells_per_chunk=${cells_per_chunk:-24}   # cells mapped per gmap job (1 genome load each)
maxpar=${maxpar:-64}                     # %N concurrency cap on the per-cell arrays
maxpar_gmap=${maxpar_gmap:-8}            # %N cap on the gmap array (each holds ~56G of shm)

### START -- first stage to submit (2=trim, 3=ribo, 4=gmap, 5=b2bs/b2bm).
### Stages below it are assumed already done and present in $folder; the first
### stage submitted simply carries no dependency. Use this to re-run the tail of
### the pipeline after changing something that only affects that tail -- e.g.
### START=3 with a corrected rRNA reference, which does not touch trimming.
START=${START:-2}

### RIBOREF -- override the rRNA reference chosen by the MOUSE/HUMAN/MIXED
### branch below. The branch defaults are the v1 Ensembl-only references; see
### own_version/build_rrna_reference*.sh for why those are wrong and what
### replaced them.
RIBOREF=${RIBOREF:-}

### REFBED -- override the annotation BED chosen by the branch below. Same
### situation as RIBOREF: the MIXED default has 0 tRNA rows despite being named
### "...IntronExonTrna.bed", so every *_tRNA.*Counts.tsv step 7 writes is empty.
### own_version/build_annotation_bed_mixed.sh builds the replacement. Only
### stages 5-7 read this, so pair it with START=5.
REFBED=${REFBED:-}

EBROOT=/camp/apps/eb/software
p2trimgalore=${EBROOT}/Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
p2cutadapt=${EBROOT}/cutadapt/1.18-foss-2018b-Python-3.6.6/bin
p2bwa=${EBROOT}/BWA/0.7.17-GCC-10.3.0/bin
p2samtools=${EBROOT}/SAMtools/1.11-GCC-10.2.0/bin
p2star=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin
p2bedtools=${EBROOT}/BEDTools/2.30.0-GCC-11.2.0/bin

### per-stage module loads (never all at once -- see submit_vasaplate_map.sh)
ml_init="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=${EBROOT%/software}/modules/all"
ml_trim="${ml_init}; module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6"
ml_gmap="${ml_init}; module load STAR/2.7.7a-GCC-10.2.0 SAMtools/1.11-GCC-10.2.0"
ml_b2b="${ml_init}; module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0"
ca="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa"
ml_ribo="${ml_init}; module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0; ${ca}"

VASA_REFS=${VASA_REFS:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq}

### check input parameters
if [ $# -ne 6 ]
then
    echo "Please, give:"
    echo "1) library name (prefix of the per-cell *_cbc.fastq.gz files)"
    echo "2) genome: MOUSE / HUMAN / MIXED"
    echo "3) read length (must match the STAR index suffix)"
    echo "4) prefix for output files (unused, kept for argument parity with submit_vasaplate_map.sh)"
    echo "5) folder holding the *_cbc.fastq.gz files (also the output folder)"
    echo "6) cellID from filename or readname (f/r)"
    exit 1
fi

lib=$1
ref=$2
n=$3
out=$4
folder=$5
cellidori=$6

### set references (same three branches as submit_vasaplate_map.sh)
if [[ $ref == "MOUSE" ]]
then
    riboref=${VASA_REFS}/mouse/rRNA_mouse_Rn45S.fa
    genome=${VASA_REFS}/mouse/star_index_$n
    refBED=${VASA_REFS}/mouse/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
elif [[ $ref == "HUMAN" ]]
then
    riboref=${VASA_REFS}/human/unique_rRNA_human.fa
    genome=${VASA_REFS}/human/star_index_$n
    refBED=${VASA_REFS}/human/Homosapines_ensemble99.homemade_IntronExonTrna.bed
elif [[ $ref == "MIXED" ]]
then
    riboref=${VASA_REFS}/mixed/unique_rRNA_human_mouse.fa
    genome=${VASA_REFS}/mixed/star_index_$n
    refBED=${VASA_REFS}/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.bed
else
    echo "unknown reference '$ref' (expected MOUSE / HUMAN / MIXED)"
    exit 1
fi
if [ -n "$RIBOREF" ]; then riboref=$RIBOREF; fi
if [ -n "$REFBED"  ]; then refBED=$REFBED;   fi

### preflight
missing=0
for tool in ${p2star}/STAR ${p2bwa}/bwa ${p2samtools}/samtools ${p2bedtools}/bedtools ${p2trimgalore}/trim_galore ${p2cutadapt}/cutadapt
do
    if [ ! -x "$tool" ]; then echo "MISSING tool:      $tool"; missing=1; fi
done
for f in "$riboref" "$refBED"
do
    if [ ! -s "$f" ]; then echo "MISSING reference: $f"; missing=1; fi
done
if [ ! -d "$genome" ]; then echo "MISSING STAR index: $genome"; missing=1; fi
for s in trim.sh ribo-bwamem.sh gmap_chunk.sh deal_with_singlemappers.sh deal_with_multimappers.sh countTables_2pickle_cellsSpliced.py countTables_fromPickle.py riboread-selection.py
do
    if [ ! -x "${p2s}/$s" ]; then echo "MISSING/NOT-EXEC script: ${p2s}/$s"; missing=1; fi
done
if [ ! -d "$folder" ]; then echo "MISSING folder: $folder"; missing=1; fi
if [ $missing -ne 0 ]
then
    echo ""
    echo "Refusing to submit until the above exist."
    exit 1
fi

### build the cell manifest
# One cell basename per line (e.g. SRR14783059_001), sorted, taken from the
# STAGE 1 output. Array task i handles line i, which is what makes the
# per-cell arrays reproducible and resumable.
manifest=${folder}/.cells.manifest
ls ${folder}/${lib}*_cbc.fastq.gz 2>/dev/null \
    | sed 's#.*/##; s#_cbc\.fastq\.gz$##' \
    | sort > "$manifest"
ncell=$(wc -l < "$manifest")
if [ "$ncell" -eq 0 ]
then
    echo "no ${lib}*_cbc.fastq.gz found in ${folder} -- has STAGE 1 been run?"
    exit 1
fi
nchunk=$(( (ncell + cells_per_chunk - 1) / cells_per_chunk ))

echo "library      : ${lib}"
echo "reference    : ${ref} (readlen ${n})"
echo "folder       : ${folder}"
echo "cells        : ${ncell}  (manifest: ${manifest})"
echo "gmap chunks  : ${nchunk} x ${cells_per_chunk} cells"
echo "concurrency  : ${maxpar} (per-cell stages), ${maxpar_gmap} (gmap)"
echo ""

SB=sbatch
if [ "${DRYRUN:-0}" != "0" ]; then SB="echo [dryrun] sbatch"; fi

# every array task resolves its own cell from the manifest
cellsel="cell=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" ${manifest})"

### STAGE 2 (trim) ----------------------------------------------------------
jtrim=""; jribo=""; jgmap=""; jb2bs=""; jb2bm=""
if [ "$START" -le 2 ]; then
jtrim=$(${SB} --parsable --export=All -N 1 -c 1 -t 05:00:00 --mem=10G \
    -a 1-${ncell}%${maxpar} -J trim-${lib} \
    -e ${folder}/trim-%a.err -o ${folder}/trim-%a.out \
    --wrap="${ml_trim}; ${cellsel}; ${p2s}/trim.sh ${folder}/\${cell}_cbc.fastq.gz ${folder} ${p2trimgalore} ${p2cutadapt}")
echo "STAGE 2 trim  : ${jtrim}  (array 1-${ncell}%${maxpar})"
else
echo "STAGE 2 trim  : SKIPPED (START=${START})"
fi

### STAGE 3 (ribo) ----------------------------------------------------------
# aftercorr: task i starts as soon as task i of trim is done, rather than
# waiting for the whole trim array -- the stages pipeline per cell.
if [ "$START" -le 3 ]; then
dep=""; [ -n "$jtrim" ] && dep="--dependency=aftercorr:${jtrim}"
jribo=$(${SB} --parsable --export=All -N 1 -c 8 -t 10:00:00 --mem=40G \
    -a 1-${ncell}%${maxpar} -J ribo-${lib} ${dep} \
    -e ${folder}/ribo-%a.err -o ${folder}/ribo-%a.out \
    --wrap="${ml_ribo}; ${cellsel}; ${p2s}/ribo-bwamem.sh ${riboref} ${folder}/\${cell}_cbc_trimmed_homoATCG.fq.gz ${folder}/\${cell}_cbc_trimmed_homoATCG ${p2bwa} ${p2samtools} y ${p2s}")
echo "STAGE 3 ribo  : ${jribo}  (array 1-${ncell}%${maxpar}) ${dep:-(no dependency)}"
else
echo "STAGE 3 ribo  : SKIPPED (START=${START})"
fi

### STAGE 4 (gmap, batched) --------------------------------------------------
# --mem=80G: the index is 53 GB and lives in shared memory that is charged to
# this job's cgroup, plus STAR's own working set. The upstream 40G is what
# would have OOM'd here.
if [ "$START" -le 4 ]; then
dep=""; [ -n "$jribo" ] && dep="--dependency=afterany:${jribo}"
jgmap=$(${SB} --parsable --export=All -N 1 -c 8 -t 24:00:00 --mem=80G \
    -a 1-${nchunk}%${maxpar_gmap} -J gmap-${lib} ${dep} \
    -e ${folder}/gmap-%a.err -o ${folder}/gmap-%a.out \
    --wrap="${ml_gmap}; ${p2s}/gmap_chunk.sh ${manifest} \${SLURM_ARRAY_TASK_ID} ${cells_per_chunk} ${folder} ${genome} ${p2star} ${p2samtools}")
echo "STAGE 4 gmap  : ${jgmap}  (array 1-${nchunk}%${maxpar_gmap}) ${dep:-(no dependency)}"
else
echo "STAGE 4 gmap  : SKIPPED (START=${START})"
fi

### STAGE 5a (b2bs) / 5b (b2bm) ---------------------------------------------
dep=""; [ -n "$jgmap" ] && dep="--dependency=afterany:${jgmap}"
jb2bs=$(${SB} --parsable --export=All -N 1 -c 1 -t 12:00:00 --mem=40G \
    -a 1-${ncell}%${maxpar} -J b2bs-${lib} ${dep} \
    -e ${folder}/b2bs-%a.err -o ${folder}/b2bs-%a.out \
    --wrap="${ml_b2b}; ${cellsel}; ${p2s}/deal_with_singlemappers.sh ${folder}/\${cell}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
echo "STAGE 5a b2bs : ${jb2bs}  (array 1-${ncell}%${maxpar}, afterany:${jgmap})"

jb2bm=$(${SB} --parsable --export=All -N 1 -c 1 -t 10:00:00 --mem=40G \
    -a 1-${ncell}%${maxpar} -J b2bm-${lib} ${dep} \
    -e ${folder}/b2bm-%a.err -o ${folder}/b2bm-%a.out \
    --wrap="${ml_b2b}; ${cellsel}; ${p2s}/deal_with_multimappers.sh ${folder}/\${cell}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
echo "STAGE 5b b2bm : ${jb2bm}  (array 1-${ncell}%${maxpar}, afterany:${jgmap})"

### STAGE 6 (cout) -----------------------------------------------------------
# Run once over the whole folder, after every cell's b2bs AND b2bm are done.
jcout=$(${SB} --parsable --export=All -c 8 -t 60:00:00 --mem=160G \
    --dependency=afterany:${jb2bs}:${jb2bm} -J cnt${folder} \
    -e cnt${folder}.err -o cnt${folder}.out --mail-type=END --mail-user=${email} \
    --wrap="${ca}; ${p2s}/countTables_2pickle_cellsSpliced.py ${folder} ${folder} vasa ${cellidori}")
echo "STAGE 6 cout  : ${jcout}  (afterany:${jb2bs}:${jb2bm})"

### STAGE 7 (pick) -----------------------------------------------------------
jpick=$(${SB} --parsable --export=All -c 1 -t 15:00:00 --mem=160G \
    --dependency=afterany:${jcout} -J pick${folder} \
    -e pick${folder}.err -o pick${folder}.out --mail-type=END --mail-user=${email} \
    --wrap="${ca}; ${p2s}/countTables_fromPickle.py ${folder}.pickle.gz ${folder} vasa y")
echo "STAGE 7 pick  : ${jpick}  (afterany:${jcout})"

echo ""
echo "submitted. cancel everything with:  scancel $(echo ${jtrim} ${jribo} ${jgmap} ${jb2bs} ${jb2bm} ${jcout} ${jpick} | tr -s ' ')"
