#!/bin/bash

################################################################################
# submit_vasaplate_map.sh
#
# PURPOSE
#   SLURM driver script for the VASA-plate scRNA-seq mapping pipeline.
#   It does not do any mapping itself; instead it submits a chain of sbatch
#   jobs (one per pipeline stage) that call helper scripts living in
#   ${p2s} (vasaplate_split), wiring them together with
#   `sbatch --dependency=afterany:<jobid>` so each stage only starts once
#   the previous one has finished (successfully or not).
#
# PIPELINE STAGES (submitted in order, one sbatch job per stage per lane
# except where noted):
#   1. extract  - extractBC.sh            : demultiplex reads / extract cell
#                                            barcodes (CBCs) from raw fastqs.
#                                            Only runs if fqext=y; the script
#                                            EXITS after submitting this job,
#                                            so you must re-run it afterwards
#                                            with fqext=n to continue the rest
#                                            of the pipeline once extraction
#                                            has finished.
#   2. trim     - trim.sh                 : adapter/quality trimming
#                                            (TrimGalore + cutadapt).
#   3. ribo     - ribo-bwamem.sh          : BWA-align reads against a
#                                            ribosomal RNA reference and
#                                            discard/flag rRNA reads.
#   4. gmap     - map_star.sh             : STAR-align non-ribosomal reads to
#                                            the genome.
#   5. b2bs/b2bm- deal_with_singlemappers.sh / deal_with_multimappers.sh
#                                          : assign uniquely-/multi-mapping
#                                            reads to exon/intron/tRNA
#                                            features using the BED
#                                            annotation (refBED). These two
#                                            run in parallel off the same
#                                            STAR output.
#   6. cout     - countTables_2pickle_cellsSpliced.py
#                                          : build per-cell, spliced/unspliced
#                                            -aware count tables (pickle),
#                                            after ALL lanes' b2bs/b2bm jobs
#                                            have finished.
#   7. pick     - countTables_fromPickle.py
#                                          : convert the pickle into the
#                                            final output count table(s).
#
# USAGE
#   ./submit_vasaplate_map.sh <lib> <genome> <readlen> <out> <folder> <fqext> <cellidori>
#   See the argument-count check below for the meaning of each of the 7
#   positional arguments.
#
#   Typical workflow: run once with fqext=y to submit barcode extraction and
#   let it finish, then run again with fqext=n (same other arguments) to
#   submit trim -> ribo -> gmap -> b2bs/b2bm -> cout -> pick for every lane
#   found under ${folder}.
#
# ENVIRONMENT / PATHS
#   All tool paths below (p2s, p2trimgalore, p2cutadapt, p2bwa, p2samtools,
#   p2star, p2bedtools) and the notification email are hardcoded for the
#   original Alemany-lab HPC environment and will need to be updated to
#   match wherever this is being run (e.g. the Crick HPC / your modified
#   pipeline setup).
#
# REFERENCES
#   Reference files (ribosomal fasta, STAR genome index, annotation BED) are
#   selected below based on the `ref` argument (MOUSE or HUMAN). Note the
#   usage message also advertises a MIXED option, but no MIXED branch is
#   actually implemented in the reference-selection block, so passing MIXED
#   will fall through without setting riboref/genome/refBED and will fail
#   the "genome not found" check.
################################################################################

### input paths (to modify by user)
p2s=/exports/ana-scarlab/aalemany/bin/vasaplate_split   # path to mapping scripts in your computer/HPC
p2trimgalore=/exports/ana-scarlab/bin/TrimGalore-0.6.6  # path to TrimGalore
p2cutadapt=/home/aalemany/anaconda3/bin/                # path to cutadapt
p2bwa=/exports/ana-scarlab/bin/bwa-0.7.17               # path to BWA
p2samtools=/exports/ana-scarlab/bin/samtools-1.11       #  path to samtools
p2star=/exports/ana-scarlab/bin/STAR-2.7.7a/bin/Linux_x86_64    # path to STAR
p2bedtools=/exports/ana-scarlab/bin/bedtools2/bin          # path to bedtools
email=a.alemany@lumc.nl                                 # email

### check input parameters
# Requires exactly 7 positional arguments (see PIPELINE STAGES / USAGE header
# above for the full workflow). If the count is wrong, print usage and exit
# without submitting anything.
if [ $# -ne 7 ]
then
    echo "Please, give:"
    echo "1) library name (prefix of the fastq files, name before _R1.fastq.gz and _R2.fastq.gz)"
    echo "2) genome: MOUSE /  HUMAN / MIXED "
    echo "3) read length (for MOUSE: 59, 74, 96, 246; for HUMAN: 74, 136; for MIXED: 74, 91, 135)"
    echo "4) prefix for output files"
    echo "5) folder for output files"
    echo "6) fastqfile extraction (y/n)"
    echo "7) cellID from filename or readname (f/r)"
    exit
fi

lib=$1
ref=$2
n=$3
out=$4
folder=$5
fqext=$6
cellidori=$7

### check existence of input fastq files
# Looks in the CURRENT working directory (not $folder) for files matching
# ${lib}*_R1*.fastq.gz / ${lib}*_R2*.fastq.gz. If a glob matches nothing, the
# `ls` call errors to stderr and $r1/$r2 stay empty, which is what the
# length checks below detect. Run this script from the directory holding the
# raw fastqs.
r1=$(ls ${lib}*_R1*.fastq.gz)
r2=$(ls ${lib}*_R2*.fastq.gz)
echo $r1 $r2
if [ ${#r1} == 0 ]
then
    echo "R1 fastq files not found"
    exit
fi
if [ ${#r2} == 0 ]
then
    echo "R2 fastq files not found"
    exit
fi

### check python version (we want version 3)
# `python` (not `python3`) must resolve to a Python 3 interpreter, since the
# downstream count-table scripts (countTables_2pickle_cellsSpliced.py,
# countTables_fromPickle.py) require it. Exits if the major version isn't 3.
v=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))' | awk -F "." '{print $1}')
if [ $v -ne "3" ]
then
    echo "python needs to be 3"
    exit
fi

### set references
# Picks riboref (rRNA fasta for the ribo-depletion alignment step),
# genome (STAR index directory, suffixed with the read-length argument `n`
# since a separate index is built per read length), and refBED (intron/
# exon/tRNA annotation used to assign reads to features) based on `ref`.
# NOTE: only MOUSE and HUMAN are implemented here, even though the usage
# message above also lists MIXED as a valid option — passing MIXED (or
# anything else) leaves riboref/genome/refBED unset and will fail the
# "genome not found" check just below.
if [[ $ref == "MOUSE" ]]
then
    riboref=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/rRNA_mouse_Rn45S.fa
    genome=/hpc/hub_oudenaarden/group_references/ensembl/99/mus_musculus/star_v273a_NOMASK_NOERCC_index_$n
    refBED=/hpc/hub_oudenaarden/aalemany/vasaseq/ref_seqs/Mus_musculus.GRCm38.99.homemade_IntronExonTrna.bed
elif [[ $ref == "HUMAN" ]]
then
    riboref=/exports/ana-scarlab/group_references/ensembl/human/99/unique_rRNA_human.fa
    genome=/exports/ana-scarlab/group_references/ensembl/human/99/star_v277a_index_$n
    refBED=/exports/ana-scarlab/group_references/ensembl/human/99/Homosapines_ensemble99.homemade_IntronExonTrna.bed
fi 

if [ ! -d $genome ]
then
    echo "genome not found"
    exit
fi

### extract cell barcodes
# STAGE 1 (extract): submits extractBC.sh via sbatch to pull cell barcodes
# out of the raw reads and write per-cell/per-lane *_cbc.fastq.gz files into
# ${folder}. `jcbc` captures the submitted job's numeric ID (parsed out of
# sbatch's "Submitted batch job <id>" output via `awk '{print $NF}'`) so
# later stages can depend on it.
# IMPORTANT: when fqext=y, the script submits this one job and then EXITS
# immediately — it does NOT continue on to trim/ribo/gmap/etc. in the same
# invocation. Wait for the extract job to finish, then re-run this script
# with fqext=n (and the same other arguments) to submit the rest of the
# pipeline against the resulting *_cbc.fastq.gz files.
# jcbc defaults to the placeholder value 1 for the fqext=n path below (no
# real extract job is submitted in that case, so the later
# --dependency=afterany:$jcbc on the trim job is effectively a no-op).
jcbc=1
if [ $fqext == "y" ]
then
    jobid=extract-${lib}
    jcbc=$(sbatch --export=All -c 1 -N 1 -J ${jobid} -e ${jobid}.err -o ${jobid}.out -t 48:00:00 --mem=10G --mail-type=END --mail-user=${email} --wrap="${p2s}/extractBC.sh ${lib} vasaplate ${p2s} ${folder}")
    jcbc=$(echo $jcbc | awk '{print $NF}')
    exit
fi

### main per-lane loop (only reached when fqext=n)
# Iterates over every ${folder}/${lib}*cbc.fastq.gz produced by stage 1
# (one file per lane) and, for each one, submits the trim -> ribo -> gmap ->
# {b2bs,b2bm} chain below. `lib` is reassigned inside the loop to the
# per-lane basename (strips the _cbc.fastq.gz suffix and any leading
# directory path), so job names/log files are per-lane. `jbeds[]` collects
# the b2bs/b2bm job IDs from every lane so the final count-table step can
# wait on all of them via a single --dependency=afterany:<id1>:<id2>:...
# string built further down.
jbeds=(); lane=0
for file in ${folder}/${lib}*cbc.fastq.gz
do
    lib=${file%_cbc.fastq.gz}
    lib=${lib#*/}
    echo $lib
    ### trim
    # STAGE 2 (trim): TrimGalore/cutadapt adapter+quality trimming of the
    # per-lane barcoded fastq; depends on the extract job (jcbc).
    jobid=trim-${lib}
    jtrim=2
    jtrim=$(sbatch --export=All -N 1 -J ${jobid} -e ${folder}/${jobid}.err -o ${folder}/${jobid}.out --dependency=afterany:$jcbc -t 05:00:00 --mem=10G --wrap="${p2s}/trim.sh ${folder}/${lib}_cbc.fastq.gz ${folder} ${p2trimgalore} ${p2cutadapt}")
    jtrim=$(echo $jtrim | awk '{print $NF}')

    ### ribo-map
    # STAGE 3 (ribo): BWA-align trimmed reads against the rRNA reference
    # (riboref) and split out non-ribosomal reads (*.nonRibo.fastq.gz) for
    # genome mapping; depends on the trim job (jtrim).
    jobid=ribo-${lib}
    jribo=3
    jribo=$(sbatch --export=All -c 8 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jtrim --wrap="${p2s}/ribo-bwamem.sh $riboref ${folder}/${lib}_cbc_trimmed_homoATCG.fq.gz ${folder}/${lib}_cbc_trimmed_homoATCG $p2bwa $p2samtools y $p2s")
    jribo=$(echo $jribo | awk '{print $NF}')

    ### map to genome
    # STAGE 4 (gmap): STAR-align the non-ribosomal reads to the genome
    # index; depends on the ribo job (jribo). Output BAM is
    # ..._nonRibo_E99_Aligned.out.bam, consumed by both b2bs and b2bm below.
    jobid=gmap-$lib
    jgmap=4
    jgmap=$(sbatch --export=All -c 8 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jribo --wrap="${p2s}/map_star.sh ${p2star} ${p2samtools} ${genome} ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo.fastq.gz ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_")
    jgmap=$(echo $jgmap | awk '{print $NF}')

    ### STAGE 5a (b2bs): assign UNIQUELY-mapping reads to exon/intron/tRNA
    # features via refBED; depends on the gmap job (jgmap). Its job ID is
    # appended to jbeds[] for the final count-table dependency.
    jobid=b2bs-$lib
    jb2bs=5
    jb2bs=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${folder}/${jobid}.err -t 12:00:00 --mem=40G --dependency=afterany:$jgmap --wrap="${p2s}/deal_with_singlemappers.sh ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
    jb2bs=$(echo $jb2bs | awk '{print $NF}')
    lane=$((lane+1))
    jbeds[$lane]=$jb2bs

    ### STAGE 5b (b2bm): assign MULTI-mapping reads to exon/intron/tRNA
    # features via refBED, run in parallel with b2bs off the same gmap
    # output; its job ID is also appended to jbeds[].
    jobid=b2bm-$lib
    jb2bm=6
    jb2bm=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jgmap --wrap="${p2s}/deal_with_multimappers.sh ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
    jb2bm=$(echo $jb2bm | awk '{print $NF}')
    lane=$((lane+1))
    jbeds[$lane]=$jb2bm
done

### count table
# STAGE 6 (cout): builds the ":jobid1:jobid2:..." dependency string `j` from
# every b2bs/b2bm job ID collected in jbeds[] across all lanes, then submits
# countTables_2pickle_cellsSpliced.py with --dependency=afterany${j} so it
# only starts once every lane's feature-assignment jobs are done. This
# produces a per-cell, spliced/unspliced-aware count table pickle
# (${folder}.pickle.gz), using $cellidori to decide whether cell IDs come
# from the filename or the read name.
jobid=cout-$lib
j=''
for k in ${jbeds[@]}
do
    j=${j}:$k
done
jcout=7
jcout=$(sbatch --export=All -c 8 -t 60:00:00 --mem=160G --dependency=afterany${j} -J cnt${folder} -e cnt${folder}.err -o cnt${folder}.out --mail-type=END --mail-user=${email} --wrap="${p2s}/countTables_2pickle_cellsSpliced.py ${folder} ${folder} vasa $cellidori";)
jcout=$(echo $jcout | awk '{print $NF}')

### STAGE 7 (pick): final step. Depends on the cout job (jcout); unpacks the
# pickle produced above into the final, human-usable count table file(s)
# via countTables_fromPickle.py. This is the last job in the chain — no
# further stage depends on jpick.
jobid=pick-$lib
jpick=8
jpick=$(sbatch --export=All -c 1 -t 15:00:00 --mem=160G --dependency=afterany:$jcout -J pick${folder} -e pick${folder}.err -o pick${folder}.out --mail-type=END --mail-user=${email} --wrap="${p2s}/countTables_fromPickle.py ${folder}.pickle.gz ${folder} vasa y")



