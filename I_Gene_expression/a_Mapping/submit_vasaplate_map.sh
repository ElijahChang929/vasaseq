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
#   Localised to the Crick HPC (NEMO/CAMP) account of zhangg on 2026-07-17.
#   The tool paths (p2trimgalore, p2cutadapt, p2bwa, p2samtools, p2star,
#   p2bedtools) now point at the cluster's EasyBuild module tree, and p2s /
#   email point at this checkout and this user. Each stage loads its own
#   modules (ml_trim / ml_ribo / ml_gmap / ml_b2b) — see the note by those
#   definitions for why they must NOT be loaded all together.
#
# REFERENCES
#   Reference files (ribosomal fasta, STAR genome index, annotation BED) are
#   selected below based on the `ref` argument (MOUSE / HUMAN / MIXED) and
#   rooted at $VASA_REFS. The MIXED branch that the usage text advertised but
#   upstream never implemented is now present — it is the branch this account's
#   only VASA-plate library (SRR14783059, the HEK293T-mESC species-mixing
#   control) actually needs. NONE of the reference files exist yet; they must
#   be built. A preflight block below checks every tool, reference and helper
#   script and refuses to submit if anything is missing.
################################################################################

### input paths (to modify by user)
# Adapted for the Crick HPC (NEMO/CAMP) account of zhangg on 2026-07-17.
# Every tool below is provided by the cluster's EasyBuild module tree — nothing
# needs installing. Versions match the published pipeline almost exactly:
# STAR 2.7.7a, bwa 0.7.17 and samtools 1.11 are the exact upstream versions;
# bedtools is 2.30.0 (upstream unpinned) and TrimGalore is 0.6.2 (upstream 0.6.6
# is not on the cluster; 0.6.2 is the closest available).
p2s=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/a_Mapping  # this dir: holds extractBC.sh, extractBC_parallel.sh, concatenator.py, bc_celseq2.tsv, countTables_*.py
email=zhangg@crick.ac.uk                      # email

# STAGE 1 (extract) parallelism. The demultiplex step is single-core in
# concatenator.py, so a full plate takes hours. extractBC_parallel.sh runs the
# SAME concatenator.py over read-shards, ncbc-way in parallel, and merges the
# per-cell outputs — output is equivalent to serial extractBC.sh (verified by
# extractBC_parallel.verify.sh; no change to any python). Set ncbc=1 to fall
# back to the serial path. reads-per-shard tunes load balancing / shard count.
ncbc=16                                                 # concatenators to run concurrently in stage 1
ncbc_shard=2000000                                      # reads per shard (~112 shards for a 224M-read plate; keeps all 16 cores fed)

EBROOT=/camp/apps/eb/software
p2trimgalore=${EBROOT}/Trim_Galore/0.6.2-foss-2018b-Python-3.6.6   # path to TrimGalore (trim_galore sits in the dir root, not bin/)
p2cutadapt=${EBROOT}/cutadapt/1.18-foss-2018b-Python-3.6.6/bin     # path to cutadapt (the one Trim_Galore/0.6.2 is built against)
p2bwa=${EBROOT}/BWA/0.7.17-GCC-10.3.0/bin                          # path to BWA
p2samtools=${EBROOT}/SAMtools/1.11-GCC-10.2.0/bin                  #  path to samtools
p2star=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin                        # path to STAR
p2bedtools=${EBROOT}/BEDTools/2.30.0-GCC-11.2.0/bin                # path to bedtools

### per-stage module loads
# These binaries need their EasyBuild runtime libraries on LD_LIBRARY_PATH, so
# each sbatch job below loads its own modules via ${ml_*} before running.
#
# They are deliberately loaded PER STAGE and never all at once: Trim_Galore/0.6.2
# is a foss-2018b (GCC 7.3.0) build, and loading it alongside the others drags
# libstdc++ back to 7.3.0, which breaks STAR and bedtools at runtime with
# "GLIBCXX_3.4.26/3.4.29 not found" and silently downgrades cutadapt. Each of the
# four groupings below was verified to work in a clean shell; the union does not.
ml_init="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=${EBROOT%/software}/modules/all"
ml_trim="${ml_init}; module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6"
ml_gmap="${ml_init}; module load STAR/2.7.7a-GCC-10.2.0 SAMtools/1.11-GCC-10.2.0"
ml_b2b="${ml_init}; module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0"

### python env for the four Python-dependent stages
# The pipeline's own scripts need packages the system python3 lacks:
#   extract (concatenator.py)            -> numpy, pandas
#   ribo    (riboread-selection.py)      -> numpy, pandas, PYSAM
#   cout    (countTables_2pickle...py)   -> numpy, pandas, MULTIPROCESS
#   pick    (countTables_fromPickle.py)  -> numpy, pandas
# All four live in a dedicated conda env (built 2026-07-17, not a shared env):
#   python 3.10, numpy 2.2, pandas 2.3 (read_csv import verified), pysam 0.24,
#   multiprocess 0.70. The scripts use `#!/usr/bin/env python3`, so `ca` must be
#   sourced before they run to put this env's python3 first on PATH.
ca="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate /nemo/lab/turnerj/working/guangxin/envs/vasa"
# ribo runs bwa/samtools (by absolute path, from modules) AND riboread-selection.py
# (needs pysam), so it loads BOTH the modules and the conda env.
ml_ribo="${ml_init}; module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0; ${ca}"

### reference root
# The rRNA fasta, STAR index and IntronExonTrna BED are NOT in this account and
# are not shipped by the repo — they must be built (see a_Mapping/README.md for
# the BED format). Only the mouse GRCm39 genome fasta exists, under
# /nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39.
# NB the published pipeline is built on Ensembl 99 = GRCm38, not GRCm39; using
# GRCm39 here would silently disagree with the published count tables.
VASA_REFS=${VASA_REFS:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq}

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
# There is no bare `python` on this cluster (only `python3`), so the upstream
# check errored out with "command not found" and an empty $v before it could
# ever report a version. Use whichever of python3/python resolves.
PY=$(command -v python3 || command -v python)
if [ -z "$PY" ]
then
    echo "no python interpreter found"
    exit 1
fi
v=$($PY -c 'import sys; print(sys.version_info[0])')
if [ "$v" -ne 3 ]
then
    echo "python needs to be 3 (found $v at $PY)"
    exit 1
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
    # MIXED was advertised in the usage text but never implemented upstream.
    # It is the branch this account's only VASA-plate library actually needs:
    # SRR14783059 is the HEK293T-mESC species-mixing control, so it must map
    # against a concatenated human+mouse reference with species-prefixed
    # contigs, NOT against either species alone.
    riboref=${VASA_REFS}/mixed/unique_rRNA_human_mouse.fa
    genome=${VASA_REFS}/mixed/star_index_$n
    refBED=${VASA_REFS}/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.bed
fi

if [ -z "$genome" ]
then
    echo "unknown reference '$ref' (expected MOUSE / HUMAN / MIXED)"
    exit 1
fi

### preflight: refuse to submit against paths that do not exist
# The upstream script only checked $genome, so a missing rRNA fasta, BED or tool
# surfaced hours later as a failed sbatch job. Check everything up front instead.
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
extract_script=extractBC.sh; [ "${ncbc}" -gt 1 ] && extract_script=extractBC_parallel.sh
for s in ${extract_script} trim.sh ribo-bwamem.sh map_star.sh deal_with_singlemappers.sh deal_with_multimappers.sh countTables_2pickle_cellsSpliced.py countTables_fromPickle.py
do
    if [ ! -f "${p2s}/$s" ]; then echo "MISSING script:    ${p2s}/$s"; missing=1; fi
done
if [ $missing -ne 0 ]
then
    echo ""
    echo "Refusing to submit until the above exist."
    echo "  tools come from the module tree at ${EBROOT} (override: \$EBROOT)"
    echo "  refs  are rooted at \$VASA_REFS = ${VASA_REFS}"
    exit 1
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
    if [ "${ncbc}" -gt 1 ]
    then
        # PARALLEL demultiplex: allocate ncbc cpus and run extractBC_parallel.sh,
        # which fans concatenator.py out across read-shards ncbc-way. Same per-cell
        # output as the serial path, produced in a fraction of the wall time. Mem is
        # bumped because each concurrent concatenator holds its own barcode tables.
        excmd="${p2s}/extractBC_parallel.sh ${lib} vasaplate ${p2s} ${folder} ${ncbc} ${ncbc_shard}"
        jcbc=$(sbatch --export=All -c ${ncbc} -N 1 -J ${jobid} -e ${jobid}.err -o ${jobid}.out -t 48:00:00 --mem=64G --mail-type=END --mail-user=${email} --wrap="${ca}; ${excmd}")
    else
        # SERIAL fallback (ncbc=1): the original single-core extractBC.sh.
        jcbc=$(sbatch --export=All -c 1 -N 1 -J ${jobid} -e ${jobid}.err -o ${jobid}.out -t 48:00:00 --mem=10G --mail-type=END --mail-user=${email} --wrap="${ca}; ${p2s}/extractBC.sh ${lib} vasaplate ${p2s} ${folder}")
    fi
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
    jtrim=$(sbatch --export=All -N 1 -J ${jobid} -e ${folder}/${jobid}.err -o ${folder}/${jobid}.out --dependency=afterany:$jcbc -t 05:00:00 --mem=10G --wrap="${ml_trim}; ${p2s}/trim.sh ${folder}/${lib}_cbc.fastq.gz ${folder} ${p2trimgalore} ${p2cutadapt}")
    jtrim=$(echo $jtrim | awk '{print $NF}')

    ### ribo-map
    # STAGE 3 (ribo): BWA-align trimmed reads against the rRNA reference
    # (riboref) and split out non-ribosomal reads (*.nonRibo.fastq.gz) for
    # genome mapping; depends on the trim job (jtrim).
    jobid=ribo-${lib}
    jribo=3
    jribo=$(sbatch --export=All -c 8 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jtrim --wrap="${ml_ribo}; ${p2s}/ribo-bwamem.sh $riboref ${folder}/${lib}_cbc_trimmed_homoATCG.fq.gz ${folder}/${lib}_cbc_trimmed_homoATCG $p2bwa $p2samtools y $p2s")
    jribo=$(echo $jribo | awk '{print $NF}')

    ### map to genome
    # STAGE 4 (gmap): STAR-align the non-ribosomal reads to the genome
    # index; depends on the ribo job (jribo). Output BAM is
    # ..._nonRibo_E99_Aligned.out.bam, consumed by both b2bs and b2bm below.
    jobid=gmap-$lib
    jgmap=4
    jgmap=$(sbatch --export=All -c 8 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jribo --wrap="${ml_gmap}; ${p2s}/map_star.sh ${p2star} ${p2samtools} ${genome} ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo.fastq.gz ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_")
    jgmap=$(echo $jgmap | awk '{print $NF}')

    ### STAGE 5a (b2bs): assign UNIQUELY-mapping reads to exon/intron/tRNA
    # features via refBED; depends on the gmap job (jgmap). Its job ID is
    # appended to jbeds[] for the final count-table dependency.
    jobid=b2bs-$lib
    jb2bs=5
    jb2bs=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${folder}/${jobid}.err -t 12:00:00 --mem=40G --dependency=afterany:$jgmap --wrap="${ml_b2b}; ${p2s}/deal_with_singlemappers.sh ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
    jb2bs=$(echo $jb2bs | awk '{print $NF}')
    lane=$((lane+1))
    jbeds[$lane]=$jb2bs

    ### STAGE 5b (b2bm): assign MULTI-mapping reads to exon/intron/tRNA
    # features via refBED, run in parallel with b2bs off the same gmap
    # output; its job ID is also appended to jbeds[].
    jobid=b2bm-$lib
    jb2bm=6
    jb2bm=$(sbatch --export=All -c 1 -N 1 -J $jobid -o ${folder}/${jobid}.err -t 10:00:00 --mem=40G --dependency=afterany:$jgmap --wrap="${ml_b2b}; ${p2s}/deal_with_multimappers.sh ${folder}/${lib}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam ${refBED} y ${p2samtools} ${p2bedtools}")
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
jcout=$(sbatch --export=All -c 8 -t 60:00:00 --mem=160G --dependency=afterany${j} -J cnt${folder} -e cnt${folder}.err -o cnt${folder}.out --mail-type=END --mail-user=${email} --wrap="${ca}; ${p2s}/countTables_2pickle_cellsSpliced.py ${folder} ${folder} vasa $cellidori";)
jcout=$(echo $jcout | awk '{print $NF}')

### STAGE 7 (pick): final step. Depends on the cout job (jcout); unpacks the
# pickle produced above into the final, human-usable count table file(s)
# via countTables_fromPickle.py. This is the last job in the chain — no
# further stage depends on jpick.
jobid=pick-$lib
jpick=8
jpick=$(sbatch --export=All -c 1 -t 15:00:00 --mem=160G --dependency=afterany:$jcout -J pick${folder} -e pick${folder}.err -o pick${folder}.out --mail-type=END --mail-user=${email} --wrap="${ca}; ${p2s}/countTables_fromPickle.py ${folder}.pickle.gz ${folder} vasa y")



