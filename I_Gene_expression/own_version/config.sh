#!/bin/bash
###############################################################################
# config.sh -- THE ONLY FILE YOU SHOULD NEED TO EDIT
#
# Everything in pipeline.sh reads its settings from here. Edit the values in
# the "EDIT ME" block, then run:
#
#     ./pipeline.sh check      # verify every path/tool exists (do this first!)
#     ./pipeline.sh step1      # ... and so on, one step at a time
#
# Nothing here runs any analysis; this file is only variable definitions.
###############################################################################

# =============================================================================
# EDIT ME -- your sample
# =============================================================================

# Prefix of your fastq files, i.e. the part BEFORE _R1/_R2.
# Files must be named:  ${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz
#                       ${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz
# If yours are named _1.fastq.gz / _2.fastq.gz, `./pipeline.sh check` will tell
# you and print the two `ln -s` commands that fix it.
SAMPLE="${SAMPLE:-ZHA9292A1}"

# Directory holding those two fastq files.
#
# The delivered files are ZHA9292A1_S181_L007_R{1,2}_001.fastq.gz, which neither
# this script nor concatenator.py's glob will match, and the delivery directory
# is not writable. FASTQ_DIR therefore holds _R1/_R2 symlinks to them:
#   /nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/PM26037/
#       20260720_LH00442_0273_B23TM55LT4/fastq/
FASTQ_DIR="${FASTQ_DIR:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037}"

# Where all output goes. Created if missing. Expect it to get large
# (roughly 3-4x your input fastq size once every stage has run) -- the two
# input fastqs are 30 GB, so budget ~120 GB here.
OUTDIR="${OUTDIR:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out}"

# Throwaway files that are NOT results. Currently just step4: STAR's
# --genomeLoad LoadAndExit / Remove calls still insist on writing a Log.out and
# Log.progress.out, which we discard. Keeping them out of OUTDIR means a run
# killed mid-step4 leaves no .star.XXXXXX dirs sitting next to the count tables.
# Anything here is safe to delete at any time; nothing downstream reads it.
#
# Default is the lab scratch, which is the same directory under either mount
# (/nemo/lab/turnerj/scratch/zhangg and /flask/scratch/turnerj/zhangg are one
# inode, verified 2026-07-26) and is writable via ACL, not group.
SCRATCH="${SCRATCH:-/nemo/lab/turnerj/scratch/zhangg/vasaseq}"

# =============================================================================
# EDIT ME -- read structure
# =============================================================================

# Uninformative 5' prefix present on BOTH mates, removed before anything else.
# This is the thing that makes your data different from the published VASA data.
# Set to 0 to disable.
#
# Verified against ZHA9292A1 by per-position base composition (100k reads):
#   R1 pos 1-21  = GAGTTCTACAGTCCGACGATC   (3' end of the RA5 adapter), 98-99%
#                  fixed at every position; pos 22 goes random.
#   R2 pos 1-21  = CCTTGGCACCCGAGAATTCCA   (revcomp of RA3, TGGAATTCTCGGGTGCC
#                  AAGG), same 98-99% fixed; pos 22 goes random.
# So 21 is right for both mates and SKIP5R1/SKIP5R2 are not needed.
SKIP5="${SKIP5:-21}"

# Per-mate overrides, if only ONE mate carries the prefix. Leave empty to use
# SKIP5 for both. e.g. SKIP5R1=0 SKIP5R2=21  -> strip R2 only.
SKIP5R1="${SKIP5R1:-}"
SKIP5R2="${SKIP5R2:-}"

# Barcode read layout, AFTER the SKIP5 prefix has been removed.
#   R1 = [SKIP5 junk][UMI  LEN_UMI][CBC  LEN_CBC][...]      when UMI_FIRST=yes
#   R1 = [SKIP5 junk][CBC  LEN_CBC][UMI  LEN_UMI][...]      when UMI_FIRST=no
#
# ZHA9292A1 is 6+6, NOT the CEL-seq2 6+8. Measured on R1 (absolute positions):
#   22-27  all 4096 6-mers observed, near-flat composition  -> UMI, 6 nt
#   28-33  strongly structured, 16 sequences dominate       -> CBC, 6 nt
#   34+    93% T, decaying over ~24 nt                      -> polyT
# The 8-mer at 28-35 is always <barcode>TT, i.e. the barcode plus the first two
# bases of the polyT -- that is what rules out an 8 nt barcode.
LEN_UMI="${LEN_UMI:-6}"
LEN_CBC="${LEN_CBC:-6}"
UMI_FIRST="${UMI_FIRST:-yes}"   # yes = UMI comes before the cell barcode (CEL-seq2 / VASA)

# Which mate holds the barcode, and which holds the cDNA.
BC_READ="${BC_READ:-R1}"
BIO_READ="${BIO_READ:-R2}"

# Cell-barcode whitelist: a 2-column TAB-separated file, "barcode<TAB>cellID".
# bc_celseq2.tsv does NOT apply here -- it is 384 x 8 nt. This library uses 16
# x 6 nt barcodes, read straight off the data (see README, "Barcodes"). They
# become cells 001-016 in the count tables, in the alphabetical order of the
# file, which is not abundance order -- the README table maps id -> read count.
BC_WHITELIST="${BC_WHITELIST:-/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/own_version/bc_PM26037_6nt.tsv}"

# Allowed mismatches when matching an observed barcode to the whitelist.
# Keep at 0: several barcode pairs in this set differ only in their last two
# bases (AGCTAG/AGCTCA, CATGAG/CATGCA, TGTCAC/TGTCGA), so a 1-mismatch variant
# would be ambiguous. concatenator.py drops ambiguous variants anyway, so
# raising this buys almost nothing and costs runtime.
BC_HAMMING="${BC_HAMMING:-0}"

# Stranded protocol? VASA is stranded -> y
STRANDED="${STRANDED:-y}"

# =============================================================================
# EDIT ME -- step 2 trimming
# =============================================================================
# Upstream's second trimming pass cuts at the first run of 6 identical bases
# and discards everything after it. On this library that truncates 49% of the
# reads that have no poly-A tail at all (median loss 120 of 130 bases), while
# never removing the adapter that is actually in the way. own_version/trim.sh
# replaces it; the benchmark behind these values is in
# ../../../data/PM26037/trimtest/ and the numbers are in README.md.
#
#   vasa    read-through adapter + real poly-A/poly-G runs  (default)
#   legacy  byte-for-byte upstream, for comparison
TRIM_MODE="${TRIM_MODE:-vasa}"

# The 3' read-through, MEASURED not assumed: 68,250 reads anchored on the 16 nt
# that is unambiguously revcomp(R1's 5' prefix), then a per-position consensus
# (>=97% agreement to position 55). It is revcomp(R1 5' prefix) followed by the
# Nextera mosaic end:
#     GATCGTCGGACTGTAGAACTC CTGTCTCTTATACACATCT
# TrimGalore auto-detects only the Nextera half, which starts 21 nt too late.
# Re-derive it for a new library before trusting it -- see README.md.
TRIM_ADAPTER3="${TRIM_ADAPTER3:-GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT}"

# Poly-A / poly-G removed only as REAL runs (20 nt, or 10 nt running off the
# read end), not as the 6-mers upstream used. Set either to "" to disable.
# Dropping TRIM_POLYA gains ~28% more uniquely mapped reads and is tempting --
# don't. Those extra reads are poly-A aligning to genomic A-tracts: the
# poly-A-only fraction of uniquely mapped reads goes 3.7% -> 16.6%, piled onto
# a few 100-kb bins.
# NB: written as an if, not "${TRIM_POLYA:-A{20}}". The `}` inside the default
# closes the parameter expansion early and you get a stray brace in the adapter
# ("A{20}}"), which cutadapt then reads as a repeat count it cannot parse. The
# `+x` test also lets TRIM_POLYA="" mean "disabled" rather than "use default".
[ -n "${TRIM_POLYA+x}" ] || TRIM_POLYA='A{20}'
[ -n "${TRIM_POLYG+x}" ] || TRIM_POLYG='G{20}'

# 5' poly-T: a read in the reverse orientation reads the poly-A tail as poly-T
# at its START, so this is applied as a 5' adapter (removes the match and
# everything before it). It deletes junk almost exclusively -- protein-coding
# exonic counts stay flat while ~3,400 non-exonic uniquely mapped reads per 300k
# go away -- and it is what makes the blank barcodes look blank again.
[ -n "${TRIM_POLYT5+x}" ] || TRIM_POLYT5='T{20}'

# Reads shorter than this after trimming are dropped. Upstream used 15; 20 is
# about the shortest that still maps somewhere believable.
TRIM_MINLEN="${TRIM_MINLEN:-20}"

# Anchor the trim on the cell barcode (pass 0, trim_bc_anchor.py). The 3' tail
# does not have to be recognised by its shape: step 1 already put the barcode
# and the UMI on the read name, so for any given read the 12 nt after the
# poly-A are a known literal string. Find it, drop everything from there to the
# 3' end, walk back over the poly-A. Nothing to configure, and no threshold to
# tune -- set this to "no" to switch pass 0 off entirely.
#
# Effect: barcode remnant in the output drops 13.6% -> ~2%, poly-A tails
# shorter than 20 nt get cleaned (which TRIM_POLYA structurally cannot do), and
# purity rises from 54.2% to 55.5% protein-coding exonic at equal yield.
TRIM_ANCHOR_BC="${TRIM_ANCHOR_BC:-yes}"


# Where cell IDs in the final count table come from:
#   f = from the filename  (ids look like "cells/MYSAMPLE_001")
#   r = from the read name (ids look like "001" -- cleaner for a single sample)
CELLID_FROM="${CELLID_FROM:-f}"

# =============================================================================
# EDIT ME -- references (these must already exist; the pipeline does NOT build them)
# =============================================================================

# ZHA9292A1 is mouse (fastq_screen on R2: 7.8% unique to MOUSE vs 0.6% HUMAN),
# so the human+mouse MIXED reference built for the published mixing control is
# the wrong reference here -- wrong species set AND wrong read length. All three
# paths below are GRCm39 + Ensembl 116, matching the nf-core rnaseq run so the
# two analyses stay comparable.
#
# Only ONE of the three still has a build script: the rRNA fasta, via
# ./build_rrna_reference.sh. The STAR index and the BED were made by
# build_mouse_reference.sh, which was deleted once its outputs existed -- so
# they exist but cannot be reproduced. Do not repeat that: new reference
# builders go in this directory, in git.
VASA_REF="${VASA_REF:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116}"

# STAR index. sjdbOverhang must be >= (biological read length after SKIP5) - 1:
# R2 is 151 nt, SKIP5 takes 21 -> 130 nt biological -> needs at least 129.
#
# We reuse the index the nf-core rnaseq FLASH-seq run already built rather than
# keeping a second 27 GB copy. It is the same genome.fa and the same Ensembl 116
# GTF, with a byte-identical junction set (562,855 junctions); only the overhang
# differs -- 150, for that run's 151 nt reads. That is fine here: an overhang
# LARGER than readLength-1 is harmless (STAR just stores flanking sequence it
# never uses), whereas one too small silently costs junction-spanning
# sensitivity. pipeline.sh check enforces the >= rule.
#
# This index was built by the nf-core rnaseq run, then moved (2026-07-24) out of
# its fragile results/ output dir into the curated reference library so it now
# survives cleaning results/ or re-running nf-core with a different --outdir.
# Named star_index_151_r116: read length 151 (sjdbOverhang 150), Ensembl rel 116.
#
# NOTE: build_mouse_reference.sh, which this line used to point at, was deleted
# once its outputs existed and is gone. If this index goes missing it has to be
# rebuilt by hand (STAR --runMode genomeGenerate over the same genome.fa +
# Ensembl 116 GTF), or via the nf-core rnaseq run that originally made it.
# The rRNA fasta is the one piece with a tracked builder --
# see build_rrna_reference.sh, and keep new reference builders in the repo.
STAR_INDEX="${STAR_INDEX:-/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116}"

# rRNA fasta, bwa-indexed (needs the .amb/.ann/.bwt/.pac/.sa files beside it).
#
# v2 (2026-07-26), built by ./build_rrna_reference.sh -- read that script's
# header for the full rationale. Short version: v1 was Ensembl-only, and because
# the rDNA array is collapsed in the GRCm39 primary assembly, Ensembl 116 has no
# Rn28s / Rn45s / Rn5-8s gene at all. Measured on 130 nt reads tiled across the
# true subunits, v1 caught 0 of 71 28S reads and 1 of 60 5'ETS+ITS1 reads; 28S
# alone is 4,730 of the 13,400 nt transcript and the most abundant rRNA by mass.
# v2 = v1 (all 356 Ensembl seqs, unchanged) + the NCBI curated 47S transcript
# (BK000964.3:1-13403), which is what the paper's Methods actually used
# ("mouse or human rRNA (National Center for Biotechnology Information)").
#
# v2 catches 159/159 of those same simulated subunit reads. Specificity is
# unchanged: 20,000 simulated protein-coding-exon reads give 0 false positives
# under BOTH v1 and v2, at 130 nt and at 50 nt -- so the added 13.4 kb costs
# nothing. v1 is kept beside it as unique_rRNA_mouse.fa for comparison.
RRNA_FASTA="${RRNA_FASTA:-${VASA_REF}/unique_rRNA_mouse.v2.fa}"

# Intron/exon/tRNA BED with biotype embedded in the gene name.
REF_BED="${REF_BED:-${VASA_REF}/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.bed}"

# =============================================================================
# EDIT ME -- how much machine to use
# =============================================================================

# This library is 16 cells, not 384, so the two halves of the pipeline want
# different widths and NCORES is a compromise between them:
#   step1     shards the fastq ~100 ways, so it scales past 16.
#   steps 2/3/5 run one worker per cell and CANNOT use more than 16.
# 16 fits both. On a wider node, override for step1 only:
#   NCORES=32 ./pipeline.sh step1   &&   ./pipeline.sh step2
NCORES="${NCORES:-16}"              # parallel workers (step1 shards, steps 2/3/5 cells)
SHARD_READS="${SHARD_READS:-2000000}"   # reads per shard in step1. More shards = better balance.
# step4 loops over cells serially with the genome held in shared memory, so
# STAR_THREADS is independent of NCORES and should just be the core count.
# Each of the 16 cells here carries ~12M reads (a 384-well plate cell carries
# ~0.5M), so this matters much more than it did upstream.
STAR_THREADS="${STAR_THREADS:-16}"      # threads for one STAR run

# =============================================================================
# You should not normally need to change anything below this line.
# =============================================================================

# --- where the unchanged upstream helper scripts live ------------------------
# own_version deliberately does NOT copy these: the awk in deal_with_*mappers.sh
# and the count-table python are long, fragile, and identical to upstream. We
# call them in place. Only concatenator.py was forked (it lives next to this
# file) because the SKIP5 change had to go inside it.
VASA_SCRIPTS="${VASA_SCRIPTS:-/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/a_Mapping}"

# our own forked scripts
OWN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONCATENATOR="${OWN_DIR}/concatenator.py"
TRIM_SH="${TRIM_SH:-${OWN_DIR}/trim.sh}"
TRIM_ANCHOR_PY="${TRIM_ANCHOR_PY:-${OWN_DIR}/trim_bc_anchor.py}"
STEP2_REPORT="${STEP2_REPORT:-${OWN_DIR}/step2_report.py}"
STEP3_REPORT="${STEP3_REPORT:-${OWN_DIR}/step3_report.py}"

# --- cluster tools (EasyBuild module tree) -----------------------------------
EBROOT=/camp/apps/eb/software
P2TRIMGALORE=${EBROOT}/Trim_Galore/0.6.2-foss-2018b-Python-3.6.6
P2CUTADAPT=${EBROOT}/cutadapt/1.18-foss-2018b-Python-3.6.6/bin
P2BWA=${EBROOT}/BWA/0.7.17-GCC-10.3.0/bin
P2SAMTOOLS=${EBROOT}/SAMtools/1.11-GCC-10.2.0/bin
P2STAR=${EBROOT}/STAR/2.7.7a-GCC-10.2.0/bin
P2BEDTOOLS=${EBROOT}/BEDTools/2.30.0-GCC-11.2.0/bin

# Modules are loaded PER STEP, never all together: Trim_Galore/0.6.2 is a
# foss-2018b (GCC 7.3.0) build, and loading it alongside the others drags
# libstdc++ backwards and breaks STAR and bedtools with "GLIBCXX not found".
ML_INIT="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=${EBROOT%/software}/modules/all"
ML_TRIM="${ML_INIT}; module load Trim_Galore/0.6.2-foss-2018b-Python-3.6.6"
ML_STAR="${ML_INIT}; module load STAR/2.7.7a-GCC-10.2.0 SAMtools/1.11-GCC-10.2.0"
ML_BED="${ML_INIT}; module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0"

# python env holding numpy/pandas/pysam/multiprocess (the system python3 has none)
CONDA_ENV="${CONDA_ENV:-/nemo/lab/turnerj/working/guangxin/envs/vasa}"

# cutadapt for step 2's SECOND pass. The module tree stops at 1.18 (2018),
# which has neither per-adapter `;min_overlap=` nor `--poly-a`, so cutadapt 5.1
# was pip-installed into the conda env above (2026-07-26). This is a capability
# upgrade, not a behaviour change: given the upstream parameters, 1.18 and 5.1
# produce byte-identical output (checked by md5 on two cells). It is called by
# absolute path so it never has to win a PATH fight with the Trim_Galore
# module -- TrimGalore's own pass still drives the module's cutadapt 1.18.
TRIM_CUTADAPT="${TRIM_CUTADAPT:-${CONDA_ENV}/bin/cutadapt}"
TRIM_PYTHON="${TRIM_PYTHON:-${CONDA_ENV}/bin/python3}"   # runs pass 0
CONDA_ACTIVATE="source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate ${CONDA_ENV}"
# bwa/samtools are called by absolute path, so loading modules + conda together
# is safe here: conda's python leads on PATH, the binaries resolve absolutely.
ML_RIBO="${ML_INIT}; module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0; ${CONDA_ACTIVATE}"

# --- derived paths (do not edit) ---------------------------------------------
CELLDIR="${OUTDIR}/cells"          # one fastq per cell lives here, and all
                                   # per-cell intermediates after that
LOGDIR="${OUTDIR}/logs"            # per-step logs
