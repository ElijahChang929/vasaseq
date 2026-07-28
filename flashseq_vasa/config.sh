#!/bin/bash
###############################################################################
# code/flashseq_vasa/config.sh -- THE ONLY FILE YOU SHOULD NEED TO EDIT
#
# Settings for pipeline_fs.sh, which pushes the ten FLASH-seq libraries through
# the VASA count-table path (steps 4-7) so that FLASH-seq and VASA-seq are
# quantified by the SAME code and can be compared as one measurement rather than
# two.
#
# Same arrangement as code/I_Gene_expression/own_version/config.sh: nothing here
# runs anything, it is only variable definitions, and pipeline_fs.sh reads every
# path from here.
#
#     ./pipeline_fs.sh check      # verify every path/tool first
#     FSV_ARM=native  ./pipeline_fs.sh prep map assign pickle tables
#     FSV_ARM=vasalen ./pipeline_fs.sh prep map assign pickle tables
#
# WHAT IS DIFFERENT FROM THE VASA SIDE, AND WHY -- all seven of these were
# established by the A9 dry run and are written up in NOUMI_PATH.md. They are
# settings, not discoveries, and every one of them is load-bearing:
#
#   1. R1 ONLY, SINGLE-END. Paired-end input cannot reach step 6 at all:
#      `bedtools bamtobed` appends the mate suffix /1,/2 to the read name, and
#      deal_with_*.sh appends ";CG:<cigar>;nM:<n>" to the QNAME BEFORE bamtobed
#      runs, so the suffix lands on the end of the nM value and step 6's
#      int('0/2') raises. Measured: 1,033,528 of 1,033,528 PE rows unparseable,
#      0 of 534,978 SE rows. a_Mapping/ is published code and this is its
#      contract, so the fix is to map R1 only -- which is also what keeps the
#      unit of measurement (one read, one verdict) identical to VASA's, and is
#      the same choice code/flashseq/05_rrna_bwa.sh made for the rRNA leg.
#   2. NO STEP 1, NO STEP 2's barcode anchor, NO STEP 3. There is no cell
#      barcode and no UMI to extract, and rRNA is deliberately NOT depleted
#      here: the rRNA leg of this comparison is measured by the bwa route in
#      res/flashseq/rrna_bwa.tsv, and the annotation route gives 0.80% against
#      that route's 3.50-6.44% because 99.7% of rRNA features are under 151 nt
#      and Ensembl lacks the 47S unit (NOUMI_PATH.md 5c).
#   3. protocol = smartseq_noUMI on steps 6 and 7.
#   4. stranded = n. FLASH-seq is unstranded -- its forward strand carries
#      49.1-50.5% of ribosomal reads, measured -- so y would halve every
#      biotype figure. VASA runs y, correctly, on its own side.
#   5. filt_unigenes = n on step 7. At y the threshold is
#      max(5, round(0.01*ncols)), written for a 384-well plate; at 10 columns it
#      demands 5 of them, i.e. 50%.
#   6. CELLID_FROM = f, and the per-library stem therefore contains the literal
#      '_cbc' -- step 6 mode 'f' does cellfile[:cellfile.index('_cbc')] and
#      raises ValueError without it. Mode 'r' would need ';SM:<lib>' injected
#      into ~29 M read names per library. The token is in the FILENAME, so
#      nothing in the data is touched. Cost: the column name comes out as
#      'cells/<LIB>' and the prefix is stripped by build_analysis_tables.py.
#   7. ReadCounts ONLY. With one literal UMI key 'A', UFICounts degenerates to a
#      0/1 detection mask and TranscriptCounts == UFICounts elementwise
#      (bc2trans(1) = 1.0 at K = 4**len('A') = 4). Neither carries abundance
#      information. There is also NO UMI CEILING on this path, so VASA's
#      ceiling filter has no analogue here.
###############################################################################

# =============================================================================
# EDIT ME -- which arm
# =============================================================================
# TWO ARMS OF THE SAME TEN LIBRARIES. Identical downstream settings; the ONLY
# difference is the length of the reads that go into STAR.
#
#   native   adapter-trimmed R1 at its natural length (up to 151 nt). This is
#            what the FLASH-seq run actually yields, and it is the arm the
#            nf-core cross-check must be compared against -- same reads, same
#            length, different quantifier, so the pipeline-vs-protocol contrast
#            is isolated.
#
#   vasalen  the native arm's reads hard-trimmed, per read, to a length drawn
#            from VASA's own STAR-input length distribution. This is the
#            READ-LENGTH CONTROL, and it exists because of a measured
#            confound: step 5 tags a read jS:IN only if the read is CONTAINED in
#            the feature (readstart >= refstart && readend <= refend), and in
#            the v2 BED 98.5% of tRNA features, 99.2% of miRNA, 96.5% of snoRNA
#            and 84.6% of snRNA are shorter than one 151 nt read. Step 6 then
#            keeps a non-spliceable biotype ONLY when jS == IN. So at 151 nt the
#            short non-poly-A species are suppressed by read length alone,
#            independently of any biology -- in the A9 dry run all six
#            tRNA-overlapping reads were discarded by containment and the tRNA
#            table came out empty, while VASA cell 005 had 229 tRNA rows.
#            Comparing native to vasalen separates the structural effect from
#            genuine poly-A depletion. Nothing else changes between the arms.
FSV_ARM="${FSV_ARM:-native}"

# =============================================================================
# EDIT ME -- libraries and where things live
# =============================================================================
FSV_ROOT="${FSV_ROOT:-/nemo/lab/turnerj/working/guangxin/vasaseq}"
FSV_CODE="${FSV_CODE:-$FSV_ROOT/code/flashseq_vasa}"

# Read-only delivery directory, run RN26038. Files are <LIB>_S<n>_L007_R1_001.fastq.gz.
FSV_FASTQ="${FSV_FASTQ:-/nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/RN26038/20260325_LH00442_0237_B23GT7GLT3/fastq}"

# The ten libraries, in titration order rather than lexical order (A10 sorts
# between A1 and A2). This is the SAME order and the same set as
# code/flashseq/sample_metadata.tsv, which is the only place the input amounts
# live: A1/A2 = 30 ng, A3/A4 = 3 ng, A5/A6 = 1.5 ng, A7/A8 = 60 pg,
# A9/A10 = 30 pg. A9/A10 is the single-cell-equivalent rung and the comparison
# point that matters. A8 is qc_verdict=exclude (18.3% human CALB1, a well
# effect at H:1) and is STILL PROCESSED here, exactly like the others -- the
# verdict filters interpretation, not QC.
FSV_LIBS="${FSV_LIBS:-ZHA8833A1 ZHA8833A2 ZHA8833A3 ZHA8833A4 ZHA8833A5 ZHA8833A6 ZHA8833A7 ZHA8833A8 ZHA8833A9 ZHA8833A10}"

# Intermediates (trimmed fastq, BAM, BED) -- tens of GB per arm, so they live on
# scratch, not in /nemo/lab (104 T at 92% full). Same scratch as
# own_version/config.sh: /nemo/lab/turnerj/scratch/zhangg and
# /flask/scratch/turnerj/zhangg are one inode.
FSV_SCRATCH="${FSV_SCRATCH:-/nemo/lab/turnerj/scratch/zhangg/vasaseq/flashseq_vasa}"

# Per-arm working directory. Step 6 is run from here with the relative folder
# name 'cells', so the column name is 'cells/<LIB>' and does not carry an
# absolute path (the same reason own_version/pipeline.sh cds to OUTDIR).
FSV_OUTDIR="${FSV_OUTDIR:-$FSV_SCRATCH/$FSV_ARM}"
FSV_CELLDIR="${FSV_CELLDIR:-$FSV_OUTDIR/cells}"
FSV_LOGDIR="${FSV_LOGDIR:-$FSV_OUTDIR/logs}"

# Where the count tables and reports that OTHER work reads are published. Small
# enough for /nemo/lab; this is the only thing that leaves scratch.
FSV_RES="${FSV_RES:-$FSV_ROOT/res/flashseq_vasa}"
FSV_TABLES="${FSV_TABLES:-$FSV_ROOT/data/flashseq_vasa}"

# Sample name step 6/7 prefix their outputs with. One per arm, so the two arms'
# tables cannot be confused for each other.
FSV_SAMPLE="${FSV_SAMPLE:-FS_$FSV_ARM}"

# =============================================================================
# EDIT ME -- the two arms' read handling
# =============================================================================
# Adapter trim, BOTH arms. This is the cutadapt call TrimGalore itself issued in
# the nf-core run, taken verbatim from results/trimgalore/*_trimming_report.txt,
# exactly as code/flashseq/05_rrna_bwa.sh already reuses it -- one recipe, not
# two. Trimming is not optional: mapping the raw fastq in the dry run returned
# 42.17% "unmapped: too short", consistent with the 55.1% adapter read-through
# measured over the whole ZHA8833A1 library.
FSV_TRIM_ADAPTER="${FSV_TRIM_ADAPTER:-CTGTCTCTTATA}"   # Nextera mosaic end
FSV_TRIM_Q="${FSV_TRIM_Q:-20}"
FSV_TRIM_OVERLAP="${FSV_TRIM_OVERLAP:-1}"
FSV_TRIM_ERR="${FSV_TRIM_ERR:-0.1}"
FSV_TRIM_MINLEN="${FSV_TRIM_MINLEN:-20}"
# ONE DELIBERATE DIFFERENCE FROM TrimGalore, stated because it changes the
# denominator: TrimGalore ran --paired and dropped a PAIR when EITHER mate fell
# under 20 bp. We map R1 only, so -m applies to R1 alone and a read whose mate
# was too short is kept. That is the right choice for an R1-only measurement,
# but it means our read count is >= TrimGalore's; the reconciliation table
# reports both.
#
# The reproduction is CHECKED, not assumed: prep compares its own R1
# "Reads with adapter" rate against the saved full-library TrimGalore report and
# fails if they differ by more than FSV_TRIM_TOL percentage points. Both run
# over the WHOLE library, so this should agree to well under a point.
FSV_TRIM_TOL="${FSV_TRIM_TOL:-0.5}"

# --- the vasalen arm's hard trim --------------------------------------------
# HOW THE LENGTHS ARE MATCHED, and why it is per-read rather than one number.
#
# A single fixed length at VASA's median would NOT be a control. VASA's median
# STAR-input read is ~110-130 nt, which still exceeds nearly every tRNA, miRNA
# and snoRNA feature, so trimming to the median would recover nothing and the
# arm would prove nothing. What lets VASA satisfy jS:IN on a short feature is
# the SHORT TAIL of its distribution. Reproducing only the median would drop
# exactly the part of the distribution that matters.
#
# So the whole distribution is matched: measure_vasa_readlen.py builds a
# 10,000-entry quantile lookup table from the pooled STAR-input read lengths of
# VASA's 12 REAL cells (blanks 001/014/015/016 excluded -- they are not
# biology), and trim_to_vasalen.sh draws one length per FLASH-seq read from that
# table and truncates the read to it. The result reproduces VASA's read-length
# composition, short tail included.
#
# Truncation is from the 3' end (the first L bases are kept), so the read's 5'
# start -- the only end STAR anchors on for an unspliced read -- is unchanged.
# A read already shorter than its drawn L is kept whole and counted; the trim
# script reports how often that happened.
#
# Reproducible rather than random: awk's srand(FSV_VASALEN_SEED) over a fixed
# input file in fixed order gives the same output every run. Re-running on the
# same fastq reproduces the same trimmed fastq bit for bit; the seed is stated
# so the draw is auditable.
FSV_VASALEN_LUT="${FSV_VASALEN_LUT:-$FSV_RES/vasa_starinput_len_lut.txt}"
FSV_VASALEN_SEED="${FSV_VASALEN_SEED:-20260728}"

# The VASA run whose length distribution is being matched, and the cells to
# pool. 001/014/015/016 are the four confirmed blanks.
FSV_VASA_OUT="${FSV_VASA_OUT:-$FSV_ROOT/data/PM26037/out}"
FSV_VASA_SAMPLE="${FSV_VASA_SAMPLE:-ZHA9292A1}"
FSV_VASA_REALCELLS="${FSV_VASA_REALCELLS:-002 003 004 005 006 007 008 009 010 011 012 013}"
FSV_VASA_BLANKS="${FSV_VASA_BLANKS:-001 014 015 016}"

# =============================================================================
# EDIT ME -- protocol switches passed to steps 5-7
# =============================================================================
FSV_STRANDED="${FSV_STRANDED:-n}"                 # note 4 above
FSV_PROTOCOL="${FSV_PROTOCOL:-smartseq_noUMI}"    # note 3
FSV_CELLID_FROM="${FSV_CELLID_FROM:-f}"           # note 6
FSV_FILT_UNIGENES="${FSV_FILT_UNIGENES:-n}"       # note 5

# =============================================================================
# EDIT ME -- references. Identical to the VASA side; that is the whole point.
# =============================================================================
FSV_VASA_REF="${FSV_VASA_REF:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116}"

# The v2 annotation BED, byte-identical to the one the VASA tables were built
# against. v1 and v2 tables cannot be mixed (v1 had no tRNA at all and its
# coordinates were 1 bp high), so this must not be changed independently of the
# VASA side.
FSV_REF_BED="${FSV_REF_BED:-$FSV_VASA_REF/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed}"

# STAR index. THIS IS THE ONE COMPONENT THAT DIFFERS BETWEEN THE TWO SIDES, and
# it is recorded in provenance.tsv rather than hidden: the VASA run used the same
# genome.fa and the same Ensembl 116 GTF, with a byte-identical junction set
# (562,855 junctions), but sjdbOverhang 150 vs the shorter-read index. 151 nt
# reads need >= 150, so this index is the correct one for FLASH-seq; an overhang
# larger than readLength-1 is harmless, one too small silently costs
# junction-spanning sensitivity. The vasalen arm's reads are SHORTER than 151,
# so 150 is comfortably large enough for it too -- both arms use this index, so
# the index is not a difference BETWEEN the arms.
FSV_STAR_INDEX="${FSV_STAR_INDEX:-/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/star_index_151_r116}"

# =============================================================================
# EDIT ME -- how much machine to use
# =============================================================================
FSV_STAR_THREADS="${FSV_STAR_THREADS:-16}"   # threads for one STAR run
FSV_NCORES="${FSV_NCORES:-10}"               # parallel workers over libraries
# Step 6's own pool width is hardcoded at ncores=8 inside
# countTables_2pickle_cellsSpliced.py and is NOT settable from here.

# =============================================================================
# You should not normally need to change anything below this line.
# =============================================================================

# --- upstream, untouched -----------------------------------------------------
# a_Mapping/ is published code: it gets explanatory comments, never logic
# changes. Steps 6 and 7 are called from here exactly as they ship.
FSV_VASA_SCRIPTS="${FSV_VASA_SCRIPTS:-$FSV_ROOT/code/I_Gene_expression/a_Mapping}"

# --- the step-5 forks, SHARED WITH THE VASA SIDE ----------------------------
# own_version/deal_with_{single,multi}mappers.sh, i.e. the SAME two scripts that
# produced the VASA tables. Using upstream's here instead would silently drop
# NH:i:10..19 on one side of the comparison only (upstream matches NH as text
# with /NH:i:1\tHI:i:1\t/ and /NH:i:[2-9]/, which between them miss 10-19
# entirely -- 894,536 reads, 4.8% of every multimapping read, measured on the
# VASA library). Both scripts are taken, not just the single one: step 6's
# per-cell glob is cell + '*_genes.bed.gz', which matches the multimapper BED
# too, so the multimapper script is part of the measurement.
#
# NB the A9 dry run used own_version for singles and UPSTREAM for multis, and
# its comment claimed that was "exactly as own_version/pipeline.sh does". That
# claim was wrong -- pipeline.sh points ASSIGN_MULTI_SH at own_version as well.
# The ten-library run uses own_version for BOTH, so that the two sides of the
# comparison share one multimapper rule. Point these at $FSV_VASA_SCRIPTS to
# get upstream behaviour.
FSV_VASA_OWN="${FSV_VASA_OWN:-$FSV_ROOT/code/I_Gene_expression/own_version}"
FSV_ASSIGN_SINGLE_SH="${FSV_ASSIGN_SINGLE_SH:-$FSV_VASA_OWN/deal_with_singlemappers.sh}"
FSV_ASSIGN_MULTI_SH="${FSV_ASSIGN_MULTI_SH:-$FSV_VASA_OWN/deal_with_multimappers.sh}"

# --- our own scripts ---------------------------------------------------------
FSV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FSV_TRIM_VASALEN="${FSV_TRIM_VASALEN:-$FSV_DIR/trim_to_vasalen.sh}"
FSV_MEASURE_LEN="${FSV_MEASURE_LEN:-$FSV_DIR/measure_vasa_readlen.py}"
FSV_BUILD_TABLES="${FSV_BUILD_TABLES:-$FSV_DIR/build_analysis_tables.py}"
FSV_RECON="${FSV_RECON:-$FSV_DIR/reconcile.py}"
# Combines the per-library step-6 dicts into the single frame step 7 expects.
# Does only what upstream step 6's parent tail does; see pipeline_fs.sh's
# `pickle` stage header for why the split is exact.
FSV_MERGE_PICKLES="${FSV_MERGE_PICKLES:-$FSV_DIR/merge_pickles.py}"

# --- cluster tools (EasyBuild module tree) -----------------------------------
# Same module builds and same absolute-path convention as own_version/config.sh
# and code/flashseq/config.sh. Trim_Galore is deliberately NOT loaded anywhere:
# it is a foss-2018b build and dragging its libstdc++ in alongside
# STAR/SAMtools/BEDTools breaks them with "GLIBCXX not found". That is why the
# adapter trim is cutadapt-direct.
FSV_EBROOT="${FSV_EBROOT:-/camp/apps/eb/software}"
FSV_P2STAR="${FSV_P2STAR:-$FSV_EBROOT/STAR/2.7.7a-GCC-10.2.0/bin}"
FSV_P2SAMTOOLS="${FSV_P2SAMTOOLS:-$FSV_EBROOT/SAMtools/1.11-GCC-10.2.0/bin}"
FSV_P2BEDTOOLS="${FSV_P2BEDTOOLS:-$FSV_EBROOT/BEDTools/2.30.0-GCC-11.2.0/bin}"

FSV_ML_INIT="source /usr/share/lmod/lmod/init/bash; export MODULEPATH=${FSV_EBROOT%/software}/modules/all"
FSV_ML_STAR="${FSV_ML_INIT}; module load STAR/2.7.7a-GCC-10.2.0 SAMtools/1.11-GCC-10.2.0"
FSV_ML_BED="${FSV_ML_INIT}; module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0"

# --- python ------------------------------------------------------------------
# envs/vasa, the env the VASA tables were built in. PANDAS MUST STAY 2.x:
# countTables_fromPickle.py calls DataFrame.applymap, deprecated in pandas 2.1
# and removed in 3.0. Never install into this env without --freeze-installed.
FSV_CONDA_ENV="${FSV_CONDA_ENV:-/nemo/lab/turnerj/working/guangxin/envs/vasa}"
FSV_CONDA_ACTIVATE="source ${FSV_EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate ${FSV_CONDA_ENV}"
FSV_PYTHON="${FSV_PYTHON:-$FSV_CONDA_ENV/bin/python3}"
# cutadapt 5.1, pip-installed into envs/vasa (the module tree stops at 1.18).
# Called by absolute path so it never has to win a PATH fight with a module.
FSV_CUTADAPT="${FSV_CUTADAPT:-$FSV_CONDA_ENV/bin/cutadapt}"

# The saved TrimGalore reports the adapter-rate check compares against.
FSV_TG_REPORTS="${FSV_TG_REPORTS:-$FSV_ROOT/data/flashseq/results/trimgalore}"
