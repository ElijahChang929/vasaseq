#!/bin/bash
################################################################################
# run_mapping_stepwise.sh
#
# WHAT THIS IS
#   A "flattened" / unrolled version of submit_vasaplate_map.sh, written so you
#   can SEE and RUN each real command by hand, one step at a time, and inspect
#   the output of every stage before moving on.
#
#   The original submit_vasaplate_map.sh hides the real work: it only submits
#   sbatch jobs that call helper scripts (extractBC.sh, trim.sh, ribo-bwamem.sh,
#   map_star.sh, deal_with_*.sh, countTables_*.py). This file INLINES the actual
#   commands those helper scripts run (bwa, STAR, cutadapt, samtools, bedtools,
#   ...), so nothing is hidden and every parameter is editable in one place.
#
# HOW TO USE
#   1. Edit the CONFIG block below (paths, sample name, species, read length).
#   2. Run the steps ONE AT A TIME. The recommended way is to open this file and
#      copy/paste each STEP into your terminal, OR set RUN_STEP at the top and
#      run `bash run_mapping_stepwise.sh` once per step.
#   3. After each step, look at the "OUTPUTS" listed in its header to check the
#      result before continuing.
#
#   This is deliberately NOT a one-shot pipeline. It processes ONE fastq file
#   (one sample / one lane). If you have several lanes or samples, run it once
#   per file, then do STEP 6/7 (counting) over the whole folder at the end.
#
# NOTE ON SLURM
#   None of the steps below submit sbatch jobs; they run directly. If a step is
#   heavy (STAR, bwa) and you are on a login node, wrap the single command in
#   `srun ... ` or `sbatch --wrap="..."` yourself. The resource hints from the
#   original script are noted in each step header (time / mem / cpus).
################################################################################


############################## CONFIG (EDIT ME) ################################
# ---- your sample -----------------------------------------------------------
LIB=mysample                 # base name of your fastq (before _R1/_R2). EDIT.
FOLDER=/path/to/output       # where all intermediate + output files go. EDIT.
PROTOCOL=vasaplate           # keep vasaplate for VASA-plate. Others: celseq1/celseq2/vasadrop/10x
STRANDED=y                   # VASA is stranded -> y
CELLIDORI=f                  # cell ID from filename (f) or read name (r); used only in STEP 6

# ---- species / references (THIS is what you change per species) ------------
# Only these three variables define the species. Point them at YOUR references.
SPECIES=MOUSE                # label only, for your own bookkeeping
READLEN=74                   # read length; must match the STAR index you built (see GENOME)

# rRNA fasta (used to remove ribosomal reads in STEP 3)
RIBOREF=/path/to/refs/rRNA_mouse_Rn45S.fa
# STAR genome index DIRECTORY. NB: the original builds a separate index per read
# length and appends the length, e.g. .../star_..._index_74 . Your index dir must
# match READLEN, or just give the exact directory here.
GENOME=/path/to/refs/star_index_${READLEN}
# intron/exon/tRNA annotation BED (used to assign reads to genes in STEP 5)
REFBED=/path/to/refs/annotation.homemade_IntronExonTrna.bed

# ---- tool locations (edit to your HPC / modules) ---------------------------
P2S=$(cd "$(dirname "$0")" && pwd)     # dir holding the helper scripts (this dir)
P2TRIMGALORE=/path/to/TrimGalore-0.6.6
P2CUTADAPT=/path/to/cutadapt/bin       # dir containing the `cutadapt` executable
P2BWA=/path/to/bwa-0.7.17              # dir containing the `bwa` executable
P2SAMTOOLS=/path/to/samtools-1.11      # dir containing the `samtools` executable
P2STAR=/path/to/STAR/bin/Linux_x86_64  # dir containing the `STAR` executable
P2BEDTOOLS=/path/to/bedtools2/bin      # dir containing the `bedtools` executable

# ---- which step to run (used only if you run this file as a script) --------
# Set to 0 to run nothing (copy/paste steps manually instead), or 1..7.
RUN_STEP=0
###############################################################################


mkdir -p "$FOLDER" 2>/dev/null || true

# Convenience name used from STEP 2 onward: the per-cell barcoded fastq that
# STEP 1 produces. After STEP 1 you will have one such file PER CELL/LANE named
# <something>_cbc.fastq.gz inside $FOLDER. Set CBC to the specific one you want
# to process through steps 2-5. (Do steps 2-5 once per _cbc.fastq.gz file.)
CBC=${FOLDER}/${LIB}_cbc.fastq.gz      # EDIT to the exact _cbc file you are mapping
BASE=${CBC%.fastq.gz}                  # e.g. .../<lib>_cbc   -> used to build all later names


###############################################################################
# STEP 1 - EXTRACT CELL BARCODES  (was: extractBC.sh -> concatenator.py)
#   Reads raw R1/R2, moves the cell-barcode + UMI onto the read name, and
#   demultiplexes into one fastq per cell.
#   For vasaplate: CELseq2 barcodes, UMI length 6, UMI first, demultiplex.
#   SLURM hint in original: 1 cpu, 10G, up to 48h.
#
#   INPUTS : ${LIB}_R1*.fastq.gz and ${LIB}_R2*.fastq.gz  (in current dir)
#   OUTPUTS: ${FOLDER}/*_cbc.fastq.gz   (one per cell)  + a demux log
#   TWEAK  : --lenumi (UMI length), --cbcfile (barcode list), --umifirst, --demux
#
# -----------------------------------------------------------------------------
# HOW concatenator.py WORKS (the "concatenator" = it CONCATENATES the barcode
# read onto the biological read's name, then demultiplexes):
#
#   1. FIND FILES. --fqf is a PREFIX, not a filename. The script globs
#      <prefix>*_R1*.fastq.gz and <prefix>*_R2*.fastq.gz and pairs them up
#      (so multiple lanes are handled together). R1/R2 must stay in sync.
#
#   2. READ LAYOUT. By default the barcode lives in R1 (--bcread R1) and the
#      biological cDNA sequence in R2 (--bioread R2). Only the FIRST
#      (lenumi + lencbc) bases of the barcode read are the barcode block; the
#      rest of that read is discarded. For VASA-plate: lencbc=8 (default),
#      lenumi=6, and --umifirst means the block is ordered [UMI | cell-barcode]:
#          bases 0..lenumi-1     -> UMI          (random molecular tag)
#          bases lenumi..end     -> cell barcode (identifies the well/cell)
#      Without --umifirst the order is [cell-barcode | UMI] instead.
#
#   3. BARCODE WHITELIST + ERROR TOLERANCE. --cbcfile (here bc_celseq2.tsv) is
#      the list of the 384 valid CEL-seq2 cell barcodes -> cell IDs. With
#      --cbchd 0 only exact matches (plus a single-N substitution) are allowed;
#      a higher Hamming distance would also accept barcodes within that many
#      mismatches, but ONLY if the corrected barcode maps to exactly one cell
#      (ambiguous ones are dropped).
#
#   4. PER-READ ASSIGNMENT. For every read the cell barcode is looked up in the
#      whitelist. If it matches a cell, the read is KEPT; if not, it is skipped.
#      The barcode/UMI info is written into the read NAME as ';'-separated tags:
#          SS = observed cell-barcode sequence
#          CB = matched/corrected whitelist barcode
#          QT = cell-barcode quality
#          RX = UMI sequence        <- used later for UMI dedup in STEP 6
#          RQ = UMI quality
#          SM = cell ID (zero-padded, e.g. 001..384)
#      The OUTPUT read is the biological (R2) sequence+quality under this new
#      name. This is why later steps can recover cell + UMI straight from the
#      read name -- the barcode read itself is thrown away after this step.
#      (Barcode/UMI phred chars are shifted +32 so they can't be confused with
#      the biological read's quality string downstream.)
#
#   5. OUTPUT MODE. With --demux (used for VASA-plate) it writes ONE fastq per
#      cell: <prefix>_<cellID>_cbc.fastq.gz. Without --demux it writes a single
#      combined <prefix>_cbc.fastq.gz. A <prefix>.log records total reads and
#      the fraction carrying a valid barcode (a key QC number).
#
#   => To change chemistry you edit the barcode geometry here: --lencbc,
#      --lenumi, --umifirst, --cbcfile, and --cbchd. Nothing else in STEP 1
#      needs touching.
# -----------------------------------------------------------------------------
###############################################################################
step1_extract () {
  python3 ${P2S}/concatenator.py \
      --fqf ${LIB} \
      --cbcfile ${P2S}/bc_celseq2.tsv \
      --cbchd 0 \
      --lenumi 6 \
      --umifirst \
      --demux \
      --outdir ${FOLDER}
  echo ">> STEP 1 done. Barcoded fastqs:"; ls -1 ${FOLDER}/*_cbc.fastq.gz
}


###############################################################################
# STEP 2 - TRIM  (was: trim.sh -> trim_galore + cutadapt)
#   (a) TrimGalore removes sequencing adapters.
#   (b) cutadapt removes homopolymer stretches (polyG/C/T/A) and short reads.
#   SLURM hint: 1 node, 10G, ~5h.
#
#   INPUTS : ${CBC}                         (a *_cbc.fastq.gz from STEP 1)
#   OUTPUTS: ${BASE}_trimmed.fq.gz          (adapter-trimmed)
#            ${BASE}_trimmed_homoATCG.fq.gz (adapter + homopolymer trimmed)  <-- used next
#   TWEAK  : cutadapt -m 15 (min length), the -a polyX patterns
###############################################################################
step2_trim () {
  # (a) adapter trimming
  ${P2TRIMGALORE}/trim_galore --path_to_cutadapt ${P2CUTADAPT}/cutadapt ${CBC} -o ${FOLDER}
  mv ${CBC}_trimming_report.txt ${BASE}_trimming_report.txt 2>/dev/null

  # (b) homopolymer trimming
  ${P2CUTADAPT}/cutadapt -m 15 --trim-n \
      -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" \
      -o ${BASE}_trimmed_homoATCG.fq.gz \
         ${BASE}_trimmed.fq.gz
  echo ">> STEP 2 done -> ${BASE}_trimmed_homoATCG.fq.gz"
}


###############################################################################
# STEP 3 - REMOVE RIBOSOMAL READS  (was: ribo-bwamem.sh -> bwa + riboread-selection.py)
#   Map reads to the rRNA reference with BWA (two strategies: aln and mem),
#   merge, sort by name, then riboread-selection.py splits reads into:
#     - ribosomal  -> kept in a .Ribo.bam
#     - non-ribosomal -> written to a new fastq for genome mapping
#   SLURM hint: 8 cpus, 40G, ~10h.
#
#   INPUTS : ${BASE}_trimmed_homoATCG.fq.gz
#            ${RIBOREF}   <-- SPECIES-SPECIFIC
#   OUTPUTS: ${BASE}_trimmed_homoATCG.nonRibo.fastq.gz  <-- used next
#            ${BASE}_trimmed_homoATCG.Ribo.bam
#            ${BASE}_trimmed_homoATCG.ribo-map.log       (read counts)
#   TWEAK  : STRANDED (y/n), bwa mem -h 15
###############################################################################
step3_ribo () {
  local FQ=${BASE}_trimmed_homoATCG.fq.gz
  local OUT=${BASE}_trimmed_homoATCG
  local SAI=${FOLDER}/aln_$(basename ${FQ%.f*q}).sai

  # map with bwa aln (short-read) and bwa mem (in parallel), both -> bam
  ${P2BWA}/bwa aln ${RIBOREF} ${FQ} > ${SAI}
  ${P2BWA}/bwa samse ${RIBOREF} ${SAI} ${FQ} | ${P2SAMTOOLS}/samtools view -Sb > ${OUT}.aln-ribo.bam &
  ${P2BWA}/bwa mem -t 8 -h 15 ${RIBOREF} ${FQ} | ${P2SAMTOOLS}/samtools view -Sb > ${OUT}.mem-ribo.bam &
  wait

  # merge, name-sort, then classify ribo vs non-ribo
  ${P2SAMTOOLS}/samtools merge -n -r -h ${OUT}.aln-ribo.bam --threads 8 ${OUT}.all-ribo.bam ${OUT}.aln-ribo.bam ${OUT}.mem-ribo.bam
  rm ${OUT}.aln-ribo.bam ${OUT}.mem-ribo.bam ${SAI}
  ${P2SAMTOOLS}/samtools sort -n --threads 8 ${OUT}.all-ribo.bam -O BAM -o ${OUT}.nsorted.all-ribo.bam
  rm ${OUT}.all-ribo.bam
  ${P2S}/riboread-selection.py ${OUT}.nsorted.all-ribo.bam ${STRANDED} ${OUT}
  echo ">> STEP 3 done -> ${OUT}.nonRibo.fastq.gz  (see ${OUT}.ribo-map.log)"
}


###############################################################################
# STEP 4 - MAP TO GENOME  (was: map_star.sh -> STAR)
#   Align the non-ribosomal reads to the genome with STAR.
#   Keeps unmapped reads inside the BAM; allows up to 20 multimap positions.
#   SLURM hint: 8 cpus, 40G, ~10h.
#
#   INPUTS : ${BASE}_trimmed_homoATCG.nonRibo.fastq.gz
#            ${GENOME}   <-- SPECIES-SPECIFIC (must match READLEN)
#   OUTPUTS: ${BASE}_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam  <-- used next
#            ${BASE}_trimmed_homoATCG.nonRibo_E99_Log.final.txt    (mapping stats)
#   TWEAK  : --outFilterMultimapNmax 20 , --outSAMattributes
###############################################################################
step4_genome () {
  local INFQ=${BASE}_trimmed_homoATCG.nonRibo.fastq.gz
  local OUTPREF=${BASE}_trimmed_homoATCG.nonRibo_E99_

  ${P2STAR}/STAR --runThreadN 8 --genomeDir ${GENOME} \
      --readFilesIn ${INFQ} --readFilesCommand zcat \
      --outFilterMultimapNmax 20 --outSAMunmapped Within \
      --outSAMtype BAM Unsorted --outSAMattributes All \
      --outFileNamePrefix ${OUTPREF}

  rm -r ${OUTPREF}_STARtmp 2>/dev/null
  rm ${OUTPREF}Log.progress.out 2>/dev/null
  mv ${OUTPREF}Log.out ${OUTPREF}Log.txt 2>/dev/null
  mv ${OUTPREF}Log.final.out ${OUTPREF}Log.final.txt 2>/dev/null
  echo ">> STEP 4 done -> ${OUTPREF}Aligned.out.bam  (stats in ${OUTPREF}Log.final.txt)"
}


###############################################################################
# STEP 5 - ASSIGN READS TO GENES  (was: deal_with_singlemappers.sh + _multimappers.sh)
#   For each aligned read, intersect with the intron/exon/tRNA BED to decide
#   which gene/feature it belongs to and whether it is intronic/exonic (spliced
#   vs unspliced). Uniquely-mapping (NH:i:1) and multi-mapping (NH:i:>=2) reads
#   are handled by two separate scripts and produce two BED outputs.
#   The heavy awk logic is left inside the helper scripts (you would not normally
#   tweak it); each script internally runs:
#       samtools view | awk (tag reads) | bedtools bamtobed | bedtools intersect -b REFBED
#   SLURM hint: 1 cpu, 40G, ~10-12h each. These two are independent (run in parallel).
#
#   INPUTS : ${BASE}_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam
#            ${REFBED}   <-- SPECIES-SPECIFIC
#   OUTPUTS: *.singlemappers_genes.bed.gz
#            *.nsorted.multimappers_genes.bed.gz
#   TWEAK  : STRANDED (y/n)
###############################################################################
step5_assign () {
  local BAM=${BASE}_trimmed_homoATCG.nonRibo_E99_Aligned.out.bam
  ${P2S}/deal_with_singlemappers.sh ${BAM} ${REFBED} ${STRANDED} ${P2SAMTOOLS} ${P2BEDTOOLS} &
  ${P2S}/deal_with_multimappers.sh  ${BAM} ${REFBED} ${STRANDED} ${P2SAMTOOLS} ${P2BEDTOOLS} &
  wait
  echo ">> STEP 5 done -> *singlemappers_genes.bed.gz and *multimappers_genes.bed.gz"
}


###############################################################################
# STEP 6 - BUILD COUNT PICKLE  (was: countTables_2pickle_cellsSpliced.py)
#   Collapses the per-read gene assignments across ALL cells in $FOLDER into a
#   single pickled count structure, UMI-aware and spliced/unspliced-aware.
#   RUN THIS ONCE, after steps 2-5 have been done for every cell/lane fastq.
#   SLURM hint: 8 cpus, 160G, up to 60h.
#
#   INPUTS : $FOLDER (all the *_genes.bed.gz produced in STEP 5)
#   OUTPUTS: ${FOLDER}.pickle.gz
#   ARGS   : <cellFolder> <outputPrefix> <protocol> <cellid f|r>
###############################################################################
step6_pickle () {
  ${P2S}/countTables_2pickle_cellsSpliced.py ${FOLDER} ${FOLDER} vasa ${CELLIDORI}
  echo ">> STEP 6 done -> ${FOLDER}.pickle.gz"
}


###############################################################################
# STEP 7 - FINAL COUNT TABLES  (was: countTables_fromPickle.py)
#   Turns the pickle into the final, human-readable count tables
#   (genes x cells), separating spliced / unspliced etc.
#   SLURM hint: 1 cpu, 160G, ~15h.
#
#   INPUTS : ${FOLDER}.pickle.gz
#   OUTPUTS: final count-table files prefixed with ${FOLDER}
#   ARGS   : <input.pickle.gz> <outputPrefix> <protocol> <filter uniq genes y|n>
###############################################################################
step7_tables () {
  ${P2S}/countTables_fromPickle.py ${FOLDER}.pickle.gz ${FOLDER} vasa y
  echo ">> STEP 7 done -> count tables prefixed ${FOLDER}"
}


############################## dispatcher ######################################
# If you run `bash run_mapping_stepwise.sh` with RUN_STEP set to 1..7 above,
# it will run just that step. Otherwise it prints this help and does nothing,
# so you can copy/paste the step bodies by hand.
case "$RUN_STEP" in
  1) step1_extract ;;
  2) step2_trim ;;
  3) step3_ribo ;;
  4) step4_genome ;;
  5) step5_assign ;;
  6) step6_pickle ;;
  *) echo "Nothing run. Edit CONFIG, then set RUN_STEP=1..7 (or copy/paste a stepN_ function body).";;
esac
