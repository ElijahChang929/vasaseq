#!/bin/bash
#SBATCH -J vasa_rrnaref
#SBATCH -p ncpu
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/build/rrna_v2.%j.out
#SBATCH -e /nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/build/rrna_v2.%j.err
###############################################################################
# build_rrna_reference.sh -- mouse rRNA reference for the in-silico ribosomal
# depletion stage (step3 / ribo-bwamem.sh).
#
# THIS SCRIPT LIVES IN THE GIT REPO ON PURPOSE. The previous reference-build
# script lived next to its outputs in the (untracked) reference tree and was
# deleted once the outputs existed, which left three reference files on disk
# with no record of how they were made. Keep build scripts in git; the outputs
# are reproducible, the provenance is not.
#
# ---------------------------------------------------------------------------
# WHAT THE PAPER'S AUTHORS USED, AND WHY v1 WAS NOT THAT
# ---------------------------------------------------------------------------
# Salmen & De Jonghe et al. 2022, Methods ("Mapping data"):
#
#     "In silico ribosomal depletion was performed by mapping the trimmed reads
#      to mouse or human rRNA (National Center for Biotechnology Information)
#      using bwa mem and bwa aln"
#
# and a_Mapping/README.md names the actual sequences:
#
#     mouse: Rn45s, Rn6s, 12s, 16s, 47s
#     human: 12S, 16S, 45SN1-5, 45S9, 45S1-17
#
# So the authors used a HANDFUL OF FULL-LENGTH NCBI PRE-rRNA TRANSCRIPTS
# (45S/47S) plus the two mitochondrial rRNAs. They did NOT use Ensembl.
#
# v1 of this reference (unique_rRNA_mouse.fa, built 2026-07-23) was made by
# grepping the Ensembl 116 GTF for gene_biotype in {rRNA, Mt_rRNA} and pulling
# those sequences out of the genome. That yields 356 sequences whose median
# length is 116 nt, because the rDNA repeat array is collapsed/absent in the
# GRCm39 primary assembly -- Ensembl can only annotate the dispersed 5S/5.8S
# copies and fragments. Ensembl 116 mouse has no Rn28s, no Rn45s and no Rn5-8s
# gene at all; the single 18S entry (Rn18s-rs5, ENSMUSG00000119584) is a
# dispersed related-sequence copy, not the rDNA locus.
#
# Measured consequence (130 nt simulated reads tiled across the true subunits,
# run through the real ribo-bwamem.sh + riboread-selection.py):
#
#     subunit        caught by v1
#     18S              27 / 27      ok  (Rn18s-rs5 is close enough)
#     5.8S              1 / 1       ok
#     28S               0 / 71      ALL LEAK
#     5'ETS + ITS1      1 / 60      ALL LEAK
#
# 28S is 4,730 of the 13,400 nt of the transcript and is the most abundant rRNA
# species by mass. That is why the VASA-plate run only called 0.3% of reads
# ribosomal (490,012 / 181,287,059).
#
# ---------------------------------------------------------------------------
# WHAT v2 IS
# ---------------------------------------------------------------------------
# v2 = v1 (Ensembl component, kept) + the NCBI curated pre-rRNA unit (added).
#
#   component A  Ensembl 116 GRCm39, gene_biotype in {rRNA, Mt_rRNA,
#                rRNA_pseudogene}          -> 356 seqs (354 + 2 + 0)
#                Keeps the dispersed 5S/5.8S/18S copies and, importantly, the
#                two mitochondrial rRNAs (mt-Rnr1 12S, mt-Rnr2 16S) at correct
#                GRCm39 coordinates. This is the "12s, 16s" of the paper's list.
#
#   component B  BK000964.3:1-13403 -- the TRANSCRIBED portion of the mouse rDNA
#                repeating unit. This is the "Rn45s / 47s" of the paper's list.
#                Structure (GenBank feature table, 1-based inclusive):
#                    1     - 4007    5' ETS
#                    4008  - 5877    18S     (1,870 nt)
#                    5878  - 6877    ITS1
#                    6878  - 7034    5.8S    (157 nt)
#                    7035  - 8122    ITS2
#                    8123  - 12852   28S     (4,730 nt)
#                    12852 - 13403   3' ETS
#
# Two deliberate choices, both of which are the difference between this being
# right and being over-aggressive:
#
#   (1) ONE SEQUENCE, NOT SEVEN. The 13.4 kb unit goes in whole rather than
#       split into named subunits. Reads straddling a subunit boundary (18S/ITS1
#       etc.) still align cleanly, which they would not against split
#       references. VASA is a total-RNA protocol and genuinely does sequence
#       pre-rRNA across those junctions. QC decomposition is recovered by
#       POSITION instead -- see rrna_intervals.tsv, written below.
#
#   (2) IGS IS EXCLUDED. BK000964.3 is 45,306 bp; 13,404-45,305 is intergenic
#       spacer, never transcribed, and dense with SINE/LINE repeats. Including
#       it would let any repeat-containing mRNA be called ribosomal and silently
#       deleted -- riboread-selection.py ignores MAPQ, so a single low-quality
#       hit is fatal. We take only the transcript.
#
# NR_046233.2 (RefSeq "Rn45s", 13,400 nt) is the same molecule: aligned against
# BK000964.3:1-13403 it is NM=17 over 13.4 kb (99.87% identical). BK000964.3 is
# used because its subunit coordinates are authoritatively annotated, which is
# what makes choice (1) workable.
#
# ---------------------------------------------------------------------------
# OUTPUTS  (in $OUT)
# ---------------------------------------------------------------------------
#   unique_rRNA_mouse.v2.fa (+ .amb .ann .bwt .pac .sa)   <- point RRNA_FASTA here
#   rrna_intervals.tsv                                    subunit coords for QC
#   rrna_v2.provenance.txt                                accessions, dates, counts
#
# v1 (unique_rRNA_mouse.fa) is left untouched so the two can be compared.
#
# Idempotent: every step skips if its output already exists. Delete the output
# to force a rebuild.
#
# Usage:  ./build_rrna_reference.sh          (runs fine on a login node, ~2 min)
#     or  sbatch build_rrna_reference.sh
###############################################################################
set -euo pipefail

OUT=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116
B=${OUT}/build
EBROOT=/camp/apps/eb/software

GENOME=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/genome.fa
GTF=/nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/annotation/release-116/gtf/Mus_musculus.GRCm39.116.gtf

# Pinned. Bump the version suffix deliberately, never silently.
RDNA_ACC=BK000964.3
RDNA_TX_START=1
RDNA_TX_END=13403
RDNA_EXPECTED_LEN=45306        # full record length, verified after download

FA_OUT=${OUT}/unique_rRNA_mouse.v2.fa

mkdir -p "$B"
cd "$B"

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=${EBROOT%/software}/modules/all

# Modules are loaded PER STEP, never all at once. BWA/0.7.17 is a GCC-10.3.0
# build and BEDTools/2.30.0 is GCC-11.2.0; loading them together puts the older
# libstdc++ first on the library path and bedtools dies with
#   "GLIBCXX_3.4.29 not found (required by bedtools)".
# SAMtools/1.11 (GCC-10.2.0) + BEDTools coexist fine, so steps 1-2 share a load
# and step 3 purges before pulling in BWA. Same discipline as ml_trim/ml_ribo/
# ml_gmap/ml_b2b in submit_vasaplate_map.sh -- see the module gotchas in
# CLAUDE.md.
module purge 2>/dev/null || true
module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0

for f in "$GENOME" "$GTF"; do
    [ -s "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
done
[ -s "${GENOME}.fai" ] || samtools faidx "$GENOME"

###############################################################################
echo "[$(date)] === STEP 1: component A -- Ensembl 116 rRNA/Mt_rRNA genes ==="
###############################################################################
if [ ! -s ensembl_rRNA.fa ]; then
    # gene-level features whose biotype is in the rRNA family -> BED (0-based)
    awk 'BEGIN{FS="\t"} !/^#/ && $3=="gene" {
            bt=""; if (match($9, /gene_biotype "([^"]*)"/, m)) bt=m[1]
            if (bt=="rRNA" || bt=="Mt_rRNA" || bt=="rRNA_pseudogene") {
                gid="NA"; if (match($9, /gene_id "([^"]*)"/, g))   gid=g[1]
                gn="";    if (match($9, /gene_name "([^"]*)"/, n)) gn="_"n[1]
                print $1"\t"($4-1)"\t"$5"\t"gid gn"_"bt"\t.\t"$7
            }
         }' "$GTF" > ensembl_rRNA.bed
    # -s: respect strand, so minus-strand genes are stored as the sense sequence.
    # This matters because riboread-selection.py only counts FORWARD-strand hits
    # when stranded=y; a reverse-complemented reference entry would be invisible.
    bedtools getfasta -s -nameOnly -fi "$GENOME" -bed ensembl_rRNA.bed > ensembl_rRNA.fa
fi
n_ens=$(grep -c '^>' ensembl_rRNA.fa)
echo "  ensembl component: ${n_ens} seqs"

###############################################################################
echo "[$(date)] === STEP 2: component B -- NCBI ${RDNA_ACC} transcribed unit ==="
###############################################################################
if [ ! -s "${RDNA_ACC}.fa" ]; then
    wget -q -O "${RDNA_ACC}.fa.tmp" \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${RDNA_ACC}&rettype=fasta&retmode=text"
    mv "${RDNA_ACC}.fa.tmp" "${RDNA_ACC}.fa"
fi
# Guard against a truncated download or a silently re-versioned accession.
got_len=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${RDNA_ACC}.fa")
if [ "$got_len" -ne "$RDNA_EXPECTED_LEN" ]; then
    echo "  FATAL: ${RDNA_ACC} is ${got_len} bp, expected ${RDNA_EXPECTED_LEN}" >&2
    echo "  The accession changed or the download truncated. Re-check the" >&2
    echo "  subunit coordinates in the header of this script before proceeding." >&2
    exit 1
fi
echo "  ${RDNA_ACC}: ${got_len} bp (as expected)"

if [ ! -s rdna_47S.fa ]; then
    samtools faidx "${RDNA_ACC}.fa"
    samtools faidx "${RDNA_ACC}.fa" "${RDNA_ACC}:${RDNA_TX_START}-${RDNA_TX_END}" \
      | sed "1s|.*|>mouse_rDNA_47S_${RDNA_ACC}_${RDNA_TX_START}-${RDNA_TX_END}|" > rdna_47S.fa
fi
echo "  47S transcript: $(awk '/^>/{next}{n+=length($0)}END{print n}' rdna_47S.fa) bp"

###############################################################################
echo "[$(date)] === STEP 3: concatenate + bwa index ==="
###############################################################################
# See the module comment at the top: BWA cannot be co-loaded with BEDTools.
module purge 2>/dev/null || true
module load BWA/0.7.17-GCC-10.3.0

if [ ! -s "$FA_OUT" ]; then
    cat rdna_47S.fa ensembl_rRNA.fa > "${FA_OUT}.tmp"
    mv "${FA_OUT}.tmp" "$FA_OUT"
fi
[ -s "${FA_OUT}.bwt" ] || bwa index "$FA_OUT"
echo "  ${FA_OUT}: $(grep -c '^>' "$FA_OUT") seqs"

###############################################################################
echo "[$(date)] === STEP 4: QC interval table ==="
###############################################################################
# Lets you decompose the .Ribo.bam by subunit without splitting the reference:
#   samtools view cell.Ribo.bam | awk '$3 ~ /rDNA_47S/ {print $4}' \
#     | ... bin against these intervals
cat > "${OUT}/rrna_intervals.tsv" <<EOF
# subunit intervals within mouse_rDNA_47S_${RDNA_ACC}_${RDNA_TX_START}-${RDNA_TX_END}
# 1-based inclusive, from the ${RDNA_ACC} GenBank feature table
#subunit	start	end	length
5ETS	1	4007	4007
18S	4008	5877	1870
ITS1	5878	6877	1000
5.8S	6878	7034	157
ITS2	7035	8122	1088
28S	8123	12852	4730
3ETS	12853	13403	551
EOF
echo "  wrote ${OUT}/rrna_intervals.tsv"

###############################################################################
echo "[$(date)] === STEP 5: provenance ==="
###############################################################################
cat > "${OUT}/rrna_v2.provenance.txt" <<EOF
unique_rRNA_mouse.v2.fa
built    : $(date -Iseconds)
by       : $(basename "$0") (tracked in code/I_Gene_expression/own_version/)
host     : $(hostname)

component A -- Ensembl
  genome : ${GENOME}
  gtf    : ${GTF}
  filter : gene feature, gene_biotype in {rRNA, Mt_rRNA, rRNA_pseudogene}
  seqs   : ${n_ens}

component B -- NCBI
  accession : ${RDNA_ACC} (${got_len} bp, full repeating unit)
  extracted : ${RDNA_TX_START}-${RDNA_TX_END} (transcribed 47S; IGS excluded)
  seqs      : 1

total seqs  : $(grep -c '^>' "$FA_OUT")
total bp    : $(awk '/^>/{next}{n+=length($0)}END{print n}' "$FA_OUT")

consumed by : ribo-bwamem.sh (bwa aln + bwa mem), via RRNA_FASTA in config.sh
EOF
echo "  wrote ${OUT}/rrna_v2.provenance.txt"

echo "[$(date)] === DONE ==="
echo "RRNA_FASTA=${FA_OUT}"
