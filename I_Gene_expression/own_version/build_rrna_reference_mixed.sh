#!/bin/bash
#SBATCH -J vasa_rrnaref_mixed
#SBATCH -p ncpu
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH -o /nemo/lab/turnerj/working/guangxin/reference/vasaseq/mixed/build/rrna_v2.%j.out
#SBATCH -e /nemo/lab/turnerj/working/guangxin/reference/vasaseq/mixed/build/rrna_v2.%j.err
###############################################################################
# build_rrna_reference_mixed.sh -- HUMAN+MOUSE rRNA reference for step 3.
#
# Deliberately a SEPARATE script from build_rrna_reference.sh (mouse-only, for
# the PM26037 data). They share a diagnosis but almost nothing else: different
# assemblies (GRCh38+GRCm38 / Ensembl 99, not GRCm39 / 116), two species, and a
# different NCBI component on the human side. Merging them would mean a species
# switch threaded through every step. Read build_rrna_reference.sh first -- its
# header carries the full argument for why the Ensembl-only reference is wrong;
# this header only covers what differs.
#
# WHICH REFERENCE THIS IS, AND WHY IT MATTERS
# -------------------------------------------
# $VASA_REFS/mixed/ is the reference for SRR14783059, the paper's own
# vasaplate-HEK293T-mESC species-mixing control -- i.e. the run we use to check
# our pipeline against the authors'. Its rRNA fasta (v1, 2026-07-17) was built
# the same Ensembl-only way as the mouse v1, and has the same hole.
#
# MEASURED on that run, aggregated over all 384 cells' *.ribo-map.log:
#
#     reads into step 3 : 181,287,059
#     called ribosomal  :     490,012   (0.27%)
#
# 0.27% for a total-RNA protocol whose entire point is that it captures rRNA is
# not a result, it is a broken reference. 28S alone is 4.7 kb of the 13.4 kb
# transcript and is the most abundant species by mass, and neither GRCh38 nor
# GRCm38 has the rDNA array in the primary assembly, so Ensembl cannot annotate
# it.
#
# WHAT THE AUTHORS ACTUALLY USED  (a_Mapping/README.md, step 3)
# -------------------------------------------------------------
#     "...a fasta file containing the reference sequences for the ribosomal RNA
#      (e.g. Rn45s, Rn6s, 12s, 16s, 47s for mouse, 12S, 16S, 45SN1-5, 45S9,
#      45S1-17 for human)"
#
# So: full-length NCBI pre-rRNA transcripts plus the two mitochondrial rRNAs,
# per species. This script reproduces that.
#
# COMPONENTS
# ----------
#   A  Ensembl 99 rRNA/Mt_rRNA/rRNA_pseudogene genes, both species (915 seqs:
#      559 human + 356 mouse). Kept from v1 -- it is what supplies the "12S,
#      16S" mitochondrial rRNAs at correct assembly coordinates, plus the
#      dispersed 5S/5.8S copies. Reused from mixed/build/{human,mouse}.rRNA.fa
#      if present, rebuilt from the Ensembl 99 GTFs otherwise.
#
#   B  human: RefSeq RNA45SN1-N5, the paper's "45SN1-5". Five near-identical
#      13.3 kb 45S pre-rRNA transcripts. VERIFIED by fetching each and checking
#      its DEFINITION line -- note NR_046235.3 is N5, not N1, which is exactly
#      the sort of thing that is wrong if you go from memory.
#
#      Unlike the mouse case there is no IGS to trim: these are RefSeq
#      TRANSCRIPTS, so the record is already the transcribed unit. That is why
#      this side uses RefSeq rather than the U13369.1 complete repeating unit --
#      no coordinates to pick, nothing to get wrong.
#
#      All five go in even though they are ~99% identical. That is what the
#      paper lists, and it costs nothing: riboread-selection.py discards a read
#      that maps ANYWHERE in this reference, multi-mappers included, so extra
#      near-duplicate copies change no decision.
#
#   B  mouse: BK000964.3:1-13403, the transcribed 47S. Identical to the mouse
#      script's component B -- see it for why the IGS (13,404-45,305) is
#      excluded and why the unit goes in whole rather than split by subunit.
#      rDNA is not in either mouse primary assembly, so GRCm38 vs GRCm39 makes
#      no difference here.
#
# OUTPUTS (in $OUT)
# -----------------
#   unique_rRNA_human_mouse.v2.fa (+ bwa index)   <- point riboref here
#   rrna_intervals_mixed.tsv                      subunit coords, both species
#   unique_rRNA_human_mouse.v3.provenance.txt     accessions, dates, counts
#                                                 (named from FA_OUT, so a v2 and
#                                                  a v3 build never clobber each other)
#
# v1 (unique_rRNA_human_mouse.fa) is left untouched so the two can be compared;
# validate_rrna_reference.sh does exactly that.
#
# Idempotent: every step skips if its output exists. Delete it to force a rebuild.
#
# Usage:  ./build_rrna_reference_mixed.sh      (login node is fine, ~3 min)
#     or  sbatch build_rrna_reference_mixed.sh
###############################################################################
set -euo pipefail

VASA_REFS=${VASA_REFS:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq}
OUT=${VASA_REFS}/mixed
B=${OUT}/build
EBROOT=/camp/apps/eb/software

HUMAN_GTF=${B}/Homo_sapiens.GRCh38.99.gtf.gz
MOUSE_GTF=${B}/Mus_musculus.GRCm38.99.gtf.gz
HUMAN_GENOME=${B}/human.genome.fa
MOUSE_GENOME=${B}/mouse.genome.fa

# Pinned. Bump deliberately, never silently. Lengths are asserted after download
# so a re-versioned or truncated accession fails loudly instead of quietly
# shrinking the reference.
MOUSE_RDNA_ACC=BK000964.3
MOUSE_RDNA_TX_START=1
MOUSE_RDNA_TX_END=13403
MOUSE_RDNA_EXPECTED_LEN=45306

# human RNA45SN1-N5, accession -> expected length. Verified 2026-07-26 by
# fetching each and reading its DEFINITION line.
HUMAN_45S_ACC=(NR_145819.1 NR_146144.1 NR_146151.1 NR_146117.1 NR_046235.3)
HUMAN_45S_NAME=(RNA45SN1   RNA45SN2   RNA45SN3   RNA45SN4   RNA45SN5)
HUMAN_45S_LEN=(13351      13315      13309      13373      13357)

# v3 = v2 + STEP 4 (orientation). v2 is left on disk untouched so the two can be
# compared and so the tables already built against it keep a valid reference.
FA_OUT=${FA_OUT:-${OUT}/unique_rRNA_human_mouse.v3.fa}

# ORIENT_TO_UNITS -- see the STEP 4 header. Set to `no` (with FA_OUT pointed at
# the v2 filename) to reproduce v2 byte-for-byte.
ORIENT_TO_UNITS=${ORIENT_TO_UNITS:-yes}

mkdir -p "$B"
cd "$B"

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=${EBROOT%/software}/modules/all

# Per-step module loads, never all at once: BWA/0.7.17 is GCC-10.3.0 and
# BEDTools/2.30.0 is GCC-11.2.0, and co-loading them puts the older libstdc++
# first so bedtools dies with "GLIBCXX_3.4.29 not found". Same discipline as
# build_rrna_reference.sh and the ml_* helpers in submit_vasaplate_map.sh.
module purge 2>/dev/null || true
module load SAMtools/1.11-GCC-10.2.0 BEDTools/2.30.0-GCC-11.2.0

###############################################################################
echo "[$(date)] === STEP 1: component A -- Ensembl 99 rRNA genes, both species ==="
###############################################################################
# Same extraction as the mouse script. -s is load-bearing: it stores minus-strand
# genes as the sense sequence, and riboread-selection.py only counts FORWARD
# hits when stranded=y, so a reverse-complemented entry would be invisible.
#
# BUT -s IS NOT SUFFICIENT, which is what STEP 4 exists to repair. It gives the
# sense strand OF THE ANNOTATED GENE, which is transcript-sense only when Ensembl
# put the gene on the same strand as the rRNA transcript. Where it did not
# (ENSMUSG00000106106 in Ensembl 99), -s faithfully produces a backwards entry.
# Do not trust this step alone to get orientation right.
extract_ensembl_rrna() {   # $1=gtf(.gz) $2=genome $3=species_tag $4=out.fa
    local gtf=$1 genome=$2 tag=$3 out=$4
    local bed="${tag}.rRNA.bed"
    [ -s "${genome}.fai" ] || samtools faidx "$genome"
    if [ ! -s "$bed" ]; then
        local cat=cat; case "$gtf" in *.gz) cat=zcat ;; esac
        $cat "$gtf" | awk -v tag="$tag" 'BEGIN{FS="\t"} !/^#/ && $3=="gene" {
                bt=""; if (match($9, /gene_biotype "([^"]*)"/, m)) bt=m[1]
                if (bt=="rRNA" || bt=="Mt_rRNA" || bt=="rRNA_pseudogene") {
                    gid="NA"; if (match($9, /gene_id "([^"]*)"/, g))   gid=g[1]
                    gn="";    if (match($9, /gene_name "([^"]*)"/, n)) gn="_"n[1]
                    print $1"\t"($4-1)"\t"$5"\t"tag"_"gid gn"_"bt"\t.\t"$7
                }
             }' > "$bed"
    fi
    [ -s "$out" ] || bedtools getfasta -s -nameOnly -fi "$genome" -bed "$bed" > "$out"
    echo "  ${tag}: $(grep -c '^>' "$out") seqs"
}

if [ -s human.rRNA.fa ] && [ -s mouse.rRNA.fa ]; then
    echo "  reusing existing human.rRNA.fa / mouse.rRNA.fa from the v1 build"
    echo "  human: $(grep -c '^>' human.rRNA.fa) seqs, mouse: $(grep -c '^>' mouse.rRNA.fa) seqs"
else
    for f in "$HUMAN_GTF" "$MOUSE_GTF" "$HUMAN_GENOME" "$MOUSE_GENOME"; do
        [ -s "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
    done
    extract_ensembl_rrna "$HUMAN_GTF" "$HUMAN_GENOME" human human.rRNA.fa
    extract_ensembl_rrna "$MOUSE_GTF" "$MOUSE_GENOME" mouse mouse.rRNA.fa
fi
n_ens=$(( $(grep -c '^>' human.rRNA.fa) + $(grep -c '^>' mouse.rRNA.fa) ))

###############################################################################
echo "[$(date)] === STEP 2: component B -- human RNA45SN1-N5 (RefSeq) ==="
###############################################################################
fetch_ncbi() {   # $1=accession $2=outfile
    [ -s "$2" ] && return 0
    wget -q -O "${2}.tmp" \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${1}&rettype=fasta&retmode=text"
    [ -s "${2}.tmp" ] || { echo "  FATAL: empty download for $1" >&2; exit 1; }
    mv "${2}.tmp" "$2"
}

: > human_45S.fa.tmp
for i in "${!HUMAN_45S_ACC[@]}"; do
    acc=${HUMAN_45S_ACC[$i]}; nm=${HUMAN_45S_NAME[$i]}; want=${HUMAN_45S_LEN[$i]}
    fetch_ncbi "$acc" "${acc}.fa"
    got=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${acc}.fa")
    if [ "$got" -ne "$want" ]; then
        echo "  FATAL: ${acc} is ${got} bp, expected ${want}" >&2
        echo "  The accession was re-versioned or the download truncated." >&2
        exit 1
    fi
    # Assert the record really is the gene we think it is. NR_046235.3 is N5,
    # not N1 -- getting this pairing wrong silently is entirely possible.
    grep -q "(${nm})" "${acc}.fa" || {
        echo "  FATAL: ${acc} DEFINITION does not mention ${nm}:" >&2
        head -1 "${acc}.fa" >&2; exit 1; }
    sed "1s|.*|>human_45S_${nm}_${acc}|" "${acc}.fa" >> human_45S.fa.tmp
    echo "  ${nm} ${acc}: ${got} bp  ok"
done
mv human_45S.fa.tmp human_45S.fa

###############################################################################
echo "[$(date)] === STEP 3: component B -- mouse ${MOUSE_RDNA_ACC} 47S unit ==="
###############################################################################
fetch_ncbi "$MOUSE_RDNA_ACC" "${MOUSE_RDNA_ACC}.fa"
got=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${MOUSE_RDNA_ACC}.fa")
if [ "$got" -ne "$MOUSE_RDNA_EXPECTED_LEN" ]; then
    echo "  FATAL: ${MOUSE_RDNA_ACC} is ${got} bp, expected ${MOUSE_RDNA_EXPECTED_LEN}" >&2
    echo "  Re-check the subunit coordinates in this script before proceeding." >&2
    exit 1
fi
echo "  ${MOUSE_RDNA_ACC}: ${got} bp (as expected)"
if [ ! -s mouse_rdna_47S.fa ]; then
    samtools faidx "${MOUSE_RDNA_ACC}.fa"
    samtools faidx "${MOUSE_RDNA_ACC}.fa" \
        "${MOUSE_RDNA_ACC}:${MOUSE_RDNA_TX_START}-${MOUSE_RDNA_TX_END}" \
      | sed "1s|.*|>mouse_rDNA_47S_${MOUSE_RDNA_ACC}_${MOUSE_RDNA_TX_START}-${MOUSE_RDNA_TX_END}|" \
      > mouse_rdna_47S.fa
fi
echo "  47S transcript: $(awk '/^>/{next}{n+=length($0)}END{print n}' mouse_rdna_47S.fa) bp"

###############################################################################
echo "[$(date)] === STEP 4: orient every Ensembl entry to the NCBI units ==="
###############################################################################
# WHY THIS STEP EXISTS
#
# riboread-selection.py decides what is ribosomal by STRAND:
#
#     mapreads = [x for x in reads if (not x.is_unmapped) and (not x.is_reverse)]
#
# Only FORWARD hits count; a read whose hits are all reverse is kept, on the
# reasoning that VASA is stranded so a genuine rRNA read must be sense. That is
# correct -- but it silently assumes EVERY SEQUENCE IN THIS FASTA IS STORED IN
# TRANSCRIPT-SENSE ORIENTATION. Orientation is the decision rule, so a single
# backwards entry does not merely fail to help, it actively shields reads from
# depletion.
#
# Component A cannot guarantee that on its own. `bedtools getfasta -s` yields the
# sense strand OF THE ANNOTATED GENE, which equals transcript sense only when
# Ensembl annotates the gene on the same strand as the rRNA transcript. Measured
# on the v2 build: 914 of 915 entries were fine and exactly ONE was not --
# ENSMUSG00000106106 (CT010467.1), annotated opposite the transcript, so its
# entry was the 18S region stored backwards.
#
# What that one entry cost, measured on the 384-cell control:
#   - it and the 47S unit capture the SAME reads (30/30 on the test set), but it
#     scores HIGHER (mean AS 67.3 vs 64.7) because it is the exact GRCm38 locus
#     the reads came from while BK000964.3 is a curated consensus;
#   - bwa mem reports one primary alignment, so the backwards entry wins it and
#     returns REVERSE;
#   - the read is therefore spared, then written reverse-complemented by
#     riboread-selection.py, then mapped by STAR sense to the minus-strand gene;
#   - result: 95,823 spurious UFIs on that one locus, against 72 in the
#     published table, and step-3 depletion understated by ~0.28 points
#     (5.06% -> 5.34%).
#
# So this step fixes the CLASS, not the instance: align every Ensembl entry to
# the NCBI full-length units and reverse-complement any that come back flag 16.
# Entries with no homology to the units (5S, dispersed fragments) are left alone
# -- there is nothing to orient them against, and they are not 45S-derived, so
# they cannot be captured by the units either way.
module purge 2>/dev/null || true
module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0

ENS_ALL=ensembl_all.fa
cat human.rRNA.fa mouse.rRNA.fa > "$ENS_ALL"

if [ "$ORIENT_TO_UNITS" = "yes" ]; then
    cat human_45S.fa mouse_rdna_47S.fa > units.fa
    [ -s units.fa.bwt ] || bwa index units.fa 2>/dev/null
    bwa mem -t 4 units.fa "$ENS_ALL" 2>/dev/null | samtools view - \
      | awk '$2==16 {print $1}' | sort -u > antisense_names.txt
    n_anti=$(wc -l < antisense_names.txt)
    echo "  Ensembl entries: $(grep -c '^>' "$ENS_ALL"); stored antisense: ${n_anti}"
    [ "$n_anti" -eq 0 ] || sed 's/^/    flipping: /' antisense_names.txt

    # Reverse-complement exactly those entries; leave every other byte alone.
    awk 'function rc(s,   i,c,o) {
             o = ""
             for (i = length(s); i > 0; i--) {
                 c = substr(s, i, 1)
                 o = o ( c=="A"?"T": c=="T"?"A": c=="C"?"G": c=="G"?"C":
                         c=="a"?"t": c=="t"?"a": c=="c"?"g": c=="g"?"c": "N" )
             }
             return o
         }
         function flush_rec() {
             if (name == "") return
             print name
             print (hit ? rc(seq) : seq)
         }
         NR == FNR { bad[$1] = 1; next }
         /^>/ {
             flush_rec()
             name = $0
             # the FASTA id is the header up to the first whitespace, minus ">"
             id = substr($1, 2)
             hit = (id in bad)
             seq = ""
             next
         }
         { seq = seq $0 }
         END { flush_rec() }
        ' antisense_names.txt "$ENS_ALL" > ensembl_oriented.fa
    ENS_ALL=ensembl_oriented.fa

    # Prove it: re-align and require zero reverse hits.
    n_left=$(bwa mem -t 4 units.fa "$ENS_ALL" 2>/dev/null | samtools view - \
             | awk '$2==16' | wc -l)
    if [ "$n_left" -ne 0 ]; then
        echo "  FATAL: ${n_left} entries still antisense after the flip" >&2
        exit 1
    fi
    echo "  verified: 0 entries remain antisense to the units"
else
    n_anti=NA
    echo "  SKIPPED (ORIENT_TO_UNITS=no) -- reproduces v2, backwards entry included"
fi

###############################################################################
echo "[$(date)] === STEP 5: concatenate + bwa index ==="
###############################################################################
if [ ! -s "$FA_OUT" ]; then
    cat human_45S.fa mouse_rdna_47S.fa "$ENS_ALL" > "${FA_OUT}.tmp"
    mv "${FA_OUT}.tmp" "$FA_OUT"
fi
[ -s "${FA_OUT}.bwt" ] || bwa index "$FA_OUT"
echo "  ${FA_OUT}: $(grep -c '^>' "$FA_OUT") seqs"

###############################################################################
echo "[$(date)] === STEP 6: QC interval table ==="
###############################################################################
# Human subunit coordinates differ per RNA45SN copy by a few nt, so only the
# mouse unit gets exact intervals here; for human, decompose by which of the
# five transcripts a read hit and use the RefSeq feature table if you need
# subunit resolution.
cat > "${OUT}/rrna_intervals_mixed.tsv" <<EOF
# subunit intervals within mouse_rDNA_47S_${MOUSE_RDNA_ACC}_${MOUSE_RDNA_TX_START}-${MOUSE_RDNA_TX_END}
# 1-based inclusive, from the ${MOUSE_RDNA_ACC} GenBank feature table
#subunit	start	end	length
5ETS	1	4007	4007
18S	4008	5877	1870
ITS1	5878	6877	1000
5.8S	6878	7034	157
ITS2	7035	8122	1088
28S	8123	12852	4730
3ETS	12853	13403	551
EOF
echo "  wrote ${OUT}/rrna_intervals_mixed.tsv"

###############################################################################
echo "[$(date)] === STEP 7: provenance ==="
###############################################################################
{
  echo "unique_rRNA_human_mouse.v2.fa"
  echo "built    : $(date -Iseconds)"
  echo "by       : $(basename "$0") (tracked in code/I_Gene_expression/own_version/)"
  echo "host     : $(hostname)"
  echo
  echo "component A -- Ensembl 99 (GRCh38 + GRCm38)"
  echo "  filter : gene feature, gene_biotype in {rRNA, Mt_rRNA, rRNA_pseudogene}"
  echo "  human  : $(grep -c '^>' human.rRNA.fa) seqs"
  echo "  mouse  : $(grep -c '^>' mouse.rRNA.fa) seqs"
  echo
  echo "component B -- NCBI, human"
  for i in "${!HUMAN_45S_ACC[@]}"; do
    echo "  ${HUMAN_45S_NAME[$i]} : ${HUMAN_45S_ACC[$i]} (${HUMAN_45S_LEN[$i]} bp, RefSeq transcript, no IGS)"
  done
  echo
  echo "component B -- NCBI, mouse"
  echo "  accession : ${MOUSE_RDNA_ACC} (${got} bp, full repeating unit)"
  echo "  extracted : ${MOUSE_RDNA_TX_START}-${MOUSE_RDNA_TX_END} (transcribed 47S; IGS excluded)"
  echo
  echo "orientation : ORIENT_TO_UNITS=${ORIENT_TO_UNITS}"
  if [ "$ORIENT_TO_UNITS" = "yes" ]; then
      echo "  every Ensembl entry aligned to the NCBI 45S/47S units; any coming"
      echo "  back flag 16 (antisense) was reverse-complemented, then re-checked"
      echo "  to 0 remaining. This matters because riboread-selection.py decides"
      echo "  what is ribosomal BY STRAND, so a backwards entry shields reads from"
      echo "  depletion instead of removing them."
      echo "  entries flipped : ${n_anti}"
      [ "${n_anti}" = "0" ] || sed 's/^/    /' "${B}/antisense_names.txt"
  else
      echo "  NOT APPLIED -- this file may contain backwards entries"
  fi
  echo "total seqs  : $(grep -c '^>' "$FA_OUT")"
  echo "total bp    : $(awk '/^>/{next}{n+=length($0)}END{print n}' "$FA_OUT")"
  echo
  echo "supersedes  : unique_rRNA_human_mouse.v2.fa, in which ENSMUSG00000106106"
  echo "              was stored antisense, sparing 95,823 UFIs worth of genuine"
  echo "              18S reads on that locus (published table: 72) and"
  echo "              understating step-3 depletion by ~0.28 points"
  echo "replaces    : unique_rRNA_human_mouse.fa (v1, Ensembl only), which called"
  echo "              0.27% of the vasaplate control's reads ribosomal"
  echo "consumed by : ribo-bwamem.sh (bwa aln + bwa mem)"
} > "${FA_OUT%.fa}.provenance.txt"
echo "  wrote ${FA_OUT%.fa}.provenance.txt"

echo "[$(date)] === DONE ==="
echo "riboref=${FA_OUT}"
