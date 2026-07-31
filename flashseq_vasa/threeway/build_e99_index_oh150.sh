#!/bin/bash
#SBATCH -J e99_oh150_index
#SBATCH -p ncpu
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 8:00:00
###############################################################################
# build_e99_index_oh150.sh
#
# WHY THIS EXISTS
# ---------------
# The published VASA plate (SRR14783059) was mapped to a species-prefixed
# GRCh38+GRCm38 Ensembl 99 reference; the own plate (ZHA9292A1) was mapped to
# GRCm39 / Ensembl 116. Every published-vs-own gap therefore mixes protocol
# biology with annotation release. To remove the release term instead of
# correcting for it, the own plate has to be re-quantified against the SAME
# Ensembl 99 reference.
#
# The existing mixed index reference/vasaseq/mixed/star_index_74 CANNOT be
# reused for that: it was built with sjdbOverhang 73, correct for the published
# library's 75 nt reads and much too small for the own plate's 130 nt reads.
# sjdbOverhang is baked into STAR's Genome/SA files (junction sequences are
# appended to the genome), so it cannot be overridden at mapping time.
# pipeline.sh step_check enforces sjdbOverhang >= readlen-skip-1 = 129 and would
# refuse the oh73 index. Mapping 130 nt reads against oh73 does not error -- it
# silently loses junction-spanning alignments, which is exactly the artefact
# this control is supposed to eliminate.
#
# WHAT IS HELD CONSTANT
# ---------------------
# genomeFastaFiles and sjdbGTFfile are the SAME FILES the published index was
# built from (combined_genome.fa, combined.gtf, checksummed below), so the
# genome sequence and the gene models are identical, not merely the same
# release. STAR version is pinned to 2.7.7a, as used for both existing indices.
# The ONLY parameter that differs from star_index_74 is sjdbOverhang: 150 here
# vs 73 there, each sized correctly for the read length it will serve. That is
# a per-library correctness setting, not a comparison variable -- using 73 for
# 130 nt reads would be the error.
#
# Naming follows the existing convention (star_index_<readlen>): READLEN 151 ->
# OVERHANG 150, mirroring GRCm39/star_index_151_r116 which the own plate's
# Ensembl 116 arm used, so the own-E99 and own-E116 arms carry the same
# junction-overhang budget and differ only in genome build + annotation.
#
# Idempotent: re-running skips the build if SA is already present.
###############################################################################
set -euo pipefail

REF=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mixed
BUILD="${REF}/build"
READLEN=151
OVERHANG=$(( READLEN - 1 ))
IDX="${REF}/star_index_${READLEN}"

echo "[$(date)] === precheck ==="
for f in "${BUILD}/combined_genome.fa" "${BUILD}/combined.gtf"; do
    [ -s "$f" ] || { echo "FATAL: missing $f"; exit 1; }
    printf '  ok %s  %s bytes\n' "$(basename "$f")" "$(stat -c%s "$f")"
done
# Disk: the oh73 index is ~56 GB; oh150 appends more junction sequence, so
# budget 80 GB and refuse rather than fill a 92%-full filesystem.
# NB `df -P --output=` is rejected (mutually exclusive on this coreutils); take
# the POSIX 4th column instead.
avail_kb=$(df -P "$REF" | awk 'NR==2{print $4}')
echo "  avail on $REF: $(( avail_kb / 1024 / 1024 )) GB"
[ "$avail_kb" -gt $(( 80 * 1024 * 1024 )) ] || { echo "FATAL: <80 GB free"; exit 1; }

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/camp/apps/eb/modules/all
module purge 2>/dev/null || true
module load STAR/2.7.7a-GCC-10.2.0
echo "  STAR $(STAR --version)"

if [ -s "${IDX}/SA" ] && [ -s "${IDX}/Genome" ] && [ -s "${IDX}/SAindex" ]; then
    echo "[$(date)] index already present, skipping build: ${IDX}"
else
    echo "[$(date)] === STAR genomeGenerate, sjdbOverhang ${OVERHANG} ==="
    mkdir -p "$IDX"
    cd "$BUILD"
    STAR --runMode genomeGenerate \
         --runThreadN 32 \
         --genomeDir "$IDX" \
         --genomeFastaFiles combined_genome.fa \
         --sjdbGTFfile combined.gtf \
         --sjdbOverhang "$OVERHANG" \
         --limitGenomeGenerateRAM 190000000000 \
         --outFileNamePrefix "${IDX}/"
fi

echo "[$(date)] === verify ==="
# The index must agree with the Ensembl 99 BED that step 5 will intersect
# against, or every read silently fails to assign. Both must be species-
# prefixed and must carry the same 260 contigs.
n_chr=$(wc -l < "${IDX}/chrName.txt")
echo "  contigs: ${n_chr}"
[ "$n_chr" -eq 260 ] || { echo "FATAL: expected 260 contigs, got ${n_chr}"; exit 1; }
n_m=$(grep -c '^GRCm38_' "${IDX}/chrName.txt" || true)
n_h=$(grep -c '^GRCh38_' "${IDX}/chrName.txt" || true)
echo "  GRCm38_ contigs: ${n_m}   GRCh38_ contigs: ${n_h}"
[ "$n_m" -gt 0 ] && [ "$n_h" -gt 0 ] || { echo "FATAL: species prefixes absent"; exit 1; }
got_ovh=$(awk '$1=="sjdbOverhang"{print $2; exit}' "${IDX}/genomeParameters.txt")
echo "  sjdbOverhang recorded: ${got_ovh}"
[ "$got_ovh" -eq "$OVERHANG" ] || { echo "FATAL: overhang ${got_ovh} != ${OVERHANG}"; exit 1; }
# Chromosome-name identity against the BED step 5 uses -- not just "both are
# prefixed", but the same set. A silent mismatch here is the single most
# expensive failure mode in this whole exercise.
BED="${REF}/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed"
cut -f1 "$BED" | sort -u > /tmp/bed_chr.$$
sort -u "${IDX}/chrName.txt" > /tmp/idx_chr.$$
missing=$(comm -23 /tmp/bed_chr.$$ /tmp/idx_chr.$$ | wc -l)
echo "  BED contigs absent from index: ${missing}"
rm -f /tmp/bed_chr.$$ /tmp/idx_chr.$$
[ "$missing" -eq 0 ] || { echo "FATAL: BED references contigs the index lacks"; exit 1; }

{
  echo "star_index_${READLEN}  (Ensembl 99, GRCh38+GRCm38 species-prefixed)"
  echo "built              : $(date -Iseconds)"
  echo "built_by           : flashseq_vasa/threeway/build_e99_index_oh150.sh"
  echo "STAR               : $(STAR --version)"
  echo "sjdbOverhang       : ${OVERHANG}   (star_index_74 uses 73)"
  echo "genomeFastaFiles   : ${BUILD}/combined_genome.fa"
  echo "sjdbGTFfile        : ${BUILD}/combined.gtf"
  echo "contigs            : ${n_chr}  (GRCm38_ ${n_m}, GRCh38_ ${n_h})"
  echo "purpose            : re-quantify own plate ZHA9292A1 under the published"
  echo "                     plate's annotation, to remove the release confound"
  echo "sha256 combined_genome.fa : $(sha256sum "${BUILD}/combined_genome.fa" | cut -c1-32)"
  echo "sha256 combined.gtf       : $(sha256sum "${BUILD}/combined.gtf"       | cut -c1-32)"
} > "${IDX}/PROVENANCE.txt"
cat "${IDX}/PROVENANCE.txt"
cp "${IDX}/PROVENANCE.txt" ./e99_index_provenance.txt
du -sh "$IDX"
echo "[$(date)] === DONE ==="
