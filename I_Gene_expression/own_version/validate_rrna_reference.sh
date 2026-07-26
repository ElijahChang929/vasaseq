#!/bin/bash
#SBATCH -J vasa_rrnaval
#SBATCH -p ncpu
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 4:00:00
###############################################################################
# validate_rrna_reference.sh -- does the new rRNA reference actually catch more?
#
# Separate from the two build scripts on purpose: this one runs REAL READS
# through the REAL step 3 (ribo-bwamem.sh + riboread-selection.py, unmodified),
# once per reference, and diffs the result. A reference that looks right on
# paper is not evidence; the depletion rate on the authors' own control is.
#
# Default subject is the paper's vasaplate-HEK293T-mESC run, whose 384 cells
# were depleted against the Ensembl-only v1 and came out at 0.27% ribosomal --
# the number that started this.
#
# Nothing is overwritten: each reference's output goes to its own scratch dir,
# and the pipeline's own vasaplate_out/ is read but never written.
#
# Usage:
#   ./validate_rrna_reference.sh [NCELLS] [WORKDIR]
#     NCELLS   how many cells to test (default 8; 384 = the whole run)
#
#   REFS="/path/a.fa /path/b.fa" ./validate_rrna_reference.sh
#     compare an arbitrary list instead of the mixed v1/v2 pair
###############################################################################
set -uo pipefail

NCELLS=${1:-8}
WORK=${2:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/ref/fastq_vasaplate/rrna_validation}

EBROOT=/camp/apps/eb/software
VASA_REFS=${VASA_REFS:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq}
P2S=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/a_Mapping
CELLS_DIR=${CELLS_DIR:-/nemo/lab/turnerj/working/guangxin/vasaseq/data/ref/fastq_vasaplate/vasaplate_out}

REFS=${REFS:-"${VASA_REFS}/mixed/unique_rRNA_human_mouse.fa ${VASA_REFS}/mixed/unique_rRNA_human_mouse.v2.fa"}

P2BWA=${EBROOT}/BWA/0.7.17-GCC-10.3.0/bin
P2SAMTOOLS=${EBROOT}/SAMtools/1.11-GCC-10.2.0/bin
CONDA_ENV=${CONDA_ENV:-/nemo/lab/turnerj/working/guangxin/envs/vasa}

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=${EBROOT%/software}/modules/all
# ribo needs BOTH: bwa/samtools from modules (called by absolute path) and the
# conda env for riboread-selection.py's pysam. Conda's python leads on PATH, the
# binaries resolve absolutely, so this pairing is safe -- see CLAUDE.md.
module purge 2>/dev/null || true
module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0
source ${EBROOT}/Anaconda3/2024.10-1/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

mapfile -t CELLFILES < <(ls "${CELLS_DIR}"/*_cbc_trimmed_homoATCG.fq.gz 2>/dev/null | head -n "$NCELLS")
[ ${#CELLFILES[@]} -gt 0 ] || { echo "no trimmed fastqs under ${CELLS_DIR}" >&2; exit 1; }
echo "subject : ${CELLS_DIR}"
echo "cells   : ${#CELLFILES[@]}"
echo

for ref in $REFS; do
    [ -s "$ref" ] || { echo "MISSING reference: $ref" >&2; exit 1; }
    [ -s "${ref}.bwt" ] || { echo "MISSING bwa index for: $ref" >&2; exit 1; }
    tag=$(basename "$ref" .fa)
    d="${WORK}/${tag}"
    mkdir -p "$d"
    echo "=== ${tag} ($(grep -c '^>' "$ref") seqs) ==="
    for fq in "${CELLFILES[@]}"; do
        cell=$(basename "$fq" _cbc_trimmed_homoATCG.fq.gz)
        [ -s "${d}/${cell}.ribo-map.log" ] && continue
        "${P2S}/ribo-bwamem.sh" "$ref" "$fq" "${d}/${cell}" \
            "$P2BWA" "$P2SAMTOOLS" y "$P2S" > "${d}/${cell}.stdout" 2>&1 \
            || echo "  FAILED: $cell"
        # ribo-bwamem.sh writes <prefix>.ribo-map.log next to its other outputs
        [ -s "${d}/${cell}.ribo-map.log" ] || mv "${d}/${cell}"*ribo-map.log \
            "${d}/${cell}.ribo-map.log" 2>/dev/null || true
    done
    python - "$d" <<'PY'
import glob, re, sys
d = sys.argv[1]
tot = unm = 0
for f in glob.glob(f"{d}/*ribo-map.log"):
    t = open(f, errors="replace").read()
    a = re.search(r"Number of reads:\s*(\d+)", t)
    b = re.search(r"Number of unmapped reads:\s*(\d+)", t)
    if a and b:
        tot += int(a.group(1)); unm += int(b.group(1))
if tot:
    print(f"  reads in        : {tot:,}")
    print(f"  kept (non-ribo) : {unm:,}")
    print(f"  called ribosomal: {tot-unm:,}  ({100*(tot-unm)/tot:.2f}%)")
else:
    print("  no logs parsed")
PY
    echo
done
echo "outputs under ${WORK}/ -- delete it to re-run"
