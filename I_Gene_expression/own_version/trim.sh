#!/bin/bash
###############################################################################
# trim.sh -- step 2, forked from ../a_Mapping/trim.sh
#
# Same two passes as upstream (TrimGalore, then cutadapt) and the SAME output
# filename, so steps 3-7 are untouched. Only the second pass differs.
#
# WHY IT WAS FORKED
# -----------------
# Upstream's second pass is
#
#     cutadapt -m 15 --trim-n -a polyG1=GG{5} -a polyC1=CC{5} \
#                             -a polyT1=TT{5} -a polyA1=AA{5}
#
# i.e. "cut at the first run of 6 identical bases, and throw away everything
# after it". As a 3' adapter a 6-mer matches the FIRST time it occurs, so a
# perfectly good 130 nt read with an internal AAAAAA at position 10 is cut to
# 10 nt and then dropped by -m 15. Measured on this library: of the reads that
# carry no poly-A tail at all, that rule truncates 49% of them, median loss 120
# of 130 bases.
#
# It also never removed the thing that is actually in the way. This library's
# 3' read-through is
#
#     [insert][poly-A][12 nt = revcomp(CBC+UMI)][revcomp(R1 5' prefix) + Nextera]
#
# TrimGalore auto-detects only the Nextera part, which sits ~26 nt too far
# downstream, so every short-insert read kept ~40 nt of adapter and failed
# STAR's outFilterMatchNminOverLread. Naming that adapter explicitly is where
# nearly all of the gain comes from.
#
# See ../../../data/PM26037/trimtest/ for the benchmark this is based on, and
# README.md "Step 2: what is trimmed and why" for the numbers.
#
# Usage (positional args are upstream's; the rest come from config.sh):
#   trim.sh <in.fastq.gz> <outdir> <path2trimgalore> <path2cutadapt>
###############################################################################

set -uo pipefail

if [ $# -ne 4 ]; then
    echo "Please, give (1) input fastq file; (2) output directory; (3) path2trimgalore; (4) path2cutadapt;"
    exit 1
fi

file2trim=$1
outdir=$2
path2trimgalore=$3
path2cutadapt=$4

# --- settings, normally exported by config.sh --------------------------------
TRIM_MODE="${TRIM_MODE:-vasa}"                 # vasa | legacy
TRIM_CUTADAPT="${TRIM_CUTADAPT:-${path2cutadapt}/cutadapt}"   # pass 2 binary
TRIM_ADAPTER3="${TRIM_ADAPTER3:-}"
TRIM_MINLEN="${TRIM_MINLEN:-20}"
# an if, not "${VAR:-A{20}}" -- the `}` in the default would close the expansion
# early and leave a stray brace in the adapter. See config.sh for the long note.
[ -n "${TRIM_POLYA+x}" ] || TRIM_POLYA='A{20}'
[ -n "${TRIM_POLYG+x}" ] || TRIM_POLYG='G{20}'

stem=${file2trim%.fastq.gz}

# --- pass 1: adapters + Phred<20, unchanged from upstream --------------------
${path2trimgalore}/trim_galore --path_to_cutadapt "${path2cutadapt}/cutadapt" \
    "${file2trim}" -o "${outdir}" || exit 1
mv "${file2trim}_trimming_report.txt" "${stem}_trimming_report.txt"

# --- pass 2 --------------------------------------------------------------
if [ "$TRIM_MODE" = "legacy" ]; then
    # byte-for-byte the upstream command, kept so the two can be compared
    "${path2cutadapt}/cutadapt" -m 15 --trim-n \
        -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" \
        -o "${stem}_trimmed_homoATCG.fq.gz" "${stem}_trimmed.fq.gz"
    exit $?
fi

# vasa mode. Every adapter is a literal sequence with an explicit min_overlap:
#   rt     the measured read-through construct. min_overlap 8 is high enough
#          that it will not fire on a chance 3-mer at the read end, which is
#          cutadapt's default and is what makes short adapters dangerous.
#   polyA  requires a real 20 nt run (or 10 nt running off the read end), so it
#          removes the protocol's poly-A tail without touching the 6-mers that
#          occur in ordinary mRNA. This is the one line that stops poly-A reads
#          aligning to genomic A-tracts -- without it 16.6% of uniquely mapped
#          reads are poly-A stuck to a handful of loci, with it 3.7%.
#   polyG  the two-colour "no signal" artefact. Costs nothing on this run
#          (18 reads in 300k); kept because it is free insurance on a NovaSeq X.
# -n 2 so the poly-A can still be reached after the read-through is removed.
opts=(-m "${TRIM_MINLEN}" --trim-n -n 2)
[ -n "$TRIM_ADAPTER3" ] && opts+=(-a "rt=${TRIM_ADAPTER3};min_overlap=8")
[ -n "$TRIM_POLYA"    ] && opts+=(-a "polyA=${TRIM_POLYA};min_overlap=10")
[ -n "$TRIM_POLYG"    ] && opts+=(-a "polyG=${TRIM_POLYG};min_overlap=10")

# `env -u PYTHONPATH` is load-bearing. Loading the Trim_Galore module puts
# cutadapt 1.18's site-packages on PYTHONPATH, so the conda env's cutadapt 5.1
# launcher imports the 1.18 package instead of its own and dies with
# "No module named 'cutadapt.cli'". Unsetting it lets the launcher's interpreter
# find its own package. Pass 1 above is unaffected -- it wants 1.18.
env -u PYTHONPATH "${TRIM_CUTADAPT}" "${opts[@]}" \
    -o "${stem}_trimmed_homoATCG.fq.gz" "${stem}_trimmed.fq.gz"
