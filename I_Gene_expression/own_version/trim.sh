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
TRIM_ANCHOR_BC="${TRIM_ANCHOR_BC:-yes}"        # anchor on the cell barcode
TRIM_LEN_UMI="${TRIM_LEN_UMI:-${LEN_UMI:-6}}"  # wildcards between CBC and adapter
TRIM_ANCHOR_ADLEN="${TRIM_ANCHOR_ADLEN:-21}"   # adapter nt carried by the anchor
# an if, not "${VAR:-A{20}}" -- the `}` in the default would close the expansion
# early and leave a stray brace in the adapter. See config.sh for the long note.
[ -n "${TRIM_POLYA+x}" ] || TRIM_POLYA='A{20}'
[ -n "${TRIM_POLYG+x}" ] || TRIM_POLYG='G{20}'
[ -n "${TRIM_POLYT5+x}" ] || TRIM_POLYT5='T{20}'

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
#   bcumi  the barcode anchor, built below. This is the one that knows where
#          the tail starts rather than guessing from its shape.
# -n 3 so poly-A can still be reached after the anchor and read-through go.
opts=(-m "${TRIM_MINLEN}" --trim-n -n 3)

# --- the barcode anchor ------------------------------------------------------
# The tail does not have to be recognised by shape -- we know what follows it.
# concatenator.py put the cell barcode and UMI on the read name at step 1, and
# within one cell's fastq the barcode is CONSTANT (verified: a per-cell file has
# exactly one distinct CB tag). So per cell the whole tail is a fixed pattern
# with only the UMI unknown:
#
#   revcomp(CBC)  N x LEN_UMI  <first TRIM_ANCHOR_ADLEN nt of the adapter>
#    6 specific     6 any               21 specific       = 27 specific of 33
#
# TRIM_ANCHOR_ADLEN defaults to 21 because TRIM_ADAPTER3[:21] is EXACTLY
# revcomp(the 21 nt prefix stripped from R1) -- everything past 21 is the
# Nextera mosaic end. So it is the same 21 nt at both ends of the read and there
# is no magic number to remember: it is SKIP5. Measured plateau (round 9):
# 12, 16 and 21 all give 40,678-40,681 protein-coding exonic reads, i.e. equal
# to within 3 reads, so the principled boundary costs nothing. Outside that
# window the anchor stops working -- at 26 the pattern is too long to fit inside
# most reads (fires 12,760 times vs 80,697 at 21), and at 8 it is short enough
# that the 40 nt `rt` adapter outscores it and wins the match instead (12,080).
#
# min_overlap is set to the FULL length, which forbids partial 3'-end matching.
# That matters: partial matching plus wildcards is what made an earlier attempt
# eat the end of every read (see README, "What did not work"). With full-length
# matching required, 27 specific bases make a chance hit impossible.
#
# What this unlocks is the poly-A. Once the anchor is removed the poly-A is at
# the very 3' END of the read, which is the only place --poly-a looks, so tails
# SHORTER than 20 nt -- which A{20} structurally cannot catch -- are cleaned
# too. Barcode remnant in the output: 13.6% without this, 2.2% with it.
if [ "$TRIM_ANCHOR_BC" = "yes" ] && [ -n "$TRIM_ADAPTER3" ]; then
    # take the barcode off the first read's CB: tag rather than the filename or
    # the whitelist order, so this cannot silently pair a cell with the wrong
    # barcode
    cb=$(zcat "$file2trim" 2>/dev/null | head -1 \
         | tr ';' '\n' | awk -F: '$1=="CB"{print $2; exit}')
    if [ -n "$cb" ]; then
        cbrc=$(printf '%s' "$cb" | tr ACGTacgt TGCAtgca | rev)
        umiN=$(printf 'N%.0s' $(seq 1 "$TRIM_LEN_UMI"))
        anchor="${cbrc}${umiN}${TRIM_ADAPTER3:0:$TRIM_ANCHOR_ADLEN}"
        opts+=(-a "bcumi=${anchor};min_overlap=${#anchor}")
        opts+=(--poly-a)
    else
        echo "trim.sh: no CB: tag on the first read of ${file2trim}," \
             "skipping the barcode anchor" >&2
    fi
fi

[ -n "$TRIM_ADAPTER3" ] && opts+=(-a "rt=${TRIM_ADAPTER3};min_overlap=8")
[ -n "$TRIM_POLYA"    ] && opts+=(-a "polyA=${TRIM_POLYA};min_overlap=10")
[ -n "$TRIM_POLYG"    ] && opts+=(-a "polyG=${TRIM_POLYG};min_overlap=10")

# 5' poly-T. A read in the reverse orientation reads the poly-A tail as poly-T
# at its START, so this is a -g (5') adapter: it removes the match and
# everything BEFORE it. It removes ~3,400 uniquely mapped reads per 300k from a
# real cell of which ~98% were non-exonic, leaving the protein-coding exonic
# count flat -- i.e. it deletes junk and nothing else. On the blank barcode 016
# it is the difference between 47,284 and 12,396 uniquely mapped reads: the
# blank finally looks like a blank.
#
# The mirror-image barcode anchor for reverse-orientation reads is deliberately
# NOT here. It was built and measured (trimtest/bench_trim8.sh) and fired 13
# times in 300,000 reads: this library has no meaningful reverse population, and
# all of the gain above comes from the poly-T alone.
[ -n "$TRIM_POLYT5"   ] && opts+=(-g "polyT5=${TRIM_POLYT5};min_overlap=10")

# `env -u PYTHONPATH` is load-bearing. Loading the Trim_Galore module puts
# cutadapt 1.18's site-packages on PYTHONPATH, so the conda env's cutadapt 5.1
# launcher imports the 1.18 package instead of its own and dies with
# "No module named 'cutadapt.cli'". Unsetting it lets the launcher's interpreter
# find its own package. Pass 1 above is unaffected -- it wants 1.18.
env -u PYTHONPATH "${TRIM_CUTADAPT}" "${opts[@]}" \
    -o "${stem}_trimmed_homoATCG.fq.gz" "${stem}_trimmed.fq.gz"
