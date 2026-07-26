#!/bin/bash
###############################################################################
# trim.sh -- step 2, forked from ../a_Mapping/trim.sh
#
# Three passes now, and the SAME output filename as upstream, so steps 3-7 are
# untouched.
#
#   pass 0  trim_bc_anchor.py -- cut the 3' technical tail by FINDING it
#   pass 1  TrimGalore        -- adapters + Phred<20   (upstream's, unchanged)
#   pass 2  cutadapt          -- whatever pass 0 could not anchor, + poly-A/G/T
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
# STAR's outFilterMatchNminOverLread.
#
# PASS 0 IS THE INTERESTING ONE
# -----------------------------
# The tail never had to be recognised by its shape. Step 1 put the cell barcode
# and UMI on the read name, so for any given read the 12 nt after the poly-A
# are a known literal string. trim_bc_anchor.py finds it, drops everything from
# there to the 3' end, and walks back over the poly-A. Nothing to tune.
#
# Earlier revisions did this inside cutadapt, which cannot take a per-read
# pattern: the UMI became 6 wildcards, wildcards had to run at min_overlap =
# the whole pattern length or they ate the end of every read, and that forced
# 21 nt of adapter into the pattern so it only fired when the read had gone
# that far. All of that machinery is gone. Pass 0 matches 29.6% of reads
# against the wildcard version's 26.9%, needs no adapter at all, and reaches
# the case the wildcard version structurally could not -- a read that ENDS
# partway through the barcode. Equal protein-coding exonic yield (40,684 vs
# 40,678) at higher purity (55.5% vs 54.2% exonic, 87.1% vs 85.1% annotated).
#
# See trimtest/ for the benchmark, and README.md "Step 2" for the numbers.
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
TRIM_ANCHOR_BC="${TRIM_ANCHOR_BC:-yes}"        # run pass 0
TRIM_ANCHOR_PY="${TRIM_ANCHOR_PY:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/trim_bc_anchor.py}"
TRIM_PYTHON="${TRIM_PYTHON:-${TRIM_CUTADAPT%/*}/python3}"
# an if, not "${VAR:-A{20}}" -- the `}` in the default would close the expansion
# early and leave a stray brace in the adapter. See config.sh for the long note.
[ -n "${TRIM_POLYA+x}" ]  || TRIM_POLYA='A{20}'
[ -n "${TRIM_POLYG+x}" ]  || TRIM_POLYG='G{20}'
[ -n "${TRIM_POLYT5+x}" ] || TRIM_POLYT5='T{20}'

stem=${file2trim%.fastq.gz}

# =============================================================================
# legacy -- byte-for-byte upstream, kept so the two can be compared
# =============================================================================
if [ "$TRIM_MODE" = "legacy" ]; then
    ${path2trimgalore}/trim_galore --path_to_cutadapt "${path2cutadapt}/cutadapt" \
        "${file2trim}" -o "${outdir}" || exit 1
    mv "${file2trim}_trimming_report.txt" "${stem}_trimming_report.txt"
    "${path2cutadapt}/cutadapt" -m 15 --trim-n \
        -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" \
        -o "${stem}_trimmed_homoATCG.fq.gz" "${stem}_trimmed.fq.gz"
    exit $?
fi

# =============================================================================
# pass 0 -- cut the tail at the barcode
# =============================================================================
tg_in="$file2trim"
p0=""
if [ "$TRIM_ANCHOR_BC" = "yes" ]; then
    [ -x "$TRIM_ANCHOR_PY" ] || { echo "trim.sh: ${TRIM_ANCHOR_PY} missing or not executable" >&2; exit 1; }
    p0="${stem}.bcanchor.fastq.gz"
    # `env -u PYTHONPATH` for the same reason as pass 2 below -- the conda
    # interpreter must not inherit the Trim_Galore module's PYTHONPATH.
    env -u PYTHONPATH "$TRIM_PYTHON" "$TRIM_ANCHOR_PY" \
        "$file2trim" "$p0" --log "${stem}_bcanchor.log" || exit 1
    tg_in="$p0"
fi

# =============================================================================
# pass 1 -- TrimGalore, unchanged from upstream
# =============================================================================
${path2trimgalore}/trim_galore --path_to_cutadapt "${path2cutadapt}/cutadapt" \
    "${tg_in}" -o "${outdir}" || exit 1
# TrimGalore names its outputs after the input basename, so when pass 0 ran they
# come out as *.bcanchor_trimmed.fq.gz. Put them back on the expected stem and
# drop the intermediate -- everything downstream keys off ${stem}_trimmed*.
if [ -n "$p0" ]; then
    mv "${stem}.bcanchor_trimmed.fq.gz" "${stem}_trimmed.fq.gz"
    mv "${p0}_trimming_report.txt"      "${stem}_trimming_report.txt"
    rm -f "$p0"
else
    mv "${file2trim}_trimming_report.txt" "${stem}_trimming_report.txt"
fi

# =============================================================================
# pass 2 -- cutadapt, for what pass 0 could not anchor
# =============================================================================
# Pass 0 handles the ~30% of reads that reached their own barcode. The rest
# either never reached the tail (long insert -- nothing to trim) or carry too
# many errors in the barcode to match, and those are what these are for:
#
#   rt     the measured read-through construct. min_overlap 8 is high enough
#          that it will not fire on a chance 3-mer at the read end, which is
#          cutadapt's default and is what makes short adapters dangerous.
#   polyA  requires a real 20 nt run (or 10 nt running off the read end), so it
#          removes the protocol's poly-A tail without touching the 6-mers that
#          occur in ordinary mRNA. This is what stops poly-A reads aligning to
#          genomic A-tracts -- without it 16.6% of uniquely mapped reads are
#          poly-A stuck to a handful of loci, with it 3.7%.
#   polyG  the two-colour "no signal" artefact. Costs nothing on this run
#          (18 reads in 300k); kept as free insurance on a NovaSeq X.
#   polyT5 a read in the REVERSE orientation reads the poly-A as poly-T at its
#          START, so this is a -g adapter: it removes the match and everything
#          BEFORE it. It deletes ~3,400 uniquely mapped reads per 300k of which
#          ~98% are non-exonic, leaving the exonic count flat -- and it is what
#          makes the blank barcodes look blank. The matching reverse-orientation
#          barcode anchor was built and measured (trimtest/bench_trim8.sh) and
#          fired 13 times in 300,000 reads, so it is deliberately absent.
# -n 3 so poly-A can still be reached after the read-through is removed.
opts=(-m "${TRIM_MINLEN}" --trim-n -n 3 --poly-a)
[ -n "$TRIM_ADAPTER3" ] && opts+=(-a "rt=${TRIM_ADAPTER3};min_overlap=8")
[ -n "$TRIM_POLYA"    ] && opts+=(-a "polyA=${TRIM_POLYA};min_overlap=10")
[ -n "$TRIM_POLYG"    ] && opts+=(-a "polyG=${TRIM_POLYG};min_overlap=10")
[ -n "$TRIM_POLYT5"   ] && opts+=(-g "polyT5=${TRIM_POLYT5};min_overlap=10")

# `env -u PYTHONPATH` is load-bearing. Loading the Trim_Galore module puts
# cutadapt 1.18's site-packages on PYTHONPATH, so the conda env's cutadapt 5.1
# launcher imports the 1.18 package instead of its own and dies with
# "No module named 'cutadapt.cli'". Unsetting it lets the launcher's interpreter
# find its own package. Pass 1 above is unaffected -- it wants 1.18.
env -u PYTHONPATH "${TRIM_CUTADAPT}" "${opts[@]}" \
    -o "${stem}_trimmed_homoATCG.fq.gz" "${stem}_trimmed.fq.gz"
