#!/bin/bash
###############################################################################
# 05_rrna_bwa.sh -- measure FLASH-seq rRNA the way the VASA run measured it.
#
# WHY THIS EXISTS
# ---------------
# 01_rrna_kmer_screen.py gives 3.6-6.3% rRNA across these libraries, but that
# number is an exact-31-mer-containment LOWER BOUND, and the VASA figure it is
# meant to be compared against (21.39% over ZHA9292A1, own_version/README.md)
# came from bwa. Two different methods is not a comparison. This script runs the
# VASA pipeline's OWN rRNA stage -- a_Mapping/ribo-bwamem.sh and
# riboread-selection.py, unmodified, against the same unique_rRNA_mouse.v2.fa --
# over the FLASH-seq reads, so the two percentages are the same measurement.
#
# Nothing about rRNA detection is reimplemented here. This script only prepares
# input for those scripts and tabulates what they wrote.
#
# THE THREE CHOICES THAT ARE NOT OBVIOUS
# --------------------------------------
# 1. stranded = n, not y. riboread-selection.py's `stranded` flag decides
#    whether a reverse-mapping read may be called ribosomal: with `y` only
#    forward-mapping reads count. VASA is a stranded protocol so `y` is right
#    there. FLASH-seq is unstranded -- nf-core's own salmon inference said
#    `unstranded` on all ten libraries -- so `y` would discard roughly half the
#    rRNA reads for no reason. We run `n`, and then re-run riboread-selection.py
#    with `y` over the SAME merged BAM (no realignment) so the strand split is
#    measured rather than asserted. If the library really is unstranded the two
#    should differ by ~2x; that is a free internal check on the flag choice.
#
# 2. Two arms: raw and trimmed, over the same sampled reads.
#      raw     -- untrimmed, so it is comparable with 01_rrna_kmer_screen.py and
#                 the k-mer/bwa gap is a difference of method alone.
#      trimmed -- comparable DENOMINATOR to VASA, whose 21.39% is over trimmed
#                 reads. 16-32% of raw FLASH-seq R1 carries adapter read-through and
#                 poly-G (README.md), and those reads cannot map to rRNA, so
#                 leaving them in deflates the percentage.
#    The trimmed arm is the headline number. Both are reported: the two
#    protocols' trimming is not the same operation and pretending otherwise
#    would be the dishonest version of "comparable".
#
# 3. R1 only, single-end. VASA has one biological read per fragment, so giving
#    FLASH-seq two chances per fragment would bias in its favour. R1 alone keeps
#    the unit of measurement -- one read, one verdict -- identical on both sides.
#
# TRIMMING: why cutadapt directly and not TrimGalore
# --------------------------------------------------
# Trim_Galore/0.6.2 is the only build on this cluster and it is foss-2018b; per
# the repo CLAUDE.md, loading it alongside BWA/SAMtools drags libstdc++ backwards
# and breaks them. This stage needs bwa and samtools loaded, so TrimGalore
# cannot be. Instead we issue the cutadapt call TrimGalore itself issued, taken
# verbatim from the saved report
# (results/trimgalore/*_trimming_report.txt):
#
#     cutadapt -j 8 -e 0.1 -q 20 -O 1 -a CTGTCTCTTATA <mate>
#
# plus TrimGalore's --paired post-filter, which drops a pair when either mate
# falls under 20 bp. In one cutadapt call that is `-A` for the second mate and
# `-m 20:20 --pair-filter=any`. The reproduction is CHECKED, not assumed: the
# script compares its own "Reads with adapter" rate against the rate in the
# saved full-library TrimGalore report and fails if they disagree by more than
# FS_TRIM_TOL percentage points.
#
# OUTPUT
# ------
# Per library, under $FS_RRNA_BWA_DIR/<lib>/ :
#     <lib>.<arm>.ribo-map.log        counts, written by riboread-selection.py
#     <lib>.<arm>.Ribo.bam            the reads it called ribosomal
#     <lib>.<arm>.stranded_y.ribo-map.log   the same BAM re-selected with y
#     <lib>.<arm>.riboloci.tsv        flag/ref/pos per ribosomal read
#     <lib>.<arm>.readlen.tsv         read-length histogram of the arm's input
#     <lib>.cutadapt.log, <lib>.trimcheck.txt
# 05_rrna_bwa_report.py turns those into res/flashseq/rrna_bwa.tsv.
#
# RUNNING
#     sbatch code/flashseq/05_rrna_bwa.sbatch          # all ten, ~40 min
#     FS_LIBS=ZHA8833A1 bash code/flashseq/05_rrna_bwa.sh   # one, for testing
#
# Re-runnable: each library's directory is cleared before it is rebuilt.
# riboread-selection.py ends in a bare `gzip`, which refuses to clobber, so a
# leftover .nonRibo.fastq.gz from a previous run would make the rerun look like
# it worked while writing nothing -- the same trap own_version/pipeline.sh
# guards against in step3.
###############################################################################
set -eu

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$HERE/config.sh"

NCORES="${FS_NCORES:-8}"
WORK="$FS_RRNA_BWA_DIR"
# TrimGalore's own defaults, echoed here so the values are visible at the call.
TRIM_ADAPTER="CTGTCTCTTATA"   # Nextera transposase mosaic end, auto-detected
TRIM_Q=20
TRIM_OVERLAP=1
TRIM_ERR=0.1
TRIM_MINLEN=20
FS_TRIM_TOL="${FS_TRIM_TOL:-3.0}"   # percentage points; see "TRIMMING" above

RIBO_SH="$FS_VASA_SCRIPTS/ribo-bwamem.sh"
RIBO_SEL="$FS_VASA_SCRIPTS/riboread-selection.py"

say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'FATAL: %s\n' "$*" >&2; exit 1; }

###############################################################################
# preflight -- fail in seconds, not an hour in
###############################################################################
for f in "$RIBO_SH" "$RIBO_SEL" "$FS_RRNA_FA" "$FS_CUTADAPT"; do
    [ -e "$f" ] || die "missing $f"
done
[ -x "$RIBO_SH" ]  || die "$RIBO_SH is not executable (chmod +x it)"
[ -x "$RIBO_SEL" ] || die "$RIBO_SEL is not executable (chmod +x it)"
# v1 has no 47S unit and would reproduce exactly the bug this whole line of work
# exists to correct, so refuse to run against it however it was reached.
grep -q '^>mouse_rDNA_47S' "$FS_RRNA_FA" \
    || die "$FS_RRNA_FA has no mouse_rDNA_47S record -- that is v1, not v2"
[ -e "$FS_RRNA_FA.bwt" ] || die "no bwa index beside $FS_RRNA_FA"

eval "$FS_ML_RIBO"
[ -x "$FS_P2BWA/bwa" ]           || die "no bwa at $FS_P2BWA"
[ -x "$FS_P2SAMTOOLS/samtools" ] || die "no samtools at $FS_P2SAMTOOLS"
# riboread-selection.py is `#!/usr/bin/env python3`, so it takes whatever python3
# leads PATH. After FS_ML_RIBO that must be the conda env's, which has pysam.
python3 -c 'import pysam' 2>/dev/null \
    || die "python3 on PATH has no pysam -- FS_ML_RIBO did not activate $FS_VASA_ENV"
say "preflight ok: bwa $("$FS_P2BWA/bwa" 2>&1 | awk '/^Version/{print $2}'), \
samtools $("$FS_P2SAMTOOLS/samtools" --version | head -1 | awk '{print $2}'), \
cutadapt $("$FS_CUTADAPT" --version), python3 $(python3 -V 2>&1 | awk '{print $2}')"

###############################################################################
# libraries, ordered A1..A10 rather than lexically
###############################################################################
if [ -n "${FS_LIBS:-}" ]; then
    LIBS="$FS_LIBS"
else
    LIBS=$(ls "$FS_FASTQ"/ZHA8833A*_R1_001.fastq.gz \
           | sed 's#.*/##; s/_S[0-9]*_.*//' \
           | sort -t A -k3 -n)
fi
[ -n "$LIBS" ] || die "no libraries found under $FS_FASTQ"
say "libraries: $(echo "$LIBS" | tr '\n' ' ')"

###############################################################################
# per library
###############################################################################
for lib in $LIBS; do
    r1=$(ls "$FS_FASTQ/${lib}_S"*"_R1_001.fastq.gz")
    r2=$(ls "$FS_FASTQ/${lib}_S"*"_R2_001.fastq.gz")
    [ -s "$r1" ] && [ -s "$r2" ] || die "$lib: missing R1/R2 under $FS_FASTQ"

    d="$WORK/$lib"
    rm -rf "$d"; mkdir -p "$d"
    say "=== $lib ==="

    # --- sample -------------------------------------------------------------
    # Every STRIDE-th read pair across the WHOLE file, not the first N. Head
    # sampling is measurably biased -- see the FS_STRIDE comment in config.sh --
    # and this is where that bias would enter. R1 and R2 are in the same order,
    # so the same arithmetic on each keeps pairs together; the trim check below
    # would fail if it did not.
    # pipefail is set only here. It is what makes a zcat that dies halfway --
    # the one failure that would silently shrink the sample -- an error rather
    # than a short file; it is NOT set globally because the version probes in
    # the preflight pipe into `head`, which raises SIGPIPE by design.
    for m in 1 2; do
        src=$([ "$m" = 1 ] && echo "$r1" || echo "$r2")
        ( set -o pipefail
          zcat "$src" | awk -v s="$FS_STRIDE" 'int((NR-1)/4)%s==0' \
              | gzip -1 > "$d/$lib.sub_R$m.fq.gz" ) \
            || die "$lib: sampling R$m failed (truncated or unreadable $src)"
    done
    got=$(( $(zcat "$d/$lib.sub_R1.fq.gz" | wc -l) / 4 ))
    got2=$(( $(zcat "$d/$lib.sub_R2.fq.gz" | wc -l) / 4 ))
    [ "$got" -eq "$got2" ] || die "$lib: R1/R2 sampled to $got vs $got2 reads"
    [ "$got" -gt 1000 ] || die "$lib: only $got reads sampled at stride $FS_STRIDE"
    say "  sampled $got read pairs (every ${FS_STRIDE}th)"

    # --- trim ---------------------------------------------------------------
    "$FS_CUTADAPT" -j "$NCORES" -e "$TRIM_ERR" -q "$TRIM_Q" -O "$TRIM_OVERLAP" \
        -a "$TRIM_ADAPTER" -A "$TRIM_ADAPTER" \
        -m "${TRIM_MINLEN}:${TRIM_MINLEN}" --pair-filter=any \
        -o "$d/$lib.trim_R1.fq.gz" -p "$d/$lib.trim_R2.fq.gz" \
        "$d/$lib.sub_R1.fq.gz" "$d/$lib.sub_R2.fq.gz" > "$d/$lib.cutadapt.log"

    # Did we actually reproduce TrimGalore? Compare R1 "Reads with adapters"
    # against the saved full-library report. A mismatch means the reproduction
    # is wrong and every trimmed-arm number below would be wrong with it.
    report="$FS_RESULTS/trimgalore/${lib}_trimmed_1.fastq.gz_trimming_report.txt"
    if [ -s "$report" ]; then
        ref_pct=$(awk -F'[()%]' '/^Reads with adapters:/{print $2; exit}' "$report")
        our_pct=$(awk -F'[()%]' '/^Reads with adapters:/{print $2; exit}' "$d/$lib.cutadapt.log")
        # cutadapt >=2 prints the pair-aware "Read 1 with adapter:" instead
        [ -n "$our_pct" ] || our_pct=$(awk -F'[()%]' '/Read 1 with adapter:/{print $2; exit}' "$d/$lib.cutadapt.log")
        printf 'trimgalore_report_pct\t%s\nour_pct\t%s\ntol\t%s\n' \
            "$ref_pct" "$our_pct" "$FS_TRIM_TOL" > "$d/$lib.trimcheck.txt"
        awk -v a="$ref_pct" -v b="$our_pct" -v t="$FS_TRIM_TOL" \
            'BEGIN{d=a-b; if(d<0)d=-d; exit !(a!="" && b!="" && d<=t)}' \
            || die "$lib: cutadapt reproduction disagrees with the TrimGalore report ($ref_pct% vs $our_pct%, tol $FS_TRIM_TOL) -- see $d/$lib.cutadapt.log"
        say "  trim check: TrimGalore ${ref_pct}% vs ours ${our_pct}% adapter-containing R1"
    else
        say "  trim check SKIPPED: no $report"
    fi

    # --- the two arms -------------------------------------------------------
    for arm in raw trimmed; do
        fq="$d/$lib.$([ "$arm" = raw ] && echo sub || echo trim)_R1.fq.gz"
        out="$d/$lib.$arm"
        n_in=$(( $(zcat "$fq" | wc -l) / 4 ))

        # read-length histogram of exactly what bwa is about to see; this is
        # what answers "is the VASA/FLASH-seq gap just read length?"
        zcat "$fq" | awk 'NR%4==2{h[length($0)]++} END{for(l in h) print l"\t"h[l]}' \
            | sort -n > "$out.readlen.tsv"

        # gzip inside riboread-selection.py will not clobber; clear first.
        rm -f "$out.nonRibo.fastq.gz" "$out.nonRibo.fastq" "$out.Ribo.bam"

        say "  $arm: $n_in reads -> bwa aln + bwa mem"
        "$RIBO_SH" "$FS_RRNA_FA" "$fq" "$out" \
            "$FS_P2BWA" "$FS_P2SAMTOOLS" n "$FS_VASA_SCRIPTS" \
            > "$out.ribo-bwamem.log" 2>&1 \
            || die "$lib/$arm: ribo-bwamem.sh failed -- see $out.ribo-bwamem.log"
        [ -s "$out.ribo-map.log" ] || die "$lib/$arm: no ribo-map.log"

        # The same merged BAM, re-selected with stranded=y. No realignment, so
        # the only difference is the flag -- which is the point.
        rm -f "${out}.stranded_y.nonRibo.fastq.gz" "${out}.stranded_y.nonRibo.fastq"
        "$RIBO_SEL" "$out.nsorted.all-ribo.bam" y "${out}.stranded_y" \
            > "$out.stranded_y.log" 2>&1 \
            || die "$lib/$arm: riboread-selection.py (y) failed"

        # flag / reference / position of every ribosomal read, for the 47S
        # composition and strand tables. Small: ~5% of $n_in rows.
        "$FS_P2SAMTOOLS/samtools" view "$out.Ribo.bam" \
            | awk -v OFS='\t' 'BEGIN{print "flag","ref","pos"} {print $2,$3,$4}' \
            > "$out.riboloci.tsv"

        printf 'reads_in\t%s\n' "$n_in" > "$out.n_in.tsv"

        # nonRibo fastqs are not used downstream here -- FLASH-seq is already
        # mapped by nf-core. Dropped so the work dir stays a few hundred MB.
        rm -f "$out.nonRibo.fastq.gz" "${out}.stranded_y.nonRibo.fastq.gz" \
              "${out}.stranded_y.Ribo.bam"
    done

    rm -f "$d/$lib".sub_R?.fq.gz "$d/$lib".trim_R?.fq.gz
    say "  done $lib"
done

say "all libraries done -> $WORK"
say "now run: fs_python code/flashseq/05_rrna_bwa_report.py"
