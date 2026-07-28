#!/usr/bin/env python3
"""
measure_vasa_readlen.py -- measure the read-length distribution VASA actually
gave STAR, and turn it into a lookup table the FLASH-seq vasalen arm draws from.

WHY THIS EXISTS
---------------
The comparison's second axis is "which RNA species are detected", and the short
non-poly-A species (tRNA, snoRNA, snRNA, miRNA) are exactly where VASA is
expected to win. But read length alone suppresses them on the FLASH-seq side,
independently of any biology:

  * step 5 tags a read jS:IN only if the read is CONTAINED in the feature
    (readstart >= refstart && readend <= refend);
  * step 6 keeps a non-spliceable biotype ONLY when jS == IN;
  * in the v2 BED, 98.5% of tRNA features, 99.2% of miRNA, 96.5% of snoRNA and
    84.6% of snRNA are shorter than one 151 nt read.

So a 151 nt read physically cannot be contained in almost any of them. In the A9
dry run all six tRNA-overlapping reads were discarded by containment and the
tRNA table came out empty, while VASA cell 005 had 229 tRNA rows.

The control the user chose is to hard-trim FLASH-seq R1 to VASA's own length
distribution before mapping, mirroring what code/flashseq/05_rrna_bwa.sh did for
the rRNA leg -- one measurement, both sides. This script measures the target.

WHAT IS MEASURED, AND WHY THAT QUANTITY
---------------------------------------
Two different lengths are involved and they are NOT interchangeable:

  (a) STAR-INPUT READ LENGTH -- the length of the sequence in the fastq that
      STAR is handed. This is the ONLY one we can control by trimming, so it is
      what the lookup table is built from. Source: the step-3 output
      *_cbc_trimmed_homoATCG.nonRibo.fastq.gz, i.e. the exact files step 4 read.

  (b) ALIGNED SPAN -- (End - Start) of the step-5 BED interval, which is what
      the containment test actually compares against a feature. It is a
      consequence of (a) plus soft-clipping and splicing, not something a
      trimmer can set.

Matching (a) is the honest, controllable intervention; (b) is measured on both
sides afterwards and reported, so how well the containment-relevant quantity
ended up matching is visible rather than assumed. Both are written here.

WHY A DISTRIBUTION AND NOT A SINGLE NUMBER
------------------------------------------
Trimming to VASA's MEDIAN would not be a control. VASA's median STAR-input read
still exceeds nearly every tRNA/miRNA/snoRNA feature, so a median-trimmed arm
would recover nothing and would prove nothing. What lets VASA satisfy jS:IN on a
short feature is the SHORT TAIL of its distribution. Reproducing only the median
discards exactly the part of the distribution that carries the effect. So the
whole distribution is reproduced.

WHICH CELLS
-----------
The 12 real cells only. The four confirmed blanks (001, 014, 015, 016) are not
biology -- their STAR average mapped lengths are 78.3, 57.3, 71.9 and 51.6 nt
against 108-115 for the real cells, so including them would drag the target
distribution short for the wrong reason. They are measured and reported anyway,
so the exclusion is auditable rather than silent.

OUTPUT
------
  <lut>            10,000 integers, one per line. Draw a uniform index and read
                   off a length: that reproduces the pooled distribution to
                   0.01% resolution. Consumed by trim_to_vasalen.sh.
  <report>         the full histograms, percentiles, the per-cell breakdown, the
                   blanks for comparison, and the aligned-span distribution.
  <lut>.hist.tsv   length<TAB>count, the full pooled STAR-input histogram.
  <lut>.span.tsv   the aligned-span histogram (BED End-Start), stride-sampled.

Usage:
  measure_vasa_readlen.py <celldir> <sample> <realcells> <blanks> <lut> <report>
                          [bed_stride]

  <realcells>, <blanks>  space-separated cell ids, e.g. "002 003 ..."
  [bed_stride]           sample every Nth BED row for the aligned-span
                         histogram (default 64; the BEDs are ~2.6 GB gz)
"""
import os
import subprocess
import sys
from collections import Counter

LUT_N = 10000


def die(msg):
    sys.exit('FATAL: %s' % msg)


def fastq_len_hist(path):
    """Histogram of sequence lengths over the WHOLE fastq.

    Shelled out to awk rather than looped in python: this is ~70 M reads across
    the 12 real cells and awk does it in a few minutes where python would take
    an hour. Note NR%4==2 -- the sequence line. Not a sample: the whole file, so
    there is no stride bias to argue about.
    """
    cmd = "zcat %s | awk 'NR%%4==2{h[length($0)]++} END{for(l in h) print l\"\\t\"h[l]}'" % path
    p = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd],
                       capture_output=True, text=True)
    if p.returncode != 0:
        die('reading %s failed: %s' % (path, p.stderr.strip()[:400]))
    h = Counter()
    for line in p.stdout.splitlines():
        if not line.strip():
            continue
        a, b = line.split('\t')
        h[int(a)] += int(b)
    return h


def bed_span_hist(path, stride):
    """Histogram of BED interval widths (End - Start), every Nth row.

    A stride, not a head: the BED is coordinate-sorted, so the head of the file
    is chromosome 1 and its length composition is not the library's.
    """
    cmd = ("zcat %s | awk -v s=%d 'NR%%s==1{h[$3-$2]++} "
           "END{for(l in h) print l\"\\t\"h[l]}'" % (path, stride))
    p = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return None
    h = Counter()
    for line in p.stdout.splitlines():
        if not line.strip():
            continue
        a, b = line.split('\t')
        h[int(a)] += int(b)
    return h


def pct(h, q):
    """Percentile of a histogram, by cumulative count."""
    n = sum(h.values())
    if n == 0:
        return None
    target = q / 100.0 * n
    c = 0
    for k in sorted(h):
        c += h[k]
        if c >= target:
            return k
    return max(h)


def describe(h):
    n = sum(h.values())
    if n == 0:
        return 'empty'
    mean = sum(k * v for k, v in h.items()) / float(n)
    return ('n=%d  mean=%.2f  min=%d  p1=%s p5=%s p10=%s p25=%s p50=%s '
            'p75=%s p90=%s p95=%s p99=%s  max=%d'
            % (n, mean, min(h), pct(h, 1), pct(h, 5), pct(h, 10), pct(h, 25),
               pct(h, 50), pct(h, 75), pct(h, 90), pct(h, 95), pct(h, 99),
               max(h)))


def build_lut(h, n_slots):
    """Expand a histogram into n_slots length values, proportional to frequency.

    Largest-remainder allocation: every length gets floor(share) slots, and the
    leftover slots go to the largest fractional remainders. That keeps rare short
    lengths -- the ones that carry the whole point of this arm -- from being
    rounded out of existence wholesale, while never inventing a length that was
    not observed.
    """
    total = float(sum(h.values()))
    keys = sorted(h)
    exact = [(k, h[k] / total * n_slots) for k in keys]
    slots = {k: int(v) for k, v in exact}
    used = sum(slots.values())
    rema = sorted(((v - int(v)), k) for k, v in exact)
    for _, k in reversed(rema):
        if used >= n_slots:
            break
        slots[k] += 1
        used += 1
    lut = []
    for k in keys:
        lut.extend([k] * slots[k])
    return lut


def main():
    if len(sys.argv) < 7:
        sys.exit(__doc__)
    celldir, sample, realcells_s, blanks_s, lutpath, reportpath = sys.argv[1:7]
    bed_stride = int(sys.argv[7]) if len(sys.argv) > 7 else 64
    realcells = realcells_s.split()
    blanks = blanks_s.split()

    out = []

    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 78)
        emit(t)
        emit('=' * 78)

    section('WHAT IS BEING MEASURED')
    emit('target quantity : STAR-input read length (the fastq step 4 read)')
    emit('source files    : %s/<cell>_cbc_trimmed_homoATCG.nonRibo.fastq.gz' % celldir)
    emit('sample          : %s' % sample)
    emit('real cells (pooled into the LUT) : %s' % ' '.join(realcells))
    emit('blanks (measured, NOT pooled)    : %s' % ' '.join(blanks))
    emit('whole files, not a sample -- no stride bias in the target distribution')

    section('PER-CELL STAR-INPUT READ LENGTHS')
    pooled = Counter()
    percell = {}
    for grp, cells in (('real', realcells), ('blank', blanks)):
        for cell in cells:
            fq = os.path.join(celldir, '%s_%s_cbc_trimmed_homoATCG.nonRibo.fastq.gz'
                              % (sample, cell))
            if not os.path.exists(fq):
                die('missing %s' % fq)
            h = fastq_len_hist(fq)
            percell[cell] = (grp, h)
            if grp == 'real':
                pooled += h
            emit('%-6s %-5s %s' % (grp, cell, describe(h)))

    section('POOLED TARGET DISTRIBUTION (12 real cells)')
    emit(describe(pooled))
    n_pool = sum(pooled.values())
    emit()
    emit('This is the distribution the vasalen arm reproduces. Note the SHORT')
    emit('TAIL: it is what lets a VASA read be CONTAINED in a sub-151 nt feature')
    emit('and satisfy jS:IN, and it is the reason a single fixed trim length')
    emit('would not be a control.')
    emit()
    emit('fraction of the pooled distribution at or below each length:')
    cum = 0
    for thr in (20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 151):
        cum = sum(v for k, v in pooled.items() if k <= thr)
        emit('  <= %3d nt : %10d  %6.2f%%' % (thr, cum, 100.0 * cum / n_pool))

    section('BLANKS, FOR COMPARISON -- reported, deliberately excluded')
    bl = Counter()
    for cell in blanks:
        bl += percell[cell][1]
    if sum(bl.values()):
        emit('pooled blanks: %s' % describe(bl))
        emit()
        emit('The blanks are shorter than the real cells (their STAR average')
        emit('mapped lengths are 51.6-78.3 nt against 108-115 for the real')
        emit('cells). Pooling them would drag the target distribution short for')
        emit('a reason that has nothing to do with the biology being compared,')
        emit('which is why they are out.')

    section('ALIGNED SPAN (BED End-Start) -- the containment-relevant quantity')
    emit('This is NOT what the LUT is built from -- it cannot be set by a')
    emit('trimmer. It is measured here so that after the vasalen arm is mapped,')
    emit('how well the containment-relevant quantity actually matched can be')
    emit('checked rather than assumed. Stride %d over the step-5 singlemapper BEDs.'
         % bed_stride)
    span = Counter()
    for cell in realcells:
        bed = os.path.join(
            celldir,
            '%s_%s_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.singlemappers_genes.bed.gz'
            % (sample, cell))
        if not os.path.exists(bed):
            emit('  (missing %s -- skipped)' % os.path.basename(bed))
            continue
        h = bed_span_hist(bed, bed_stride)
        if h is None:
            emit('  (unreadable %s -- skipped)' % os.path.basename(bed))
            continue
        span += h
        emit('  %-5s %s' % (cell, describe(h)))
    if sum(span.values()):
        emit()
        emit('pooled real cells: %s' % describe(span))
        n_span = sum(span.values())
        emit('rows spanning >= 140 nt : %d (%.2f%%)'
             % (sum(v for k, v in span.items() if k >= 140),
                100.0 * sum(v for k, v in span.items() if k >= 140) / n_span))

    section('LOOKUP TABLE')
    lut = build_lut(pooled, LUT_N)
    if len(lut) != LUT_N:
        die('LUT has %d entries, expected %d' % (len(lut), LUT_N))
    lo = min(lut)
    emit('%d entries, min=%d max=%d' % (len(lut), lo, max(lut)))
    emit('distinct lengths represented: %d of %d observed'
         % (len(set(lut)), len(pooled)))
    # VASA's own step 2 floors reads at TRIM_MINLEN=20, so nothing below 20 is
    # expected. Assert rather than trust: a shorter entry would mean the LUT can
    # ask the trimmer for a read length VASA never produced.
    if lo < 20:
        die('LUT contains a length below 20 nt (%d) -- VASA step 2 floors at 20, '
            'so this means the source fastqs are not what this script assumes' % lo)
    # How faithfully does the 10,000-slot LUT reproduce the pooled distribution?
    lh = Counter(lut)
    worst = max(abs(lh[k] / float(LUT_N) - pooled[k] / float(n_pool))
                for k in pooled)
    emit('largest per-length probability error from the %d-slot quantisation: %.5f'
         % (LUT_N, worst))
    emit('LUT percentiles: %s' % describe(lh))

    with open(lutpath, 'w') as fh:
        fh.write('\n'.join(str(x) for x in lut) + '\n')
    with open(lutpath + '.hist.tsv', 'w') as fh:
        fh.write('length\tcount\n')
        for k in sorted(pooled):
            fh.write('%d\t%d\n' % (k, pooled[k]))
    if sum(span.values()):
        with open(lutpath + '.span.tsv', 'w') as fh:
            fh.write('span\tcount_sampled\tstride\n')
            for k in sorted(span):
                fh.write('%d\t%d\t%d\n' % (k, span[k], bed_stride))
    emit()
    emit('written: %s' % lutpath)
    emit('written: %s.hist.tsv' % lutpath)

    d = os.path.dirname(reportpath)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(reportpath, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % reportpath, flush=True)


if __name__ == '__main__':
    main()
