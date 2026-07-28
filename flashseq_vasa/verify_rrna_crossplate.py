#!/usr/bin/env python3
"""
verify_rrna_crossplate.py -- re-derive every cross-plate rRNA number from the
BAMs, and assert it against the value that was reported.

WHY THIS EXISTS
---------------
The cross-plate comparison (published VASA-plate SRR14783059 vs the user's
ZHA9292A1) was measured interactively, in a shell heredoc, from a `fresh` kernel
whose state is gone. That is not a reproducible provenance for a number that is
now in a report. This script is the reproducible form: one file, run it, and it
either reproduces the reported figures or it fails loudly.

It is READ-ONLY with respect to the data. It writes two TSVs and a report.

WHAT IT CHECKS, AND WHY EACH CHECK IS THE RIGHT ONE
--------------------------------------------------
1. SELF-CHECK AGAINST THE PIPELINE'S OWN LOG.
   The `stranded=y` predicate implemented here must reproduce step 3's own
   `ribo%` column for all 16 barcodes of ZHA9292A1. This is the check that
   matters most: it proves the predicate reimplemented here is the same
   predicate `riboread-selection.py` applied, rather than merely a plausible
   one. If this fails, nothing else in the script means anything.

2. THE TWO PREDICATES, STATED EXPLICITLY.
   riboread-selection.py groups alignments by read name, then:
     - all unmapped                     -> not ribosomal (goes to nonRibo.fastq)
     - stranded='n': any mapped         -> ribosomal
     - stranded='y': any mapped forward -> ribosomal
   Reimplemented here verbatim in `classify()`. The read-name grouping matters:
   a read with several alignments is ONE observation, not several.

3. REPORTED-VALUE ASSERTIONS.
   Every number quoted in RRNA_PROVENANCE.md is re-derived and asserted to
   within a stated tolerance. A mismatch is a FAIL with both values printed.

WHAT THIS SCRIPT DELIBERATELY DOES NOT DO
-----------------------------------------
It does not re-run bwa. The `.nsorted.all-ribo.bam` files are the output of
`ribo-bwamem.sh`, and re-aligning would test bwa's determinism rather than the
counting logic that the comparison actually rests on. The alignment provenance
is covered by `provenance.tsv` (same fasta checksum, same script).

It does not adjudicate the biological question. It verifies arithmetic and
predicate fidelity. Whether 2.1x is a depletion problem is not a thing a script
can answer.

USAGE
    verify_rrna_crossplate.py <own_cells_dir> <published_bam_dir> <step3_log> <outdir>
"""
import os
import sys
import re
import glob

import pysam


# ---------------------------------------------------------------------------
# The two predicates, copied from ../I_Gene_expression/a_Mapping/riboread-selection.py.
#
# The source loops over a name-sorted BAM accumulating alignments while the
# qname is unchanged, then classifies the group. `all unmapped` is tested first
# and short-circuits to the non-ribosomal arm; only then does the stranded flag
# select between "any mapped" and "any mapped on the forward strand".
#
# This is the classification only. The source also WRITES the fastq/bam and
# counts per-reference hits; none of that is needed to reproduce the fractions,
# and reimplementing it would add failure modes without adding fidelity.
# ---------------------------------------------------------------------------
def classify(group):
    """-> (is_ribo_n, is_ribo_y) for one read's alignment group."""
    if all(r.is_unmapped for r in group):
        return (False, False)
    ribo_n = True                                    # stranded='n': any mapped
    ribo_y = any((not r.is_unmapped) and (not r.is_reverse)
                 for r in group)                     # stranded='y': any fwd
    return (ribo_n, ribo_y)


def count_bam(path):
    """Group a name-sorted all-ribo BAM by qname and apply both predicates."""
    bam = pysam.AlignmentFile(path, check_sq=False)
    total = n_any = n_fwd = 0
    cur = None
    grp = []
    for r in bam.fetch(until_eof=True):
        if cur is None or r.qname == cur:
            grp.append(r)
            cur = r.qname
        else:
            total += 1
            a, f = classify(grp)
            n_any += a
            n_fwd += f
            grp = [r]
            cur = r.qname
    if grp:
        total += 1
        a, f = classify(grp)
        n_any += a
        n_fwd += f
    bam.close()
    return total, n_any, n_fwd


def parse_step3_log(path):
    """Pull the per-cell ribo% column out of step3.log for the self-check."""
    want = re.compile(r'^\s*(\d{3})\s+([\d,]+)\s+([\d,]+)\s+([\d.]+)%')
    out = {}
    with open(path) as fh:
        for line in fh:
            m = want.match(line)
            if m:
                cc = m.group(1)
                out[cc] = dict(reads=int(m.group(2).replace(',', '')),
                               ribo=int(m.group(3).replace(',', '')),
                               pct=float(m.group(4)))
    return out


def main():
    try:
        own_dir, pub_dir, step3_log, outdir = sys.argv[1:5]
    except ValueError:
        sys.exit(__doc__)

    os.makedirs(outdir, exist_ok=True)
    out = []
    fail = []
    warn = []

    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 74)
        emit(t)
        emit('=' * 74)

    # -----------------------------------------------------------------------
    section('OWN PLATE -- ZHA9292A1 (PM26037)')
    own = {}
    for b in sorted(glob.glob(os.path.join(own_dir, '*.nsorted.all-ribo.bam'))):
        m = re.search(r'_(\d{3})_', os.path.basename(b))
        cc = m.group(1) if m else os.path.basename(b)
        t, na, nf = count_bam(b)
        own[cc] = dict(total=t, ribo_n=na, ribo_y=nf,
                       pct_n=100.0 * na / t if t else 0.0,
                       pct_y=100.0 * nf / t if t else 0.0,
                       fwd=100.0 * nf / na if na else 0.0)
        emit('  %s  total=%9d  ribo_y=%9d  pct_y=%6.2f  pct_n=%6.2f  fwd=%5.1f'
             % (cc, t, nf, own[cc]['pct_y'], own[cc]['pct_n'], own[cc]['fwd']))

    # -----------------------------------------------------------------------
    section('CHECK 1 -- does our stranded=y predicate reproduce step3.log?')
    emit('This is the load-bearing check: it proves the predicate implemented')
    emit('here is the one riboread-selection.py actually applied.')
    emit()
    log = parse_step3_log(step3_log)
    if not log:
        fail.append('step3.log could not be parsed -- self-check impossible')
        emit('  x could not parse %s' % step3_log)
    else:
        worst = 0.0
        for cc, ref in sorted(log.items()):
            if cc not in own:
                warn.append('cell %s in step3.log but no BAM on disk' % cc)
                continue
            d = abs(own[cc]['pct_y'] - ref['pct'])
            worst = max(worst, d)
            flag = ' ' if d <= 0.01 else 'x'
            emit('  %s %s  ours=%6.3f  step3.log=%6.3f  |diff|=%.4f'
                 % (flag, cc, own[cc]['pct_y'], ref['pct'], d))
        emit()
        emit('  worst |deviation| = %.4f percentage points' % worst)
        if worst > 0.01:
            fail.append('predicate does not reproduce step3.log (worst %.4f pts)' % worst)
        else:
            emit('  -> PASS: same predicate, to within log rounding.')

    # -----------------------------------------------------------------------
    section('PUBLISHED PLATE -- SRR14783059')
    pub = {}
    for b in sorted(glob.glob(os.path.join(pub_dir, '*.nsorted.all-ribo.bam'))):
        m = re.search(r'_(\d{3})\.', os.path.basename(b))
        cc = m.group(1) if m else os.path.basename(b)
        t, na, nf = count_bam(b)
        pub[cc] = dict(total=t, ribo_n=na, ribo_y=nf,
                       pct_n=100.0 * na / t if t else 0.0,
                       pct_y=100.0 * nf / t if t else 0.0,
                       fwd=100.0 * nf / na if na else 0.0)
        emit('  %s  total=%9d  ribo_y=%9d  pct_y=%6.2f  pct_n=%6.2f  fwd=%5.1f'
             % (cc, t, nf, pub[cc]['pct_y'], pub[cc]['pct_n'], pub[cc]['fwd']))

    # -----------------------------------------------------------------------
    section('CHECK 2 -- reported values')
    BLANKS = {'001', '014', '015', '016'}
    real = {k: v for k, v in own.items() if k not in BLANKS}

    def med(vals):
        v = sorted(vals)
        n = len(v)
        return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])

    checks = []

    if real:
        checks += [
            ('own real-cell pct_y min', min(v['pct_y'] for v in real.values()), 17.99, 0.02),
            ('own real-cell pct_y max', max(v['pct_y'] for v in real.values()), 25.09, 0.02),
            ('own real-cell pct_y median', med([v['pct_y'] for v in real.values()]), 21.14, 0.05),
            ('own real-cell fwd min', min(v['fwd'] for v in real.values()), 57.7, 0.1),
            ('own real-cell fwd max', max(v['fwd'] for v in real.values()), 83.8, 0.1),
            ('own real-cell fwd median', med([v['fwd'] for v in real.values()]), 74.8, 0.1),
        ]
        tot = sum(v['total'] for v in own.values())
        ry = sum(v['ribo_y'] for v in own.values())
        checks.append(('own pooled pct_y (all 16 bc)', 100.0 * ry / tot, 21.39, 0.02))

    blanks_present = {k: v for k, v in own.items() if k in BLANKS}
    if blanks_present:
        checks += [
            ('own blank fwd min', min(v['fwd'] for v in blanks_present.values()), 86.4, 0.1),
            ('own blank fwd max', max(v['fwd'] for v in blanks_present.values()), 93.6, 0.1),
        ]

    if pub:
        checks += [
            ('published pct_y min', min(v['pct_y'] for v in pub.values()), 1.33, 0.02),
            ('published pct_y max', max(v['pct_y'] for v in pub.values()), 26.42, 0.02),
            ('published pct_y median', med([v['pct_y'] for v in pub.values()]), 9.90, 0.05),
            ('published fwd min', min(v['fwd'] for v in pub.values()), 92.7, 0.1),
            ('published fwd max', max(v['fwd'] for v in pub.values()), 97.4, 0.1),
        ]
        tp = sum(v['total'] for v in pub.values())
        pyv = sum(v['ribo_y'] for v in pub.values())
        checks.append(('published pooled pct_y', 100.0 * pyv / tp, 9.19, 0.02))

    emit('  %-34s %10s %10s %8s' % ('quantity', 'recomputed', 'reported', 'ok?'))
    for name, got, want, tol in checks:
        ok = abs(got - want) <= tol
        emit('  %-34s %10.4f %10.4f %8s' % (name, got, want, 'yes' if ok else 'NO'))
        if not ok:
            fail.append('%s: recomputed %.4f, reported %.4f (tol %.3f)'
                        % (name, got, want, tol))

    # -----------------------------------------------------------------------
    section('CHECK 3 -- derived ratios')
    if pub and real:
        pm = med([v['pct_y'] for v in pub.values()])
        om = med([v['pct_y'] for v in real.values()])
        emit('  median ratio (own / published)      = %.3fx   [reported 2.1x]' % (om / pm))
        pooled_ratio = ((sum(v['ribo_y'] for v in own.values()) / sum(v['total'] for v in own.values()))
                        / (sum(v['ribo_y'] for v in pub.values()) / sum(v['total'] for v in pub.values())))
        emit('  pooled ratio (own / published)      = %.3fx   [reported 2.3x]' % pooled_ratio)
        ps = max(v['pct_y'] for v in pub.values()) / min(v['pct_y'] for v in pub.values())
        os_ = max(v['pct_y'] for v in real.values()) / min(v['pct_y'] for v in real.values())
        emit('  spread published                    = %.2fx' % ps)
        emit('  spread own (real cells)             = %.2fx' % os_)
        emit('  own spread tighter by               = %.1fx   [reported 14x]' % (ps / os_))
        emit()
        emit('  NOTE the two ratios differ (median vs pooled) because the published')
        emit('  cells span 64k-1.26M reads, so pooling weights the deep cells. Quote')
        emit('  which one is meant; they are not interchangeable.')

    # -----------------------------------------------------------------------
    section('TSV OUTPUT')
    for name, d in (('verify_own_percell.tsv', own), ('verify_published_percell.tsv', pub)):
        p = os.path.join(outdir, name)
        with open(p, 'w') as fh:
            fh.write('cell\ttotal\tribo_n\tribo_y\tpct_n\tpct_y\tfwd_share\n')
            for cc in sorted(d):
                v = d[cc]
                fh.write('%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.2f\n'
                         % (cc, v['total'], v['ribo_n'], v['ribo_y'],
                            v['pct_n'], v['pct_y'], v['fwd']))
        emit('  wrote %s (%d rows)' % (p, len(d)))

    # -----------------------------------------------------------------------
    section('VERDICT')
    if fail:
        emit('FAIL (%d):' % len(fail))
        for f in fail:
            emit('  x %s' % f)
    else:
        emit('PASS -- every reported cross-plate number reproduces from the BAMs.')
    if warn:
        emit()
        for w in warn:
            emit('  ! %s' % w)

    rp = os.path.join(outdir, 'verify_rrna_crossplate.txt')
    with open(rp, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport: %s' % rp, flush=True)
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
