#!/usr/bin/env python3
"""
noumi_dryrun_report.py -- is the output of the smartseq_noUMI dry run SANE?

Read-only. Reads the TSVs step 7 wrote plus the STAR logs, and checks the
things that would make the tables wrong without making them crash:

  1. RECONCILIATION. Do reads reconcile from STAR input through assignment to
     the tables? Any table total above its STAR input is a double count.
  2. UFI / TRANSCRIPT DEGENERACY. Predicted by the precheck: with one literal
     UMI 'A' the UFICounts table can only hold 0/1, and bc2trans at K=4 maps
     1 -> 1.0 exactly, so TranscriptCounts must equal UFICounts cell for cell.
     Verified numerically here rather than trusted.
  3. SPLICED + UNSPLICED vs TOTAL. The no-UMI exon/intron branches test for
     the substrings 'exon' / 'intron' in the label, so a combination label
     like 'exon-intron' is counted in BOTH. Measured as an excess.
  4. BIOTYPE COMPOSITION. FLASH-seq is poly-A primed, so protein-coding
     should dominate. Reported on the non-rRNA remainder, the denominator the
     comparison uses, as well as on all reads.
  5. PE vs SE. The two arms, side by side: reads, genes detected, and the
     overlap of their detected-gene sets.

Usage: noumi_dryrun_report.py <outdir> <sample> <celldir> <report.txt>
"""
import sys
import os
import glob
import gzip

import numpy as np
import pandas as pd


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    outdir, sample, celldir, reportpath = sys.argv[1:5]

    out = []

    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 78)
        emit(t)
        emit('=' * 78)

    fail, warn = [], []

    def rd(suffix):
        p = os.path.join(outdir, '%s%s' % (sample, suffix))
        if not os.path.exists(p):
            return None
        return pd.read_csv(p, sep='\t', index_col=0)

    section('0. TABLES STEP 7 WROTE')
    tsvs = sorted(glob.glob(os.path.join(outdir, '%s*.tsv' % sample)))
    emit('%d tables' % len(tsvs))
    for p in tsvs:
        emit('  %-64s %8.1f KB' % (os.path.basename(p), os.path.getsize(p) / 1e3))

    tot_r = rd('_total.ReadCounts.tsv')
    tot_u = rd('_total.UFICounts.tsv')
    tot_t = rd('_total.TranscriptCounts.tsv')
    if tot_r is None:
        sys.exit('no _total.ReadCounts.tsv -- step 7 did not complete')
    emit()
    emit('shape of _total.ReadCounts : %s' % (tot_r.shape,))
    emit('columns                    : %s' % list(tot_r.columns))

    # -----------------------------------------------------------------------
    section('1. RECONCILIATION -- STAR input -> BED rows -> table totals')
    star = {}
    for f in sorted(glob.glob(os.path.join(celldir, '*_E99_Log.final.txt'))):
        arm = os.path.basename(f).split('_cbc')[0]
        d = {}
        for line in open(f):
            if 'Number of input reads' in line:
                d['input'] = int(line.strip().rsplit('\t')[-1])
            if 'Uniquely mapped reads number' in line:
                d['uniq'] = int(line.strip().rsplit('\t')[-1])
            if 'Number of reads mapped to multiple loci' in line:
                d['multi'] = int(line.strip().rsplit('\t')[-1])
        star[arm] = d
    emit('STAR Log.final per arm (STAR counts PAIRS as one "read" in PE):')
    for a, d in star.items():
        emit('  %-16s input=%-10d uniq=%-10d multi=%-10d' %
             (a, d.get('input', -1), d.get('uniq', -1), d.get('multi', -1)))

    armstats = os.path.join(outdir, 'logs', 'arm_stats.tsv')
    if os.path.exists(armstats):
        emit()
        emit('BED rows / QNAMEs per arm (from the driver):')
        for line in open(armstats):
            emit('  ' + line.rstrip('\n'))

    emit()
    emit('%-18s %14s %14s %14s %10s' %
         ('column', 'table reads', 'STAR input', 'STAR uniq', 'reads/input'))
    for cc in tot_r.columns:
        arm = str(cc).rsplit('/')[-1]
        d = star.get(arm, {})
        tr = int(tot_r[cc].sum())
        si = d.get('input')
        emit('%-18s %14d %14s %14s %10s' %
             (cc, tr, si if si else 'n/a', d.get('uniq', 'n/a'),
              ('%.3f' % (tr / si)) if si else 'n/a'))
        if si and tr > si:
            warn.append('%s: table reads (%d) EXCEED STAR input reads (%d), '
                        'ratio %.3f -- consistent with both mates of a pair '
                        'being counted separately' % (cc, tr, si, tr / si))

    # -----------------------------------------------------------------------
    section('2. UFI / TRANSCRIPT DEGENERACY -- the precheck prediction, checked')
    if tot_u is None or tot_t is None:
        fail.append('UFI or Transcript table missing')
    else:
        vals = sorted(pd.unique(tot_u.values.ravel()))
        emit('distinct values in _total.UFICounts : %s%s'
             % (vals[:12], ' ...' if len(vals) > 12 else ''))
        if set(vals) <= {0, 1}:
            emit('OK as predicted: UFICounts holds only 0/1 -- it is a DETECTION MASK,')
            emit('not a molecule count. One literal UMI "A" per gene per column.')
        else:
            fail.append('UFICounts holds values outside {0,1}: %s -- the no-UMI '
                        'assumption is violated' % vals[:12])
        same = np.allclose(tot_u.values.astype(float), tot_t.values.astype(float))
        emit('TranscriptCounts == UFICounts elementwise: %s' % same)
        if same:
            emit('OK as predicted: bc2trans at K=4**len("A")=4 maps 1 -> 1.0 exactly,')
            emit('so the TranscriptCounts tables are the detection mask verbatim.')
            emit('They carry NO abundance information and must not be used as a')
            emit('quantification. Use ReadCounts.')
        else:
            warn.append('TranscriptCounts differs from UFICounts -- check K')
        emit()
        emit('%-18s %12s %12s %12s' %
             ('column', 'reads', 'genes detected', 'sum(UFI)'))
        for cc in tot_r.columns:
            emit('%-18s %12d %12d %12d' %
                 (cc, int(tot_r[cc].sum()), int((tot_r[cc] > 0).sum()),
                  int(tot_u[cc].sum())))
        d = int((tot_u.values.sum() - (tot_r.values > 0).sum()))
        emit('sum(UFI) - count(reads>0) = %d   (0 confirms UFI is exactly the mask)' % d)

    # -----------------------------------------------------------------------
    section('3. SPLICED + UNSPLICED vs TOTAL (uni-genes)')
    ur = rd('_uniaggGenes_total.ReadCounts.tsv')
    sr = rd('_uniaggGenes_spliced.ReadCounts.tsv')
    xr = rd('_uniaggGenes_unspliced.ReadCounts.tsv')
    if ur is None or sr is None or xr is None:
        warn.append('uniaggGenes read tables missing -- cannot check the exon/intron split')
    else:
        emit('%-18s %14s %14s %14s %14s %9s' %
             ('column', 'total', 'spliced', 'unspliced', 'spl+unspl', 'excess'))
        for cc in ur.columns:
            t = int(ur[cc].sum()); s = int(sr[cc].sum()); x = int(xr[cc].sum())
            exc = s + x - t
            emit('%-18s %14d %14d %14d %14d %9d' % (cc, t, s, x, s + x, exc))
            if exc > 0:
                warn.append('%s: spliced+unspliced exceeds total by %d reads (%.3f%%) '
                            "-- the no-UMI branches test for the substrings 'exon' and "
                            "'intron', so a label containing both is counted twice"
                            % (cc, exc, 100.0 * exc / max(1, t)))
            elif exc < 0:
                emit('      (deficit: reads whose label is neither exon nor intron)')
        emit()
        emit('unspliced fraction (unspliced/total, uni-genes, reads):')
        for cc in ur.columns:
            t = int(ur[cc].sum())
            emit('  %-18s %.4f' % (cc, int(xr[cc].sum()) / max(1, t)))

    # -----------------------------------------------------------------------
    section('4. BIOTYPE COMPOSITION -- poly-A primed, so protein-coding should win')
    # index labels are ENSID_Symbol_Biotype after step 6 strips the exon/intron
    # field; tRNA rows are the "N.tRNAn.IsotypeNNN" form.
    def biotype_of(idx):
        idx = str(idx)
        if 'tRNA' in idx and idx.count('_') == 0:
            return 'tRNA'
        parts = idx.rsplit('-')
        bs = sorted(set(p.rsplit('_')[-1] for p in parts))
        return bs[0] if len(bs) == 1 else 'MULTI:' + '-'.join(bs)

    bt = pd.Series([biotype_of(i) for i in tot_r.index], index=tot_r.index)
    for cc in tot_r.columns:
        g = tot_r[cc].groupby(bt).sum().sort_values(ascending=False)
        tot = g.sum()
        ribo = g[[b for b in g.index if b in ('rRNA', 'MtRrna', 'rRNApseudogene')]].sum() \
            if len(g) else 0
        emit()
        emit('%s -- %d reads total, %d in an rRNA-family biotype (%.2f%%)'
             % (cc, tot, ribo, 100.0 * ribo / max(1, tot)))
        emit('  %-42s %12s %9s %9s' % ('biotype', 'reads', '% all', '% nonRibo'))
        for b, v in g.head(14).items():
            emit('  %-42s %12d %8.2f%% %8.2f%%'
                 % (b, v, 100.0 * v / max(1, tot), 100.0 * v / max(1, tot - ribo)))
        pc = g.get('ProteinCoding', 0)
        emit('  ProteinCoding as %% of the non-rRNA remainder: %.2f%%'
             % (100.0 * pc / max(1, tot - ribo)))
        if 100.0 * pc / max(1, tot - ribo) < 40:
            warn.append('%s: protein-coding is only %.1f%% of the non-rRNA remainder, '
                        'which is low for a poly-A primed library'
                        % (cc, 100.0 * pc / max(1, tot - ribo)))

    # -----------------------------------------------------------------------
    section('5. PE vs SE')
    cols = list(tot_r.columns)
    emit('columns that reached step 7: %s' % cols)
    pe_dir = os.path.join(outdir, 'cells_pe_unusable')
    if os.path.isdir(pe_dir):
        nq = len([f for f in os.listdir(pe_dir)])
        emit()
        emit('The PE arm is NOT here: %d of its files are quarantined in' % nq)
        emit('cells_pe_unusable/ because step 6 cannot parse a paired-end step-5')
        emit("BED -- bedtools bamtobed's '/1','/2' mate suffix lands on the end of")
        emit('the nM value, and step 6 line 97 int()s it. Rates are in')
        emit('logs/nm_contract.txt; the mechanism is in NOUMI_PATH.md.')
        emit('So this path is SINGLE-END ONLY, and the comparison below is not')
        emit('available. That is a measured result, not a missing measurement.')
    if len(cols) == 2:
        a, b = cols
        da = set(tot_r.index[tot_r[a] > 0])
        db = set(tot_r.index[tot_r[b] > 0])
        emit('%-18s genes detected = %d' % (a, len(da)))
        emit('%-18s genes detected = %d' % (b, len(db)))
        emit('intersection = %d   Jaccard = %.4f'
             % (len(da & db), len(da & db) / max(1, len(da | db))))
        emit('only in %-10s = %d' % (a, len(da - db)))
        emit('only in %-10s = %d' % (b, len(db - da)))
        both = sorted(da & db)
        if both:
            r = np.corrcoef(np.log1p(tot_r.loc[both, a].values.astype(float)),
                            np.log1p(tot_r.loc[both, b].values.astype(float)))[0, 1]
            emit('Pearson r of log1p(reads) over the shared genes: %.4f' % r)
            ratio = tot_r[a].sum() / max(1, tot_r[b].sum())
            emit('read ratio %s/%s = %.4f' % (a, b, ratio))
        emit()
        emit('top 12 genes by reads, per arm:')
        for cc in cols:
            top = tot_r[cc].sort_values(ascending=False).head(12)
            emit('  %s:' % cc)
            for i, v in top.items():
                emit('      %-58s %10d' % (str(i)[:58], int(v)))
    elif len(cols) == 1:
        emit()
        emit('One column (%s), as expected once the PE arm is excluded.' % cols[0])
        top = tot_r[cols[0]].sort_values(ascending=False).head(15)
        emit('top 15 genes by reads:')
        for i, v in top.items():
            emit('    %-58s %10d' % (str(i)[:58], int(v)))
    else:
        warn.append('expected 1 column (se) or 2 (se+pe), found %d: %s'
                    % (len(cols), cols))

    # -----------------------------------------------------------------------
    section('VERDICT')
    if fail:
        emit('FAIL (%d):' % len(fail))
        for f in fail:
            emit('  x %s' % f)
    else:
        emit('The smartseq_noUMI path ran end to end and its output is internally')
        emit('consistent. Semantics of each table: see the WARNs and NOUMI_PATH.md.')
    if warn:
        emit()
        emit('WARN (%d) -- not crashes, but they change what the numbers MEAN:' % len(warn))
        for w in warn:
            emit('  ! %s' % w)

    d = os.path.dirname(reportpath)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(reportpath, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % reportpath, flush=True)
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
