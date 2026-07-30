#!/usr/bin/env python3
"""verify_coverage_threeway.py -- re-derive and ASSERT every number claimed in
the three-way coverage writeup, from the tables on disk.

This project's rule is that a number entering a report is verified by a
committed script rather than transcribed from a log. Each check below states
the claim, recomputes it from `res/threeway/coverage_threeway*.tsv`, and
asserts. A failure here means the writeup is wrong, not that the check is
mis-tuned -- so nothing is given a tolerance wider than the printing precision
of the value it tests.

Two things are deliberately NOT re-derived:
  * the base-vs-mid discrepancy's CAUSE. Edge clipping is the leading untested
    hypothesis; the withdrawn refutation is not resurrected. The loss
    decomposition is asserted to be internally exact, which is a statement
    about arithmetic, not mechanism.
  * qualimap's 3'/5' = 0.93 on the FLASH-seq BAMs. That is nf-core's number
    from another run; it is quoted as external corroboration for choosing
    `mid`, and this script cannot recompute it.

Usage: verify_coverage_threeway.py <res_dir> [report.txt]
"""
import os
import sys

import numpy as np
import pandas as pd

BC = ['b%02d' % i for i in range(100)]


def main(res, report=None):
    out = []

    def say(s=''):
        print(s, flush=True)
        out.append(s)

    n_ok = n_bad = 0

    def check(name, got, want, tol=0.0005, note=''):
        nonlocal n_ok, n_bad
        if got is None or want is None:
            # a value this script could not obtain is a FAILED check, not a
            # skipped one -- silently passing on None is how an unchecked
            # number gets mistaken for a checked one
            n_bad += 1
            say('  [FAIL] %-52s got %-12s claim %-12s %s'
                % (name, repr(got), repr(want), note or 'value unavailable'))
            return False
        ok = (abs(float(got) - float(want)) <= tol
              if not isinstance(want, str) else got == want)
        n_ok += ok
        n_bad += (not ok)
        say('  [%s] %-52s got %-12s claim %-12s %s'
            % ('OK' if ok else 'FAIL', name,
               ('%.4f' % got) if not isinstance(got, str) else got,
               ('%.4f' % want) if not isinstance(want, str) else want, note))
        return ok

    per = pd.read_csv(os.path.join(res, 'coverage_threeway.tsv'), sep='\t')
    prof = pd.read_csv(os.path.join(res, 'coverage_threeway_profile.tsv'), sep='\t')
    aln = pd.read_csv(os.path.join(res, 'coverage_threeway_alnlen.tsv'), sep='\t')
    rob = pd.read_csv(os.path.join(res, 'coverage_threeway_robust.tsv'), sep='\t')
    cells = pd.read_csv(os.path.join(res, 'coverage_threeway_pubcells.tsv'), sep='\t')
    genes = pd.read_csv(os.path.join(res, 'coverage_threeway_geneset.genes.tsv'),
                        sep='\t')
    est = pd.read_csv(os.path.join(res, 'threeway_published_cellcalls.tsv'),
                      sep='\t', dtype={'unit': str})

    say('=' * 78)
    say('VERIFY  three-way transcript-body coverage')
    say('=' * 78)

    # ---------------------------------------------------------------- scope
    say()
    say('1. SCOPE AND DENOMINATORS (Rule 5: state the denominator)')
    mid = per[per.metric == 'mid']
    for g, n in (('VASA_published', 12), ('VASA_own', 6),
                 ('FLASHseq_native', 4), ('FLASHseq_vasalen', 4)):
        check('n units, %s' % g, len(mid[mid.group == g]), n, 0)
    check('gene set size', len(genes), 4000, 0,
          'one longest transcript per gene, >=1000 nt in BOTH releases')
    check('bins', len(BC), 100, 0)

    # every profile row must be a normalised distribution, or the fractions
    # below are not fractions
    worst = max(abs(prof[BC].sum(axis=1) - 1))
    check('max |sum(profile) - 1| over all rows', worst, 0.0, 1e-6)

    # ---------------------------------------------------------------- calls
    say()
    say('2. PUBLISHED-PLATE CELL SELECTION')
    est['unit'] = est.unit.astype(str).str.zfill(3)
    check('barcodes on the plate', len(est), 384, 0)
    check('mouse by Fig.1d rule (UFI fraction)',
          int((est.call_fig1d == 'mouse').sum()), 173, 0, 'the paper\'s Fig.1d')
    check('mouse by Methods rule (gene fraction)',
          int((est.call_methods == 'mouse').sum()), 141, 0, 'the paper\'s Methods')
    both = int(((est.call_fig1d == 'mouse') & (est.call_methods == 'mouse')).sum())
    check('mouse by BOTH rules', both, 141, 0, 'the set selected from')
    check('rules disagree on', int((est.call_fig1d != est.call_methods).sum()),
          32, 0)
    sel = cells[cells.selected]
    check('cells profiled', len(sel), 12, 0, 'the 12 deepest by mouse UFI')
    assert sel.mouse_both_rules.all(), 'a selected cell is not mouse under both rules'
    say('       every selected cell is mouse under BOTH rules: asserted')
    # the selection must be the deepest, not an arbitrary 12
    ok_order = (cells[cells.mouse_both_rules]
                .sort_values('ufi_mouse', ascending=False)
                .head(12).cell.astype(str).str.zfill(3).tolist())
    got_order = sel.cell.astype(str).str.zfill(3).tolist()
    check('selected == 12 deepest mouse-pure', str(sorted(got_order) == sorted(ok_order)),
          'True', note='by mouse UFI')

    # ---------------------------------------------------------------- the answer
    say()
    say('3. THE HEADLINE CLAIM: the 3\' rise is the LIBRARY, not the PROTOCOL')
    say('   rise = mean(bins 90-99) / mean(bins 40-59) of the MIDPOINT profile')
    r = {g: mid[mid.group == g].rise.values.astype(float)
         for g in mid.group.unique()}
    # the writeup's headline table, cell for cell (rise, 3'/5', n, p50)
    WRITEUP_HEAD = {
        # group:             n   rise   ratio  p50
        'VASA_published':   (12, 0.649, 0.790,  74),
        'VASA_own':         (6,  1.392, 1.188, 127),
        'FLASHseq_native':  (4,  0.994, 0.882, 150),
        'FLASHseq_vasalen': (4,  0.911, 0.734, 107),
    }
    for g, (cn, crise, cratio, cp50) in WRITEUP_HEAD.items():
        s = mid[mid.group == g]
        check('%s n' % g, len(s), cn, 0)
        check('%s rise' % g, s.rise.mean(), crise, 0.0005)
        check("%s 3'/5'" % g, s.ratio.mean(), cratio, 0.0005)
        check('%s p50 (headline table)' % g,
              float(aln[aln.group == g].p50.iloc[0]), cp50, 0)
    check('published rise max', r['VASA_published'].max(), 0.7626, 0.0005)
    check('own rise min', r['VASA_own'].min(), 1.3442, 0.0005)
    disjoint = r['VASA_published'].max() < r['VASA_own'].min()
    check('per-unit ranges DISJOINT', str(disjoint), 'True',
          note='no published cell reaches any own cell')
    # the direction is what the claim rests on: published is BELOW flat, own ABOVE
    check('published rise < 1 (3\' depleted)',
          str(bool(r['VASA_published'].max() < 1.0)), 'True')
    check('own rise > 1 (3\' enriched)', str(bool(r['VASA_own'].min() > 1.0)), 'True')

    # exact permutation, recomputed here rather than read from the log
    from itertools import combinations
    from math import comb
    a, b = r['VASA_published'], r['VASA_own']
    pool = np.concatenate([a, b])
    na, n = len(a), len(pool)
    obs = abs(a.mean() - b.mean())
    cnt = 0
    for cmb in combinations(range(n), na):
        m = np.zeros(n, bool)
        m[list(cmb)] = True
        cnt += abs(pool[m].mean() - pool[~m].mean()) >= obs - 1e-12
    tot = comb(n, na)
    check('exact permutation splits total', tot, 18564, 0)
    check('splits as extreme as observed', cnt, 1, 0)
    check('exact two-sided p', cnt / tot, 5.387e-05, 1e-7)
    psd = np.sqrt(((na - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                  / (na + len(b) - 2))
    check('gap (own - published)', b.mean() - a.mean(), 0.7424, 0.0005)
    check("Cohen's d", (b.mean() - a.mean()) / psd, 10.5, 0.05,
          'within-group sd: pub %.4f, own %.4f' % (a.std(ddof=1), b.std(ddof=1)))
    say('       n=12 vs n=6, so the smallest attainable p is 1/18,564 = %.2e'
        % (1 / tot))
    say('       and the observed split IS that minimum -- the test cannot')
    say('       resolve smaller, so this is a floor, not a small p-value.')

    # ---------------------------------------------------------------- controls
    say()
    say('4. THE READ-LENGTH CONTROL (the first objection to a shape claim)')
    sh = per[per.metric == 'mid_short']
    rs = {g: sh[sh.group == g].rise.values.astype(float) for g in sh.group.unique()}
    # EVERY cell of the writeup's aligned-length table, not just p50. A review
    # caught a transcribed p75 of 130 for FLASHseq_vasalen where the run
    # measured 129 (130 is that group's p95/p99) -- and the earlier version of
    # this section asserted only p50 and frac_le_80nt, so the check could not
    # catch it. The claims below are the table as published, cell for cell.
    WRITEUP_ALNLEN = {
        # group:            p05  p25  p50  p75   mean   frac_le_80nt(%)
        'VASA_published':   (46,  71,  74,  75,  69.8, 100.0),
        'VASA_own':         (26, 100, 127, 129, 109.7,  15.3),
        'FLASHseq_native':  (51, 107, 150, 151, 127.1,  15.3),
        'FLASHseq_vasalen': (28,  71, 107, 129,  96.7,  30.8),
    }
    for g, (c5, c25, c50, c75, cmean, cle80) in WRITEUP_ALNLEN.items():
        row = aln[aln.group == g]
        assert len(row) == 1, (g, len(row))
        row = row.iloc[0]
        for nm, got, want, tol in (('p05', row.p05, c5, 0), ('p25', row.p25, c25, 0),
                                   ('p50', row.p50, c50, 0), ('p75', row.p75, c75, 0),
                                   ('mean', row['mean'], cmean, 0.05),
                                   ('<=80nt %', 100 * row.frac_le_80nt, cle80, 0.05)):
            check('%s %s' % (g, nm), float(got), want, tol)
    check('published frac <= 80 nt',
          float(aln[aln.group == 'VASA_published'].frac_le_80nt.iloc[0]), 1.0, 1e-6,
          'the whole published library is inside the shared regime')
    check('published rise, reads <=80nt', rs['VASA_published'].mean(), 0.649, 0.0005)
    check('own rise, reads <=80nt', rs['VASA_own'].mean(), 1.431, 0.0005)
    disj_s = rs['VASA_published'].max() < rs['VASA_own'].min()
    check('length-matched ranges DISJOINT', str(disj_s), 'True',
          'the gap SURVIVES matching on read length')
    say('       own rise RISES under length matching (%.3f -> %.3f), so the gap'
        % (r['VASA_own'].mean(), rs['VASA_own'].mean()))
    say('       is not produced by the own plate\'s longer reads.')

    say()
    say('5. ROBUSTNESS TO THE PUBLISHED CELLS\' ~10x LOWER DEPTH')
    rm = rob[rob.metric == 'mid']
    for g, want in (('VASA_published', 0.616), ('VASA_own', 1.398)):
        s = rm[rm.group == g]
        check('%s rise, tx with >=20 reads' % g, s.rise.mean(), want, 0.0005,
              'n_tx=%.0f' % s.n_tx.mean())
    pub_rb = rm[rm.group == 'VASA_published'].rise.values
    own_rb = rm[rm.group == 'VASA_own'].rise.values
    check('robust ranges still DISJOINT',
          str(bool(pub_rb.max() < own_rb.min())), 'True')

    say()
    say('6. THE BASE METRIC IS DIAGNOSTIC ONLY -- and still disagrees by ~2x')
    ba = per[per.metric == 'base']
    for g, want in (('VASA_published', 0.537), ('VASA_own', 0.560),
                    ('FLASHseq_native', 0.455), ('FLASHseq_vasalen', 0.460)):
        check("%s base 3'/5'" % g, ba[ba.group == g].ratio.mean(), want, 0.0005)
    say('       Under `base`, published (0.537) and own (0.560) are nearly EQUAL')
    say('       -- the metric that disagrees with qualimap also erases the')
    say('       library difference that `mid` and the paper-scale profile show.')
    say('       Cause of the base-vs-mid disagreement: STILL UNKNOWN. Edge')
    say('       clipping remains the leading UNTESTED hypothesis; the earlier')
    say('       refutation is withdrawn and is not re-derived here.')

    say()
    say('7. LOSS DECOMPOSITION IS ARITHMETICALLY EXACT (not a mechanism claim)')
    m = per[per.metric == 'mid']
    tot_b = (m.bases_binned + m.lost_txend + m.lost_txstart + m.lost_internal
             + m.lost_modelerr + m.lost_unclassified)
    check('max |sum(classes) - bases_aligned|',
          float(np.max(np.abs(tot_b - m.bases_aligned))), 0.0, 0.0)
    check('lost_modelerr total (tripwire, must be 0)',
          float(m.lost_modelerr.sum()), 0.0, 0.0)
    say('       reads overhanging the annotated 3\' end carry only %.3f%% of'
        % (100 * m.lost_txend.sum() / m.bases_aligned.sum()))
    say('       aligned bases pooled over all units -- too small to be the')
    say('       whole base-vs-mid gap, but this bounds the loss, it does not')
    say('       test the hypothesis.')

    say()
    say('7b. THE REPORTED METRIC\'S OWN LIMITATION, MEASURED NOT ASSUMED')
    say('    `mid` votes at the GENOMIC-SPAN midpoint, so a junction-spanning')
    say('    read whose midpoint lands in an intron casts NO vote while still')
    say('    counting in reads_placed. The dropout is large and is coupled to')
    say('    read length -- the variable `mid` exists to neutralise.')
    dp = os.path.join(res, 'coverage_threeway_middrop.tsv')
    if os.path.exists(dp):
        dd = pd.read_csv(dp, sep='\t')
        # the writeup's dropout table: BOTH columns, not only the percentage
        for g, w_placed, want in (('VASA_published', 2463966, 30.554),
                                  ('VASA_own', 14402031, 41.896),
                                  ('FLASHseq_native', 12000000, 47.798),
                                  ('FLASHseq_vasalen', 12000000, 40.591)):
            s = dd[dd.group == g]
            check('%s reads placed' % g, float(s.reads_placed.sum()), w_placed, 0)
            check('%s midpoint-vote dropout (%%)' % g,
                  100 * s.dropped.sum() / s.reads_placed.sum(), want, 0.001)
        rr = float(np.corrcoef(dd.frac_dropped, dd.alnlen_p50)[0, 1])
        check('corr(dropout, p50 aligned length)', rr, 0.982, 0.001,
              'so mid is NOT read-length-neutral')
        # the dropout cannot invert the claim: base has no such dropout
        pbase = prof[(prof.group == 'VASA_published') & (prof.metric == 'base')]
        obase = prof[(prof.group == 'VASA_own') & (prof.metric == 'base')]

        def rise_of(row):
            v = row[BC].values.ravel().astype(float)
            return float(v[90:].mean() / v[40:60].mean())

        rb = rise_of(obase) / rise_of(pbase)
        pmid = prof[(prof.group == 'VASA_published') & (prof.metric == 'mid')]
        omid = prof[(prof.group == 'VASA_own') & (prof.metric == 'mid')]
        rm2 = rise_of(omid) / rise_of(pmid)
        check('own/published rise ratio, by base', rb, 1.850, 0.001)
        check('own/published rise ratio, by mid', rm2, 2.144, 0.001)
        check('own is 3\'-heavier under BOTH metrics',
              str(bool(rb > 1 and rm2 > 1)), 'True',
              'the dropout cannot invert the claim\'s direction')
    else:
        say('  FAIL %s absent -- the dropout is unquantified' % dp)
        n_bad += 1

    say()
    say('8. THE ANNOTATION-RELEASE CONFOUND IS PRESENT AND UNREMOVED')
    lr = genes.len_ratio_116_99.values
    check('median longest-transcript length ratio E116/E99',
          float(np.median(lr)), 1.1129, 0.0005)
    check('genes with identical length (%)', 100 * float((lr == 1).mean()), 13.8,
          0.05)
    check('genes within +-5%% (%)', 100 * float((abs(lr - 1) <= 0.05).mean()),
          37.6, 0.05)
    check('strand disagreements', int((~genes.same_strand).sum()), 4, 0)
    ident = ((genes.len_ratio_116_99 == 1) & genes.same_strand
             & (genes.nexon_E99 == genes.nexon_E116))
    check('genes with a fully identical model', int(ident.sum()), 552, 0,
          '%.1f%% of the gene set' % (100 * ident.mean()))

    say()
    say('8b. THE ONE MECHANISM BY WHICH THE RELEASE COULD FAKE THIS, TESTED')
    say('    E99 models are the SHORTER ones (median 11.3%). If that truncates')
    say('    the published arm\'s annotated 3\' end, reads at the true transcript')
    say('    end would fall BEYOND the model and be lost instead of binned,')
    say('    depleting published\'s last bins artefactually. That mechanism')
    say('    predicts published lost_txend > own lost_txend.')
    pubm, ownm = m[m.group == 'VASA_published'], m[m.group == 'VASA_own']
    p_pct = 100 * pubm.lost_txend.sum() / pubm.bases_aligned.sum()
    o_pct = 100 * ownm.lost_txend.sum() / ownm.bases_aligned.sum()
    check('published lost_txend (% of aligned bases)', p_pct, 0.212, 0.0005)
    check('own lost_txend (% of aligned bases)', o_pct, 0.297, 0.0005)
    check('prediction FALSIFIED (published is LOWER, not higher)',
          str(bool(p_pct < o_pct)), 'True',
          'so 3\'-end model truncation is not producing the gap')

    say()
    say('8c. AND THE GAP SURVIVES ON READS THAT LOSE NO BASES AT ALL')
    say('    mid_fullexon restricts to reads every base of which landed in the')
    say('    model, so no model-edge effect of any kind can touch them.')
    mfe = per[per.metric == 'mid_fullexon']
    pfe = mfe[mfe.group == 'VASA_published'].rise.values.astype(float)
    ofe = mfe[mfe.group == 'VASA_own'].rise.values.astype(float)
    check('published rise, fullexon reads only', pfe.mean(), 0.6517, 0.0005)
    check('own rise, fullexon reads only', ofe.mean(), 1.3848, 0.0005)
    check('fullexon-only ranges STILL DISJOINT',
          str(bool(pfe.max() < ofe.min())), 'True',
          'pub max %.4f < own min %.4f' % (pfe.max(), ofe.min()))

    say()
    say('       The published arm is quantified on E99/GRCm38 models and the')
    say('       other two on E116/GRCm39; the models differ (median 11.3%')
    say('       longer in E116). This is NOT controlled here. It is why the')
    say('       claim is stated on the DIRECTION of the rise -- published')
    say('       BELOW flat, own ABOVE -- which no rescaling of a transcript')
    say('       axis can invert, rather than on the size of the difference.')

    say()
    say('9. REGRESSION: the already-reported two-way numbers survive')
    up = os.path.join(os.path.dirname(res.rstrip('/')), 'flashseq_vasa',
                      'gene_coverage_profile.tsv')
    if os.path.exists(up):
        u = pd.read_csv(up, sep='\t')
        ub = [c for c in u.columns if c.startswith('b')]
        mp = {'VASA': 'VASA_own', 'FLASH-seq native': 'FLASHseq_native',
              'FLASH-seq vasalen': 'FLASHseq_vasalen'}
        for _, row in u[u.metric == 'mid'].iterrows():
            g = mp.get(row.grp)
            if g is None:
                continue
            v = row[ub].values.astype(float)
            v = v / v.sum()
            old = v[80:].sum() / v[:20].sum()
            new = mid[mid.group == g].ratio.mean()
            check('%s: same side of 1.0 as upstream' % g,
                  str(bool((old - 1) * (new - 1) > 0)), 'True',
                  'upstream %.3f, here %.3f' % (old, new))
    else:
        say('  upstream table absent: %s' % up)

    say()
    say('10. THE WRITEUP\'S OWN TABLES, PARSED FROM THE MARKDOWN')
    say('    Sections 3, 4 and 7b assert hardcoded claim values that a human')
    say('    typed into BOTH this script and THREEWAY_COVERAGE.md. That catches')
    say('    a drifting pipeline but not a transcription slip made in both')
    say('    places -- and a review did catch one (p75 130 for FLASHseq_vasalen')
    say('    where the run measured 129; 130 is that group\'s p95). So the')
    say('    markdown is parsed here and every numeric cell of its two')
    say('    length/rise tables is re-derived from the TSVs directly.')
    # res is <project>/res/threeway, so the repo root is two levels up. Try the
    # committed location first, then res/ itself, then this script's own dir --
    # and if none resolves, FAIL rather than skip: an unfindable writeup means
    # its tables are unchecked, which is exactly the hole this section closes.
    root = os.path.dirname(os.path.dirname(os.path.abspath(res.rstrip('/'))))
    cands = [os.path.join(root, 'code', 'flashseq_vasa', 'threeway',
                          'THREEWAY_COVERAGE.md'),
             os.path.join(res, 'THREEWAY_COVERAGE.md'),
             os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'THREEWAY_COVERAGE.md')]
    md = next((p for p in cands if os.path.exists(p)), cands[0])
    if os.path.exists(md):
        say('    reading %s' % md)
        LAB = {'VASA-seq, published plate': 'VASA_published',
               'VASA-seq, own plate': 'VASA_own',
               'FLASH-seq, native': 'FLASHseq_native',
               'FLASH-seq, VASA-trimmed': 'FLASHseq_vasalen'}
        rows = {}
        for line in open(md):
            if not line.lstrip().startswith('|'):
                continue
            cell = [c.strip() for c in line.strip().strip('|').split('|')]
            if cell and cell[0] in LAB:
                rows.setdefault(LAB[cell[0]], []).append(cell[1:])
        check('groups found in the markdown tables', len(rows), 4, 0)

        import re as _re

        def num(x):
            """First number in a markdown cell, tolerating units and emphasis.

            Cells look like '12 cells', '74 nt', '**0.649**', '30.554%',
            '2,463,966'. Returns None only when the cell holds no number at
            all -- and a None reaching a check is treated as a FAILURE below,
            never skipped, because a cell this parser cannot read is a cell
            nobody is checking.
            """
            m = _re.search(r'-?\d[\d,]*\.?\d*',
                           x.replace('**', '').replace(',', ''))
            return float(m.group()) if m else None

        for g, blocks in rows.items():
            srow = aln[aln.group == g].iloc[0]
            mrow = mid[mid.group == g]
            for cells in blocks:
                v = []
                for c in cells:
                    try:
                        v.append(num(c))
                    except ValueError:
                        v.append(None)
                if len(cells) == 4 and 'cell' in cells[0] or 'librar' in cells[0]:
                    # headline table: n | rise | 3'/5' | p50
                    check('md headline %s rise' % g, mrow.rise.mean(), v[1], 0.0005)
                    check("md headline %s 3'/5'" % g, mrow.ratio.mean(), v[2],
                          0.0005)
                    check('md headline %s p50' % g, float(srow.p50), v[3], 0)
                elif len(cells) == 6:
                    # alnlen table: p05 p25 p50 p75 mean <=80nt
                    for nm, got, want, tol in (
                            ('p05', srow.p05, v[0], 0), ('p25', srow.p25, v[1], 0),
                            ('p50', srow.p50, v[2], 0), ('p75', srow.p75, v[3], 0),
                            ('mean', srow['mean'], v[4], 0.05),
                            ('<=80nt', 100 * srow.frac_le_80nt, v[5], 0.05)):
                        check('md alnlen %s %s' % (g, nm), float(got), want, tol)
                elif len(cells) == 2 and dp and os.path.exists(dp):
                    # dropout table: reads placed | dropout %
                    s = dd[dd.group == g]
                    check('md dropout %s reads placed' % g,
                          float(s.reads_placed.sum()), v[0], 0)
                    check('md dropout %s pct' % g,
                          100 * s.dropped.sum() / s.reads_placed.sum(), v[1],
                          0.001)
    else:
        say('  FAIL THREEWAY_COVERAGE.md not found -- its tables are unchecked')
        n_bad += 1

    say()
    say('=' * 78)
    say('%d checks passed, %d FAILED' % (n_ok, n_bad))
    say('=' * 78)
    if report:
        with open(report, 'w') as fh:
            fh.write('\n'.join(out) + '\n')
        print('wrote %s' % report)
    return 1 if n_bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None))
