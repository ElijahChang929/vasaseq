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
    check("published rise (mean of 12 cells)", r['VASA_published'].mean(), 0.649,
          0.0005)
    check("own rise (mean of 6 cells)", r['VASA_own'].mean(), 1.392, 0.0005)
    check("FLASH-seq native rise", r['FLASHseq_native'].mean(), 0.994, 0.0005)
    check("FLASHseq_vasalen rise", r['FLASHseq_vasalen'].mean(), 0.911, 0.0005)
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
    check('published p50 aligned length',
          float(aln[aln.group == 'VASA_published'].p50.iloc[0]), 74, 0)
    check('own p50 aligned length',
          float(aln[aln.group == 'VASA_own'].p50.iloc[0]), 127, 0)
    check('FLASH-seq native p50', float(aln[aln.group == 'FLASHseq_native'].p50.iloc[0]),
          150, 0)
    check('FLASH-seq vasalen p50',
          float(aln[aln.group == 'FLASHseq_vasalen'].p50.iloc[0]), 107, 0)
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
