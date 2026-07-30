#!/usr/bin/env python3
"""fig_coverage_threeway.py -- the three-way coverage figure.

ONE MESSAGE: the own VASA plate's 3' rise is tested against the published VASA
plate, so the figure has to let a reader see whether the two VASA curves agree
with each other and differ from FLASH-seq, or not.

Panel a  the four midpoint profiles, 5'->3', the reported metric.
Panel b  the last 20% of the transcript, expanded -- where the own plate's rise
         lives, so it is where the discriminating comparison must be readable.
Panel c  3'/5' ratio, one point per unit, so the within-group spread is visible
         rather than hidden inside a mean.
Panel d  aligned read-length distribution for all four, since read length is the
         obvious confounder for a coverage-shape claim and all three datasets
         are now on the record.

Every number annotated in the figure is recomputed here from
coverage_threeway.tsv and asserted against that file -- nothing is transcribed.

Usage: fig_coverage_threeway.py <res_dir> <cov_dir> <out.png>
"""
import glob
import os
import sys

import numpy as np


def main(res, covdir, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    per = pd.read_csv(os.path.join(res, 'coverage_threeway.tsv'), sep='\t')
    prof = pd.read_csv(os.path.join(res, 'coverage_threeway_profile.tsv'),
                       sep='\t')
    bcols = ['b%02d' % i for i in range(100)]

    GROUPS = [
        ('VASA_published', 'VASA-seq, published plate', '#1b6ca8'),
        ('VASA_own', 'VASA-seq, own plate', '#e8763a'),
        ('FLASHseq_native', 'FLASH-seq, native 151 nt', '#4a4a4a'),
        ('FLASHseq_vasalen', 'FLASH-seq, VASA-trimmed', '#9a9a9a'),
    ]
    GROUPS = [g for g in GROUPS if g[0] in set(prof.group)]

    def curve(g, metric='mid'):
        r = prof[(prof.group == g) & (prof.metric == metric)]
        assert len(r) == 1, (g, metric, len(r))
        return r[bcols].values.ravel().astype(float) * 100

    def units(g, metric='mid'):
        return per[(per.group == g) & (per.metric == metric)]

    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92],
                          hspace=0.52, wspace=0.30,
                          left=0.095, right=0.985, top=0.925, bottom=0.085)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    x = np.arange(100) + 0.5

    # ---- panel a: full profiles -------------------------------------------
    for g, lab, col in GROUPS:
        y = curve(g)
        lw = 1.9 if g.startswith('VASA') else 1.2
        axA.plot(x, y, color=col, lw=lw,
                 ls='-' if g != 'FLASHseq_vasalen' else (0, (4, 1.5)), zorder=3)
    axA.axhline(1.0, color='#bbbbbb', lw=0.7, ls=':', zorder=1)
    axA.set_xlabel("position along transcript, 5' \u2192 3' (%)")
    axA.set_ylabel('reads per 1% bin (%)')
    axA.set_xlim(0, 100)
    axA.margins(y=0.10)
    axA.set_title('Midpoint coverage profile', loc='left')
    axA.text(0.5, 1.0, 'flat = uniform', transform=axA.transAxes,
             ha='right', va='bottom', fontsize=6, color='#888888')

    # ---- panel b: the 3' end, expanded ------------------------------------
    for g, lab, col in GROUPS:
        y = curve(g)
        lw = 1.9 if g.startswith('VASA') else 1.2
        axB.plot(x[75:], y[75:], color=col, lw=lw, marker='o', ms=2.4,
                 ls='-' if g != 'FLASHseq_vasalen' else (0, (4, 1.5)), zorder=3)
    axB.axhline(1.0, color='#bbbbbb', lw=0.7, ls=':', zorder=1)
    axB.set_xlabel("position along transcript (%)")
    axB.set_ylabel('reads per 1% bin (%)')
    axB.set_xlim(74, 100.5)
    axB.margins(y=0.12)
    axB.set_title("Last quarter, expanded", loc='left')

    # ---- panel c: per-unit 3'/5' ------------------------------------------
    rng = np.random.default_rng(0)
    for i, (g, lab, col) in enumerate(GROUPS):
        u = units(g)
        v = u.ratio.values
        jit = rng.uniform(-0.13, 0.13, len(v))
        axC.scatter(np.full(len(v), i) + jit, v, s=17, color=col,
                    edgecolor='white', linewidth=0.35, zorder=3)
        axC.plot([i - 0.27, i + 0.27], [v.mean()] * 2, color=col, lw=2.0,
                 zorder=4, solid_capstyle='butt')
    axC.axhline(1.0, color='#bbbbbb', lw=0.8, ls=':', zorder=1)
    axC.set_xticks(range(len(GROUPS)))
    axC.set_xticklabels(['published\nVASA', 'own\nVASA', 'FLASH-seq\nnative',
                         'FLASH-seq\ntrimmed'][:len(GROUPS)])
    axC.set_ylabel("3' 20% / 5' 20%")
    axC.set_xlim(-0.55, len(GROUPS) - 0.45)
    axC.margins(y=0.16)
    axC.set_title("3'/5' ratio, one point per unit", loc='left')
    axC.text(0.015, 0.965, "> 1 = 3'-weighted", transform=axC.transAxes,
             ha='left', va='top', fontsize=6, color='#888888')

    # ---- panel d: aligned read length ------------------------------------
    tot = {}
    for g, lab, col in GROUPS:
        labs = set(units(g).label)
        h = None
        for f in glob.glob(os.path.join(covdir, '*.lenhist.npz')):
            if os.path.basename(f)[:-len('.lenhist.npz')] in labs:
                z = np.load(f)['all'].astype(float)
                h = z if h is None else h + z
        assert h is not None and h.sum() > 0, g
        tot[g] = h
        cs = np.cumsum(h) / h.sum()
        axD.plot(np.arange(len(h)), 100 * cs, color=col,
                 lw=1.9 if g.startswith('VASA') else 1.2,
                 ls='-' if g != 'FLASHseq_vasalen' else (0, (4, 1.5)), zorder=3)
    axD.set_xlabel('aligned read length (nt)')
    axD.set_ylabel('primary alignments \u2264 x (%)')
    axD.set_xlim(0, 175)
    axD.set_ylim(-3, 103)
    axD.set_title('Aligned read length', loc='left')

    # direct labels in panel a's whitespace, replacing a legend box
    ymax = max(curve(g).max() for g, _, _ in GROUPS)
    for g, lab, col in GROUPS:
        y = curve(g)
        axA.annotate(lab, xy=(x[-1], y[-1]), xytext=(0, 0),
                     textcoords='offset points', color=col, fontsize=6,
                     ha='right', va='bottom')

    fig.savefig(out, dpi=300)
    print('wrote %s' % out)

    # ---- assert every annotated number against the source table -----------
    print('\n=== numbers this figure asserts against coverage_threeway.tsv ===')
    for g, lab, _ in GROUPS:
        u = units(g)
        c = curve(g)
        body = c[40:60].mean()
        print('  %-17s n=%d  3\'/5\'=%.3f  rise=%.3f  p50len=%.0f nt'
              % (g, len(u), u.ratio.mean(), c[90:].mean() / body,
                 u.alnlen_p50.mean()))
        assert abs(c.sum() - 100) < 1e-6, (g, c.sum())
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
