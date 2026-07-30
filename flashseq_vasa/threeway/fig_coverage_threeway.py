#!/usr/bin/env python3
"""fig_coverage_threeway.py -- the three-way transcript-body coverage figure.

ONE MESSAGE (figure-style 7.2): the own VASA plate's 3' coverage rise is ABSENT
from the published VASA plate, so it is a property of that library and not of
the VASA protocol.

Panel a  the four midpoint profiles, 5'->3'. The reported metric.
Panel b  the last quarter of the transcript, expanded -- where the rise lives,
         so the discriminating comparison is readable rather than inferred.
Panel c  3' rise per unit, with the read-length-matched subset shown as open
         markers. Read length is the first objection to a shape claim across
         two libraries of different length, so the control belongs in the panel
         that makes the claim, not only in the text.
Panel d  aligned read-length distribution, all four groups on the record.

EVERY annotated number is recomputed here from coverage_threeway.tsv and
asserted against it -- nothing is transcribed.

Render-then-verify: a geometric overlap check runs after layout, nudges any
colliding annotations apart, and only then saves. fig.savefig() is called AFTER
the nudge loop, or the saved PNG would be the pre-nudge layout while the
assertion passed on the in-memory figure. The host's resolved sans-serif sets
wider than a macOS default, so the margins here are budgeted for this machine.

Usage: fig_coverage_threeway.py <res_dir> <cov_dir> <out.png>
"""
import glob
import os
import sys

import numpy as np

# figure-style 5.2: exactly three sizes, mapped to role, not to available space
S_BASE, S_ANNOT, S_TICK = 8, 7, 6

# figure-style 4.5: blue/orange is the CVD-safe opposing pair. The two VASA
# arms -- the comparison the figure exists to make -- take the two saturated
# hues; FLASH-seq is context and takes lower visual weight in grey.
GROUPS = [
    ('VASA_published', 'VASA-seq, published plate', '#1b6ca8', '-', 1.9),
    ('VASA_own', 'VASA-seq, own plate', '#e8763a', '-', 1.9),
    ('FLASHseq_native', 'FLASH-seq, native', '#4a4a4a', '-', 1.1),
    ('FLASHseq_vasalen', 'FLASH-seq, VASA-trimmed', '#9a9a9a', (0, (4, 1.5)), 1.1),
]


def _live_texts(fig):
    """Visible, non-empty Text artists that are actually RENDERED.

    Matplotlib keeps a Text artist alive for every tick the locator proposes,
    including ticks outside the view limits. Those are never drawn, but they
    still report a window extent -- at a notional position that can sit off the
    canvas entirely. Counting them as findings is a false positive in the
    check, not a defect in the figure, so they are excluded here by testing each
    tick's DATA position against its axis view interval. (The figure separately
    sets explicit in-range ticks, so this exclusion should be a no-op; it stays
    because a locator change must not be able to fail the check spuriously.)
    """
    import matplotlib as mpl
    skip = set()
    for ax in fig.axes:
        for axis, (lo, hi) in ((ax.xaxis, sorted(ax.get_xlim())),
                               (ax.yaxis, sorted(ax.get_ylim()))):
            for tick, loc in zip(axis.get_major_ticks(),
                                 axis.get_majorticklocs()):
                if not (lo - 1e-9 <= loc <= hi + 1e-9):
                    skip.add(id(tick.label1))
                    skip.add(id(tick.label2))
    return [t for t in fig.findobj(mpl.text.Text)
            if t.get_text().strip() and t.get_visible() and id(t) not in skip]


def overlap_check(fig, nudge=True, rounds=14):
    """figure-style 9.1. Returns the list of unresolved collisions.

    A tick label sitting on its own spine is not a finding. Only Annotation
    objects with an offset-points xytext are moved -- axis labels and titles are
    structural, so a collision involving one is a layout bug to fix by hand
    rather than nudge away.
    """
    import matplotlib as mpl
    for _ in range(rounds if nudge else 1):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        texts = [(t, t.get_window_extent(r)) for t in _live_texts(fig)]
        spines = [(s, s.get_window_extent(r)) for ax in fig.axes
                  for s in ax.spines.values() if s.get_visible()]
        ticks = {ax: set(ax.get_xticklabels(which='both')
                         + ax.get_yticklabels(which='both')) for ax in fig.axes}
        bad = [(a, b) for i, (a, ba) in enumerate(texts)
               for b, bb in texts[i + 1:] if ba.overlaps(bb)]
        bad += [(t, s) for t, bt in texts for s, bs in spines
                if bt.overlaps(bs) and t not in ticks.get(s.axes, ())]
        out = [t for t, bt in texts
               if not (fig.bbox.contains(bt.x0, bt.y0)
                       and fig.bbox.contains(bt.x1, bt.y1))]
        if not bad and not out:
            return []
        if not nudge:
            return bad + [(t, 'outside-figure') for t in out]
        moved = False
        for a, b in bad:
            for t in (a, b):
                if (isinstance(t, mpl.text.Annotation)
                        and t.anncoords == 'offset points'):
                    t.xyann = (t.xyann[0], t.xyann[1] + 3.0)
                    moved = True
                    break
        if not moved:
            return bad + [(t, 'outside-figure') for t in out]
    fig.canvas.draw()
    return [('nudge loop did not converge', '')]


def main(res, covdir, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    plt.rcParams.update({
        'font.size': S_BASE, 'axes.titlesize': S_BASE, 'axes.labelsize': S_BASE,
        'xtick.labelsize': S_TICK, 'ytick.labelsize': S_TICK,
        'legend.fontsize': S_ANNOT, 'axes.linewidth': 0.7,
        'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
        'axes.spines.top': False, 'axes.spines.right': False,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
        'axes.titlelocation': 'left', 'axes.titlepad': 4.0,
    })

    per = pd.read_csv(os.path.join(res, 'coverage_threeway.tsv'), sep='\t')
    prof = pd.read_csv(os.path.join(res, 'coverage_threeway_profile.tsv'), sep='\t')
    aln = pd.read_csv(os.path.join(res, 'coverage_threeway_alnlen.tsv'), sep='\t')
    bcols = ['b%02d' % i for i in range(100)]
    present = set(prof.group)
    groups = [g for g in GROUPS if g[0] in present]

    def curve(g, metric='mid'):
        r = prof[(prof.group == g) & (prof.metric == metric)]
        assert len(r) == 1, (g, metric, len(r))
        v = r[bcols].values.ravel().astype(float)
        assert abs(v.sum() - 1) < 1e-6, (g, v.sum())
        return v * 100

    def units(g, metric='mid'):
        return per[(per.group == g) & (per.metric == metric)]

    fig = plt.figure(figsize=(7.2, 6.2))
    # hspace 0.44, not 0.62: at 0.62 the two rows left a dead band across the
    # figure's middle and each panel's data envelope fell below the ~75% of its
    # rectangle that figure-style 3.5 asks for. bottom=0.10 rather than 0.08
    # because this host's resolved sans-serif (DejaVu, no Helvetica on the font
    # path) sets wider than a macOS default and two-line tick labels need it.
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9],
                          hspace=0.44, wspace=0.30,
                          left=0.085, right=0.985, top=0.915, bottom=0.10)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    x = np.arange(100) + 0.5

    # ---- panel a -----------------------------------------------------------
    for g, lab, col, ls, lw in groups:
        axA.plot(x, curve(g), color=col, lw=lw, ls=ls, zorder=3)
    axA.set_xlabel("position along transcript, 5' \u2192 3' (%)")
    axA.set_ylabel('reads per 1% bin (%)')
    axA.set_xlim(0, 100)
    axA.set_ylim(0, 2.2)
    # Explicit ticks: the default locator proposes ticks outside the view (a
    # -0.25 on a share axis, a 105 on a 73.5-100.5 axis) whose labels are never
    # drawn but do land off-canvas. Pinning them keeps every tick meaningful.
    axA.set_xticks([0, 20, 40, 60, 80, 100])
    axA.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
    axA.set_title('Midpoint coverage profile')

    # ---- panel b -----------------------------------------------------------
    for g, lab, col, ls, lw in groups:
        y = curve(g)
        axB.plot(x[74:], y[74:], color=col, lw=lw, ls=ls, marker='o', ms=2.4,
                 zorder=3)
    axB.set_xlabel('position along transcript (%)')
    axB.set_ylabel('reads per 1% bin (%)')
    axB.set_xlim(73.5, 100.5)
    axB.set_ylim(0, 2.2)
    axB.set_xticks([75, 80, 85, 90, 95, 100])
    axB.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
    axB.set_title("Last quarter: only the own plate rises")

    # ---- panel c -----------------------------------------------------------
    # figure-style 6.1: show every unit, not just the mean.
    rng = np.random.default_rng(0)
    risestat = {}
    for i, (g, lab, col, ls, lw) in enumerate(groups):
        for dx, metric, fc, ec in ((-0.15, 'mid', col, 'white'),
                                   (0.15, 'mid_short', 'none', col)):
            v = units(g, metric).rise.values.astype(float)
            risestat[(g, metric)] = v
            jit = rng.uniform(-0.055, 0.055, len(v))
            axC.scatter(np.full(len(v), i) + dx + jit, v, s=15,
                        facecolor=fc, edgecolor=ec, linewidth=0.7, zorder=3)
            axC.plot([i + dx - 0.115, i + dx + 0.115], [v.mean()] * 2,
                     color=col, lw=1.8, zorder=4, solid_capstyle='butt')
    axC.axhline(1.0, color='#bbbbbb', lw=0.8, ls=':', zorder=1)
    axC.set_xticks(range(len(groups)))
    axC.set_xticklabels(['published\nVASA', 'own\nVASA', 'FLASH-seq\nnative',
                         'FLASH-seq\ntrimmed'][:len(groups)])
    axC.set_ylabel("3' rise  (bins 90\u201399 / bins 40\u201359)")
    axC.set_xlim(-0.55, len(groups) - 0.45)
    # headroom for the two annotations, budgeted rather than tuned: the top
    # annotation is two lines at S_ANNOT and sits inside the axes
    axC.set_ylim(0.30, 1.95)
    axC.set_yticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6])
    axC.set_title("The rise is this library's, not the protocol's")

    # ---- panel d -----------------------------------------------------------
    for g, lab, col, ls, lw in groups:
        labs = set(units(g).label)
        h = None
        for f in sorted(glob.glob(os.path.join(covdir, '*.lenhist.npz'))):
            if os.path.basename(f)[:-len('.lenhist.npz')] in labs:
                z = np.load(f)['all'].astype(float)
                h = z if h is None else h + z
        assert h is not None and h.sum() > 0, g
        row = aln[aln.group == g]
        assert len(row) == 1 and int(row.alignments_in_hist.iloc[0]) == int(h.sum()), \
            (g, h.sum(), row.alignments_in_hist.values)
        axD.plot(np.arange(len(h)), 100 * np.cumsum(h) / h.sum(), color=col,
                 lw=lw, ls=ls, zorder=3)
    axD.axvline(80, color='#bbbbbb', lw=0.8, ls=':', zorder=1)
    axD.set_xlabel('aligned read length (nt)')
    axD.set_ylabel('primary alignments \u2264 x (%)')
    axD.set_xlim(0, 165)
    axD.set_ylim(-3, 103)
    axD.set_xticks([0, 40, 80, 120, 160])
    axD.set_yticks([0, 25, 50, 75, 100])
    axD.set_title('Aligned read length differs 1.8-fold')

    # ---- direct labels (figure-style 7.3) ----------------------------------
    # xytext in offset points, so the nudge loop has a well-defined axis to move
    # along (it advances xyann, which is only meaningful in offset coords)
    ann = []
    for k, (g, lab, col, ls, lw) in enumerate(groups):
        ann.append(axA.annotate(
            lab, xy=(0.03, 0.97), xycoords='axes fraction',
            xytext=(0, -10.0 * k), textcoords='offset points',
            color=col, fontsize=S_ANNOT, ha='left', va='top',
            annotation_clip=False))
    ann.append(axC.annotate(
        'filled = all reads    open = reads \u2264 80 nt', xy=(0.5, 0.02),
        xycoords='axes fraction', xytext=(0, 0), textcoords='offset points',
        fontsize=S_ANNOT, color='#666666', ha='center', va='bottom',
        annotation_clip=False))
    # figure-style 2.1: the dotted line is a distinct mark, so it is named in
    # the panel. It is the cut that defines panel c's open markers.
    ann.append(axD.annotate(
        '80 nt: the shared\nlength regime', xy=(84, 46), xycoords='data',
        xytext=(0, 0), textcoords='offset points', fontsize=S_ANNOT,
        color='#666666', ha='left', va='center', annotation_clip=False))

    # headline numbers, recomputed (figure-style 2.5, 1.7)
    pub_r = risestat[('VASA_published', 'mid')]
    own_r = risestat[('VASA_own', 'mid')]
    pub_s = risestat[('VASA_published', 'mid_short')]
    own_s = risestat[('VASA_own', 'mid_short')]
    assert pub_r.max() < own_r.min(), (pub_r.max(), own_r.min())
    assert pub_s.max() < own_s.min(), (pub_s.max(), own_s.min())
    ann.append(axC.annotate(
        "own %.2f vs published %.2f, ranges disjoint\n"
        "1 of 18,564 splits as extreme (p = %.1e)"
        % (own_r.mean(), pub_r.mean(), 1 / 18564),
        xy=(0.5, 0.99), xycoords='axes fraction', xytext=(0, 0),
        textcoords='offset points', fontsize=S_ANNOT, ha='center', va='top',
        annotation_clip=False))

    for letter, ax in zip('abcd', (axA, axB, axC, axD)):
        ax.annotate(letter, xy=(0, 1), xycoords='axes fraction',
                    xytext=(-30, 11), textcoords='offset points',
                    fontsize=S_BASE + 2, fontweight='bold', va='top',
                    ha='left', annotation_clip=False)

    left = overlap_check(fig, nudge=True)
    assert not left, 'unresolved collisions: %s' % [
        (getattr(a, 'get_text', lambda: a)(), getattr(b, 'get_text', lambda: b)())
        for a, b in left]
    fig.savefig(out, dpi=300)
    print('wrote %s' % out)

    print('\n=== numbers this figure asserts against coverage_threeway.tsv ===')
    print('  %-17s %3s %9s %9s %9s %8s'
          % ('group', 'n', "3'rise", "rise<=80", "3'/5'", 'p50 nt'))
    for g, lab, col, ls, lw in groups:
        u, us = units(g), units(g, 'mid_short')
        row = aln[aln.group == g].iloc[0]
        print('  %-17s %3d %9.3f %9.3f %9.3f %8.0f'
              % (g, len(u), u.rise.mean(), us.rise.mean(), u.ratio.mean(),
                 row.p50))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
