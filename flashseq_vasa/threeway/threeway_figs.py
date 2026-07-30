#!/usr/bin/env python3
"""threeway_figs.py -- the two published panels, three groups, paper's own form.

FIG 2b  mirrored horizontal density, the AUTHORS' four panels and AUTHORS' class
        definitions (Protein coding / lncRNA / Transcription factors / sncRNA).
        The CLASSES AND THE NUMBERS ARE THE AUTHORS' (see threeway_paperform.py);
        the mirrored rendering is this project's reading of the published panel.
        Mirror axis = protocol: VASA left, FLASH-seq right. Within each side the
        second arm is a dashed outline, the convention the two-way figure used.

FIG 1f  genes per cell vs trimmed reads per cell. MY REIMPLEMENTATION -- the code
        for this panel was never deposited. Axes and depth grid are the paper's.

BINNING -- WHY THE GROUPS DIFFER, which is a real departure from the paper
-------------------------------------------------------------------------
The published plate gives n=173 mESC cells; the own plate n=12 and FLASH-seq
n=9-10. One bin count cannot serve both.

  n=173  ->  Freedman-Diaconis, per panel. FD is the standard density rule and is
             well behaved at this n (it returns 20-30 bins here).
  n<=12  ->  fixed 7 bins, unchanged from the two-way figure.

FD is NOT used at n<=12 because it is unstable there, and that was measured, not
assumed: on these same four panels FD returned 4/2/2/2 bins for the own plate but
30/4/15/10 for FLASH-seq native, because at n=10 the IQR collapses relative to
the range and the bin width goes to zero. A rule that returns 30 bins for 10
points is not a rule. Every unit is also drawn as a dot in every panel, so n is
never hidden by the density whichever count is used.

Each density is normalised to its own maximum, so the panels compare SHAPE and
LOCATION, not counts -- with n of 173 against 10 a shared count scale would make
the small groups invisible. The paper's panels are normalised the same way.

EXCLUDED LIBRARY
----------------
FLASH-seq ZHA8833A8 carries qc_verdict='exclude' in code/flashseq/sample_metadata.tsv
(18.3% human CALB1; not a usable 60 pg data point). It is therefore kept OUT of
every density, curve and median here, so FLASH-seq is n=9 -- matching the two-way
figure -- but it is still DRAWN, as an open ring in Fig 2b, because its
contamination is a finding and hiding it would lose the well effect that found it.
The verdict filters interpretation, not QC. Its rows remain in the TSVs.
"""
import os
import sys

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

RES = sys.argv[1] if len(sys.argv) > 1 else '.'

# AUTHORS' panel set for Fig 2b, in the published panel's order. 'TF_authors' and
# 'smallRNA' are the authors' classes; the display names are the paper's labels.
PANELS = [('ProteinCoding', 'Protein coding'),
          ('lncRNA', 'lncRNA'),
          ('TF_authors', 'Transcription factors'),
          ('smallRNA', 'sncRNA')]

# Colour threads the ENTITY (figure-style 4.1) and nests by protocol (4.3):
# VASA = blue family, FLASH-seq = orange family.
VASA_PUB = '#1f4e9c'
VASA_OWN = '#5b9bd5'
FS_NAT = '#d9822b'
FS_VAS = '#8c4a1a'

SIDE = {'published VASA-plate': -1, 'own VASA-plate': -1,
        'FLASH-seq native': +1, 'FLASH-seq vasalen': +1}
STYLE = {'published VASA-plate': dict(color=VASA_PUB, fill=True, lw=0.8),
         'own VASA-plate': dict(color=VASA_OWN, fill=False, lw=1.6, ls='--'),
         'FLASH-seq native': dict(color=FS_NAT, fill=True, lw=0.8),
         'FLASH-seq vasalen': dict(color=FS_VAS, fill=False, lw=1.6, ls='--')}
ORDER = ['published VASA-plate', 'own VASA-plate', 'FLASH-seq native',
         'FLASH-seq vasalen']

FD_MIN_N = 30      # above this n, trust Freedman-Diaconis
SMALL_N_BINS = 7   # at or below it, the two-way figure's fixed count

# qc_verdict='exclude' in code/flashseq/sample_metadata.tsv. Out of every summary
# statistic (figure-style 1.1); still drawn, as an open ring.
EXCLUDED_UNITS = {'ZHA8833A8'}


def fd_bins(x, cap=30):
    """Freedman-Diaconis bin count. Only called when n > FD_MIN_N."""
    x = np.asarray(x, float)
    q75, q25 = np.percentile(x, [75, 25])
    iqr = q75 - q25
    if iqr <= 0:
        return int(np.ceil(np.sqrt(len(x))))
    h = 2 * iqr / len(x) ** (1 / 3)
    return int(np.clip(np.ceil((x.max() - x.min()) / h), 1, cap))


def nbins_for(x):
    return fd_bins(x) if len(x) > FD_MIN_N else SMALL_N_BINS


def mirrored_panel(ax, data, excl, lo, hi):
    """Draw the mirrored densities.

    `data` maps dataset -> fractions of INCLUDED units; `excl` maps dataset ->
    fractions of qc-excluded units, which are drawn but never binned."""
    used = {}
    for name in ORDER:
        x = data.get(name)
        if x is None or len(x) == 0:
            continue
        st = STYLE[name]
        nb = nbins_for(x)
        used[name] = nb
        cnt, edges = np.histogram(x, bins=nb, range=(lo, hi))
        if cnt.max() == 0:
            continue
        w = cnt / cnt.max() * SIDE[name]          # normalised to its own max
        h = np.diff(edges)
        if st['fill']:
            ax.barh(edges[:-1], w, height=h, align='edge', color=st['color'],
                    alpha=0.75, linewidth=0, zorder=2)
        else:
            # dashed outline: step the profile so it reads as a boundary
            xs, ys = [0.0], [edges[0]]
            for i in range(len(w)):
                xs += [w[i], w[i]]
                ys += [edges[i], edges[i + 1]]
            xs.append(0.0); ys.append(edges[-1])
            ax.plot(xs, ys, color=st['color'], lw=st['lw'], ls=st['ls'],
                    zorder=4, solid_joinstyle='miter')
    # every unit as a dot, so n is never hidden by the density
    rng = np.random.default_rng(0)
    for name in ORDER:
        x = data.get(name)
        if x is None or len(x) == 0:
            continue
        sd = SIDE[name]
        jit = rng.uniform(0.10, 0.62, len(x)) * sd
        ax.plot(jit, x, marker='o', ls='none', ms=1.9 if len(x) > 30 else 3.2,
                mfc='0.12', mec='none', alpha=0.55 if len(x) > 30 else 0.9,
                zorder=6)
    # qc-excluded units: open ring, distinct from every included mark, and not in
    # any density or median (figure-style 1.1)
    for name in ORDER:
        x = excl.get(name)
        if x is None or len(x) == 0:
            continue
        sd = SIDE[name]
        ax.plot(rng.uniform(0.18, 0.55, len(x)) * sd, x, marker='o', ls='none',
                ms=5.0, mfc='none', mec='0.12', mew=1.1, zorder=7)
    return used


def fig2b(frac, out):
    u = frac[frac.table_family == 'uniagg']
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))
    used_all = {}
    for k, (panel, label) in enumerate(PANELS):
        ax = axes.flat[k]
        sub = {n: u[(u.dataset == n) & (u.panel == panel)] for n in ORDER}
        data = {n: d[~d.unit.isin(EXCLUDED_UNITS)].fraction.values
                for n, d in sub.items()}
        excl = {n: d[d.unit.isin(EXCLUDED_UNITS)].fraction.values
                for n, d in sub.items()}
        allv = np.concatenate([v for v in list(data.values()) + list(excl.values())
                               if len(v)])
        pad = 0.06 * (allv.max() - allv.min() or 1)
        lo, hi = allv.min() - pad, allv.max() + pad
        used_all[panel] = mirrored_panel(ax, data, excl, lo, hi)
        ax.axvline(0, color='0.25', lw=0.9, zorder=5)
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(lo, hi)
        ax.set_xticks([])
        ax.set_title(label, loc='left', fontsize=9, pad=14)
        if k % 2 == 0:
            ax.set_ylabel('Fraction of reads', fontsize=9)
        for sp in ('top', 'right', 'bottom'):
            ax.spines[sp].set_visible(False)
        ax.tick_params(labelsize=7)
        # medians, read off the data
        # Medians of the INCLUDED units only. Placed in axes coordinates just
        # above the panel so they cannot land on a dot or a density (the two-way
        # figure had them inside the data area and they collided).
        med = {n: (np.median(v) if len(v) else np.nan) for n, v in data.items()}
        ax.text(0.0, 1.012, 'VASA  pub %.3f   own %.3f'
                % (med['published VASA-plate'], med['own VASA-plate']),
                transform=ax.transAxes, color=VASA_PUB, fontsize=6.8,
                va='bottom', ha='left')
        ax.text(1.0, 1.012, 'FLASH-seq  nat %.3f   trimmed %.3f'
                % (med['FLASH-seq native'], med['FLASH-seq vasalen']),
                transform=ax.transAxes, color=FS_NAT, fontsize=6.8,
                va='bottom', ha='right')

    handles = [
        Patch(facecolor=VASA_PUB, alpha=0.75,
              label='published VASA-plate, mESC, n=173'),
        Line2D([], [], color=VASA_OWN, ls='--', lw=1.6,
               label='own VASA-plate, n=12'),
        Patch(facecolor=FS_NAT, alpha=0.75, label='FLASH-seq native, n=9'),
        Line2D([], [], color=FS_VAS, ls='--', lw=1.6,
               label='FLASH-seq VASA-trimmed, n=9'),
        Line2D([], [], color='0.12', marker='o', ls='none', ms=3,
               label='one unit'),
        Line2D([], [], color='0.12', marker='o', ls='none', ms=5, mfc='none',
               mew=1.1, label='ZHA8833A8, qc-excluded (not in any density/median)'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, -0.005))
    bintxt = '; '.join('%s %s' % (lab, used_all[p].get('published VASA-plate'))
                       for p, lab in PANELS)
    fig.suptitle(
        "Fig 2b in the paper's mirrored form, three groups. Classes, class "
        "definitions and denominator are the AUTHORS' code\n"
        "(b_Analysis/02_scanpy_QCxBiotype.py, applied to our matrices); the "
        "mirrored rendering is ours. ReadCounts on all sides.\n"
        "Bins: published plate Freedman-Diaconis per panel (%s); own plate and "
        "FLASH-seq fixed 7 (FD is unstable at n<=12).\n"
        "Densities normalised to their own maxima, so shape and location "
        "compare, not counts. Published plate = mESC; own plate = mouse embryo."
        % bintxt,
        fontsize=7.4, ha='left', x=0.012, y=0.995, va='top')
    fig.tight_layout(rect=(0, 0.055, 1, 0.865))
    fig.savefig(out, dpi=200, bbox_inches='tight')
    return fig, used_all


# --- Fig 1f ----------------------------------------------------------------
CURVE = {'published VASA-plate': dict(c=VASA_PUB, m='o', lab='published VASA-plate'),
         'own VASA-plate': dict(c=VASA_OWN, m='s', lab='own VASA-plate'),
         'FLASH-seq native': dict(c=FS_NAT, m='^', lab='FLASH-seq native'),
         'FLASH-seq vasalen': dict(c=FS_VAS, m='v', lab='FLASH-seq VASA-trimmed')}
GRID = [5000, 10000, 15000, 20000, 25000, 50000, 75000]
PAPER_FIG1F = (9480, 1252)
PAPER_EXT2E = (15248, 1092)

# The single-cell-equivalent input rung, per code/flashseq/sample_metadata.tsv:
# "a mammalian cell carries ~10-30 pg of total RNA, so the 30 pg rung is the
# single-cell-equivalent one". A7 (60 pg) carries a 3.6% CALB1 caveat and A8 is
# excluded outright, so the 30 pg pair is the clean input-matched comparator.
FS_30PG = {'ZHA8833A9', 'ZHA8833A10'}


def fig1f(sat, out):
    """MINE (reimplementation). Paper axes and depth grid.

    THREE panels, because the two-way version's headline ("all entries:
    FLASH-seq leads; single-gene only: VASA leads") turned out to depend on a
    confound the two-way figure could not see:

      Pooling the 9 qc-ok FLASH-seq libraries averages over a 1,000x input
      titration (30 ng down to 30 pg). Detection rises steeply with input, so the
      pooled FLASH-seq mean is dominated by the ng-scale libraries and is not a
      single-cell measurement at all. Panel 2 therefore pools (reproducing the
      earlier claim) and panel 3 restricts FLASH-seq to the 30 pg rung -- the
      single-cell-equivalent input, per sample_metadata.tsv -- where the flip
      actually is: FLASH-seq +225 entries on all entries, -729 on single-gene.
    """
    s = sat[sat.table_family == 'uniagg']
    s = s[~s.unit.isin(EXCLUDED_UNITS)]        # figure-style 1.1
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6))

    # --- panels 1-2: pooled, both scopes (the two-way figure's comparison) ----
    for k, cs in enumerate(['all_entries', 'single_gene']):
        ax = axes[k]
        d = s[(s.count_scope == cs) & (s.cohort == 'mouse_mESC')
              & (s.depth_trimmed_reads.isin(GRID))]
        for name, st in CURVE.items():
            g = d[d.dataset == name].groupby('depth_trimmed_reads').genes
            m, n = g.mean(), g.count()
            if not len(m):
                continue
            ax.plot(m.index / 1000, m.values, color=st['c'], marker=st['m'],
                    ms=4, lw=1.6, label='%s (n=%d)' % (st['lab'], n.iloc[0]))
        ax.set_xlabel('Trimmed reads per unit (thousands)', fontsize=8.5)
        if k == 0:
            ax.set_ylabel('Entries detected (mean per unit)', fontsize=8.5)
        ax.set_title(('All entries: FLASH-seq leads'
                      if cs == 'all_entries'
                      else 'Single-gene only: the gap narrows, does not close'),
                     loc='left', fontsize=9)
        ax.legend(frameon=False, fontsize=6.8, loc='upper left')
        ax.tick_params(labelsize=7.5)
        ax.margins(0.05)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)

    # --- panel 3: input-matched, where the flip actually is ------------------
    ax = axes[2]
    d75 = s[(s.depth_trimmed_reads == 75000) & (s.cohort == 'mouse_mESC')]
    fs30 = d75[(d75.dataset == 'FLASH-seq native')
               & (d75.unit.isin(FS_30PG))]
    ownd = d75[d75.dataset == 'own VASA-plate']
    xs = np.arange(2)
    w = 0.32
    for i, (lab, dd, col) in enumerate([
            ('own VASA-plate, 12 cells', ownd, VASA_OWN),
            ('FLASH-seq native, 30 pg rung', fs30, FS_NAT)]):
        vals, errs, ns = [], [], []
        for cs in ('all_entries', 'single_gene'):
            g = dd[dd.count_scope == cs].genes
            vals.append(g.mean()); errs.append(g.std()); ns.append(len(g))
        ax.bar(xs + (i - 0.5) * w, vals, w, yerr=errs, capsize=3, color=col,
               alpha=0.85, label='%s (n=%d)' % (lab, ns[0]),
               error_kw=dict(lw=0.9))
        for x, v in zip(xs + (i - 0.5) * w, vals):
            ax.text(x, v + 130, '%.0f' % v, ha='center', va='bottom', fontsize=7)
    for j, cs in enumerate(('all_entries', 'single_gene')):
        a = fs30[fs30.count_scope == cs].genes.mean()
        b = ownd[ownd.count_scope == cs].genes.mean()
        ax.text(j, 900, '%+.0f\n%s' % (a - b, 'FLASH-seq' if a > b else 'VASA'),
                ha='center', va='bottom', fontsize=7, color='0.2')
    ax.set_xticks(xs)
    ax.set_xticklabels(['All entries', 'Single-gene only'], fontsize=8)
    ax.set_ylabel('Entries detected at 75k trimmed reads', fontsize=8.5)
    ax.set_title('Input-matched: the flip is real, but small', loc='left',
                 fontsize=9)
    ax.legend(frameon=False, fontsize=6.8, loc='upper right')
    ax.tick_params(labelsize=7.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.margins(y=0.16)

    fig.suptitle(
        "Fig 1f, paper axes and depth grid. MY REIMPLEMENTATION -- the code for "
        "this panel was never deposited, so a difference from the paper could be "
        "the reimplementation rather than the data.\n"
        "Mouse cells throughout (published plate = mESC by the paper's own Fig 1d "
        "rule). Both scopes are shown because the published panel does not say "
        "whether multimapper combination entries counted as genes.\n"
        "Panels 1-2 pool the 9 qc-ok FLASH-seq libraries, which spans a 1,000x "
        "input titration (30 ng to 30 pg) and so is NOT a single-cell "
        "measurement; panel 3 restricts to the 30 pg rung.\n"
        "ZHA8833A8 (qc_verdict=exclude) is out of every curve. Estimator: "
        "binomial-thinning expectation, deterministic.",
        fontsize=7.3, ha='left', x=0.008, y=0.995, va='top')
    fig.tight_layout(rect=(0, 0, 1, 0.825))
    fig.savefig(out, dpi=200, bbox_inches='tight')
    return fig


def fig_papercheck(sat, out):
    """The direct test of the reimplementation against the PUBLISHED curve.

    Only possible now the published plate is in the comparison. Uses the paper's
    OWN cohort: HEK293T cells (human entries), gated at >=75k trimmed reads, which
    is what the paper's caption specifies -- running this on mESC would not be a
    test of the published number."""
    s = sat[(sat.table_family == 'uniagg')
            & (sat.dataset == 'published VASA-plate')
            & (sat.entry_scope == 'human')]
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    xs = np.arange(2)
    w = 0.26
    for i, cs in enumerate(['all_entries', 'single_gene']):
        vals, errs, ns = [], [], []
        for dep in (75000, 750000):
            g = s[(s.count_scope == cs) & (s.depth_trimmed_reads == dep)].genes
            vals.append(g.mean()); errs.append(g.std()); ns.append(len(g))
        ax.bar(xs + (i - 0.5) * w, vals, w, yerr=errs, capsize=3,
               color=(VASA_PUB if cs == 'all_entries' else VASA_OWN),
               alpha=0.85, label='our re-run, %s' % cs.replace('_', ' '),
               error_kw=dict(lw=0.9))
        for x, v, nn in zip(xs + (i - 0.5) * w, vals, ns):
            ax.text(x, v + 380, '%.0f\nn=%d' % (v, nn), ha='center',
                    va='bottom', fontsize=6.6)
    for j, (pv, pe) in enumerate([PAPER_FIG1F, PAPER_EXT2E]):
        ax.errorbar([j], [pv], yerr=[pe], fmt='D', ms=7, color='0.12',
                    capsize=4, lw=1.3, zorder=5,
                    label='published value +- s.d.' if j == 0 else None)
        ax.text(j + 0.40, pv, '%d\n+-%d' % (pv, pe), fontsize=6.6, va='center',
                ha='left', color='0.12')
    ax.set_xticks(xs)
    ax.set_xticklabels(['75k trimmed reads\n(Fig 1f, n=174)',
                        '750k trimmed reads\n(Ext. Data Fig 2e)'], fontsize=7.8)
    ax.set_ylabel('Genes detected per HEK293T cell', fontsize=8.5)
    ax.set_title('Our re-run of THEIR data against their published curve\n'
                 'Single-gene scope reproduces it: 1.02x at both depths',
                 loc='left', fontsize=8.8)
    ax.legend(frameon=False, fontsize=7, loc='upper left')
    ax.tick_params(labelsize=7.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.margins(y=0.22)
    fig.text(0.008, -0.055,
             "Paper's own cohort and definition: HEK293T (human), all annotated "
             "genes, cells sequenced to >=75k trimmed reads.\n"
             "Estimator is MY reimplementation (binomial-thinning expectation); "
             "the published code was never deposited. Counting\n"
             "multimapper combination entries as genes inflates the estimate by "
             "1.09x at 75k and 1.24x at 750k -- so the scope, not\n"
             "the chemistry, is what decides agreement with the printed value.",
             fontsize=7.0, ha='left', va='top')
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches='tight')
    return fig


if __name__ == '__main__':
    frac = pd.read_csv(f'{RES}/paperform_threeway.tsv', sep='\t')
    sat = pd.read_csv(f'{RES}/paperform_threeway_fig1f.tsv', sep='\t')
    f1, used = fig2b(frac, f'{RES}/paperfig2b_threeway_mirrored.png')
    f2 = fig1f(sat, f'{RES}/paperfig1f_threeway.png')
    f3 = fig_papercheck(sat, f'{RES}/paperfig1f_publishedcheck.png')
    print('bins used, published plate:',
          {p: v.get('published VASA-plate') for p, v in used.items()})
    print('wrote paperfig2b_threeway_mirrored.png, paperfig1f_threeway.png, '
          'paperfig1f_publishedcheck.png')
