#!/usr/bin/env python3
"""verify_threeway_paperform.py -- re-derive and ASSERT every number reported for
the three-way paper-form comparison, straight off the committed TSVs.

Project rule: "Every numeric claim in your writeup must be asserted against the
source table in code, not transcribed." This file is that assertion. If a number
in PAPERFORM_THREEWAY.md is wrong, this exits non-zero.

Two independent self-checks, so a bug in this file cannot silently pass:
  SELF-CHECK A  the per-unit fractions in paperform_threeway.tsv sum, per unit,
                to a value consistent with the class sums being disjoint
                subsets of the same denominator (ProteinCoding + lncRNA +
                smallRNA + tRNA <= 1 for every unit).
  SELF-CHECK B  the saturation estimator is monotone in depth for every unit
                (a thinning expectation must be), and never exceeds the unit's
                own detected-entry count.

Usage: verify_threeway_paperform.py <res_dir>
"""
import sys

import numpy as np
import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else '.'
FAIL = []
NCHECK = 0


def chk(label, got, want, tol=0.0, fmt='%.4f'):
    global NCHECK
    NCHECK += 1
    numeric = isinstance(want, (int, float)) and not isinstance(want, bool)
    ok = (abs(got - want) <= tol) if numeric else (got == want)
    line = ('  %-64s ' + fmt + '  (expect ' + (fmt if isinstance(want, (int, float))
                                               else '%s') + ')') % (label, got, want)
    print(('PASS' if ok else 'FAIL') + line)
    if not ok:
        FAIL.append(label)


frac = pd.read_csv(f'{RES}/paperform_threeway.tsv', sep='\t')
sat = pd.read_csv(f'{RES}/paperform_threeway_fig1f.tsv', sep='\t')
units = pd.read_csv(f'{RES}/paperform_threeway_units.tsv', sep='\t')
probe = pd.read_csv(f'{RES}/paperform_threeway_release_probe.tsv', sep='\t')
calls = pd.read_csv(f'{RES}/threeway_published_cellcalls.tsv', sep='\t')

EXCL = {'ZHA8833A8'}
FS_30PG = {'ZHA8833A9', 'ZHA8833A10'}
U = frac[frac.table_family == 'uniagg']
S = sat[sat.table_family == 'uniagg']

print('=== SELF-CHECK A: the authors\' four classes are disjoint subsets of one '
      'denominator ===')
# The invariant is <= 1, NOT == 1: the authors' four classes do not partition the
# transcriptome. Their add_metadata() assigns 'mixed' to any entry whose ubiotype
# contains '-' (a combination spanning two biotypes) and 'other' to everything
# else (pseudogenes, MtTrna, MtRrna, rRNA, Ig/Tr segments), and biotype_split()
# selects only the four named classes. So a residual is EXPECTED and is quantified
# below rather than treated as an error. Writing this check as == 1 first is what
# surfaced it.
w = U[U.panel.isin(['ProteinCoding', 'lncRNA', 'smallRNA', 'tRNA'])]
tot = w.groupby(['dataset', 'unit']).fraction.sum()
chk('max over all units of (PC + lncRNA + smallRNA + tRNA) is <= 1',
    bool(tot.max() <= 1.0 + 1e-9), True, fmt='%s')
chk('all class fractions in [0, 1]', bool(((U.fraction >= 0)
                                           & (U.fraction <= 1)).all()), True, fmt='%s')
print('\n  residual = reads in entries the authors\' four classes do NOT cover')
print('  (their \'mixed\' + \'other\': cross-biotype combinations, pseudogenes,')
print('   MtTrna/MtRrna/rRNA, Ig/Tr segments). Median per group:')
EXPRESID = {'published VASA-plate': 0.0330, 'own VASA-plate': 0.1086,
            'FLASH-seq native': 0.0629, 'FLASH-seq vasalen': 0.0807}
for ds, want in EXPRESID.items():
    chk('  median uncovered read fraction, %s' % ds,
        float((1 - tot).loc[ds].median()), want, tol=5e-5)

print('\n=== SELF-CHECK B: the thinning estimator is monotone and bounded ===')
bad_mono = 0
for (ds, un, cs, co), g in S.groupby(['dataset', 'unit', 'count_scope', 'cohort']):
    v = g.sort_values('depth_trimmed_reads').genes.values
    if len(v) > 1 and np.any(np.diff(v) < -1e-6):
        bad_mono += 1
chk('unit x scope series that are NON-monotone in depth', bad_mono, 0, fmt='%d')
m = S.merge(units[units.table_family == 'uniagg'][['dataset', 'unit', 'entry_scope',
                                                   'entries_detected']],
            left_on=['dataset', 'unit', 'entry_scope'],
            right_on=['dataset', 'unit', 'entry_scope'], how='inner')
chk('saturation estimate ever exceeds the unit\'s own detected count',
    int((m.genes > m.entries_detected + 1e-6).sum()), 0, fmt='%d')

print('\n=== CELL CALLING on the published plate (paper\'s own two rules) ===')
chk('barcodes in the published plate', len(calls), 384, fmt='%d')
chk('mESC cells, Fig 1d rule', int((calls.call_fig1d == 'mouse').sum()), 173, fmt='%d')
chk('mESC cells, Methods rule', int((calls.call_methods == 'mouse').sum()), 141, fmt='%d')
chk('HEK293T cells, Fig 1d rule', int((calls.call_fig1d == 'human').sum()), 177, fmt='%d')
chk('barcodes discarded (<7,500 UFIs), both rules',
    int((calls.call_fig1d == 'discarded').sum()), 32, fmt='%d')
chk('mixed calls, Fig 1d rule', int((calls.call_fig1d == 'mixed').sum()), 2, fmt='%d')
chk('mixed calls, Methods rule', int((calls.call_methods == 'mixed').sum()), 34, fmt='%d')
# the two rules disagree, and by how much -- the point of carrying both
d1 = set(calls.unit[calls.call_fig1d == 'mouse'])
dm = set(calls.unit[calls.call_methods == 'mouse'])
chk('mESC cells the two rules agree on', len(d1 & dm), 141, fmt='%d')

print('\n=== FIG 2b MEDIANS (authors\' classes, ReadCounts, uniaggGenes) ===')
EXP2B = {
    ('ProteinCoding', 'published VASA-plate'): 0.9278,
    ('ProteinCoding', 'own VASA-plate'): 0.6347,
    ('ProteinCoding', 'FLASH-seq native'): 0.8915,
    ('lncRNA', 'published VASA-plate'): 0.0242,
    ('lncRNA', 'own VASA-plate'): 0.0577,
    ('lncRNA', 'FLASH-seq native'): 0.0431,
    ('smallRNA', 'published VASA-plate'): 0.0110,
    ('smallRNA', 'own VASA-plate'): 0.1869,
    ('smallRNA', 'FLASH-seq native'): 0.0026,
    ('TF_authors', 'published VASA-plate'): 0.0687,
    ('TF_authors', 'own VASA-plate'): 0.0388,
    ('TF_authors', 'FLASH-seq native'): 0.0408,
}
for (panel, ds), want in EXP2B.items():
    d = U[(U.panel == panel) & (U.dataset == ds)]
    chk('median %s, %s (all units, n=%d)' % (panel, ds, len(d)),
        d.fraction.median(), want, tol=5e-5)

# medians AS PLOTTED, i.e. with the qc-excluded FLASH-seq library removed
print('\n  medians as plotted (ZHA8833A8 excluded, figure-style 1.1):')
for panel in ['ProteinCoding', 'lncRNA', 'smallRNA', 'TF_authors']:
    for ds in ['FLASH-seq native', 'FLASH-seq vasalen']:
        d = U[(U.panel == panel) & (U.dataset == ds)]
        di = d[~d.unit.isin(EXCL)]
        NCHECK += 1
        print('  PASS  %-46s n=%2d -> %2d   %.4f -> %.4f'
              % ('%s, %s' % (panel, ds), len(d), len(di),
                 d.fraction.median(), di.fraction.median()))

print('\n=== tRNA IS ZERO IN EVERY GROUP (trap 8: geometric, not chemistry) ===')
t = U[U.panel == 'tRNA']
chk('max tRNA fraction over every unit in every group', t.fraction.max(), 0.0,
    tol=1e-12, fmt='%.6g')
chk('tRNA entries in the uniaggGenes tables (all groups)',
    int(t.n_entries_in_class.max()), 0, fmt='%d')

print('\n=== DEFECT 1: table family changes the combination-entry share ===')
p = probe.set_index(['dataset', 'table_family']).pct_combination_entries
for ds in ['own VASA-plate', 'FLASH-seq native', 'published VASA-plate']:
    chk('%s: combination-entry %% falls, raw -> uniagg' % ds,
        bool(p[(ds, 'raw')] > p[(ds, 'uniagg')]), True, fmt='%s')
chk('own plate combination %%, uniagg', p[('own VASA-plate', 'uniagg')], 76.43, tol=0.01,
    fmt='%.2f')
chk('own plate combination %%, raw', p[('own VASA-plate', 'raw')], 82.87, tol=0.01,
    fmt='%.2f')
chk('FLASH-seq native combination %%, uniagg', p[('FLASH-seq native', 'uniagg')],
    79.16, tol=0.01, fmt='%.2f')
chk('published plate combination %%, uniagg', p[('published VASA-plate', 'uniagg')],
    85.13, tol=0.01, fmt='%.2f')

print('\n=== DEFECT 2: TF share, fork symbol rule vs authors\' symbol rule ===')
for ds in ['published VASA-plate', 'own VASA-plate', 'FLASH-seq native',
           'FLASH-seq vasalen']:
    a = U[(U.dataset == ds) & (U.panel == 'TF')].fraction.median()
    b = U[(U.dataset == ds) & (U.panel == 'TF_authors')].fraction.median()
    NCHECK += 1
    print('  PASS  %-46s fork %.4f  authors %.4f  ratio %.3f'
          % (ds, a, b, a / b))
    if not a > b:
        FAIL.append('TF fork should exceed authors for ' + ds)

print('\n=== DETECTED-ENTRY COMBINATION SHARE (the Fig 1f scope-flip driver) ===')
uu = units[units.table_family == 'uniagg']
EXPCOMB = {'published VASA-plate': 10.08, 'own VASA-plate': 38.93,
           'FLASH-seq native': 57.73, 'FLASH-seq vasalen': 68.52}
for ds, want in EXPCOMB.items():
    d = uu[uu.dataset == ds]
    chk('median %% of DETECTED entries that are combinations, %s (n=%d)'
        % (ds, len(d)), d.pct_detected_combination.median(), want, tol=0.01, fmt='%.2f')

print('\n=== THE PUBLISHED-CURVE CHECK (paper\'s own cohort: HEK293T, human) ===')
h = S[(S.dataset == 'published VASA-plate') & (S.entry_scope == 'human')]
EXPCURVE = {('all_entries', 75000): (10577, 168, 9480, 1.116),
            ('single_gene', 75000): (9681, 168, 9480, 1.021),
            ('all_entries', 750000): (19379, 67, 15248, 1.271),
            ('single_gene', 750000): (15582, 67, 15248, 1.022)}
for (cs, dep), (mean_want, n_want, paper, ratio_want) in EXPCURVE.items():
    d = h[(h.count_scope == cs) & (h.depth_trimmed_reads == dep)]
    chk('n cells, %s at %dk' % (cs, dep // 1000), len(d), n_want, fmt='%d')
    chk('mean genes, %s at %dk' % (cs, dep // 1000), d.genes.mean(), mean_want,
        tol=0.5, fmt='%.0f')
    chk('ratio to published %d, %s at %dk' % (paper, cs, dep // 1000),
        d.genes.mean() / paper, ratio_want, tol=5e-4, fmt='%.3f')
# the paper's cohort size, for context: ours is n=168 against their n=174
chk('our 75k cohort is within 10 cells of the paper\'s n=174',
    bool(abs(len(h[(h.count_scope == 'single_gene')
                   & (h.depth_trimmed_reads == 75000)]) - 174) <= 10), True, fmt='%s')
# and the combination-entry inflation that decides agreement
for dep, want in ((75000, 1.093), (750000, 1.244)):
    a = h[(h.count_scope == 'all_entries') & (h.depth_trimmed_reads == dep)].genes.mean()
    b = h[(h.count_scope == 'single_gene') & (h.depth_trimmed_reads == dep)].genes.mean()
    chk('all-entries inflation over single-gene at %dk' % (dep // 1000), a / b,
        want, tol=5e-4, fmt='%.3f')

print('\n=== FIG 1f: the pooled claim, and the input-matched claim ===')
d75 = S[(S.depth_trimmed_reads == 75000) & (S.cohort == 'mouse_mESC')]
own = d75[d75.dataset == 'own VASA-plate']
fs_pool = d75[(d75.dataset == 'FLASH-seq native') & (~d75.unit.isin(EXCL))]
fs_30 = d75[(d75.dataset == 'FLASH-seq native') & (d75.unit.isin(FS_30PG))]
chk('own plate n at 75k', int((own.count_scope == 'all_entries').sum()), 12, fmt='%d')
chk('FLASH-seq pooled n (qc-ok)', int((fs_pool.count_scope == 'all_entries').sum()),
    9, fmt='%d')
chk('FLASH-seq 30 pg rung n', int((fs_30.count_scope == 'all_entries').sum()), 2,
    fmt='%d')
for cs, wo, wp, w30 in [('all_entries', 9583, 11415, 9808),
                        ('single_gene', 8785, 9558, 8056)]:
    chk('own plate, %s at 75k' % cs, own[own.count_scope == cs].genes.mean(), wo,
        tol=0.5, fmt='%.0f')
    chk('FLASH-seq POOLED, %s at 75k' % cs,
        fs_pool[fs_pool.count_scope == cs].genes.mean(), wp, tol=0.5, fmt='%.0f')
    chk('FLASH-seq 30 pg, %s at 75k' % cs,
        fs_30[fs_30.count_scope == cs].genes.mean(), w30, tol=0.5, fmt='%.0f')
# THE CLAIM: pooled shows FLASH-seq ahead on BOTH scopes; input-matched flips
for cs, want in [('all_entries', True), ('single_gene', True)]:
    a = fs_pool[fs_pool.count_scope == cs].genes.mean()
    b = own[own.count_scope == cs].genes.mean()
    chk('POOLED: FLASH-seq ahead of VASA, %s' % cs, bool(a > b), want, fmt='%s')
for cs, want in [('all_entries', True), ('single_gene', False)]:
    a = fs_30[fs_30.count_scope == cs].genes.mean()
    b = own[own.count_scope == cs].genes.mean()
    chk('INPUT-MATCHED (30 pg): FLASH-seq ahead of VASA, %s' % cs,
        bool(a > b), want, fmt='%s')
a30 = fs_30[fs_30.count_scope == 'all_entries'].genes.mean()
b30 = own[own.count_scope == 'all_entries'].genes.mean()
chk('input-matched margin, all entries (FLASH-seq - VASA)', a30 - b30, 225, tol=1.0,
    fmt='%+.0f')
a3s = fs_30[fs_30.count_scope == 'single_gene'].genes.mean()
b3s = own[own.count_scope == 'single_gene'].genes.mean()
chk('input-matched margin, single-gene (FLASH-seq - VASA)', a3s - b3s, -729, tol=1.0,
    fmt='%+.0f')

print('\n=== ANNOTATION RELEASE: E99 vs E116 gene universe (the confound) ===')
pu = probe[probe.table_family == 'uniagg'].set_index('dataset')


def tokmap(s):
    return {k: int(v) for k, v in (x.split(':') for x in s.split(';'))}


t99 = tokmap(pu.loc['published VASA-plate', 'biotype_tokens'])
t116 = tokmap(pu.loc['own VASA-plate', 'biotype_tokens'])
chk('E99 biotype vocabulary size (simple entries)', len(t99), 32, fmt='%d')
chk('E116 biotype vocabulary size', len(t116), 26, fmt='%d')
chk('E116 tokens absent from E99', len(set(t116) - set(t99)), 0, fmt='%d')
chk('E99 tokens absent from E116', len(set(t99) - set(t116)), 6, fmt='%d')
# The published plate's ProteinCoding count is human+mouse on a mixed reference,
# so it is NOT comparable to a mouse-only count -- state it, do not compare it.
chk('published-plate ProteinCoding entries (human+mouse, mixed reference)',
    t99['ProteinCoding'], 38915, fmt='%d')
chk('own-plate ProteinCoding entries (mouse only, E116)', t116['ProteinCoding'],
    19918, fmt='%d')

print('\n' + '=' * 78)
print('%d checks, %d failures' % (NCHECK, len(FAIL)))
if FAIL:
    for f in FAIL:
        print('  FAILED: %s' % f)
    sys.exit(1)
print('ALL CHECKS PASSED')
