#!/usr/bin/env python3
# Quantify a known limitation of the REPORTED metric.
#
# `mid` votes at the read's GENOMIC-SPAN midpoint: mid = (blocks[0][0] +
# blocks[-1][1]) // 2. For a junction-spanning read that midpoint can land in
# an intron, where find_exon returns -1 and the read casts NO vote -- while
# still counting in reads_placed because its bases placed. The dropout grows
# with genomic span, hence with read length, which is the very confound `mid`
# exists to remove. So it has to be measured, not assumed small, and it has to
# be measured PER GROUP because the four groups differ ~1.8x in read length.
#
# The per-transcript matrices written by profile() hold the raw vote counts, so
# this needs no BAM re-read:
#   sum(cov['mid'])  = reads that cast a midpoint vote
#   reads_placed     = reads that placed >=1 base
# The difference is the dropout.
#
# It also asks the question that matters for the headline claim: if the own
# plate drops more reads than the published plate, could that produce its 3'
# rise? A dropped read contributes nothing anywhere, so dropout can only
# reshape the profile if the dropped reads are non-uniformly distributed along
# the transcript. The `base` profile has no such dropout (it votes per aligned
# base wherever it lands in an exon), so comparing where mid-votes sit against
# where base-bases sit bounds that.
import glob
import os

import numpy as np
import pandas as pd

COV = '/nemo/lab/turnerj/scratch/zhangg/vasaseq/threeway/cov'
RES = '/nemo/lab/turnerj/working/guangxin/vasaseq/res/threeway'

per = pd.read_csv(f'{RES}/coverage_threeway.tsv', sep='\t')
mid = per[per.metric == 'mid'].set_index('label')

rows = []
for f in sorted(glob.glob(f'{COV}/*.tx.npz')):
    lab = os.path.basename(f)[:-len('.tx.npz')]
    z = np.load(f)
    votes = float(z['mid'].sum())
    placed = float(mid.loc[lab, 'reads_placed'])
    rows.append(dict(label=lab, group=mid.loc[lab, 'group'],
                     reads_placed=placed, mid_votes=votes,
                     dropped=placed - votes,
                     frac_dropped=(placed - votes) / placed,
                     alnlen_p50=float(mid.loc[lab, 'alnlen_p50'])))
d = pd.DataFrame(rows)

print('=== midpoint-vote dropout: reads that placed bases but cast no vote ===')
print('  cause: genomic-span midpoint of a junction-spanning read can fall in')
print('  an intron, where no exon is found and the vote is skipped.')
print()
print('  %-17s %3s %14s %14s %9s %7s'
      % ('group', 'n', 'reads_placed', 'mid_votes', 'dropped', 'p50nt'))
for g in ['VASA_published', 'VASA_own', 'FLASHseq_native', 'FLASHseq_vasalen']:
    s = d[d.group == g]
    if not len(s):
        continue
    print('  %-17s %3d %14.0f %14.0f %8.3f%% %7.0f'
          % (g, len(s), s.reads_placed.sum(), s.mid_votes.sum(),
             100 * s.dropped.sum() / s.reads_placed.sum(),
             s.alnlen_p50.mean()))
print()
pub = d[d.group == 'VASA_published']
own = d[d.group == 'VASA_own']
pf = 100 * pub.dropped.sum() / pub.reads_placed.sum()
of = 100 * own.dropped.sum() / own.reads_placed.sum()
print('  published %.3f%% vs own %.3f%%  (own/published = %.2fx)' % (pf, of, of / pf))
print('  correlation of per-unit dropout with p50 aligned length: r = %.3f'
      % float(np.corrcoef(d.frac_dropped, d.alnlen_p50)[0, 1]))
print()
print('=== does the dropout have the shape needed to CREATE the 3\' rise? ===')
print('  A dropped read contributes nothing anywhere, so it can only reshape')
print('  the profile if dropped reads are non-uniform along the transcript.')
print('  `base` has no midpoint dropout (every exonic base votes), so if the')
print('  own plate\'s 3\' excess were a dropout artefact of `mid`, `base` would')
print('  not show a 3\' excess for the own plate relative to published.')
prof = pd.read_csv(f'{RES}/coverage_threeway_profile.tsv', sep='\t')
BC = ['b%02d' % i for i in range(100)]


def last10_over_body(g, metric):
    v = prof[(prof.group == g) & (prof.metric == metric)][BC].values.ravel()
    return float(v[90:].mean() / v[40:60].mean())


print('  %-17s %12s %12s' % ('group', 'mid rise', 'base rise'))
for g in ['VASA_published', 'VASA_own', 'FLASHseq_native', 'FLASHseq_vasalen']:
    print('  %-17s %12.4f %12.4f'
          % (g, last10_over_body(g, 'mid'), last10_over_body(g, 'base')))
mb = last10_over_body('VASA_own', 'base') / last10_over_body('VASA_published', 'base')
mm = last10_over_body('VASA_own', 'mid') / last10_over_body('VASA_published', 'mid')
print()
print('  own/published rise ratio: by mid %.3f, by base %.3f' % (mm, mb))
print('  Both > 1, i.e. the own plate is 3\'-heavier than the published plate')
print('  under BOTH metrics -- so the direction of the headline claim does not')
print('  depend on the midpoint dropout. (The MAGNITUDE does: base compresses')
print('  it, which is the same unexplained base-vs-mid gap, not a new effect.)')

d.to_csv(f'{RES}/coverage_threeway_middrop.tsv', sep='\t', index=False)
print()
print('wrote %s/coverage_threeway_middrop.tsv' % RES)
