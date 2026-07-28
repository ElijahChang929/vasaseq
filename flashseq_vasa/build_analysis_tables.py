#!/usr/bin/env python3
"""
build_analysis_tables.py -- turn one arm's raw step-7 output into the tables the
comparison actually reads, applying the analysis-set filters and recording
exactly which of VASA's filters apply here and which do not.

WHAT THIS DOES, AND WHAT IT DELIBERATELY DOES NOT
-------------------------------------------------
It does three things:

  1. Strips the 'cells/' prefix step 6 leaves on every column name under
     CELLID_FROM=f, so columns are bare library ids (ZHA8833A1 ...) and are
     ordered by input amount rather than lexically (A10 sorts between A1 and A2).
  2. Applies the filters that apply, and says which of VASA's do not.
  3. Writes gene metadata for the filtering it does NOT do, in the same shape as
     the VASA side's gene_metadata.tsv, so downstream choices stay the user's.

It does NOT drop libraries. A8 is qc_verdict=exclude (18.3% human CALB1, a well
effect at H:1) and A7 is caveat (3.6%). Both are carried through with the verdict
in a metadata column: the verdict filters interpretation, not QC, and hiding A8
would lose the well effect that identified it. This mirrors the VASA side, where
the four blank barcodes were dropped only because they are confirmed blanks, not
because they looked wrong.

HOW THE FILTERS DIFFER FROM THE VASA SIDE -- verified, not assumed
------------------------------------------------------------------
VASA's analysis set (data/PM26037/out/analysis/) applied three filters:

  VASA filter 1, the unnamed row  -> APPLIES HERE, checked. reduceGeneName can
      produce an empty label, which pandas writes as an unnamed index entry. On
      the VASA side there was exactly 1, carrying 42 UFIs. Whether it appears on
      this path is an empirical question; this script counts it and reports.

  VASA filter 2, the UMI ceiling  -> DOES NOT APPLY, and cannot. bc2trans clamps
      at x >= K; K = 4**len(umi) and the only UMI on this path is the literal
      'A', so K = 4, while UFI per gene per column is at most 1. The clamp is
      unreachable. VASA dropped 8 genes at its K = 4096 ceiling; there is no
      analogue. This script asserts the precondition (max UFI <= 1) rather than
      just asserting the conclusion, and records the result.

  VASA filter 3, blank barcodes   -> DOES NOT APPLY. There are no barcodes and
      no blanks: each column is a whole library. The qc_verdict field plays a
      different role and is carried as metadata, not applied.

One filter is NEW here and is a consequence of the arm design:

  All-zero rows. Step 6 builds the row index from the union over libraries, so a
  gene detected in no library can exist in the frame. Dropped, and counted.

WHAT IS WRITTEN
---------------
  <tables>/FS_<arm>_analysis_total.ReadCounts.tsv        genes x libraries
  <tables>/FS_<arm>_analysis_spliced.ReadCounts.tsv
  <tables>/FS_<arm>_analysis_unspliced.ReadCounts.tsv
  <tables>/FS_<arm>_analysis_uniagg_total.ReadCounts.tsv the aggregated uni-gene set
  <tables>/FS_<arm>_analysis_tRNA.ReadCounts.tsv
  <tables>/FS_<arm>_gene_metadata.tsv                    per-row biotype etc.
  <tables>/FS_<arm>_library_metadata.tsv                 per-column, with qc_verdict
  <tables>/manifest.json                                 which file is which arm

ReadCounts ONLY. UFICounts is a 0/1 detection mask on this path and
TranscriptCounts equals it elementwise, so neither carries abundance
information. The detection mask is preserved as a boolean derived from
ReadCounts > 0 -- identical information, one file instead of two, and no reader
can mistake it for a molecule count.

Usage: build_analysis_tables.py <outdir> <sample> <arm> <tablesdir> <report.txt>
"""
import json
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

# Titration order, from code/flashseq/sample_metadata.tsv. A10 sorts between A1
# and A2 lexically, which is why the order is written out rather than sorted.
LIB_ORDER = ['ZHA8833A1', 'ZHA8833A2', 'ZHA8833A3', 'ZHA8833A4', 'ZHA8833A5',
             'ZHA8833A6', 'ZHA8833A7', 'ZHA8833A8', 'ZHA8833A9', 'ZHA8833A10']
LIB_META = {
    'ZHA8833A1':  ('30 ng',  30000, 1, 'A:1', 'ok',      ''),
    'ZHA8833A2':  ('30 ng',  30000, 2, 'B:1', 'ok',      ''),
    'ZHA8833A3':  ('3 ng',    3000, 1, 'C:1', 'ok',      ''),
    'ZHA8833A4':  ('3 ng',    3000, 2, 'D:1', 'ok',      ''),
    'ZHA8833A5':  ('1.5 ng',  1500, 1, 'E:1', 'ok',      ''),
    'ZHA8833A6':  ('1.5 ng',  1500, 2, 'F:1', 'ok',      ''),
    'ZHA8833A7':  ('60 pg',     60, 1, 'G:1', 'caveat',
                   '3.6% human CALB1 (well G:1); reads are identifiable and can be filtered'),
    'ZHA8833A8':  ('60 pg',     60, 2, 'H:1', 'exclude',
                   '18.3% human CALB1 (well H:1); not a usable 60 pg data point'),
    'ZHA8833A9':  ('30 pg',     30, 1, 'A:2', 'ok',      ''),
    'ZHA8833A10': ('30 pg',     30, 2, 'B:2', 'ok',      ''),
}


def biotype_of(idx):
    idx = str(idx)
    parts = idx.rsplit('-')
    bs = sorted(set(p.rsplit('_')[-1] for p in parts))
    return bs[0] if len(bs) == 1 else 'MULTI:' + '-'.join(bs)


def symbol_of(idx):
    """Gene symbol from an ENSID_Symbol_Biotype label; '' if not that shape."""
    idx = str(idx)
    if '-' in idx:
        return ''
    p = idx.rsplit('_')
    return p[1] if len(p) >= 3 else ''


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    outdir, sample, arm, tablesdir, reportpath = sys.argv[1:6]
    os.makedirs(tablesdir, exist_ok=True)

    out = []

    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 78)
        emit(t)
        emit('=' * 78)

    def rd(suffix):
        p = os.path.join(outdir, '%s%s' % (sample, suffix))
        return pd.read_csv(p, sep='\t', index_col=0) if os.path.exists(p) else None

    section('PROVENANCE')
    emit('generated   : %s' % datetime.now().isoformat(timespec='seconds'))
    emit('script      : %s' % os.path.abspath(__file__))
    emit('arm         : %s' % arm)
    emit('sample      : %s' % sample)
    emit('source      : %s' % outdir)
    emit('output      : %s' % tablesdir)
    emit('pandas %s / numpy %s / python %s'
         % (pd.__version__, np.__version__, sys.version.split()[0]))

    section('LOAD')
    want = {
        'total.ReadCounts': '_total.ReadCounts.tsv',
        'total.UFICounts': '_total.UFICounts.tsv',
        'total.TranscriptCounts': '_total.TranscriptCounts.tsv',
        'uniagg_total.ReadCounts': '_uniaggGenes_total.ReadCounts.tsv',
        'uniagg_spliced.ReadCounts': '_uniaggGenes_spliced.ReadCounts.tsv',
        'uniagg_unspliced.ReadCounts': '_uniaggGenes_unspliced.ReadCounts.tsv',
        'tRNA.ReadCounts': '_tRNA.ReadCounts.tsv',
    }
    tabs = {}
    for name, suf in want.items():
        t = rd(suf)
        tabs[name] = t
        emit('  %-30s %s' % (name, (t.shape,) if t is not None else 'MISSING'))
    if tabs['total.ReadCounts'] is None:
        sys.exit('no _total.ReadCounts.tsv -- step 7 did not complete')

    # ------------------------------------------------------------- columns
    section("COLUMNS -- strip step 6's 'cells/' prefix and order by input amount")
    raw_cols = list(tabs['total.ReadCounts'].columns)
    emit('as written by step 6 : %s' % raw_cols)
    ren = {c: str(c).rsplit('/')[-1] for c in raw_cols}
    libs_present = [ren[c] for c in raw_cols]
    unknown = [l for l in libs_present if l not in LIB_META]
    if unknown:
        sys.exit('columns not in the known library set: %s' % unknown)
    order = [l for l in LIB_ORDER if l in libs_present]
    emit('bare library ids     : %s' % libs_present)
    emit('written in order     : %s' % order)
    emit('(titration order, not lexical -- A10 sorts between A1 and A2)')
    for name, t in tabs.items():
        if t is not None:
            t.rename(columns={c: str(c).rsplit('/')[-1] for c in t.columns}, inplace=True)
            tabs[name] = t[[l for l in order if l in t.columns]]

    tr = tabs['total.ReadCounts']

    # -------------------------------------------------------------- filters
    section('FILTER 1 -- the unnamed row (VASA filter 1: APPLIES)')
    emit('reduceGeneName can produce an empty label, which pandas writes as an')
    emit('unnamed index entry. VASA had exactly 1, carrying 42 UFIs.')
    unnamed = [i for i in tr.index if pd.isna(i) or str(i).strip() == ''
               or str(i).lower() == 'nan']
    emit('unnamed rows found: %d' % len(unnamed))
    for i in unnamed:
        emit('  carrying %d reads across %d libraries'
             % (int(tr.loc[[i]].values.sum()), int((tr.loc[[i]].sum(axis=0) > 0).sum())))

    section('FILTER 2 -- the UMI ceiling (VASA filter 2: DOES NOT APPLY)')
    emit('bc2trans clamps at x >= K, K = 4**len(umi). The only UMI on this path is')
    emit("the literal 'A', so K = 4, while UFI per gene per column is at most 1.")
    emit('The clamp is unreachable. Precondition checked below rather than assumed.')
    tu = tabs['total.UFICounts']
    tt = tabs['total.TranscriptCounts']
    ceiling_ok = None
    if tu is not None:
        mx = int(np.nanmax(tu.values))
        emit('max UFI observed anywhere : %d   (K = 4; clamp needs >= 4)' % mx)
        vals = sorted(pd.unique(tu.values.ravel()))
        emit('distinct UFI values       : %s%s' % (vals[:12], ' ...' if len(vals) > 12 else ''))
        ceiling_ok = mx < 4 and set(vals) <= {0, 1}
        emit('CONFIRMED: no gene can reach the ceiling -> filter has no analogue here'
             if ceiling_ok else
             'UNEXPECTED: UFI is not a 0/1 mask -- re-derive before trusting anything')
        if tt is not None:
            same = np.allclose(tu.values.astype(float), tt.values.astype(float))
            emit('TranscriptCounts == UFICounts elementwise : %s' % same)
            emit('-> neither carries abundance information; ReadCounts only.')

    section('FILTER 3 -- blank barcodes (VASA filter 3: DOES NOT APPLY)')
    emit('There are no barcodes and no blanks: each column is a whole library.')
    emit('qc_verdict plays a different role and is carried as metadata, NOT applied:')
    for l in order:
        v = LIB_META[l]
        emit('  %-11s %-7s %-8s well %-4s %s' % (l, v[0], v[4], v[3], v[5][:52]))
    emit('A8 (exclude) and A7 (caveat) are KEPT in the tables. The verdict filters')
    emit('interpretation, not QC -- dropping A8 would lose the well effect that')
    emit('identified its contamination.')

    section('FILTER 4 -- all-zero rows (NEW here, not a VASA filter)')
    emit("Step 6 builds the row index from the union over libraries, so a gene")
    emit('detected in no library can exist in the frame.')
    allzero = tr.index[(tr.fillna(0).values.sum(axis=1) == 0)]
    allzero = [i for i in allzero if i not in unnamed]
    emit('all-zero rows found: %d of %d' % (len(allzero), len(tr)))

    # ---------------------------------------------------------------- apply
    section('APPLY')
    drop = list(unnamed) + list(allzero)
    keep = [i for i in tr.index if i not in set(drop)]
    emit('rows: %d -> %d  (dropped %d: %d unnamed, %d all-zero)'
         % (len(tr), len(keep), len(drop), len(unnamed), len(allzero)))
    emit('cols: %d -> %d  (no library dropped, by design)' % (len(raw_cols), len(order)))
    dropped_reads = int(tr.loc[drop].values.sum()) if drop else 0
    emit('dropped rows carry %d reads of %d (%.6f%%)'
         % (dropped_reads, int(tr.values.sum()),
            100.0 * dropped_reads / max(1, int(tr.values.sum()))))

    written = {}

    def write(t, fname, what):
        p = os.path.join(tablesdir, fname)
        t.to_csv(p, sep='\t')
        written[fname] = what
        emit('  wrote %-52s %s' % (fname, (t.shape,)))

    pre = 'FS_%s_analysis' % arm
    write(tr.loc[keep], '%s_total.ReadCounts.tsv' % pre,
          'genes x libraries, total reads per gene, arm=%s -- THE quantification' % arm)

    for nm, suf, what in (
        ('uniagg_total.ReadCounts', 'uniagg_total.ReadCounts.tsv',
         'aggregated uni-gene set, total reads'),
        ('uniagg_spliced.ReadCounts', 'uniagg_spliced.ReadCounts.tsv',
         'aggregated uni-gene set, exon-labelled reads'),
        ('uniagg_unspliced.ReadCounts', 'uniagg_unspliced.ReadCounts.tsv',
         'aggregated uni-gene set, intron-labelled reads'),
        ('tRNA.ReadCounts', 'tRNA.ReadCounts.tsv',
         'tRNA reads collapsed by isotype (may be empty at native read length)'),
    ):
        t = tabs.get(nm)
        if t is None:
            emit('  SKIP %s -- not written by step 7' % nm)
            continue
        k = [i for i in t.index if not (pd.isna(i) or str(i).strip() == '')]
        write(t.loc[k], '%s_%s' % (pre, suf), what + ', arm=%s' % arm)

    # Detection mask, derived rather than copied: identical information to
    # UFICounts, but named so it cannot be mistaken for a molecule count.
    mask = (tr.loc[keep] > 0).astype('int8')
    write(mask, '%s_detected.mask.tsv' % pre,
          'boolean detection mask derived from ReadCounts>0 (this is what '
          'UFICounts degenerates to on the no-UMI path), arm=%s' % arm)

    # ------------------------------------------------------------- metadata
    section('GENE METADATA -- for the filtering this script does NOT do')
    kept = tr.loc[keep]
    gm = pd.DataFrame(index=kept.index)
    gm['entry'] = kept.index
    gm['symbol'] = [symbol_of(i) for i in kept.index]
    gm['biotype'] = [biotype_of(i) for i in kept.index]
    gm['n_genes_in_entry'] = [str(i).count('-') + 1 for i in kept.index]
    gm['is_combination'] = gm.n_genes_in_entry > 1
    gm['n_libraries_detected'] = (kept > 0).sum(axis=1)
    gm['total_reads'] = kept.sum(axis=1)
    ut, us = tabs.get('uniagg_total.ReadCounts'), tabs.get('uniagg_unspliced.ReadCounts')
    if ut is not None and us is not None:
        tsum, usum = ut.sum(axis=1), us.sum(axis=1)
        pctu = (100.0 * usum / tsum.replace(0, np.nan)).reindex(kept.index)
        gm['pct_unspliced'] = pctu.round(2)
        # >=95% unspliced: the same threshold the VASA side used to FLAG (not
        # drop) intron-only entries.
        gm['intron_only'] = gm.pct_unspliced >= 95
    p = os.path.join(tablesdir, 'FS_%s_gene_metadata.tsv' % arm)
    gm.to_csv(p, sep='\t', index=False)
    written['FS_%s_gene_metadata.tsv' % arm] = \
        'per-row biotype/symbol/detection, for filters left to the user'
    emit('  wrote FS_%s_gene_metadata.tsv %s' % (arm, (gm.shape,)))
    emit()
    emit('entries: %d   single-gene: %d   combination: %d'
         % (len(gm), int((~gm.is_combination).sum()), int(gm.is_combination.sum())))
    emit()
    emit('single-gene entries by biotype (top 12):')
    sg = gm[~gm.is_combination]
    agg = sg.groupby('biotype').agg(entries=('entry', 'size'), reads=('total_reads', 'sum'))
    for b, r in agg.sort_values('entries', ascending=False).head(12).iterrows():
        emit('  %-34s %8d entries %14d reads' % (b, r.entries, r.reads))

    emit()
    emit('what the OPTIONAL downstream filters would leave (NONE applied here):')
    tot_reads = float(gm.total_reads.sum())
    for label, sel in (
        ('as written (unnamed + all-zero only)', pd.Series(True, index=gm.index)),
        ('+ single-gene entries only', ~gm.is_combination),
        ('+ single-gene AND protein-coding', (~gm.is_combination) & (gm.biotype == 'ProteinCoding')),
        ('+ ... AND detected in >=3 libraries',
         (~gm.is_combination) & (gm.biotype == 'ProteinCoding') & (gm.n_libraries_detected >= 3)),
    ):
        emit('  %-52s %8d entries  %5.1f%% of reads'
             % (label, int(sel.sum()), 100.0 * gm.loc[sel, 'total_reads'].sum() / max(1.0, tot_reads)))

    lm = pd.DataFrame([{
        'library': l, 'arm': arm, 'input_amount': LIB_META[l][0],
        'input_pg': LIB_META[l][1], 'replicate': LIB_META[l][2],
        'well': LIB_META[l][3], 'qc_verdict': LIB_META[l][4],
        'qc_note': LIB_META[l][5],
        'reads_in_tables': int(kept[l].sum()),
        'genes_detected': int((kept[l] > 0).sum()),
    } for l in order])
    p = os.path.join(tablesdir, 'FS_%s_library_metadata.tsv' % arm)
    lm.to_csv(p, sep='\t', index=False)
    written['FS_%s_library_metadata.tsv' % arm] = \
        'per-library input amount, well, qc_verdict, reads, genes detected'
    emit()
    emit('%-11s %-8s %-8s %14s %10s' % ('library', 'input', 'verdict', 'reads', 'genes'))
    for _, r in lm.iterrows():
        emit('%-11s %-8s %-8s %14d %10d'
             % (r.library, r.input_amount, r.qc_verdict, r.reads_in_tables, r.genes_detected))

    # -------------------------------------------------------------- manifest
    mpath = os.path.join(tablesdir, 'manifest.json')
    man = json.load(open(mpath)) if os.path.exists(mpath) else {'arms': {}}
    man.setdefault('arms', {})[arm] = {
        'generated': datetime.now().isoformat(timespec='seconds'),
        'sample': sample,
        'source_dir': outdir,
        'purpose': ('native = adapter-trimmed R1 at its natural length; what the '
                    'FLASH-seq run yields, and the arm the nf-core cross-check '
                    'compares against'
                    if arm == 'native' else
                    "vasalen = additionally hard-trimmed per read to a draw from "
                    "VASA's own STAR-input length distribution; the read-length "
                    'control for the short-biotype axis'),
        'libraries': order,
        'filters_applied': ['unnamed_row', 'all_zero_row'],
        'filters_not_applicable': {
            'umi_ceiling': 'K = 4**len("A") = 4 and UFI <= 1, so the bc2trans clamp '
                           'is unreachable; VASA dropped 8 genes at K = 4096',
            'blank_barcodes': 'no barcodes and no blanks; each column is a library',
        },
        'unnamed_rows_dropped': len(unnamed),
        'all_zero_rows_dropped': len(allzero),
        'shape_before': list(tr.shape),
        'shape_after': [len(keep), len(order)],
        'ufi_is_01_mask': bool(ceiling_ok) if ceiling_ok is not None else None,
        'counting_currency': 'reads (not molecules) -- compare against VASA ReadCounts',
        'files': written,
        'pandas': pd.__version__, 'numpy': np.__version__,
    }
    with open(mpath, 'w') as fh:
        json.dump(man, fh, indent=1)
    emit()
    emit('manifest updated: %s (arms recorded: %s)' % (mpath, sorted(man['arms'])))

    with open(reportpath, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % reportpath, flush=True)


if __name__ == '__main__':
    main()
