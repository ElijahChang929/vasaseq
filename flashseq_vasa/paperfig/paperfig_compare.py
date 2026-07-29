#!/usr/bin/env python3
"""paperfig_compare.py -- FLASH-seq vs VASA in the published figures' own terms.

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
Three panels of the VASA-seq paper are reproduced here with our two datasets in
place of the paper's method panel:

  Fig 2b  fraction of transcripts per biotype class
          ** THE AUTHORS' OWN CODE **. The grouping is lifted character-for-
          character from b_Analysis/02_scanpy_QCxBiotype.py (see PAPER_* below),
          including their sncRNA membership list and their TF/cofactor lists.
          Nothing here is my judgement.

  Fig 1e  gene-body coverage, 5' -> 3'
  Fig 1f  genes per cell vs reads per cell after trimming
          ** MY REIMPLEMENTATION **. Grepped every .py/.R in the repo: the code
          for these two panels was never deposited. The AXES, BINS and UNITS are
          taken from the published panels (1e: 0-1 gene body, 'Coverage' summing
          to 1 over 100 bins; 1f: 5k/10k/15k/20k/25k/50k/75k reads per cell), so
          the plots are directly comparable -- but the code is mine. Any
          difference from the paper could be my reimplementation rather than a
          real protocol difference, and that cannot be removed without the
          authors' code.

This file is a FORK. It does not import from, write to, or modify anything in
b_Analysis/ or a_Mapping/. The originals are untouched.

WHY NOT JUST RUN 02_scanpy_QCxBiotype.py
----------------------------------------
It expects the paper's merged feather tables, a four-timepoint sample sheet, and
scanpy objects keyed by their cell naming. We have per-library TSVs from
countTables_fromPickle.py. So the DEFINITIONS are copied verbatim and applied to
our matrices, rather than the driver being run on inputs it was not written for.
The functions below are marked PAPER (copied) or MINE (reimplemented) so the
provenance of every number is visible.

Usage:
  paperfig_compare.py fig2b <vasa_tsv> <fs_native_tsv> <fs_vasalen_tsv> \\
                            <tf_txt> <cof_txt> <out_dir>
  paperfig_compare.py fig1f <vasa_tsv> <fs_native_tsv> <fs_vasalen_tsv> <out_dir>
"""
import os
import sys

import numpy as np
import pandas as pd

BLANKS = {'001', '014', '015', '016'}

# ---------------------------------------------------------------------------
# PAPER -- copied verbatim from b_Analysis/02_scanpy_QCxBiotype.py.
#
# add_metadata() builds `ubiotype` by splitting a combination entry on '-',
# taking the biotype suffix of each part, de-duplicating and sorting. Then
# biotype_split() selects on membership of these exact lists. Both are
# reproduced here unchanged -- in particular the sncRNA/smallRNA membership,
# which is the authors' choice and not mine:
#
#   ('ProteinCoding', ["ProteinCoding"])
#   ('lncRNA',        ["lncRNA"])
#   ('smallRNA',      ["snRNA","snoRNA","MiscRna","scaRNA",'ribozyme','miRNA'])
#   ('tRNA',          ['tRNA'])
#
# and the TF / cofactor assignment, which keys on the gene SYMBOL against
# Mus_musculus_TF.txt / Mus_musculus_TF_cofactors.txt.
# ---------------------------------------------------------------------------
PAPER_SMALLRNA = ["snRNA", "snoRNA", "MiscRna", "scaRNA", 'ribozyme', 'miRNA']
PAPER_CLASSES = [('ProteinCoding', ["ProteinCoding"]),
                 ('lncRNA', ["lncRNA"]),
                 ('smallRNA', PAPER_SMALLRNA),
                 ('tRNA', ['tRNA'])]


def paper_ubiotype(entry):
    """PAPER. adata.var['biotype'] then ['ubiotype'], for a non-s/u index label.

    Original:
        biotype  = '-'.join([k.rsplit('_')[-1] for k in idx.rsplit('-')])
        ubiotype = '-'.join(sorted(set(b.rsplit('-'))))
    """
    biotype = '-'.join([k.rsplit('_')[-1] for k in entry.rsplit('-')])
    return '-'.join(sorted(set(biotype.rsplit('-'))))


def paper_id(entry):
    """PAPER. adata.var['id'] -- the gene id(s) of an entry.

        id = '-'.join([k.rsplit('_')[0] for k in idx.rsplit('-')])
    """
    return '-'.join([k.rsplit('_')[0] for k in entry.rsplit('-')])


def paper_symbol(entry):
    """MINE, but forced by a format difference, and it is the one place where the
    paper's code cannot be applied unchanged.

    adata.var['reg'] tests `id in list(tfs['Symbol'])` -- i.e. the paper's `id`
    field held gene SYMBOLS. Our entries are ENSMUSG..._Symbol_Biotype, so
    paper_id() returns the ENSMUSG accession and the symbol is the middle field.
    Extracting it here keeps the paper's TF test (`symbol in TF list`) intact
    rather than silently matching nothing.
    """
    parts = entry.rsplit('_')
    return parts[1] if len(parts) >= 3 else parts[0]


# ---------------------------------------------------------------------------
def load(path, drop_blanks=False):
    df = pd.read_csv(path, sep='\t', index_col=0)
    df.columns = [str(c).split('/')[-1] for c in df.columns]
    if drop_blanks:
        df.columns = [c.zfill(3) for c in df.columns]
        df = df[[c for c in df.columns if c not in BLANKS]]
    df = df[~df.index.isna()]
    return df


def fig2b(vasa_tsv, fsn_tsv, fsv_tsv, tf_txt, cof_txt, outdir):
    tfs = set(pd.read_csv(tf_txt, sep='\t')['Symbol'].astype(str))
    cofs = set(pd.read_csv(cof_txt, sep='\t')['Symbol'].astype(str))
    print('PAPER TF list: %d symbols; cofactor list: %d symbols' % (len(tfs), len(cofs)))

    mats = {'VASA': load(vasa_tsv, True),
            'FLASH-seq native': load(fsn_tsv),
            'FLASH-seq vasalen': load(fsv_tsv)}

    rows = []
    for name, m in mats.items():
        ent = m.index.astype(str)
        ub = np.array([paper_ubiotype(e) for e in ent])
        sym = np.array([paper_symbol(e) for e in ent])
        reg = np.array(['TF' if s in tfs else ('Cof' if s in cofs else '-') for s in sym])
        # PAPER: n_counts is the FULL matrix sum; each class fraction is
        # n_counts_<class> / n_counts. rRNA is NOT excluded by the paper's code,
        # so it is left in -- otherwise the denominator would not be theirs.
        tot = m.sum(axis=0)
        print('%-19s %6d entries, %2d units, %s counts' % (
            name, m.shape[0], m.shape[1], format(int(tot.sum()), ',')))
        for cls, bts in PAPER_CLASSES:
            sel = np.isin(ub, bts)
            frac = m[sel].sum(axis=0) / tot
            for u, f in frac.items():
                rows.append(dict(dataset=name, unit=u, panel=cls, fraction=float(f),
                                 n_entries=int(sel.sum()), source='PAPER code'))
        for cls, tag in [('TF', 'TF'), ('Cofactor', 'Cof')]:
            sel = reg == tag
            frac = m[sel].sum(axis=0) / tot
            for u, f in frac.items():
                rows.append(dict(dataset=name, unit=u, panel=cls, fraction=float(f),
                                 n_entries=int(sel.sum()), source='PAPER code'))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(outdir, 'paperfig2b_fractions.tsv'), sep='\t', index=False)

    print('\n=== Fig 2b, the paper\'s own classes and denominator ===')
    print('  %-14s %-19s %8s %8s %8s  %s' % ('panel', 'dataset', 'min', 'median', 'max', 'n_entries'))
    for cls in ['ProteinCoding', 'lncRNA', 'smallRNA', 'tRNA', 'TF', 'Cofactor']:
        for name in mats:
            d = out[(out.panel == cls) & (out.dataset == name)]
            if not len(d):
                continue
            print('  %-14s %-19s %8.4f %8.4f %8.4f  %d' % (
                cls, name, d.fraction.min(), d.fraction.median(), d.fraction.max(),
                d.n_entries.iloc[0]))
    print('\nwrote paperfig2b_fractions.tsv')


def fig1f(vasa_tsv, fsn_tsv, fsv_tsv, outdir):
    """MINE (reimplementation). Axes and depth grid taken from the published panel:
    x = reads per cell after trimming at 5k, 10k, 15k, 20k, 25k, 50k, 75k;
    y = average number of genes per cell.

    The paper does not say how it downsampled. This uses the binomial-thinning
    expectation, E[genes] = sum_g [1 - (1-p)^{c_g}], which is exact in expectation
    and deterministic -- the same estimator used elsewhere in this project. A
    resampling implementation would give the same curve to within its own noise.
    """
    GRID = [5000, 10000, 15000, 20000, 25000, 50000, 75000]
    mats = {'VASA': load(vasa_tsv, True),
            'FLASH-seq native': load(fsn_tsv),
            'FLASH-seq vasalen': load(fsv_tsv)}
    rows = []
    for name, m in mats.items():
        for u in m.columns:
            c = m[u].values.astype(float)
            c = c[c > 0]
            tot = c.sum()
            for g in GRID:
                p = min(1.0, g / tot)
                rows.append(dict(dataset=name, unit=u, depth=g,
                                 genes=float(np.sum(1.0 - np.power(1.0 - p, c))),
                                 total_reads=int(tot), source='MINE (reimplementation)'))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(outdir, 'paperfig1f_saturation.tsv'), sep='\t', index=False)
    print('\n=== Fig 1f, paper axes, my reimplementation ===')
    piv = out.pivot_table(index='depth', columns='dataset', values='genes', aggfunc='mean')
    print(piv.round(0).to_string())
    print('\nwrote paperfig1f_saturation.tsv')


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'fig2b':
        fig2b(*sys.argv[2:8])
    elif mode == 'fig1f':
        fig1f(*sys.argv[2:6])
    else:
        sys.exit(__doc__)
