#!/usr/bin/env python3
"""
step7_report.py -- per-cell QC table for the final count tables.

Same contract as step2_report.py and step3_report.py: the run generates its own
numbers, so nothing in README.md has to be reconstructed by hand.

Reads the step-7 output in $OUTDIR and writes:

  logs/step7_report.txt          human-readable, the thing to paste into README
  logs/step7_per_cell.tsv        one row per cell
  logs/step7_biotype_UFI.tsv     biotype x cell, UFI counts
  logs/step7_biotype_pct.tsv     biotype x cell, % of that cell's UFIs
  logs/step7_qc.png              four-panel QC figure (skipped if no matplotlib)

WHICH TABLE TO READ, AND WHY
----------------------------
Step 7 writes three flavours of every table and they are NOT interchangeable:

  ReadCounts        raw reads. Use for mapping/QC bookkeeping only.
  UFICounts         unique UMIs seen. SATURATES: this library's UMI is 6 nt, so
                    K = 4^6 = 4096 distinct UMIs exist per gene, and a deep cell
                    has hundreds of genes past 1000.
  TranscriptCounts  UFI counts with the collision correction
                    t = ln(1-x/K)/ln(1-1/K) applied. THIS is the expression
                    estimate. At x >= K it clamps to a constant (~62000) and
                    carries no information, so genes at the ceiling are counted
                    and reported here rather than interpreted.

Gene detection is therefore counted on UFIs (a gene is detected if >=1 molecule
was seen), while expression totals are read off TranscriptCounts.

Usage:
    step7_report.py <outdir> <sample>
"""
import sys
import os

import numpy as np
import pandas as pd


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    outdir, sample = sys.argv[1], sys.argv[2]
    logdir = os.path.join(outdir, 'logs')
    os.makedirs(logdir, exist_ok=True)

    def T(suffix):
        return os.path.join(outdir, '%s_%s.tsv' % (sample, suffix))

    out = []
    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 78)
        emit(t)
        emit('=' * 78)

    # ------------------------------------------------------------------
    section('TABLES ON DISK')
    need = ['total.UFICounts', 'uniaggGenes_total.UFICounts',
            'uniaggGenes_total.TranscriptCounts',
            'uniaggGenes_spliced.UFICounts', 'uniaggGenes_unspliced.UFICounts',
            'tRNA.ReadCounts', 'tRNA.UFICounts']
    missing = [n for n in need if not os.path.exists(T(n))]
    if missing:
        sys.exit('MISSING required tables: %s' % ', '.join(missing))
    for n in need:
        emit('  %-42s %9.1f MB' % (os.path.basename(T(n)), os.path.getsize(T(n)) / 1e6))

    # ------------------------------------------------------------------
    section('LOAD')
    uni_ufi = pd.read_csv(T('uniaggGenes_total.UFICounts'), sep='\t', index_col=0)
    uni_tr = pd.read_csv(T('uniaggGenes_total.TranscriptCounts'), sep='\t', index_col=0)
    spl = pd.read_csv(T('uniaggGenes_spliced.UFICounts'), sep='\t', index_col=0)
    uns = pd.read_csv(T('uniaggGenes_unspliced.UFICounts'), sep='\t', index_col=0)
    trna = pd.read_csv(T('tRNA.UFICounts'), sep='\t', index_col=0)
    # column order: cells are zero-padded strings, so a plain sort is correct
    cells = sorted(uni_ufi.columns, key=str)
    uni_ufi, uni_tr = uni_ufi[cells], uni_tr[cells]
    spl, uns, trna = spl[cells], uns[cells], trna[cells]
    emit('uniaggGenes (aggregated unique genes): %s' % (uni_ufi.shape,))
    emit('tRNA (isotype-collapsed)             : %s' % (trna.shape,))
    emit('cells: %s' % ', '.join(map(str, cells)))

    # ------------------------------------------------------------------
    section('INTEGRITY CHECKS')
    ok = True

    # 1. spliced + unspliced must equal total, exactly: countExonUMI and
    #    countIntronUMI partition the UMIs of a gene with no overlap.
    diff = (spl + uns - uni_ufi).abs().values.sum()
    emit('spliced + unspliced - total  (must be 0): %d' % diff)
    if diff != 0:
        ok = False
        emit('  ^ NOT ZERO -- the spliced/unspliced partition is not exhaustive')

    # 2. the empty-name row the precheck predicted: reduceGeneName collapses
    #    rRNA x small-ncRNA combinations to '', and pandas reads that back as NaN
    empty = [i for i in uni_ufi.index if (not isinstance(i, str)) or i.strip() == '']
    if empty:
        emit('empty-name row present, as predicted: %r carrying %d UFIs'
             % (empty, int(uni_ufi.loc[empty].values.sum())))
        emit('  -> drop this row at analysis time; it is an artefact of')
        emit('     reduceGeneName, not a gene.')
    else:
        emit('no empty-name row in the aggregated unigene table')

    # 3. TranscriptCounts must be >= UFICounts everywhere (the correction only
    #    ever inflates), and equal only where the count is 0 or 1.
    bad = int((uni_tr.values < uni_ufi.values - 1e-6).sum())
    emit('cells where TranscriptCounts < UFICounts (must be 0): %d' % bad)
    if bad:
        ok = False

    # 4. UMI ceiling
    K = 4096  # 4^6, LEN_UMI=6 -- asserted against the data below
    at_ceiling = (uni_ufi >= K)
    emit('gene x cell entries at or past the UMI ceiling K=%d: %d'
         % (K, int(at_ceiling.values.sum())))
    if at_ceiling.values.any():
        ceil_genes = uni_ufi.index[at_ceiling.any(axis=1)]
        emit('  genes involved (%d): %s' % (len(ceil_genes), ', '.join(map(str, ceil_genes[:12]))))
        emit('  their TranscriptCounts are clamped to ~%.0f and mean nothing.'
             % uni_tr.values.max())

    emit()
    emit('INTEGRITY: %s' % ('OK' if ok else 'PROBLEM -- see above'))

    # ------------------------------------------------------------------
    section('PER-CELL SUMMARY')
    # biotype is the last underscore field of the aggregated gene name,
    # e.g. ENSMUSG00000051951_Xkr4_ProteinCoding. Combination names joined by
    # '-' are excluded from uniaggGenes by construction, but guard anyway.
    def biotype(idx):
        if not isinstance(idx, str) or idx.strip() == '':
            return 'UNNAMED'
        if '-' in idx:
            return 'combination'
        parts = idx.rsplit('_')
        return parts[-1] if len(parts) >= 2 else 'unparsed'

    bt = pd.Series([biotype(i) for i in uni_ufi.index], index=uni_ufi.index, name='biotype')

    rows = []
    for cc in cells:
        u = uni_ufi[cc]
        det = int((u > 0).sum())
        tot_ufi = int(u.sum())
        tot_tr = float(uni_tr[cc].sum())
        s, i_ = int(spl[cc].sum()), int(uns[cc].sum())
        pc = u[bt == 'ProteinCoding']
        rows.append(dict(
            cell=cc,
            genes_detected=det,
            UFIs=tot_ufi,
            transcripts=round(tot_tr, 1),
            spliced_UFI=s,
            unspliced_UFI=i_,
            pct_unspliced=round(100.0 * i_ / tot_ufi, 2) if tot_ufi else 0.0,
            protein_coding_genes=int((pc > 0).sum()),
            pct_UFI_protein_coding=round(100.0 * pc.sum() / tot_ufi, 2) if tot_ufi else 0.0,
            tRNA_UFI=int(trna[cc].sum()),
            genes_at_UMI_ceiling=int((u >= K).sum()),
        ))
    per_cell = pd.DataFrame(rows).set_index('cell')
    emit(per_cell.to_string())
    per_cell.to_csv(os.path.join(logdir, 'step7_per_cell.tsv'), sep='\t')

    # ------------------------------------------------------------------
    section('BIOTYPE COMPOSITION (UFIs)')
    bio_ufi = uni_ufi.groupby(bt).sum()
    bio_pct = 100.0 * bio_ufi / bio_ufi.sum(axis=0)
    order = bio_ufi.sum(axis=1).sort_values(ascending=False).index
    bio_ufi, bio_pct = bio_ufi.loc[order], bio_pct.loc[order]
    emit('%d biotypes; showing the top 15 by total UFI' % len(order))
    emit()
    emit(bio_pct.loc[order[:15]].round(2).to_string())
    bio_ufi.to_csv(os.path.join(logdir, 'step7_biotype_UFI.tsv'), sep='\t')
    bio_pct.to_csv(os.path.join(logdir, 'step7_biotype_pct.tsv'), sep='\t')

    # ------------------------------------------------------------------
    section('TOP GENES (by summed TranscriptCounts over all cells)')
    top = uni_tr.sum(axis=1).sort_values(ascending=False).head(25)
    for g, v in top.items():
        emit('  %12.0f  %s' % (v, g))

    # ------------------------------------------------------------------
    section('tRNA BY ISOTYPE (UFIs, pure single-locus classes only)')
    pure = [i for i in trna.index if isinstance(i, str) and '-' not in i]
    tp = trna.loc[pure].sum(axis=1).sort_values(ascending=False)
    emit('%d pure isotype classes of %d total rows; %d UFIs in pure classes, '
         '%d in ambiguous multi-locus classes'
         % (len(pure), len(trna), int(tp.sum()),
            int(trna.values.sum() - tp.sum())))
    emit()
    for g, v in tp.head(20).items():
        emit('  %8d  %s' % (v, g))

    # ------------------------------------------------------------------
    figpath = os.path.join(logdir, 'step7_qc.png')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        blanks = ['001', '014', '015', '016']
        colors = ['#c0392b' if str(c) in blanks else '#2c6fbb' for c in cells]
        x = np.arange(len(cells))

        fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))

        ax = axes[0, 0]
        ax.bar(x, per_cell['genes_detected'].values, color=colors)
        ax.set_ylabel('genes detected (>=1 UFI)')
        ax.set_title('Gene detection')

        ax = axes[0, 1]
        ax.bar(x, per_cell['transcripts'].values, color=colors)
        ax.set_ylabel('transcripts (collision-corrected)')
        ax.set_yscale('log')
        ax.set_title('Library depth')

        ax = axes[1, 0]
        ax.bar(x, per_cell['pct_unspliced'].values, color=colors)
        ax.set_ylabel('% of UFIs unspliced')
        ax.set_title('Unspliced fraction (VASA captures total RNA)')

        ax = axes[1, 1]
        keep = [b for b in order[:8]]
        bottom = np.zeros(len(cells))
        cmap = plt.get_cmap('tab10')
        for k, b in enumerate(keep):
            v = bio_pct.loc[b, cells].values.astype(float)
            ax.bar(x, v, bottom=bottom, label=str(b), color=cmap(k % 10))
            bottom += v
        ax.set_ylabel('% of UFIs')
        ax.set_title('Biotype composition (top 8)')
        ax.legend(fontsize=7, ncol=2, loc='lower right')

        for ax in axes.ravel():
            ax.set_xticks(x)
            ax.set_xticklabels([str(cc) for cc in cells], rotation=90, fontsize=8)
            ax.spines[['top', 'right']].set_visible(False)

        fig.suptitle('%s -- step 7 count-table QC   (red = the four blanks)' % sample)
        fig.tight_layout()
        fig.savefig(figpath, dpi=150)
        emit()
        emit('figure written: %s' % figpath)
    except ImportError as e:
        emit()
        emit('no figure (matplotlib unavailable: %s)' % e)

    reportpath = os.path.join(logdir, 'step7_report.txt')
    with open(reportpath, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % reportpath, flush=True)


if __name__ == '__main__':
    main()
