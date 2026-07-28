#!/usr/bin/env python3
"""
make_analysis_matrix.py -- turn step 7's raw tables into an analysis-ready set.

Applies exactly THREE filters, each one a decision the user made explicitly
(2026-07-28), and nothing else:

  1. DROP the unnamed row. `reduceGeneName` collapses six rRNA x small-ncRNA
     combinations to the empty string, so step 7 emits one row with no name
     (pandas reads it back as NaN). It is an artefact of the name reduction, not
     a gene.
  2. DROP genes that hit the UMI ceiling. The UMI is 6 nt, so K = 4^6 = 4096
     distinct UMIs exist per gene. `bc2trans` clamps at x >= K to a constant, so
     those genes carry a ceiling value rather than an estimate.
  3. DROP the four blank barcodes (001, 014, 015, 016). They were processed
     identically to real cells on purpose -- they are the negative control --
     and they behave like one at every stage. Confirmed blanks by the user.

WHAT THIS SCRIPT DELIBERATELY DOES *NOT* DO
-------------------------------------------
It does not filter on biotype, and it does not drop multi-gene ("-"-joined)
combination entries, even though after the three filters above the top of the
ranking is still dominated by structural RNA and combination rows (ranks 9-14 of
the raw ranking: Rn7s6-Rn7s2-Rn7s1, Rn18s.rs5, a 12-way snRNA combination, ...).
That is a real analysis decision with a defensible answer either way -- a
total-RNA method is *supposed* to see structural RNA -- so instead of choosing,
this script writes `gene_metadata.tsv` with `biotype`, `is_combination`,
`n_genes_in_entry`, `pct_unspliced` and `intron_only` columns. Filtering is then
one line downstream, e.g.:

    md = pd.read_csv('gene_metadata.tsv', sep='\\t', index_col=0)
    keep = md.index[(md.biotype == 'ProteinCoding') & (~md.is_combination)
                    & (~md.intron_only)]
    mat = mat.loc[keep]

Nor does it drop **intron-only** entries -- 2,104 single-gene protein-coding
entries are >=95% unspliced, against a median of ~34%, and two of them (`Gphn`
99.8%, `Cmss1` 99.5%) rank in the top 15 by transcripts on pure intronic signal.
Same reasoning: unspliced transcription is real signal in a total-RNA method and
is the point of the spliced/unspliced tables, but such an entry is not evidence
the protein is expressed. Flagged, measured in the report, left in.

Nothing is deleted silently: every dropped row is written to
`dropped_rows.tsv` with the reason.

INPUT  : $OUTDIR/<sample>_*.tsv                (step 7 output, untouched)
OUTPUT : $OUTDIR/analysis/                     (new, self-contained)

Usage:
    make_analysis_matrix.py <outdir> <sample> [blank_cells_comma_separated]

    make_analysis_matrix.py .../out ZHA9292A1 001,014,015,016
"""
import sys
import os
import json
import datetime

import numpy as np
import pandas as pd

LEN_UMI = 6                 # config.sh LEN_UMI. Verified against the data below.
K = 4 ** LEN_UMI            # 4096 distinct UMIs per gene


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    outdir = sys.argv[1]
    sample = sys.argv[2]
    blanks = (sys.argv[3].split(',') if len(sys.argv) > 3
              else ['001', '014', '015', '016'])
    blanks = [b.strip() for b in blanks if b.strip()]

    adir = os.path.join(outdir, 'analysis')
    os.makedirs(adir, exist_ok=True)

    out = []
    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 78)
        emit(t)
        emit('=' * 78)

    def T(suffix):
        return os.path.join(outdir, '%s_%s.tsv' % (sample, suffix))

    section('PROVENANCE')
    emit('generated      : %s' % datetime.datetime.now().isoformat(timespec='seconds'))
    emit('script         : %s' % os.path.abspath(__file__))
    emit('sample         : %s' % sample)
    emit('source tables  : %s' % outdir)
    emit('output         : %s' % adir)
    emit('pandas %s / numpy %s / python %s'
         % (pd.__version__, np.__version__, sys.version.split()[0]))
    emit('blank barcodes to exclude: %s' % ', '.join(blanks))

    # ------------------------------------------------------------------
    # The tables that get filtered. All are gene x cell and must stay
    # row-aligned with each other, so the same row mask is applied to all.
    # ------------------------------------------------------------------
    GENE_TABLES = {
        'total.TranscriptCounts':     'uniaggGenes_total.TranscriptCounts',
        'spliced.TranscriptCounts':   'uniaggGenes_spliced.TranscriptCounts',
        'unspliced.TranscriptCounts': 'uniaggGenes_unspliced.TranscriptCounts',
        'total.UFICounts':            'uniaggGenes_total.UFICounts',
        'spliced.UFICounts':          'uniaggGenes_spliced.UFICounts',
        'unspliced.UFICounts':        'uniaggGenes_unspliced.UFICounts',
        'total.ReadCounts':           'uniaggGenes_total.ReadCounts',
    }
    TRNA_TABLES = {'tRNA.UFICounts': 'tRNA.UFICounts',
                   'tRNA.ReadCounts': 'tRNA.ReadCounts'}

    section('LOAD')
    missing = [s for s in list(GENE_TABLES.values()) + list(TRNA_TABLES.values())
               if not os.path.exists(T(s))]
    if missing:
        sys.exit('MISSING input tables: %s' % ', '.join(missing))

    tabs = {}
    for label, suffix in GENE_TABLES.items():
        tabs[label] = pd.read_csv(T(suffix), sep='\t', index_col=0)
        emit('  %-28s %s' % (label, tabs[label].shape))

    # every gene table must have identical row and column labels, or the mask
    # cannot be shared -- check rather than assume
    ref = tabs['total.UFICounts']
    for label, df in tabs.items():
        if not df.index.equals(ref.index):
            sys.exit('ROW MISMATCH: %s does not share an index with total.UFICounts' % label)
        if list(df.columns) != list(ref.columns):
            sys.exit('COL MISMATCH: %s' % label)
    emit('OK: all %d gene tables share an index and column set' % len(tabs))

    # cell ids are zero-padded strings ('001'); read_csv may hand back ints
    cells_raw = list(ref.columns)
    cells = [str(cc).zfill(3) for cc in cells_raw]
    ren = dict(zip(cells_raw, cells))
    for label in tabs:
        tabs[label] = tabs[label].rename(columns=ren)
    emit('cells: %s' % ', '.join(cells))

    # ------------------------------------------------------------------
    section('FILTER 1 -- the unnamed row')
    idx = tabs['total.UFICounts'].index
    unnamed = [i for i in idx if (not isinstance(i, str)) or i.strip() == '']
    if unnamed:
        n = int(tabs['total.UFICounts'].loc[unnamed].values.sum())
        emit('found %d unnamed row(s) carrying %d UFIs -> DROP' % (len(unnamed), n))
    else:
        emit('no unnamed row found (nothing to drop)')

    # ------------------------------------------------------------------
    section('FILTER 2 -- UMI-ceiling genes')
    emit('LEN_UMI = %d  ->  K = %d' % (LEN_UMI, K))
    ufi = tabs['total.UFICounts']
    tr = tabs['total.TranscriptCounts']

    # sanity: the clamp value bc2trans produces at x >= K
    clamp = float(np.log(1. - (float(K) - 1e-3) / K) / np.log(1. - 1. / K))
    emit('bc2trans clamp at x >= K: %.1f transcripts' % clamp)
    obs_max = float(tr.values.max())
    emit('max TranscriptCounts observed: %.1f' % obs_max)
    if abs(obs_max - clamp) > 1.0:
        emit('NOTE: observed max is not the clamp -- either no gene is saturated,')
        emit('      or LEN_UMI does not match this data. Check config.sh.')

    real_cells = [cc for cc in cells if cc not in blanks]
    ceil_any = ufi.index[(ufi >= K).any(axis=1)]
    ceil_real = ufi.index[(ufi[real_cells] >= K).any(axis=1)]
    emit()
    emit('genes at ceiling over ALL %d cells      : %d' % (len(cells), len(ceil_any)))
    emit('genes at ceiling over the %d real cells : %d' % (len(real_cells), len(ceil_real)))
    if set(ceil_any) != set(ceil_real):
        emit('(the two sets differ -- using the ALL-cells set, the stricter one)')
    ceiling = list(ceil_any)
    for g in ceiling:
        ncell = int((ufi.loc[g] >= K).sum())
        emit('  DROP  saturated in %2d/%d cells, max %5d UFIs  %s'
             % (ncell, len(cells), int(ufi.loc[g].max()), g))

    # ------------------------------------------------------------------
    section('FILTER 3 -- blank barcodes')
    present = [b for b in blanks if b in cells]
    absent = [b for b in blanks if b not in cells]
    if absent:
        emit('WARNING: requested blanks not present in the tables: %s' % ', '.join(absent))
    emit('dropping %d blank cells: %s' % (len(present), ', '.join(present)))
    emit('keeping  %d real  cells: %s' % (len(real_cells), ', '.join(real_cells)))

    # ------------------------------------------------------------------
    section('APPLY')
    drop_rows = list(dict.fromkeys(list(unnamed) + ceiling))
    keep_rows = [i for i in idx if i not in set(drop_rows)]
    emit('rows: %d -> %d  (dropped %d)' % (len(idx), len(keep_rows), len(drop_rows)))
    emit('cols: %d -> %d  (dropped %d)' % (len(cells), len(real_cells), len(present)))

    # accountability: what was removed, and how much signal went with it
    rec = []
    for i in drop_rows:
        reason = 'unnamed_row' if i in set(unnamed) else 'umi_ceiling'
        rec.append(dict(entry=('' if not isinstance(i, str) else i),
                        reason=reason,
                        UFIs_all_cells=int(ufi.loc[i].sum()),
                        UFIs_real_cells=int(ufi.loc[i, real_cells].sum()),
                        max_UFIs_in_one_cell=int(ufi.loc[i].max())))
    dropped = pd.DataFrame(rec)
    dropped.to_csv(os.path.join(adir, 'dropped_rows.tsv'), sep='\t', index=False)
    emit()
    emit('dropped rows carry %d UFIs of %d in the real cells (%.4f%%)'
         % (int(ufi.loc[drop_rows, real_cells].values.sum()),
            int(ufi[real_cells].values.sum()),
            100.0 * ufi.loc[drop_rows, real_cells].values.sum()
            / max(1, ufi[real_cells].values.sum())))

    written = []
    for label in GENE_TABLES:
        sub = tabs[label].loc[keep_rows, real_cells]
        p = os.path.join(adir, '%s_analysis_%s.tsv' % (sample, label))
        sub.to_csv(p, sep='\t')
        written.append((os.path.basename(p), sub.shape))

    # tRNA tables: only the cell filter applies (isotype rows, no gene names)
    for label, suffix in TRNA_TABLES.items():
        t = pd.read_csv(T(suffix), sep='\t', index_col=0)
        t = t.rename(columns={cc: str(cc).zfill(3) for cc in t.columns})
        t = t[[cc for cc in real_cells if cc in t.columns]]
        p = os.path.join(adir, '%s_analysis_%s.tsv' % (sample, label))
        t.to_csv(p, sep='\t')
        written.append((os.path.basename(p), t.shape))

    # ------------------------------------------------------------------
    section('GENE METADATA (for the filtering this script does not do)')

    def parse(entry):
        """ENSMUSG..._Xkr4_ProteinCoding -> (symbol, biotype). Combination
        entries are '-'-joined; report the set of biotypes."""
        parts = entry.rsplit('-')
        bios, syms = [], []
        for p_ in parts:
            f = p_.rsplit('_')
            bios.append(f[-1] if len(f) >= 2 else 'unparsed')
            syms.append('_'.join(f[1:-1]) if len(f) >= 3 else (f[0] if f else ''))
        ubio = sorted(set(bios))
        return ('-'.join(syms), '-'.join(ubio) if len(ubio) > 1 else ubio[0], len(parts))

    kept_ufi = tabs['total.UFICounts'].loc[keep_rows, real_cells]
    kept_tr = tabs['total.TranscriptCounts'].loc[keep_rows, real_cells]
    kept_sp = tabs['spliced.UFICounts'].loc[keep_rows, real_cells]
    kept_un = tabs['unspliced.UFICounts'].loc[keep_rows, real_cells]

    # pct_unspliced per entry. This is NOT a filter -- see the block below the
    # table for why it is reported instead.
    s_sum = kept_sp.sum(axis=1).astype(float)
    u_sum = kept_un.sum(axis=1).astype(float)
    denom = (s_sum + u_sum).replace(0, np.nan)
    pct_unsp = (u_sum / denom * 100.0)

    meta = []
    for e in keep_rows:
        sym, bio, n = parse(e)
        pu = pct_unsp.loc[e]
        meta.append(dict(entry=e, symbol=sym, biotype=bio, n_genes_in_entry=n,
                         is_combination=(n > 1),
                         n_cells_detected=int((kept_ufi.loc[e] > 0).sum()),
                         total_UFIs=int(kept_ufi.loc[e].sum()),
                         total_transcripts=round(float(kept_tr.loc[e].sum()), 1),
                         pct_unspliced=(round(float(pu), 2) if pd.notna(pu) else np.nan),
                         intron_only=(bool(pd.notna(pu) and pu >= 95.0))))
    md = pd.DataFrame(meta).set_index('entry')
    md.to_csv(os.path.join(adir, 'gene_metadata.tsv'), sep='\t')
    written.append(('gene_metadata.tsv', md.shape))

    emit('entries: %d   single-gene: %d   combination: %d'
         % (len(md), int((~md.is_combination).sum()), int(md.is_combination.sum())))
    emit()
    emit('single-gene entries by biotype (top 12):')
    sg = md[~md.is_combination]
    for b, n in sg.biotype.value_counts().head(12).items():
        emit('  %-34s %7d entries  %12d UFIs'
             % (b, n, int(sg.loc[sg.biotype == b, 'total_UFIs'].sum())))

    emit()
    emit('what the OPTIONAL downstream filters would leave (none applied here):')
    _sg = ~md.is_combination
    _pc = md.biotype == 'ProteinCoding'
    _io = md.intron_only.fillna(False).astype(bool)
    for desc, sel in [
        ('as written (3 agreed filters only)', md.index),
        ('+ single-gene entries only', md.index[_sg]),
        ('+ single-gene AND protein-coding', md.index[_sg & _pc]),
        ('+ ... AND not intron-only (>=95% unspliced)', md.index[_sg & _pc & ~_io]),
        ('+ ... AND detected in >=3 cells',
         md.index[_sg & _pc & ~_io & (md.n_cells_detected >= 3)]),
    ]:
        emit('  %-52s %7d entries  %5.1f%% of UFIs'
             % (desc, len(sel), 100.0 * kept_ufi.loc[sel].values.sum()
                / max(1, kept_ufi.values.sum())))

    # ------------------------------------------------------------------
    section('INTRON-ONLY ENTRIES -- reported, deliberately NOT filtered')
    pc = md[(~md.is_combination) & (md.biotype == 'ProteinCoding')]
    pcv = pc.dropna(subset=['pct_unspliced'])
    emit('Distribution of pct_unspliced over %d single-gene protein-coding entries'
         % len(pcv))
    emit('with any signal:')
    for q in (5, 25, 50, 75, 90, 95, 99):
        emit('    p%-3d %6.1f%%' % (q, float(np.percentile(pcv.pct_unspliced, q))))
    emit()
    for thr in (90, 95, 99):
        m = pcv.pct_unspliced >= thr
        emit('  >= %d%% unspliced : %5d entries, %9d UFIs (%.2f%% of protein-coding UFIs)'
             % (thr, int(m.sum()), int(pcv.total_UFIs[m].sum()),
                100.0 * pcv.total_UFIs[m].sum() / max(1, pcv.total_UFIs.sum())))
    emit()
    emit('The MEDIAN entry is ~34% unspliced, matching the library-wide ~30%; the')
    emit('>=95% class is therefore a distinct population, not the tail of one')
    emit('distribution. Entries here have essentially NO exonic signal: their')
    emit('counts come from intronic reads in a long locus.')
    emit()
    emit('This is flagged (`intron_only`, `pct_unspliced` in gene_metadata.tsv) and')
    emit('NOT filtered, because the honest answer depends on the question:')
    emit('  - VASA captures total RNA, so nascent/unspliced transcription is real')
    emit('    signal and is the whole point of the spliced/unspliced tables;')
    emit('  - but an entry at 99.8% unspliced is not evidence the PROTEIN is')
    emit('    expressed, and it will behave like a highly-expressed gene in any')
    emit('    analysis that reads `total` as expression.')
    emit('Decide per analysis. For steady-state mRNA abundance, exclude them or use')
    emit('the spliced table; for velocity/nascent transcription, keep them.')
    emit()
    hi = pcv[(pcv.pct_unspliced >= 95) & (pcv.total_UFIs >= 1000)] \
        .sort_values('total_UFIs', ascending=False)
    emit('Intron-only AND abundant (>=1000 UFIs): %d entries. Top 20:' % len(hi))
    emit('  %-16s %10s %9s' % ('symbol', 'UFIs', 'pct_unsp'))
    for e, row in hi.head(20).iterrows():
        emit('  %-16s %10d %8.1f%%' % (str(row.symbol), int(row.total_UFIs),
                                       float(row.pct_unspliced)))

    t200 = pc.sort_values('total_transcripts', ascending=False).head(200)
    n200 = int((t200.pct_unspliced >= 95).sum())
    emit()
    emit('Of the top 200 protein-coding entries by transcripts, %d are intron-only:'
         % n200)
    emit('  %s' % ', '.join(str(x) for x in t200.symbol[t200.pct_unspliced >= 95]))

    # ------------------------------------------------------------------
    section('TOP 20 SINGLE-GENE PROTEIN-CODING ENTRIES (transcripts, %d cells)'
            % len(real_cells))
    emit('  %10s %8s  %-14s %s' % ('transcr', 'pct_uns', 'symbol', 'entry'))
    for e, row in pc.sort_values('total_transcripts', ascending=False).head(20).iterrows():
        flag = '  <- intron-only' if (pd.notna(row.pct_unspliced) and row.pct_unspliced >= 95) else ''
        emit('  %10.0f %7.1f%%  %-14s %s%s'
             % (row.total_transcripts,
                float(row.pct_unspliced) if pd.notna(row.pct_unspliced) else -1,
                row.symbol, e, flag))

    # ------------------------------------------------------------------
    section('PER-CELL, AFTER FILTERING')
    emit('%-6s %10s %14s %14s %12s %8s' %
         ('cell', 'genes', 'UFIs', 'transcripts', 'unspliced', '%unspl'))
    rows = []
    for cc in real_cells:
        u = kept_ufi[cc]
        uns = tabs['unspliced.UFICounts'].loc[keep_rows, cc]
        pct = 100.0 * uns.sum() / u.sum() if u.sum() else 0.0
        emit('%-6s %10d %14d %14.0f %12d %8.2f'
             % (cc, int((u > 0).sum()), int(u.sum()), float(kept_tr[cc].sum()),
                int(uns.sum()), pct))
        rows.append(dict(cell=cc, genes_detected=int((u > 0).sum()),
                         UFIs=int(u.sum()), transcripts=round(float(kept_tr[cc].sum()), 1),
                         unspliced_UFI=int(uns.sum()), pct_unspliced=round(pct, 2)))
    pc_df = pd.DataFrame(rows).set_index('cell')
    pc_df.to_csv(os.path.join(adir, 'per_cell_after_filter.tsv'), sep='\t')
    written.append(('per_cell_after_filter.tsv', pc_df.shape))

    # ------------------------------------------------------------------
    section('WRITTEN')
    for n, shp in written:
        emit('  %-52s %s' % (n, shp))

    manifest = dict(
        generated=datetime.datetime.now().isoformat(timespec='seconds'),
        sample=sample, source_dir=os.path.abspath(outdir),
        filters=dict(
            unnamed_rows_dropped=len(unnamed),
            umi_ceiling_genes_dropped=[str(g) for g in ceiling],
            K=K, len_umi=LEN_UMI,
            blank_cells_dropped=present),
        shape_before=[len(idx), len(cells)],
        shape_after=[len(keep_rows), len(real_cells)],
        pandas=pd.__version__, numpy=np.__version__)
    with open(os.path.join(adir, 'manifest.json'), 'w') as fh:
        json.dump(manifest, fh, indent=1)

    rp = os.path.join(adir, 'filter_report.txt')
    with open(rp, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % rp, flush=True)


if __name__ == '__main__':
    main()
