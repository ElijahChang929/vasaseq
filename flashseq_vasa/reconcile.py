#!/usr/bin/env python3
"""
reconcile.py -- per-library read reconciliation for one FLASH-seq arm, plus the
three checks that decide what the count tables MEAN.

The chain, per library:

  delivered FASTQ records
    -> after adapter trim (cutadapt)
      -> after the vasalen hard trim, if this is that arm
        -> STAR "Number of input reads"
          -> STAR uniquely mapped
            -> step-5 singlemapper BED rows
              -> reads present in the step-7 _total.ReadCounts column

A break anywhere is a finding, not something to smooth over, so every link is
printed with its ratio and the ones that CANNOT hold as equalities are named
explicitly rather than left for the reader to wonder about:

  * BED rows > uniquely mapped reads is EXPECTED and is not a double count. One
    read overlapping k annotated features produces k rows; step 6 groups the BED
    by read name (df.groupby('Name')) and assigns each read once. The ratio is
    reported as features-per-read.
  * table reads < STAR input is EXPECTED. It drops unmapped reads,
    multimapper-only reads that lost their tie-break, and reads overlapping no
    annotated feature.
  * table reads > STAR input would be a REAL double count and is the one thing
    that fails this script.

The three semantic checks (findings 3f, 3h, and the tRNA/short-species axis):

  1. UFICounts must be a 0/1 detection mask, and TranscriptCounts must equal it
     elementwise. With one literal UMI key 'A', countTotalUMI returns len({'A':..})
     == 1 for any detected gene, and bc2trans at K = 4**len('A') = 4 maps 1 -> 1.0
     exactly. Verified numerically, not trusted: this is what justifies carrying
     ReadCounts and only ReadCounts forward.
  2. spliced + unspliced vs total, per library. The no-UMI exon/intron branches
     test SUBSTRING containment ('exon' in k / 'intron' in k) where the vasa
     branch tests exact membership, so a combination label 'exon-intron' is
     counted in BOTH and spliced+unspliced can exceed total. It was 0 in the A9
     dry run; that is not a guarantee, so it is measured here and reported per
     library.
  3. The short-biotype tables, per library: tRNA, miRNA, snoRNA, snRNA and the
     aligned-span distribution. This is what the native-vs-vasalen comparison
     turns on, so the numbers are produced here for BOTH arms in the same format.

Usage: reconcile.py <outdir> <sample> <celldir> <logdir> <arm> <out.tsv> <report.txt>
"""
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd

SHORT_BIOTYPES = ['tRNA', 'miRNA', 'snoRNA', 'snRNA', 'MiscRna', 'scaRNA',
                  'ribozyme', 'MtTrna', 'MtRrna', 'rRNA']


def kv(path):
    """Read a two-column key<TAB>value file written by the pipeline stages."""
    d = {}
    if os.path.exists(path):
        for line in open(path):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                d[p[0]] = p[1]
    return d


def star_log(path):
    d = {}
    if not os.path.exists(path):
        return d
    for line in open(path):
        if 'Number of input reads' in line:
            d['star_input'] = int(line.strip().rsplit('\t')[-1])
        elif 'Uniquely mapped reads number' in line:
            d['star_uniq'] = int(line.strip().rsplit('\t')[-1])
        elif 'Number of reads mapped to multiple loci' in line:
            d['star_multi'] = int(line.strip().rsplit('\t')[-1])
        elif 'too short' in line and '%' in line:
            d['star_tooshort_pct'] = line.strip().rsplit('\t')[-1]
        elif 'Average input read length' in line:
            d['star_avg_input_len'] = int(line.strip().rsplit('\t')[-1])
        elif 'Average mapped length' in line:
            d['star_avg_mapped_len'] = float(line.strip().rsplit('\t')[-1])
    return d


def biotype_of(idx):
    """Biotype of a step-7 row label.

    After step 6 strips the exon/intron field, a single-gene label is
    ENSID_Symbol_Biotype; tRNA rows are the 'N.tRNAn.IsotypeNNN_tRNA' form; and a
    combination row is those joined by '-'. A row whose components disagree is
    reported as MULTI:<sorted set> rather than being assigned to one of them.
    """
    idx = str(idx)
    parts = idx.rsplit('-')
    bs = sorted(set(p.rsplit('_')[-1] for p in parts))
    return bs[0] if len(bs) == 1 else 'MULTI:' + '-'.join(bs)


def main():
    if len(sys.argv) < 8:
        sys.exit(__doc__)
    outdir, sample, celldir, logdir, arm, tsvpath, reportpath = sys.argv[1:8]

    out = []

    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 78)
        emit(t)
        emit('=' * 78)

    fail, warn = [], []

    def rd(suffix):
        p = os.path.join(outdir, '%s%s' % (sample, suffix))
        return pd.read_csv(p, sep='\t', index_col=0) if os.path.exists(p) else None

    tot_r = rd('_total.ReadCounts.tsv')
    if tot_r is None:
        sys.exit('no %s_total.ReadCounts.tsv in %s -- step 7 did not complete'
                 % (sample, outdir))
    tot_u = rd('_total.UFICounts.tsv')
    tot_t = rd('_total.TranscriptCounts.tsv')

    # Step 6 with CELLID_FROM=f keeps the folder string it globbed, so the column
    # comes out as 'cells/<LIB>'. Map back to the bare library id.
    col2lib = {c: str(c).rsplit('/')[-1] for c in tot_r.columns}
    libs = [col2lib[c] for c in tot_r.columns]

    section('ARM AND SHAPE')
    emit('arm            : %s' % arm)
    emit('sample         : %s' % sample)
    emit('outdir         : %s' % outdir)
    emit('_total.ReadCounts shape : %s' % (tot_r.shape,))
    emit('columns as written by step 6 : %s' % list(tot_r.columns))
    emit('libraries                    : %s' % libs)

    # ------------------------------------------------------------------ chain
    section('1. RECONCILIATION CHAIN, PER LIBRARY')
    emit('FASTQ -> trimmed -> STAR input -> uniquely mapped -> BED rows -> assigned')
    emit()
    rows = []
    for col in tot_r.columns:
        lib = col2lib[col]
        p = kv(os.path.join(logdir, 'prep_%s.tsv' % lib))
        a = kv(os.path.join(logdir, 'assign_%s.tsv' % lib))
        v = kv(os.path.join(logdir, 'vasalen_%s.tsv' % lib))
        s = star_log(os.path.join(celldir, '%s_cbc_noumi_E99_Log.final.txt' % lib))
        r = {
            'arm': arm,
            'library': lib,
            'column_in_tables': col,
            'fastq_records': int(p.get('fastq_records', 0) or 0),
            'prep_records': int(p.get('prep_records', 0) or 0),
            'star_input': s.get('star_input', 0),
            'star_uniq': s.get('star_uniq', 0),
            'star_multi': s.get('star_multi', 0),
            'star_avg_input_len': s.get('star_avg_input_len', 0),
            'star_avg_mapped_len': s.get('star_avg_mapped_len', 0.0),
            'star_tooshort_pct': s.get('star_tooshort_pct', 'NA'),
            'single_bed_rows': int(a.get('single_bed_rows', 0) or 0),
            'multi_bed_rows': int(a.get('multi_bed_rows', 0) or 0),
            'table_reads': int(tot_r[col].sum()),
            'genes_detected': int((tot_r[col] > 0).sum()),
            'vasalen_truncated': int(v.get('truncated', 0) or 0),
            'vasalen_already_shorter': int(v.get('already_shorter_than_draw', 0) or 0),
            'vasalen_mean_len_out': float(v.get('mean_len_out', 0) or 0),
        }
        rows.append(r)
    rec = pd.DataFrame(rows)

    # Derived ratios, each one a named question rather than a bare number.
    rec['trim_kept_frac'] = rec.prep_records / rec.fastq_records.replace(0, np.nan)
    rec['star_input_eq_prep'] = rec.star_input == rec.prep_records
    rec['uniq_frac_of_input'] = rec.star_uniq / rec.star_input.replace(0, np.nan)
    rec['bed_rows_per_uniq_read'] = rec.single_bed_rows / rec.star_uniq.replace(0, np.nan)
    rec['assigned_frac_of_input'] = rec.table_reads / rec.star_input.replace(0, np.nan)
    rec['assigned_frac_of_uniq'] = rec.table_reads / rec.star_uniq.replace(0, np.nan)

    emit('%-11s %12s %12s %12s %12s %13s %12s' %
         ('library', 'FASTQ', 'trimmed', 'STAR in', 'uniq', 'BED rows', 'assigned'))
    for _, r in rec.iterrows():
        emit('%-11s %12d %12d %12d %12d %13d %12d' %
             (r.library, r.fastq_records, r.prep_records, r.star_input,
              r.star_uniq, r.single_bed_rows, r.table_reads))
    emit()
    emit('%-11s %10s %10s %10s %12s %10s' %
         ('library', 'kept', 'in==prep', 'uniq/in', 'rows/uniq', 'asgn/in'))
    for _, r in rec.iterrows():
        emit('%-11s %9.4f %10s %10.4f %12.4f %10.4f' %
             (r.library, r.trim_kept_frac, str(r.star_input_eq_prep),
              r.uniq_frac_of_input, r.bed_rows_per_uniq_read,
              r.assigned_frac_of_input))

    # The links that must be exact.
    for _, r in rec.iterrows():
        if not r.star_input_eq_prep:
            fail.append('%s: STAR input (%d) != prep records (%d) -- STAR did not read '
                         'the whole fastq' % (r.library, r.star_input, r.prep_records))
        if r.table_reads > r.star_input:
            fail.append('%s: table reads (%d) EXCEED STAR input (%d) -- real double count'
                        % (r.library, r.table_reads, r.star_input))
        if r.fastq_records and r.prep_records > r.fastq_records:
            fail.append('%s: trimmed records (%d) exceed delivered records (%d)'
                        % (r.library, r.prep_records, r.fastq_records))
    emit()
    emit('EXPECTED inequalities, stated so they are not read as breaks:')
    emit('  BED rows > uniquely mapped reads -- one read overlapping k annotated')
    emit('    features gives k rows; step 6 groups the BED by read name and assigns')
    emit('    each read ONCE. rows/uniq is features per read, not a double count.')
    emit('  assigned < STAR input -- drops unmapped reads, multimapper-only reads')
    emit('    that lost their tie-break, and reads over no annotated feature.')
    emit('  trimmed <= delivered -- cutadapt -m %s drops reads that fall under it.'
         % os.environ.get('FSV_TRIM_MINLEN', '20'))

    if arm == 'vasalen':
        emit()
        emit('vasalen hard trim, per library:')
        emit('%-11s %14s %14s %14s' %
             ('library', 'truncated', 'already shorter', 'mean len out'))
        for _, r in rec.iterrows():
            emit('%-11s %14d %14d %14.2f' %
                 (r.library, r.vasalen_truncated, r.vasalen_already_shorter,
                  r.vasalen_mean_len_out))
        emit('"already shorter" = the read was shorter than its drawn target, so it')
        emit('was kept whole. The output distribution is therefore the pointwise')
        emit('MINIMUM of input and target: it can only be shorter than VASA, never')
        emit('longer, which is the conservative direction for this control.')

    # ------------------------------------------------- UFI / transcript checks
    section('2. UFI IS A MASK, AND TranscriptCounts CARRIES NO ABUNDANCE (finding 3f)')
    if tot_u is None or tot_t is None:
        fail.append('UFI or Transcript table missing')
    else:
        vals = sorted(pd.unique(tot_u.values.ravel()))
        emit('distinct values in _total.UFICounts : %s%s'
             % (vals[:12], ' ...' if len(vals) > 12 else ''))
        if set(vals) <= {0, 1}:
            emit('CONFIRMED: UFICounts holds only 0/1 -- a DETECTION MASK, not a')
            emit('molecule count. One literal UMI key "A" per gene per column.')
        else:
            fail.append('UFICounts holds values outside {0,1}: %s' % vals[:12])
        same = np.allclose(tot_u.values.astype(float), tot_t.values.astype(float))
        emit('TranscriptCounts == UFICounts elementwise : %s' % same)
        if same:
            emit('CONFIRMED: bc2trans at K = 4**len("A") = 4 maps 1 -> 1.0 exactly, so')
            emit('the TranscriptCounts tables ARE the mask. Neither carries abundance')
            emit('information. ReadCounts is the quantification; carry only that.')
        else:
            fail.append('TranscriptCounts != UFICounts -- check K and len(umi)')
        d = int(tot_u.values.sum() - (tot_r.values > 0).sum())
        emit('sum(UFI) - count(reads>0) = %d   (0 confirms UFI is exactly the mask)' % d)
        if d != 0:
            warn.append('sum(UFI) - count(reads>0) = %d, expected 0' % d)
        emit()
        emit('NO UMI CEILING EXISTS ON THIS PATH. bc2trans clamps at x >= K; K = 4')
        emit('here, and UFI per gene per column is at most 1, so the clamp is')
        emit('unreachable. VASA dropped 8 genes at its K = 4096 ceiling; that filter')
        emit('has no analogue here and must not be carried over.')
        emit('max UFI observed: %d (clamp would need %d)' % (int(tot_u.values.max()), 4))

    # ------------------------------------------------ spliced/unspliced excess
    section('3. SPLICED + UNSPLICED vs TOTAL (finding 3h)')
    emit('The no-UMI branches test SUBSTRING containment ("exon" in k / "intron" in')
    emit('k) where the vasa branch tests exact membership, so a combination label')
    emit('"exon-intron" is counted in BOTH. Excess was 0 in the A9 dry run; that is')
    emit('not a guarantee, so it is measured here per library.')
    ur, sr, xr = (rd('_uniaggGenes_total.ReadCounts.tsv'),
                  rd('_uniaggGenes_spliced.ReadCounts.tsv'),
                  rd('_uniaggGenes_unspliced.ReadCounts.tsv'))
    excess = {}
    if ur is None or sr is None or xr is None:
        warn.append('uniaggGenes read tables missing -- cannot check the exon/intron split')
    else:
        emit()
        emit('%-11s %14s %14s %14s %14s %10s %9s' %
             ('library', 'total', 'spliced', 'unspliced', 'spl+unspl', 'excess', 'excess%'))
        for col in ur.columns:
            lib = col2lib.get(col, str(col).rsplit('/')[-1])
            t, s_, x = int(ur[col].sum()), int(sr[col].sum()), int(xr[col].sum())
            exc = s_ + x - t
            excess[lib] = exc
            emit('%-11s %14d %14d %14d %14d %10d %8.4f%%' %
                 (lib, t, s_, x, s_ + x, exc, 100.0 * exc / max(1, t)))
            if exc > 0:
                warn.append('%s: spliced+unspliced exceeds total by %d reads (%.4f%%) '
                            '-- the substring double count is live in this library'
                            % (lib, exc, 100.0 * exc / max(1, t)))
        emit()
        emit('unspliced fraction (unspliced/total, uni-genes, reads):')
        for col in ur.columns:
            t = int(ur[col].sum())
            emit('  %-11s %.4f' % (col2lib.get(col, col), int(xr[col].sum()) / max(1, t)))
        # mapStats names the mechanism directly: a non-zero multi-label line is
        # where the excess can come from.
        ms = os.path.join(outdir, '%s_mapStats.log' % sample)
        if os.path.exists(ms):
            for line in open(ms):
                if 'multiple labels' in line:
                    emit('  mapStats: %s' % line.strip())

    # -------------------------------------------------------- short biotypes
    section('4. SHORT / NON-POLY-A BIOTYPES -- the axis the two arms exist to separate')
    bt = pd.Series([biotype_of(i) for i in tot_r.index], index=tot_r.index)
    emit('Rows whose components disagree on biotype are reported as MULTI:<set>,')
    emit('not assigned to one of them. Counts below are SINGLE-biotype rows only,')
    emit('so they are a lower bound on each species and are comparable between arms.')
    emit()
    bio_tbl = []
    for b in SHORT_BIOTYPES:
        sel = bt == b
        if not sel.any():
            bio_tbl.append((b, 0, 0, 0))
            continue
        sub = tot_r.loc[sel]
        bio_tbl.append((b, int(sel.sum()), int(sub.values.sum()),
                        int((sub.sum(axis=1) > 0).sum())))
    emit('%-12s %10s %14s %16s' % ('biotype', 'rows', 'reads (all libs)', 'rows detected'))
    for b, nrow, nread, ndet in bio_tbl:
        emit('%-12s %10d %14d %16d' % (b, nrow, nread, ndet))

    emit()
    emit('per library, reads in each short biotype:')
    hdr = '%-11s' % 'library' + ''.join('%11s' % b for b in SHORT_BIOTYPES)
    emit(hdr)
    perlib_bio = {}
    for col in tot_r.columns:
        lib = col2lib[col]
        vals = []
        for b in SHORT_BIOTYPES:
            sel = bt == b
            vals.append(int(tot_r.loc[sel, col].sum()) if sel.any() else 0)
        perlib_bio[lib] = vals
        emit('%-11s' % lib + ''.join('%11d' % v for v in vals))

    # The dedicated tRNA table step 7 writes, which is collapsed by isotype and
    # is the table the VASA side's tRNA claim uses.
    trna = rd('_tRNA.ReadCounts.tsv')
    emit()
    if trna is None:
        warn.append('no _tRNA.ReadCounts.tsv')
    else:
        emit('_tRNA.ReadCounts.tsv (step 7, collapsed by isotype): %d rows' % len(trna))
        if len(trna) == 0:
            emit('EMPTY. At native length that is the expected structural result:')
            emit('jS:IN requires containment and 98.5%% of tRNA features are under')
            emit('151 nt. Whether the vasalen arm recovers it is the whole question.')
        else:
            emit('total reads: %d' % int(trna.values.sum()))
            emit('%-28s %12s' % ('isotype', 'reads'))
            for i, v in trna.sum(axis=1).sort_values(ascending=False).head(15).items():
                emit('%-28s %12d' % (str(i)[:28], int(v)))

    # ----------------------------------------------------- biotype composition
    section('5. BIOTYPE COMPOSITION -- poly-A primed, so protein-coding should win')
    ribo_fams = ('rRNA', 'MtRrna', 'rRNApseudogene')
    for col in tot_r.columns:
        g = tot_r[col].groupby(bt).sum().sort_values(ascending=False)
        tot = int(g.sum())
        ribo = int(g[[b for b in g.index if b in ribo_fams]].sum()) if len(g) else 0
        pc = int(g.get('ProteinCoding', 0))
        emit()
        emit('%s -- %d reads, %d rRNA-family (%.2f%%), ProteinCoding %.2f%% of the '
             'non-rRNA remainder' % (col2lib[col], tot, ribo,
                                     100.0 * ribo / max(1, tot),
                                     100.0 * pc / max(1, tot - ribo)))
        for b, v in g.head(8).items():
            emit('    %-44s %12d %7.2f%%' % (str(b)[:44], int(v), 100.0 * v / max(1, tot)))
        if 100.0 * pc / max(1, tot - ribo) < 40:
            warn.append('%s: protein-coding only %.1f%% of the non-rRNA remainder'
                        % (col2lib[col], 100.0 * pc / max(1, tot - ribo)))
    emit()
    emit('NOTE: rRNA here is the ANNOTATION route and is NOT the rRNA figure for')
    emit('this comparison. Use res/flashseq/rrna_bwa.tsv -- the annotation route')
    emit('reads low (0.80%% in the dry run vs 3.50-6.44% by bwa) because 99.7% of')
    emit('rRNA features are under 151 nt and Ensembl lacks the 47S unit.')
    emit()
    emit('ASYMMETRY WORTH STATING when these rows are compared with the VASA side:')
    emit('the FLASH-seq rows above are clean of one specific artefact that the VASA')
    emit('rows are not. FLASH-seq runs stranded=n, so an antisense ribosomal read is')
    emit('caught by step 3 on its own side. VASA runs stranded=y (correctly, for a')
    emit('stranded protocol), and this library is measurably less strand-specific')
    emit('than the published one -- 57.7-83.8% of its ribosomal reads are forward,')
    emit('median 74.8%, against 92.7-97.4% published (measured on both plates with')
    emit('riboread-selection.py\'s own two predicates). So ~24% of VASA ribosomal')
    emit('reads are sent to the non-ribosomal arm and go on to be gene-assigned;')
    emit('vasaplate_check/README.md documents that failure mode inflating one')
    emit('antisense rRNA locus ~600x. Treat any rRNA / MiscRna / MtRrna row on the')
    emit('VASA side as possibly carrying antisense rRNA leakage rather than as a')
    emit('clean measurement. It does not affect the FLASH-seq rows above.')

    # ------------------------------------------------------------------ write
    for b in SHORT_BIOTYPES:
        rec['reads_' + b] = [perlib_bio[l][SHORT_BIOTYPES.index(b)] for l in rec.library]
    rec['spliced_unspliced_excess'] = [excess.get(l, np.nan) for l in rec.library]
    os.makedirs(os.path.dirname(tsvpath) or '.', exist_ok=True)
    rec.to_csv(tsvpath, sep='\t', index=False)
    emit()
    emit('written: %s' % tsvpath)

    section('VERDICT')
    if fail:
        emit('FAIL (%d):' % len(fail))
        for f in fail:
            emit('  x %s' % f)
    else:
        emit('Reconciliation holds for all %d libraries and the no-UMI semantics are'
             % len(rec))
        emit('as predicted. ReadCounts is the quantification.')
    if warn:
        emit()
        emit('WARN (%d) -- not crashes, but they change what the numbers MEAN:' % len(warn))
        for w in warn:
            emit('  ! %s' % w)

    with open(reportpath, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % reportpath, flush=True)
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
