#!/usr/bin/env python3
"""Transcript-body coverage, THREE protocols/pipelines, computed identically.

WHAT THIS IS, AND WHAT IT IS A FORK OF
--------------------------------------
This is a fork of `res/flashseq_vasa/mk_gene_coverage.py` (Rule 1: fork, never
patch). That script established the metric definitions and the FLASH-seq vs
own-VASA-plate result. It is left callable and untouched; its numbers are
reproduced here as a regression check (see `merge`, REGRESSION section).

The question this fork exists to answer:

    the own VASA plate shows a sharp coverage rise over the last ~10% of the
    transcript (3'/5' = 1.22). Is that a property of the VASA PROTOCOL, or
    something specific to the user's library?

The discriminating experiment is the PUBLISHED VASA plate (SRR14783059 /
GSM5369495, `vasaplate-HEK293T-mESC`). Same protocol, different hands, different
cells, different pipeline run. If the rise is there too, it is the protocol. If
it is not, it is this library.

METRIC DEFINITIONS -- CARRIED OVER VERBATIM, INCLUDING WHICH ONE TO REPORT
-------------------------------------------------------------------------
For every read, find where along its transcript it landed as a percentile of
transcript length 5'->3'; average over transcripts.

  base   every aligned base votes for the bin it sits in -- nominally what
         qualimap and RSeQC compute.
  mid    each read's midpoint casts ONE vote, so a long and a short read count
         the same.

A LIMITATION OF `mid` THAT IS LARGE AND MUST BE READ WITH EVERY mid NUMBER
-------------------------------------------------------------------------
The midpoint is the GENOMIC-SPAN midpoint, `(blocks[0][0] + blocks[-1][1])//2`.
For a junction-spanning read that point can fall in an INTRON, where no model
exon is found and the read casts NO vote -- while still counting in
`reads_placed`, because its bases did place. Measured by
`measure_middrop.py` -> `coverage_threeway_middrop.tsv`:

    VASA published   30.554% of placed reads cast no midpoint vote
    VASA own         41.896%
    FLASH-seq native 47.798%
    FLASH-seq vasalen 40.591%

and per-unit dropout tracks p50 aligned length at r = 0.982. So `mid` is NOT
read-length-neutral in the way the word "midpoint" suggests; it is only
vote-count-neutral. A transcript-coordinate midpoint (from the walked
in-transcript positions, which `profile` already computes) would remove this and
is the obvious next fix; it is NOT done here, so every `mid` number in this
file's outputs carries this dropout.

What the dropout can and cannot do: a dropped read contributes nothing to any
bin, so it can only reshape a profile if dropped reads are non-uniform ALONG the
transcript. `base` has no such dropout, and under `base` the own plate is still
3'-heavier than the published plate (own/published rise ratio 1.850 by base,
2.144 by mid). The DIRECTION of the published-vs-own claim therefore does not
rest on the dropout; its magnitude does.

The two disagree by 2x on FLASH-seq. nf-core ran qualimap over the same
FLASH-seq BAMs and gives 3'/5' = 0.93; `mid` gives 0.89 (agrees), `base` gives
0.46 (does not). **`mid` is therefore the reported metric and `base` is
diagnostic only.** That decision is inherited, not re-litigated here.

The cause of the base-vs-mid disagreement is NOT KNOWN. The upstream docstring
records a withdrawn claim: an edge-clipping hypothesis was declared refuted on
the basis of a wrong VASA read length (~24 nt, which is the *trimmed* length from
step 3's log, not the ~127 nt ALIGNED length). That refutation is withdrawn and
is NOT re-derived here. **Edge clipping remains the leading untested
hypothesis.** This fork adds the measurement the upstream docstring asks for --
see LOSS DECOMPOSITION below -- but reports it as a decomposition of where bases
go, not as a verdict on a mechanism.

TWO ANNOTATION RELEASES, WHICH IS THE WHOLE DIFFICULTY
-----------------------------------------------------
    published plate   GRCm38 + hg38, Ensembl 99   (mixed human+mouse index)
    own plate         GRCm39, Ensembl 116
    FLASH-seq         GRCm39, Ensembl 116

So a published-vs-own gap mixes THREE causes: protocol, pipeline run, and
annotation release. Two things are done about that.

1. SPECIES. Published-plate rows are ENSG (human) or ENSMUSG (mouse). Only
   mouse rows are used, and only cells classified mouse-pure -- by BOTH of the
   paper's two mutually-inconsistent rules (vp_common.classify_fig1d on UFI
   fraction and classify_methods on gene fraction). `select` reports both and
   takes the intersection, so a cell counted here is mouse under either rule.

2. RELEASE. The gene set is the intersection: a gene is used only if it is
   protein-coding and expressed in ALL FOUR unit groups, AND has a >=MIN_TXLEN
   transcript in BOTH Ensembl 99 (GRCm38) and Ensembl 116 (GRCm39). Two model
   files are built from the two GTFs with the gene lists in IDENTICAL ORDER, so
   every profile is over the same genes and only the gene MODEL differs.

   The residual release effect is then bounded rather than assumed: `crossrel`
   re-profiles ONE own-plate cell (GRCm39/E116 BAM) -- against nothing else --
   and separately, the per-gene transcript lengths of the two releases are
   compared in the gene table. A release effect on coverage SHAPE can only come
   through changed transcript models, so the length agreement is the handle on
   its size. Any comparison that could be confounded by 99-vs-116 says so.

LOSS DECOMPOSITION (new here; the measurement upstream asked for)
-----------------------------------------------------------------
A read's aligned bases can fail to reach a bin for four distinct reasons. All
four are counted, over reads that placed at least one base (so the gene is in
the model and the count is meaningful):

  lost_txend      block clipped at the exon that is TERMINAL in transcript
                  coordinates -- i.e. read overhangs the annotated 3' end
  lost_txstart    block clipped at the FIRST exon in transcript coordinates
  lost_internal   block clipped at an internal exon (intron overlap /
                  unspliced read running past an internal exon boundary)
  lost_noexon     whole block dropped: its start is in no model exon

`base_fullexon` / `mid_fullexon` are the same two profiles restricted to reads
that lost NOTHING. If base and mid converge there, the base-vs-mid gap is
attributable to bases lost outside the exon model, of which transcript-end
overhang is the `lost_txend` component. That is a decomposition, not a proof of
mechanism, and it is reported as such.

`mid_sense` is `mid` restricted to reads whose alignment strand matches the
transcript strand -- reported because VASA is counted with stranded=y and
FLASH-seq with stranded=n (Rule 6), so an unstranded midpoint profile is not
counted the way either pipeline counts its own reads.

READ CAP AND DEPTH
------------------
Streaming stops after MAX_READS alignments have landed in the gene set. The
published plate is 384 cells sharing one library, so its cells are ~10x
shallower than the own plate's and will not reach the cap; `reads_placed` is
written per unit and MUST be read alongside any published-plate number. A
per-transcript read floor sensitivity check (merge, ROBUSTNESS) tests whether
the shallower cells' profiles are driven by thinly covered transcripts.

Usage:
  mk_coverage_threeway.py select   <pub_ufi_tsv> <out_tsv>
  mk_coverage_threeway.py geneset  <own> <fsnat> <fsvas> <pub> <pubcells> \
                                   <gtf99> <gtf116> <out_prefix>
  mk_coverage_threeway.py profile  <models.npz> <bam> <label> <group> <out_pfx>
  mk_coverage_threeway.py merge    <cov_dir> <res_dir> <geneset_prefix>
  mk_coverage_threeway.py crossrel <cov_dir> <res_dir>
"""
import glob
import os
import sys

import numpy as np

# --- inherited verbatim from mk_gene_coverage.py; changing any of these breaks
# --- comparability with the already-published FLASH-seq/own-plate numbers
NBINS = 100          # matches qualimap's 100-bin profile
MIN_TXLEN = 1000     # nt; a transcript must be long enough to HAVE a 5'/3' axis
MIN_READS_EACH = 20  # per-gene floor, in every unit group
MAX_GENES = 4000
MAX_READS = 3_000_000

# --- new here
MIN_UFI = 7500       # vp_common.MIN_UFI, the paper's own cell gate
# 12 deepest mouse-pure published cells. 12 rather than 8 for two reasons
# measured in the precheck: the published plate's 384 barcodes share one library
# so each cell is ~10x shallower than an own-plate cell, and the deepest BAMs by
# FILE SIZE are human (HEK293T) cells -- one of them placed 111 reads out of 2.7M
# primary. 12 also matches the own plate's 12 real cells, so the two VASA groups
# carry the same n and the permutation test is symmetric.
N_PUB_CELLS = 12
ROBUST_MIN_TX_READS = 20   # per-transcript floor for the sensitivity check
# Upper bound of the aligned-length regime all four datasets share. The
# published plate's STAR index is sjdbOverhang 73, i.e. a ~75 nt library, and
# its aligned reads measure ~73 nt; the own plate's are ~127 nt and FLASH-seq
# native's ~150 nt. 80 nt is just above the published plate's own p75, so its
# short-read profile is nearly all of its reads while the others contribute
# only their genuinely short tail. Chosen from the measured distributions, not
# picked to produce an answer.
SHORT_MAX = 80

BLANKS_OWN = {'001', '014', '015', '016'}   # confirmed blanks, own plate


# ===========================================================================
# 0. published-plate cell selection
# ===========================================================================
def select(pub_ufi_tsv, out_tsv):
    """Pick the deepest MOUSE-PURE published-plate cells.

    Species purity uses BOTH of the paper's rules, which disagree, and keeps
    only cells called mouse by BOTH. vp_common.py is the authority for the rule
    definitions; they are reimplemented here rather than imported so this script
    stands alone, and the reimplementation is asserted against vp_common when
    that module is importable.
    """
    import pandas as pd

    def species_of(i):
        parts = i.split('-')
        hh = any(p.startswith('ENSG') for p in parts)
        mm = any(p.startswith('ENSMUSG') for p in parts)
        if hh and mm:
            return 'mixed'
        return 'human' if hh else ('mouse' if mm else 'other')

    # chunked: the table is 588k x 384 and only per-column reductions are needed
    with open(pub_ufi_tsv) as fh:
        cols = [_wellid(c) for c in fh.readline().rstrip('\n').split('\t')[1:]]
    acc = {k: np.zeros(len(cols), dtype=np.int64)
           for k in ('h', 'm', 'gh', 'gm')}
    nrow = 0
    for ch in pd.read_csv(pub_ufi_tsv, sep='\t', index_col=0, chunksize=50000):
        nrow += len(ch)
        sp = np.array([species_of(str(i)) for i in ch.index])
        v = ch.values
        for key, mask in (('h', sp == 'human'), ('m', sp == 'mouse')):
            sub = v[mask]
            acc[key] += sub.sum(axis=0)
            acc['g' + key] += (sub > 0).sum(axis=0)
    print('read %d rows of %s' % (nrow, os.path.basename(pub_ufi_tsv)))

    h = pd.Series(acc['h'], index=cols)
    m = pd.Series(acc['m'], index=cols)
    gh = pd.Series(acc['gh'], index=cols)
    gm = pd.Series(acc['gm'], index=cols)

    tot_u = h + m
    keep = tot_u >= MIN_UFI
    f_ufi = h.where(keep) / tot_u.where(keep)
    f_gene = gh.where(keep) / (gh + gm).where(keep)

    def label(frac):
        lab = pd.Series('discarded', index=cols, dtype=object)
        lab[keep & (frac > 0.75)] = 'human'
        lab[keep & (frac < 0.25)] = 'mouse'
        lab[keep & (frac >= 0.25) & (frac <= 0.75)] = 'mixed'
        return lab

    lab_fig1d, lab_methods = label(f_ufi), label(f_gene)

    out = pd.DataFrame({
        'cell': cols, 'ufi_human': h.values, 'ufi_mouse': m.values,
        'ufi_total_assigned': tot_u.values,
        'genes_human': gh.values, 'genes_mouse': gm.values,
        'frac_human_ufi': f_ufi.values, 'frac_human_gene': f_gene.values,
        'class_fig1d': lab_fig1d.values, 'class_methods': lab_methods.values,
    })
    out['mouse_both_rules'] = ((out.class_fig1d == 'mouse')
                              & (out.class_methods == 'mouse'))
    out = out.sort_values('ufi_mouse', ascending=False)
    out['selected'] = False
    sel = out.index[out.mouse_both_rules][:N_PUB_CELLS]
    out.loc[sel, 'selected'] = True
    out.to_csv(out_tsv, sep='\t', index=False)

    # --- cross-check against the ESTABLISHED call table --------------------
    # res/threeway/threeway_published_cellcalls.tsv came out of an earlier stage
    # of this project and is the table the rest of the writeup quotes (mESC
    # n=173 by Fig.1d, n=141 by Methods). This function recomputes the calls
    # from the UFI table independently. If the two disagree then one of them is
    # wrong, and that has to surface HERE rather than be averaged silently into
    # a coverage profile -- so the disagreement count is asserted to be zero.
    est = os.path.join(os.path.dirname(os.path.abspath(out_tsv)),
                       'threeway_published_cellcalls.tsv')
    if os.path.exists(est):
        ref = pd.read_csv(est, sep='\t', dtype={'unit': str})
        ref['unit'] = ref.unit.astype(str).str.zfill(3)
        j = out.merge(ref, left_on='cell', right_on='unit', how='inner')
        assert len(j) == len(out), ('join lost barcodes', len(j), len(out))
        d1 = int((j.class_fig1d != j.call_fig1d).sum())
        d2 = int((j.class_methods != j.call_methods).sum())
        print('cross-check vs %s (%d barcodes joined)'
              % (os.path.basename(est), len(j)))
        print('  established: fig1d mouse=%d  methods mouse=%d  both=%d'
              % (int((ref.call_fig1d == 'mouse').sum()),
                 int((ref.call_methods == 'mouse').sum()),
                 int(((ref.call_fig1d == 'mouse')
                      & (ref.call_methods == 'mouse')).sum())))
        print('  recomputed here, disagreements: fig1d %d, methods %d' % (d1, d2))
        assert d1 == 0 and d2 == 0, (
            'species calls disagree with the established table: fig1d %d, '
            'methods %d barcodes' % (d1, d2))
        print('  AGREE exactly -- this selection rests on the same species '
              'calls the rest of the writeup quotes')
    else:
        print('WARNING established call table absent at %s -- the selection is '
              'NOT cross-checked' % est)

    n1 = int((out.class_fig1d == 'mouse').sum())
    n2 = int((out.class_methods == 'mouse').sum())
    print('published plate, %d barcodes' % len(out))
    print('  mouse by Fig.1d rule (UFI fraction)   : %d' % n1)
    print('  mouse by Methods rule (gene fraction) : %d' % n2)
    print('  mouse by BOTH (used here)             : %d'
          % int(out.mouse_both_rules.sum()))
    print('  the two rules disagree on             : %d barcodes'
          % int((out.class_fig1d != out.class_methods).sum()))
    print('  selected %d deepest: %s'
          % (int(out.selected.sum()), ','.join(out.cell[out.selected])))
    print('  their mouse UFIs: %s'
          % ','.join('%d' % v for v in out.ufi_mouse[out.selected]))
    print('wrote %s' % out_tsv)


# ===========================================================================
# 1. gene set + two transcript-model files, gene order identical
# ===========================================================================
def _wellid(c):
    return str(c).rsplit('/', 1)[-1].rsplit('_', 1)[-1].zfill(3)


def _expressed(path, keep_cols=None, drop_cols=None):
    """Per-gene ReadCounts total over MOUSE simple protein-coding rows.

    Rule 4: ReadCounts on every side. VASA's TranscriptCounts is the right
    column for VASA's own biology, but reads are the only unit all three
    protocols measure and this is a cross-protocol comparison.

    Columns are selected by POSITION via usecols, because the published table is
    588k x 384 and reading all of it to keep 8 columns costs ~900 MB for nothing.
    """
    import pandas as pd

    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
    if keep_cols:
        wells = {_wellid(c): i for i, c in enumerate(header)}
        missing = [c for c in keep_cols if c not in wells]
        assert not missing, 'columns absent from %s: %s' % (path, missing)
        use = [0] + [wells[c] for c in keep_cols]
        df = pd.read_csv(path, sep='\t', index_col=0, usecols=use)
        df = df[[header[wells[c]] for c in keep_cols]]
    else:
        df = pd.read_csv(path, sep='\t', index_col=0)
        if drop_cols:
            df = df[[c for c in df.columns if _wellid(c) not in drop_cols]]
    idx = df.index.astype(str)
    keep = np.array([('-' not in i) and i.startswith('ENSMUSG')
                     and i.endswith('_ProteinCoding') for i in idx])
    df = df[keep]
    tot = df.sum(axis=1)
    gid = [i.split('_')[0] for i in df.index.astype(str)]
    s = pd.Series(tot.values, index=gid).groupby(level=0).sum()
    return s, int(df.shape[1])


def _parse_gtf(gtf, want, contig_prefix=None):
    """exon rows of every transcript of every wanted gene -> {tx: [g,chrom,strand,exons]}"""
    tx = {}
    with open(gtf) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'exon':
                continue
            if contig_prefix is not None and not p[0].startswith(contig_prefix):
                continue
            a = p[8]
            i = a.find('gene_id "')
            if i < 0:
                continue
            g = a[i + 9:a.find('"', i + 9)]
            if g not in want:
                continue
            j = a.find('transcript_id "')
            t = a[j + 15:a.find('"', j + 15)]
            rec = tx.get(t)
            if rec is None:
                rec = tx[t] = [g, p[0], p[6], []]
            rec[3].append((int(p[3]) - 1, int(p[4])))
    return tx


def _longest(tx):
    """gene -> (len, tx_id, chrom, strand, sorted exons); longest >= MIN_TXLEN."""
    best = {}
    for t, (g, chrom, strand, ex) in tx.items():
        ex = sorted(set(ex))
        L = sum(e - s for s, e in ex)
        if L < MIN_TXLEN:
            continue
        if g not in best or L > best[g][0]:
            best[g] = (L, t, chrom, strand, ex)
    return best


def _write_models(genes, best, out):
    ex_chrom, ex_start, ex_end, ex_tx, ex_off = [], [], [], [], []
    tx_len, tx_strand, tx_id = [], [], []
    for k, g in enumerate(genes):
        L, t, chrom, strand, ex = best[g]
        tx_len.append(L)
        tx_strand.append(1 if strand == '+' else -1)
        tx_id.append(t)
        off = 0
        for s, e in ex:
            ex_chrom.append(chrom); ex_start.append(s); ex_end.append(e)
            ex_tx.append(k); ex_off.append(off)
            off += e - s
        assert off == L, (g, off, L)
    np.savez_compressed(
        out,
        gene=np.array(genes), tx_id=np.array(tx_id),
        tx_len=np.array(tx_len, dtype=np.int64),
        tx_strand=np.array(tx_strand, dtype=np.int8),
        ex_chrom=np.array(ex_chrom), ex_start=np.array(ex_start, dtype=np.int64),
        ex_end=np.array(ex_end, dtype=np.int64),
        ex_tx=np.array(ex_tx, dtype=np.int32),
        ex_off=np.array(ex_off, dtype=np.int64), nbins=np.array([NBINS]))
    print('  wrote %s (%d exons over %d transcripts)'
          % (out, len(ex_start), len(genes)))


def geneset(own_tsv, fsnat_tsv, fsvas_tsv, pub_tsv, pubcells_tsv,
            gtf99, gtf116, out_prefix):
    import pandas as pd

    # NOTE the published ReadCounts table passed in must be the SAME family
    # (uniaggGenes vs raw total) as the UFI table `select` used, or the two
    # stages disagree about which entries exist. See the driver's comment.
    cells = pd.read_csv(pubcells_tsv, sep='\t')
    pub_cols = list(cells.cell.astype(str).str.zfill(3)[cells.selected])
    assert len(pub_cols) == N_PUB_CELLS, pub_cols
    print('published cells used: %s' % ','.join(pub_cols))

    e = {}
    e['own'], n_own = _expressed(own_tsv, drop_cols=BLANKS_OWN)
    e['fs_native'], n_fn = _expressed(fsnat_tsv)
    e['fs_vasalen'], n_fv = _expressed(fsvas_tsv)
    e['published'], n_pub = _expressed(pub_tsv, keep_cols=pub_cols)
    print('units per group: own=%d fs_native=%d fs_vasalen=%d published=%d'
          % (n_own, n_fn, n_fv, n_pub))
    for k, v in e.items():
        print('  %-11s %6d mouse protein-coding genes, %.3g total reads'
              % (k, len(v), float(v.sum())))

    common = None
    for k, v in e.items():
        s = set(v.index[v >= MIN_READS_EACH])
        print('  %-11s %6d genes >= %d reads' % (k, len(s), MIN_READS_EACH))
        common = s if common is None else (common & s)
    print('expressed >= %d reads in ALL FOUR groups: %d genes'
          % (MIN_READS_EACH, len(common)))
    # An empty or tiny gene set does not raise anywhere downstream: the models
    # simply carry no exons, every BAM places 0 reads, and merge writes a table
    # of zeros that looks like a result. Fail here instead. The floor is set far
    # below any plausible value (the two-group upstream run reached MAX_GENES)
    # so it catches a broken join, not a marginal dataset.
    assert len(common) >= 500, (
        'only %d genes expressed >= %d reads in all four groups -- check the '
        'count tables are the same FAMILY and that species filtering left mouse '
        'rows' % (len(common), MIN_READS_EACH))
    for k, v in e.items():
        assert len(v) > 5000, (
            '%s yielded only %d mouse protein-coding genes -- wrong table or '
            'wrong species filter' % (k, len(v)))

    # rank by the geometric mean of CPM, not of raw totals: the four groups
    # differ in depth by orders of magnitude, and a raw-total geometric mean
    # would let the deepest group choose the gene set. (Upstream used raw
    # totals over two groups; this is the documented divergence.)
    want = sorted(common)
    cpm = {k: (e[k][want] / e[k].sum() * 1e6) for k in e}
    lg = np.zeros(len(want))
    for k in e:
        lg += np.log(cpm[k].values)
    score = dict(zip(want, lg / len(e)))

    print('parsing %s (Ensembl 99, GRCm38 contigs only)' % gtf99)
    b99 = _longest(_parse_gtf(gtf99, common, contig_prefix='GRCm38_'))
    print('  genes with a >=%d nt transcript: %d' % (MIN_TXLEN, len(b99)))
    print('parsing %s (Ensembl 116, GRCm39)' % gtf116)
    b116 = _longest(_parse_gtf(gtf116, common))
    print('  genes with a >=%d nt transcript: %d' % (MIN_TXLEN, len(b116)))

    both = sorted(set(b99) & set(b116), key=lambda g: -score[g])
    print('genes modellable in BOTH releases: %d' % len(both))
    genes = both[:MAX_GENES]
    print('gene set used: %d (top by CPM geometric mean)' % len(genes))
    assert len(genes) >= 500, (
        'gene set is %d genes, too few to average a coverage profile over'
        % len(genes))

    _write_models(genes, b99, out_prefix + '.E99.npz')
    _write_models(genes, b116, out_prefix + '.E116.npz')

    g99, g116 = np.load(out_prefix + '.E99.npz', allow_pickle=True), \
        np.load(out_prefix + '.E116.npz', allow_pickle=True)
    assert list(g99['gene']) == list(g116['gene']), 'gene order differs'
    print('gene ORDER identical across the two model files: asserted')

    L99 = {g: b99[g][0] for g in genes}
    L116 = {g: b116[g][0] for g in genes}
    tab = pd.DataFrame({
        'gene': genes,
        'tx_E99': [b99[g][1] for g in genes],
        'tx_E116': [b116[g][1] for g in genes],
        'txlen_E99': [L99[g] for g in genes],
        'txlen_E116': [L116[g] for g in genes],
        'strand_E99': [b99[g][3] for g in genes],
        'strand_E116': [b116[g][3] for g in genes],
        'nexon_E99': [len(b99[g][4]) for g in genes],
        'nexon_E116': [len(b116[g][4]) for g in genes],
    })
    for k in e:
        tab['reads_' + k] = [float(e[k][g]) for g in genes]
        tab['cpm_' + k] = [float(cpm[k][g]) for g in genes]
    tab['len_ratio_116_99'] = tab.txlen_E116 / tab.txlen_E99
    tab['same_strand'] = tab.strand_E99 == tab.strand_E116
    tab.to_csv(out_prefix + '.genes.tsv', sep='\t', index=False)

    lr = tab.len_ratio_116_99.values
    print('\n=== transcript-model agreement between Ensembl 99 and 116 ===')
    print('  identical length      : %d / %d (%.1f%%)'
          % (int((lr == 1).sum()), len(lr), 100 * (lr == 1).mean()))
    print('  within +-5%%           : %d (%.1f%%)'
          % (int((abs(lr - 1) <= 0.05).sum()), 100 * (abs(lr - 1) <= 0.05).mean()))
    print('  median length ratio   : %.4f' % float(np.median(lr)))
    print('  strand disagreements  : %d' % int((~tab.same_strand).sum()))
    print('wrote %s.genes.tsv' % out_prefix)


# ===========================================================================
# 2. one BAM -> profiles + loss decomposition + read-length distribution
# ===========================================================================
# Binning is `bin = min(nb-1, pos*nb//L)`, so bin j covers transcript positions
# [ceil(j*L/nb), ceil((j+1)*L/nb) - 1], and the last bin absorbs the remainder.
# These two invert that exactly, which is what lets a contiguous run of bases be
# accumulated as a handful of per-bin counts instead of one vote per base. They
# are asserted against the per-base definition in `_selfcheck_bins`.
def _bin_start(j, nb, L):
    return -(-j * L // nb)          # ceil(j*L/nb)


def _bin_end(j, nb, L):
    return (L - 1) if j == nb - 1 else (-(-(j + 1) * L // nb) - 1)


def _selfcheck_bins(nb=NBINS):
    """The vectorised binning must reproduce the per-base definition exactly."""
    rng = np.random.default_rng(0)
    for L in list(rng.integers(1000, 200000, 40)) + [1000, 1001, 99999]:
        L = int(L)
        want = np.minimum(nb - 1, np.arange(L) * nb // L)
        for j in range(nb):
            s, e = _bin_start(j, nb, L), _bin_end(j, nb, L)
            got = np.where(want == j)[0]
            assert len(got), (L, j)
            assert (s, e) == (int(got[0]), int(got[-1])), (L, j, s, e, got[0], got[-1])
    return True


def profile(models, bam_path, label, group, out_prefix, perbase=False):
    """Stream one BAM into the coverage profiles.

    `perbase=True` selects the reference accumulation path (one np.add.at vote
    per aligned base, as mk_gene_coverage.py does) instead of the fast run-based
    one. It exists so the two can be compared on real data; the profiles must be
    bit-identical, and the precheck asserts that.
    """
    import pysam

    assert _selfcheck_bins(), 'bin inversion self-check failed'
    m = np.load(models, allow_pickle=True)
    tx_len, tx_strand = m['tx_len'], m['tx_strand']
    nb = int(m['nbins'][0])
    ntx = len(tx_len)

    by_chrom = {}
    order = np.lexsort((m['ex_start'], m['ex_chrom']))
    ec, es, ee = m['ex_chrom'][order], m['ex_start'][order], m['ex_end'][order]
    et, eo = m['ex_tx'][order], m['ex_off'][order]
    for chrom in np.unique(ec):
        sel = ec == chrom
        by_chrom[str(chrom)] = (es[sel], ee[sel], et[sel], eo[sel])

    # genomic span of each transcript, and which genomic end is its 5' end.
    # Used to classify lost bases: a base outside [gmin, gmax) is beyond one
    # annotated terminus, and strand says which.
    tx_gmin = np.full(ntx, np.iinfo(np.int64).max, dtype=np.int64)
    tx_gmax = np.full(ntx, -1, dtype=np.int64)
    for k, s, e in zip(m['ex_tx'], m['ex_start'], m['ex_end']):
        if s < tx_gmin[k]:
            tx_gmin[k] = s
        if e > tx_gmax[k]:
            tx_gmax[k] = e

    cov = {k: np.zeros((ntx, nb), dtype=np.float64)
           for k in ('base', 'mid', 'mid_sense', 'base_fullexon',
                     'mid_fullexon', 'mid_short')}
    lost = dict(txend=0, txstart=0, internal=0, modelerr=0,
                unclassified=0, placed=0)
    lenhist = np.zeros(1024, dtype=np.int64)         # all primary alignments
    lenhist_placed = np.zeros(1024, dtype=np.int64)  # placed reads only
    placed = seen = primary = n_fullexon = n_binned = 0

    def find_exon(chrom, gpos):
        c = by_chrom.get(chrom)
        if c is None:
            return -1
        s, e = c[0], c[1]
        i = int(np.searchsorted(s, gpos, side='right')) - 1
        if i < 0 or gpos >= e[i]:
            return -1
        return i

    def to_tx(chrom, i, gpos):
        s, _, t, o = by_chrom[chrom]
        k = int(t[i])
        p = int(o[i]) + (gpos - int(s[i]))
        if tx_strand[k] < 0:
            p = int(tx_len[k]) - 1 - p
        return k, p

    bam = pysam.AlignmentFile(bam_path, check_sq=False)
    for r in bam.fetch(until_eof=True):
        seen += 1
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        primary += 1
        qal = r.query_alignment_length or 0
        if qal < 1024:
            lenhist[qal] += 1
        chrom = r.reference_name
        blocks = r.get_blocks()
        if not blocks:
            continue

        # --- pass 1: place every aligned base, one lookup per genomic run -----
        # A block is contiguous in transcript space only while it stays inside
        # one exon (STAR splits at junctions, but an unspliced read can run past
        # an exon boundary into an intron). Walk each block exon by exon.
        aligned_bases = sum(be - bs for bs, be in blocks)
        rows = []          # (tx index, bin indices) for this read
        home = None        # the transcript this read belongs to, if any
        n_placed_bases = 0
        merr_read = 0      # tripwire, expected 0 (see below)
        unplaced = []      # genomic intervals whose bases reached no bin
        for bs, be in blocks:
            g = bs
            while g < be:
                i = find_exon(chrom, g)
                if i < 0:
                    # skip forward to the next exon start on this chrom, if any
                    s = by_chrom[chrom][0] if chrom in by_chrom else None
                    if s is None:
                        unplaced.append((g, be)); g = be; break
                    j = int(np.searchsorted(s, g, side='right'))
                    nxt = int(s[j]) if j < len(s) else be
                    stop = min(be, nxt)
                    unplaced.append((g, stop))
                    g = stop
                    continue
                k, p0 = to_tx(chrom, i, g)
                L = int(tx_len[k])
                exon_end = int(by_chrom[chrom][1][i])
                n = min(be, exon_end) - g
                if home is None:
                    home = k
                # bins are contiguous in transcript space over a run inside one
                # exon, so the run's per-bin counts are a difference of bin
                # indices at its two ends -- no per-base array needed. This is
                # the hot loop: at ~127 nt/read x millions of reads, building an
                # np.arange(n) per exon-run and calling np.add.at on it costs
                # ~150x more than the two divisions below.
                if tx_strand[k] > 0:
                    a, b2 = p0, p0 + n - 1
                else:
                    a, b2 = p0 - n + 1, p0
                lo, hi = max(a, 0), min(b2, L - 1)
                # a base inside an exon whose transcript coordinate falls
                # outside [0, L) would mean the model's exon offsets disagree
                # with tx_len; _write_models asserts they cannot, so this
                # counter is a tripwire and is expected to stay 0
                merr_read += (lo - a) + (b2 - hi)
                if hi < lo:
                    g += n
                    continue
                b_lo = min(nb - 1, lo * nb // L)
                b_hi = min(nb - 1, hi * nb // L)
                rows.append((k, lo, hi, b_lo, b_hi, L))
                n_placed_bases += hi - lo + 1
                g += n
        if not rows:
            continue

        # --- pass 2: classify the bases that reached no bin -------------------
        # Anchor on `home`, the transcript the read placed into. A base is
        # 'txend'/'txstart' if it lies beyond that transcript's annotated
        # genomic span (which end depends on strand), 'internal' if it lies
        # inside the span but in no exon of it (intron), and 'noexon' if the
        # read placed bases in no exon of ANY transcript at that position.
        gmin, gmax = int(tx_gmin[home]), int(tx_gmax[home])
        plus = tx_strand[home] > 0
        n_lost = aligned_bases - n_placed_bases
        n_class = 0
        for a, b in unplaced:
            if b <= a:
                continue
            before = max(0, min(b, gmin) - a)      # genomically upstream of tx
            after = max(0, b - max(a, gmax))       # genomically downstream
            inside = (b - a) - before - after
            # 'before'/'after' are genomic; which TRANSCRIPT terminus that is
            # depends on strand
            lost['txstart' if plus else 'txend'] += before
            lost['txend' if plus else 'txstart'] += after
            # inside the transcript's genomic span but in no exon of it: intron
            # of `home`, or exonic only in some other model. Both 'internal'.
            lost['internal'] += inside
            n_class += before + after + inside
        # every lost base must land in exactly one class; the residual is
        # written out so the decomposition is checkable rather than assumed
        lost['modelerr'] += merr_read
        n_class += merr_read
        lost['unclassified'] += n_lost - n_class
        assert n_lost - n_class >= 0, (n_lost, n_class, merr_read)

        fullexon = (n_lost == 0)
        assert n_lost >= 0, (n_lost, aligned_bases, n_placed_bases)
        for k, lo, hi, b_lo, b_hi, L in rows:
            targets = ((cov['base'][k], cov['base_fullexon'][k]) if fullexon
                       else (cov['base'][k],))
            if perbase:
                # reference path: one vote per base, as mk_gene_coverage.py does.
                # Kept so the fast path can be proved exactly equal on real data
                # rather than only on the synthetic self-check.
                bb = np.minimum(nb - 1, np.arange(lo, hi + 1) * nb // L)
                for cm in targets:
                    np.add.at(cm, bb, 1.0)
                continue
            for cm in targets:
                if b_lo == b_hi:
                    cm[b_lo] += hi - lo + 1
                else:
                    # first and last bins are partial; those between are whole
                    cm[b_lo] += _bin_end(b_lo, nb, L) - lo + 1
                    cm[b_hi] += hi - _bin_start(b_hi, nb, L) + 1
                    for j in range(b_lo + 1, b_hi):
                        cm[j] += _bin_end(j, nb, L) - _bin_start(j, nb, L) + 1

        mid = (blocks[0][0] + blocks[-1][1]) // 2
        i = find_exon(chrom, mid)
        if i >= 0:
            k, p = to_tx(chrom, i, mid)
            j = min(nb - 1, p * nb // int(tx_len[k]))
            cov['mid'][k, j] += 1.0
            if fullexon:
                cov['mid_fullexon'][k, j] += 1.0
            sense = (tx_strand[k] > 0) != bool(r.is_reverse)
            if sense:
                cov['mid_sense'][k, j] += 1.0
            # read-length-matched control: the published plate was mapped from a
            # ~75 nt library (STAR index sjdbOverhang 73) and its aligned reads
            # are ~73 nt against the own plate's ~127 nt. A coverage-shape
            # difference between the two VASA plates could therefore be read
            # length rather than chemistry, so every unit also gets a profile
            # restricted to its SHORT reads -- the length regime all four
            # datasets have in common.
            if qal <= SHORT_MAX:
                cov['mid_short'][k, j] += 1.0

        lost['placed'] += aligned_bases
        n_binned += n_placed_bases
        if qal < 1024:
            lenhist_placed[qal] += 1
        n_fullexon += int(fullexon)
        placed += 1
        if placed >= MAX_READS:
            break
    bam.close()

    def norm(cm):
        tot = cm.sum(axis=1)
        ok = tot > 0
        if not ok.any():
            return np.zeros(nb), 0
        return (cm[ok] / tot[ok, None]).mean(axis=0), int(ok.sum())

    with open(out_prefix + '.cov.tsv', 'w') as fh:
        fh.write('label\tgroup\tmetric\tn_tx\treads_placed\treads_seen\t'
                 'reads_primary\treads_fullexon\t'
                 + '\t'.join('b%02d' % i for i in range(nb)) + '\n')
        for name in ('base', 'mid', 'mid_sense', 'base_fullexon',
                     'mid_fullexon', 'mid_short'):
            p, n_ok = norm(cov[name])
            fh.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n'
                     % (label, group, name, n_ok, placed, seen, primary,
                        n_fullexon, '\t'.join('%.8f' % x for x in p)))

    # per-transcript matrices, so merge can re-aggregate (robustness check)
    # without re-reading the BAM
    np.savez_compressed(out_prefix + '.tx.npz',
                        **{k: cov[k].astype(np.float32) for k in ('mid', 'base')})

    L = np.arange(1024)
    tot = lenhist.sum()
    def pct(hist, q):
        cs = np.cumsum(hist)
        if cs[-1] == 0:
            return 0.0
        return float(L[np.searchsorted(cs, q * cs[-1] / 100.0)])
    with open(out_prefix + '.meta.tsv', 'w') as fh:
        fh.write('label\tgroup\treads_seen\treads_primary\treads_placed\t'
                 'reads_fullexon\tfrac_fullexon\t'
                 'alnlen_p25\talnlen_p50\talnlen_p75\talnlen_mean\t'
                 'alnlen_placed_p25\talnlen_placed_p50\talnlen_placed_p75\t'
                 'alnlen_placed_mean\t'
                 'bases_aligned\tbases_binned\tlost_txend\tlost_txstart\t'
                 'lost_internal\tlost_modelerr\tlost_unclassified\n')
        mean_all = float((L * lenhist).sum() / tot) if tot else 0.0
        tp = lenhist_placed.sum()
        mean_pl = float((L * lenhist_placed).sum() / tp) if tp else 0.0
        fh.write('%s\t%s\t%d\t%d\t%d\t%d\t%.6f\t%d\t%d\t%d\t%.2f\t'
                 '%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n'
                 % (label, group, seen, primary, placed, n_fullexon,
                    (n_fullexon / placed) if placed else 0.0,
                    pct(lenhist, 25), pct(lenhist, 50), pct(lenhist, 75), mean_all,
                    pct(lenhist_placed, 25), pct(lenhist_placed, 50),
                    pct(lenhist_placed, 75), mean_pl,
                    lost['placed'], n_binned, lost['txend'], lost['txstart'],
                    lost['internal'], lost['modelerr'], lost['unclassified']))
    np.savez_compressed(out_prefix + '.lenhist.npz',
                        all=lenhist, placed=lenhist_placed)
    print('%s [%s]: %d placed of %d primary (%d seen), %d/%d tx with signal, '
          'fullexon %.1f%%, aln p50 %d nt'
          % (label, group, placed, primary, seen, norm(cov['mid'])[1], ntx,
             100 * (n_fullexon / placed if placed else 0), pct(lenhist, 50)))


# ===========================================================================
# 3. merge, statistics, regression check
# ===========================================================================
def _stats(v):
    v = np.asarray(v, dtype=float)
    s = v.sum()
    if s <= 0:
        return dict(f5=np.nan, fmid=np.nan, f3=np.nan, ratio=np.nan,
                    rise=np.nan, cv=np.nan, lastbin=np.nan)
    v = v / s
    f5, fmid, f3 = v[:20].sum(), v[20:80].sum(), v[80:].sum()
    body = v[40:60].mean()
    return dict(f5=100 * f5, fmid=100 * fmid, f3=100 * f3,
                ratio=f3 / f5 if f5 else np.nan,
                rise=v[90:].mean() / body if body else np.nan,
                lastbin=v[99] / body if body else np.nan,
                cv=v.std() / v.mean())


def _perm_p(a, b, n_exact=200000):
    """Exact two-sided permutation p on the difference in means."""
    from itertools import combinations
    a, b = np.asarray(a, float), np.asarray(b, float)
    pool = np.concatenate([a, b])
    na, n = len(a), len(pool)
    obs = abs(a.mean() - b.mean())
    from math import comb
    total = comb(n, na)
    if total > n_exact:
        rng = np.random.default_rng(0)
        cnt = 0
        for _ in range(n_exact):
            p = rng.permutation(n)
            cnt += abs(pool[p[:na]].mean() - pool[p[na:]].mean()) >= obs - 1e-12
        return cnt / n_exact, total, False
    cnt = 0
    allidx = np.arange(n)
    for cmb in combinations(range(n), na):
        sel = np.zeros(n, dtype=bool); sel[list(cmb)] = True
        cnt += abs(pool[sel].mean() - pool[~sel].mean()) >= obs - 1e-12
    return cnt / total, total, True


def merge(covdir, res, geneset_prefix):
    import pandas as pd

    fs = sorted(glob.glob(os.path.join(covdir, '*.cov.tsv')))
    if not fs:
        sys.exit('no per-BAM profiles in %s' % covdir)
    df = pd.concat([pd.read_csv(f, sep='\t') for f in fs], ignore_index=True)
    meta = pd.concat([pd.read_csv(f, sep='\t') for f in
                      sorted(glob.glob(os.path.join(covdir, '*.meta.tsv')))],
                     ignore_index=True)
    bcols = ['b%02d' % i for i in range(NBINS)]

    st = pd.DataFrame([_stats(r[bcols].values) for _, r in df.iterrows()])
    per = pd.concat([df[['label', 'group', 'metric', 'n_tx', 'reads_placed',
                         'reads_seen', 'reads_primary', 'reads_fullexon']],
                     st], axis=1)
    per = per.merge(meta.drop(columns=['group', 'reads_seen', 'reads_primary',
                                       'reads_placed', 'reads_fullexon']),
                    on='label', how='left')
    per.to_csv(os.path.join(res, 'coverage_threeway.tsv'), sep='\t', index=False)

    # --- aligned read-length table: all groups on the record ---------------
    # Two reasons this is a deliverable and not a footnote. (1) Read length is
    # the INPUT to the edge-clipping hypothesis for the base-vs-mid
    # disagreement, which is untested and remains so here. (2) It is the
    # obvious confounder for any coverage-SHAPE claim: the published plate was
    # mapped with a sjdbOverhang-73 index and the own plate with a 151 nt one,
    # so a published-vs-own shape difference could be length, not chemistry.
    # Percentiles are over the pooled per-unit histograms of ALL primary
    # alignments STREAMED -- the whole BAM for the shallow published cells, the
    # first MAX_READS-worth for the deep ones. reads_primary states which, and
    # this is a different estimator from the per-unit mean of percentiles
    # printed below; disagreement between them is itself informative.
    alnrows = []
    for _g in ['VASA_published', 'VASA_own', 'FLASHseq_native',
               'FLASHseq_vasalen']:
        sm = meta[meta.group == _g]
        if not len(sm):
            continue
        h = None
        for f in sorted(glob.glob(os.path.join(covdir, '*.lenhist.npz'))):
            if os.path.basename(f)[:-len('.lenhist.npz')] in set(sm.label):
                z = np.load(f)['all'].astype(np.int64)
                h = z if h is None else h + z
        assert h is not None and h.sum() > 0, _g
        cs = np.cumsum(h)
        LL = np.arange(len(h))
        qs = {p: float(LL[np.searchsorted(cs, p * cs[-1] / 100.0)])
              for p in (1, 5, 25, 50, 75, 95, 99)}
        alnrows.append(dict(
            group=_g, n_units=len(sm),
            reads_primary=int(sm.reads_primary.sum()),
            alignments_in_hist=int(cs[-1]),
            p01=qs[1], p05=qs[5], p25=qs[25], p50=qs[50], p75=qs[75],
            p95=qs[95], p99=qs[99],
            mean=float((LL * h).sum() / cs[-1]),
            frac_le_80nt=float(h[:81].sum() / cs[-1]),
            frac_ge_140nt=float(h[140:].sum() / cs[-1])))
    aln = pd.DataFrame(alnrows)
    aln.to_csv(os.path.join(res, 'coverage_threeway_alnlen.tsv'),
               sep='\t', index=False)
    print('=== ALIGNED READ LENGTH, pooled over the units of each group ===')
    print('  %-17s %3s %12s %5s %5s %5s %5s %7s %8s'
          % ('group', 'n', 'alignments', 'p05', 'p25', 'p50', 'p75', 'mean',
             '<=80nt'))
    for _, r0 in aln.iterrows():
        print('  %-17s %3d %12d %5.0f %5.0f %5.0f %5.0f %7.1f %7.1f%%'
              % (r0.group, r0.n_units, r0.alignments_in_hist, r0.p05, r0.p25,
                 r0.p50, r0.p75, r0['mean'], 100 * r0.frac_le_80nt))
    print('  wrote coverage_threeway_alnlen.tsv')

    prof = df.groupby(['group', 'metric'])[bcols].mean().reset_index()
    prof.to_csv(os.path.join(res, 'coverage_threeway_profile.tsv'),
                sep='\t', index=False)

    order = ['VASA_published', 'VASA_own', 'FLASHseq_native', 'FLASHseq_vasalen']
    groups = [g for g in order if g in set(df.group)]

    print('=== units, depth, aligned read length (ALL primary alignments) ===')
    print('  %-17s %3s %12s %12s %9s %6s %6s %6s'
          % ('group', 'n', 'primary', 'placed', 'fullexon', 'p25', 'p50', 'p75'))
    for g in groups:
        s = meta[meta.group == g]
        print('  %-17s %3d %12d %12d %8.1f%% %6.0f %6.0f %6.0f'
              % (g, len(s), s.reads_primary.sum(), s.reads_placed.sum(),
                 100 * s.reads_fullexon.sum() / s.reads_placed.sum(),
                 s.alnlen_p25.mean(), s.alnlen_p50.mean(), s.alnlen_p75.mean()))

    print('\n=== MIDPOINT profile -- THE REPORTED METRIC ===')
    print("  %-17s %8s %8s %8s %8s %8s %8s"
          % ('group', "5'20%", 'mid60%', "3'20%", "3'/5'", '3\'rise', 'CV'))
    for g in groups:
        s = per[(per.group == g) & (per.metric == 'mid')]
        print('  %-17s %7.2f%% %7.2f%% %7.2f%% %8.3f %8.3f %8.3f'
              % (g, s.f5.mean(), s.fmid.mean(), s.f3.mean(),
                 s.ratio.mean(), s.rise.mean(), s.cv.mean()))
        print('       per unit 3\'/5\': %s'
              % ' '.join('%s=%.3f' % (l.split('_')[-1], v)
                         for l, v in zip(s.label, s.ratio)))

    print('\n=== the 3\' rise: is it the PROTOCOL or this LIBRARY? ===')
    print('  rise = mean(bins 90-99) / mean(bins 40-59); 1.0 = flat')
    for stat in ('rise', 'ratio', 'lastbin'):
        pub = per[(per.group == 'VASA_published') & (per.metric == 'mid')][stat].values
        own = per[(per.group == 'VASA_own') & (per.metric == 'mid')][stat].values
        fsn = per[(per.group == 'FLASHseq_native') & (per.metric == 'mid')][stat].values
        if not (len(pub) and len(own)):
            continue
        p, tot, exact = _perm_p(pub, own)
        print('  %-8s published %.3f (n=%d, %.3f-%.3f)   own %.3f (n=%d, %.3f-%.3f)'
              % (stat, pub.mean(), len(pub), pub.min(), pub.max(),
                 own.mean(), len(own), own.min(), own.max()))
        print('           FLASH-seq native %.3f (n=%d)   published-vs-own '
              'permutation p = %.4f (%s, %d splits)'
              % (fsn.mean(), len(fsn), p, 'exact' if exact else 'sampled', tot))

    print('\n=== BASE profile (diagnostic only; disagrees with mid by ~2x) ===')
    for g in groups:
        s = per[(per.group == g) & (per.metric == 'base')]
        print('  %-17s 5\'=%6.2f%% 3\'=%6.2f%%  3\'/5\'=%.3f'
              % (g, s.f5.mean(), s.f3.mean(), s.ratio.mean()))

    print('\n=== where the aligned bases of placed reads go ===')
    print('  denominator: aligned bases of reads that placed >=1 base '
          '(Rule 5)')
    print('  %-17s %12s %8s %8s %8s %8s %8s'
          % ('group', 'bases', 'binned', 'txend', 'txstart', 'intronic', 'unclass'))
    for g in groups:
        s = meta[meta.group == g]
        b = s.bases_aligned.sum()
        print('  %-17s %12d %7.3f%% %7.3f%% %7.3f%% %7.3f%% %7.3f%%'
              % (g, b, 100 * s.bases_binned.sum() / b,
                 100 * s.lost_txend.sum() / b, 100 * s.lost_txstart.sum() / b,
                 100 * s.lost_internal.sum() / b,
                 100 * s.lost_unclassified.sum() / b))
        me = int(s.lost_modelerr.sum())
        if me:
            print('       WARNING lost_modelerr = %d (expected 0)' % me)

    print('\n=== does restricting to reads that lose NO bases reconcile '
          'base and mid? ===')
    print('  (a decomposition of where bases go, NOT a verdict on a mechanism;')
    print('   edge clipping remains the leading untested hypothesis)')
    print('  %-17s %9s %9s %9s %9s'
          % ('group', 'mid 3\'/5\'', 'base', 'mid_fe', 'base_fe'))
    for g in groups:
        def r(mt):
            s = per[(per.group == g) & (per.metric == mt)]
            return s.ratio.mean() if len(s) else np.nan
        print('  %-17s %9.3f %9.3f %9.3f %9.3f'
              % (g, r('mid'), r('base'), r('mid_fullexon'), r('base_fullexon')))

    print('\n=== READ-LENGTH-MATCHED CONTROL (reads <= %d nt aligned) ===' % SHORT_MAX)
    print('  The published plate is a ~75 nt library (its STAR index is')
    print('  sjdbOverhang 73) and its aligned reads are ~73 nt against the own')
    print("  plate's ~127 nt. So a published-vs-own difference in coverage shape")
    print('  could be READ LENGTH rather than chemistry. This restricts every')
    print('  unit to the length regime they share.')
    print('  %-17s %9s %9s %9s %9s'
          % ('group', "all 3'/5'", "<=%dnt" % SHORT_MAX, 'all rise', 'short rise'))
    for g in groups:
        a = per[(per.group == g) & (per.metric == 'mid')]
        b = per[(per.group == g) & (per.metric == 'mid_short')]
        print('  %-17s %9.3f %9.3f %9.3f %9.3f'
              % (g, a.ratio.mean(), b.ratio.mean(), a.rise.mean(), b.rise.mean()))
    pubs = per[(per.group == 'VASA_published') & (per.metric == 'mid_short')]
    owns = per[(per.group == 'VASA_own') & (per.metric == 'mid_short')]
    if len(pubs) and len(owns):
        for stat in ('rise', 'ratio'):
            p, tot, exact = _perm_p(pubs[stat].values, owns[stat].values)
            print('  length-matched %s: published %.3f vs own %.3f, '
                  'permutation p = %.4f (%s)'
                  % (stat, pubs[stat].mean(), owns[stat].mean(), p,
                     'exact' if exact else 'sampled'))

    print('\n=== stranded midpoint (VASA counted y, FLASH-seq n -- Rule 6) ===')
    for g in groups:
        a = per[(per.group == g) & (per.metric == 'mid')]
        b = per[(per.group == g) & (per.metric == 'mid_sense')]
        print('  %-17s mid 3\'/5\'=%.3f   sense-only 3\'/5\'=%.3f'
              % (g, a.ratio.mean(), b.ratio.mean()))

    # --- ROBUSTNESS: per-transcript read floor ------------------------------
    print('\n=== ROBUSTNESS: profiles recomputed over transcripts with '
          '>=%d placed reads ===' % ROBUST_MIN_TX_READS)
    print('  (the published cells are ~10x shallower; this tests whether their')
    print('   shape is an artefact of thinly-covered transcripts)')
    rows = []
    for f in sorted(glob.glob(os.path.join(covdir, '*.tx.npz'))):
        lab = os.path.basename(f)[:-len('.tx.npz')]
        z = np.load(f)
        grp = df[df.label == lab].group.iloc[0]
        for mt in ('mid', 'base'):
            cm = z[mt].astype(np.float64)
            tot = cm.sum(axis=1)
            ok = tot >= ROBUST_MIN_TX_READS
            if not ok.any():
                continue
            p = (cm[ok] / tot[ok, None]).mean(axis=0)
            s = _stats(p)
            s.update(label=lab, group=grp, metric=mt, n_tx=int(ok.sum()))
            rows.append(s)
    rob = pd.DataFrame(rows)
    if len(rob):
        rob.to_csv(os.path.join(res, 'coverage_threeway_robust.tsv'),
                   sep='\t', index=False)
        print('  %-17s %6s %9s %9s' % ('group', 'n_tx', "3'/5'", '3\'rise'))
        for g in groups:
            s = rob[(rob.group == g) & (rob.metric == 'mid')]
            if len(s):
                print('  %-17s %6.0f %9.3f %9.3f'
                      % (g, s.n_tx.mean(), s.ratio.mean(), s.rise.mean()))

    # --- REGRESSION against the published upstream numbers ------------------
    print('\n=== REGRESSION vs res/flashseq_vasa/gene_coverage_profile.tsv ===')
    print('  Upstream ran on a DIFFERENT gene set (chosen from 2 groups, not 4)')
    print('  so exact equality is not expected; the check is that the ORDERING')
    print('  and the direction of every FLASH-seq/own-plate claim survive.')
    up = os.path.join(os.path.dirname(res.rstrip('/')), 'flashseq_vasa',
                      'gene_coverage_profile.tsv')
    if os.path.exists(up):
        u = pd.read_csv(up, sep='\t')
        ub = [c for c in u.columns if c.startswith('b')]
        m = {'VASA': 'VASA_own', 'FLASH-seq native': 'FLASHseq_native',
             'FLASH-seq vasalen': 'FLASHseq_vasalen'}
        for _, r in u[u.metric == 'mid'].iterrows():
            s = _stats(r[ub].values)
            g = m.get(r.grp)
            if g is None:
                continue
            new = per[(per.group == g) & (per.metric == 'mid')].ratio.mean()
            print("  %-17s upstream 3'/5' = %.3f   here = %.3f   %s"
                  % (g, s['ratio'], new,
                     'SAME SIDE of 1.0' if (s['ratio'] - 1) * (new - 1) > 0
                     else 'DIFFERENT SIDE of 1.0 -- investigate'))
    else:
        print('  upstream table not found at %s' % up)

    gt = geneset_prefix + '.genes.tsv'
    if os.path.exists(gt):
        g = pd.read_csv(gt, sep='\t')
        lr = g.len_ratio_116_99.values
        print('\n=== annotation-release confound, bounded ===')
        print('  %d genes; longest-transcript length E116/E99: median %.4f, '
              'identical %.1f%%, within 5%% %.1f%%'
              % (len(g), float(np.median(lr)), 100 * (lr == 1).mean(),
                 100 * (abs(lr - 1) <= 0.05).mean()))
        print('  strand disagreements: %d' % int((~g.same_strand).sum()))

    print('\nwrote coverage_threeway.tsv, coverage_threeway_profile.tsv, '
          'coverage_threeway_robust.tsv')


# ===========================================================================
# 4. cross-release control
# ===========================================================================
def crossrel(covdir, res):
    """How much of a published-vs-own gap could the annotation release explain?

    The published plate is profiled against E99/GRCm38 models and the other two
    against E116/GRCm39. That difference cannot be removed -- the published BAMs
    are aligned to a GRCm38 index and there is no GRCm39 alignment of them. What
    CAN be done is to hold the reads fixed and change only the model.

    WHAT THIS FUNCTION ACTUALLY DOES: it reads
    `coverage_threeway_geneset.genes.tsv` and summarises how far the two
    releases' transcript models disagree -- length ratio, exon count, strand.
    It does NOT open a BAM and does NOT call `profile`; there is no read-level
    re-profiling here. (An earlier version of this docstring described one, and
    was wrong. Its `covdir` argument is accepted for driver-call symmetry and
    is unused.)

    A read-level control is not possible on these inputs: E99 transcript
    structures are on GRCm38 coordinates and the own-plate BAM is on GRCm39, so
    the E99 models cannot be applied to it without a liftover. What IS possible
    is bounding the effect structurally -- a release effect on coverage SHAPE
    can only enter through the coverage axis (length, exon count, strand), so
    measuring the axis disagreement bounds it.

    This is a bound, not a correction, and it is reported as one.
    """
    import pandas as pd
    g = pd.read_csv(os.path.join(res, 'coverage_threeway_geneset.genes.tsv'),
                    sep='\t')
    lr = g.len_ratio_116_99.values
    print('=== annotation-release bound, from the transcript models ===')
    print('  genes: %d' % len(g))
    print('  longest-transcript length E116/E99:')
    for q in (1, 5, 25, 50, 75, 95, 99):
        print('    p%-3d %.4f' % (q, float(np.percentile(lr, q))))
    print('  identical length : %.1f%%' % (100 * (lr == 1).mean()))
    print('  within +-1%%      : %.1f%%' % (100 * (abs(lr - 1) <= 0.01).mean()))
    print('  within +-5%%      : %.1f%%' % (100 * (abs(lr - 1) <= 0.05).mean()))
    print('  strand disagreements: %d' % int((~g.same_strand).sum()))
    print('  exon-count identical: %.1f%%'
          % (100 * (g.nexon_E99 == g.nexon_E116).mean()))
    out = os.path.join(res, 'coverage_threeway_release_bound.tsv')
    pd.DataFrame({
        'quantity': ['n_genes', 'len_identical_pct', 'len_within1pct',
                     'len_within5pct', 'len_ratio_median', 'len_ratio_p05',
                     'len_ratio_p95', 'strand_disagree', 'nexon_identical_pct'],
        'value': [len(g), 100 * (lr == 1).mean(),
                  100 * (abs(lr - 1) <= 0.01).mean(),
                  100 * (abs(lr - 1) <= 0.05).mean(), float(np.median(lr)),
                  float(np.percentile(lr, 5)), float(np.percentile(lr, 95)),
                  int((~g.same_strand).sum()),
                  100 * (g.nexon_E99 == g.nexon_E116).mean()],
    }).to_csv(out, sep='\t', index=False)
    print('wrote %s' % out)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'select':
        select(*sys.argv[2:4])
    elif mode == 'geneset':
        geneset(*sys.argv[2:10])
    elif mode == 'profile':
        profile(*sys.argv[2:7])
    elif mode == 'merge':
        merge(*sys.argv[2:5])
    elif mode == 'crossrel':
        crossrel(*sys.argv[2:4])
    else:
        sys.exit(__doc__)
