#!/usr/bin/env python3
"""
step7_precheck.py -- verify step 6's pickle, and prove step 7 can run on it.

WHY THIS EXISTS
---------------
Step 7 (`../a_Mapping/countTables_fromPickle.py`) is a long single-threaded
pandas job over the whole 16-cell structure, and it has already failed once in
this repo for a reason that was invisible beforehand: on the VASA-plate mixed
run (job 50542441) `reduceGeneName` raised IndexError because the annotation
BED carried gene names that did not satisfy the pipeline's 4-field naming
contract. That cost the whole run. This script replays every string operation
step 7 performs on an index label, over every label, before any of it is
submitted.

It is READ-ONLY. It writes one report file and touches nothing else.

It answers four questions:

  1. Is the step-6 pickle complete and fresh? Cells in the pickle == cells on
     disk, and each cell's step-5 BEDs are older than the pickle. ("Counting
     files is not proof" -- see own_version/README.md, the step-4 staleness
     incident.)
  2. Did the v2 BED actually deliver? Under the v1 BED every tRNA table step 7
     writes was empty. tRNA rows must now be present in the pickle itself.
  3. Would step 7 crash, or silently produce a garbage row? Every index label
     is pushed through `reduceGeneName` exactly as step 7 calls it, and both
     failure modes are counted: an exception, and an empty-string result.
  4. How big is it, and what does the unigene filter do at n=16 cells? Step 7's
     `ncells = max(5, round(0.01*ncols))` was written for a 384-cell plate;
     at 16 cells it demands a gene be seen in 5/16 = 31% of cells rather than
     5/384 = 1.3%. That is a dataset-size sensitivity, not a bug, but it
     changes which combination names collapse, so it is measured here.

Usage:
    step7_precheck.py <pickle.gz> <celldir> <outdir> [report.txt]
"""
import sys
import os
import glob
import gzip
import time
import pickle
import resource
from collections import Counter

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# these definitions are COPIED VERBATIM from
# ../a_Mapping/countTables_fromPickle.py so that what is tested here is exactly
# what step 7 will execute. Do not "clean them up" -- the point is fidelity.
# ---------------------------------------------------------------------------
def reduceGeneName(gene, uni_genes):
    rg = gene
    if gene.count('-') == 0:
        rg = gene
    else:
        bios = set([x.rsplit('_')[-1] for x in gene.rsplit('-')])
        shortlived = ['miRNA', 'tRNA', 'MtTrna']
        longstuff = ['lncRNA']
        shortstuff = ['snRNA', 'snoRNA', 'MiscRna', 'scaRNA']
        ribos = ['rRNA', 'ribozyme']
        if any([b in ribos for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in ribos])
            rg = gene
        if any([b not in shortlived for b in bios]) and any([b in shortlived for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] not in shortlived])
            rg = gene
        if any([b in shortstuff for b in bios]) and any([b not in shortstuff for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in shortstuff])
            rg = gene
        if sum([g in uni_genes for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g in uni_genes][0]
            gene = rg
        if gene.count('-') >= 1 and sum([g.rsplit('_')[1][:2] != "Gm" for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g.rsplit('_')[1][:2] != "Gm"][0]
    return rg


def countTotalReads(x, protocol='vasa'):
    if protocol in ['vasa', '10x', 'smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x]) if type(x) == dict else 0
    return y


def countTotalUMI(x):
    return len(x) if type(x) == dict else 0


def remove_ENSandGm(gene):
    rg = sorted(set(['_'.join(x.rsplit("_")[1:]) for x in gene.rsplit('-')]))
    xg = [g for g in rg if g[:2] != 'Gm']
    if len(xg) == 0:
        xg = rg
    xg = '-'.join(xg)
    return xg


# ---------------------------------------------------------------------------
def rss_gb():
    # ru_maxrss is KiB on Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    picklegz, celldir, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    reportpath = sys.argv[4] if len(sys.argv) > 4 else os.path.join(outdir, 'logs', 'step7_precheck.txt')

    out = []
    def emit(s=''):
        print(s, flush=True)
        out.append(s)

    def section(t):
        emit()
        emit('=' * 78)
        emit(t)
        emit('=' * 78)

    fail = []
    warn = []

    section('0. INPUTS')
    emit('pickle : %s' % picklegz)
    emit('celldir: %s' % celldir)
    for p in (picklegz, celldir):
        if not os.path.exists(p):
            sys.exit('MISSING: %s' % p)
    pkl_mtime = os.path.getmtime(picklegz)
    emit('pickle mtime : %s' % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(pkl_mtime)))
    emit('pickle size  : %.1f MB' % (os.path.getsize(picklegz) / 1e6))

    # -----------------------------------------------------------------------
    section('1. FRESHNESS -- is the pickle newer than every input it consumed?')
    # step 6 takes its cell list from the singlemapper BEDs alone; a cell whose
    # single file is missing disappears from the pickle silently.
    singles = sorted(glob.glob(os.path.join(celldir, '*.singlemappers_genes.bed.gz')))
    multis = sorted(glob.glob(os.path.join(celldir, '*multimappers_genes.bed.gz')))
    emit('step5 singlemapper BEDs on disk: %d' % len(singles))
    emit('step5 multimapper  BEDs on disk: %d' % len(multis))
    stale = [os.path.basename(f) for f in singles + multis if os.path.getmtime(f) > pkl_mtime]
    if stale:
        fail.append('%d step-5 BED(s) are NEWER than the pickle -- step 6 output is stale' % len(stale))
        for s in stale[:10]:
            emit('  STALE: %s' % s)
    else:
        emit('OK: all %d step-5 BEDs predate the pickle' % (len(singles) + len(multis)))
    # a leftover uncompressed file is the tell for the bare-gzip re-run hazard
    leftovers = glob.glob(os.path.join(celldir, '*_genes.bed')) + \
                [f for f in glob.glob(os.path.join(outdir, '*.pickle')) if not f.endswith('dict.pickle')]
    if leftovers:
        warn.append('uncompressed leftovers beside the outputs: %s'
                    % ', '.join(os.path.basename(f) for f in leftovers))

    # -----------------------------------------------------------------------
    section('2. LOAD -- this is exactly what step 7 does first')
    t0 = time.time()
    cntdf = pickle.load(gzip.open(picklegz, 'rb'))
    cntdf = cntdf[sorted(cntdf.columns)]
    t_load = time.time() - t0
    emit('load time     : %.1f s' % t_load)
    emit('peak RSS      : %.2f GB' % rss_gb())
    emit('cntdf.shape   : %s   (rows = gene entries, cols = cells)' % (cntdf.shape,))
    emit('columns       : %s' % ', '.join(map(str, cntdf.columns)))

    disk_cells = set(os.path.basename(f).split('_cbc')[0] for f in singles)
    emit()
    emit('cells on disk (from BED filenames): %d' % len(disk_cells))
    emit('cells in pickle                   : %d' % len(cntdf.columns))
    if len(cntdf.columns) != len(disk_cells):
        fail.append('cell count mismatch: %d in pickle vs %d on disk'
                    % (len(cntdf.columns), len(disk_cells)))

    # CELLID_FROM=r should give bare '001'..'016'; =f gives 'cells/SAMPLE_001'
    if any('/' in str(cc) for cc in cntdf.columns):
        warn.append('column names contain a path separator -- CELLID_FROM was f, not r')

    # -----------------------------------------------------------------------
    section('3. tRNA -- did the v2 BED actually deliver?')
    t0 = time.time()
    idx = list(cntdf.index)
    tRNAs = [i for i in idx if 'tRNA' in i]
    tset = set(tRNAs)
    genes = [i for i in idx if i not in tset]
    emit('tRNA index labels   : %d' % len(tRNAs))
    emit('non-tRNA labels     : %d' % len(genes))
    emit('(label split took %.1f s)' % (time.time() - t0))
    if len(tRNAs) == 0:
        fail.append('ZERO tRNA labels in the pickle -- the v1 BED (no cytoplasmic tRNA) '
                    'was used, or REF_BED did not take effect')
    else:
        emit('examples: %s' % ', '.join(tRNAs[:5]))
        iso = sorted(set('-'.join(sorted(set([t.rsplit('.')[-1] for t in i.rsplit('-')]))) for i in tRNAs))
        pure = [i for i in iso if '-' not in i]
        emit('isotype classes after collapse: %d total, %d pure (single-locus)' % (len(iso), len(pure)))
        emit('pure isotypes: %s' % ', '.join(pure[:70]))

    mt = [i for i in genes if 'MtTrna' in i]
    emit('MtTrna labels landing in the GENE tables (case-sensitivity, expected): %d' % len(mt))

    # -----------------------------------------------------------------------
    section('4. NAMING CONTRACT -- would reduceGeneName crash or empty a name?')
    uni_genes = [g for g in genes if '-' not in g]
    emit('uni (single-gene) labels : %d' % len(uni_genes))
    emit('combination labels       : %d' % (len(genes) - len(uni_genes)))

    # the exact IndexError from run 1: g.rsplit('_')[1] on a component with no '_'
    bad_fields = []
    for g in genes:
        for comp in g.rsplit('-'):
            if len(comp.rsplit('_')) < 2:
                bad_fields.append((g, comp))
                break
    if bad_fields:
        fail.append('%d label(s) have a component with <2 underscore fields -- '
                    'reduceGeneName will raise IndexError (this is how run 1 died)'
                    % len(bad_fields))
        for g, comp in bad_fields[:10]:
            emit('  BAD: %-60s component=%r' % (g, comp))
    else:
        emit('OK: every component of every label has >=2 underscore fields')

    # -----------------------------------------------------------------------
    section('5. UNIGENE FILTER -- what ncells=max(5, round(0.01*n)) does here')
    t0 = time.time()
    total_reads_df = cntdf.apply(lambda col: col.map(lambda x: countTotalReads(x, 'vasa')))
    emit('reads-per-cell pass: %.1f s, peak RSS %.2f GB' % (time.time() - t0, rss_gb()))
    ncells = max(5, round(0.01 * len(cntdf.columns)))
    nreads = 1
    emit()
    emit('n cells = %d  ->  ncells threshold = %d  (%.1f%% of cells)'
         % (len(cntdf.columns), ncells, 100.0 * ncells / len(cntdf.columns)))
    emit('  for reference, the published 384-cell plate: threshold 5 = 1.3% of cells')
    det = (total_reads_df.loc[uni_genes] >= nreads).sum(axis=1)
    uni_genes_filt = np.array(uni_genes)[det >= ncells]
    emit('unigenes passing (>=1 read in >=%d cells): %d / %d  (%.1f%%)'
         % (ncells, len(uni_genes_filt), len(uni_genes),
            100.0 * len(uni_genes_filt) / max(1, len(uni_genes))))
    for thr in (1, 2, 3, 5, 8, 12, 16):
        emit('    >= %2d cells : %7d unigenes' % (thr, int((det >= thr).sum())))

    # -----------------------------------------------------------------------
    section('6. reduceGeneName REPLAY over every label')
    t0 = time.time()
    errs = []
    empties = []
    newnames = []
    for g in genes:
        try:
            r = reduceGeneName(g, uni_genes_filt)
        except Exception as e:
            errs.append((g, '%s: %s' % (type(e).__name__, e)))
            newnames.append(g)
            continue
        if r == '':
            empties.append(g)
        newnames.append(r)
    emit('replay time: %.1f s' % (time.time() - t0))
    if errs:
        fail.append('reduceGeneName RAISED on %d label(s) -- step 7 would abort' % len(errs))
        for g, e in errs[:10]:
            emit('  RAISE: %-60s %s' % (g, e))
    else:
        emit('OK: no exception on any of %d labels' % len(genes))

    if empties:
        warn.append('reduceGeneName returns an EMPTY name for %d label(s); step 7 will '
                    'merge them all into one row named "" (upstream behaviour, silent)'
                    % len(empties))
        emit()
        emit('EMPTY-NAME labels (first 10):')
        for g in empties[:10]:
            emit('  %s' % g)
        er = int(total_reads_df.loc[empties].sum().sum())
        gr = int(total_reads_df.loc[genes].sum().sum())
        emit('reads carried by those labels: %d' % er)
        emit('  (as a fraction of all gene reads: %.4f%%)' % (100.0 * er / max(1, gr)))
    else:
        emit('OK: no label reduces to an empty name')

    emit()
    emit('distinct names after reduction: %d  (from %d labels)' % (len(set(newnames)), len(genes)))

    # -----------------------------------------------------------------------
    section('7. SCALE -- per-cell totals and the STAR cross-check')
    t0 = time.time()
    total_umi_df = cntdf.apply(lambda col: col.map(countTotalUMI))
    emit('umi pass: %.1f s, peak RSS %.2f GB' % (time.time() - t0, rss_gb()))
    star = {}
    for f in glob.glob(os.path.join(celldir, '*_E99_Log.final.txt')):
        b = os.path.basename(f)
        cell = b.split('_cbc')[0].rsplit('_', 1)[-1]
        with open(f) as fh:
            for line in fh:
                if 'Uniquely mapped reads number' in line:
                    star.setdefault(cell, {})['uniq'] = int(line.strip().rsplit('\t')[-1])
                if 'Number of input reads' in line:
                    star.setdefault(cell, {})['input'] = int(line.strip().rsplit('\t')[-1])
    emit()
    emit('%-6s %14s %14s %10s %12s %8s' % ('cell', 'reads', 'molecules', 'entries', 'STAR input', 'assign%'))
    for cc in cntdf.columns:
        s = star.get(str(cc), {})
        rd = int(total_reads_df[cc].sum())
        si = s.get('input')
        pct = ('%.1f' % (100.0 * rd / si)) if si else 'n/a'
        emit('%-6s %14d %14d %10d %12s %8s'
             % (cc, rd, int(total_umi_df[cc].sum()),
                int((total_umi_df[cc] > 0).sum()), si if si else 'n/a', pct))
    emit()
    emit('TOTAL  %14d %14d' % (int(total_reads_df.values.sum()), int(total_umi_df.values.sum())))

    # -----------------------------------------------------------------------
    section('8. UMI SATURATION -- 6 nt UMI means K=4096 per gene')
    umi0 = sorted([x for x in cntdf[cntdf.columns[0]] if type(x) == dict][0].keys())[0]
    K = 4 ** len(umi0)
    emit('inferred UMI length: %d  ->  K = %d' % (len(umi0), K))
    for cc in cntdf.columns:
        col = total_umi_df[cc]
        emit('  %-5s  >1000 UMIs: %5d   >2048: %5d   >=K: %4d   share of molecules in >1000: %5.1f%%'
             % (cc, int((col > 1000).sum()), int((col > 2048).sum()), int((col >= K).sum()),
                100.0 * col[col > 1000].sum() / max(1, col.sum())))

    # -----------------------------------------------------------------------
    section('9. SIZING for the step-7 sbatch')
    emit('peak RSS of this precheck (load + 2 full passes): %.2f GB' % rss_gb())
    emit('step 7 does ~12 such passes plus a groupby-aggregate over object dtype,')
    emit('and fixGeneLabels rewrites Counters in place for every reduced label.')
    emit('Request at least 3x this precheck peak, and give it a generous walltime.')

    # -----------------------------------------------------------------------
    section('VERDICT')
    if fail:
        emit('FAIL (%d):' % len(fail))
        for f in fail:
            emit('  x %s' % f)
    else:
        emit('No blocking problem found. Step 7 can be submitted.')
    if warn:
        emit()
        emit('WARN (%d) -- not blocking, but know about them:' % len(warn))
        for w in warn:
            emit('  ! %s' % w)

    os.makedirs(os.path.dirname(reportpath), exist_ok=True)
    with open(reportpath, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % reportpath, flush=True)
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
