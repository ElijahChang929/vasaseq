#!/usr/bin/env python3
"""
noumi_precheck.py -- prove VASA's protocol='smartseq_noUMI' branch can run on
FLASH-seq data, BEFORE any of it is submitted.

WHY THIS EXISTS
---------------
The `smartseq_noUMI` branch exists in upstream `a_Mapping/` but has never been
exercised in this repo. Step 7 has already destroyed one run here for a reason
that was invisible beforehand (`reduceGeneName` IndexError on gene names that
broke the 4-field naming contract, job 50542441), and the response was
`own_version/step7_precheck.py`. This script is the same idea aimed at the
no-UMI branch: it replays the branch's own string operations over real
FLASH-seq read names and over the real v2 annotation BED, and reports what
would break.

It is READ-ONLY. It writes one report file and touches nothing else.

It answers six questions:

  1. READ NAMES. `get_UMI(..., 'smartseq_noUMI')` builds a dict out of the
     read name before throwing the result away. Do FLASH-seq read names
     survive that? In which of the three forms they can reach it (raw FASTQ
     header, STAR-truncated, SM-tagged)?
  2. CELL ID. Step 6 gets the column name either from the BED filename
     (`cellidFROMfilename='f'`) or from the read name's SM tag (`='r'`).
     Replay both over real strings and report the literal column name each
     would produce.
  3. NAMING CONTRACT. Every gene name in the v2 BED, after the transforms
     step 6 applies to it, must satisfy the >=2-underscore-field contract or
     `reduceGeneName` raises -- and must not contain a '-', which is the
     separator step 6 uses to build combination names.
  4. NO-UMI SEMANTICS. What the counts MEAN once every read is filed under the
     single literal UMI 'A'. Demonstrated numerically, not asserted.
  5. STEP-5 BED CONTRACT. Can step 6 `int()` the nM field of every row of a
     real step-5 BED? This check was ADDED after the first dry run, where it
     is exactly what broke: on paired-end input `bedtools bamtobed` appends the
     '/1','/2' mate suffix, which lands on the end of the nM value, and step 6
     dies with `ValueError: invalid literal for int() with base 10: '0/2'`.
     Give it the BEDs and it reports the rate rather than discovering it in
     an hour-long job.
  6. FILTER SENSITIVITY. `ncells = max(5, round(0.01*ncols))` at the library
     counts this comparison will actually use (1 and 10), not 384.

PROVENANCE OF THE HELPERS BELOW
-------------------------------
The point is fidelity, not tidiness. What is tested here must be what will
execute.

  VERBATIM, character-for-character identical to the upstream source
  (`a_Mapping/countTables_2pickle_cellsSpliced.py`, `countTables_fromPickle.py`
  at the commit this repo has):
    - get_UMI                (2pickle, lines 36-45)
    - addCount               (2pickle, lines 26-34)
    - countTotalReads        (fromPickle, lines 47-54)  -- BOTH branches kept
    - countTotalUMI          (fromPickle, lines 55-56)
    - countExonReads         (fromPickle, lines 59-68)  -- BOTH branches kept
    - countIntronReads       (fromPickle, lines 73-82)  -- BOTH branches kept
    - countExonUMI           (fromPickle, lines 69-71)
    - countIntronUMI         (fromPickle, lines 83-85)
    - reduceGeneName         (fromPickle, lines 142-165)
    - remove_ENSandGm        (fromPickle, lines 265-271)

  TRIMMED, and why:
    - bc2trans (fromPickle, lines 89-96) is reproduced as `bc2trans_K(x, K)`
      with K passed in rather than read off a module-level global. Upstream
      computes K once at line 88 from the first UMI it finds in the pickle
      (line 87); this script needs to evaluate it at several K to show what
      the no-UMI case does, so K became an argument. The arithmetic in the
      three branches is unchanged.
    - `cellid` extraction (2pickle, lines 89-92) is lifted out of
      `get_cellDict` into `cellid_from_filename` / `cellid_from_readname`.
      The expressions are copied unchanged; only the surrounding loop, which
      needs a real BED file, is dropped.

  NOT copied: `gene_assignment` and `get_cellDict` need real BED dataframes,
  so they are exercised by the end-to-end dry run (`dryrun_a9_noumi.sh`),
  not here.

Usage:
    noumi_precheck.py <R1.fastq.gz> <annotation.bed> <report.txt>
                      [stride] [libid] [maxnames] [step5.bed.gz ...]

    stride    sample every Nth read record over the WHOLE file (default 64,
              matching FS_STRIDE). The head of a fastq is not a sample of it.
    libid     library id used to build the synthetic SM tag (default ZHA8833A9)
    maxnames  cap on kept read names (default 50000)
    step5     any number of *_genes.bed.gz from step 5, checked in section 5.
              Optional: without them sections 1-4 and 6 still run.
"""
import sys
import os
import gzip
import time
from collections import Counter

# ---------------------------------------------------------------------------
# VERBATIM from a_Mapping/countTables_2pickle_cellsSpliced.py, lines 26-45.
# ---------------------------------------------------------------------------
def addCount(cnt, gene, umi, label):
    try:
        cnt[gene][umi].update([label])
    except:
        try:
            cnt[gene][umi] =  Counter([label])
        except:
            cnt[gene] = {umi: Counter([label])}
    return cnt


def get_UMI(name, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';')[1:]}
        umi = name_info['RX']
    elif protocol == 'ramda':
        umi = 'A'
    elif protocol == 'smartseq_noUMI':
        name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';')[1:]}
        umi = 'A'
    return umi


# ---------------------------------------------------------------------------
# VERBATIM from a_Mapping/countTables_2pickle_cellsSpliced.py, lines 89-92,
# lifted out of get_cellDict (see PROVENANCE above).
# ---------------------------------------------------------------------------
def cellid_from_filename(cellfile):
    return cellfile[:cellfile.index('_cbc')]


def cellid_from_readname(name):
    return {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';') if ':' in x}['SM']


# ---------------------------------------------------------------------------
# VERBATIM from a_Mapping/countTables_fromPickle.py, lines 47-85, 142-165,
# 265-271. Both protocol branches are kept in every counter -- unlike
# step7_precheck.py, which dropped the no-UMI branch because it did not need it.
# ---------------------------------------------------------------------------
def countTotalReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']:
        umi = 'A'
        y = sum(x[umi].values()) if type(x) == dict else 0
    return y


def countTotalUMI(x):
    return len(x) if type(x) == dict else 0


def countExonReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x','smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x if 'intron' not in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']:
        umi = 'A'
        y = 0
        if type(x) == dict:
            y = sum([x[umi][k] if 'exon' in k else 0 for k in x[umi]])
    return y


def countExonUMI(x):
    return len([u for u in x if 'intron' not in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0


def countIntronReads(x, protocol = 'vasa'):
    if protocol in ['vasa','10x', 'smartseq_UMI']:
        y = sum([sum(x[u].values()) for u in x if 'intron' in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0
    elif protocol in ['ramda','smartseq_noUMI']:
        umi = 'A'
        y = 0
        if type(x) == dict:
            y = sum([x[umi][k] if 'intron' in k else 0 for k in x[umi]])
    return y


def countIntronUMI(x):
    return len([u for u in x if 'intron' in ['-'.join(set(k.rsplit('-'))) for k in x[u]]]) if type(x) == dict else 0


def reduceGeneName(gene, uni_genes):
    rg = gene
    if gene.count('-') == 0:
        rg = gene
    else:
        bios = set([x.rsplit('_')[-1] for x in gene.rsplit('-')])
        shortlived = ['miRNA', 'tRNA','MtTrna']
        longstuff = ['lncRNA']
        shortstuff = ['snRNA','snoRNA','MiscRna','scaRNA']
        ribos = ['rRNA','ribozyme']
        if any([b in ribos for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in ribos])
            rg = gene
        if any([b not in shortlived for b in bios])  and any([b in shortlived for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] not in shortlived])
            rg = gene
        if any([b in shortstuff for b in bios]) and any([b not in shortstuff for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in shortstuff])
            rg = gene
        if sum([g in uni_genes for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g in uni_genes][0]
            gene = rg
        if gene.count('-') >= 1 and sum([g.rsplit('_')[1][:2]!="Gm" for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g.rsplit('_')[1][:2]!="Gm"][0]
    return rg


def remove_ENSandGm(gene):
    rg = sorted(set(['_'.join(x.rsplit("_")[1:]) for x in gene.rsplit('-')]))
    xg = [g for g in rg if g[:2] != 'Gm']
    if len(xg) == 0:
        xg = rg
    xg = '-'.join(xg)
    return xg


# ---------------------------------------------------------------------------
# TRIMMED from a_Mapping/countTables_fromPickle.py, lines 89-96: K is an
# argument here instead of a module global, so several K can be evaluated.
# ---------------------------------------------------------------------------
def bc2trans_K(x, K):
    import math
    if x >= K:
        t = math.log(1.-(float(K)-1e-3)/K)/math.log(1.-1./K)
    elif x > 0 and x < K:
        t = math.log(1.-float(x)/K)/math.log(1.-1./K)
    elif x == 0:
        t = 0
    return t


# ---------------------------------------------------------------------------
# The transforms step 6 applies to the BED's gene-name column before any of
# the fromPickle helpers ever see it. VERBATIM expressions from
# countTables_2pickle_cellsSpliced.py lines 93-97 and 99, restated as a
# function over one name (upstream applies them column-wise with df.apply).
# ---------------------------------------------------------------------------
def step6_name_transforms(raw):
    """raw = column 5 of the annotation BED. Returns (index_label, biotype,
    label) exactly as step 6 derives them, or raises."""
    g = raw.replace('-', '.') + '_tRNA' if 'tRNA' in raw else raw   # line 93
    label = g.rsplit('_')[-1]                                       # line 94
    biotype = g.rsplit('_')[-2]                                     # line 95
    index_label = '_'.join(g.rsplit('_')[:-1])                       # line 99
    return index_label, biotype, label


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    fq = sys.argv[1]
    bed = sys.argv[2]
    reportpath = sys.argv[3]
    stride = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    libid = sys.argv[5] if len(sys.argv) > 5 else 'ZHA8833A9'
    maxnames = int(sys.argv[6]) if len(sys.argv) > 6 else 50000

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
    emit('fastq      : %s' % fq)
    emit('annotation : %s' % bed)
    emit('stride     : %d   (every Nth read record over the WHOLE file)' % stride)
    emit('libid      : %s' % libid)
    for p in (fq, bed):
        if not os.path.exists(p):
            sys.exit('MISSING: %s' % p)
    emit('fastq size : %.2f GB' % (os.path.getsize(fq) / 1e9))
    emit('bed size   : %.1f MB' % (os.path.getsize(bed) / 1e6))

    # -----------------------------------------------------------------------
    section('1. READ NAMES -- what get_UMI requires, and what FLASH-seq has')
    emit("get_UMI under smartseq_noUMI (2pickle lines 42-44) is:")
    emit("    elif protocol == 'smartseq_noUMI':")
    emit("        name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] "
         "for x in name.rsplit(';')[1:]}")
    emit("        umi = 'A'")
    emit()
    emit("name_info is BUILT AND THEN DISCARDED -- umi is the literal 'A' either")
    emit("way. But the dict comprehension still runs, so the name must not make it")
    emit("raise. It raises IndexError iff some ';'-separated field AFTER the first")
    emit("contains no ':'. A name with NO ';' at all gives rsplit(';')[1:] == [] ,")
    emit("the comprehension is empty, and it cannot raise.")
    emit()

    t0 = time.time()
    names_raw = []
    nrec = 0
    with gzip.open(fq, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            if nrec % stride == 0 and len(names_raw) < maxnames:
                names_raw.append(line.rstrip('\n'))
            nrec += 1
    emit('read records scanned : %d' % nrec)
    emit('read names kept      : %d   (%.1f s)' % (len(names_raw), time.time() - t0))
    if not names_raw:
        sys.exit('no read names recovered from %s' % fq)
    emit()
    emit('first 3 raw FASTQ header lines, verbatim:')
    for n in names_raw[:3]:
        emit('  %s' % n)

    # The three forms a name can reach get_UMI in.
    forms = {
        'raw FASTQ header (with the " 1:N:0:..." suffix)':
            [n.lstrip('@') for n in names_raw],
        'STAR-truncated (first whitespace token -- what a BAM QNAME is)':
            [n.lstrip('@').split()[0] for n in names_raw],
        'SM-tagged (STAR-truncated + ";SM:%s")' % libid:
            [n.lstrip('@').split()[0] + ';SM:' + libid for n in names_raw],
    }
    emit()
    emit('%-62s %8s %8s' % ('form', 'raises', 'ok'))
    form_results = {}
    for label, ns in forms.items():
        nerr = 0
        firsterr = None
        for n in ns:
            try:
                u = get_UMI(n, 'smartseq_noUMI')
                if u != 'A':
                    nerr += 1
            except Exception as e:
                nerr += 1
                if firsterr is None:
                    firsterr = '%s: %s' % (type(e).__name__, e)
        form_results[label] = (nerr, firsterr)
        emit('%-62s %8d %8d' % (label[:62], nerr, len(ns) - nerr))
        if firsterr:
            emit('      first exception: %s' % firsterr)
    emit()
    emit('example, the STAR-truncated form as get_UMI decomposes it:')
    ex = forms['STAR-truncated (first whitespace token -- what a BAM QNAME is)'][0]
    emit('  name              = %r' % ex)
    emit('  name.rsplit(";")  = %r' % ex.rsplit(';'))
    emit('  [1:] (parsed)     = %r' % ex.rsplit(';')[1:])
    emit('  umi returned      = %r' % get_UMI(ex, 'smartseq_noUMI'))
    ext = forms['SM-tagged (STAR-truncated + ";SM:%s")' % libid][0]
    emit('  SM-tagged name    = %r' % ext)
    emit('  [1:] (parsed)     = %r' % ext.rsplit(';')[1:])
    emit('  umi returned      = %r' % get_UMI(ext, 'smartseq_noUMI'))

    # And the same names down the 'vasa' path, to show what the branch buys.
    nv = 0
    for n in forms['STAR-truncated (first whitespace token -- what a BAM QNAME is)'][:2000]:
        try:
            get_UMI(n, 'vasa')
        except Exception:
            nv += 1
    emit()
    emit("for contrast, the SAME names through protocol='vasa': %d/%d raise "
         "(KeyError 'RX')" % (nv, min(2000, len(names_raw))))

    star_form = 'STAR-truncated (first whitespace token -- what a BAM QNAME is)'
    if form_results[star_form][0] == 0:
        emit()
        emit('OK: untransformed FLASH-seq read names satisfy get_UMI under '
             'smartseq_noUMI. No transformation is required for the UMI field.')
    else:
        fail.append('untransformed FLASH-seq read names break get_UMI under '
                    'smartseq_noUMI (%d/%d)' % (form_results[star_form][0], len(names_raw)))

    # a real hazard: a ';' already inside the instrument name
    semis = [n for n in names_raw if ';' in n.split()[0]]
    if semis:
        fail.append('%d read name(s) already contain a ";" -- get_UMI would try to '
                    'parse them as VASA tag fields' % len(semis))
        emit('  OFFENDING: %s' % semis[0])
    else:
        emit('OK: no read name contains a ";" of its own (%d checked)' % len(names_raw))

    # -----------------------------------------------------------------------
    section('2. CELL ID -- what column name each cellidFROMfilename mode gives')
    emit("step 6 (2pickle lines 89-92):")
    emit("    if cellidFROMfilename == 'f':  cellid = cellfile[:cellfile.index('_cbc')]")
    emit("    else:                          cellid = {...}['SM']   # from the read name")
    emit()
    # mode 'f' -- the filename must contain '_cbc' or .index() raises ValueError
    for cand in ['%s_cbc_noumi_E99_Aligned.out.singlemappers_genes.bed.gz' % libid,
                 './%s_cbc_noumi_E99_Aligned.out.singlemappers_genes.bed.gz' % libid,
                 '%s_E99_Aligned.out.singlemappers_genes.bed.gz' % libid]:
        try:
            emit("  mode 'f'  %-64s -> %r" % (cand, cellid_from_filename(cand)))
        except Exception as e:
            emit("  mode 'f'  %-64s -> %s: %s" % (cand, type(e).__name__, e))
    emit()
    emit("  NB step 6 globs the folder you hand it and keeps the path it globbed,")
    emit("  so mode 'f' puts that folder string into the column name. There is no")
    emit("  folder argument that yields a bare library id.")
    emit()
    # mode 'r' -- needs an SM tag
    for cand in [names_raw[0].lstrip('@').split()[0],
                 names_raw[0].lstrip('@').split()[0] + ';SM:' + libid]:
        try:
            emit("  mode 'r'  %-64s -> %r" % (cand[:64], cellid_from_readname(cand)))
        except Exception as e:
            emit("  mode 'r'  %-64s -> %s: %s" % (cand[:64], type(e).__name__, e))
    emit()
    emit("  An untagged Illumina name HAS colons, so the dict builds; it simply has")
    emit("  no 'SM' key. Mode 'r' therefore needs the SM tag; mode 'f' does not.")

    # -----------------------------------------------------------------------
    section('3. NAMING CONTRACT -- the v2 BED through step 6 and reduceGeneName')
    t0 = time.time()
    raw_names = set()
    ncols = Counter()
    nlines = 0
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            ncols[len(f)] += 1
            nlines += 1
            if len(f) >= 5:
                raw_names.add(f[4])
    emit('BED lines            : %d   (%.1f s)' % (nlines, time.time() - t0))
    emit('column counts        : %s' % dict(ncols))
    emit('distinct raw names   : %d' % len(raw_names))
    if set(ncols) != {8}:
        warn.append('annotation BED does not have a uniform 8 columns: %s' % dict(ncols))

    idx_labels = set()
    biotypes = Counter()
    xform_err = []
    for raw in raw_names:
        try:
            il, bt, lb = step6_name_transforms(raw)
        except Exception as e:
            xform_err.append((raw, '%s: %s' % (type(e).__name__, e)))
            continue
        idx_labels.add(il)
        biotypes[bt] += 1
    if xform_err:
        fail.append("step 6's own name transforms raise on %d BED name(s)" % len(xform_err))
        for raw, e in xform_err[:10]:
            emit('  RAISE: %-60s %s' % (raw, e))
    else:
        emit("OK: step 6's name transforms succeed on all %d raw names" % len(raw_names))
    emit('distinct index labels: %d' % len(idx_labels))
    emit('biotype field values : %d distinct; commonest: %s'
         % (len(biotypes), ', '.join('%s=%d' % kv for kv in biotypes.most_common(8))))

    # the exact IndexError from the VASA run-1 failure: g.rsplit('_')[1]
    trna_labels = set(i for i in idx_labels if 'tRNA' in i)
    gene_labels = sorted(i for i in idx_labels if i not in trna_labels)
    emit()
    emit("index labels routed to the tRNA tables ('tRNA' in label): %d" % len(trna_labels))
    emit('index labels routed to the GENE tables                  : %d' % len(gene_labels))
    bad_fields = [g for g in gene_labels if len(g.rsplit('_')) < 2]
    if bad_fields:
        fail.append('%d gene label(s) have <2 underscore fields -- reduceGeneName '
                    'raises IndexError on any combination containing one '
                    '(this is how the VASA run died)' % len(bad_fields))
        for g in bad_fields[:10]:
            emit('  BAD: %r' % g)
    else:
        emit('OK: every gene-table label has >=2 underscore fields')

    # '-' in a single-gene label poisons the combination separator
    dashed = [g for g in gene_labels if '-' in g]
    if dashed:
        fail.append("%d gene label(s) contain '-', which is the separator step 6 "
                    "uses to join multi-gene names -- such a label is "
                    "indistinguishable from a combination" % len(dashed))
        for g in dashed[:10]:
            emit('  DASHED: %r' % g)
    else:
        emit("OK: no gene-table label contains a '-' (the combination separator)")

    # replay reduceGeneName over synthetic combinations built the way step 6
    # builds them: '-'.join(sorted set of genes hit by one read)
    emit()
    emit('reduceGeneName replay over synthetic combinations '
         '(step 6 builds them as "-".join(genes)):')
    import random
    random.seed(0)
    pool = gene_labels
    uni = set(pool)
    combos = []
    for k in (2, 3, 4):
        for _ in range(4000):
            combos.append('-'.join(random.sample(pool, k)))
    errs, empties = [], []
    for cb in combos:
        try:
            r = reduceGeneName(cb, uni)
        except Exception as e:
            errs.append((cb, '%s: %s' % (type(e).__name__, e)))
            continue
        if r == '':
            empties.append(cb)
    emit('  combinations tested : %d (seed 0, k=2,3,4)' % len(combos))
    if errs:
        fail.append('reduceGeneName raised on %d/%d synthetic combinations'
                    % (len(errs), len(combos)))
        for cb, e in errs[:5]:
            emit('  RAISE: %-70s %s' % (cb[:70], e))
    else:
        emit('  OK: no exception')
    if empties:
        warn.append('reduceGeneName returns an EMPTY name for %d/%d synthetic '
                    'combinations (upstream behaviour: they all merge into one row '
                    'named "")' % (len(empties), len(combos)))
        for cb in empties[:3]:
            emit('  EMPTY: %s' % cb[:70])
    else:
        emit('  OK: no combination reduces to an empty name')
    # remove_ENSandGm runs on every surviving label at the end of step 7
    rerr = [g for g in gene_labels[:20000] if _raises(remove_ENSandGm, g)]
    if rerr:
        fail.append('remove_ENSandGm raises on %d label(s)' % len(rerr))
    else:
        emit('  OK: remove_ENSandGm raises on none of %d single labels'
             % min(20000, len(gene_labels)))

    # -----------------------------------------------------------------------
    section('4. NO-UMI SEMANTICS -- what the counts mean when every read is "A"')
    emit('Built with addCount exactly as step 6 calls it: one gene, 7 exon reads')
    emit('and 3 intron reads, all with umi="A" because get_UMI returned "A".')
    cnt = {}
    for _ in range(7):
        cnt = addCount(cnt, 'ENSMUSG00000000001_Gnai3_ProteinCoding', 'A', 'exon')
    for _ in range(3):
        cnt = addCount(cnt, 'ENSMUSG00000000001_Gnai3_ProteinCoding', 'A', 'intron')
    x = cnt['ENSMUSG00000000001_Gnai3_ProteinCoding']
    emit('  cell entry            : %r' % x)
    emit()
    emit('%-34s %14s %14s' % ('step-7 counter', "protocol='vasa'", "'smartseq_noUMI'"))
    rows = [
        ('countTotalReads', countTotalReads(x, 'vasa'), countTotalReads(x, 'smartseq_noUMI')),
        ('countExonReads', countExonReads(x, 'vasa'), countExonReads(x, 'smartseq_noUMI')),
        ('countIntronReads', countIntronReads(x, 'vasa'), countIntronReads(x, 'smartseq_noUMI')),
    ]
    for nm, a, b in rows:
        emit('%-34s %14s %14s' % (nm, a, b))
    emit('%-34s %14s %14s' % ('countTotalUMI  (no protocol arg)', countTotalUMI(x), countTotalUMI(x)))
    emit('%-34s %14s %14s' % ('countExonUMI   (no protocol arg)', countExonUMI(x), countExonUMI(x)))
    emit('%-34s %14s %14s' % ('countIntronUMI (no protocol arg)', countIntronUMI(x), countIntronUMI(x)))
    emit()
    emit('countTotalUMI/countExonUMI/countIntronUMI have NO protocol argument')
    emit('(fromPickle lines 55, 69, 83) -- they count keys of the per-gene dict.')
    emit('With one literal key "A" that is 1 for every detected gene and 0')
    emit('otherwise, i.e. the UFICounts table degenerates to a DETECTION MASK.')
    emit()
    emit('The vasa/noUMI columns agree above only because there is a single UMI.')
    emit('The no-UMI branch is not a deduplication step: it reads x["A"] directly')
    emit('instead of summing over UMIs. It would raise KeyError on a VASA pickle')
    emit('(no "A" key); the vasa branch on a no-UMI pickle gives the same number.')
    emit()
    emit('BUT the exon/intron branches are NOT equivalent, and this matters:')
    emit("  vasa   : a UMI is exonic iff the string 'intron' is not an ELEMENT of")
    emit("           its collapsed label list -- exact membership, not substring")
    emit("           (fromPickle lines 61, 75: \"'intron' not in ['-'.join(...)]\")")
    emit("  noUMI  : a read is exonic iff 'exon' is a SUBSTRING of its label, and")
    emit("           intronic iff 'intron' is a substring (lines 65, 79)")
    emit("A combination label such as 'exon-intron' (one read spanning two genes,")
    emit("one exonically one intronically) collapses to 'exon-intron', which is not")
    emit("the literal string 'intron', so the VASA branch counts it as EXON only.")
    emit("The no-UMI branch tests substrings, so both 'exon' and 'intron' match and")
    emit('the read is counted in exon AND intron reads. Demonstrated:')
    cnt2 = {}
    for _ in range(5):
        cnt2 = addCount(cnt2, 'G1-G2', 'A', 'exon-intron')
    x2 = cnt2['G1-G2']
    emit('  entry %r' % x2)
    emit('  vasa  : total=%d exon=%d intron=%d  (exon+intron=%d)'
         % (countTotalReads(x2, 'vasa'), countExonReads(x2, 'vasa'),
            countIntronReads(x2, 'vasa'),
            countExonReads(x2, 'vasa') + countIntronReads(x2, 'vasa')))
    emit('  noUMI : total=%d exon=%d intron=%d  (exon+intron=%d)  <- DOUBLE COUNTS'
         % (countTotalReads(x2, 'smartseq_noUMI'), countExonReads(x2, 'smartseq_noUMI'),
            countIntronReads(x2, 'smartseq_noUMI'),
            countExonReads(x2, 'smartseq_noUMI') + countIntronReads(x2, 'smartseq_noUMI')))
    emit('  (step 7 writes the spliced/unspliced tables only for uni-genes plus')
    emit('   single-label multi-genes, so how much this bites is measured by the')
    emit('   dry run, not asserted here.)')

    emit()
    emit('-- the UMI ceiling, and what bc2trans does to a no-UMI pickle --')
    emit('fromPickle line 87 reads the UMI length off the FIRST UMI in the pickle:')
    emit('    umi = sorted([x for x in cntdf[cntdf.columns[0]] '
         'if type(x)==dict][0].keys())[0]')
    emit('    K = 4**len(umi)')
    emit("With the literal UMI 'A', len(umi)==1, so K = 4**1 = 4.")
    K_noumi = 4 ** len('A')
    K_vasa = 4 ** 6
    emit('  K (no-UMI pickle) = %d        K (VASA 6 nt UMI) = %d' % (K_noumi, K_vasa))
    emit()
    emit('bc2trans is then applied to the UFI table (lines 102, 239, 247-249).')
    emit('%8s %22s %22s' % ('UFI', 'bc2trans @ K=4', 'bc2trans @ K=4096'))
    for v in (0, 1, 2, 3, 4, 10, 4096, 100000):
        emit('%8d %22.6f %22.6f' % (v, bc2trans_K(v, K_noumi), bc2trans_K(v, K_vasa)))
    emit()
    emit('The only UFI values a no-UMI pickle can produce are 0 and 1, and')
    emit('bc2trans(0)=0, bc2trans(1)=%.6f at K=4. So the TranscriptCounts tables'
         % bc2trans_K(1, K_noumi))
    emit('become the detection mask scaled by a constant -- they carry no')
    emit('abundance information at all. They are NOT read counts corrected for')
    emit('collisions; they must not be used as a quantification.')
    emit()
    emit('Consequences for the analysis-set filter:')
    emit('  - the UMI-saturation ceiling CANNOT arise. Saturation needs UFI to')
    emit('    approach K; here UFI<=1 always, and the x>=K clamp branch')
    emit('    (bc2trans line 90, value %.4f at K=4) is unreachable.'
         % bc2trans_K(K_noumi, K_noumi))
    emit('  - so the 8 genes dropped on the VASA side for clamping at a constant')
    emit('    62356 have no FLASH-seq counterpart. The FLASH-seq analysis set must')
    emit('    be filtered on ReadCounts, and only ReadCounts.')

    # -----------------------------------------------------------------------
    section('5. STEP-5 BED CONTRACT -- the nM field step 6 must int()')
    # ADDED after the first dry run: this is the check that would have caught
    # the paired-end break. 2pickle line 97 does
    #     int(x['Info'].rsplit(';nM:')[1].rsplit(';jS:')[0])
    # on column 7 of the step-5 BED. `bedtools bamtobed` appends '/1' or '/2'
    # to the read name of a PAIRED-END mate, and deal_with_*.sh appends
    # ';CG:<cigar>;nM:<n>' to the QNAME BEFORE bamtobed runs -- so the mate
    # suffix lands on the END of the nM value ('nM:0/2') and int() raises.
    beds = sys.argv[7:] if len(sys.argv) > 7 else []
    if not beds:
        emit('no step-5 BED passed (argv[7:]) -- skipping.')
        emit('Pass the *_genes.bed.gz files to check this contract on real rows.')
    for bp in beds:
        emit()
        emit('BED: %s' % os.path.basename(bp))
        if not os.path.exists(bp):
            emit('  MISSING -- skipped')
            continue
        n = 0
        bad = 0
        firstbad = None
        firstok = None
        opener = gzip.open if bp.endswith('.gz') else open
        with opener(bp, 'rt') as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 7:
                    continue
                n += 1
                info = f[6]
                try:
                    int(info.rsplit(';nM:')[1].rsplit(';jS:')[0])
                    if firstok is None:
                        firstok = info
                except Exception as e:
                    bad += 1
                    if firstbad is None:
                        firstbad = (info, '%s: %s' % (type(e).__name__, e))
        emit('  rows                       : %d' % n)
        emit('  rows whose nM step 6 CANNOT int(): %d  (%.2f%%)'
             % (bad, 100.0 * bad / max(1, n)))
        if firstok:
            emit('  first parseable Info       : %r' % firstok)
        if firstbad:
            emit('  first UNPARSEABLE Info     : %r' % firstbad[0])
            emit('  the exception step 6 raises: %s' % firstbad[1])
            fail.append('%s: %d/%d rows (%.1f%%) have an nM field step 6 cannot '
                        'parse -- step 6 dies with ValueError on the first one. '
                        'This is the bedtools bamtobed /1,/2 mate suffix; it '
                        'means PAIRED-END input is not supported by this path.'
                        % (os.path.basename(bp), bad, n, 100.0 * bad / max(1, n)))
        else:
            emit('  OK: every row satisfies the nM contract')

    # -----------------------------------------------------------------------
    section('6. FILTER SENSITIVITY -- ncells = max(5, round(0.01*ncols))')
    emit('step 7 line ~223, reached only when argv[4] == "y":')
    emit('    ncells = max(5, round(0.01*len(cntdf.columns))); nreads = 1')
    emit()
    emit('%8s %8s %10s' % ('ncols', 'ncells', 'as % cols'))
    for n in (1, 2, 10, 12, 16, 384):
        nc = max(5, round(0.01 * n))
        emit('%8d %8d %9.1f%%' % (n, nc, 100.0 * nc / n))
    emit()
    emit('At ncols=1 (one library at a time) the threshold is 5 cells out of 1 --')
    emit('no unigene can pass, uni_genes_filt is empty, and reduceGeneName loses')
    emit('its "exactly one component is a known unigene" rule entirely. At')
    emit('ncols=10 it is 5/10 = 50%. Pass "n" for the dry run, and decide the')
    emit('filter deliberately downstream rather than inheriting a 384-cell default.')

    # -----------------------------------------------------------------------
    section('VERDICT')
    if fail:
        emit('FAIL (%d):' % len(fail))
        for f in fail:
            emit('  x %s' % f)
    else:
        emit('No blocking problem found for the smartseq_noUMI path.')
    if warn:
        emit()
        emit('WARN (%d) -- not blocking, but know about them:' % len(warn))
        for w in warn:
            emit('  ! %s' % w)

    d = os.path.dirname(reportpath)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(reportpath, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\nreport written: %s' % reportpath, flush=True)
    sys.exit(1 if fail else 0)


def _raises(fn, *a):
    try:
        fn(*a)
        return False
    except Exception:
        return True


if __name__ == '__main__':
    main()
