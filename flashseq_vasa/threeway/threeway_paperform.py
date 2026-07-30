#!/usr/bin/env python3
"""threeway_paperform.py -- published VASA-plate / own VASA-plate / FLASH-seq,
rendered in the published figures' own terms.

WHAT IS THE AUTHORS' CODE AND WHAT IS MINE
==========================================
This is the whole point of the file, so it is stated first and repeated at every
function.

** THE AUTHORS' CODE ** (character-for-character from
   I_Gene_expression/b_Analysis/02_scanpy_QCxBiotype.py, which is READ-ONLY here
   and is neither imported nor modified):

     paper_ubiotype()      = add_metadata()'s adata.var['biotype'] then ['ubiotype']
     paper_id()            = add_metadata()'s adata.var['id']
     PAPER_CLASSES         = biotype_split()'s class list, verbatim
     PAPER_SMALLRNA        = ["snRNA","snoRNA","MiscRna","scaRNA",'ribozyme','miRNA']
     the TF / cofactor test = add_metadata()'s adata.var['reg'], keyed on their
                              Mus_musculus_TF.txt (1,623 symbols) and
                              Mus_musculus_TF_cofactors.txt (970 symbols)
     the Fig 2b denominator = add_metadata()'s n_counts, i.e. the FULL matrix sum
                              per unit. rRNA is NOT excluded by their code and is
                              not excluded here.

** THE AUTHORS' CODE, SECOND HAND ** (copied from
   I_Gene_expression/vasaplate_check/vp_common.py, which is this project's own
   transcription of the paper's two stated cell-calling rules -- the rules are the
   authors', the implementation is this project's):

     species_of(), classify_fig1d(), classify_methods(), MIN_UFI = 7500

** MY REIMPLEMENTATION ** (the code for these panels was never deposited; grepped
   every .py/.R in the repo. The AXES, DEPTH GRID and UNITS are taken from the
   published panels so the plots are directly comparable, but the estimator is
   mine and any difference from the paper could be mine rather than real):

     fig1f saturation -- E[genes] = sum_g [1 - (1-p)^{c_g}], binomial-thinning
     expectation, exact in expectation and deterministic (no seed to argue about).

TWO PROVENANCE DEFECTS IN THE TWO-WAY VERSION, FIXED HERE AND QUANTIFIED
=======================================================================
Found while reading code/flashseq_vasa/paperfig/paperfig_compare.py. Neither is
silently corrected: both variants are computed and both are in the output TSV.

DEFECT 1 -- the two-way figure compared two DIFFERENT TABLE FAMILIES.
  It loaded the own plate's `ZHA9292A1_uniaggGenes_total.ReadCounts.tsv`
  (222,421 rows) against FLASH-seq's `FSall10_native_total.ReadCounts.tsv`
  (270,217 rows). The first is post-gene-aggregation, the second is pre-. Read
  off its own log: "VASA 222420 entries" vs "FLASH-seq native 270217 entries".
  Aggregation is exactly what collapses multimapper combination entries, so
  mixing the families biases the combination-entry rate -- the statistic the
  Fig 1f scope flip rests on. Here ALL THREE groups use uniaggGenes.
  `table_family` records which, and --also-raw recomputes the raw variant.

DEFECT 2 -- the TF class silently included combination entries.
  The authors' `id` field held gene SYMBOLS, because their tables are
  `shortGeneNames_*` (verified: their row labels are `0610005C13Rik_lncRNA`, so
  id = idx.rsplit('_')[0] = the symbol). For a combination entry their id is
  `Sym1-Sym2`, which cannot match a single symbol in the TF list -- so THE
  AUTHORS' TF CLASS CONTAINS SINGLE-GENE ENTRIES ONLY. paperfig_compare.py's
  paper_symbol() takes parts[1] of the full ENSMUSG_Sym_Biotype label, which for
  a combination returns the FIRST gene's symbol and therefore TF-assigns
  combination entries the authors' code would have scored '-'. Both are computed:
  panel 'TF'/'Cofactor' = the fork's behaviour (unchanged, so the two-way numbers
  stay reproducible), panel 'TF_authors'/'Cofactor_authors' = the authors' exact
  semantics.

THE CONFOUND THIS FILE CANNOT REMOVE
====================================
The published plate is Ensembl 99 on GRCm38; the own plate and FLASH-seq are
Ensembl 116 on GRCm39. Gene sets, gene models and biotype assignments all differ
between releases, so a published-vs-others gap mixes protocol, pipeline AND
annotation release. Two consequences, both handled explicitly:

  1. SPECIES. The published library is HEK293T + mESC on a concatenated
     GRCh38+GRCm38 reference, so entries are ENSG* (human) or ENSMUSG* (mouse).
     For the three-way comparison only mouse cells and mouse entries are used.
     Cells are called with the authors' Fig 1d rule; the Methods rule is
     tabulated alongside because the two disagree.
  2. RELEASE. Measured, not just flagged: `release_probe` in the TSV reports the
     biotype vocabulary and per-class entry counts for each group, and
     E99-only / E116-only biotype tokens are listed in the report.

Usage:
  threeway_paperform.py --stage precheck --out DIR
  threeway_paperform.py --stage main     --out DIR [--also-raw]
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

W = '/nemo/lab/turnerj/working/guangxin/vasaseq'
FSDIR = '/nemo/lab/turnerj/scratch/zhangg/vasaseq/flashseq_vasa'

TF_TXT = f'{W}/code/I_Gene_expression/b_Analysis/Mus_musculus_TF.txt'
COF_TXT = f'{W}/code/I_Gene_expression/b_Analysis/Mus_musculus_TF_cofactors.txt'

# 4 of the own plate's 16 barcodes are blanks (no cell). paperfig_compare.py's
# BLANKS, unchanged.
BLANKS = {'001', '014', '015', '016'}

# The three groups. `species_scope` says which entries are in scope: 'mouse' for
# the published plate (mouse cells on a mixed reference), 'all' for the other two
# (mouse-only reference, so every entry is already mouse).
GROUPS = {
    'published VASA-plate': dict(
        uniagg=f'{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_uniaggGenes_total.ReadCounts.tsv',
        raw=f'{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_total.ReadCounts.tsv',
        ufi=f'{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_uniaggGenes_total.UFICounts.tsv',
        protocol='vasa', genome='GRCm38+GRCh38', annotation='Ensembl 99',
        mixed_species=True, drop_blanks=False, unit_style='srr_well'),
    'own VASA-plate': dict(
        uniagg=f'{W}/data/PM26037/out/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv',
        raw=f'{W}/data/PM26037/out/ZHA9292A1_total.ReadCounts.tsv',
        ufi=None,
        protocol='vasa', genome='GRCm39', annotation='Ensembl 116',
        mixed_species=False, drop_blanks=True, unit_style='well'),
    'FLASH-seq native': dict(
        uniagg=f'{FSDIR}/native/FSall10_native_uniaggGenes_total.ReadCounts.tsv',
        raw=f'{FSDIR}/native/FSall10_native_total.ReadCounts.tsv',
        ufi=None,
        protocol='smartseq_noUMI', genome='GRCm39', annotation='Ensembl 116',
        mixed_species=False, drop_blanks=False, unit_style='library'),
    'FLASH-seq vasalen': dict(
        uniagg=f'{FSDIR}/vasalen/FSall10_vasalen_uniaggGenes_total.ReadCounts.tsv',
        raw=f'{FSDIR}/vasalen/FSall10_vasalen_total.ReadCounts.tsv',
        ufi=None,
        protocol='smartseq_noUMI', genome='GRCm39', annotation='Ensembl 116',
        mixed_species=False, drop_blanks=False, unit_style='library'),
}

# Rule 4: ReadCounts on every side. VASA's TranscriptCounts is the better quantity
# for VASA's own biology, but FLASH-seq has no UMI, so reads are the only unit both
# protocols measure. Cross-protocol comparison => ReadCounts everywhere.
COUNT_COLUMN = 'ReadCounts'

# ---------------------------------------------------------------------------
# ** THE AUTHORS' CODE **, from b_Analysis/02_scanpy_QCxBiotype.py.
# ---------------------------------------------------------------------------
PAPER_SMALLRNA = ["snRNA", "snoRNA", "MiscRna", "scaRNA", 'ribozyme', 'miRNA']
PAPER_CLASSES = [('ProteinCoding', ["ProteinCoding"]),
                 ('lncRNA', ["lncRNA"]),
                 ('smallRNA', PAPER_SMALLRNA),
                 ('tRNA', ['tRNA'])]


def paper_ubiotype(entry):
    """** AUTHORS' CODE. ** adata.var['biotype'] then ['ubiotype'].

        biotype  = '-'.join([k.rsplit('_')[-1] for k in idx.rsplit('-')])
        ubiotype = '-'.join(sorted(set(b.rsplit('-'))))

    Note their set() de-duplicates, so a ProteinCoding+ProteinCoding combination
    entry has ubiotype 'ProteinCoding' and DOES enter the ProteinCoding class.
    That is their behaviour and it is preserved.
    """
    biotype = '-'.join([k.rsplit('_')[-1] for k in entry.rsplit('-')])
    return '-'.join(sorted(set(biotype.rsplit('-'))))


def paper_id(entry):
    """** AUTHORS' CODE. ** adata.var['id'].

        id = '-'.join([k.rsplit('_')[0] for k in idx.rsplit('-')])
    """
    return '-'.join([k.rsplit('_')[0] for k in entry.rsplit('-')])


def paper_symbol(entry):
    """FORK, from paperfig_compare.py, unchanged so its numbers stay reproducible.

    Our labels are ENSMUSG..._Symbol_Biotype, so paper_id() returns the accession
    and the symbol is the middle field. Consequence, quantified in the output:
    for a COMBINATION entry this returns the first gene's symbol, so combination
    entries can be TF-assigned. See paper_reg_authors() for the authors' semantics.
    """
    parts = entry.rsplit('_')
    return parts[1] if len(parts) >= 3 else parts[0]


def paper_symbol_authors(entry):
    """** AUTHORS' SEMANTICS, restored. ** Their tables are `shortGeneNames_*`
    (verified: row labels look like `0610005C13Rik_lncRNA`), so their
        id = '-'.join([k.rsplit('_')[0] for k in idx.rsplit('-')])
    is the SYMBOL for a single-gene entry and 'Sym1-Sym2' for a combination --
    which cannot match a single symbol in the TF list. Returning the '-'-joined
    symbol string reproduces that exactly: single genes can match, combinations
    never do.
    """
    return '-'.join([k.rsplit('_')[1] if len(k.rsplit('_')) >= 3 else k.rsplit('_')[0]
                     for k in entry.rsplit('-')])


# ---------------------------------------------------------------------------
# From vasaplate_check/vp_common.py -- this project's transcription of the
# paper's two stated cell-calling rules. The RULES are the authors' (Fig. 1d
# caption p.1781 and Methods p.18); this implementation is the project's.
# ---------------------------------------------------------------------------
MIN_UFI = 7500


def species_of(idx):
    """'human' | 'mouse' | 'mixed' | 'trna' | 'other' for one row label."""
    parts = idx.split("-")
    h = any(p.startswith("ENSG") for p in parts)
    m = any(p.startswith("ENSMUSG") for p in parts)
    if h and m:
        return "mixed"
    if h:
        return "human"
    if m:
        return "mouse"
    if "tRNA" in idx:
        return "trna"
    return "other"


def classify_fig1d(h, m, min_ufi=MIN_UFI):
    """Fig. 1d caption rule: >25% of detected UFIs from the other species = mixed."""
    tot = h + m
    keep = tot >= min_ufi
    frac_h = h.where(keep) / tot.where(keep)
    lab = pd.Series("discarded", index=h.index, dtype=object)
    lab[keep & (frac_h > 0.75)] = "human"
    lab[keep & (frac_h < 0.25)] = "mouse"
    lab[keep & (frac_h >= 0.25) & (frac_h <= 0.75)] = "mixed"
    return lab, frac_h


def classify_methods(gh, gm, h, m, min_ufi=MIN_UFI):
    """Methods p.18 rule: same 7,500-UFI gate, purity measured on GENES detected."""
    tot_u = h + m
    keep = tot_u >= min_ufi
    gtot = gh + gm
    frac_h = gh.where(keep) / gtot.where(keep)
    lab = pd.Series("discarded", index=gh.index, dtype=object)
    lab[keep & (frac_h > 0.75)] = "human"
    lab[keep & (frac_h < 0.25)] = "mouse"
    lab[keep & (frac_h >= 0.25) & (frac_h <= 0.75)] = "mixed"
    return lab, frac_h


# ---------------------------------------------------------------------------
# MINE. Depth grid and axis from the published panel; estimator mine.
# ---------------------------------------------------------------------------
FIG1F_GRID = [5000, 10000, 15000, 20000, 25000, 50000, 75000]
# Two extra rungs so the paper's two printed VASA-plate values are both testable:
# Fig 1f 9,480 +- 1,252 at 75k, Ext.Data Fig 2e 15,248 +- 1,092 at 750k.
PAPER_CHECK_GRID = [75000, 750000]

# The paper's two printed VASA-plate saturation values. BOTH ARE HEK293T, i.e.
# HUMAN, and both are "all annotated genes" -- read off p.1781 and the Fig 1
# caption:
#
#   "the HEK293T datasets for each method were downsampled to determine the gene
#    detection sensitivity ... for all annotated genes. VASA-drop showed the
#    highest sensitivity, followed by VASA-plate, with 9,825+-280 and
#    9,480+-1,252 (mean +- s.d.) detected genes per cell, respectively, at a
#    sequencing depth of 75,000 trimmed reads per cell"
#
#   caption f: "The number of detected annotated genes in HEK293T cells ...
#    plotted against the number of reads (after quality filtering, adapter
#    removal and homopolymer trimming) per cell ... Only cells that were
#    sequenced to at least 75,000 reads were used (VASA-plate: n = 174 ...)"
#
# So the reproduction check MUST use HEK293T cells and human entries. Running it
# on mESC would not be a test of the published curve. Both are reported.
PAPER_FIG1F = (9480, 1252, 174,
               'Fig 1f p.1781, VASA-plate, HEK293T (human), all annotated genes, '
               '75,000 trimmed reads/cell, n=174 cells with >=75k trimmed reads')
PAPER_EXTFIG2E = (15248, 1092, None,
                  'Ext.Data Fig 2e, VASA-plate, HEK293T (human), '
                  '750,000 trimmed reads/cell')

# "reads (after quality filtering, adapter removal and homopolymer trimming)"
# (Fig 1 caption) and Methods p.17: "trimmed with TrimGalore ... homopolymers at
# the end of the read were removed with cutadapt. In silico ribosomal depletion
# was performed by mapping the trimmed reads to ... rRNA". So the paper's
# trimmed-read count is the input to the ribosomal-depletion step -- which is
# exactly what trimmed_reads() reads out of ribo-map.log. Same quantity, not an
# approximation of it.

CHUNK = 20000


def header_of(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        return fh.readline().rstrip('\n').split('\t')


def unit_labels(cols, g):
    """Column labels -> bare unit ids.

    Three different header conventions have to land on the same key space as
    trimmed_reads(), or the Fig 1f x-axis silently goes missing:

      published  'vasaplate_out_v3/SRR14783059_001'  -> '001'
      own        '001'                                -> '001'
      FLASH-seq  'cells/ZHA8833A1'                    -> 'ZHA8833A1'

    The published rule is vp_common.normalise_columns()'s (rsplit on '_', then
    zfill(3)); the own-plate rule is paperfig_compare.load()'s (zfill(3)). Getting
    this wrong is not loud -- it just makes every trimmed_reads() lookup return
    NaN and empties the saturation panel -- so unit_style is explicit per group
    and asserted against the trimmed-read keys in the precheck."""
    out = [str(cc).split('/')[-1] for cc in cols]
    if g['unit_style'] == 'srr_well':
        out = [u.rsplit('_', 1)[-1].zfill(3) for u in out]
    elif g['unit_style'] == 'well':
        out = [u.zfill(3) for u in out]
    return out


def read_chunks(path, cols):
    dt = {cc: np.int32 for cc in cols}
    for ch in pd.read_csv(path, sep='\t', index_col=0, dtype=dt, chunksize=CHUNK):
        ch = ch[~ch.index.isna()]                      # paperfig_compare.load()
        yield ch.index.astype(str).to_numpy(), ch.to_numpy(dtype=np.float64)


def thin_expect(counts, p):
    """MINE. E[# entries with >=1 read] after binomial thinning at rate p.
    Additive over entries, which is what makes a chunked pass exact."""
    if p >= 1.0:
        return float((counts > 0).sum())
    c = counts[counts > 0]
    return float(np.sum(1.0 - np.power(1.0 - p, c)))


# ---------------------------------------------------------------------------
def load_tf_lists():
    tfs = set(pd.read_csv(TF_TXT, sep='\t')['Symbol'].astype(str))
    cofs = set(pd.read_csv(COF_TXT, sep='\t')['Symbol'].astype(str))
    return tfs, cofs


def precheck(outdir):
    """Rule 2: replay the real operations over the real data, read-only, and
    report what would break -- before committing the expensive pass."""
    log = []

    def say(s):
        print(s)
        log.append(s)

    say('=== PRECHECK: contracts the main pass depends on ===')
    tfs, cofs = load_tf_lists()
    assert len(tfs) == 1623, f'TF list is {len(tfs)}, expected 1,623'
    assert len(cofs) == 970, f'cofactor list is {len(cofs)}, expected 970'
    say(f'TF list {len(tfs)} symbols, cofactor list {len(cofs)} symbols  [authors\' files]')

    # the authors' functions must behave as documented on real labels
    e_simple = 'ENSMUSG00000060938_Rpl26_ProteinCoding'
    e_combo = ('ENSMUSG00000004263_Atn1_ProteinCoding-'
               'ENSMUSG00000107478_Gm45234_ProteinCoding')
    e_mixbio = 'ENSMUSG00000000031_H19_lncRNA-ENSMUSG00000060938_Rpl26_ProteinCoding'
    assert paper_ubiotype(e_simple) == 'ProteinCoding'
    assert paper_ubiotype(e_combo) == 'ProteinCoding', paper_ubiotype(e_combo)
    assert paper_ubiotype(e_mixbio) == 'ProteinCoding-lncRNA', paper_ubiotype(e_mixbio)
    assert paper_id(e_simple) == 'ENSMUSG00000060938'
    assert paper_symbol(e_simple) == 'Rpl26'
    assert paper_symbol(e_combo) == 'Atn1', paper_symbol(e_combo)
    assert paper_symbol_authors(e_simple) == 'Rpl26'
    assert paper_symbol_authors(e_combo) == 'Atn1-Gm45234', paper_symbol_authors(e_combo)
    say('authors\' paper_ubiotype/paper_id verified on real labels; the de-duplicating')
    say('  set() puts ProteinCoding+ProteinCoding combinations INTO ProteinCoding')
    say('DEFECT 2 confirmed live: fork symbol for the combination = %r (TF-assignable),'
        % paper_symbol(e_combo))
    say('  authors\' symbol = %r (cannot match a TF list entry)'
        % paper_symbol_authors(e_combo))

    assert species_of('ENSG00000000003_TSPAN6_ProteinCoding') == 'human'
    assert species_of(e_combo) == 'mouse'
    assert species_of('ENSG00000000003_TSPAN6_ProteinCoding-' + e_simple) == 'mixed'
    say('species_of verified: human / mouse / mixed')

    # every input resolves, and the first chunk parses
    say('')
    say('%-22s %10s %8s %9s  %s' % ('group', 'MB', 'units', 'chunk1', 'first label'))
    for name, g in GROUPS.items():
        for fam in ('uniagg', 'raw'):
            p = g[fam]
            assert os.path.exists(p), f'{name} {fam}: missing {p}'
            mb = os.path.getsize(p) / 1e6
            cols = header_of(p)[1:]
            units = unit_labels(cols, g)
            if g['drop_blanks']:
                units = [u for u in units if u not in BLANKS]
            idx, X = next(read_chunks(p, cols))
            assert X.shape[1] == len(cols), f'{name} {fam}: width {X.shape[1]} != {len(cols)}'
            say('%-22s %10.1f %8d %9d  %s' % (f'{name}/{fam}', mb, len(units),
                                              len(idx), idx[0][:52]))
        if g['ufi']:
            assert os.path.exists(g['ufi']), f'{name}: missing UFI table {g["ufi"]}'

    # DEFECT 1, measured: the two families disagree on row count and combo share
    say('')
    say('DEFECT 1, measured on chunk-free row counts is expensive; measured on the')
    say('  full tables in the main pass. Header widths agree, so the mismatch the')
    say('  two-way figure had was table FAMILY, not cell set.')

    # the published plate's trimmed-read logs, which the paper's x-axis needs
    tr = trimmed_reads()
    for name in GROUPS:
        have = sum(1 for u in tr.get(name, {}) if tr[name][u] > 0)
        say('%-22s trimmed-read counts available for %d units' % (name, have))
    assert len(tr['published VASA-plate']) == 384, len(tr['published VASA-plate'])

    # THE SILENT FAILURE THIS GUARDS. Three groups use three different header
    # conventions; if unit_labels() and trimmed_reads() disagree on the key space
    # the join returns NaN for every unit, the Fig 1f depth gate excludes
    # everything, and the panel comes out EMPTY rather than wrong. Assert the
    # join, do not trust it.
    say('')
    say('unit-key join, count-table columns vs trimmed-read logs:')
    for name, g in GROUPS.items():
        cols = header_of(g['uniagg'])[1:]
        units = [u for u in unit_labels(cols, g)
                 if not (g['drop_blanks'] and u in BLANKS)]
        keys = set(tr[name])
        matched = [u for u in units if u in keys]
        say('  %-22s %3d units, %3d matched a trimmed-read log  (e.g. %s -> %s)'
            % (name, len(units), len(matched), str(cols[0])[:40], units[0]))
        assert len(matched) == len(units), (
            f'{name}: only {len(matched)}/{len(units)} units join to a '
            f'trimmed-read log. units={units[:4]} keys={sorted(keys)[:4]}')

    # the paper's own Fig 1f cohort gate must be satisfiable on the published plate
    n75 = sum(1 for u, v in tr['published VASA-plate'].items() if v >= 75000)
    n750 = sum(1 for u, v in tr['published VASA-plate'].items() if v >= 750000)
    say('published plate: %d/384 barcodes >= 75,000 trimmed reads, %d >= 750,000'
        % (n75, n750))
    say('  (the paper used n=174 HEK293T cells at the 75k gate; ours is over ALL')
    say('   384 barcodes of both species before cell calling, so it is not yet')
    say('   comparable -- the species split happens in the main pass)')
    assert n75 > 0 and n750 > 0
    say('')
    say('PRECHECK PASSED')
    with open(f'{outdir}/threeway_precheck.txt', 'w') as fh:
        fh.write('\n'.join(log) + '\n')


def trimmed_reads():
    """Reads entering the pipeline after trimming, per unit -- the paper's Fig 1f
    x-axis ("trimmed reads per cell").

    VASA: stage 3's ribo-map.log "Number of reads:" is the trimmed fastq it was
    handed, so it is the trimmed-read count. FLASH-seq has no ribosomal-depletion
    stage, so its trimmed fastq goes straight to STAR and STAR's "Number of input
    reads" IS the trimmed-read count. Different files, same quantity."""
    out = {}
    vp = f'{W}/data/ref/fastq_vasaplate/vasaplate_out_v3'
    d = {}
    for i in range(1, 385):
        u = '%03d' % i
        f = f'{vp}/SRR14783059_{u}_cbc_trimmed_homoATCG.ribo-map.log'
        if os.path.exists(f):
            for ln in open(f):
                if ln.startswith('Number of reads:'):
                    d[u] = int(ln.split(':')[1]); break
    out['published VASA-plate'] = d

    own = f'{W}/data/PM26037/out/cells'
    d = {}
    for i in range(1, 17):
        u = '%03d' % i
        f = f'{own}/ZHA9292A1_{u}_cbc_trimmed_homoATCG.ribo-map.log'
        if os.path.exists(f):
            for ln in open(f):
                if ln.startswith('Number of reads:'):
                    d[u] = int(ln.split(':')[1]); break
    out['own VASA-plate'] = d

    for arm in ('native', 'vasalen'):
        d = {}
        cd = f'{FSDIR}/{arm}/cells'
        for f in sorted(os.listdir(cd)):
            if f.endswith('_cbc_noumi_E99_Log.final.txt'):
                u = f.replace('_cbc_noumi_E99_Log.final.txt', '')
                for ln in open(f'{cd}/{f}'):
                    if 'Number of input reads' in ln:
                        d[u] = int(ln.split('|')[1].strip()); break
        out[f'FLASH-seq {arm}'] = d
    return out


# ---------------------------------------------------------------------------
def pass1(path, cols, tfs, cofs, mixed_species):
    """Per-unit class sums, totals, detection and combination-entry share.

    Everything accumulated here is additive over entries, so a chunked pass gives
    the same answer as loading the whole table. Two entry scopes are carried:
    'all' (every row) and 'mouse' (species_of == 'mouse'; identical to 'all' for
    the two GRCm39 groups, so it is only computed when mixed_species)."""
    n = len(cols)
    scopes = ['all', 'mouse', 'human'] if mixed_species else ['all']
    acc = {s: dict(total=np.zeros(n), det=np.zeros(n), det_single=np.zeros(n),
                   det_combo=np.zeros(n), total_single=np.zeros(n),
                   cls={}, cls_n={}) for s in scopes}
    panels = [c for c, _ in PAPER_CLASSES] + ['TF', 'Cofactor',
                                              'TF_authors', 'Cofactor_authors']
    for s in scopes:
        for p in panels:
            acc[s]['cls'][p] = np.zeros(n)
            acc[s]['cls_n'][p] = 0
    biotype_tokens = {}
    nrows = 0
    nrows_combo = 0

    for idx, X in read_chunks(path, cols):
        nrows += len(idx)
        ub = np.array([paper_ubiotype(e) for e in idx])
        sym = np.array([paper_symbol(e) for e in idx])
        sym_a = np.array([paper_symbol_authors(e) for e in idx])
        is_combo = np.array(['-' in e for e in idx])
        nrows_combo += int(is_combo.sum())
        sp = (np.array([species_of(e) for e in idx]) if mixed_species
              else np.full(len(idx), 'mouse'))
        reg = np.array(['TF' if s in tfs else ('Cof' if s in cofs else '-') for s in sym])
        reg_a = np.array(['TF' if s in tfs else ('Cof' if s in cofs else '-') for s in sym_a])
        for u in ub[~is_combo]:
            biotype_tokens[u] = biotype_tokens.get(u, 0) + 1

        for s in scopes:
            m = np.ones(len(idx), bool) if s == 'all' else (sp == s)
            if not m.any():
                continue
            Xs = X[m]
            acc[s]['total'] += Xs.sum(axis=0)
            acc[s]['det'] += (Xs > 0).sum(axis=0)
            sing = ~is_combo[m]
            acc[s]['det_single'] += (Xs[sing] > 0).sum(axis=0)
            acc[s]['det_combo'] += (Xs[~sing] > 0).sum(axis=0)
            acc[s]['total_single'] += Xs[sing].sum(axis=0)
            for cls, bts in PAPER_CLASSES:
                sel = np.isin(ub[m], bts)
                acc[s]['cls'][cls] += Xs[sel].sum(axis=0)
                acc[s]['cls_n'][cls] += int(sel.sum())
            for cls, tag, rv in [('TF', 'TF', reg), ('Cofactor', 'Cof', reg),
                                 ('TF_authors', 'TF', reg_a),
                                 ('Cofactor_authors', 'Cof', reg_a)]:
                sel = rv[m] == tag
                acc[s]['cls'][cls] += Xs[sel].sum(axis=0)
                acc[s]['cls_n'][cls] += int(sel.sum())

    return acc, biotype_tokens, nrows, nrows_combo


def pass2(path, cols, mixed_species, p_by_unit):
    """MINE. Saturation: E[entries detected >=1 read] after binomial thinning,
    accumulated chunk by chunk.

    p_by_unit maps (entry_scope, count_scope, axis, depth) -> per-unit thinning
    rate array (NaN where that unit cannot reach that depth). Vectorised form of
    thin_expect over a whole chunk at once:

        sum_i [1 - (1-p_j)^{C_ij}]  =  nrows - sum_i exp(C_ij * log(1-p_j))

    An entry with C_ij = 0 contributes exp(0) = 1 and hence 0 to the sum, so the
    'counts > 0' filter in thin_expect is implicit and the identity is exact.
    p_j = 1 is handled separately because log(0) = -inf."""
    n = len(cols)
    acc = {k: np.zeros(n) for k in p_by_unit}
    for idx, X in read_chunks(path, cols):
        is_combo = np.array(['-' in e for e in idx])
        sp = (np.array([species_of(e) for e in idx]) if mixed_species
              else np.full(len(idx), 'mouse'))
        masks = {}
        for (es, cs, _, _) in p_by_unit:
            if (es, cs) in masks:
                continue
            m = np.ones(len(idx), bool) if es == 'all' else (sp == es)
            if cs == 'single_gene':
                m = m & ~is_combo
            masks[(es, cs)] = m
        sub = {k: np.ascontiguousarray(X[m]) for k, m in masks.items() if m.any()}
        for key, pv in p_by_unit.items():
            es, cs = key[0], key[1]
            Xs = sub.get((es, cs))
            if Xs is None:
                continue
            ok = np.isfinite(pv)
            if not ok.any():
                continue
            p = np.where(ok, pv, 0.0)
            full = ok & (p >= 1.0)
            part = ok & (p < 1.0)
            if part.any():
                L = np.log1p(-p[part])                      # log(1-p), p<1
                acc[key][part] += Xs.shape[0] - np.exp(Xs[:, part] * L[None, :]).sum(axis=0)
            if full.any():
                acc[key][full] += (Xs[:, full] > 0).sum(axis=0)
    return acc


# ---------------------------------------------------------------------------
def call_published_cells(outdir):
    """Assign the 384 published-plate barcodes to a species, using the paper's own
    two rules. The rules are the authors'; the implementation is vp_common.py's.

    The published library is HEK293T (human) + mESC (mouse) on a concatenated
    GRCh38+GRCm38 reference. UFICounts is used here because both of the paper's
    rules are stated in UFIs ("<7,500 UFIs were filtered out"), and this is a
    cell-calling gate, not an abundance comparison -- Rule 4's ReadCounts
    requirement is about cross-protocol quantities, which this is not."""
    g = GROUPS['published VASA-plate']
    path = g['ufi']
    cols = header_of(path)[1:]
    units = unit_labels(cols, g)
    n = len(cols)
    h = np.zeros(n); m = np.zeros(n); gh = np.zeros(n); gm = np.zeros(n)
    for idx, X in read_chunks(path, cols):
        sp = np.array([species_of(e) for e in idx])
        for tag, u_acc, g_acc in (('human', h, gh), ('mouse', m, gm)):
            sel = sp == tag
            if sel.any():
                u_acc += X[sel].sum(axis=0)
                g_acc += (X[sel] > 0).sum(axis=0)
    h = pd.Series(h, index=units); m = pd.Series(m, index=units)
    gh = pd.Series(gh, index=units); gm = pd.Series(gm, index=units)

    lab_1d, frac_u = classify_fig1d(h, m)
    lab_me, frac_g = classify_methods(gh, gm, h, m)
    calls = pd.DataFrame(dict(unit=units, ufi_human=h.values, ufi_mouse=m.values,
                              genes_human=gh.values, genes_mouse=gm.values,
                              frac_human_ufi=frac_u.values,
                              frac_human_genes=frac_g.values,
                              call_fig1d=lab_1d.values, call_methods=lab_me.values))
    calls.to_csv(f'{outdir}/threeway_published_cellcalls.tsv', sep='\t', index=False)
    print('published plate, %d barcodes:' % len(calls))
    for rule in ('call_fig1d', 'call_methods'):
        vc = calls[rule].value_counts().to_dict()
        print('  %-13s %s' % (rule.replace('call_', ''),
                              ', '.join('%s=%d' % (k, vc.get(k, 0))
                                        for k in ('mouse', 'human', 'mixed', 'discarded'))))
    return calls


def main(outdir, also_raw):
    tfs, cofs = load_tf_lists()
    assert len(tfs) == 1623 and len(cofs) == 970
    print('AUTHORS\' lists: TF %d symbols, cofactor %d symbols' % (len(tfs), len(cofs)))

    calls = call_published_cells(outdir)
    # Fig 1d is the paper's own figure rule, so it is the primary. The Methods
    # rule is carried through so the disagreement stays visible (Rule 3).
    mouse_1d = set(calls.unit[calls.call_fig1d == 'mouse'])
    mouse_me = set(calls.unit[calls.call_methods == 'mouse'])
    human_1d = set(calls.unit[calls.call_fig1d == 'human'])
    print('mESC cells: Fig1d rule n=%d, Methods rule n=%d, agree on %d'
          % (len(mouse_1d), len(mouse_me), len(mouse_1d & mouse_me)))

    tr = trimmed_reads()
    families = ['uniagg'] + (['raw'] if also_raw else [])
    frac_rows, sat_rows, probe_rows, unit_rows = [], [], [], []

    for name, g in GROUPS.items():
        for fam in families:
            path = g[fam]
            cols = header_of(path)[1:]
            units = unit_labels(cols, g)
            keep_col = np.array([True] * len(cols))
            if g['drop_blanks']:
                keep_col = np.array([u not in BLANKS for u in units])
            print('\n[%s / %s] %d columns, %d in scope'
                  % (name, fam, len(cols), keep_col.sum()))

            acc, btok, nrows, nrows_combo = pass1(path, cols, tfs, cofs,
                                                  g['mixed_species'])
            probe_rows.append(dict(
                dataset=name, table_family=fam, count_column=COUNT_COLUMN,
                annotation=g['annotation'], genome=g['genome'],
                protocol=g['protocol'], n_entries=nrows,
                n_combination_entries=nrows_combo,
                pct_combination_entries=round(100 * nrows_combo / nrows, 2),
                n_biotype_tokens=len(btok),
                biotype_tokens=';'.join('%s:%d' % (k, v)
                                        for k, v in sorted(btok.items(),
                                                           key=lambda x: -x[1]))))

            # --- Fig 2b: the authors' classes over the authors' denominator ---
            # entry scope: mouse only for the published plate (mixed reference),
            # everything for the two GRCm39 groups (already mouse-only).
            escope = 'mouse' if g['mixed_species'] else 'all'
            A = acc[escope]
            for ci, u in enumerate(units):
                if not keep_col[ci]:
                    continue
                # published plate: mESC cells only, by the paper's Fig 1d rule
                if g['mixed_species'] and u not in mouse_1d:
                    continue
                tot = A['total'][ci]
                if tot <= 0:
                    continue
                for panel in list(A['cls']):
                    frac_rows.append(dict(
                        dataset=name, table_family=fam, unit=u, panel=panel,
                        fraction=float(A['cls'][panel][ci] / tot),
                        n_entries_in_class=A['cls_n'][panel],
                        denominator_counts=float(tot),
                        entry_scope=escope, count_column=COUNT_COLUMN,
                        source=('AUTHORS code' if panel in
                                ('ProteinCoding', 'lncRNA', 'smallRNA', 'tRNA',
                                 'TF_authors', 'Cofactor_authors')
                                else 'AUTHORS code, fork symbol rule')))
                unit_rows.append(dict(
                    dataset=name, table_family=fam, unit=u,
                    assigned_reads=float(tot),
                    trimmed_reads=tr.get(name, {}).get(u, np.nan),
                    entries_detected=float(A['det'][ci]),
                    entries_detected_single_gene=float(A['det_single'][ci]),
                    entries_detected_combination=float(A['det_combo'][ci]),
                    pct_detected_combination=(
                        round(100 * A['det_combo'][ci] / A['det'][ci], 2)
                        if A['det'][ci] > 0 else np.nan),
                    entry_scope=escope,
                    published_call_fig1d=(calls.set_index('unit').call_fig1d.get(u)
                                          if g['mixed_species'] else ''),
                    published_call_methods=(calls.set_index('unit').call_methods.get(u)
                                            if g['mixed_species'] else '')))

            # --- Fig 1f: MY reimplementation, paper axes ---
            # x = trimmed reads per cell; a unit enters a depth only if it was
            # actually sequenced that deep (the paper's own gate).
            p_by_unit = {}
            axes = [('mouse_mESC', escope, FIG1F_GRID)]
            if g['mixed_species']:
                # the paper's own cohort: HEK293T, human entries
                axes.append(('human_HEK293T', 'human', FIG1F_GRID + PAPER_CHECK_GRID))
                axes.append(('mouse_mESC_papercheck', 'mouse', PAPER_CHECK_GRID))
            for axname, es, grid in axes:
                for dep in grid:
                    pv = np.full(len(cols), np.nan)
                    for ci, u in enumerate(units):
                        if not keep_col[ci]:
                            continue
                        if g['mixed_species']:
                            want = human_1d if es == 'human' else mouse_1d
                            if u not in want:
                                continue
                        t = tr.get(name, {}).get(u, np.nan)
                        if not np.isfinite(t) or t < dep:
                            continue
                        pv[ci] = dep / t
                    if np.isfinite(pv).any():
                        for cs in ('all_entries', 'single_gene'):
                            p_by_unit[(es, cs, axname, dep)] = pv
            if p_by_unit:
                sat = pass2(path, cols, g['mixed_species'], p_by_unit)
                for (es, cs, axname, dep), vec in sat.items():
                    pv = p_by_unit[(es, cs, axname, dep)]
                    for ci, u in enumerate(units):
                        if not np.isfinite(pv[ci]):
                            continue
                        sat_rows.append(dict(
                            dataset=name, table_family=fam, unit=u,
                            cohort=axname, entry_scope=es, count_scope=cs,
                            depth_trimmed_reads=dep,
                            genes=float(vec[ci]),
                            trimmed_reads=tr.get(name, {}).get(u, np.nan),
                            assigned_reads=float(acc[es]['total'][ci]),
                            thin_fraction=round(float(pv[ci]), 6),
                            source='MINE (reimplementation)'))

    frac = pd.DataFrame(frac_rows)
    sat = pd.DataFrame(sat_rows)
    probe = pd.DataFrame(probe_rows)
    units_df = pd.DataFrame(unit_rows)
    frac.to_csv(f'{outdir}/paperform_threeway.tsv', sep='\t', index=False)
    sat.to_csv(f'{outdir}/paperform_threeway_fig1f.tsv', sep='\t', index=False)
    probe.to_csv(f'{outdir}/paperform_threeway_release_probe.tsv', sep='\t', index=False)
    units_df.to_csv(f'{outdir}/paperform_threeway_units.tsv', sep='\t', index=False)

    # ---- report -----------------------------------------------------------
    print('\n=== Fig 2b, AUTHORS\' classes and AUTHORS\' denominator '
          '(uniaggGenes, ReadCounts) ===')
    u = frac[frac.table_family == 'uniagg']
    print('  %-18s %-22s %5s %8s %8s %8s' % ('panel', 'dataset', 'n', 'min',
                                             'median', 'max'))
    for panel in ['ProteinCoding', 'lncRNA', 'smallRNA', 'tRNA', 'TF',
                  'TF_authors', 'Cofactor', 'Cofactor_authors']:
        for name in GROUPS:
            d = u[(u.panel == panel) & (u.dataset == name)]
            if not len(d):
                continue
            print('  %-18s %-22s %5d %8.4f %8.4f %8.4f'
                  % (panel, name, len(d), d.fraction.min(), d.fraction.median(),
                     d.fraction.max()))

    print('\n=== DEFECT 2, quantified: TF share, fork rule vs authors\' rule ===')
    for name in GROUPS:
        a = u[(u.dataset == name) & (u.panel == 'TF')].fraction.median()
        b = u[(u.dataset == name) & (u.panel == 'TF_authors')].fraction.median()
        if np.isfinite(a) and np.isfinite(b):
            print('  %-22s fork %.4f   authors %.4f   ratio %.2fx'
                  % (name, a, b, a / b if b else np.nan))

    print('\n=== Fig 1f, paper axes, MY reimplementation ===')
    for cs in ('all_entries', 'single_gene'):
        print('\n  scope = %s' % cs)
        d = sat[(sat.table_family == 'uniagg') & (sat.count_scope == cs)
                & (sat.cohort.isin(['mouse_mESC', 'human_HEK293T']))]
        piv = d.pivot_table(index='depth_trimmed_reads',
                            columns=['dataset', 'cohort'], values='genes',
                            aggfunc='mean')
        nn = d.pivot_table(index='depth_trimmed_reads',
                           columns=['dataset', 'cohort'], values='genes',
                           aggfunc='count')
        print(piv.round(0).to_string())
        print('  n per cell of the above:')
        print(nn.to_string())

    print('\n=== THE PUBLISHED-CURVE CHECK ===')
    print('  paper: %s' % PAPER_FIG1F[3])
    print('         %d +- %d genes, n=%d' % PAPER_FIG1F[:3])
    for cs in ('all_entries', 'single_gene'):
        for dep, ref in ((75000, PAPER_FIG1F), (750000, PAPER_EXTFIG2E)):
            d = sat[(sat.dataset == 'published VASA-plate')
                    & (sat.table_family == 'uniagg') & (sat.count_scope == cs)
                    & (sat.entry_scope == 'human')
                    & (sat.depth_trimmed_reads == dep)]
            if not len(d):
                continue
            print('  ours, HEK293T, %-11s at %6dk: %8.0f +- %6.0f  n=%3d   '
                  'paper %5d +- %4d  ratio %.3f'
                  % (cs, dep // 1000, d.genes.mean(), d.genes.std(), len(d),
                     ref[0], ref[1], d.genes.mean() / ref[0]))
    print('\nwrote paperform_threeway.tsv, _fig1f.tsv, _units.tsv, _release_probe.tsv')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', required=True, choices=['precheck', 'main'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--also-raw', action='store_true',
                    help='also compute the pre-aggregation table family, to '
                         'quantify DEFECT 1')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.stage == 'precheck':
        precheck(a.out)
    else:
        main(a.out, a.also_raw)
