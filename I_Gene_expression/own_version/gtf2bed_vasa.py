#!/usr/bin/env python3
"""
gtf2bed_vasa.py -- Ensembl GTF -> the 8-column exon/intron BED this pipeline's
step 5 consumes.

    chrom  start  end  strand  GENEID_SYMBOL_BIOTYPE_LABEL  genelen  genestart  geneend

Called by build_annotation_bed.sh. Written to reproduce the rules of the
existing (unscripted) v1 BED exactly; those rules were reverse-engineered from
v1 itself and each one is verified by the caller's --validate mode, which
rebuilds v1 and compares it as a sorted set of rows.

THE RULES, and how each was established
---------------------------------------
1. One entry per GENE, never per transcript. Confirmed: v1 has exactly 78,348
   distinct gene ids and the GTF has exactly 78,348 `gene` features.

2. Exons are the MERGED UNION of every transcript's exons. Confirmed on Xkr4,
   whose GTF exons 3283662-3285855, 3283832-3286567 and 3284705-3287191 appear
   in v1 as the single block 3283662-3287191.

3. Introns fill the gaps between merged exons: (prevExonEnd+1, nextExonStart-1)
   in GTF 1-based coordinates. Confirmed on Xkr4 row for row.

4. genestart/geneend are the GTF `gene` feature's own span, repeated on every
   row of that gene; genelen is their difference.

5. Name is GENEID_SYMBOL_BIOTYPE_LABEL where
     SYMBOL  = gene_name with '-' replaced by '.', or the gene id when the GTF
               has no gene_name (551 such genes). '-' must go because
               countTables_2pickle_cellsSpliced.py joins combination genes with
               '-'. No gene_name contains '_' (checked), so the 4-field split
               downstream is safe.
     BIOTYPE = gene_biotype, CamelCased per underscore-separated token if it
               contains an underscore, otherwise left alone:
                   protein_coding -> ProteinCoding    lncRNA -> lncRNA
                   misc_RNA       -> MiscRna          miRNA  -> miRNA
                   Mt_tRNA        -> MtTrna           TEC    -> TEC
                   IG_V_gene      -> IgVGene
               Verified to reproduce all 37 biotypes present in v1 with no
               misses and no extras.
     LABEL   = exon | intron

COORDINATES
-----------
--coord asis  emit GTF's 1-based inclusive numbers verbatim, which is what v1
              does -- and which is wrong, because BED is 0-based half-open, so
              every feature ends up 1 bp too short at its 5' end. Only for
              reproducing v1.
--coord fix   subtract 1 from start and genestart, leaving end/geneend alone,
              giving true 0-based half-open. This is the default in the caller.
              See build_annotation_bed.sh for the measured cost of not doing it.
"""

import sys, re, argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument('gtf')
ap.add_argument('out')
ap.add_argument('--coord', choices=['fix', 'asis'], default='fix')
a = ap.parse_args()

SHIFT = 1 if a.coord == 'fix' else 0

attr_re = {k: re.compile(r'%s "([^"]*)"' % k) for k in ('gene_id', 'gene_name', 'gene_biotype')}


def attr(s, key):
    m = attr_re[key].search(s)
    return m.group(1) if m else None


def mangle(bt):
    """gene_biotype -> the form v1 embeds in the name."""
    if '_' not in bt:
        return bt
    return ''.join(t[:1].upper() + t[1:].lower() for t in bt.split('_'))


genes = {}                       # gid -> (chrom, start, end, strand, name)
exons = defaultdict(list)        # gid -> [(start, end)]

with open(a.gtf) as fh:
    for line in fh:
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        feat = f[2]
        if feat == 'gene':
            gid = attr(f[8], 'gene_id')
            sym = attr(f[8], 'gene_name') or gid
            bt = mangle(attr(f[8], 'gene_biotype') or 'NA')
            genes[gid] = (f[0], int(f[3]), int(f[4]), f[6],
                          '%s_%s_%s' % (gid, sym.replace('-', '.'), bt))
        elif feat == 'exon':
            exons[attr(f[8], 'gene_id')].append((int(f[3]), int(f[4])))

sys.stderr.write('genes: %d, genes with exons: %d\n' % (len(genes), len(exons)))

nrow = 0
noexon = 0
with open(a.out, 'w') as out:
    # Gene order follows the GTF, which is how v1 is laid out (it is not
    # globally coordinate-sorted -- overlapping genes interleave). Order is
    # irrelevant to bedtools intersect; the caller sorts anyway.
    for gid, (chrom, gstart, gend, strand, name) in genes.items():
        ex = exons.get(gid)
        if not ex:
            noexon += 1
            continue
        ex.sort()
        merged = []
        cs, ce = ex[0]
        for s, e in ex[1:]:
            if s <= ce + 1:          # overlapping OR exactly adjacent
                ce = max(ce, e)
            else:
                merged.append((cs, ce))
                cs, ce = s, e
        merged.append((cs, ce))

        gs = gstart - SHIFT
        glen = gend - gs
        rows = []
        prev_end = None
        for s, e in merged:
            if prev_end is not None:
                rows.append((prev_end + 1, s - 1, 'intron'))
            rows.append((s, e, 'exon'))
            prev_end = e
        for s, e, label in rows:
            out.write('%s\t%d\t%d\t%s\t%s_%s\t%d\t%d\t%d\n'
                      % (chrom, s - SHIFT, e, strand, name, label, glen, gs, gend))
            nrow += 1

sys.stderr.write('rows: %d (genes with no exon feature: %d)\n' % (nrow, noexon))
