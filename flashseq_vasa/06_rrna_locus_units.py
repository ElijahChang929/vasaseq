#!/usr/bin/env python3
"""
06_rrna_locus_units.py -- redo the rRNA-biotype locus attribution in the RIGHT unit.

Script 05 counted FRAGMENTS (-f 0x40) and divided by the featureCounts biotype
total, giving ~50% for Rn18s-rs5. That 50% is a unit artefact: the ratio
(featureCounts rRNA total)/(my fragment count) came out 2.0009 for A9, which is
the tell. nf-core's biotype table counts each MATE, i.e. it is in reads, so the
comparison has to be reads-to-reads. Both units are reported here so the
arithmetic is checkable, and Mt_rRNA is included as an independent confirmation
of the unit (its two loci should reconcile with the Mt_rRNA row).

A region query is an UPPER BOUND on what featureCounts would assign: samtools
counts any overlap, while featureCounts -B -C requires exon overlap with both
mates assigned and no biotype ambiguity. Stated on every row.
"""
import gzip, os, re, subprocess, sys
import pandas as pd, numpy as np

W   = "/nemo/lab/turnerj/working/guangxin/vasaseq"
FSR = f"{W}/data/flashseq/results"
SRC_GTF = "/nemo/lab/turnerj/working/guangxin/reference/gtf/Mus_musculus.GRCm39.116.gtf.gz"
OUT = os.environ.get("OUTDIR", ".")
ST  = os.environ.get("SAMTOOLS", "samtools")
SC  = ["ZHA8833A9", "ZHA8833A10"]
log = []
def note(m): log.append(str(m)); print(m, file=sys.stderr, flush=True)

gid = re.compile(r'gene_id "([^"]+)"'); gbt = re.compile(r'gene_biotype "([^"]+)"')
gnm = re.compile(r'gene_name "([^"]+)"')
loci, src_bt = [], {}
with gzip.open(SRC_GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene": continue
        g = gid.search(f[8]).group(1); b = gbt.search(f[8]).group(1)
        src_bt[g] = b
        if b in ("rRNA", "Mt_rRNA"):
            n = gnm.search(f[8]); n = n.group(1) if n else g
            loci.append((g, n, b, f[0], int(f[3]), int(f[4])))

def q(bam, regs, extra):
    tot = 0
    for i in range(0, len(regs), 300):
        tot += int(subprocess.run([ST, "view", "-c", "-F", "0x100", "-q", "255"] + extra
                                  + [bam] + regs[i:i + 300],
                                  capture_output=True, text=True, check=True).stdout.strip())
    return tot

rows = []
for lib in SC:
    bam = f"{FSR}/star_rsem/{lib}.markdup.sorted.bam"
    fcd = {}
    for line in open(f"{FSR}/star_rsem/featurecounts/{lib}.biotype_counts_mqc.tsv"):
        if line.startswith("#"): continue
        k, v = line.rstrip("\n").split("\t"); fcd[k] = float(v)
    g = pd.read_csv(f"{FSR}/star_rsem/{lib}.genes.results", sep="\t").set_index("gene_id")
    g["biotype"] = pd.Series(src_bt)
    named = [(gi, n, b, c, s, e) for gi, n, b, c, s, e in loci if (e - s + 1) >= 400]
    for gi, n, b, c, s, e in named:
        reg = [f"{c}:{s}-{e}"]
        rows.append(dict(library=lib, gene_name=n, gene_id=gi, biotype=b, contig=c,
                         start=s, end=e, gene_len_bp=e - s + 1,
                         unique_reads=q(bam, reg, []),
                         unique_fragments=q(bam, reg, ["-f", "0x40"]),
                         rsem_expected_count=float(g.expected_count.get(gi, np.nan)),
                         featurecounts_biotype_row_reads=fcd.get(b, np.nan)))
    for b in ("rRNA", "Mt_rRNA"):
        short = [f"{c}:{s}-{e}" for gi, n, bb, c, s, e in loci
                 if bb == b and (e - s + 1) < 400]
        if not short: continue
        rows.append(dict(library=lib, gene_name=f"[{len(short)} {b} loci <400 bp, pooled]",
                         gene_id="-", biotype=b, contig="-", start=-1, end=-1, gene_len_bp=-1,
                         unique_reads=q(bam, short, []),
                         unique_fragments=q(bam, short, ["-f", "0x40"]),
                         rsem_expected_count=np.nan,
                         featurecounts_biotype_row_reads=fcd.get(b, np.nan)))

d = pd.DataFrame(rows)
d["reads_per_fragment"] = d.unique_reads / d.unique_fragments.replace(0, np.nan)
d["pct_of_featurecounts_row_READS"] = 100 * d.unique_reads / d.featurecounts_biotype_row_reads
d["note"] = ("region overlap = UPPER BOUND on featureCounts assignment "
             "(featureCounts -B -C also requires exon overlap, both mates assigned, "
             "no biotype ambiguity)")
d.to_csv(f"{OUT}/rrna_locus_attribution.tsv", sep="\t", index=False)
note("rRNA/Mt_rRNA locus attribution, READS vs FRAGMENTS:\n" + d[[
    "library","gene_name","biotype","gene_len_bp","unique_reads","unique_fragments",
    "reads_per_fragment","featurecounts_biotype_row_reads",
    "pct_of_featurecounts_row_READS","rsem_expected_count"]].round(3).to_string(index=False))
for lib in SC:
    s = d[d.library == lib]
    for b in ("rRNA", "Mt_rRNA"):
        sub = s[s.biotype == b]
        note(f"{lib} {b}: sum of locus reads {sub.unique_reads.sum():,.0f} vs "
             f"featureCounts row {sub.featurecounts_biotype_row_reads.iloc[0]:,.0f} "
             f"(ratio {sub.unique_reads.sum()/sub.featurecounts_biotype_row_reads.iloc[0]:.4f}) "
             f"-- >1 is the expected direction for an overlap upper bound")
with open(f"{OUT}/locus_units_log.txt", "w") as fh: fh.write("\n".join(log) + "\n")
