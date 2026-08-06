#!/usr/bin/env python3
###############################################################################
# 07_genebody_coverage.py -- 5'->3' coverage along the mature transcript, on all
# four datasets.
#
# WHAT THIS MEASURES, AND WHY IT SEPARATES THESE FOUR PROTOCOLS
# ------------------------------------------------------------
# Where along a transcript the reads land. A polyA-primed, oligo-dT protocol
# piles up at the 3' end when input is degraded; a random-primed / fragmented
# total-RNA protocol should be flat. FLASH-seq and VASA prime differently, so
# this is one of the few figures here where a protocol difference is expected
# to be visible directly rather than through a confound.
#
# BUILT FROM THE STEP-5 BEDs, NOT FROM BAMs
# -----------------------------------------
# RSeQC's geneBody_coverage.py wants BED12 transcript models and a *sorted,
# indexed* BAM. These BAMs are `--outSAMtype BAM Unsorted` and nothing
# downstream needs them indexed, so using RSeQC would mean sorting+indexing 215
# units first. The step-5 BEDs already carry, per read: genomic span, strand,
# gene, and the CIGAR. That is everything needed, it is one streaming pass, and
# it scores exactly the reads the rest of this folder scores.
#
# EXON-ONLY, AND THE EXCLUDED FRACTION IS REPORTED
# ------------------------------------------------
# Gene body coverage is defined on the MATURE transcript, so intronic reads have
# no transcript coordinate and are dropped. That is the standard definition and
# also the only one that makes a 5'/3' axis meaningful.
#
# But VASA deliberately captures unspliced RNA -- dropping introns silently
# would hide a real protocol difference behind a QC plot. So the intronic
# fraction is written to genebody_qc.tsv per dataset and belongs in the caption.
#
# EACH DATASET IS SCORED AGAINST THE ANNOTATION IT WAS MAPPED TO
# --------------------------------------------------------------
# own130/own75/fs -> GRCm39 E116; plate -> the mixed GRCh38+GRCm38 E99, mouse
# contigs only. Building one shared model would mean re-deriving exon structure
# in coordinates the plate's reads were never aligned in. Since only the
# RELATIVE position (0-100%) is compared, using each dataset's own model is both
# correct and the only honest option. Genes are matched across the two releases
# by Ensembl gene id, which is stable.
#
# UNION-EXON MODEL, NOT PER-ISOFORM
# ---------------------------------
# One model per gene: the union of its exons. A gene with alternative first or
# last exons smears across the axis. Per-isoform assignment is not recoverable
# from a step-5 BED (the read was assigned to a GENE), so this is a limit of the
# input, not a shortcut -- and it is what a union-exon QC plot always means.
#
# Output: tables/cross/genebody_coverage.tsv   dataset, bin(1-100), coverage
#         tables/cross/genebody_qc.tsv         per-dataset genes/reads/intronic
###############################################################################
import gzip
import os
import re
import subprocess
import sys
from bisect import bisect_right
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Products may go somewhere other than the checkout -- see datasets.sh's
# OUTROOT. Same split for the same reason: ROOT still has to find datasets.sh.
OUTROOT = os.environ.get("FOURWAY_OUT", ROOT)
NBIN = 100
MIN_EXONIC = 1000      # shorter genes cannot fill 100 bins meaningfully
MAX_EXONIC = 15000     # very long genes are dominated by a few isoforms
MIN_READS_PER_GENE = 50    # per dataset, before a gene is used at all

# CIGAR ops that consume the REFERENCE (so they place bases on the genome).
# N is a splice gap: it consumes reference but carries no read base, so it is
# excluded here -- that is exactly what keeps a spliced read from painting
# coverage across its intron.
REF_CONSUMING = set("MD=X")


def sh(cmd):
    # NOT subprocess.run(capture_output=...): that keyword is Python 3.7+, and
    # the default python3 on the compute nodes is 3.6. Both scanners died on it
    # at 00:00 on 2026-08-06 (jobs 51392134/51392135) despite an identical
    # earlier run passing, because the interpreter the shebang resolves to
    # depends on which PATH the job inherited. stdout=PIPE works on both.
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
    out, _ = p.communicate()
    return out


def load_models(bed, contig_prefix=None):
    """gene -> (strand, [(start,end)...] sorted, [cumulative offsets], total_len)

    Union of exons per gene. Overlapping exons from different isoforms are
    merged, so a base is counted once however many isoforms contain it.
    """
    ex = defaultdict(list)
    strand = {}
    opener = gzip.open if bed.endswith(".gz") else open
    with opener(bed, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or not f[4].endswith("_exon"):
                continue
            if contig_prefix and not f[0].startswith(contig_prefix):
                continue
            name = f[4]
            parts = name.split("_")
            if len(parts) < 4 or parts[2] != "ProteinCoding":
                continue
            gid = parts[0]
            ex[gid].append((int(f[1]), int(f[2])))
            strand[gid] = f[3]

    models = {}
    for gid, iv in ex.items():
        iv.sort()
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]:
                if e > merged[-1][1]:
                    merged[-1][1] = e
            else:
                merged.append([s, e])
        total = sum(e - s for s, e in merged)
        if not (MIN_EXONIC <= total <= MAX_EXONIC):
            continue
        cum, run = [], 0
        for s, e in merged:
            cum.append(run)
            run += e - s
        # `starts` is precomputed, not derived per lookup: to_transcript runs
        # once per aligned block over tens of millions of reads, and rebuilding
        # it there turned an O(log n) bisect into an O(n) list comprehension.
        starts = [s for s, _ in merged]
        models[gid] = (strand[gid], merged, cum, total, starts)
    return models


def blocks_from_cigar(start, cig):
    """Reference-consuming blocks of the alignment, as [(s,e)...].

    `CG:37M124327N44M1D59M;nM:5;jS:3` -> the field is stripped to the CIGAR
    first. N advances the reference without emitting a block, which is what
    separates the two exonic halves of a spliced read.
    """
    cig = cig.split(";")[0]
    if cig.startswith("CG:"):
        cig = cig[3:]
    out, pos, num = [], start, ""
    for ch in cig:
        if ch.isdigit():
            num += ch
        else:
            n = int(num) if num else 0
            num = ""
            if ch in REF_CONSUMING:
                out.append((pos, pos + n))
                pos += n
            elif ch == "N":
                pos += n
            # S/H/I/P consume no reference
    # merge adjacent blocks (M1D59M is one contiguous stretch on the reference)
    merged = []
    for s, e in out:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def to_transcript(model, gs, ge):
    """Genomic [gs,ge) -> list of transcript-coordinate intervals."""
    _, exons, cum, _, starts = model
    out = []
    # exons are sorted and non-overlapping; find the first that could overlap
    i = bisect_right(starts, gs) - 1
    if i < 0:
        i = 0
    for j in range(i, len(exons)):
        s, e = exons[j]
        if s >= ge:
            break
        lo, hi = max(s, gs), min(e, ge)
        if lo < hi:
            out.append((cum[j] + lo - s, cum[j] + hi - s))
    return out


def scan_dataset(label, units, models, cov, qc):
    for unit, map_stem in units:
        bed = map_stem + "_E99_Aligned.out.singlemappers_genes.bed.gz"
        if not os.path.exists(bed):
            bed = map_stem + ".singlemappers_genes.bed.gz"
        if not os.path.exists(bed):
            sys.stderr.write("  missing singlemappers BED for %s\n" % unit)
            return False
        prev_key = None
        with gzip.open(bed, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 7:
                    continue
                name = f[5]
                parts = name.split("_")
                if len(parts) < 4 or parts[2] != "ProteinCoding":
                    continue
                gid = parts[0]
                label_kind = parts[3]
                # bedtools emits one row per read x overlapping feature, and the
                # rows for one read are contiguous (script 04 asserts this). The
                # same read must not be counted once per exon it touches.
                key = (f[3], gid)
                if key == prev_key:
                    continue
                prev_key = key
                if label_kind != "exon":
                    if label_kind == "intron":
                        qc[label]["intronic"] += 1
                    continue
                m = models.get(gid)
                if m is None:
                    continue
                qc[label]["exonic"] += 1
                strand, _, _, total, _ = m
                for bs, be in blocks_from_cigar(int(f[1]), f[6]):
                    for ts, te in to_transcript(m, bs, be):
                        if strand == "-":
                            ts, te = total - te, total - ts
                        b0 = ts * NBIN // total
                        b1 = (te - 1) * NBIN // total
                        for b in range(b0, min(b1, NBIN - 1) + 1):
                            cov[label][gid][b] += 1
    return True


def main():
    ds_sh = os.path.join(ROOT, "scripts", "datasets.sh")
    keys = sh("source %s && echo $DS_KEYS" % ds_sh).split()
    # Each dataset's own annotation -- the one its reads were actually aligned
    # against. Overridable, but the defaults are the references config.sh and
    # submit_vasaplate_map_array.sh really used.
    VR = "/nemo/lab/turnerj/working/guangxin/reference/vasaseq"
    e116 = os.environ.get(
        "E116_BED",
        VR + "/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed")
    e99 = os.environ.get(
        "E99_BED",
        VR + "/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed")
    for p in (e116, e99):
        if not os.path.exists(p):
            sys.stderr.write("missing annotation BED: %s\n" % p)
            sys.exit(1)
    refs = {
        "own130": (e116, None),
        "own75": (e116, None),
        "fs": (e116, None),
        "plate": (e99, "GRCm38_"),   # mouse contigs only, this is a mixed ref
    }

    cov = {}
    qc = {}
    labels = {}
    cache = {}
    for k in keys:
        label = sh("source %s && ds_label %s" % (ds_sh, k)).strip()
        labels[k] = label
        rows = sh("source %s && ds_units %s" % (ds_sh, k)).strip().split("\n")
        units = [(r.split("\t")[0], r.split("\t")[1]) for r in rows if r]
        bedref, prefix = refs[k]
        if (bedref, prefix) not in cache:
            sys.stderr.write("building models from %s\n" % os.path.basename(bedref))
            cache[(bedref, prefix)] = load_models(bedref, prefix)
        models = cache[(bedref, prefix)]
        sys.stderr.write("%-24s %d units, %d gene models\n" % (label, len(units), len(models)))
        cov[label] = defaultdict(lambda: [0] * NBIN)
        qc[label] = defaultdict(int)
        if not scan_dataset(label, units, models, cov, qc):
            sys.exit(1)
        qc[label]["genes_covered"] = len(cov[label])

    # Only genes with real depth in EVERY dataset. Comparing curves built from
    # different gene sets would put gene composition on the same axis as
    # positional bias, which is the whole thing this figure is trying to isolate.
    shared = None
    for label in cov:
        ok = {g for g, b in cov[label].items() if sum(b) >= MIN_READS_PER_GENE}
        shared = ok if shared is None else (shared & ok)
    sys.stderr.write("\n%d genes pass >=%d reads in all four datasets\n"
                     % (len(shared), MIN_READS_PER_GENE))
    if len(shared) < 100:
        sys.stderr.write("too few shared genes -- refusing to write a curve\n")
        sys.exit(1)

    tab = os.path.join(OUTROOT, "tables", "cross")
    os.makedirs(tab, exist_ok=True)
    with open(os.path.join(tab, "genebody_coverage.tsv"), "w") as out:
        out.write("dataset\tbin\tcoverage\tgenes\n")
        for k in keys:
            label = labels[k]
            # Per-gene normalisation BEFORE averaging: without it the curve is
            # whatever the few most-expressed genes do.
            acc = [0.0] * NBIN
            for g in shared:
                b = cov[label][g]
                t = sum(b)
                if t:
                    for i in range(NBIN):
                        acc[i] += b[i] / t
            n = len(shared)
            for i in range(NBIN):
                out.write("%s\t%d\t%.8f\t%d\n" % (label, i + 1, acc[i] / n, n))
    with open(os.path.join(tab, "genebody_qc.tsv"), "w") as out:
        out.write("dataset\texonic_reads\tintronic_reads\tpct_intronic\tgenes_covered\tgenes_used\n")
        for k in keys:
            label = labels[k]
            e, i = qc[label]["exonic"], qc[label]["intronic"]
            out.write("%s\t%d\t%d\t%.2f\t%d\t%d\n"
                      % (label, e, i, 100.0 * i / (e + i) if (e + i) else 0,
                         qc[label]["genes_covered"], len(shared)))
    sys.stderr.write("wrote tables/cross/genebody_coverage.tsv + genebody_qc.tsv\n")


if __name__ == "__main__":
    main()
