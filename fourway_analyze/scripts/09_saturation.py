#!/usr/bin/env python3
###############################################################################
# 09_saturation.py -- gene-detection saturation, four datasets on one axis.
#
# THE FIGURE THE PAPER USES TO COMPARE SENSITIVITY (its Fig 1f)
# -------------------------------------------------------------
# x = reads per cell, y = genes detected per cell. It exists because a bare
# gene count is not a property of a protocol -- it is a property of how deeply
# you sequenced. The paper never quotes one without a depth attached:
#
#   "VASA-drop ... 9,825 +/- 280 and VASA-plate 9,480 +/- 1,252 detected genes
#    per cell, at a sequencing depth of 75,000 trimmed reads per cell"
#
# and states outright that the curve has not flattened there:
#
#   "Curvature of gene detection indicated that full complexity was not reached
#    for the method when 75,000 reads were allocated to each cell."
#
# That was confirmed on this data on 2026-08-06: a real cell (011, 2.44M UFI,
# 45,800 genes) downsampled to a blank well's 37,680 molecules detects 11,037
# genes against the blank's 10,992 -- 0.4% apart. In the unsaturated regime the
# gene count IS the depth, for a good cell and a bad one alike.
#
# NO RANDOM DRAWS -- THE EXPECTATION IS CLOSED-FORM
# --------------------------------------------------
# Subsampling D of a cell's N reads without replacement, a gene seen c times is
# missed with probability C(N-c,D)/C(N,D), so
#
#   E[genes at depth D] = sum_g [ 1 - C(N-c_g, D) / C(N, D) ]
#
# evaluated through lgamma: O(1) per gene per depth, exact, and identical on
# every re-run. Monte Carlo would add sampling noise to a quantity that has an
# exact answer, and at 215 units x 10 depths that noise is what a reader would
# mistake for structure.
#
# WHICH READS, AND WHY THE ABSOLUTE LEVEL IS NOT THE COUNT-TABLE NUMBER
# ---------------------------------------------------------------------
# Built from step-5 singlemapper BEDs: one read -> one gene, deduplicated,
# because bedtools emits a row per read x feature. That is READ-level detection,
# before UMI collapse and before step 6's multimapper aggregation, so the
# absolute y is ~1.4x below the step-7 tables (own75 cell 011: 33,357 here,
# 45,800 there). The SHAPE -- where the curve bends -- is what the figure is
# for, and it is unaffected: UMI collapse and multimapper rescue both act on
# reads already counted here.
#
# x is reads ASSIGNED TO A GENE, not trimmed reads as in the paper's axis.
# Assigned/trimmed is written to the QC table per dataset so the two can be
# converted; using it directly would push a per-dataset ratio into the axis.
#
# Output: tables/cross/saturation.tsv       dataset, depth, cells, mean_genes, sd
#         tables/cross/saturation_qc.tsv    per-dataset read totals + depth range
###############################################################################
import gzip
import os
import subprocess
import sys
from math import lgamma, exp, sqrt
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTROOT = os.environ.get("FOURWAY_OUT", ROOT)

# Log-spaced, and deliberately reaching 75,000 exactly -- that is the depth
# every number in the paper's sensitivity comparison is quoted at.
DEPTHS = [5000, 10000, 25000, 50000, 75000, 100000, 200000, 500000, 1000000]
MIN_CELLS = 3          # below this a mean is not worth drawing


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout


def gene_counts(bed):
    """gene -> reads, from one cell's singlemapper BED.

    One read may overlap several features of the same gene; bedtools emits a
    row for each and the rows for a read are contiguous (04_step5_biotype.sh
    checks this at runtime and exits 3 if not), so (read, gene) dedup is a
    single comparison against the previous row.
    """
    c = defaultdict(int)
    prev = None
    with gzip.open(bed, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            gid = f[5].split("_")[0]
            key = (f[3], gid)
            if key == prev:
                continue
            prev = key
            c[gid] += 1
    return c


def expected_genes(counts, N, D):
    """E[distinct genes] when D of the cell's N reads are kept."""
    if D >= N:
        return float(len(counts))
    # log C(N-c,D) - log C(N,D), regrouped so only the c-dependent terms move
    base = lgamma(N - D + 1) - lgamma(N + 1)
    tot = 0.0
    for c in counts.values():
        if N - D - c + 1 <= 0:          # gene cannot survive being missed
            tot += 1.0
            continue
        tot += 1.0 - exp(base - lgamma(N - D - c + 1) + lgamma(N - c + 1))
    return tot


def main():
    ds_sh = os.path.join(ROOT, "scripts", "datasets.sh")
    keys = sh("source %s && echo $DS_KEYS" % ds_sh).split()
    rows, qc = [], []
    for k in keys:
        label = sh("source %s && ds_label %s" % (ds_sh, k)).strip()
        units = [r.split("\t") for r in
                 sh("source %s && ds_units %s" % (ds_sh, k)).strip().split("\n") if r]
        per_cell = []            # (N, counts)
        for u in units:
            unit, stem = u[0], u[1]
            bed = stem + "_E99_Aligned.out.singlemappers_genes.bed.gz"
            if not os.path.exists(bed):
                bed = stem + ".singlemappers_genes.bed.gz"
            if not os.path.exists(bed):
                sys.stderr.write("  missing BED for %s/%s\n" % (label, unit))
                return 1
            c = gene_counts(bed)
            per_cell.append((sum(c.values()), c))
        tot = sum(n for n, _ in per_cell)
        sys.stderr.write("%-24s %3d units, %d assigned reads, median %d/cell\n"
                         % (label, len(per_cell), tot,
                            sorted(n for n, _ in per_cell)[len(per_cell) // 2]))
        deep = [n for n, _ in per_cell]
        qc.append((label, len(per_cell), tot, min(deep), max(deep),
                   sum(1 for n in deep if n >= 75000)))
        for D in DEPTHS:
            vals = [expected_genes(c, n, D) for n, c in per_cell if n >= D]
            if len(vals) < MIN_CELLS:
                continue
            m = sum(vals) / len(vals)
            sd = sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) if len(vals) > 1 else 0.0
            rows.append((label, D, len(vals), m, sd))

    tab = os.path.join(OUTROOT, "tables", "cross")
    os.makedirs(tab, exist_ok=True)
    with open(os.path.join(tab, "saturation.tsv"), "w") as out:
        out.write("dataset\tdepth\tcells\tmean_genes\tsd_genes\n")
        for r in rows:
            out.write("%s\t%d\t%d\t%.2f\t%.2f\n" % r)
    with open(os.path.join(tab, "saturation_qc.tsv"), "w") as out:
        out.write("dataset\tunits\tassigned_reads\tmin_per_cell\tmax_per_cell\tcells_ge_75k\n")
        for r in qc:
            out.write("%s\t%d\t%d\t%d\t%d\t%d\n" % r)
    sys.stderr.write("\nwrote tables/cross/saturation.tsv + saturation_qc.tsv\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
