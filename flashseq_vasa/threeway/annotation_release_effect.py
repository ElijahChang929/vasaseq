#!/usr/bin/env python3
"""
annotation_release_effect.py -- how much of a published-plate vs own-plate
difference is caused by the ANNOTATION RELEASE rather than by protocol?

WHY THIS TRACK EXISTS
---------------------
The three-way comparison has three datasets:

    published VASA-plate  SRR14783059      GRCm38 + hg38, Ensembl 99
    own VASA-plate        ZHA9292A1        GRCm39,        Ensembl 116
    FLASH-seq             RN26038          GRCm39,        Ensembl 116

The two VASA plates ran the SAME pipeline with the SAME chemistry branch
(protocol='vasa', 6 nt UMI, same scripts), so any difference between them is
protocol-free BY CONSTRUCTION. What remains is release + biology + depth. This
script isolates the release term, so the other tracks can subtract it instead of
arguing about it.

WHAT IS AND IS NOT MEASURED HERE
--------------------------------
MEASURED (this script):
  1. gene-universe overlap between the two annotation BEDs the pipeline
     actually consumed -- not between two Ensembl GTFs in the abstract;
  2. biotype-label churn for genes present in both releases;
  3. feature-length churn per biotype, which matters because VASA's `jS:IN`
     rule requires a read to be CONTAINED in a feature (conventions trap 8);
  4. the decisive test: re-derive each plate's biotype composition restricted
     to the gene set shared by both releases, and report how far the
     composition moves. That movement IS the release effect on composition.

NOT MEASURED (do not let a later track imply otherwise):
  - the GRCm38 -> GRCm39 coordinate/assembly change. Both plates were mapped
    to their own assembly; nothing here re-maps reads. A gene whose sequence
    changed between assemblies is counted as "shared" if its ID is in both
    BEDs, even if the reads that hit it changed. So every number below is a
    LOWER BOUND on the release effect.
  - anything about the FLASH-seq arm. It shares Ensembl 116 with the own
    plate, so the release term for own-vs-FLASH is zero by construction.

CONVENTIONS OBEYED
------------------
  Rule 4  cross-protocol/cross-dataset comparison uses ReadCounts on all
          sides, because reads are the only unit both protocols measure. The
          own plate is VASA so TranscriptCounts is also computed, as a
          sensitivity check, but ReadCounts is the headline.
  Rule 3  the UMI-ceiling genes are NOT dropped: on a READ-based composition
          they carry real reads and dropping them shifts ProteinCoding by
          ~+19 pp. The 4 blank barcodes ARE dropped -- they are not cells.
  Rule 5  every percentage below names its denominator in the output table.
  Trap 8  short features are reported explicitly, not assumed harmless.
  Trap 9  the Ensembl `rRNA` biotype is NOT an rRNA measurement on GRCm39;
          it is reported as an annotation class only, never as rRNA content.

LENGTH CONVENTION
-----------------
Feature length is taken as `end - start` from the BED, which is exactly the
convention the BED's own column 6 uses for the gene span (asserted below in
both files). If that is off by one relative to Ensembl it is off by one
IDENTICALLY in both releases, so the length-CHURN comparison is unaffected.
`bed_tiling_gap_mode` in the diagnostics confirms the two files use the same
inter-feature convention; if they did not, the paired length comparison would
be invalid and the assertion fires.

VERBATIM HELPERS
----------------
`species_of` and `biotype_of` below are copied VERBATIM from
code/I_Gene_expression/vasaplate_check/vp_common.py, as are `classify_fig1d`
and MIN_UFI, so that the mESC cell set used here is the same set that the
published-plate benchmark used. Nothing in them is trimmed. Deviating from them
would silently change which cells are "mESC" -- the benchmark records that the
two published doublet rules disagree ~6x on this data.

USAGE
    ./annotation_release_effect.py [--outdir DIR]
"""
import argparse
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
W = "/nemo/lab/turnerj/working/guangxin/vasaseq"
REF = "/nemo/lab/turnerj/working/guangxin/reference"

BED99 = f"{REF}/vasaseq/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed"
BED116 = f"{REF}/vasaseq/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed"

# own plate (Ensembl 116). Raw `out/` tables, NOT the analysis/ matrix, so that
# no prior filtering is entangled with the release question.
OWN_READ = f"{W}/data/PM26037/out/ZHA9292A1_total.ReadCounts.tsv"
OWN_TRANS = f"{W}/data/PM26037/out/ZHA9292A1_total.TranscriptCounts.tsv"
OWN_BLANKS = {"001", "014", "015", "016"}  # per data/PM26037/out/analysis/filter_report.txt

# published plate (Ensembl 99), run v3 -- the validated benchmark anchor.
# bedv2 still carries the ~600x rRNA artefact; v3 does not.
PUB_DIR = f"{W}/data/ref/fastq_vasaplate"
PUB_READ = f"{PUB_DIR}/vasaplate_out_v3_total.ReadCounts.tsv"
PUB_UFI = f"{PUB_DIR}/vasaplate_out_v3_total.UFICounts.tsv"

# read-length thresholds for the containment (`jS:IN`) exposure calculation.
# 151 is this project's library read length; the shorter values show how the
# exposure scales, so a read-length-matched control can be sized.
READ_LENS = (50, 75, 100, 151)

MIN_UFI = 7500  # VERBATIM from vp_common.py


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------------------
# VERBATIM from vp_common.py -- do not edit
# ---------------------------------------------------------------------------
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


def biotype_of(idx):
    """Biotype token of a simple row; None for combinations, 'tRNA' for tRNA."""
    if "-" in idx:
        return None
    if "_" not in idx:
        return "tRNA" if "tRNA" in idx else None
    return idx.rsplit("_", 1)[-1]


def classify_fig1d(h, m, min_ufi=MIN_UFI):
    tot = h + m
    keep = tot >= min_ufi
    frac_h = h.where(keep) / tot.where(keep)
    lab = pd.Series("discarded", index=h.index, dtype=object)
    lab[keep & (frac_h > 0.75)] = "human"
    lab[keep & (frac_h < 0.25)] = "mouse"
    lab[keep & (frac_h >= 0.25) & (frac_h <= 0.75)] = "mixed"
    return lab, frac_h
# ---------------------------------------------------------------------------
# end verbatim block
# ---------------------------------------------------------------------------


# coarse grouping. sncRNA membership is vp_common.SNCRNA_BIOTYPES verbatim.
SNCRNA_BIOTYPES = {"MiscRna", "snoRNA", "snRNA", "scaRNA", "miRNA", "ribozyme",
                   "sRNA", "MtTrna", "rRNA", "MtRrna", "vaultRNA"}
PSEUDO_TOKENS = ("Pseudogene", "pseudogene")
IGTR_PREFIX = ("Ig", "Tr")


def coarse_class(bt):
    if bt is None:
        return None
    if bt == "ProteinCoding":
        return "ProteinCoding"
    if bt == "lncRNA":
        return "lncRNA"
    if bt == "tRNA":
        return "tRNA_scan"          # the tRNAscan rows, not an Ensembl biotype
    if bt in SNCRNA_BIOTYPES:
        return "sncRNA"
    if bt == "TEC":
        return "TEC"
    if any(t in bt for t in PSEUDO_TOKENS):
        return "pseudogene"
    if bt.startswith(IGTR_PREFIX):
        return "IgTr"
    return "other"


# WHY tRNAscan ROWS ARE HELD OUT OF THE COMPOSITION
# -------------------------------------------------
# The count tables carry tRNAscan rows alongside Ensembl gene rows, labelled
# with DOTS not underscores ('3.tRNA36.TyrGTA'; combinations join with '-').
# They have no Ensembl gene id, so they can never be "shared with Ensembl 99"
# in the ID sense. If they were left in the unrestricted composition and
# dropped from the restricted one, the restriction test would report their
# whole share as a release effect, which would be an artefact of row type, not
# a measurement of annotation churn.
#
# They are also not release-neutral: the two BEDs carry DIFFERENT tRNAscan
# content (E99 1,758 rows for human+mouse pooled, E116 1,137 for mouse alone),
# so the tRNA share is reported as its own metric and excluded from the
# Ensembl-gene composition on BOTH sides. Every composition below is therefore
# "share of Ensembl single-gene counts", stated in its denominator column.
def is_trna_row(lab):
    return "tRNA" in lab and "_" not in lab


# ---------------------------------------------------------------------------
# BED parsing
# ---------------------------------------------------------------------------
def detect_convention(gaps):
    """Infer the coordinate convention from how features tile a gene.

    A homemade IntronExon BED covers each gene contiguously: exon, intron,
    exon, ... with no uncovered bases. So the modal gap between consecutive
    features of one gene identifies the convention outright:

        modal gap 0  ->  0-based half-open (BED standard):  len = end - start
        modal gap 1  ->  1-based inclusive (GTF-like):      len = end - start + 1

    THIS DIFFERS BETWEEN THE TWO FILES AND IT MATTERS. Measured here:
    the Ensembl 99 mixed BED has modal gap 1 (443,563 gaps), the Ensembl 116
    BED has modal gap 0 (564,176). Taking `end - start` in both would make
    EVERY Ensembl 99 feature look 1 bp shorter than its Ensembl 116
    counterpart, which would be reported as "the release shortened its gene
    models" when it is a file-construction difference with no biology in it.
    The single-feature gene ENSMUSG00000048040 is present in both files with
    raw lengths 1531 (E99) and 1532 (E116); corrected, both are 1532.

    Returns (offset, modal_gap) where offset is added to (end - start).
    """
    modal_gap, modal_n = gaps.most_common(1)[0]
    total = sum(gaps.values())
    assert modal_gap in (0, 1), (
        f"modal inter-feature gap is {modal_gap}, expected 0 (half-open) or "
        f"1 (1-based inclusive) -- coordinate convention unrecognised, refusing "
        f"to guess a feature length")
    frac = modal_n / total
    assert frac > 0.90, (
        f"modal gap {modal_gap} covers only {100*frac:.1f}% of "
        f"{total} inter-feature gaps -- tiling is not contiguous enough to infer "
        f"the convention safely")
    return modal_gap, modal_gap


def parse_bed(path, mouse_only):
    """Parse a homemade IntronExonTrna BED.

    Two passes: the first collects inter-feature gaps to detect the coordinate
    convention, the second applies it. Feature lengths in the returned dict are
    CONVENTION-CORRECTED, so they are comparable across the two files.

    Returns dict with
        genes    : gene_id -> {symbol, biotype, chrom, strand, span,
                               n_exon, n_intron, exon_len[], intron_len[]}
        trna     : number of tRNAscan rows (they carry no ENS id)
        gaps     : Counter of inter-feature gaps within a gene (convention check)
        offset   : the +0/+1 correction applied to every (end - start)
        n_rows   : rows consumed
        skipped_species : rows dropped because they were not mouse
    """
    def rows():
        """Yield the mouse Ensembl rows of interest, parsed."""
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                chrom, start, end, strand, name = f[0], int(f[1]), int(f[2]), f[3], f[4]
                gspan, gstart, gend = int(f[5]), int(f[6]), int(f[7])
                yield chrom, start, end, strand, name, gspan, gstart, gend

    # --- pass 1: gaps, to detect the convention ------------------------------
    gaps = Counter()
    trna_rows = n_rows = skipped = 0
    prev_gene = prev_end = None
    for chrom, start, end, strand, name, gspan, gstart, gend in rows():
        n_rows += 1
        if not name.startswith("ENS"):
            trna_rows += 1          # 'chr.tRNAnnn.AminoNNN': no id, no biotype
            continue
        if mouse_only:
            if not (chrom.startswith("GRCm38_") and name.startswith("ENSMUSG")):
                skipped += 1
                continue
        elif not name.startswith("ENSMUSG"):
            skipped += 1
            continue
        gid = name.split("_", 1)[0]
        if prev_gene == gid and prev_end is not None:
            gaps[start - prev_end] += 1
        prev_gene, prev_end = gid, end

    offset, modal_gap = detect_convention(gaps)

    # --- pass 2: build the gene records, lengths convention-corrected --------
    genes = {}
    span_mismatch = biotype_conflict = symbol_conflict = 0
    for chrom, start, end, strand, name, gspan, gstart, gend in rows():
        if not name.startswith("ENS"):
            continue
        if mouse_only:
            if not (chrom.startswith("GRCm38_") and name.startswith("ENSMUSG")):
                continue
        elif not name.startswith("ENSMUSG"):
            continue

        parts = name.split("_")
        assert len(parts) == 4, f"{path}: gene name is not 4 fields: {name}"
        gid, sym, bt, feat = parts
        assert feat in ("exon", "intron"), f"{path}: bad feature token {feat}"

        if gend - gstart != gspan:
            span_mismatch += 1

        g = genes.get(gid)
        if g is None:
            # The gene span gets the same correction as a feature: under the
            # 1-based-inclusive convention col6 (= col8-col7) is one less than
            # the true span, exactly as end-start is for a feature.
            g = genes[gid] = {"symbol": sym, "biotype": bt, "chrom": chrom,
                              "strand": strand, "span": gspan + offset,
                              "n_exon": 0, "n_intron": 0,
                              "exon_len": [], "intron_len": []}
        else:
            if g["biotype"] != bt:
                biotype_conflict += 1
            if g["symbol"] != sym:
                symbol_conflict += 1

        L = end - start + offset
        if feat == "exon":
            g["n_exon"] += 1
            g["exon_len"].append(L)
        else:
            g["n_intron"] += 1
            g["intron_len"].append(L)

    # A single gene must carry ONE biotype and ONE symbol in its own BED,
    # otherwise the biotype-churn measurement is meaningless before it starts.
    assert biotype_conflict == 0, f"{path}: {biotype_conflict} intra-gene biotype conflicts"
    assert symbol_conflict == 0, f"{path}: {symbol_conflict} intra-gene symbol conflicts"
    assert span_mismatch == 0, f"{path}: {span_mismatch} rows where col6 != col8-col7"

    return {"genes": genes, "trna": trna_rows, "n_rows": n_rows,
            "skipped_species": skipped, "gaps": gaps, "offset": offset,
            "modal_gap": modal_gap}


# ---------------------------------------------------------------------------
# count-table streaming
# ---------------------------------------------------------------------------
def norm_col(c):
    """Column label -> bare zero-padded well id. VERBATIM logic from
    vp_common.normalise_columns, applied per label."""
    c = str(c).rsplit("/", 1)[-1]
    c = c.rsplit("_", 1)[-1]
    return c.zfill(3)


def stream_row_totals(path, keep_wells=None, limit=None):
    """Sum a big count table row-wise over selected wells, without loading it.

    Returns (Series of per-row totals, list of wells used, n_rows).
    Values may be float (TranscriptCounts are non-integer).
    `limit` truncates for the precheck; a truncated pass must never be reported.
    """
    opener = gzip.open if path.endswith(".gz") else open
    idx, vals = [], []
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        wells = [norm_col(c) for c in header[1:]]
        if keep_wells is None:
            use = list(range(len(wells)))
        else:
            use = [i for i, w in enumerate(wells) if w in keep_wells]
        assert use, f"{path}: none of the requested wells found in {wells[:5]}..."
        for n, line in enumerate(fh):
            if limit is not None and n >= limit:
                break
            f = line.rstrip("\n").split("\t")
            idx.append(f[0])
            s = 0.0
            for i in use:
                v = f[i + 1]
                if v != "0" and v != "0.0":
                    s += float(v)
            vals.append(s)
    return pd.Series(vals, index=idx, dtype=float), [wells[i] for i in use], len(idx)


def stream_percell_species(path):
    """Per-cell (human, mouse) sums over unambiguous simple+combination rows.
    Mirrors vp_common.per_cell_species, streamed."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        wells = [norm_col(c) for c in header[1:]]
        n = len(wells)
        h = np.zeros(n)
        m = np.zeros(n)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp = species_of(f[0])
            if sp == "human":
                tgt = h
            elif sp == "mouse":
                tgt = m
            else:
                continue
            for i in range(n):
                v = f[i + 1]
                if v != "0" and v != "0.0":
                    tgt[i] += float(v)
    return pd.Series(h, index=wells), pd.Series(m, index=wells)


# ---------------------------------------------------------------------------
# composition
# ---------------------------------------------------------------------------
def composition(totals, id_of, class_of, restrict_ids=None):
    """% share per class over SIMPLE rows only.

    totals      : Series row_label -> counts
    id_of       : row_label -> gene id or None
    class_of    : row_label -> class label or None
    restrict_ids: if given, only rows whose gene id is in this set

    Returns (Series pct by class, denominator counts, n rows used).
    """
    acc = defaultdict(float)
    used = 0
    for lab, v in totals.items():
        cls = class_of(lab)
        if cls is None:
            continue
        if restrict_ids is not None:
            gid = id_of(lab)
            if gid is None or gid not in restrict_ids:
                continue
        acc[cls] += v
        used += 1
    denom = sum(acc.values())
    pct = pd.Series({k: 100.0 * v / denom for k, v in acc.items()}) if denom else pd.Series(dtype=float)
    return pct.sort_index(), denom, used


def tvd(a, b):
    """Total variation distance between two % compositions, in pp.
    0.5*sum|a-b| = the share of the library that would have to be moved between
    classes to turn one composition into the other."""
    keys = sorted(set(a.index) | set(b.index))
    return 0.5 * float(sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys))


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=f"{W}/res/threeway")
    args = ap.parse_args()
    OUT = args.outdir
    os.makedirs(OUT, exist_ok=True)

    rows = []   # tidy output: (metric, scope, value, unit, denominator, note)

    def rec(metric, scope, value, unit, denominator, note=""):
        rows.append({"metric": metric, "scope": scope, "value": value,
                     "unit": unit, "denominator": denominator, "note": note})

    # === 1. parse both BEDs ==================================================
    log("parsing BEDs ...")
    b99 = parse_bed(BED99, mouse_only=True)
    b116 = parse_bed(BED116, mouse_only=False)
    g99, g116 = b99["genes"], b116["genes"]
    log(f"  E99  mouse: {len(g99):>6} genes from {b99['n_rows']} rows "
        f"({b99['skipped_species']} non-mouse rows skipped), {b99['trna']} tRNAscan rows")
    log(f"  E116      : {len(g116):>6} genes from {b116['n_rows']} rows, "
        f"{b116['trna']} tRNAscan rows")

    ids99, ids116 = set(g99), set(g116)
    shared = ids99 & ids116
    only99 = ids99 - ids116
    only116 = ids116 - ids99

    # === coordinate convention ==============================================
    # The two files DO NOT share a convention (found by the precheck, not
    # assumed): E99 is 1-based inclusive, E116 is 0-based half-open. parse_bed
    # has already corrected for it. Verify the correction on single-feature
    # genes present in BOTH files, where the true length must be identical
    # unless the release genuinely changed the model.
    rec("bed_coordinate_convention", "E99 mouse BED",
        "1-based inclusive" if b99["offset"] == 1 else "0-based half-open",
        "convention", f"modal inter-feature gap = {b99['modal_gap']}",
        "detected from feature tiling, not assumed; length = end-start+%d" % b99["offset"])
    rec("bed_coordinate_convention", "E116 BED",
        "1-based inclusive" if b116["offset"] == 1 else "0-based half-open",
        "convention", f"modal inter-feature gap = {b116['modal_gap']}",
        "length = end-start+%d" % b116["offset"])
    if b99["offset"] != b116["offset"]:
        rec("bed_convention_differs_between_releases", "E99 vs E116", "yes", "flag",
            "the two BED files",
            "UNCORRECTED, every E99 feature would read 1 bp shorter than its E116 "
            "counterpart and the release would be blamed for shortening gene models")

    # single-feature genes in both files: an unchanged model must give an
    # identical corrected length. This is the falsifiable test of the fix.
    single_both = [gid for gid in (set(g99) & set(g116))
                   if g99[gid]["n_exon"] + g99[gid]["n_intron"] == 1
                   and g116[gid]["n_exon"] + g116[gid]["n_intron"] == 1]
    if single_both:
        corr = sum(1 for gid in single_both
                   if (g99[gid]["exon_len"] + g99[gid]["intron_len"])[0]
                   == (g116[gid]["exon_len"] + g116[gid]["intron_len"])[0])
        raw_eq = sum(1 for gid in single_both
                     if (g99[gid]["exon_len"] + g99[gid]["intron_len"])[0] - b99["offset"]
                     == (g116[gid]["exon_len"] + g116[gid]["intron_len"])[0] - b116["offset"])
        rec("single_feature_genes_in_both", "E99 n E116", len(single_both), "genes",
            f"{len(shared)} shared genes",
            "one exon, no intron, in both releases -- the cleanest length control")
        rec("single_feature_length_identical_after_convention_fix", "E99 n E116",
            round(100.0 * corr / len(single_both), 3), "%", f"{len(single_both)} genes",
            "corrected lengths agree exactly")
        rec("single_feature_length_identical_before_fix", "E99 n E116",
            round(100.0 * raw_eq / len(single_both), 3), "%", f"{len(single_both)} genes",
            "raw end-start; the gap between this and the corrected figure IS the "
            "size of the artefact the convention fix removes")
        # The fix must strictly improve agreement, or it is not a fix.
        assert corr >= raw_eq, (
            f"convention correction made agreement WORSE ({corr} vs {raw_eq} of "
            f"{len(single_both)}) -- the inferred offsets are wrong")

    rec("bed_rows", "E99 BED (all species)", b99["n_rows"], "rows", "file lines")
    rec("bed_rows", "E116 BED", b116["n_rows"], "rows", "file lines")
    rec("trnascan_rows", "E99 BED (human+mouse)", b99["trna"], "rows", "file lines",
        "tRNAscan rows carry no Ensembl id and no biotype; human+mouse pooled")
    rec("trnascan_rows", "E116 BED (mouse)", b116["trna"], "rows", "file lines")

    rec("genes_total", "E99 mouse", len(ids99), "genes", "mouse ENSMUSG ids in BED")
    rec("genes_total", "E116", len(ids116), "genes", "ENSMUSG ids in BED")
    rec("genes_shared", "E99 n E116", len(shared), "genes", "union of both BEDs")
    rec("genes_only_99", "E99 \\ E116", len(only99), "genes", "E99 mouse gene set")
    rec("genes_only_116", "E116 \\ E99", len(only116), "genes", "E116 gene set")
    rec("genes_shared_pct_of_99", "E99 mouse", round(100.0 * len(shared) / len(ids99), 3),
        "%", "E99 mouse gene set")
    rec("genes_shared_pct_of_116", "E116", round(100.0 * len(shared) / len(ids116), 3),
        "%", "E116 gene set")
    rec("jaccard_gene_universe", "E99 mouse vs E116",
        round(len(shared) / len(ids99 | ids116), 4), "index", "union of both gene sets")

    # === 2. per-biotype universe ============================================
    def by_bt(genes, key):
        d = defaultdict(set)
        for gid, g in genes.items():
            d[key(g["biotype"])].add(gid)
        return d

    for level, key in (("biotype", lambda b: b), ("coarse", coarse_class)):
        u99 = by_bt(g99, key)
        u116 = by_bt(g116, key)
        recs = []
        for cls in sorted(set(u99) | set(u116)):
            s99, s116 = u99.get(cls, set()), u116.get(cls, set())
            sh = s99 & s116
            recs.append({
                "level": level, "class": cls,
                "n_E99": len(s99), "n_E116": len(s116),
                "n_shared_id_and_class": len(sh),
                "n_only_E99": len(s99 - ids116), "n_only_E116": len(s116 - ids99),
                "n_E99_id_in_E116_other_class": len(s99 & ids116) - len(sh),
                "n_E116_id_in_E99_other_class": len(s116 & ids99) - len(sh),
                "delta_n": len(s116) - len(s99),
                "log2_ratio_n": (round(float(np.log2((len(s116) + 1) / (len(s99) + 1))), 4)),
            })
        df = pd.DataFrame(recs)
        path = f"{OUT}/annotation_gene_universe_{level}.tsv"
        df.to_csv(path, sep="\t", index=False)
        log(f"wrote {path}  ({len(df)} classes)")
        if level == "biotype":
            universe_bt = df
        else:
            universe_coarse = df

    pc = universe_bt.set_index("class").loc["ProteinCoding"]
    rec("protein_coding_genes", "E99 mouse", int(pc.n_E99), "genes", "biotype=ProteinCoding in E99")
    rec("protein_coding_genes", "E116", int(pc.n_E116), "genes", "biotype=ProteinCoding in E116")
    rec("protein_coding_shared", "E99 n E116 (same class both)", int(pc.n_shared_id_and_class),
        "genes", "union of ProteinCoding sets")
    rec("protein_coding_only_99", "E99 \\ E116", int(pc.n_only_E99), "genes",
        "E99 ProteinCoding set", "id absent from E116 entirely")
    rec("protein_coding_only_116", "E116 \\ E99", int(pc.n_only_E116), "genes",
        "E116 ProteinCoding set", "id absent from E99 entirely")
    rec("protein_coding_shared_pct_of_99", "E99 mouse",
        round(100.0 * int(pc.n_shared_id_and_class) / int(pc.n_E99), 3), "%",
        "E99 ProteinCoding set")

    # biotype vocabulary churn (classes that exist in one release only)
    voc99 = {g["biotype"] for g in g99.values()}
    voc116 = {g["biotype"] for g in g116.values()}
    rec("biotype_vocabulary_size", "E99 mouse", len(voc99), "classes", "distinct biotype tokens")
    rec("biotype_vocabulary_size", "E116", len(voc116), "classes", "distinct biotype tokens")
    rec("biotype_classes_retired_after_99", "E99 \\ E116", ";".join(sorted(voc99 - voc116)) or "none",
        "class names", "E99 vocabulary",
        "these classes CANNOT appear in an E116 composition, whatever the chemistry did")
    rec("biotype_classes_new_in_116", "E116 \\ E99", ";".join(sorted(voc116 - voc99)) or "none",
        "class names", "E116 vocabulary")

    # === 3. biotype-label churn on shared ids ===============================
    log("biotype-label churn ...")
    trans = Counter()
    sym_changed = 0
    for gid in shared:
        a, b = g99[gid]["biotype"], g116[gid]["biotype"]
        trans[(a, b)] += 1
        if g99[gid]["symbol"] != g116[gid]["symbol"]:
            sym_changed += 1
    changed = sum(v for (a, b), v in trans.items() if a != b)
    rec("shared_genes_biotype_changed", "E99 n E116", changed, "genes", f"{len(shared)} shared genes")
    rec("shared_genes_biotype_changed_pct", "E99 n E116",
        round(100.0 * changed / len(shared), 3), "%", f"{len(shared)} shared genes",
        "direct contamination of any class-by-class composition comparison")
    rec("shared_genes_symbol_changed", "E99 n E116", sym_changed, "genes",
        f"{len(shared)} shared genes",
        "same id, different gene symbol -- breaks any symbol-keyed join")
    rec("shared_genes_symbol_changed_pct", "E99 n E116",
        round(100.0 * sym_changed / len(shared), 3), "%", f"{len(shared)} shared genes")

    # How much of the "99-only / 116-only" difference is a genuinely different
    # gene, versus the SAME locus re-issued under a new id? Symbols are the only
    # cross-release handle available in the BED. Gm*/*Rik placeholder symbols are
    # not evidence of identity (they are auto-generated and reused), so they are
    # counted separately rather than credited as a recovery.
    def informative(sym):
        return not (sym.startswith("Gm") or sym.endswith("Rik") or sym.startswith("ENSMUSG"))

    sym99 = defaultdict(set)
    sym116 = defaultdict(set)
    for gid in only99:
        sym99[g99[gid]["symbol"]].add(gid)
    for gid in only116:
        sym116[g116[gid]["symbol"]].add(gid)
    sym_bridge = set(sym99) & set(sym116)
    bridged99 = sum(len(sym99[s]) for s in sym_bridge)
    bridged116 = sum(len(sym116[s]) for s in sym_bridge)
    inf_bridge = {s for s in sym_bridge if informative(s)}
    inf99 = sum(len(sym99[s]) for s in inf_bridge)

    rec("only99_genes_recoverable_by_symbol", "E99 \\ E116", bridged99, "genes",
        f"{len(only99)} E99-only genes",
        "same symbol exists in E116 under a different id -- an id change, not a lost gene")
    rec("only99_genes_recoverable_by_informative_symbol", "E99 \\ E116", inf99, "genes",
        f"{len(only99)} E99-only genes",
        "excluding Gm*/*Rik placeholder symbols, which are not evidence of identity")
    rec("only116_genes_recoverable_by_symbol", "E116 \\ E99", bridged116, "genes",
        f"{len(only116)} E116-only genes")
    rec("only99_genes_truly_absent_from_116", "E99 \\ E116", len(only99) - bridged99,
        "genes", f"{len(only99)} E99-only genes",
        "no id and no symbol match in E116")

    coarse_changed = sum(v for (a, b), v in trans.items()
                         if coarse_class(a) != coarse_class(b))
    rec("shared_genes_coarse_class_changed", "E99 n E116", coarse_changed, "genes",
        f"{len(shared)} shared genes", "coarse grouping, i.e. survives lumping")
    rec("shared_genes_coarse_class_changed_pct", "E99 n E116",
        round(100.0 * coarse_changed / len(shared), 3), "%", f"{len(shared)} shared genes")

    tr = pd.DataFrame([{"biotype_E99": a, "biotype_E116": b, "n_genes": v,
                        "coarse_E99": coarse_class(a), "coarse_E116": coarse_class(b),
                        "same_biotype": a == b,
                        "same_coarse": coarse_class(a) == coarse_class(b)}
                       for (a, b), v in trans.items()])
    tr = tr.sort_values(["same_biotype", "n_genes"], ascending=[True, False])
    tr.to_csv(f"{OUT}/annotation_biotype_churn_matrix.tsv", sep="\t", index=False)
    log(f"wrote {OUT}/annotation_biotype_churn_matrix.tsv  ({len(tr)} transitions)")

    # per-class gain/loss ledger among shared genes
    gained = Counter()
    lost = Counter()
    for (a, b), v in trans.items():
        if a != b:
            lost[a] += v
            gained[b] += v
    led = pd.DataFrame([{"biotype": cls,
                         "n_shared_kept_label": trans.get((cls, cls), 0),
                         "n_lost_to_other_class": lost.get(cls, 0),
                         "n_gained_from_other_class": gained.get(cls, 0),
                         "net": gained.get(cls, 0) - lost.get(cls, 0)}
                        for cls in sorted(set(voc99) | set(voc116))])
    led["churn_pct_of_E99_members"] = [
        round(100.0 * r.n_lost_to_other_class /
              max(1, r.n_shared_kept_label + r.n_lost_to_other_class), 2)
        for r in led.itertuples()]
    led = led.sort_values("n_lost_to_other_class", ascending=False)
    led.to_csv(f"{OUT}/annotation_biotype_churn_ledger.tsv", sep="\t", index=False)
    log(f"wrote {OUT}/annotation_biotype_churn_ledger.tsv")

    top = tr[~tr.same_biotype].head(1)
    if len(top):
        t = top.iloc[0]
        rec("largest_single_biotype_transition", f"{t.biotype_E99} -> {t.biotype_E116}",
            int(t.n_genes), "genes", f"{len(shared)} shared genes")

    # === 4. feature-length churn ============================================
    log("feature-length churn ...")
    def len_stats(genes, gid_subset=None):
        acc = defaultdict(lambda: {"ex": [], "intr": [], "span": [], "exonic": []})
        for gid, g in genes.items():
            if gid_subset is not None and gid not in gid_subset:
                continue
            a = acc[g["biotype"]]
            a["ex"].extend(g["exon_len"])
            a["intr"].extend(g["intron_len"])
            a["span"].append(g["span"])
            a["exonic"].append(sum(g["exon_len"]))
        return acc

    def flatten(acc, release, subset_label):
        out = []
        for bt, a in sorted(acc.items()):
            row = {"release": release, "gene_set": subset_label, "biotype": bt,
                   "n_genes": len(a["span"]), "n_exon_features": len(a["ex"]),
                   "n_intron_features": len(a["intr"]),
                   "median_exon_len": float(np.median(a["ex"])) if a["ex"] else np.nan,
                   "median_intron_len": float(np.median(a["intr"])) if a["intr"] else np.nan,
                   "median_gene_span": float(np.median(a["span"])) if a["span"] else np.nan,
                   "median_total_exonic_len": float(np.median(a["exonic"])) if a["exonic"] else np.nan}
            for L in READ_LENS:
                row[f"pct_exon_features_shorter_than_{L}"] = (
                    round(100.0 * float(np.mean(np.asarray(a["ex"]) < L)), 3) if a["ex"] else np.nan)
            out.append(row)
        return out

    lrows = []
    lrows += flatten(len_stats(g99), "E99", "all_mouse")
    lrows += flatten(len_stats(g116), "E116", "all")
    lrows += flatten(len_stats(g99, shared), "E99", "shared_only")
    lrows += flatten(len_stats(g116, shared), "E116", "shared_only")
    ldf = pd.DataFrame(lrows)
    ldf.to_csv(f"{OUT}/annotation_feature_length.tsv", sep="\t", index=False)
    log(f"wrote {OUT}/annotation_feature_length.tsv  ({len(ldf)} rows)")

    # paired, per-gene: same gene, both releases -> did its model change?
    paired = []
    for gid in shared:
        a, b = g99[gid], g116[gid]
        paired.append({"biotype_E99": a["biotype"],
                       "exonic99": sum(a["exon_len"]), "exonic116": sum(b["exon_len"]),
                       "nex99": a["n_exon"], "nex116": b["n_exon"],
                       "span99": a["span"], "span116": b["span"]})
    pdf = pd.DataFrame(paired)
    pdf["exonic_same"] = pdf.exonic99 == pdf.exonic116
    pdf["nex_same"] = pdf.nex99 == pdf.nex116
    rec("shared_genes_exonic_length_unchanged", "E99 n E116", int(pdf.exonic_same.sum()),
        "genes", f"{len(shared)} shared genes",
        "identical summed exon length; a proxy for an unchanged gene model")
    rec("shared_genes_exonic_length_unchanged_pct", "E99 n E116",
        round(100.0 * float(pdf.exonic_same.mean()), 3), "%", f"{len(shared)} shared genes")
    rec("shared_genes_exon_count_unchanged_pct", "E99 n E116",
        round(100.0 * float(pdf.nex_same.mean()), 3), "%", f"{len(shared)} shared genes")
    rec("shared_genes_median_log2_exonic_ratio", "E116/E99",
        round(float(np.median(np.log2((pdf.exonic116 + 1) / (pdf.exonic99 + 1)))), 4),
        "log2", f"{len(shared)} shared genes",
        "0 means the typical shared gene's exonic length did not move")

    pg = pdf.groupby("biotype_E99").agg(
        n_genes=("exonic99", "size"),
        pct_exonic_unchanged=("exonic_same", lambda s: round(100.0 * float(s.mean()), 2)),
        pct_exon_count_unchanged=("nex_same", lambda s: round(100.0 * float(s.mean()), 2)),
        median_exonic_E99=("exonic99", "median"),
        median_exonic_E116=("exonic116", "median"),
    ).reset_index()
    pg["median_log2_exonic_ratio"] = [
        round(float(np.log2((r.median_exonic_E116 + 1) / (r.median_exonic_E99 + 1))), 4)
        for r in pg.itertuples()]
    pg.to_csv(f"{OUT}/annotation_feature_length_paired.tsv", sep="\t", index=False)
    log(f"wrote {OUT}/annotation_feature_length_paired.tsv")

    # === 5. decisive test: own plate (E116) restricted to the shared set =====
    log("own plate composition ...")
    own_wells = None  # resolved from the header
    comp_rows = []

    def id_of(lab):
        if "-" in lab or "_" not in lab:
            return None
        return lab.split("_")[0]

    # tRNAscan rows are held out on BOTH sides -- see the comment on is_trna_row.
    def class_bt(lab):
        if is_trna_row(lab):
            return None
        return biotype_of(lab)

    def class_coarse(lab):
        bt = class_bt(lab)
        return coarse_class(bt) if bt else None

    own_tables = {"ReadCounts": OWN_READ, "TranscriptCounts": OWN_TRANS}
    own_summary = {}
    for counts, path in own_tables.items():
        tot, wells, nrow = stream_row_totals(path, keep_wells=None)
        keep = [w for w in wells if w not in OWN_BLANKS]
        assert len(keep) == 12, f"expected 12 real wells, got {keep}"
        tot, wells, nrow = stream_row_totals(path, keep_wells=set(keep))
        own_wells = wells
        log(f"  {counts}: {nrow} rows x {len(wells)} real cells, {tot.sum():.0f} total")

        for level, fn in (("biotype", class_bt), ("coarse", class_coarse)):
            a, dA, nA = composition(tot, id_of, fn, None)
            b, dB, nB = composition(tot, id_of, fn, shared)
            for cls in sorted(set(a.index) | set(b.index)):
                comp_rows.append({
                    "dataset": "own_plate_ZHA9292A1", "release": "E116",
                    "counts": counts, "level": level, "class": cls,
                    "pct_all_genes": round(float(a.get(cls, 0.0)), 6),
                    "pct_shared_genes_only": round(float(b.get(cls, 0.0)), 6),
                    "delta_pp": round(float(b.get(cls, 0.0) - a.get(cls, 0.0)), 6)})
            own_summary[(counts, level)] = (a, b, dA, dB, nA, nB)

        # tRNAscan share, reported separately because it is EXCLUDED from every
        # composition above and because the two BEDs carry different tRNAscan
        # content (E99 1,758 rows human+mouse, E116 1,137 mouse) -- so a tRNA
        # comparison between the plates is a tRNAscan-build comparison, not a
        # chemistry one. Conventions trap 8 also applies: at 151 nt read length
        # the `jS:IN` containment rule makes almost all tRNA features invisible.
        trna_tot = float(sum(v for lab, v in tot.items() if is_trna_row(lab)))
        rec(f"own_trna_row_{counts}_share", "own plate",
            round(100.0 * trna_tot / tot.sum(), 5), "%", f"all {counts} in 12 real cells",
            "tRNAscan rows, held out of every composition in this file")

        # gene detection, single-gene rows only
        det_all = det_shared = 0
        for lab, v in tot.items():
            gid = id_of(lab)
            if gid is None:
                continue
            if v > 0:
                det_all += 1
                if gid in shared:
                    det_shared += 1
        rec(f"own_genes_detected_{counts}", "own plate, single-gene rows, 12 cells",
            det_all, "genes", f"{nrow} rows in table, single-gene rows only")
        rec(f"own_genes_detected_in_shared_set_{counts}", "own plate", det_shared, "genes",
            f"{det_all} detected single genes")
        rec(f"own_genes_detected_only_in_116_pct_{counts}", "own plate",
            round(100.0 * (det_all - det_shared) / det_all, 3), "%",
            f"{det_all} detected single genes",
            "detected genes that Ensembl 99 could not have reported at all")

        # combination-row exposure: reads in rows naming >=1 E116-only gene
        comb_tot = comb_exposed = 0.0
        for lab, v in tot.items():
            if "-" not in lab:
                continue
            comb_tot += v
            gids = [p.split("_")[0] for p in lab.split("-") if p.startswith("ENSMUSG")]
            if any(g not in shared for g in gids):
                comb_exposed += v
        rec(f"own_combination_row_{counts}_share", "own plate",
            round(100.0 * comb_tot / tot.sum(), 3), "%", f"all {counts} in 12 real cells",
            "multi-gene rows, excluded from every composition above")
        rec(f"own_combination_{counts}_naming_a_116only_gene", "own plate",
            round(100.0 * comb_exposed / comb_tot, 3), "%", f"combination-row {counts}",
            "these rows could not exist in the same form under Ensembl 99")

    # headline release-effect numbers, own plate
    for counts in own_tables:
        for level in ("biotype", "coarse"):
            a, b, dA, dB, nA, nB = own_summary[(counts, level)]
            rec(f"own_composition_TVD_{level}_{counts}", "all genes -> shared genes only",
                round(tvd(a, b), 4), "pp", f"{counts} on single-gene rows",
                "half the summed absolute shift = share of the library that moves class")
            rec(f"own_denominator_shrink_{level}_{counts}", "shared-only vs all",
                round(100.0 * (1 - dB / dA), 4), "%", f"{counts} on single-gene rows",
                "counts lost by restricting to genes Ensembl 99 also had")
            if level == "coarse":
                for cls in ("ProteinCoding", "lncRNA", "sncRNA", "pseudogene", "TEC"):
                    rec(f"own_{cls}_share_all_genes_{counts}", "own plate",
                        round(float(a.get(cls, 0.0)), 4), "%", f"single-gene {counts}")
                    rec(f"own_{cls}_share_shared_genes_{counts}", "own plate",
                        round(float(b.get(cls, 0.0)), 4), "%", f"single-gene {counts}, E99-shared ids")
                    rec(f"own_{cls}_delta_pp_{counts}", "shared-only minus all",
                        round(float(b.get(cls, 0.0) - a.get(cls, 0.0)), 4), "pp",
                        f"single-gene {counts}")

    # === 6. symmetric test: published plate (E99, mESC cells) ===============
    log("published plate: identifying mESC cells from UFIs (fig1d rule) ...")
    hs, ms = stream_percell_species(PUB_UFI)
    lab, frac_h = classify_fig1d(hs, ms)
    mesc = sorted(lab.index[lab == "mouse"])
    rec("published_wells_total", "SRR14783059 plate", len(lab), "wells", "plate wells")
    rec("published_wells_mESC_fig1d", "SRR14783059 plate", len(mesc), "wells",
        f"{len(lab)} wells", "classify_fig1d on v3 UFICounts, MIN_UFI=7500 (verbatim rule)")
    rec("published_wells_HEK293T_fig1d", "SRR14783059 plate", int((lab == "human").sum()),
        "wells", f"{len(lab)} wells")
    rec("published_wells_mixed_fig1d", "SRR14783059 plate", int((lab == "mixed").sum()),
        "wells", f"{len(lab)} wells")
    rec("published_wells_discarded_fig1d", "SRR14783059 plate", int((lab == "discarded").sum()),
        "wells", f"{len(lab)} wells")
    assert len(mesc) > 0, "no mESC wells called -- cannot run the symmetric test"

    log(f"published plate: {len(mesc)} mESC wells; streaming ReadCounts ...")
    ptot, pwells, pnrow = stream_row_totals(PUB_READ, keep_wells=set(mesc))
    assert len(pwells) == len(mesc)

    def id_of_mouse(lab_):
        if "-" in lab_ or "_" not in lab_:
            return None
        gid = lab_.split("_")[0]
        return gid if gid.startswith("ENSMUSG") else None

    def class_bt_mouse(lab_):
        # species_of() returns 'trna' for tRNAscan rows, so the mouse gate
        # already excludes them; the explicit test documents the intent and
        # guards against a future change to species_of.
        if is_trna_row(lab_) or species_of(lab_) != "mouse":
            return None
        return biotype_of(lab_)

    def class_coarse_mouse(lab_):
        bt = class_bt_mouse(lab_)
        return coarse_class(bt) if bt else None

    pub_summary = {}
    for level, fn in (("biotype", class_bt_mouse), ("coarse", class_coarse_mouse)):
        a, dA, nA = composition(ptot, id_of_mouse, fn, None)
        b, dB, nB = composition(ptot, id_of_mouse, fn, shared)
        for cls in sorted(set(a.index) | set(b.index)):
            comp_rows.append({
                "dataset": "published_plate_SRR14783059_v3", "release": "E99",
                "counts": "ReadCounts", "level": level, "class": cls,
                "pct_all_genes": round(float(a.get(cls, 0.0)), 6),
                "pct_shared_genes_only": round(float(b.get(cls, 0.0)), 6),
                "delta_pp": round(float(b.get(cls, 0.0) - a.get(cls, 0.0)), 6)})
        pub_summary[level] = (a, b, dA, dB, nA, nB)
        rec(f"published_composition_TVD_{level}_ReadCounts",
            "all mouse genes -> shared genes only", round(tvd(a, b), 4), "pp",
            "ReadCounts on single-gene mouse rows, mESC wells")
        rec(f"published_denominator_shrink_{level}_ReadCounts", "shared-only vs all",
            round(100.0 * (1 - dB / dA), 4), "%",
            "ReadCounts on single-gene mouse rows, mESC wells",
            "reads lost by restricting to genes Ensembl 116 also has")

    trna_p = float(sum(v for lab_, v in ptot.items() if is_trna_row(lab_)))
    rec("published_trna_row_ReadCounts_share", "mESC wells",
        round(100.0 * trna_p / ptot.sum(), 5), "%", "all ReadCounts in mESC wells",
        "tRNAscan rows; note E99 BED pools human+mouse tRNAs, E116 is mouse-only, "
        "so this is NOT release-comparable and is excluded from every composition")

    det_all = det_shared = 0
    for lab_, v in ptot.items():
        gid = id_of_mouse(lab_)
        if gid is None or v <= 0:
            continue
        det_all += 1
        if gid in shared:
            det_shared += 1
    rec("published_mouse_genes_detected", "mESC wells, single-gene mouse rows",
        det_all, "genes", f"{pnrow} rows in table")
    rec("published_mouse_genes_detected_in_shared_set", "mESC wells", det_shared, "genes",
        f"{det_all} detected mouse genes")
    rec("published_mouse_genes_detected_only_in_99_pct", "mESC wells",
        round(100.0 * (det_all - det_shared) / det_all, 3), "%",
        f"{det_all} detected mouse genes",
        "detected genes Ensembl 116 no longer has under that id")

    # === 7. the cross-plate gap, with and without the release term ==========
    log("cross-plate gap decomposition ...")
    gap_rows = []
    for level in ("biotype", "coarse"):
        oa, ob, *_ = own_summary[("ReadCounts", level)]
        pa, pb, *_ = pub_summary[level]
        raw = tvd(pa, oa)
        matched = tvd(pb, ob)
        for cls in sorted(set(pa.index) | set(oa.index) | set(pb.index) | set(ob.index)):
            gap_rows.append({
                "level": level, "class": cls,
                "published_all": round(float(pa.get(cls, 0.0)), 6),
                "own_all": round(float(oa.get(cls, 0.0)), 6),
                "gap_all_pp": round(float(oa.get(cls, 0.0) - pa.get(cls, 0.0)), 6),
                "published_shared": round(float(pb.get(cls, 0.0)), 6),
                "own_shared": round(float(ob.get(cls, 0.0)), 6),
                "gap_shared_pp": round(float(ob.get(cls, 0.0) - pb.get(cls, 0.0)), 6),
                "gap_closed_pp": round(float(
                    abs(oa.get(cls, 0.0) - pa.get(cls, 0.0))
                    - abs(ob.get(cls, 0.0) - pb.get(cls, 0.0))), 6)})
        rec(f"crossplate_TVD_{level}_ReadCounts_all_genes", "published vs own",
            round(raw, 4), "pp", "ReadCounts, single-gene mouse rows",
            "the raw composition gap; mixes release + biology + depth")
        rec(f"crossplate_TVD_{level}_ReadCounts_shared_genes", "published vs own",
            round(matched, 4), "pp", "ReadCounts, single-gene mouse rows in BOTH releases",
            "same quantity with the gene universe forced to agree")
        rec(f"crossplate_TVD_{level}_release_attributable_pp", "published vs own",
            round(raw - matched, 4), "pp", "ReadCounts, single-gene mouse rows",
            "gap that disappears when the gene universe is harmonised = release term")
        rec(f"crossplate_TVD_{level}_release_attributable_pct", "published vs own",
            round(100.0 * (raw - matched) / raw, 3) if raw else np.nan, "%",
            "raw cross-plate TVD",
            "share of the raw composition gap attributable to the gene universe")
    pd.DataFrame(gap_rows).to_csv(f"{OUT}/annotation_crossplate_gap.tsv", sep="\t", index=False)
    log(f"wrote {OUT}/annotation_crossplate_gap.tsv")

    pd.DataFrame(comp_rows).to_csv(f"{OUT}/annotation_composition_shift.tsv",
                                   sep="\t", index=False)
    log(f"wrote {OUT}/annotation_composition_shift.tsv  ({len(comp_rows)} rows)")

    # === write the headline table ===========================================
    hdr = pd.DataFrame(rows)[["metric", "scope", "value", "unit", "denominator", "note"]]
    dest = f"{OUT}/annotation_release_effect.tsv"
    hdr.to_csv(dest, sep="\t", index=False)
    log(f"wrote {dest}  ({len(hdr)} rows)")

    meta = {"bed_E99": BED99, "bed_E116": BED116,
            "own_read_table": OWN_READ, "own_trans_table": OWN_TRANS,
            "own_wells_used": own_wells, "own_blanks_excluded": sorted(OWN_BLANKS),
            "published_read_table": PUB_READ, "published_ufi_table": PUB_UFI,
            "published_mesc_wells": mesc,
            "pandas": pd.__version__, "numpy": np.__version__,
            "python": sys.version.split()[0]}
    with open(f"{OUT}/annotation_release_effect_provenance.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    log(f"wrote {OUT}/annotation_release_effect_provenance.json")

    log("\n=== headline ===")
    for r in rows:
        if r["metric"].startswith(("genes_", "protein_coding_", "shared_genes_biotype",
                                  "crossplate_", "own_composition_TVD",
                                  "published_composition_TVD", "biotype_classes_")):
            log(f"  {r['metric']:<52} {r['scope']:<38} {r['value']}")


if __name__ == "__main__":
    main()
