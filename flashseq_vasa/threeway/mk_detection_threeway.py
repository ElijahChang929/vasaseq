#!/usr/bin/env python3
"""Three-way gene detection and saturation, depth-matched:
   published VASA-plate  vs  own VASA-plate  vs  FLASH-seq (both arms).

WHAT THIS EXTENDS
-----------------
res/flashseq_vasa/mk_gene_detection.py already did own-plate vs FLASH-seq. This
script adds the published VASA-plate (SRR14783059, run `vasaplate_out_v3`) as a
third arm. The estimator is REUSED VERBATIM from that script, and its
faithfulness is asserted against res/flashseq_vasa/gene_detection_matched.tsv
(see `self-check` below) rather than assumed.

    E[genes detected at depth d] = sum_g [ 1 - (1-p)^{c_g} ],  p = d / total_reads

exact in expectation, deterministic, no seed.

A TABLE-TYPE MISMATCH IN THE COMPARISON BEING EXTENDED
------------------------------------------------------
Writing the self-check surfaced this, and it is not cosmetic, so it is measured
here rather than inherited. countTables_fromPickle.py writes TWO different gene
tables:

  * `_total.ReadCounts.tsv`          (line 107, from `total_reads_df`)
        raw, PRE-aggregation. Every row is a literal annotation entry; every
        `geneA-geneB` combination is still its own row.
  * `_uniaggGenes_total.ReadCounts.tsv` (line 251, from `agg_cntdf_genes`)
        POST-aggregation. `reduceGeneName` collapses many combination entries
        onto ONE component gene, and `groupby('new_gene')` ADDS those reads into
        that gene's row. A single-gene row therefore carries strictly more reads
        here than in `_total`.

mk_gene_detection.py read the own plate from `_uniaggGenes_total` but FLASH-seq
from `_total`. So its VASA-minus-FLASH-seq lead mixes the protocol difference
with a table difference, in the direction that favours VASA. This script
therefore carries `table` as an explicit dimension ('uniagg' and 'total') for
all three datasets, puts the headline on `uniagg` everywhere (the task's
published-plate table, and the family the rest of the three-way work treats as
primary -- NOT, as an earlier version of this note claimed, the only family that
exists for all three: `_total.ReadCounts.tsv` exists for all three too, which is
what makes the `table` dimension below possible), and reports the size of the table effect per dataset in
detection_threeway_tableeffect.tsv. The self-check reproduces the reference
numbers using the reference's OWN table choice, which is what proves the
estimator identical without adopting the mismatch.

WHY vasaplate_out_v3 AND NOT bedv2
----------------------------------
v3 is the validated benchmark anchor: it carries the corrected rRNA reference
(one antisense entry reverse-complemented), which took the rRNA-biotype excess
from 623x published down to 3.2x. bedv2 still carries the ~600x artefact. See
code/I_Gene_expression/vasaplate_check/BENCHMARK_published_plate.md.

THREE THINGS THAT ARE NOT PROTOCOL DIFFERENCES
----------------------------------------------
1. ANNOTATION RELEASE. The published plate was mapped to GRCm38 + hg38 with
   Ensembl 99; the own plate and FLASH-seq to GRCm39 with Ensembl 116. Ensembl
   116 has more protein-coding mouse genes than 99, so the two GRCm39 datasets
   can detect genes that DO NOT EXIST in the 99 universe. The headline
   comparison is therefore restricted to the shared gene universe (protein-coding
   in BOTH releases, matched on the stable ENSMUSG accession); the unrestricted
   numbers are reported alongside, and the overlap counts are derived here from
   the two annotation BEDs (universe='all' vs universe='shared').
2. SPECIES. The published library is a HEK293T + mESC mixing control on a
   concatenated reference, so its rows are ENSG* (human) or ENSMUSG* (mouse) and
   its wells are human, mouse, mixed or discarded. Only mouse rows in mouse cells
   are used. Cells are classified with vp_common's two published rules, which
   disagree; both are carried per cell and the headline states which was used.
3. DEPTH. The published plate's cells are far shallower than the own plate's, so
   the deep rungs of the grid are supported by only some datasets. No curve is
   ever extrapolated past a unit's own native depth, and every track reports the
   maximum depth at which ALL its units still contribute.

DENOMINATORS (Rule 5)
---------------------
* 'depth' throughout means READS IN THE COUNTED SCOPE, i.e. reads on single-gene
  protein-coding entries of the universe in force -- NOT library reads, NOT STAR
  input, NOT mapped reads. This is what mk_gene_detection.py thins, and it is
  kept identical.
* 'genes' means single-gene protein-coding entries of that universe with >=1 read
  ('-' not in the entry name, biotype token == ProteinCoding).
* The '-' rule is upstream's and is kept for comparability, but it also drops
  real single genes whose SYMBOL contains a hyphen (Nkx2-5). That loss is
  quantified per dataset in detection_threeway_scope.tsv rather than silently
  accepted (Rule 3).

COUNT COLUMN (Rule 4)
---------------------
ReadCounts on all three sides. FLASH-seq is smartseq_noUMI, where every read is
filed under one literal UMI key, so UFICounts degenerates to a detection mask and
TranscriptCounts == UFICounts; reads are the only unit all three protocols
measure. Species classification of the published plate is the one exception: it
uses UFICounts, because the published >=7,500 gate is defined on UFIs.

USAGE
    mk_detection_threeway.py            # paths from the CONFIG block / env
    mk_detection_threeway.py --precheck # cheap assertions only, no big reads
"""
import os
import sys

import numpy as np
import pandas as pd

# --- CONFIG ----------------------------------------------------------------
W = os.environ.get("W", "/nemo/lab/turnerj/working/guangxin/vasaseq")
REF = os.environ.get("VASA_REF", os.path.join(W, "..", "reference", "vasaseq"))
FS_DIR = os.environ.get(
    "FS_DIR", "/nemo/lab/turnerj/scratch/zhangg/vasaseq/flashseq_vasa")
VP_DIR = os.path.join(W, "code", "I_Gene_expression", "vasaplate_check")

PUB_DIR = os.path.join(W, "data/ref/fastq_vasaplate")
PUB_UNIAGG = os.environ.get("PUB_UNIAGG", os.path.join(
    PUB_DIR, "vasaplate_out_v3_uniaggGenes_total.ReadCounts.tsv"))
PUB_TOTAL = os.environ.get("PUB_TOTAL", os.path.join(
    PUB_DIR, "vasaplate_out_v3_total.ReadCounts.tsv"))
PUB_UFI = os.environ.get("PUB_UFI", os.path.join(
    PUB_DIR, "vasaplate_out_v3_total.UFICounts.tsv"))
OWN = os.environ.get("OWN", os.path.join(
    W, "data/PM26037/out/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv"))
OWN_TOTAL = os.environ.get("OWN_TOTAL", os.path.join(
    W, "data/PM26037/out/ZHA9292A1_total.ReadCounts.tsv"))
META = os.environ.get("META", os.path.join(W, "code/flashseq/sample_metadata.tsv"))
BED99 = os.environ.get("BED99", os.path.join(
    REF, "mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed"))
BED116 = os.environ.get("BED116", os.path.join(
    REF, "mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed"))
REF_MATCHED = os.environ.get("REF_MATCHED", os.path.join(
    W, "res/flashseq_vasa/gene_detection_matched.tsv"))
OUT = os.environ.get("OUT", os.path.join(W, "res/threeway"))

BLANKS = {"001", "014", "015", "016"}          # own plate: empty wells
GRID = [1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7]
SUMMARY_RUNGS = [1e5, 2e5, 5e5, 1e6, 2e6, 5e6]
PRECHECK = "--precheck" in sys.argv

sys.path.insert(0, VP_DIR)
import vp_common as vp                                          # noqa: E402

print("versions: pandas %s numpy %s python %s"
      % (pd.__version__, np.__version__, sys.version.split()[0]))


# --- estimator (verbatim from res/flashseq_vasa/mk_gene_detection.py) ------
def thinned(counts, p):
    """E[# genes with >=1 read] after binomial thinning at rate p."""
    c = counts[counts > 0].astype(float)
    if p >= 1.0:
        return float(len(c))
    return float(np.sum(1.0 - np.power(1.0 - p, c)))


# --- annotation universes --------------------------------------------------
def bed_universe(path, mouse_only):
    """{ens_id: biotype} over gene entries of an annotation BED.

    Field 5 is `<ens>_<symbol>_<biotype>_<exon|intron>`; tRNA rows do not follow
    that contract and are not protein-coding, so they are counted and skipped.
    A gene id appearing with two biotypes would break the mapping, so that is
    asserted rather than assumed.
    """
    uni, sym, n_trna, n_odd, clash = {}, {}, 0, 0, 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name = f[4]
            if not (name.startswith("ENSG") or name.startswith("ENSMUSG")):
                n_trna += 1
                continue
            if mouse_only and not name.startswith("ENSMUSG"):
                continue
            parts = name.split("_")
            if len(parts) < 4 or parts[-1] not in ("exon", "intron"):
                n_odd += 1
                continue
            ens, bio = parts[0], parts[-2]
            if ens in uni and uni[ens] != bio:
                clash += 1
            uni[ens] = bio
            sym[ens] = "_".join(parts[1:-2])
    return uni, sym, dict(non_ens_rows=n_trna, off_contract_rows=n_odd,
                          biotype_clashes=clash)


# --- count-table loading ---------------------------------------------------
def read_counts(path, big=False):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    dt = {c: np.int32 for c in header[1:]} if big else None
    df = pd.read_csv(path, sep="\t", index_col=0, dtype=dt)
    assert df.shape[1] == len(header) - 1, "%s: width mismatch" % path
    return df


def scope(df, mouse_only=False):
    """Upstream's counting scope, unchanged: single-gene protein-coding entries.

    Returns (kept_df, diagnostics). `single` is upstream's '-'-free test; the
    ENS-token count is the stricter test and the gap between them is the price
    of the hyphen rule, reported not fixed.
    """
    idx = df.index.astype(str)
    bio = pd.Series([i.rsplit("_", 1)[-1] if "_" in i else "NA" for i in idx],
                    index=df.index)
    single = ~pd.Series([("-" in i) for i in idx], index=df.index)
    ntok = pd.Series([i.count("ENSMUSG") + i.count("ENSG") for i in idx],
                     index=df.index)
    pc = bio == "ProteinCoding"
    mouse = pd.Series([i.startswith("ENSMUSG") for i in idx], index=df.index)
    keep = single & pc
    if mouse_only:
        keep = keep & mouse
    # entries the hyphen rule loses although they name exactly one gene
    lost = (~single) & pc & (ntok == 1) & (mouse if mouse_only else True)
    diag = dict(rows_total=int(df.shape[0]),
                rows_protein_coding=int(pc.sum()),
                rows_in_scope=int(keep.sum()),
                rows_pc_single_gene_lost_to_hyphen_rule=int(lost.sum()),
                reads_lost_to_hyphen_rule=int(df[lost.values].values.sum()),
                reads_in_scope=int(df[keep.values].values.sum()))
    return df[keep.values], diag


def ens_of(index):
    return pd.Series([str(i).split("_", 1)[0] for i in index], index=index)


# ===========================================================================
# 1. annotation universes and the shared gene set
# ===========================================================================
def fs_path(arm, table):
    stem = "_uniaggGenes_total" if table == "uniagg" else "_total"
    return "%s/%s/FSall10_%s%s.ReadCounts.tsv" % (FS_DIR, arm, arm, stem)


INPUTS = [PUB_UNIAGG, PUB_TOTAL, PUB_UFI, OWN, OWN_TOTAL, META, BED99, BED116,
          REF_MATCHED] + [fs_path(a, t) for a in ("native", "vasalen")
                          for t in ("uniagg", "total")]
for p in INPUTS:
    assert os.path.exists(p) and os.path.getsize(p) > 0, "missing/empty: %s" % p
os.makedirs(OUT, exist_ok=True)

u99, s99, d99 = bed_universe(BED99, mouse_only=True)
u116, s116, d116 = bed_universe(BED116, mouse_only=False)
pc99 = {e for e, b in u99.items() if b == "ProteinCoding"}
pc116 = {e for e, b in u116.items() if b == "ProteinCoding"}
shared_pc = pc99 & pc116
shared_any = set(u99) & set(u116)

pc99_only = pc99 - pc116
pc116_only = pc116 - pc99
reassigned_99pc = {e: u116[e] for e in pc99_only if e in u116}
reassigned_116pc = {e: u99[e] for e in pc116_only if e in u99}

uni_rows = [
    ("mouse gene ids in Ensembl 99 BED (any biotype)", len(u99)),
    ("mouse gene ids in Ensembl 116 BED (any biotype)", len(u116)),
    ("mouse gene ids shared, any biotype", len(shared_any)),
    ("protein-coding mouse genes, Ensembl 99", len(pc99)),
    ("protein-coding mouse genes, Ensembl 116", len(pc116)),
    ("protein-coding in BOTH releases (shared universe)", len(shared_pc)),
    ("protein-coding in 99 only", len(pc99_only)),
    ("  ... of which absent from 116 altogether", len(pc99_only) - len(reassigned_99pc)),
    ("  ... of which present in 116 with another biotype", len(reassigned_99pc)),
    ("protein-coding in 116 only", len(pc116_only)),
    ("  ... of which absent from 99 altogether", len(pc116_only) - len(reassigned_116pc)),
    ("  ... of which present in 99 with another biotype", len(reassigned_116pc)),
]
for k, v in sorted(d99.items()):
    uni_rows.append(("BED99 parse: %s" % k, v))
for k, v in sorted(d116.items()):
    uni_rows.append(("BED116 parse: %s" % k, v))
pd.DataFrame(uni_rows, columns=["quantity", "n_genes"]).to_csv(
    "%s/annotation_universe.tsv" % OUT, sep="\t", index=False)
print("\n=== ANNOTATION UNIVERSE (mouse, from the two BEDs) ===")
for k, v in uni_rows[:12]:
    print("  %-52s %8d" % (k, v))
assert len(shared_pc) > 15000, "shared protein-coding universe implausibly small"

if PRECHECK:
    print("\n[precheck] annotation universes parsed; stopping before the big reads.")
    for p in INPUTS:
        if not p.endswith(".tsv"):
            continue
        with open(p) as fh:
            h = fh.readline().rstrip("\n").split("\t")
        print("  %-62s %4d columns  %8.1f MB"
              % (os.path.basename(p), len(h) - 1, os.path.getsize(p) / 1e6))
    sys.exit(0)

# ===========================================================================
# 2. published plate: species assignment (vp_common's two published rules)
# ===========================================================================
print("\nreading %s ..." % os.path.basename(PUB_UFI))
ufi = vp.normalise_columns(read_counts(PUB_UFI, big=True))
sp = vp.species_vector(ufi.index)
h_ufi, m_ufi = vp.per_cell_species(ufi, sp)
gh = (ufi[sp.values == "human"] > 0).sum()
gm = (ufi[sp.values == "mouse"] > 0).sum()
lab_f, frac_uf = vp.classify_fig1d(h_ufi, m_ufi)
lab_m, frac_gf = vp.classify_methods(gh, gm, h_ufi, m_ufi)
cells = pd.DataFrame(dict(well=ufi.columns, ufi_human=h_ufi.values,
                          ufi_mouse=m_ufi.values, genes_human=gh.values,
                          genes_mouse=gm.values,
                          frac_human_ufi=frac_uf.values.round(4),
                          frac_human_genes=frac_gf.values.round(4),
                          label_fig1d=lab_f.values, label_methods=lab_m.values))
cells.to_csv("%s/published_cell_species.tsv" % OUT, sep="\t", index=False)
print("=== PUBLISHED PLATE, 384 wells, species calls ===")
print("  Fig.1d rule (UFI purity):   %s" % dict(lab_f.value_counts()))
print("  Methods rule (gene purity): %s" % dict(lab_m.value_counts()))
mouse_f = set(cells.well[cells.label_fig1d == "mouse"])
mouse_m = set(cells.well[cells.label_methods == "mouse"])
print("  mouse cells: fig1d n=%d, methods n=%d, intersection n=%d"
      % (len(mouse_f), len(mouse_m), len(mouse_f & mouse_m)))
del ufi, sp

# ===========================================================================
# 3. counting scope, all three datasets
# ===========================================================================
def load_scoped(path, kind, table):
    """Read one count table, normalise columns, apply upstream's scope."""
    print("  reading %-58s (%s, %s)" % (os.path.basename(path), kind, table))
    if kind == "pub":
        df = vp.normalise_columns(read_counts(path, big=True))
        return scope(df, mouse_only=True)
    df = read_counts(path)
    if kind == "own":
        df.columns = [str(c).split("/")[-1].zfill(3) for c in df.columns]
        df = df[[c for c in df.columns if c not in BLANKS]]
    else:
        df.columns = [str(c).split("/")[-1] for c in df.columns]
    return scope(df)


PATHS = {
    ("VASA_published", "-"): ("pub", {"uniagg": PUB_UNIAGG, "total": PUB_TOTAL}),
    ("VASA_own", "-"): ("own", {"uniagg": OWN, "total": OWN_TOTAL}),
    ("FLASH-seq", "native"): ("fs", {t: fs_path("native", t) for t in ("uniagg", "total")}),
    ("FLASH-seq", "vasalen"): ("fs", {t: fs_path("vasalen", t) for t in ("uniagg", "total")}),
}
print("\nloading count tables ...")
D, scope_rows = {}, []
for (ds, arm), (kind, tp) in PATHS.items():
    for table, path in tp.items():
        df, diag = load_scoped(path, kind, table)
        D[(ds, arm, table)] = df
        e = ens_of(df.index)
        keep_shared = e.isin(shared_pc)
        scope_rows.append(dict(
            dataset=ds, arm=arm, table=table, source=os.path.basename(path),
            n_units=df.shape[1], **diag,
            rows_in_shared_universe=int(keep_shared.sum()),
            rows_outside_shared_universe=int((~keep_shared).sum()),
            reads_outside_shared_universe=int(df[(~keep_shared).values].values.sum()),
            duplicate_ens_ids_in_scope=int(len(e) - e.nunique()),
            shared_universe_ids_present_as_rows=int(e[keep_shared].nunique())))
scope_df = pd.DataFrame(scope_rows)
scope_df.to_csv("%s/detection_threeway_scope.tsv" % OUT, sep="\t", index=False)
print("\n=== COUNTING SCOPE (denominators) ===")
print(scope_df[["dataset", "arm", "table", "n_units", "rows_total", "rows_in_scope",
                "rows_in_shared_universe", "reads_in_scope",
                "rows_pc_single_gene_lost_to_hyphen_rule",
                "duplicate_ens_ids_in_scope"]].to_string(index=False))

# ===========================================================================
# 4. per-unit detection, both universes
# ===========================================================================
meta = pd.read_csv(META, sep="\t", comment="#")
amt = dict(zip(meta.library, meta.input_amount.astype(str).str.strip()))
qcv = dict(zip(meta.library, meta.qc_verdict))
lab_f_d, lab_m_d = dict(zip(cells.well, cells.label_fig1d)), dict(zip(cells.well, cells.label_methods))

def annotate(ds, unit):
    """(group, qc_verdict, fig1d label, methods label) for one unit."""
    if ds == "VASA_published":
        return lab_f_d[unit], "n/a", lab_f_d[unit], lab_m_d[unit]
    if ds == "VASA_own":
        return "real cell", "ok", "n/a", "n/a"
    return amt.get(unit, ""), qcv.get(unit, "ok"), "n/a", "n/a"


def unit_label(ds, col):
    return {"VASA_published": "well_%s", "VASA_own": "cell_%s"}.get(ds, "%s") % col


# vec[(table, universe, dataset, arm, unit)] -> the count vector actually
# counted. Cached explicitly so the summary, the curves and the common-depth
# column are provably computed from the same numbers.
vec, rows, sat = {}, [], []
for (ds, arm, table), df in D.items():
    e_in = ens_of(df.index).isin(shared_pc).values
    for col in df.columns:
        unit = unit_label(ds, col)
        grp, qc, lf, lm = annotate(ds, col)
        series = df[col]
        for universe, v in (("all", series.values), ("shared", series[e_in].values)):
            key = (table, universe, ds, arm, unit)
            assert key not in vec, "duplicate unit key %s" % (key,)
            vec[key] = v
            tot = int(v.sum())
            r = dict(table=table, universe=universe, dataset=ds, arm=arm, unit=unit,
                     group=grp, qc_verdict=qc, label_fig1d=lf, label_methods=lm,
                     n_entries_in_scope=int(len(v)), total_reads=tot,
                     genes_ge1=int((v >= 1).sum()), genes_ge5=int((v >= 5).sum()))
            for g in SUMMARY_RUNGS:
                r["genes_at_%s" % format(int(g), "d")] = (
                    round(thinned(v, g / tot), 1) if tot and g <= tot else np.nan)
            rows.append(r)
            for g in GRID:
                if tot and g <= tot:
                    sat.append(dict(table=table, universe=universe, dataset=ds,
                                    arm=arm, unit=unit, group=grp, qc_verdict=qc,
                                    label_fig1d=lf, label_methods=lm,
                                    native_depth=tot, depth=int(g),
                                    genes=round(thinned(v, g / tot), 1)))
det = pd.DataFrame(rows)
sat = pd.DataFrame(sat)
# `unit` is NOT unique: the same FLASH-seq library id (ZHA8833A1..A10) appears in
# BOTH arms, so selecting a track on `unit` alone silently doubles it. Every
# selection below therefore keys on dataset|arm|unit.
for _d in (det, sat):
    _d["key"] = _d.dataset + "|" + _d.arm + "|" + _d.unit
assert det.groupby(["table", "universe"]).key.apply(lambda s: s.is_unique).all()

# Common depth over every unit that enters the headline three-way comparison.
# This is the SHALLOWEST such unit, and on this data that is a published-plate
# well -- which is exactly why the per-rung columns, not this one, carry the
# comparison.
tw = det[(det.table == "uniagg") & (det.universe == "shared") & (det.total_reads > 0)
         & ((det.dataset != "VASA_published") | (det.label_fig1d == "mouse"))]
target = int(tw.total_reads.min())
tf, gcd = [], []
for _, r in det.iterrows():
    if not r.total_reads:
        tf.append(np.nan)
        gcd.append(np.nan)
        continue
    p = min(1.0, target / r.total_reads)
    tf.append(round(p, 5))
    gcd.append(round(thinned(vec[(r.table, r.universe, r.dataset, r.arm, r.unit)], p), 1))
det["thin_fraction_to_common"] = tf
det["genes_at_common_depth"] = gcd
det.to_csv("%s/detection_threeway.tsv" % OUT, sep="\t", index=False)

# --- how big is the uniagg-vs-total table effect, per dataset? --------------
te = det[det.universe == "shared"].pivot_table(
    index=["dataset", "arm", "unit"], columns="table",
    values=["total_reads", "genes_ge1"] + ["genes_at_%d" % g for g in SUMMARY_RUNGS])
eff = []
for (ds, arm), grp in te.groupby(level=[0, 1]):
    row = dict(dataset=ds, arm=arm, n_units=len(grp),
               median_reads_uniagg=float(grp[("total_reads", "uniagg")].median()),
               median_reads_total=float(grp[("total_reads", "total")].median()))
    row["reads_uplift_pct"] = round(
        100 * (row["median_reads_uniagg"] / row["median_reads_total"] - 1), 2)
    for g in SUMMARY_RUNGS:
        d = grp[("genes_at_%d" % g, "uniagg")] - grp[("genes_at_%d" % g, "total")]
        row["genes_uniagg_minus_total_at_%d" % g] = (
            round(float(d.median()), 1) if d.notna().any() else np.nan)
    eff.append(row)
eff = pd.DataFrame(eff)
eff.to_csv("%s/detection_threeway_tableeffect.tsv" % OUT, sep="\t", index=False)
print("\n=== TABLE EFFECT: uniaggGenes minus total, shared universe ===")
print("(this is the confound in mk_gene_detection.py, which read uniagg for the")
print(" own plate but total for FLASH-seq; median over units)")
print(eff[["dataset", "arm", "n_units", "median_reads_uniagg", "median_reads_total",
           "reads_uplift_pct", "genes_uniagg_minus_total_at_1000000"]].to_string(index=False))
sat.to_csv("%s/detection_threeway_saturation.tsv" % OUT, sep="\t", index=False)
print("\ncommon depth over all three-way units (shared universe) = %s reads"
      % format(target, ","))

# --- self-check --------------------------------------------------------------
# Does this reimplementation reproduce res/flashseq_vasa/gene_detection_matched.tsv?
# It must be checked using the REFERENCE's own table choice per dataset (uniagg
# for the own plate, total for FLASH-seq), because that mismatch is the thing this
# script declines to inherit. Reproducing those numbers is what proves the
# estimator and the scope filters are identical; the headline then uses uniagg on
# all sides deliberately.
REF_TABLE = {"VASA_own": "uniagg", "FLASH-seq": "total"}
ref = pd.read_csv(REF_MATCHED, sep="\t")
mine = det[(det.universe == "all") & (det.dataset != "VASA_published")].copy()
mine = mine[[t == REF_TABLE[d] for d, t in zip(mine.dataset, mine.table)]]
mine["ds"] = mine.dataset.replace({"VASA_own": "VASA"})
j = ref.merge(mine, left_on=["dataset", "arm", "unit"],
              right_on=["ds", "arm", "unit"], suffixes=("_ref", "_new"))
assert len(j) == len(ref), "self-check: %d of %d reference rows matched" % (len(j), len(ref))
for col in ("total_reads", "genes_ge1", "genes_ge5"):
    bad = j[j["%s_ref" % col] != j["%s_new" % col]]
    assert bad.empty, "self-check FAILED on %s:\n%s" % (
        col, bad[["unit", "%s_ref" % col, "%s_new" % col]].to_string(index=False))
print("\nself-check: reproduces all %d rows of gene_detection_matched.tsv exactly "
      "(total_reads, genes_ge1, genes_ge5), using that script's own per-dataset "
      "table choice" % len(ref))

# ===========================================================================
# 5. tracks, and where they cross
# ===========================================================================
# The headline five tracks plus the published plate's alternative cell rule. Names
# are the contract between this script, fig_detection_threeway.py and
# verify_detection_threeway.py -- changing one means changing all three.
TRACKS = {
    "VASA published (mouse cells)":
        lambda d: (d.dataset == "VASA_published") & (d.label_fig1d == "mouse"),
    "VASA own plate":
        lambda d: d.dataset == "VASA_own",
    "FLASH-seq native":
        lambda d: (d.arm == "native") & (d.qc_verdict != "exclude"),
    "FLASH-seq VASA-trimmed":
        lambda d: (d.arm == "vasalen") & (d.qc_verdict != "exclude"),
    "FLASH-seq 30 pg (trimmed)":
        lambda d: (d.arm == "vasalen") & (d.group == "30 pg"),
    # Sensitivity track, not a headline one: same plate, the paper's OTHER doublet
    # rule (gene purity rather than UFI purity). It calls 144 mouse cells against
    # the Fig.1d rule's 173, and is carried so the writeup can say whether the
    # choice of rule moves the result.
    "VASA published (mouse cells, Methods rule)":
        lambda d: (d.dataset == "VASA_published") & (d.label_methods == "mouse"),
}
HEADLINE_TRACKS = [t for t in TRACKS if "Methods rule" not in t]
# SUPPORT MODE. Thinning never invents reads, so a rung deeper than a unit's own
# native depth is simply absent for that unit. Two honest ways to report a track
# above its shallowest unit, and this script gives both rather than choosing:
#
#   support='all_units'      every unit of the track contributes. The median
#                            cannot move because membership changed. This is the
#                            only mode in which a between-track difference is
#                            purely a depth-matched difference.
#   support='deepest_subset' only the units that natively reach that depth. The
#                            comparison is still depth-matched, but the cells are
#                            now a DEPTH-SELECTED subset of the track, and deeper
#                            cells are systematically better cells.
#
# `selection_bias_genes` measures the second mode's cost directly: it re-evaluates
# BOTH the full track and the surviving subset at the deepest all-units rung, and
# reports the difference. It is the number of genes the median gains purely from
# dropping the shallow units, with depth held fixed -- so a reader can see whether
# a gap at a deep rung is bigger or smaller than the selection effect that rung
# carries. This matters here: the published plate reaches 500 k on only 17 of 173
# mouse cells.
trk = []
for table in ("uniagg", "total"):
    for universe in ("shared", "all"):
        du = det[(det.universe == universe) & (det.table == table)]
        su = sat[(sat.universe == universe) & (sat.table == table)]
        for tname, fn in TRACKS.items():
            u = set(du[fn(du)].key)
            if not u:
                continue
            sub = su[su.key.isin(u)]
            nat = du[du.key.isin(u)].total_reads
            ref_depth = int(sub[sub.depth <= nat.min()].depth.max())
            base = sub[sub.depth == ref_depth]
            assert len(base) == len(u), "ref rung %d: %d of %d units" % (
                ref_depth, len(base), len(u))
            for g in GRID:
                v = sub[sub.depth == int(g)]
                if not len(v):
                    continue
                here = set(v.key)
                trk.append(dict(
                    table=table, universe=universe, track=tname, depth=int(g),
                    support="all_units" if len(here) == len(u) else "deepest_subset",
                    n_units=len(here), n_units_total=len(u),
                    median_genes=round(v.genes.median(), 1),
                    q1=round(v.genes.quantile(.25), 1),
                    q3=round(v.genes.quantile(.75), 1),
                    min_genes=round(v.genes.min(), 1),
                    max_genes=round(v.genes.max(), 1),
                    native_depth_min=int(nat.min()),
                    native_depth_median=int(nat.median()),
                    native_depth_max=int(nat.max()),
                    ref_depth=ref_depth,
                    selection_bias_genes=round(float(
                        base[base.key.isin(here)].genes.median()
                        - base.genes.median()), 1)))
trk = pd.DataFrame(trk)
trk.to_csv("%s/detection_threeway_tracks.tsv" % OUT, sep="\t", index=False)

cross = []
for table in ("uniagg", "total"):
    for universe in ("shared", "all"):
        t = trk[(trk.universe == universe) & (trk.table == table)]
        names = [n for n in TRACKS if n in set(t.track)]
        for i in range(len(names)):
            for k in range(i + 1, len(names)):
                a, b = names[i], names[k]
                A = t[t.track == a].set_index("depth").median_genes
                B = t[t.track == b].set_index("depth").median_genes
                d = sorted(set(A.index) & set(B.index))
                if len(d) < 2:
                    continue
                diff = np.array([A[x] - B[x] for x in d])
                xs = None
                for j in range(len(d) - 1):
                    if diff[j] == 0 or np.sign(diff[j]) != np.sign(diff[j + 1]):
                        l1, l2 = np.log10(d[j]), np.log10(d[j + 1])
                        f = 0.0 if diff[j] == 0 else diff[j] / (diff[j] - diff[j + 1])
                        xs = 10 ** (l1 + f * (l2 - l1))
                        break
                ta, tb = t[t.track == a].set_index("depth"), t[t.track == b].set_index("depth")
                cross.append(dict(
                    table=table, universe=universe, track_a=a, track_b=b,
                    n_a_total=int(ta.n_units_total.iloc[0]),
                    n_b_total=int(tb.n_units_total.iloc[0]),
                    min_n_a_over_range=int(ta.loc[d, "n_units"].min()),
                    min_n_b_over_range=int(tb.loc[d, "n_units"].min()),
                    support_a="+".join(sorted(set(ta.loc[d, "support"]))),
                    support_b="+".join(sorted(set(tb.loc[d, "support"]))),
                    shared_support_min_depth=int(min(d)),
                    shared_support_max_depth=int(max(d)),
                    n_shared_rungs=len(d),
                    diff_at_min_depth=round(float(diff[0]), 1),
                    diff_at_max_depth=round(float(diff[-1]), 1),
                    diff_min=round(float(diff.min()), 1),
                    diff_max=round(float(diff.max()), 1),
                    abs_diff_min=round(float(np.abs(diff).min()), 1),
                    abs_diff_max=round(float(np.abs(diff).max()), 1),
                    a_leads_throughout=("yes" if (diff > 0).all()
                                        else ("no" if (diff < 0).all() else "mixed")),
                    crosses="yes" if xs else "no",
                    crossing_depth=int(round(xs)) if xs else np.nan,
                    max_selection_bias_a=round(float(ta.loc[d, "selection_bias_genes"].abs().max()), 1),
                    max_selection_bias_b=round(float(tb.loc[d, "selection_bias_genes"].abs().max()), 1)))
cross = pd.DataFrame(cross)
cross.to_csv("%s/detection_threeway_crossings.tsv" % OUT, sep="\t", index=False)

HL = (trk.table == "uniagg") & (trk.universe == "shared")
print("\n=== DEPTH SUPPORT per track (uniagg; depth = reads on in-scope entries) ===")
print(trk[HL].groupby("track").agg(
    n_units=("n_units_total", "first"), native_min=("native_depth_min", "first"),
    native_median=("native_depth_median", "first"),
    native_max=("native_depth_max", "first"),
    all_units_to=("ref_depth", "first"), max_rung=("depth", "max")).to_string())

print("\n=== HEADLINE: median genes, uniagg tables, SHARED universe "
      "(protein-coding in E99 AND E116) ===")
print(trk[HL].pivot_table(index="depth", columns="track",
                          values="median_genes").to_string())
print("\n=== UNRESTRICTED (uniagg, each dataset's own annotation universe) ===")
print(trk[(trk.table == "uniagg") & (trk.universe == "all")].pivot_table(
    index="depth", columns="track", values="median_genes").to_string())
print("\n   annotation-release effect on the headline = unrestricted minus shared:")
sh = trk[HL].pivot_table(index="depth", columns="track", values="median_genes")
al = trk[(trk.table == "uniagg") & (trk.universe == "all")].pivot_table(
    index="depth", columns="track", values="median_genes")
print((al - sh).round(1).to_string())

print("\n=== SUPPORT / DEPTH-SELECTION per rung (published + own plate) ===")
print(trk[HL & trk.track.isin(["VASA published (mouse cells)", "VASA own plate"])][
    ["track", "depth", "support", "n_units", "n_units_total", "median_genes",
     "ref_depth", "selection_bias_genes"]].to_string(index=False))

print("\n=== CROSSINGS (uniagg, shared universe) ===")
cs = cross[(cross.table == "uniagg") & (cross.universe == "shared")]
print(cs[["track_a", "track_b", "shared_support_min_depth", "shared_support_max_depth",
          "diff_at_min_depth", "diff_at_max_depth", "a_leads_throughout", "crosses",
          "crossing_depth"]].to_string(index=False))

# --- three-way convergence / separation -----------------------------------
# The literal three-way question: one track per dataset, at every rung all three
# support. FLASH-seq is represented by its VASA-trimmed arm (read-length-matched
# to VASA, so the arm that isolates protocol rather than read length).
THREE = ["VASA published (mouse cells)", "VASA own plate", "FLASH-seq VASA-trimmed"]
assert all(t in sh.columns for t in THREE), sorted(sh.columns)
sp = sh[THREE].dropna()
conv = pd.DataFrame(dict(
    depth=sp.index, spread_genes=(sp.max(axis=1) - sp.min(axis=1)).round(1).values,
    spread_pct_of_lowest=(100 * (sp.max(axis=1) - sp.min(axis=1))
                          / sp.min(axis=1)).round(2).values,
    highest=sp.idxmax(axis=1).values, lowest=sp.idxmin(axis=1).values,
    **{t.replace(" ", "_"): sp[t].values for t in THREE}))
conv.to_csv("%s/detection_threeway_convergence.tsv" % OUT, sep="\t", index=False)
print("\n=== THREE-WAY SPREAD (published / own / FLASH-seq VASA-trimmed) ===")
print("Do they converge or separate? spread = max minus min of the three medians")
print(conv[["depth", "spread_genes", "spread_pct_of_lowest", "highest", "lowest"]]
      .to_string(index=False))

print("\nwrote to %s:" % OUT)
for f in ("detection_threeway.tsv", "detection_threeway_saturation.tsv",
          "detection_threeway_tracks.tsv", "detection_threeway_crossings.tsv",
          "detection_threeway_scope.tsv", "detection_threeway_tableeffect.tsv",
          "detection_threeway_convergence.tsv", "annotation_universe.tsv",
          "published_cell_species.tsv"):
    print("  %-42s %10d bytes" % (f, os.path.getsize("%s/%s" % (OUT, f))))
