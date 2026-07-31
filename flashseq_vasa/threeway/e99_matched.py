#!/usr/bin/env python3
"""Own VASA plate vs published VASA plate under ONE annotation (Ensembl 99).

WHY
---
`res/threeway/composition_threeway.tsv` compares the published plate (GRCm38 /
Ensembl 99, human+mouse) with the own plate (GRCm39 / Ensembl 116, mouse only).
Every published-vs-own difference in it therefore mixes protocol-free biology
with annotation release. `threeway_release_control.tsv` corrects for that on the
composition scale by restricting to genes present in both releases -- but each
dataset is still quantified under its OWN annotation, gene models and genome
build, so it is a PARTIAL control.

This script uses the complete control: the own plate re-quantified through the
published plate's own reference (`out_E99/`, built by run_own_e99.sh), so the two
VASA plates are compared with nothing left to correct for.

WHAT IS HELD CONSTANT, AND WHAT IS NOT
--------------------------------------
Constant across the two arms compared here: genome build (GRCm38), annotation
(Ensembl 99), the annotation BED, the step-5 assignment code (both use upstream
a_Mapping/), STAR version (2.7.7a), STAR parameters, `stranded=y`, and the
counting code.

NOT constant, and reported rather than hidden:
  1. sjdbOverhang 150 (own) vs 73 (published). Forced by read length -- 130 nt
     vs 75 nt reads -- and using 73 for 130 nt reads would silently drop
     junction-spanning alignments. Sized correctly per library, not equalised.
  2. Read length itself, 130 nt vs 75 nt. Unfixable without truncating reads;
     `readlen_effect_pp` in the existing three-way table is the FLASH-seq-derived
     estimate of its size, and Trap 8 (feature shorter than one read) applies.
  3. Sequencing depth per cell, ~4.4M vs ~0.5M reads. Handled by depth-matched
     detection, not by composition.
  4. n: 12 real own-plate cells vs the published plate's mouse-called cells.

Rules honoured: ReadCounts on both sides (Rule 4, and Rule 4 for cross-protocol
comparison); non-rRNA denominator stated everywhere (Rule 5); nothing silently
filtered (Rule 3).

Outputs, all under res/threeway/:
  own_plate_E99.tsv              own plate under E99 vs under E116, per biotype
  e99_matched_composition.tsv    both plates under E99, the payoff table
  e99_matched_structural.tsv     the structural-RNA gap under each control
  e99_matched_assignment.tsv     reads assigned under E99 but not E116, and back
  e99_matched_detection.tsv      depth-matched gene detection under E99
  e99_matched_report.txt         the log
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/vasaplate_check")
import vp_common as vp  # noqa: E402  species_of / classify_fig1d / classify_methods / MIN_UFI

W = "/nemo/lab/turnerj/working/guangxin/vasaseq"
REF = "/nemo/lab/turnerj/working/guangxin/reference/vasaseq"
RES = f"{W}/res/threeway"

PUB_UNIAGG = f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_uniaggGenes_total.ReadCounts.tsv"
PUB_UFI = f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_total.UFICounts.tsv"
OWN99_UNIAGG = f"{W}/data/PM26037/out_E99/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv"
OWN116_UNIAGG = f"{W}/data/PM26037/out/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv"
BED99 = f"{REF}/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed"
BED116 = f"{REF}/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed"

BLANKS = {"001", "014", "015", "016"}          # own-plate empty wells, verbatim
RRNA_PRIMARY = {"rRNA", "Mt_rRNA"}             # verbatim from the two-way scripts
STRUCTURAL = ["MiscRna", "snRNA", "snoRNA", "scaRNA", "ribozyme"]
CHUNK = 40_000

# Values this script must reproduce from the committed tables rather than trust.
# Read from res/threeway/{composition_threeway,threeway_release_control,threeway_structural}.tsv.
EXPECT = {}
log = []


def say(s=""):
    print(s, flush=True)
    log.append(s)


def sha(path, n=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(n)
            if not b:
                break
            h.update(b)
    return h.hexdigest()[:16]


def biotype_token(labels):
    """Verbatim rule from mk_vasa_composition.py / threeway_composition.py."""
    return pd.Series([str(i).rsplit("_", 1)[-1] if "_" in str(i) else "NA" for i in labels],
                     index=labels)


def norm_cols(cols):
    """Column labels -> bare unit ids, covering all three naming conventions."""
    out = []
    for c in cols:
        c = str(c).replace("cells/", "")
        c = c.rsplit("/", 1)[-1]
        if c.startswith("SRR"):
            c = c.rsplit("_", 1)[-1].zfill(3)
        elif c.isdigit():
            c = c.zfill(3)
        out.append(c)
    return out


def bed_gene_biotypes(path, mouse_only):
    """{gene_id: biotype} from a homemade_IntronExonTrna BED (col 5 = name).

    Verbatim from threeway_composition.py so the biotype vocabulary is identical.
    """
    got = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            name = f[4]
            if not name.startswith("ENS"):
                continue
            if mouse_only and not name.startswith("ENSMUSG"):
                continue
            body = name.rsplit("_", 1)[0]
            gid = body.split("_", 1)[0]
            bt = body.rsplit("_", 1)[-1]
            got[gid] = bt
    return got


class Comp:
    """Read sums per biotype on the non-rRNA denominator.

    Same accumulator contract as threeway_composition.py's Comp, plus the
    species split this arm needs: the own plate is a mouse-only library on a
    human+mouse reference, so anything landing on a human gene id is
    mismapping and has to be measured, not assumed away.
    """

    def __init__(self, name, species_filter=None):
        self.name, self.species_filter = name, species_filter
        self.bt = {}
        self.gene_reads = {}         # gene_id -> reads (simple rows only)
        self.n_entries = self.n_detected = 0
        self.reads_all = 0.0
        self.reads_rrna = 0.0
        self.reads_mtrrna = 0.0
        self.reads_dropped_species = 0.0
        self.reads_human = 0.0       # rows whose id is human, before any filter
        self.reads_combo = 0.0
        self.reads_combo_multibt = 0.0
        self.percell = {}            # unit -> reads on the denominator

    def add(self, df):
        rs = df.sum(axis=1)
        self.reads_all += float(rs.sum())
        idx = df.index.astype(str)

        spv = np.array([vp.species_of(i) for i in idx])
        self.reads_human += float(rs.values[spv == "human"].sum())
        if self.species_filter is not None:
            keep = np.isin(spv, self.species_filter)
            self.reads_dropped_species += float(rs.values[~keep].sum())
            df, rs, idx = df[keep], rs[keep], idx[keep]
            if len(df) == 0:
                return

        bt = biotype_token(list(idx))
        rr = bt.isin(RRNA_PRIMARY).values
        self.reads_rrna += float(rs.values[rr].sum())
        df, rs, idx, bt = df[~rr], rs[~rr], idx[~rr], bt[~rr]
        self.n_entries += len(df)
        self.n_detected += int((rs > 0).sum())
        self.reads_mtrrna += float(rs.values[bt.isin(["MtRrna"]).values].sum())

        for b, v in rs.groupby(bt.values).sum().items():
            self.bt[b] = self.bt.get(b, 0.0) + float(v)
        for col, v in df.sum(axis=0).items():
            self.percell[col] = self.percell.get(col, 0.0) + float(v)

        is_combo = np.array(["-" in i for i in idx])
        self.reads_combo += float(rs.values[is_combo].sum())
        multi = np.array([
            len({p.rsplit("_", 1)[-1] for p in i.split("-")}) > 1 if "-" in i else False
            for i in idx])
        self.reads_combo_multibt += float(rs.values[multi].sum())

        # simple rows -> per-gene reads, for the assignment-churn comparison
        srs, sidx = rs[~is_combo], idx[~is_combo]
        if len(srs):
            gid = [i.split("_", 1)[0] for i in sidx]
            for g, v in srs.groupby(np.array(gid)).sum().items():
                self.gene_reads[g] = self.gene_reads.get(g, 0.0) + float(v)

    def denom(self):
        return sum(self.bt.values())

    def pct(self):
        t = self.denom()
        return pd.Series({k: 100.0 * v / t for k, v in self.bt.items()}) if t else pd.Series(dtype=float)


def stream(comp, path, cols=None, chunk=CHUNK):
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=chunk):
        ch.columns = norm_cols(ch.columns)
        if cols is not None:
            ch = ch[[c for c in cols if c in ch.columns]]
        comp.add(ch)
    return comp


# ==========================================================================
say("=" * 78)
say("INPUTS")
say("=" * 78)
for p in (PUB_UNIAGG, OWN99_UNIAGG, OWN116_UNIAGG, BED99, BED116):
    if not os.path.exists(p):
        say(f"  MISSING {p}")
        sys.exit(f"required input absent: {p}")
    say(f"  {sha(p)}  {os.path.getsize(p):>12,d}  {p}")

bt99 = bed_gene_biotypes(BED99, mouse_only=True)
bt116 = bed_gene_biotypes(BED116, mouse_only=False)
shared = set(bt99) & set(bt116)
say(f"  Ensembl  99 mouse gene ids: {len(bt99):,}")
say(f"  Ensembl 116 mouse gene ids: {len(bt116):,}")
say(f"  shared                    : {len(shared):,}")

# ---- published plate: which barcodes are mouse (Fig.1d rule, as before) ----
say()
say("=" * 78)
say("PUBLISHED PLATE -- MOUSE CELL CALLING (Fig.1d rule, as in threeway_composition.py)")
say("=" * 78)
h_ufi = m_ufi = gh = gm = None
for ch in pd.read_csv(PUB_UFI, sep="\t", index_col=0, chunksize=CHUNK):
    ch.columns = norm_cols(ch.columns)
    sp = pd.Series([vp.species_of(str(i)) for i in ch.index], index=ch.index).values
    hh, mm = ch[sp == "human"], ch[sp == "mouse"]

    def acc(a, b):
        return b if a is None else a.add(b, fill_value=0)

    h_ufi, m_ufi = acc(h_ufi, hh.sum()), acc(m_ufi, mm.sum())
    gh, gm = acc(gh, (hh > 0).sum()), acc(gm, (mm > 0).sum())
lab_f, _ = vp.classify_fig1d(h_ufi, m_ufi)
MOUSE_F = sorted(lab_f[lab_f == "mouse"].index)
say(f"  mouse-called barcodes: {len(MOUSE_F)}")

# ---- the three arms -------------------------------------------------------
say()
say("=" * 78)
say("COMPOSITION -- ReadCounts, non-rRNA denominator, mouse entries only")
say("=" * 78)
head = pd.read_csv(OWN99_UNIAGG, sep="\t", index_col=0, nrows=1)
own_cols = [c for c in norm_cols(head.columns) if c not in BLANKS]
assert len(own_cols) == 12, f"expected 12 real own-plate cells, got {len(own_cols)}: {own_cols}"

ds = {}
# Published plate: species filter to mouse, as the existing table does.
ds["published_E99"] = stream(Comp("published_E99", species_filter=("mouse",)),
                             PUB_UNIAGG, cols=MOUSE_F)
# Own plate under E99: the SAME mixed reference, so the SAME species filter
# applies. Its human share is a pure artefact and is reported separately.
ds["own_E99"] = stream(Comp("own_E99", species_filter=("mouse",)),
                       OWN99_UNIAGG, cols=own_cols)
# Own plate under E116: mouse-only reference, no filter needed.
ds["own_E116"] = stream(Comp("own_E116"), OWN116_UNIAGG, cols=own_cols)

for k, cp in ds.items():
    d = cp.denom()
    say(f"  {k:14s} entries={cp.n_entries:9,d} detected={cp.n_detected:9,d}")
    say(f"  {'':14s} reads_all={cp.reads_all:14,.0f}  rRNA_token={cp.reads_rrna:11,.0f}"
        f"  denominator={d:14,.0f}")
    say(f"  {'':14s} human-id reads={cp.reads_human:12,.0f}"
        f" ({100 * cp.reads_human / cp.reads_all:.4f}% of reads_all)"
        f"  combos={100 * cp.reads_combo / d:.2f}%")

# The own plate's human share under the mixed reference is the new artefact.
say()
say("  NEW ARTEFACT -- mouse-only library on a human+mouse reference:")
own_h = 100 * ds["own_E99"].reads_human / ds["own_E99"].reads_all
pub_h_of_mouse_cells = 100 * ds["published_E99"].reads_human / ds["published_E99"].reads_all
say(f"    own plate  reads on human gene ids: {own_h:.4f}% of all rows")
say(f"    published mouse-called cells      : {pub_h_of_mouse_cells:.4f}% of all rows")
say("    The published figure is the empirical yardstick: those cells are mouse")
say("    cells measured on the same mixed reference, so their human share is")
say("    also mismapping plus ambient. A comparable own-plate figure means the")
say("    mixed reference costs the own plate no more than it cost the published")
say("    plate, and the comparison is not biased by it.")

pct = {k: cp.pct() for k, cp in ds.items()}
comp = pd.DataFrame(pct).fillna(0.0)
comp.index.name = "biotype"
comp = comp.sort_values("own_E99", ascending=False)
comp["gap_pub_minus_own_E99_pp"] = comp.published_E99 - comp.own_E99
comp["own_E99_minus_E116_pp"] = comp.own_E99 - comp.own_E116
reads = pd.DataFrame({f"reads_{k}": pd.Series(cp.bt) for k, cp in ds.items()}).fillna(0.0)
comp = comp.join(reads.astype("int64"))
comp.round(6).to_csv(f"{RES}/e99_matched_composition.tsv", sep="\t")

# own_plate_E99.tsv -- the own plate under both annotations, the deliverable the
# task names explicitly.
own = pd.DataFrame({
    "pct_under_E99": pct["own_E99"],
    "pct_under_E116": pct["own_E116"],
    "reads_under_E99": pd.Series(ds["own_E99"].bt),
    "reads_under_E116": pd.Series(ds["own_E116"].bt),
}).fillna(0.0)
own.index.name = "biotype"
own["delta_pp_E99_minus_E116"] = own.pct_under_E99 - own.pct_under_E116
own = own.sort_values("pct_under_E99", ascending=False)
own.round(6).to_csv(f"{RES}/own_plate_E99.tsv", sep="\t")

# ---- THE KEY NUMBER -------------------------------------------------------
say()
say("=" * 78)
say("THE KEY NUMBER -- structural RNA gap between the two VASA plates")
say("=" * 78)


def struct_pct(series):
    return float(sum(series.get(b, 0.0) for b in STRUCTURAL))


s_pub = struct_pct(pct["published_E99"])
s_own99 = struct_pct(pct["own_E99"])
s_own116 = struct_pct(pct["own_E116"])

# Read the two prior controls from the committed tables -- never transcribed.
prior = pd.read_csv(f"{RES}/threeway_release_control.tsv", sep="\t", index_col=0)
gap_raw = float(sum(prior.loc[b, "gap_raw_pp"] for b in STRUCTURAL if b in prior.index))
gap_shared = float(sum(prior.loc[b, "gap_shared_pp"] for b in STRUCTURAL if b in prior.index))
struct_tbl = pd.read_csv(f"{RES}/threeway_structural.tsv", sep="\t")
s_pub_prior = float(struct_tbl.set_index("dataset").loc["published_vasa", "structural_pct_of_nonrRNA"])
s_own_prior = float(struct_tbl.set_index("dataset").loc["own_vasa", "structural_pct_of_nonrRNA"])

gap_matched = s_own99 - s_pub
say(f"  published plate, E99            : {s_pub:8.4f}%   (committed: {s_pub_prior:.4f}%)")
say(f"  own plate,       E116 (as before): {s_own116:8.4f}%   (committed: {s_own_prior:.4f}%)")
say(f"  own plate,       E99 (matched)   : {s_own99:8.4f}%")
say()
say(f"  gap raw            (own E116 - pub E99): {abs(-gap_raw):8.4f} pp   [committed table]")
say(f"  gap shared-universe control            : {abs(-gap_shared):8.4f} pp   [committed table]")
say(f"  gap FULLY MATCHED  (own E99  - pub E99): {gap_matched:8.4f} pp   [this run]")

# Per-class, which is what "name it per class" asks for.
rows = []
for b in STRUCTURAL + ["ProteinCoding", "lncRNA"]:
    p, o99, o116 = (pct["published_E99"].get(b, 0.0), pct["own_E99"].get(b, 0.0),
                    pct["own_E116"].get(b, 0.0))
    rows.append({
        "biotype": b,
        "published_E99_pct": p,
        "own_E116_pct": o116,
        "own_E99_pct": o99,
        "gap_raw_pp": o116 - p,
        "gap_shared_pp": (-float(prior.loc[b, "gap_shared_pp"]) if b in prior.index else np.nan),
        "gap_matched_pp": o99 - p,
        "release_explained_pp": (o116 - p) - (o99 - p),
        "pct_of_raw_gap_explained_by_release":
            (100.0 * ((o116 - p) - (o99 - p)) / (o116 - p)) if (o116 - p) != 0 else np.nan,
    })
st = pd.DataFrame(rows)
st["structural"] = st.biotype.isin(STRUCTURAL)
st.round(6).to_csv(f"{RES}/e99_matched_structural.tsv", sep="\t", index=False)
say()
say("  per class (pp, positive = own above published):")
for _, r in st.iterrows():
    say(f"    {r.biotype:16s} raw {r.gap_raw_pp:+9.4f}  matched {r.gap_matched_pp:+9.4f}"
        f"  release explains {r.pct_of_raw_gap_explained_by_release:7.1f}%")

# ---- assignment churn ----------------------------------------------------
say()
say("=" * 78)
say("ASSIGNMENT CHURN -- own plate, E99 vs E116")
say("=" * 78)
g99, g116 = ds["own_E99"].gene_reads, ds["own_E116"].gene_reads
ids99, ids116 = set(g99), set(g116)
r99 = sum(g99.values())
r116 = sum(g116.values())
only99 = sum(g99[g] for g in ids99 - ids116)
only116 = sum(g116[g] for g in ids116 - ids99)
both = ids99 & ids116
say(f"  simple-row reads assigned under E99 : {r99:14,.0f}  ({len(ids99):,} gene ids)")
say(f"  simple-row reads assigned under E116: {r116:14,.0f}  ({len(ids116):,} gene ids)")
say(f"  on gene ids present ONLY in E99     : {only99:14,.0f}  ({100 * only99 / r99:.3f}% of E99)")
say(f"  on gene ids present ONLY in E116    : {only116:14,.0f}  ({100 * only116 / r116:.3f}% of E116)")
say(f"  shared gene ids                     : {len(both):,}")
# Total assigned differs for a second reason too: STAR mapped to a different
# genome build, so the input to step 5 is not the same read set.
say(f"  net change in assigned reads        : {r99 - r116:+14,.0f}"
    f"  ({100 * (r99 - r116) / r116:+.3f}%)")
churn = pd.DataFrame([{
    "reads_simple_E99": r99, "reads_simple_E116": r116,
    "n_gene_ids_E99": len(ids99), "n_gene_ids_E116": len(ids116),
    "n_gene_ids_shared": len(both),
    "reads_on_ids_only_in_E99": only99, "reads_on_ids_only_in_E116": only116,
    "pct_reads_on_ids_only_in_E99": 100 * only99 / r99,
    "pct_reads_on_ids_only_in_E116": 100 * only116 / r116,
    "net_read_change_E99_minus_E116": r99 - r116,
    "own_plate_human_id_reads_E99": ds["own_E99"].reads_human,
    "own_plate_pct_human_id_reads_E99": own_h,
    "published_pct_human_id_reads_mousecells": pub_h_of_mouse_cells,
}])
churn.round(6).to_csv(f"{RES}/e99_matched_assignment.tsv", sep="\t", index=False)

# ---- gene detection, depth-matched, both plates under E99 ----------------
say()
say("=" * 78)
say("GENE DETECTION -- depth-matched, both plates under Ensembl 99")
say("=" * 78)


def thinned(counts, p):
    """E[# genes with >=1 read] after binomial thinning at rate p.

    VERBATIM from mk_detection_threeway.py (itself verbatim from
    res/flashseq_vasa/mk_gene_detection.py). Exact in expectation,
    deterministic, no seed.
    """
    cc = counts[counts > 0].astype(float)
    if p >= 1.0:
        return float(len(cc))
    return float(np.sum(1.0 - np.power(1.0 - p, cc)))


def scope(df, mouse_only=False):
    """Upstream's counting scope: single-gene protein-coding entries.

    The keep MASK is byte-identical to mk_detection_threeway.py's scope()
    (`single & pc [& mouse]`), so the rows counted here are exactly the rows it
    counts. This is NOT the whole upstream function: its `ntok`/`lost`
    hyphen-rule diagnostics and its `(df, diag)` return contract are dropped,
    because the hyphen-rule loss is already reported for these same tables in
    detection_threeway_scope.tsv and is not re-derived here. Returns the kept
    frame only.
    """
    idx = df.index.astype(str)
    bio = pd.Series([i.rsplit("_", 1)[-1] if "_" in i else "NA" for i in idx], index=df.index)
    single = ~pd.Series([("-" in i) for i in idx], index=df.index)
    pc = bio == "ProteinCoding"
    mouse = pd.Series([i.startswith("ENSMUSG") for i in idx], index=df.index)
    keep = single & pc
    if mouse_only:
        keep = keep & mouse
    return df[keep.values]


# 'depth' = reads in the counted scope, as in mk_detection_threeway.py.
GRID = [1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6]

# Both plates are now under E99, so the release-driven universe difference is
# gone by construction -- no shared-universe restriction is needed here, which
# is the point. The scope is E99 protein-coding, mouse, single-gene.
det_rows = []
for arm, path, cols in (("published_E99", PUB_UNIAGG, MOUSE_F),
                        ("own_E99", OWN99_UNIAGG, own_cols),
                        ("own_E116", OWN116_UNIAGG, own_cols)):
    # Each chunk's column is a vector over the GENES IN THAT CHUNK, so chunks
    # must be CONCATENATED, never added: `prev + v` would add chunk-1 gene i to
    # chunk-2 gene i, which either raises on unequal post-scope row counts or
    # silently merges distinct genes and understates both genes_native and
    # thinned(). (This differs from the composition path above, where percell
    # accumulates a per-cell SCALAR and += is right.)
    chunks = {}
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=CHUNK):
        ch.columns = norm_cols(ch.columns)
        ch = ch[[c for c in cols if c in ch.columns]]
        sc = scope(ch, mouse_only=True)
        if not len(sc):
            continue
        for col in sc.columns:
            chunks.setdefault(col, []).append(sc[col].values.astype(float))
    percell = {col: np.concatenate(v) for col, v in chunks.items()}
    # Every cell must end up with one entry per in-scope gene, and all cells of
    # an arm read the same rows, so the vectors must all be the same length.
    lens = {len(v) for v in percell.values()}
    assert len(lens) <= 1, f"{arm}: ragged per-cell gene vectors {sorted(lens)}"
    say(f"  {arm:14s} in-scope genes {lens.pop() if lens else 0:,} x {len(percell)} units")
    for col, v in percell.items():
        tot = float(v.sum())
        det_rows.append(dict(arm=arm, unit=col, native_depth=int(tot),
                             genes_native=int((v > 0).sum())))
        for g in GRID:
            if g <= tot:
                det_rows.append(dict(arm=arm, unit=col, depth=int(g),
                                     native_depth=int(tot),
                                     genes_expected=thinned(v, g / tot)))
det = pd.DataFrame(det_rows)
det.round(4).to_csv(f"{RES}/e99_matched_detection.tsv", sep="\t", index=False)

# Headline: the deepest rung every unit of BOTH E99 arms reaches natively.
nat = det.dropna(subset=["native_depth"]).groupby(["arm", "unit"]).native_depth.max()
e99_units = nat.loc[["published_E99", "own_E99"]]
common = max([g for g in GRID if g <= e99_units.min()], default=None)
say(f"  native depth (scope reads): published min {int(nat['published_E99'].min()):,} "
    f"median {int(nat['published_E99'].median()):,} | own min {int(nat['own_E99'].min()):,} "
    f"median {int(nat['own_E99'].median()):,}")
if common is None:
    say("  no common rung -- the two arms do not overlap in depth; detection not compared")
else:
    say(f"  deepest rung ALL units of both E99 arms reach natively: {int(common):,} reads")
    sub = det[(det.depth == common) & det.arm.isin(["published_E99", "own_E99", "own_E116"])]
    for arm, grp in sub.groupby("arm"):
        say(f"    {arm:14s} n={len(grp):3d}  genes at {int(common):,} reads: "
            f"median {grp.genes_expected.median():8.1f}  mean {grp.genes_expected.mean():8.1f}")
    pubm = sub[sub.arm == "published_E99"].genes_expected.median()
    ownm = sub[sub.arm == "own_E99"].genes_expected.median()
    o116 = sub[sub.arm == "own_E116"].genes_expected
    say(f"    own/published at matched depth, both under E99: {ownm / pubm:.3f}x")
    if len(o116):
        say(f"    own under E116 for comparison: median {o116.median():.1f} "
            f"({o116.median() / ownm:.3f}x its own E99 value)")
    # Underpowered by construction: 12 own cells. State it rather than imply
    # significance.
    say("    n=12 own cells vs "
        f"{len(sub[sub.arm == 'published_E99'])} published; a ratio at this n is "
        "descriptive, not a test.")

# ---- provenance + report ------------------------------------------------
prov = pd.DataFrame({
    "file": [PUB_UNIAGG, OWN99_UNIAGG, OWN116_UNIAGG, BED99, BED116],
    "sha256_16": [sha(p) for p in (PUB_UNIAGG, OWN99_UNIAGG, OWN116_UNIAGG, BED99, BED116)],
    "bytes": [os.path.getsize(p) for p in (PUB_UNIAGG, OWN99_UNIAGG, OWN116_UNIAGG, BED99, BED116)],
})
prov.to_csv(f"{RES}/e99_matched_provenance.tsv", sep="\t", index=False)
with open(f"{RES}/e99_matched_report.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")
say()
say("wrote: own_plate_E99.tsv, e99_matched_composition.tsv, e99_matched_structural.tsv,")
say("       e99_matched_assignment.tsv, e99_matched_detection.tsv,")
say("       e99_matched_provenance.tsv, e99_matched_report.txt")
