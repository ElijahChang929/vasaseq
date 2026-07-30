#!/usr/bin/env python3
"""Why do the two VASA plates disagree by ~9x on structural RNA classes?

`threeway_composition.py` found the published plate at 2.25% structural
(MiscRna+snRNA+snoRNA+scaRNA+ribozyme, non-rRNA read denominator) against the own
plate's 20.59% and FLASH-seq's 0.26%. So the published plate sits BETWEEN, not
with the own plate, and the two-way conclusion ("protocol difference") does not
survive as stated. This script tests the four confounds that could each
manufacture that gap, before any of it is called biology.

  A. ALLOCATION. The biotype rule both upstream scripts use is the token after
     the last '_' of the row label. On a multi-gene combination key
     ('A_ProteinCoding-B_snRNA') that silently credits the LAST member. The
     release control showed own-plate MiscRna is 9.82% raw but 0.030% on shared
     simple rows -- so almost all of that read mass is in combination keys or
     E116-only genes. Decomposed here into four disjoint buckets, and the
     structural share recomputed under three allocation rules.
  B. DUPLICATION. The own plate carries ~18x more reads per unit than the
     published plate. Reads include PCR duplicates, which pile up on the most
     abundant molecules, so a deeper library inflates abundant short RNAs in READ
     space specifically. Tested as reads-vs-molecules on BOTH VASA plates
     (TranscriptCounts exists for both) -- if it is duplication, the gap shrinks
     in molecule space.
  C. GEOMETRY (conventions trap 8). Stage 6 keeps a non-splicing biotype only if
     the read is CONTAINED in the feature (jS:IN). Read length therefore gates
     short features. Measured per plate from the per-cell annotated BEDs with a
     deterministic stride (trap 4), reporting the realised jS:IN rate rather
     than assuming it.
  D. POOLING. A share computed over pooled units can be carried by one outlier.
     Reported per unit, with the spread.

Usage: threeway_mechanism.py <RESDIR>
"""
import glob
import gzip
import os
import random
import sys
import numpy as np
import pandas as pd

W = "/nemo/lab/turnerj/working/guangxin/vasaseq"
sys.path.insert(0, f"{W}/code/I_Gene_expression/vasaplate_check")
import vp_common as vp  # noqa: E402

RES = sys.argv[1]
FSDIR = "/nemo/lab/turnerj/scratch/zhangg/vasaseq/flashseq_vasa"
PUB = f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_uniaggGenes_total"
OWN = f"{W}/data/PM26037/out/ZHA9292A1_uniaggGenes_total"
STRUCTURAL = ["MiscRna", "snRNA", "snoRNA", "scaRNA", "ribozyme"]
RRNA_PRIMARY = {"rRNA", "Mt_rRNA"}
BLANKS = {"001", "014", "015", "016"}
CHUNK = 40_000
SEED = 20260730
log = []


def say(s=""):
    print(s, flush=True)
    log.append(s)


def norm_cols(cols):
    out = []
    for c in cols:
        c = str(c).replace("cells/", "").rsplit("/", 1)[-1]
        if c.startswith("SRR"):
            c = c.rsplit("_", 1)[-1].zfill(3)
        elif c.isdigit():
            c = c.zfill(3)
        out.append(c)
    return out


def token(lbl):
    return lbl.rsplit("_", 1)[-1] if "_" in lbl else "NA"


# published mouse cells, Fig.1d rule -- recomputed, not carried over
say("=" * 78)
say("PUBLISHED-PLATE MOUSE CELLS (Fig.1d rule, recomputed here)")
say("=" * 78)
h = m = None
for ch in pd.read_csv(f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_total.UFICounts.tsv",
                      sep="\t", index_col=0, chunksize=CHUNK):
    ch.columns = norm_cols(ch.columns)
    sp = np.array([vp.species_of(str(i)) for i in ch.index])
    a, b = ch[sp == "human"].sum(), ch[sp == "mouse"].sum()
    h = a if h is None else h.add(a, fill_value=0)
    m = b if m is None else m.add(b, fill_value=0)
lab, _ = vp.classify_fig1d(h, m)
MOUSE = sorted(lab[lab == "mouse"].index)
say(f"  mouse cells: {len(MOUSE)}")
assert len(MOUSE) == 173, f"expected 173 mouse cells (composition run), got {len(MOUSE)}"

own_head = pd.read_csv(f"{OWN}.ReadCounts.tsv", sep="\t", index_col=0, nrows=1)
OWNC = [c for c in norm_cols(own_head.columns) if c not in BLANKS]
assert len(OWNC) == 12


# --------------------------------------------------------------------------
# A. allocation: decompose the structural read mass, and re-allocate 3 ways
# --------------------------------------------------------------------------
def decompose(path, cols, species_filter, name):
    """Read mass per biotype under 3 allocation rules + 4 provenance buckets."""
    last = {}          # last-token rule (what the pipeline scripts do)
    unan = {}          # combos counted only when all members agree
    frac = {}          # combo reads split equally across member biotypes
    buckets = {b: dict(simple=0.0, combo_unanimous=0.0, combo_discordant=0.0)
               for b in STRUCTURAL}
    denom = 0.0
    n_units = None
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=CHUNK):
        ch.columns = norm_cols(ch.columns)
        if cols is not None:
            ch = ch[cols]
        n_units = ch.shape[1]
        idx = ch.index.astype(str)
        if species_filter is not None:
            keep = np.isin([vp.species_of(i) for i in idx], species_filter)
            ch, idx = ch[keep], idx[keep]
            if len(ch) == 0:
                continue
        bt = np.array([token(i) for i in idx])
        keep = ~np.isin(bt, list(RRNA_PRIMARY))
        ch, idx, bt = ch[keep], idx[keep], bt[keep]
        rs = ch.sum(axis=1).values.astype(float)
        denom += rs.sum()
        for i, lbl in enumerate(idx):
            r = rs[i]
            if r == 0:
                continue
            b = bt[i]
            last[b] = last.get(b, 0.0) + r
            if "-" in lbl:
                mem = [token(p) for p in lbl.split("-")]
                uniq = set(mem)
                if len(uniq) == 1:
                    unan[b] = unan.get(b, 0.0) + r
                    if b in buckets:
                        buckets[b]["combo_unanimous"] += r
                else:
                    if b in buckets:
                        buckets[b]["combo_discordant"] += r
                for bb in uniq:
                    frac[bb] = frac.get(bb, 0.0) + r * mem.count(bb) / len(mem)
            else:
                unan[b] = unan.get(b, 0.0) + r
                frac[b] = frac.get(b, 0.0) + r
                if b in buckets:
                    buckets[b]["simple"] += r
    return dict(name=name, n_units=n_units, denom=denom,
                last=last, unan=unan, frac=frac, buckets=buckets)


say()
say("=" * 78)
say("A. ALLOCATION -- where the structural read mass actually sits")
say("=" * 78)
srcs = [("published_vasa", f"{PUB}.ReadCounts.tsv", MOUSE, ("mouse",)),
        ("own_vasa", f"{OWN}.ReadCounts.tsv", OWNC, None),
        ("flashseq_native", f"{FSDIR}/native/FSall10_native_uniaggGenes_total.ReadCounts.tsv", None, None),
        ("flashseq_vasalen", f"{FSDIR}/vasalen/FSall10_vasalen_uniaggGenes_total.ReadCounts.tsv", None, None)]
dec = {k: decompose(p, c, s, k) for k, p, c, s in srcs}

alloc_rows, bucket_rows = [], []
for k, d in dec.items():
    row = dict(dataset=k, n_units=d["n_units"], denominator_reads=d["denom"])
    for rule in ("last", "unan", "frac"):
        row[f"structural_pct_{rule}"] = 100 * sum(d[rule].get(b, 0.0) for b in STRUCTURAL) / d["denom"]
    alloc_rows.append(row)
    for b in STRUCTURAL:
        bk = d["buckets"][b]
        t = sum(bk.values())
        bucket_rows.append(dict(dataset=k, biotype=b, reads_total=t,
                                pct_of_denominator=100 * t / d["denom"],
                                pct_from_simple_rows=100 * bk["simple"] / t if t else np.nan,
                                pct_from_unanimous_combos=100 * bk["combo_unanimous"] / t if t else np.nan,
                                pct_from_discordant_combos=100 * bk["combo_discordant"] / t if t else np.nan))
al = pd.DataFrame(alloc_rows).set_index("dataset")
fsn = al.loc["flashseq_native"]
for rule in ("last", "unan", "frac"):
    al[f"fold_vs_FSnative_{rule}"] = al[f"structural_pct_{rule}"] / fsn[f"structural_pct_{rule}"]
al.round(6).to_csv(f"{RES}/threeway_allocation_rules.tsv", sep="\t")
bk = pd.DataFrame(bucket_rows)
bk.round(6).to_csv(f"{RES}/threeway_structural_provenance.tsv", sep="\t")
say("  structural % of the non-rRNA read denominator under three allocation rules")
say("    last = biotype of the LAST member (what mk_vasa_composition.py does)")
say("    unan = combination keys counted only where every member agrees")
say("    frac = combination-key reads split equally across member biotypes")
say(al.round(4).to_string())
say()
say("  where each dataset's structural reads come from (% of that class's reads)")
say(bk.round(3).to_string(index=False))


# --------------------------------------------------------------------------
# B. duplication: reads vs molecules on BOTH VASA plates
# --------------------------------------------------------------------------
say()
say("=" * 78)
say("B. DUPLICATION -- reads vs molecules, both VASA plates")
say("=" * 78)


def biotype_sums(path, cols, species_filter):
    got, denom = {}, 0.0
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=CHUNK):
        ch.columns = norm_cols(ch.columns)
        if cols is not None:
            ch = ch[cols]
        idx = ch.index.astype(str)
        if species_filter is not None:
            keep = np.isin([vp.species_of(i) for i in idx], species_filter)
            ch, idx = ch[keep], idx[keep]
            if len(ch) == 0:
                continue
        bt = pd.Series([token(i) for i in idx], index=ch.index)
        keep = ~bt.isin(RRNA_PRIMARY)
        ch, bt = ch[keep], bt[keep]
        rs = ch.sum(axis=1)
        denom += float(rs.sum())
        for b, v in rs.groupby(bt.values).sum().items():
            got[b] = got.get(b, 0.0) + float(v)
    return got, denom


dup = {}
for k, base, cols, spf in (("published_vasa", PUB, MOUSE, ("mouse",)),
                           ("own_vasa", OWN, OWNC, None)):
    r, dr = biotype_sums(f"{base}.ReadCounts.tsv", cols, spf)
    t, dt = biotype_sums(f"{base}.TranscriptCounts.tsv", cols, spf)
    dup[k] = dict(reads=r, dr=dr, mols=t, dt=dt)
    say(f"  {k:16s} reads={dr:14,.0f}  molecules={dt:14,.0f}  "
        f"reads/molecule = {dr / dt:6.2f}   per unit: "
        f"{dr / dec[k]['n_units']:12,.0f} reads, {dt / dec[k]['n_units']:12,.0f} molecules")

drows = []
for b in STRUCTURAL + ["ProteinCoding", "lncRNA"]:
    row = dict(biotype=b)
    for k in dup:
        d = dup[k]
        row[f"{k}_pct_reads"] = 100 * d["reads"].get(b, 0.0) / d["dr"]
        row[f"{k}_pct_molecules"] = 100 * d["mols"].get(b, 0.0) / d["dt"]
        row[f"{k}_reads_per_molecule"] = (d["reads"].get(b, 0.0) / d["mols"].get(b, np.nan)
                                          if d["mols"].get(b, 0) else np.nan)
    drows.append(row)
dp = pd.DataFrame(drows).set_index("biotype")
dp["own_over_pub_reads"] = dp.own_vasa_pct_reads / dp.published_vasa_pct_reads
dp["own_over_pub_molecules"] = dp.own_vasa_pct_molecules / dp.published_vasa_pct_molecules
dp.round(6).to_csv(f"{RES}/threeway_reads_vs_molecules.tsv", sep="\t")
say()
say(dp.round(4).to_string())
s_r = {k: 100 * sum(dup[k]["reads"].get(b, 0.0) for b in STRUCTURAL) / dup[k]["dr"] for k in dup}
s_m = {k: 100 * sum(dup[k]["mols"].get(b, 0.0) for b in STRUCTURAL) / dup[k]["dt"] for k in dup}
say()
say(f"  structural, READS    : own {s_r['own_vasa']:7.4f}%  published {s_r['published_vasa']:7.4f}%"
    f"   own/pub = {s_r['own_vasa'] / s_r['published_vasa']:6.2f}x")
say(f"  structural, MOLECULES: own {s_m['own_vasa']:7.4f}%  published {s_m['published_vasa']:7.4f}%"
    f"   own/pub = {s_m['own_vasa'] / s_m['published_vasa']:6.2f}x")


# --------------------------------------------------------------------------
# C. geometry: read length and the realised jS:IN rate, per plate
# --------------------------------------------------------------------------
say()
say("=" * 78)
say("C. GEOMETRY -- aligned read length and realised jS:IN rate (trap 8)")
say("=" * 78)
say("  Deterministic stride over each cell's singlemapper annotated BED (trap 4).")
say("  jS:IN = read contained in the feature, i.e. what stage 6 requires of a")
say("  non-splicing biotype. A LONGER read is HARDER to contain, so if the own")
say("  plate's reads are longer, geometry cannot be what gives it more short RNA.")

STRIDE, CAP = 37, 400_000


def bed_lengths(pattern, cells, tag):
    rows = []
    for cell in cells:
        hits = sorted(glob.glob(pattern.format(cell=cell)))
        if not hits:
            continue
        f = hits[0]
        L, jin, n, feat_short = [], 0, 0, 0
        with gzip.open(f, "rt") as fh:
            for i, line in enumerate(fh):
                if i % STRIDE:
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 12:
                    continue
                try:
                    rs, re_ = int(p[1]), int(p[2])
                    fs, fe = int(p[7]), int(p[8])
                except (ValueError, IndexError):
                    continue
                rl, fl = re_ - rs, fe - fs
                L.append(rl)
                if rs >= fs and re_ <= fe:
                    jin += 1
                if fl < rl:
                    feat_short += 1
                n += 1
                if n >= CAP:
                    break
        if n:
            a = np.array(L)
            rows.append(dict(dataset=tag, unit=cell, n_sampled=n,
                             read_len_mean=a.mean(), read_len_median=float(np.median(a)),
                             read_len_p10=float(np.percentile(a, 10)),
                             read_len_p90=float(np.percentile(a, 90)),
                             pct_jS_IN=100 * jin / n,
                             pct_feature_shorter_than_read=100 * feat_short / n,
                             bed=os.path.basename(f)))
    return rows


geo = []
geo += bed_lengths(f"{W}/data/PM26037/out/cells/ZHA9292A1_{{cell}}_*singlemappers_genes.bed.gz",
                   OWNC, "own_vasa")
random.seed(SEED)
pub_sample = sorted(random.sample(MOUSE, min(12, len(MOUSE))))
geo += bed_lengths(f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3/SRR14783059_{{cell}}_*singlemappers_genes.bed.gz",
                   pub_sample, "published_vasa")
if geo:
    g = pd.DataFrame(geo)
    g.round(4).to_csv(f"{RES}/threeway_readlength_geometry.tsv", sep="\t", index=False)
    say(f"  published-plate cells sampled (seed {SEED}): {', '.join(pub_sample)}")
    say(g.drop(columns=["bed"]).round(2).to_string(index=False))
    say()
    sm = g.groupby("dataset")[["read_len_mean", "read_len_median", "pct_jS_IN",
                               "pct_feature_shorter_than_read"]].median().round(3)
    say("  per-dataset medians across units:")
    say(sm.to_string())
else:
    say("  NO per-cell BEDs matched -- reported as not measured, NOT as zero")


# --------------------------------------------------------------------------
# D. pooling: per-unit structural share
# --------------------------------------------------------------------------
say()
say("=" * 78)
say("D. POOLING -- per-unit structural share (is the pooled number carried by one unit?)")
say("=" * 78)


def per_unit(path, cols, species_filter, tag):
    num = den = None
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=CHUNK):
        ch.columns = norm_cols(ch.columns)
        if cols is not None:
            ch = ch[cols]
        idx = ch.index.astype(str)
        if species_filter is not None:
            keep = np.isin([vp.species_of(i) for i in idx], species_filter)
            ch, idx = ch[keep], idx[keep]
            if len(ch) == 0:
                continue
        bt = np.array([token(i) for i in idx])
        keep = ~np.isin(bt, list(RRNA_PRIMARY))
        ch, bt = ch[keep], bt[keep]
        d = ch.sum()
        n = ch[np.isin(bt, STRUCTURAL)].sum()
        den = d if den is None else den.add(d, fill_value=0)
        num = n if num is None else num.add(n, fill_value=0)
    return pd.DataFrame({"dataset": tag, "unit": den.index, "denominator_reads": den.values,
                         "structural_reads": num.reindex(den.index).values,
                         "structural_pct": (100 * num.reindex(den.index) / den).values})


pu = pd.concat([per_unit(p, c, s, k) for k, p, c, s in srcs], ignore_index=True)
pu.round(6).to_csv(f"{RES}/threeway_structural_per_unit.tsv", sep="\t", index=False)
say(pu.groupby("dataset").structural_pct.describe()[["count", "min", "25%", "50%", "75%", "max"]]
    .round(4).to_string())
say()
say("  own-plate cells, individually:")
say(pu[pu.dataset == "own_vasa"].round(4).to_string(index=False))
say()
say("  Overlap test: does ANY published-plate mouse cell reach the own plate's minimum?")
lo_own = pu[pu.dataset == "own_vasa"].structural_pct.min()
hi_pub = pu[pu.dataset == "published_vasa"].structural_pct.max()
n_over = int((pu[pu.dataset == "published_vasa"].structural_pct >= lo_own).sum())
say(f"    own-plate minimum   = {lo_own:.4f}%")
say(f"    published maximum   = {hi_pub:.4f}%   ({n_over} of "
    f"{(pu.dataset == 'published_vasa').sum()} published cells at or above the own minimum)")

with open(f"{RES}/threeway_mechanism.log", "w") as fh:
    fh.write("\n".join(log) + "\n")
say()
say(f"wrote allocation / duplication / geometry / per-unit tables to {RES}")
