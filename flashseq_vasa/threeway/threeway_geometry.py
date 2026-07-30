#!/usr/bin/env python3
"""C. GEOMETRY: is the published-vs-own structural-RNA gap a read-length artefact?

Conventions trap 8: stage 6 keeps a non-splicing biotype only if the read is
CONTAINED in the feature, and 96-99% of tRNA/rRNA/miRNA/snoRNA features are
shorter than one long read -- so read length alone can suppress a whole biotype.
The two VASA plates differ in read length, so this has to be measured before the
2.25% vs 20.59% structural gap is called protocol or biology.

The first attempt reimplemented containment from BED coordinates and matched
nothing: these per-cell BEDs have NINE fields, not twelve. Corrected here, and
the reimplementation dropped entirely -- field 6 carries the pipeline's OWN
verdict (`CG:<cigar>;nM:<n>;jS:<IN|digit>`), so the measured quantity is stage
6's actual predicate rather than my model of it.

BED layout (verified by inspection, `deal_with_*mappers.sh` output):
  0 chrom  1 read_start  2 read_end  3 read_name+tags  4 strand
  5 feature name `<ID>_<NAME>_<BIOTYPE>_<exon|intron>`
  6 `CG:...;nM:...;jS:...`   7,8 gene-level annotation columns

Published-plate cells are on a mixed GRCh38+GRCm38 reference, so rows are
restricted to ENSMUSG features to match the mouse-only composition.

Also re-verifies the surprising allocation result: three allocation rules gave
IDENTICAL structural shares, which implies no biotype-discordant combination key
anywhere contains a structural biotype. Counted directly here.

Usage: threeway_geometry.py <RESDIR>
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
STRUCTURAL = {"MiscRna", "snRNA", "snoRNA", "scaRNA", "ribozyme"}
BLANKS = {"001", "014", "015", "016"}
SEED = 20260730
STRIDE, CAP = 11, 300_000
CHUNK = 40_000
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


# --- published mouse cells, Fig.1d rule ------------------------------------
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
assert len(MOUSE) == 173, f"expected 173 mouse cells, got {len(MOUSE)}"
OWNC = sorted(set(f"{i:03d}" for i in range(1, 17)) - BLANKS)

say("=" * 78)
say("C. GEOMETRY -- read length and the pipeline's OWN jS:IN verdict")
say("=" * 78)
say(f"  stride={STRIDE} (deterministic, trap 4), cap={CAP:,} rows/cell, seed={SEED}")
say("  jS:IN is read from field 6 as written by stage 5, not recomputed.")
say("  Published rows restricted to ENSMUSG features (mixed-species reference).")


def scan(path, mouse_only):
    """Per-cell read spans and jS:IN rates, split by feature biotype class."""
    n = 0
    spans, spans_s = [], []
    jin = jin_s = n_s = 0
    jin_pc = n_pc = 0
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % STRIDE:
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            feat = p[5]
            if mouse_only and not feat.startswith("ENSMUSG"):
                continue
            try:
                span = int(p[2]) - int(p[1])
            except ValueError:
                continue
            js = ""
            for kv in p[6].split(";"):
                if kv.startswith("jS:"):
                    js = kv[3:]
            is_in = js == "IN"
            body = feat.rsplit("_", 1)[0]          # drop _exon/_intron
            bt = token(body)
            n += 1
            spans.append(span)
            jin += is_in
            if bt in STRUCTURAL:
                n_s += 1
                jin_s += is_in
                spans_s.append(span)
            elif bt == "ProteinCoding":
                n_pc += 1
                jin_pc += is_in
            if n >= CAP:
                break
    if not n:
        return None
    a = np.array(spans)
    s = np.array(spans_s) if spans_s else np.array([np.nan])
    return dict(n_rows=n, span_mean=a.mean(), span_median=float(np.median(a)),
                span_p10=float(np.percentile(a, 10)), span_p90=float(np.percentile(a, 90)),
                pct_jS_IN_all=100 * jin / n,
                n_rows_structural=n_s,
                pct_jS_IN_structural=100 * jin_s / n_s if n_s else np.nan,
                span_median_structural=float(np.median(s)),
                n_rows_proteincoding=n_pc,
                pct_jS_IN_proteincoding=100 * jin_pc / n_pc if n_pc else np.nan,
                pct_rows_structural=100 * n_s / n)


rows = []
for tag, cells, pat, mo in (
        ("own_vasa", OWNC,
         f"{W}/data/PM26037/out/cells/ZHA9292A1_{{c}}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.singlemappers_genes.bed.gz",
         False),
        ("published_vasa", None,
         f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3/SRR14783059_{{c}}_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.singlemappers_genes.bed.gz",
         True)):
    if cells is None:
        random.seed(SEED)
        cells = sorted(random.sample(MOUSE, 12))
        say(f"  published-plate mouse cells sampled (seed {SEED}): {', '.join(cells)}")
    for cl in cells:
        f = pat.format(c=cl)
        if not os.path.exists(f):
            hits = sorted(glob.glob(pat.replace("{c}", cl).replace(
                "cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.", "*")))
            if not hits:
                say(f"    MISSING BED for {tag} {cl}")
                continue
            f = hits[0]
        r = scan(f, mo)
        if r is None:
            say(f"    NO USABLE ROWS for {tag} {cl}")
            continue
        r.update(dataset=tag, unit=cl, bed=os.path.basename(f))
        rows.append(r)

assert rows, "no per-cell BEDs produced usable rows -- geometry NOT measured"
g = pd.DataFrame(rows)
g.round(4).to_csv(f"{RES}/threeway_readlength_geometry.tsv", sep="\t", index=False)
say()
say(g[["dataset", "unit", "n_rows", "span_mean", "span_median", "span_p90",
       "pct_jS_IN_all", "n_rows_structural", "pct_rows_structural",
       "pct_jS_IN_structural", "span_median_structural",
       "pct_jS_IN_proteincoding"]].round(2).to_string(index=False))

say()
say("  per-dataset medians across cells:")
cols = ["span_mean", "span_median", "span_p90", "pct_jS_IN_all", "pct_rows_structural",
        "pct_jS_IN_structural", "span_median_structural", "pct_jS_IN_proteincoding"]
sm = g.groupby("dataset")[cols].median()
say(sm.round(3).to_string())
say()
o, p = sm.loc["own_vasa"], sm.loc["published_vasa"]
say(f"  own-plate reads are {o.span_median / p.span_median:.2f}x the published median span "
    f"({o.span_median:.0f} vs {p.span_median:.0f} nt).")
say("  A LONGER read is HARDER to contain in a short feature, so if the own plate")
say("  has BOTH longer reads AND more structural RNA, read-length geometry works")
say("  AGAINST the observed gap rather than producing it.")
say(f"  jS:IN on structural features: own {o.pct_jS_IN_structural:.2f}% vs "
    f"published {p.pct_jS_IN_structural:.2f}%")
say(f"  structural share of annotated rows: own {o.pct_rows_structural:.3f}% vs "
    f"published {p.pct_rows_structural:.3f}%  "
    f"({o.pct_rows_structural / p.pct_rows_structural:.1f}x)")
say("  That last line is measured upstream of every count-table decision -- no")
say("  denominator, allocation rule, annotation release or filter is involved,")
say("  only which features the reads landed in.")

# --- re-verify the allocation result --------------------------------------
say()
say("=" * 78)
say("A(verify). Do biotype-discordant combination keys ever carry a structural class?")
say("=" * 78)
say("  Three allocation rules gave identical structural shares, which implies no.")
say("  Counted directly, over every row of both VASA plates.")
vrows = []
for tag, path, cols, spf in (
        ("published_vasa", f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_uniaggGenes_total.ReadCounts.tsv",
         MOUSE, ("mouse",)),
        ("own_vasa", f"{W}/data/PM26037/out/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv", OWNC, None)):
    n_disc = n_disc_struct = 0
    r_disc = r_disc_struct = 0.0
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=CHUNK):
        ch.columns = norm_cols(ch.columns)
        ch = ch[[c for c in cols if c in ch.columns]]
        idx = ch.index.astype(str)
        if spf is not None:
            keep = np.isin([vp.species_of(i) for i in idx], spf)
            ch, idx = ch[keep], idx[keep]
            if len(ch) == 0:
                continue
        rs = ch.sum(axis=1).values.astype(float)
        for i, lbl in enumerate(idx):
            if "-" not in lbl:
                continue
            mem = [token(x) for x in lbl.split("-")]
            if len(set(mem)) == 1:
                continue
            n_disc += 1
            r_disc += rs[i]
            if STRUCTURAL & set(mem):
                n_disc_struct += 1
                r_disc_struct += rs[i]
    vrows.append(dict(dataset=tag, discordant_keys=n_disc, discordant_reads=r_disc,
                      discordant_keys_containing_structural=n_disc_struct,
                      discordant_reads_containing_structural=r_disc_struct,
                      pct_of_discordant_reads=100 * r_disc_struct / r_disc if r_disc else np.nan))
    say(f"  {tag:16s} discordant keys={n_disc:8,d} carrying {r_disc:14,.0f} reads; "
        f"of these {n_disc_struct:,d} keys / {r_disc_struct:,.0f} reads involve a structural class "
        f"({100 * r_disc_struct / r_disc if r_disc else float('nan'):.4f}%)")
pd.DataFrame(vrows).round(6).to_csv(f"{RES}/threeway_discordant_combos.tsv", sep="\t", index=False)

with open(f"{RES}/threeway_geometry.log", "w") as fh:
    fh.write("\n".join(log) + "\n")
say()
say(f"wrote {RES}/threeway_readlength_geometry.tsv and threeway_discordant_combos.tsv")
