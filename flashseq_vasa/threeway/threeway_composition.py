#!/usr/bin/env python3
"""Three-way biotype composition on the non-rRNA denominator: reads on all sides.

  published VASA-plate (SRR14783059, GRCm38 / Ensembl 99, mouse entries + mouse cells)
  own VASA-plate       (ZHA9292A1,    GRCm39 / Ensembl 116, 12 real cells)
  FLASH-seq            (FSall10,      GRCm39 / Ensembl 116, 10 libraries, both read-length arms)

Everything here is built to match `res/flashseq_vasa/mk_vasa_composition.py` and
`res/flashseq_vasa/mk_shortbiotype.py`, which produced the two-way comparison in
`res/flashseq_vasa/composition_flashseq_vs_vasa.tsv`. The three decisions carried
over verbatim, each already litigated:

  1. UNFILTERED `uniaggGenes_*.ReadCounts.tsv`, never `$OUTDIR/analysis/`. The
     UMI-ceiling filter drops the 8 most abundant small RNAs -- correct for
     molecule counting, wrong here: it shifts ProteinCoding by +19.19 pp.
  2. ReadCounts on ALL sides. VASA counts deduplicated molecules and FLASH-seq
     (smartseq_noUMI) counts reads, so reads are the only unit both protocols
     measure (conventions Rule 4).
  3. rRNA out of the denominator. Real rRNA is measured by alignment at stage 3,
     upstream of these tables; the `rRNA` biotype on GRCm39 is an annotation
     relic, not a measurement (conventions trap 9).

Biotype rule, verbatim from both upstream scripts: the token after the last '_'
of the row label. On a multi-gene combination key ('A_ProteinCoding-B_lncRNA')
that is the LAST member's biotype. Crude, but identical on all sides -- and the
read mass sitting in combination keys, and in keys whose members disagree, is
quantified below rather than assumed small.

Denominator rule, also verbatim: rows whose biotype token is in
{'rRNA', 'Mt_rRNA'} are dropped. The token the pipeline actually writes for
mitochondrial rRNA is 'MtRrna', which that set does not match, so MtRrna stays
in the numerator. Reproduced deliberately so the three-way table is comparable
with the two-way one; the effect of also dropping it is measured
(`denominator_variant` rows) rather than silently corrected (Rule 3).

Usage: threeway_composition.py <RESDIR>
"""
import os
import sys
import hashlib
import numpy as np
import pandas as pd

W = "/nemo/lab/turnerj/working/guangxin/vasaseq"
REF = "/nemo/lab/turnerj/working/guangxin/reference/vasaseq"
sys.path.insert(0, f"{W}/code/I_Gene_expression/vasaplate_check")
import vp_common as vp  # noqa: E402  species_of / classify_fig1d / classify_methods

RES = sys.argv[1]
os.makedirs(RES, exist_ok=True)

PUB_UNIAGG = f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_uniaggGenes_total.ReadCounts.tsv"
PUB_TOTAL = f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_total.ReadCounts.tsv"
PUB_UFI = f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_total.UFICounts.tsv"
PUB_TRNA = f"{W}/data/ref/fastq_vasaplate/vasaplate_out_v3_tRNA.ReadCounts.tsv"
OWN_UNIAGG = f"{W}/data/PM26037/out/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv"
OWN_TOTAL = f"{W}/data/PM26037/out/ZHA9292A1_total.ReadCounts.tsv"
OWN_TRNA = f"{W}/data/PM26037/out/ZHA9292A1_tRNA.ReadCounts.tsv"
FSDIR = "/nemo/lab/turnerj/scratch/zhangg/vasaseq/flashseq_vasa"
BED99 = f"{REF}/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed"
BED116 = f"{REF}/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed"

BLANKS = {"001", "014", "015", "016"}          # own-plate empty wells, verbatim
RRNA_PRIMARY = {"rRNA", "Mt_rRNA"}             # verbatim from the two-way scripts
RRNA_STRICT = {"rRNA", "Mt_rRNA", "MtRrna"}    # what that set was evidently meant to be
# The five classes whose ~90x protocol gap this track exists to test. Defined
# once, here, so every number downstream refers to the same set.
STRUCTURAL = ["MiscRna", "snRNA", "snoRNA", "scaRNA", "ribozyme"]
CHUNK = 40_000
log = []


def say(s=""):
    print(s, flush=True)
    log.append(s)


def biotype_token(labels):
    """Verbatim rule from mk_vasa_composition.py / mk_shortbiotype.py."""
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


def sha(path, n=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(n)
            if not b:
                break
            h.update(b)
    return h.hexdigest()[:16]


# --------------------------------------------------------------------------
# 1. annotation BEDs -> mouse gene id -> biotype, for both releases
# --------------------------------------------------------------------------
def bed_gene_biotypes(path, mouse_only):
    """{gene_id: biotype} from a homemade_IntronExonTrna BED (col 5 = name).

    Name layout is `<ID>_<NAME>_<BIOTYPE>_<exon|intron>`; tRNA rows are
    `<chr>.tRNA<n>-<aa>` with no Ensembl id and are skipped (they carry no
    gene id to intersect on, and no tRNA row reaches the gene tables anyway).
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
            body = name.rsplit("_", 1)[0]          # strip _exon / _intron
            gid = body.split("_", 1)[0]
            bt = body.rsplit("_", 1)[-1]
            got[gid] = bt
    return got


say("=" * 78)
say("ANNOTATION BEDs")
say("=" * 78)
bt99 = bed_gene_biotypes(BED99, mouse_only=True)
bt116 = bed_gene_biotypes(BED116, mouse_only=False)
shared = set(bt99) & set(bt116)
say(f"  Ensembl  99 (GRCm38), mouse gene ids   : {len(bt99):,}")
say(f"  Ensembl 116 (GRCm39), mouse gene ids   : {len(bt116):,}")
say(f"  shared gene ids                        : {len(shared):,}")
say(f"  E99-only  : {len(set(bt99) - set(bt116)):,}    E116-only : {len(set(bt116) - set(bt99)):,}")
relabelled = {g for g in shared if bt99[g] != bt116[g]}
say(f"  shared ids whose BIOTYPE changed 99->116: {len(relabelled):,} "
    f"({100 * len(relabelled) / len(shared):.3f}% of shared)")

pd.DataFrame({
    "gene_id": sorted(shared),
    "biotype_E99": [bt99[g] for g in sorted(shared)],
    "biotype_E116": [bt116[g] for g in sorted(shared)],
}).assign(changed=lambda d: d.biotype_E99 != d.biotype_E116) \
  .to_csv(f"{RES}/threeway_biotype_relabel_E99_E116.tsv", sep="\t", index=False)


# --------------------------------------------------------------------------
# 2. published plate: which cells are mouse?
# --------------------------------------------------------------------------
say()
say("=" * 78)
say("PUBLISHED PLATE -- CELL CALLING (both of the paper's rules)")
say("=" * 78)
h_ufi = m_ufi = gh = gm = all_ufi = None
for ch in pd.read_csv(PUB_UFI, sep="\t", index_col=0, chunksize=CHUNK):
    ch.columns = norm_cols(ch.columns)
    sp = pd.Series([vp.species_of(str(i)) for i in ch.index], index=ch.index).values
    hh, mm = ch[sp == "human"], ch[sp == "mouse"]
    def acc(a, b):
        return b if a is None else a.add(b, fill_value=0)
    h_ufi, m_ufi = acc(h_ufi, hh.sum()), acc(m_ufi, mm.sum())
    gh, gm = acc(gh, (hh > 0).sum()), acc(gm, (mm > 0).sum())
    all_ufi = acc(all_ufi, ch.sum())

lab_f, frac_u = vp.classify_fig1d(h_ufi, m_ufi)
lab_m, frac_g = vp.classify_methods(gh, gm, h_ufi, m_ufi)
n_unambig_ge = int(((h_ufi + m_ufi) >= vp.MIN_UFI).sum())
n_all_ge = int((all_ufi >= vp.MIN_UFI).sum())
say(f"  barcodes: {len(all_ufi)}   >=7,500 UFIs on unambiguous rows: {n_unambig_ge}"
    f"   on all rows: {n_all_ge}")
say(f"  (benchmark BENCHMARK_published_plate.md reports 353 for run 5/v3)")
for nm, lab in (("Fig.1d rule (UFI fraction)", lab_f), ("Methods rule (gene fraction)", lab_m)):
    vc = lab.value_counts()
    say(f"  {nm:32s} mouse={vc.get('mouse', 0):3d} human={vc.get('human', 0):3d} "
        f"mixed={vc.get('mixed', 0):3d} discarded={vc.get('discarded', 0):3d}"
        f"   doublet={vp.doublet_rate(lab):.2f}%")
pd.DataFrame({"call_fig1d": lab_f, "call_methods": lab_m, "human_ufi": h_ufi,
              "mouse_ufi": m_ufi, "genes_human": gh, "genes_mouse": gm,
              "frac_human_ufi": frac_u.round(5), "frac_human_genes": frac_g.round(5)}) \
    .to_csv(f"{RES}/threeway_published_cellcalls.tsv", sep="\t")

MOUSE_F = sorted(lab_f[lab_f == "mouse"].index)
MOUSE_M = sorted(lab_m[lab_m == "mouse"].index)
say(f"  mouse cells used (primary = Fig.1d rule): {len(MOUSE_F)}"
    f"   |  Methods rule would give {len(MOUSE_M)}"
    f"   |  intersection {len(set(MOUSE_F) & set(MOUSE_M))}")


# --------------------------------------------------------------------------
# 3. composition accumulators
# --------------------------------------------------------------------------
class Comp:
    """Accumulates read sums per biotype under several row restrictions."""

    def __init__(self, name, species_filter=None):
        self.name, self.species_filter = name, species_filter
        self.bt = {}            # biotype -> reads          (all kept rows)
        self.bt_shared = {}     # biotype -> reads          (simple rows, shared ids, native label)
        self.bt_relab = {}      # biotype -> reads          (simple rows, shared ids, E116 label)
        self.bt_simple = {}     # biotype -> reads          (simple rows, any id)
        self.n_entries = self.n_detected = 0
        self.reads_all = 0.0        # every row, incl. rRNA-token rows
        self.reads_rrna = 0.0       # RRNA_PRIMARY rows
        self.reads_mtrrna = 0.0     # MtRrna rows (inside the primary denominator)
        self.reads_dropped_species = 0.0   # rows removed by the species filter
        self.reads_combo = 0.0      # combination keys, inside the denominator
        self.reads_combo_multibt = 0.0     # combination keys whose members disagree on biotype
        self.reads_notshared = 0.0  # simple rows whose gene id is not in `shared`

    def add(self, df):
        df = df[[c for c in df.columns]]
        rs = df.sum(axis=1)
        self.reads_all += float(rs.sum())
        idx = df.index.astype(str)

        if self.species_filter is not None:
            spv = np.array([vp.species_of(i) for i in idx])
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

        is_combo = np.array(["-" in i for i in idx])
        self.reads_combo += float(rs.values[is_combo].sum())
        multi = np.array([
            len({p.rsplit("_", 1)[-1] for p in i.split("-")}) > 1 if "-" in i else False
            for i in idx])
        self.reads_combo_multibt += float(rs.values[multi].sum())

        # simple rows only, for the release control
        srs, sidx, sbt = rs[~is_combo], idx[~is_combo], bt[~is_combo]
        if len(srs):
            for b, v in srs.groupby(sbt.values).sum().items():
                self.bt_simple[b] = self.bt_simple.get(b, 0.0) + float(v)
            gid = np.array([i.split("_", 1)[0] for i in sidx])
            insh = np.array([g in shared for g in gid])
            self.reads_notshared += float(srs.values[~insh].sum())
            if insh.any():
                for b, v in srs[insh].groupby(sbt[insh].values).sum().items():
                    self.bt_shared[b] = self.bt_shared.get(b, 0.0) + float(v)
                relab = pd.Series([bt116[g] for g in gid[insh]], index=srs[insh].index)
                for b, v in srs[insh].groupby(relab.values).sum().items():
                    self.bt_relab[b] = self.bt_relab.get(b, 0.0) + float(v)

    def denom(self):
        return sum(self.bt.values())

    def pct(self, which="bt"):
        d = getattr(self, which)
        t = sum(d.values())
        return pd.Series({k: 100.0 * v / t for k, v in d.items()}) if t else pd.Series(dtype=float)


def stream(comp, path, cols=None, chunk=CHUNK):
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=chunk):
        ch.columns = norm_cols(ch.columns)
        if cols is not None:
            ch = ch[cols]
        comp.add(ch)
    return comp


say()
say("=" * 78)
say("READING COUNT TABLES (unfiltered uniaggGenes, ReadCounts, all three sides)")
say("=" * 78)

own_cols = None
head = pd.read_csv(OWN_UNIAGG, sep="\t", index_col=0, nrows=1)
own_cols = [c for c in norm_cols(head.columns) if c not in BLANKS]
assert len(own_cols) == 12, f"expected 12 real own-plate cells, got {len(own_cols)}: {own_cols}"

datasets = {}
datasets["published_vasa"] = stream(Comp("published_vasa", species_filter=("mouse",)),
                                   PUB_UNIAGG, cols=MOUSE_F)
datasets["own_vasa"] = stream(Comp("own_vasa"), OWN_UNIAGG, cols=own_cols)
datasets["flashseq_native"] = stream(Comp("flashseq_native"),
                                     f"{FSDIR}/native/FSall10_native_uniaggGenes_total.ReadCounts.tsv")
datasets["flashseq_vasalen"] = stream(Comp("flashseq_vasalen"),
                                      f"{FSDIR}/vasalen/FSall10_vasalen_uniaggGenes_total.ReadCounts.tsv")
NUNITS = {"published_vasa": len(MOUSE_F), "own_vasa": 12,
          "flashseq_native": 10, "flashseq_vasalen": 10}

for k, cp in datasets.items():
    d = cp.denom()
    say(f"  {k:18s} units={NUNITS[k]:3d}  entries_kept={cp.n_entries:8,d}  "
        f"detected={cp.n_detected:8,d}")
    say(f"  {'':18s} reads_all={cp.reads_all:14,.0f}  rRNA_token={cp.reads_rrna:12,.0f} "
        f"({100 * cp.reads_rrna / cp.reads_all:.4f}%)  denominator={d:14,.0f}")
    say(f"  {'':18s} dropped_by_species_filter={cp.reads_dropped_species:12,.0f} "
        f"({100 * cp.reads_dropped_species / cp.reads_all:.3f}% of reads_all)")
    say(f"  {'':18s} MtRrna inside denominator={100 * cp.reads_mtrrna / d:.4f}%   "
        f"combination keys={100 * cp.reads_combo / d:.2f}%   "
        f"of which biotype-discordant={100 * cp.reads_combo_multibt / d:.2f}%")


# --------------------------------------------------------------------------
# 4. the deliverable table
# --------------------------------------------------------------------------
pct = {k: cp.pct() for k, cp in datasets.items()}
comp = pd.DataFrame(pct).fillna(0.0)
comp.index.name = "biotype"
comp = comp.sort_values("own_vasa", ascending=False)
comp["pub_minus_own_pp"] = (comp.published_vasa - comp.own_vasa)
comp["own_minus_fsnative_pp"] = (comp.own_vasa - comp.flashseq_native)
comp["readlen_effect_pp"] = (comp.flashseq_vasalen - comp.flashseq_native)
with np.errstate(divide="ignore", invalid="ignore"):
    comp["own_over_fsnative_fold"] = (comp.own_vasa / comp.flashseq_native.replace(0, np.nan))
    comp["pub_over_fsnative_fold"] = (comp.published_vasa / comp.flashseq_native.replace(0, np.nan))
    comp["pub_over_own_fold"] = (comp.published_vasa / comp.own_vasa.replace(0, np.nan))
reads = pd.DataFrame({f"reads_{k}": pd.Series(cp.bt) for k, cp in datasets.items()}).fillna(0.0)
comp = comp.join(reads.astype("int64"))
comp.round(6).to_csv(f"{RES}/composition_threeway.tsv", sep="\t")

# Rule: verify, do not transcribe. These two constants are read from the committed
# res/flashseq_vasa/composition_flashseq_vs_vasa.tsv, and reproducing them is what
# proves this script's own-plate side is built the same way as the two-way one.
TWOWAY_OWN_PC = 64.1192          # 'vasa_reads' ProteinCoding, rounded to 4 dp there
TWOWAY_OWN_SNRNA = 7.0327
got_pc = float(comp.loc["ProteinCoding", "own_vasa"])
got_sn = float(comp.loc["snRNA", "own_vasa"])
assert abs(got_pc - TWOWAY_OWN_PC) < 5e-4, \
    f"own-plate ProteinCoding {got_pc:.6f} != committed two-way {TWOWAY_OWN_PC}"
assert abs(got_sn - TWOWAY_OWN_SNRNA) < 5e-4, \
    f"own-plate snRNA {got_sn:.6f} != committed two-way {TWOWAY_OWN_SNRNA}"
say(f"  [check] own-plate side reproduces the committed two-way table: "
    f"ProteinCoding {got_pc:.4f} == {TWOWAY_OWN_PC}, snRNA {got_sn:.4f} == {TWOWAY_OWN_SNRNA}")

meta = pd.DataFrame({
    "dataset": list(datasets),
    "protocol": ["vasa", "vasa", "smartseq_noUMI", "smartseq_noUMI"],
    "genome": ["GRCm38", "GRCm39", "GRCm39", "GRCm39"],
    "annotation": ["Ensembl 99 (human+mouse)", "Ensembl 116", "Ensembl 116", "Ensembl 116"],
    "unit": ["mouse-called cell", "real cell", "library", "library"],
    "n_units": [NUNITS[k] for k in datasets],
    "entries_kept": [datasets[k].n_entries for k in datasets],
    "entries_detected": [datasets[k].n_detected for k in datasets],
    "reads_all_rows": [datasets[k].reads_all for k in datasets],
    "reads_rRNA_token": [datasets[k].reads_rrna for k in datasets],
    "reads_dropped_species_filter": [datasets[k].reads_dropped_species for k in datasets],
    "denominator_nonrRNA_reads": [datasets[k].denom() for k in datasets],
    "pct_MtRrna_in_denominator": [100 * datasets[k].reads_mtrrna / datasets[k].denom() for k in datasets],
    "pct_reads_in_combination_keys": [100 * datasets[k].reads_combo / datasets[k].denom() for k in datasets],
    "pct_reads_in_biotype_discordant_combos": [100 * datasets[k].reads_combo_multibt / datasets[k].denom() for k in datasets],
    "pct_simple_reads_gene_id_not_shared_E99_E116": [
        100 * datasets[k].reads_notshared / sum(datasets[k].bt_simple.values()) for k in datasets],
    "table": [os.path.basename(p) for p in (PUB_UNIAGG, OWN_UNIAGG,
              f"{FSDIR}/native/FSall10_native_uniaggGenes_total.ReadCounts.tsv",
              f"{FSDIR}/vasalen/FSall10_vasalen_uniaggGenes_total.ReadCounts.tsv")],
}).set_index("dataset")
meta.round(6).to_csv(f"{RES}/threeway_denominators.tsv", sep="\t")

say()
say("=" * 78)
say("COMPOSITION ON THE NON-rRNA DENOMINATOR  (% of reads, all sides ReadCounts)")
say("=" * 78)
show = comp[["published_vasa", "own_vasa", "flashseq_native", "flashseq_vasalen",
             "pub_minus_own_pp", "own_minus_fsnative_pp", "readlen_effect_pp"]]
say(show[show.max(axis=1) >= 0.001].round(4).to_string())


# --------------------------------------------------------------------------
# 5. the headline test: structural RNA classes
# --------------------------------------------------------------------------
say()
say("=" * 78)
say("HEADLINE TEST -- structural RNA classes: " + " + ".join(STRUCTURAL))
say("=" * 78)
srows = []
for k in datasets:
    tot = float(comp[k].reindex(STRUCTURAL).fillna(0).sum())
    srows.append(dict(dataset=k, n_units=NUNITS[k], structural_pct_of_nonrRNA=tot,
                      **{b: float(comp[k].get(b, 0.0)) for b in STRUCTURAL}))
st = pd.DataFrame(srows).set_index("dataset")
st["fold_over_flashseq_native"] = st.structural_pct_of_nonrRNA / st.loc["flashseq_native", "structural_pct_of_nonrRNA"]
st.round(6).to_csv(f"{RES}/threeway_structural.tsv", sep="\t")
say(st.round(4).to_string())
say()
say(f"  own VASA / FLASH-seq native  = {st.loc['own_vasa', 'fold_over_flashseq_native']:.1f}x "
    f"(the ~90x the two-way work reported)")
say(f"  published VASA / FS native   = {st.loc['published_vasa', 'fold_over_flashseq_native']:.1f}x")
say(f"  published VASA / own VASA    = "
    f"{st.loc['published_vasa', 'structural_pct_of_nonrRNA'] / st.loc['own_vasa', 'structural_pct_of_nonrRNA']:.3f}x")
say(f"  read-length effect (vasalen - native) = "
    f"{st.loc['flashseq_vasalen', 'structural_pct_of_nonrRNA'] - st.loc['flashseq_native', 'structural_pct_of_nonrRNA']:+.4f} pp")


# --------------------------------------------------------------------------
# 6. release control -- how much of published-vs-own is Ensembl 99 vs 116?
# --------------------------------------------------------------------------
say()
say("=" * 78)
say("RELEASE CONTROL -- E99 vs E116 on published-vs-own gaps")
say("=" * 78)
say("  Track A: simple (single-gene) rows only, gene id present in BOTH BEDs, native labels.")
say("  Track B: same rows, but the PUBLISHED side relabelled with its E116 biotype.")
say("  B - A on the published side isolates biotype REASSIGNMENT; raw - A mixes")
say("  gene-set differences with the loss of combination keys, so it bounds rather")
say("  than measures. Neither track can remove the cell-type difference (mESC vs embryo).")
rc = pd.DataFrame({
    "pub_raw": comp.published_vasa,
    "own_raw": comp.own_vasa,
    "pub_shared_E99": datasets["published_vasa"].pct("bt_shared"),
    "pub_shared_E116": datasets["published_vasa"].pct("bt_relab"),
    "own_shared_E116": datasets["own_vasa"].pct("bt_shared"),
    "fsnative_shared_E116": datasets["flashseq_native"].pct("bt_shared"),
}).fillna(0.0)
rc["gap_raw_pp"] = rc.pub_raw - rc.own_raw
rc["gap_shared_pp"] = rc.pub_shared_E99 - rc.own_shared_E116
rc["gap_shared_relabelled_pp"] = rc.pub_shared_E116 - rc.own_shared_E116
rc["release_label_effect_pp"] = (rc.pub_shared_E116 - rc.pub_shared_E99).abs()
rc["gap_shrinkage_pp"] = rc.gap_raw_pp.abs() - rc.gap_shared_relabelled_pp.abs()
rc["survives_release_control"] = rc.gap_shared_relabelled_pp.abs() > 2 * rc.release_label_effect_pp.clip(lower=1e-9)
rc = rc.sort_values("gap_raw_pp", key=lambda s: s.abs(), ascending=False)
rc.round(6).to_csv(f"{RES}/threeway_release_control.tsv", sep="\t")
say(rc.head(16)[["pub_raw", "own_raw", "gap_raw_pp", "pub_shared_E99", "pub_shared_E116",
                 "own_shared_E116", "gap_shared_relabelled_pp", "release_label_effect_pp",
                 "survives_release_control"]].round(4).to_string())

say()
say("  Self-check on the BED->biotype map: for own/FLASH-seq the count-table token IS")
say("  the E116 label, so `pub_shared_E116`-style relabelling must be a no-op there.")
for k in ("own_vasa", "flashseq_native"):
    a, b = datasets[k].pct("bt_shared"), datasets[k].pct("bt_relab")
    dmax = float((a - b).abs().reindex(a.index.union(b.index)).fillna(0).max())
    say(f"    {k:16s} max |native - E116-relabelled| = {dmax:.6f} pp")
    assert dmax < 0.01, f"{k}: BED biotype map disagrees with its own count-table token by {dmax} pp"

for k in ("published_vasa", "own_vasa", "flashseq_native"):
    cp = datasets[k]
    say(f"  {k:16s} simple-row reads={sum(cp.bt_simple.values()):14,.0f} "
        f"({100 * sum(cp.bt_simple.values()) / cp.denom():5.2f}% of its denominator); "
        f"shared-id subset={sum(cp.bt_shared.values()):14,.0f} "
        f"({100 * sum(cp.bt_shared.values()) / cp.denom():5.2f}%)")


# --------------------------------------------------------------------------
# 7. sensitivity: table choice, cell-calling rule, strict denominator, tRNA
# --------------------------------------------------------------------------
say()
say("=" * 78)
say("SENSITIVITY")
say("=" * 78)
sens = []

# (a) MtRrna also out of the denominator
for k, cp in datasets.items():
    d, d2 = cp.denom(), cp.denom() - cp.reads_mtrrna
    sr = float(sum(cp.bt.get(b, 0.0) for b in STRUCTURAL))
    sens.append(dict(check="denominator_variant_drop_MtRrna", dataset=k,
                     value_pp=100 * sr / d2 - 100 * sr / d,
                     note=f"change in structural % if MtRrna leaves the denominator "
                          f"(denominator {d:,.0f} -> {d2:,.0f})"))

# (b) _total instead of uniaggGenes_total (the two-way FS side used _total)
tot_paths = {"published_vasa": (PUB_TOTAL, MOUSE_F, ("mouse",)),
             "own_vasa": (OWN_TOTAL, own_cols, None),
             "flashseq_native": (f"{FSDIR}/native/FSall10_native_total.ReadCounts.tsv", None, None),
             "flashseq_vasalen": (f"{FSDIR}/vasalen/FSall10_vasalen_total.ReadCounts.tsv", None, None)}
tot_comp = {}
for k, (p, cols, spf) in tot_paths.items():
    if not os.path.exists(p):
        say(f"  [_total] MISSING for {k}: {p}")
        continue
    cp = stream(Comp(k + "_total", species_filter=spf), p, cols=cols)
    tot_comp[k] = cp
    pv = cp.pct()
    sens.append(dict(check="table_total_vs_uniagg", dataset=k,
                     value_pp=float(pv.reindex(STRUCTURAL).fillna(0).sum())
                     - float(comp[k].reindex(STRUCTURAL).fillna(0).sum()),
                     note=f"structural pp, _total minus uniagg; denom {cp.denom():,.0f}"))
if tot_comp:
    pd.DataFrame({k: cp.pct() for k, cp in tot_comp.items()}).fillna(0.0).round(6) \
        .to_csv(f"{RES}/threeway_composition_totaltable.tsv", sep="\t")

# Which table did the committed two-way FLASH-seq side use? mk_shortbiotype.py reads
# `{TAG}_{arm}_total.ReadCounts.tsv` while mk_vasa_composition.py reads
# `uniaggGenes_total` -- so the two-way table may mix them. Settle it by matching
# its published number rather than by reading the code and guessing.
TWOWAY_FS_NATIVE_PC = 84.31994062316626
cand = {"_total": float(tot_comp["flashseq_native"].pct().get("ProteinCoding", np.nan))
        if "flashseq_native" in tot_comp else np.nan,
        "uniaggGenes_total": float(comp.loc["ProteinCoding", "flashseq_native"])}
hits = {k: v for k, v in cand.items() if abs(v - TWOWAY_FS_NATIVE_PC) < 5e-4}
say()
say(f"  [check] committed two-way flashseq_native ProteinCoding = {TWOWAY_FS_NATIVE_PC:.6f}")
for k, v in cand.items():
    say(f"          from {k:20s} -> {v:.6f}   {'MATCH' if k in hits else 'differs'}")
assert hits, ("neither FLASH-seq table reproduces the committed two-way ProteinCoding; "
              f"got {cand}")
say(f"          => the two-way FLASH-seq side was built from {list(hits)[0]}. This script "
    f"uses uniaggGenes_total on ALL sides, so any change is table choice, not protocol.")
sens.append(dict(check="twoway_fs_table_identified", dataset="flashseq_native",
                 value_pp=cand["uniaggGenes_total"] - cand["_total"],
                 note=f"ProteinCoding pp, uniagg minus _total; two-way used {list(hits)[0]}"))

# (c) Methods cell-calling rule instead of Fig.1d
cp_m = stream(Comp("published_methodsrule", species_filter=("mouse",)), PUB_UNIAGG, cols=MOUSE_M)
pm = cp_m.pct()
sens.append(dict(check="cellcall_methods_vs_fig1d", dataset="published_vasa",
                 value_pp=float(pm.reindex(STRUCTURAL).fillna(0).sum())
                 - float(comp["published_vasa"].reindex(STRUCTURAL).fillna(0).sum()),
                 note=f"structural pp, Methods rule ({len(MOUSE_M)} cells) minus Fig.1d "
                      f"({len(MOUSE_F)} cells); denom {cp_m.denom():,.0f}"))
pd.DataFrame({"fig1d": comp["published_vasa"], "methods": pm}).fillna(0.0).round(6) \
    .to_csv(f"{RES}/threeway_published_cellrule_sensitivity.tsv", sep="\t")

# (d) tRNA: absent from the gene tables entirely, lives in its own table
say()
for k, p in (("published_vasa", PUB_TRNA), ("own_vasa", OWN_TRNA),
             ("flashseq_native", f"{FSDIR}/native/FSall10_native_tRNA.ReadCounts.tsv"),
             ("flashseq_vasalen", f"{FSDIR}/vasalen/FSall10_vasalen_tRNA.ReadCounts.tsv")):
    if not os.path.exists(p):
        say(f"  [tRNA] {k:18s} table ABSENT ({os.path.basename(p)}) -- reported as absent, not as zero")
        sens.append(dict(check="tRNA_pct_of_nonrRNA_denominator", dataset=k, value_pp=float("nan"),
                         note="tRNA table absent for this dataset; not estimated"))
        continue
    t = pd.read_csv(p, sep="\t", index_col=0)
    t.columns = norm_cols(t.columns)
    cols = MOUSE_F if k == "published_vasa" else (own_cols if k == "own_vasa" else list(t.columns))
    cols = [c for c in cols if c in t.columns]
    r = float(t[cols].values.sum())
    d = datasets[k].denom()
    say(f"  [tRNA] {k:18s} rows={len(t):5d} reads={r:12,.0f}  "
        f"{100 * r / d:.5f}% of its non-rRNA denominator (separate table, never in the gene tables)")
    sens.append(dict(check="tRNA_pct_of_nonrRNA_denominator", dataset=k, value_pp=100 * r / d,
                     note=f"{os.path.basename(p)}, {len(t)} rows, {len(cols)} units"))

pd.DataFrame(sens).round(6).to_csv(f"{RES}/threeway_sensitivity.tsv", sep="\t", index=False)
say()
say(pd.DataFrame(sens).round(4).to_string(index=False))

prov = pd.DataFrame([
    dict(file=p, sha256_16=sha(p), bytes=os.path.getsize(p)) for p in
    [PUB_UNIAGG, PUB_UFI, OWN_UNIAGG, BED99, BED116,
     f"{FSDIR}/native/FSall10_native_uniaggGenes_total.ReadCounts.tsv",
     f"{FSDIR}/vasalen/FSall10_vasalen_uniaggGenes_total.ReadCounts.tsv"]])
prov.to_csv(f"{RES}/threeway_provenance.tsv", sep="\t", index=False)

with open(f"{RES}/threeway_composition.log", "w") as fh:
    fh.write("\n".join(log) + "\n")
say()
say(f"wrote {RES}/composition_threeway.tsv  ({len(comp)} biotypes x {len(datasets)} datasets)")
