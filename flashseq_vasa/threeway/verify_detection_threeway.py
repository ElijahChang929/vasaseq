#!/usr/bin/env python3
"""Re-derive and ASSERT every number reported for the three-way detection track.

Project rule: no number enters a report unless a committed script re-derives it
from the source table and asserts it. Each CLAIM below is the value as written in
the writeup; the assertion fails loudly if the table says otherwise. Run after
mk_detection_threeway.py.

    verify_detection_threeway.py [RES_DIR]
"""
import sys

import numpy as np
import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else \
    "/nemo/lab/turnerj/working/guangxin/vasaseq/res/threeway"

det = pd.read_csv("%s/detection_threeway.tsv" % RES, sep="\t")
trk = pd.read_csv("%s/detection_threeway_tracks.tsv" % RES, sep="\t")
cross = pd.read_csv("%s/detection_threeway_crossings.tsv" % RES, sep="\t")
uni = pd.read_csv("%s/annotation_universe.tsv" % RES, sep="\t")
scope = pd.read_csv("%s/detection_threeway_scope.tsv" % RES, sep="\t")
conv = pd.read_csv("%s/detection_threeway_convergence.tsv" % RES, sep="\t")
cells = pd.read_csv("%s/published_cell_species.tsv" % RES, sep="\t")
eff = pd.read_csv("%s/detection_threeway_tableeffect.tsv" % RES, sep="\t")

U = dict(zip(uni.quantity, uni.n_genes))
hl = (det.table == "uniagg") & (det.universe == "shared")
tk = trk[(trk.table == "uniagg") & (trk.universe == "shared")]
piv = tk.pivot_table(index="depth", columns="track", values="median_genes")
n = 0


def chk(label, got, want, tol=0.0):
    global n
    ok = abs(float(got) - float(want)) <= tol if isinstance(want, (int, float)) else got == want
    assert ok, "FAILED %s: table says %r, claim says %r" % (label, got, want)
    n += 1
    print("  ok  %-62s %s" % (label, got))


print("=== annotation universe (E99 vs E116, mouse) ===")
chk("protein-coding, Ensembl 99", U["protein-coding mouse genes, Ensembl 99"], 21933)
chk("protein-coding, Ensembl 116", U["protein-coding mouse genes, Ensembl 116"], 21818)
chk("protein-coding in BOTH (shared universe)",
    U["protein-coding in BOTH releases (shared universe)"], 21513)
chk("PC in 99 only", U["protein-coding in 99 only"], 420)
chk("PC in 116 only", U["protein-coding in 116 only"], 305)
chk("shared gene ids, any biotype", U["mouse gene ids shared, any biotype"], 51765)
sh = U["protein-coding in BOTH releases (shared universe)"]
chk("shared as %% of E99 PC", round(100 * sh / U["protein-coding mouse genes, Ensembl 99"], 2), 98.09)
chk("shared as %% of E116 PC", round(100 * sh / U["protein-coding mouse genes, Ensembl 116"], 2), 98.60)

print("\n=== published-plate species calls (384 wells) ===")
chk("mouse cells, Fig.1d rule", int((cells.label_fig1d == "mouse").sum()), 173)
chk("mouse cells, Methods rule", int((cells.label_methods == "mouse").sum()), 144)
chk("human cells, Fig.1d rule", int((cells.label_fig1d == "human").sum()), 178)
chk("discarded (<7,500 UFI), Fig.1d rule", int((cells.label_fig1d == "discarded").sum()), 31)
chk("Methods-rule mouse set is a subset of Fig.1d mouse set",
    bool(set(cells.well[cells.label_methods == "mouse"])
         <= set(cells.well[cells.label_fig1d == "mouse"])), True)

print("\n=== native in-scope depth ranges (uniagg, shared universe) ===")
PUBM = hl & (det.dataset == "VASA_published") & (det.label_fig1d == "mouse")
OWN = hl & (det.dataset == "VASA_own")
FSV = hl & (det.arm == "vasalen") & (det.qc_verdict != "exclude")
pub_d, own_d, fs_d = det[PUBM].total_reads, det[OWN].total_reads, det[FSV].total_reads
chk("published mouse cells, min depth", int(pub_d.min()), 16912)
chk("published mouse cells, max depth", int(pub_d.max()), 800970)
chk("published mouse cells, median depth", int(pub_d.median()), 175895)
chk("own plate, min depth", int(own_d.min()), 824608)
chk("own plate, max depth", int(own_d.max()), 5312077)
chk("own plate, median depth", int(own_d.median()), 2384584)
chk("FLASH-seq trimmed, min depth", int(fs_d.min()), 17068026)
chk("FLASH-seq trimmed, max depth", int(fs_d.max()), 23345354)
chk("published and own depth ranges are DISJOINT", bool(pub_d.max() < own_d.min()), True)
chk("own plate is deeper than published by (x, medians)",
    round(own_d.median() / pub_d.median(), 1), 13.6)
chk("FLASH-seq is deeper than own plate by (x, medians)",
    round(fs_d.median() / own_d.median(), 1), 7.9)

print("\n=== depth support ===")
chk("published: deepest all-units rung",
    int(tk[tk.track == "VASA published (mouse cells)"].ref_depth.iloc[0]), 10000)
chk("own plate: deepest all-units rung",
    int(tk[tk.track == "VASA own plate"].ref_depth.iloc[0]), 500000)
chk("published cells reaching 100k", int((pub_d >= 1e5).sum()), 127)
chk("published cells reaching 500k", int((pub_d >= 5e5).sum()), 17)
pubsel = tk[(tk.track == "VASA published (mouse cells)") & (tk.depth == 500000)]
chk("published selection bias at 500k (genes)",
    float(pubsel.selection_bias_genes.iloc[0]), 179.6, tol=0.05)
chk("published selection bias at 100k (genes)",
    float(tk[(tk.track == "VASA published (mouse cells)")
             & (tk.depth == 100000)].selection_bias_genes.iloc[0]), 48.9, tol=0.05)

print("\n=== headline medians, uniagg / shared universe ===")
CLAIM = {  # depth: (published, own, FS native, FS trimmed, FS 30 pg)
    10000:  (4066.4,  4588.0,  4139.1,  4182.7,  3875.0),
    100000: (8305.1, 10780.1,  9918.8, 10000.2,  8099.8),
    200000: (9407.6, 12191.2, 11246.8, 11350.9,  8894.3),
    500000: (10860.3, 13746.0, 12688.1, 12845.2, 9899.9),
}
ORDER = ["VASA published (mouse cells)", "VASA own plate", "FLASH-seq native",
         "FLASH-seq VASA-trimmed", "FLASH-seq 30 pg (trimmed)"]
for g, vals in CLAIM.items():
    for t, v in zip(ORDER, vals):
        chk("median genes @ %s, %s" % (format(g, ","), t), piv.loc[g, t], v, tol=0.05)

REF = "VASA published (mouse cells)"
print("\n=== gaps (published plate as reference) ===")
for g in (10000, 100000, 500000):
    chk("own minus published @ %s" % format(g, ","),
        round(piv.loc[g, "VASA own plate"] - piv.loc[g, REF], 1),
        {10000: 521.6, 100000: 2475.0, 500000: 2885.7}[g], tol=0.05)
    chk("FS-trimmed minus published @ %s" % format(g, ","),
        round(piv.loc[g, "FLASH-seq VASA-trimmed"] - piv.loc[g, REF], 1),
        {10000: 116.3, 100000: 1695.1, 500000: 1984.9}[g], tol=0.05)

print("\n=== own-plate lead over FLASH-seq, uniagg on BOTH sides ===")
# The known result (796-1,226 vs trimmed, 2,778-4,255 vs 30 pg) was computed with
# uniagg for the own plate and total for FLASH-seq. Re-measured with uniagg on both
# sides the lead SHRINKS but does not vanish or change sign.
lead = (piv["VASA own plate"] - piv["FLASH-seq VASA-trimmed"]).dropna()
lead = lead[(lead.index >= 1e5) & (lead.index <= 5e6)]
chk("own vs FS-trimmed, min over 100k-5M", round(float(lead.min()), 1), 779.9, tol=0.05)
chk("own vs FS-trimmed, max over 100k-5M", round(float(lead.max()), 1), 1129.6, tol=0.05)
chk("own vs FS-trimmed, positive at every rung 100k-5M", bool((lead > 0).all()), True)
lead30 = (piv["VASA own plate"] - piv["FLASH-seq 30 pg (trimmed)"]).dropna()
lead30 = lead30[(lead30.index >= 1e5) & (lead30.index <= 5e6)]
chk("own vs FS 30 pg, min over 100k-5M", round(float(lead30.min()), 1), 2680.3, tol=0.05)
chk("own vs FS 30 pg, max over 100k-5M", round(float(lead30.max()), 1), 4090.1, tol=0.05)

print("\n=== crossings ===")
# NB the "no pair crosses" statement is TRUE ONLY of the five three-way pairs
# below. Two FLASH-seq-internal pairs DO cross, beyond 13M reads, where the 30 pg
# rung catches up with the higher-input rungs. Stating it unqualified was a real
# error caught in review; the qualified form is asserted here and the crossings
# are asserted positively rather than being left out.
cs = cross[(cross.table == "uniagg") & (cross.universe == "shared")]
P30 = "FLASH-seq 30 pg (trimmed)"
THREEWAY = [(REF, "VASA own plate"), (REF, "FLASH-seq native"),
            (REF, "FLASH-seq VASA-trimmed"),
            ("VASA own plate", "FLASH-seq native"),
            ("VASA own plate", "FLASH-seq VASA-trimmed")]
for a, b in THREEWAY:
    row = cs[(cs.track_a == a) & (cs.track_b == b)]
    assert len(row) == 1, "no crossing row for %s vs %s" % (a, b)
    chk("%s vs %s: crosses?" % (a[:22], b[:22]), row.crosses.iloc[0], "no")
chk("published never leads own plate at any shared rung",
    cs[(cs.track_a == REF) & (cs.track_b == "VASA own plate")].a_leads_throughout.iloc[0], "no")
chk("published DOES lead FLASH-seq 30 pg throughout",
    cs[(cs.track_a == REF) & (cs.track_b == P30)].a_leads_throughout.iloc[0], "yes")
xr = cs[cs.crosses == "yes"]
chk("number of crossing pairs anywhere in the headline slice", len(xr), 2)
chk("every crossing pair is FLASH-seq-internal", sorted(set(xr.track_b)), [P30])
chk("crossing partners are the two higher-input FLASH-seq arms",
    sorted(set(xr.track_a)), ["FLASH-seq VASA-trimmed", "FLASH-seq native"])
chk("native vs 30 pg crossing depth",
    int(xr[xr.track_a == "FLASH-seq native"].crossing_depth.iloc[0]), 13370551, tol=1)
chk("trimmed vs 30 pg crossing depth",
    int(xr[xr.track_a == "FLASH-seq VASA-trimmed"].crossing_depth.iloc[0]), 19460464, tol=1)
chk("both crossings lie beyond 13M reads, far above any VASA cell's depth",
    bool(float(xr.crossing_depth.min()) > 13e6), True)
chk("30 pg is the LOWEST track at every rung where the three-way pairs are compared",
    bool((piv.loc[piv[REF].dropna().index, P30]
          < piv.loc[piv[REF].dropna().index, REF]).all()), True)

print("\n=== full ordering of ALL five tracks (not just the three-way pairs) ===")
# Review caught "own > FLASH-seq > published" in the commit message and in two
# figure titles. It is false: the FLASH-seq 30 pg rung sits BELOW the published
# plate at every rung, so "FLASH-seq" is not a single band that outranks the
# published plate. The full ordering is asserted here so no summary can restate
# the collapsed version.
FIVE = [t for t in set(tk.track) if "Methods rule" not in t]
pv5 = tk[tk.track.isin(FIVE)].pivot_table(index="depth", columns="track",
                                          values="median_genes")
rungs = pv5[REF].dropna().index          # rungs where all three datasets compare
sub = pv5.loc[rungs]
EXPECTED = ["VASA own plate", "FLASH-seq VASA-trimmed", "FLASH-seq native",
            REF, P30]
for d in rungs:
    chk("ordering @ %s" % format(int(d), ","),
        list(sub.loc[d].sort_values(ascending=False).index), EXPECTED)
chk("the published plate is NOT the lowest track -- 30 pg is below it everywhere",
    sorted(set(sub.idxmin(axis=1))), [P30])
chk("'own > FLASH-seq > published' is FALSE as stated (30 pg is a FLASH-seq arm "
    "below published)", bool((sub[P30] < sub[REF]).all()), True)

print("\n=== convergence / separation ===")
chk("three-way spread @ 10k (genes)",
    float(conv[conv.depth == 10000].spread_genes.iloc[0]), 521.6, tol=0.05)
chk("three-way spread @ 500k (genes)",
    float(conv[conv.depth == 500000].spread_genes.iloc[0]), 2885.7, tol=0.05)
chk("spread grows monotonically over 10k-500k",
    bool((conv.sort_values("depth").spread_genes.diff().dropna() > 0).all()), True)
chk("own plate is highest at every supported rung",
    sorted(set(conv.highest)), ["VASA own plate"])
# NB scoped to the THREE tracks in the convergence table (one per dataset,
# FLASH-seq represented by its VASA-trimmed arm). Across all five tracks the
# lowest is the 30 pg rung, not the published plate -- see the ordering block.
chk("published is lowest OF THE THREE compared datasets at every supported rung",
    sorted(set(conv.lowest)), ["VASA published (mouse cells)"])

print("\n=== annotation-release effect on the headline ===")
al = trk[(trk.table == "uniagg") & (trk.universe == "all")].pivot_table(
    index="depth", columns="track", values="median_genes")
delta = (al - piv).abs()
chk("largest unrestricted-minus-shared shift, any track/rung (genes)",
    round(float(delta.max().max()), 1), 107.9, tol=0.05)
chk("that shift as %% of the smallest three-way gap (521.6 @ 10k)",
    round(100 * float(delta.max().max()) / 521.6, 1), 20.7, tol=0.15)
chk("published-plate track moves least under the release control (genes)",
    round(float(delta[REF].max()), 1), 38.4, tol=0.05)
chk("release effect never reverses a pair ordering",
    bool(((al["VASA own plate"] > al["FLASH-seq VASA-trimmed"])
          == (piv["VASA own plate"] > piv["FLASH-seq VASA-trimmed"])).all()), True)

print("\n=== table-type confound in the comparison being extended ===")
for ds, arm, want in [("VASA_own", "-", 4.11), ("FLASH-seq", "native", 7.65),
                      ("FLASH-seq", "vasalen", 9.27), ("VASA_published", "-", 4.59)]:
    r = eff[(eff.dataset == ds) & (eff.arm == arm)]
    chk("read uplift uniagg vs total, %s %s (%%)" % (ds, arm),
        float(r.reads_uplift_pct.iloc[0]), want, tol=0.005)

print("\n=== does the published plate's cell rule change the answer? ===")
# The paper gives two mutually inconsistent doublet rules (UFI purity in the
# Fig.1d caption, gene purity in the Methods). They call 173 vs 144 mouse cells.
# If the conclusion depended on that choice, it would not be a conclusion.
MR = "VASA published (mouse cells, Methods rule)"
both = tk[tk.track.isin([REF, MR])].pivot_table(
    index="depth", columns="track", values="median_genes").dropna()
chk("both cell rules supported at the same 6 rungs", len(both), 6)
chk("Fig.1d minus Methods at 10k (genes)",
    round(float(both.loc[10000, REF] - both.loc[10000, MR]), 1), 16.6, tol=0.05)
chk("largest median difference between the two rules, any rung (genes)",
    round(float((both[REF] - both[MR]).abs().max()), 1), 94.0, tol=0.05)
chk("that maximum falls at 500k, the most depth-selected rung",
    int((both[REF] - both[MR]).abs().idxmax()), 500000)
chk("the Fig.1d rule reads HIGHER at every rung (its 29 extra cells are not "
    "systematically worse)", bool((both[REF] - both[MR] > 0).all()), True)
# the point of the sensitivity check: is the rule choice small against the gap?
gap500 = round(float(piv.loc[500000, "VASA own plate"] - piv.loc[500000, REF]), 1)
chk("rule effect at 500k as %% of the own-vs-published gap there",
    round(100 * 94.0 / gap500, 1), 3.3, tol=0.15)
chk("own plate still leads the published plate under the Methods rule at every rung",
    bool((piv.loc[both.index, "VASA own plate"] > both[MR]).all()), True)

print("\n=== scope / denominators ===")
s = scope[(scope.table == "uniagg")]
chk("published: in-scope entries (mouse PC single-gene)",
    int(s[s.dataset == "VASA_published"].rows_in_scope.iloc[0]), 19907)
chk("own plate: in-scope entries",
    int(s[s.dataset == "VASA_own"].rows_in_scope.iloc[0]), 19918)
chk("no dataset loses a single-gene PC entry to the hyphen rule",
    int(s.rows_pc_single_gene_lost_to_hyphen_rule.sum()), 0)
chk("no duplicate Ensembl ids within any in-scope table",
    int(s.duplicate_ens_ids_in_scope.sum()), 0)

print("\nALL %d CLAIMS VERIFIED against the tables in %s" % (n, RES))
