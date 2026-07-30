#!/usr/bin/env python3
"""
annotation_release_reid_test.py -- two targeted tests that turn suggestive
patterns from the release-effect run into measurements.

TEST 1: is the sncRNA gene-set turnover the SAME LOCI RE-ISSUED under new
        stable ids, or genuinely different annotation?

  The pattern: several sncRNA classes have IDENTICAL totals in both releases
  but nearly disjoint id sets (MiscRna 562 vs 562, 32 shared ids; rRNA 354 vs
  354, 13 shared). Symbol bridging recovers only ~30%, so symbols do not settle
  it. Coordinates cannot settle it either -- the two releases are on different
  assemblies (GRCm38 vs GRCm39) and nothing here lifts over.

  What CAN be tested from the BED alone: if a class's loci were re-issued
  unchanged, the MULTISET of feature lengths in the only-E99 set and the
  only-E116 set should be near-identical (same loci, same lengths, new ids).
  If the annotation genuinely changed, the length distributions should differ.
  This is a real test with a real null -- and it is reported as evidence about
  length-distribution identity, which is NOT the same thing as proof of locus
  identity. A shuffled-label control gives the scale of agreement expected by
  chance for a distribution of that shape.

TEST 2: WHICH genes carry the read loss when the composition is restricted to
        the shared set? The own plate loses 7.016% of single-gene reads, the
        published plate 0.935% -- a 7.5x asymmetry that needs a named cause,
        not a hand-wave.

Read-only. Writes annotation_release_reid.tsv and annotation_release_loss_top.tsv.
"""
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import annotation_release_effect as A

RNG = np.random.default_rng(20260730)   # seed fixed for the shuffle control


def log(*a):
    print(*a, flush=True)


def gene_len(g):
    """Total exonic length of a gene record (convention-corrected upstream)."""
    return sum(g["exon_len"])


def main():
    OUT = sys.argv[1] if len(sys.argv) > 1 else f"{A.W}/res/threeway"
    out = []

    def rec(metric, scope, value, unit, denominator, note=""):
        out.append({"metric": metric, "scope": scope, "value": value,
                    "unit": unit, "denominator": denominator, "note": note})

    b99 = A.parse_bed(A.BED99, mouse_only=True)
    b116 = A.parse_bed(A.BED116, mouse_only=False)
    g99, g116 = b99["genes"], b116["genes"]
    i99, i116 = set(g99), set(g116)
    shared, o99, o116 = i99 & i116, i99 - i116, i116 - i99

    # === TEST 1 ==============================================================
    log("=" * 74)
    log("TEST 1: length-multiset identity of the only-E99 vs only-E116 sets")
    log("=" * 74)
    log(f"{'class':<12} {'n99':>5} {'n116':>5} {'med99':>7} {'med116':>7} "
        f"{'exact_match_pct':>15} {'shuffled_ctrl':>13}")
    for bt in sorted(A.SNCRNA_BIOTYPES | {"lncRNA", "ProteinCoding"}):
        a = sorted(gene_len(g99[g]) for g in o99 if g99[g]["biotype"] == bt)
        b = sorted(gene_len(g116[g]) for g in o116 if g116[g]["biotype"] == bt)
        if len(a) < 10 or len(b) < 10:
            continue
        # multiset intersection size / smaller set: how much of the smaller
        # only-set is length-for-length reproduced in the other only-set
        ca, cb = Counter(a), Counter(b)
        inter = sum((ca & cb).values())
        pct = 100.0 * inter / min(len(a), len(b))
        # control: same-size samples drawn from the OTHER release's whole
        # same-biotype length pool, so the null keeps the length distribution
        # shape but breaks any locus correspondence
        pool = [gene_len(g116[g]) for g in i116 if g116[g]["biotype"] == bt]
        ctrl = []
        for _ in range(20):
            samp = Counter(RNG.choice(pool, size=len(b), replace=False))
            ctrl.append(100.0 * sum((ca & samp).values()) / min(len(a), len(b)))
        cm = float(np.mean(ctrl))
        log(f"{bt:<12} {len(a):>5} {len(b):>5} {np.median(a):>7.0f} "
            f"{np.median(b):>7.0f} {pct:>14.2f}% {cm:>12.2f}%")
        rec("only_set_length_multiset_match", bt, round(pct, 3), "%",
            f"min(only99={len(a)}, only116={len(b)}) genes of this biotype",
            f"shuffled-label control {cm:.2f}%; high match with near-disjoint ids "
            f"is consistent with the same loci re-issued under new stable ids, "
            f"but is NOT proof of locus identity (no liftover was done)")
        rec("only_set_length_multiset_match_control", bt, round(cm, 3), "%",
            f"20 shuffles, same biotype pool, seed 20260730",
            "expected agreement by chance for a distribution of this shape")

    # id blocks: are the only-E116 ncRNA ids in a later allocation block?
    log("")
    log("=== id numeric block of the only-sets (Ensembl allocates in blocks) ===")
    for nm, S, G in (("only99", o99, g99), ("only116", o116, g116)):
        snc = [int(g[7:]) for g in S if G[g]["biotype"] in A.SNCRNA_BIOTYPES]
        if snc:
            log(f"  {nm:8} sncRNA ids: n={len(snc)} min={min(snc)} "
                f"median={int(np.median(snc))} max={max(snc)}")
            rec("sncRNA_only_set_id_block", nm,
                f"min={min(snc)};median={int(np.median(snc))};max={max(snc)}",
                "numeric id", f"{len(snc)} sncRNA genes in {nm}",
                "a disjoint, later id block indicates re-issued ids rather than "
                "genes present in one release only")

    # === TEST 2 ==============================================================
    log("")
    log("=" * 74)
    log("TEST 2: which genes carry the restriction read loss on each plate")
    log("=" * 74)

    def id_of(lab):
        if "-" in lab or "_" not in lab:
            return None
        return lab.split("_")[0]

    hdr = open(A.OWN_READ).readline().rstrip("\n").split("\t")[1:]
    keep = {A.norm_col(cn) for cn in hdr} - A.OWN_BLANKS
    own, _, _ = A.stream_row_totals(A.OWN_READ, keep_wells=keep)
    hs, ms = A.stream_percell_species(A.PUB_UFI)
    labv, _ = A.classify_fig1d(hs, ms)
    mesc = set(labv.index[labv == "mouse"])
    pub, _, _ = A.stream_row_totals(A.PUB_READ, keep_wells=mesc)

    toprows = []
    for plate, tot, G, univ in (("own_E116", own, g116, i116),
                                ("published_E99", pub, g99, i99)):
        lost = {}
        denom = 0.0
        for lab2, v in tot.items():
            if A.is_trna_row(lab2):
                continue
            gid = id_of(lab2)
            if gid is None or not gid.startswith("ENSMUSG") or gid not in G:
                continue
            denom += v
            if gid not in shared:
                lost[lab2] = v
        tl = sum(lost.values())
        log(f"\n  {plate}: {tl:,.0f} of {denom:,.0f} single-gene reads "
            f"({100*tl/denom:.3f}%) sit on release-only genes")
        log(f"  {'reads':>12} {'pct_of_loss':>12}  entry")
        for lab2, v in sorted(lost.items(), key=lambda kv: -kv[1])[:15]:
            gid = id_of(lab2)
            log(f"  {v:>12,.0f} {100*v/tl:>11.2f}%  {lab2}")
            toprows.append({"plate": plate, "entry": lab2,
                            "biotype": G[gid]["biotype"], "symbol": G[gid]["symbol"],
                            "reads": v, "pct_of_restriction_loss": round(100 * v / tl, 4),
                            "pct_of_plate_single_gene_reads": round(100 * v / denom, 4)})
        # how concentrated is the loss?
        vals = sorted(lost.values(), reverse=True)
        for k in (1, 10, 50):
            if len(vals) >= k:
                rec(f"restriction_loss_top{k}_share", plate,
                    round(100.0 * sum(vals[:k]) / tl, 3), "%",
                    f"{tl:,.0f} reads lost to restriction",
                    "concentration of the loss; if a handful of genes carry it, the "
                    "release effect on composition is a few-locus effect not a diffuse one")
        rec("restriction_loss_n_entries", plate, len(lost), "entries",
            "single-gene Ensembl rows with any read on that plate")

    pd.DataFrame(toprows).to_csv(f"{OUT}/annotation_release_loss_top.tsv",
                                 sep="\t", index=False)
    pd.DataFrame(out).to_csv(f"{OUT}/annotation_release_reid.tsv", sep="\t", index=False)
    log("")
    log(f"wrote {OUT}/annotation_release_reid.tsv")
    log(f"wrote {OUT}/annotation_release_loss_top.tsv")


if __name__ == "__main__":
    main()
