#!/usr/bin/env python3
"""
annotation_release_bridge.py -- bracket the release effect from BOTH sides.

WHY THIS IS NECESSARY (and why the id-only number alone would be wrong)
----------------------------------------------------------------------
annotation_release_effect.py restricts both plates to the gene set they share
BY ENSEMBL ID, and attributes the gap that closes to the release. The re-id test
(annotation_release_reid_test.py) showed that this OVERSTATES the release effect,
and showed it with a specific, checkable observation:

  Rn7sk carries 669,205 reads on the own plate as ENSMUSG00002076161 and 4,199
  reads on the published plate as ENSMUSG00000065037. Same gene. Different
  stable id. Under an id-keyed restriction, BOTH copies are discarded as
  "release-only", so a gene that both releases annotate is counted as a release
  difference. The same happens to Ptprg (ENSMUSG00000121513 vs
  ENSMUSG00000021745), Gm25360, Gm22988 and Rn7s2.

  The length-multiset test made the mechanism concrete: for snRNA, snoRNA and
  scaRNA the only-E99 and only-E116 sets have IDENTICAL length multisets (100.00%,
  against shuffled controls of 72.25%, 56.97%, 51.18%), and their numeric id
  blocks are disjoint (only-E99 sncRNA ids span 64,393-106,670; only-E116 span
  118,674-2,076,992). That is re-issued ids for the same loci, not annotation
  gain or loss.

So the id-keyed restriction gives an UPPER bound on the release effect. This
script adds the complementary LOWER bound by additionally bridging genes that
carry the SAME SYMBOL in the two releases. Symbol bridging is imperfect in the
opposite direction -- a placeholder symbol (Gm*, *Rik) can be reused for a
different locus, so bridging on it can merge genes that are not the same -- and
that is exactly why it bounds from the other side.

The honest answer is the interval, and both endpoints are reported.

    upper bound (id-only bridge)      : release effect at most this
    lower bound (id + symbol bridge)  : release effect at least this

A third variant bridges only on INFORMATIVE symbols (excluding Gm*/*Rik/ENSMUSG
placeholders) as a sensitivity check on the lower bound.

Read-only. Writes annotation_release_bracket.tsv and
annotation_release_bracket_composition.tsv.
"""
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import annotation_release_effect as A


def log(*a):
    print(*a, flush=True)


def informative(sym):
    return not (sym.startswith("Gm") or sym.endswith("Rik")
                or sym.startswith("ENSMUSG"))


def build_bridge(g99, g116, only99, only116, require_informative, require_unique=True):
    """Map E99 id -> E116 id for genes bridged by symbol.

    require_unique: only bridge when the symbol names exactly ONE gene in each
    only-set. A symbol appearing twice on one side is ambiguous and bridging it
    would silently pick one arbitrarily.
    """
    s99, s116 = defaultdict(list), defaultdict(list)
    for gid in only99:
        s = g99[gid]["symbol"]
        if require_informative and not informative(s):
            continue
        s99[s].append(gid)
    for gid in only116:
        s = g116[gid]["symbol"]
        if require_informative and not informative(s):
            continue
        s116[s].append(gid)
    bridge = {}
    ambiguous = 0
    for s in set(s99) & set(s116):
        if require_unique and (len(s99[s]) != 1 or len(s116[s]) != 1):
            ambiguous += len(s99[s])
            continue
        bridge[s99[s][0]] = s116[s][0]
    return bridge, ambiguous


def main():
    OUT = sys.argv[1] if len(sys.argv) > 1 else f"{A.W}/res/threeway"
    out, comp_rows = [], []

    def rec(metric, scope, value, unit, denominator, note=""):
        out.append({"metric": metric, "scope": scope, "value": value,
                    "unit": unit, "denominator": denominator, "note": note})

    b99 = A.parse_bed(A.BED99, mouse_only=True)
    b116 = A.parse_bed(A.BED116, mouse_only=False)
    g99, g116 = b99["genes"], b116["genes"]
    i99, i116 = set(g99), set(g116)
    shared, o99, o116 = i99 & i116, i99 - i116, i116 - i99

    br_all, amb_all = build_bridge(g99, g116, o99, o116, require_informative=False)
    br_inf, amb_inf = build_bridge(g99, g116, o99, o116, require_informative=True)
    log(f"symbol bridge (any symbol)        : {len(br_all)} pairs "
        f"({amb_all} genes skipped as ambiguous)")
    log(f"symbol bridge (informative only)  : {len(br_inf)} pairs "
        f"({amb_inf} skipped as ambiguous)")
    rec("symbol_bridge_pairs_any", "E99-only <-> E116-only", len(br_all), "pairs",
        f"{len(o99)} E99-only and {len(o116)} E116-only genes",
        "1:1 symbol matches; ambiguous many-to-one symbols excluded")
    rec("symbol_bridge_pairs_informative", "E99-only <-> E116-only", len(br_inf), "pairs",
        f"{len(o99)} E99-only and {len(o116)} E116-only genes",
        "excluding Gm*/*Rik/ENSMUSG placeholder symbols")

    # The three universes, expressed as the set of ids ACCEPTED on each side.
    universes = {
        "id_only": (shared, shared),
        "id_plus_symbol_any": (shared | set(br_all), shared | set(br_all.values())),
        "id_plus_symbol_informative": (shared | set(br_inf), shared | set(br_inf.values())),
    }
    for nm, (a99, a116) in universes.items():
        rec("universe_size", f"{nm} / E99 side", len(a99), "genes",
            f"{len(i99)} E99 mouse genes")
        rec("universe_size", f"{nm} / E116 side", len(a116), "genes",
            f"{len(i116)} E116 genes")

    # === load both plates ====================================================
    def id_of(lab):
        if "-" in lab or "_" not in lab:
            return None
        return lab.split("_")[0]

    hdr = open(A.OWN_READ).readline().rstrip("\n").split("\t")[1:]
    keep = {A.norm_col(cn) for cn in hdr} - A.OWN_BLANKS
    own, ownw, _ = A.stream_row_totals(A.OWN_READ, keep_wells=keep)
    hs, ms = A.stream_percell_species(A.PUB_UFI)
    labv, _ = A.classify_fig1d(hs, ms)
    mesc = set(labv.index[labv == "mouse"])
    pub, pubw, _ = A.stream_row_totals(A.PUB_READ, keep_wells=mesc)
    log(f"own plate: {len(ownw)} real cells, {own.sum():,.0f} reads")
    log(f"published plate: {len(pubw)} mESC wells, {pub.sum():,.0f} reads")

    def comp(tot, G, accept, level):
        """% composition by class over single-gene mouse Ensembl rows.
        accept=None means no restriction."""
        acc = defaultdict(float)
        for lab, v in tot.items():
            if A.is_trna_row(lab):
                continue
            gid = id_of(lab)
            if gid is None or not gid.startswith("ENSMUSG") or gid not in G:
                continue
            if accept is not None and gid not in accept:
                continue
            bt = G[gid]["biotype"]
            cls = bt if level == "biotype" else A.coarse_class(bt)
            if cls is None:
                continue
            acc[cls] += v
        d = sum(acc.values())
        return pd.Series({k: 100.0 * v / d for k, v in acc.items()}).sort_index(), d

    log("")
    log(f"{'universe':<28} {'level':<8} {'raw_TVD':>8} {'restr_TVD':>10} "
        f"{'release_pp':>11} {'release_pct':>12}")
    results = {}
    for level in ("coarse", "biotype"):
        pa, da = comp(pub, g99, None, level)
        oa, doa = comp(own, g116, None, level)
        raw = A.tvd(pa, oa)
        for nm, (a99, a116) in universes.items():
            pb, db = comp(pub, g99, a99, level)
            ob, dob = comp(own, g116, a116, level)
            m = A.tvd(pb, ob)
            results[(nm, level)] = (raw, m, pa, oa, pb, ob, da, doa, db, dob)
            log(f"{nm:<28} {level:<8} {raw:>8.4f} {m:>10.4f} "
                f"{raw-m:>11.4f} {100*(raw-m)/raw:>11.2f}%")
            rec(f"crossplate_TVD_raw_{level}", nm, round(raw, 4), "pp",
                "ReadCounts, single-gene mouse Ensembl rows",
                "unrestricted composition gap; identical across universes by construction")
            rec(f"crossplate_TVD_restricted_{level}", nm, round(m, 4), "pp",
                "ReadCounts, single-gene mouse rows in the accepted universe")
            rec(f"release_attributable_pp_{level}", nm, round(raw - m, 4), "pp",
                "ReadCounts, single-gene mouse rows")
            rec(f"release_attributable_pct_{level}", nm,
                round(100.0 * (raw - m) / raw, 3), "%", "raw cross-plate TVD",
                "id_only = UPPER bound (discards re-issued ids); "
                "id_plus_symbol_any = LOWER bound (may merge distinct loci)")
            rec(f"restriction_read_loss_own_{level}", nm,
                round(100.0 * (1 - dob / doa), 4), "%",
                "own-plate single-gene Ensembl reads")
            rec(f"restriction_read_loss_published_{level}", nm,
                round(100.0 * (1 - db / da), 4), "%",
                "published-plate single-gene mouse Ensembl reads")
            for cls in sorted(set(pa.index) | set(oa.index) | set(pb.index) | set(ob.index)):
                comp_rows.append({
                    "universe": nm, "level": level, "class": cls,
                    "published_unrestricted": round(float(pa.get(cls, 0.0)), 6),
                    "own_unrestricted": round(float(oa.get(cls, 0.0)), 6),
                    "published_restricted": round(float(pb.get(cls, 0.0)), 6),
                    "own_restricted": round(float(ob.get(cls, 0.0)), 6),
                    "gap_unrestricted_pp": round(float(oa.get(cls, 0.0) - pa.get(cls, 0.0)), 6),
                    "gap_restricted_pp": round(float(ob.get(cls, 0.0) - pb.get(cls, 0.0)), 6)})

    # === how much is Rn7sk alone? ===========================================
    # A single locus carrying 25.6% of the own plate's restriction loss deserves
    # its own number: if the release effect is really one gene, downstream
    # tracks should know that rather than treating it as diffuse.
    log("")
    log("=== leave-one-gene-out: Rn7sk ===")
    rn7sk_own = [l for l in own.index if "_Rn7sk_" in l and "-" not in l]
    rn7sk_pub = [l for l in pub.index if "_Rn7sk_" in l and "-" not in l]
    log(f"  own entries : {rn7sk_own} reads={[f'{own[l]:,.0f}' for l in rn7sk_own]}")
    log(f"  pub entries : {rn7sk_pub} reads={[f'{pub[l]:,.0f}' for l in rn7sk_pub]}")
    for l in rn7sk_own:
        rec("Rn7sk_reads", f"own plate / {l}", int(own[l]), "reads",
            f"{own.sum():,.0f} total own-plate reads",
            "E116 id; absent from E99 under this id, present as ENSMUSG00000065037")
    for l in rn7sk_pub:
        rec("Rn7sk_reads", f"published plate / {l}", int(pub[l]), "reads",
            f"{pub.sum():,.0f} total published mESC reads", "E99 id")

    own_x = own.drop(index=rn7sk_own)
    pub_x = pub.drop(index=rn7sk_pub)
    for level in ("coarse",):
        pa, _ = comp(pub_x, g99, None, level)
        oa, _ = comp(own_x, g116, None, level)
        raw_x = A.tvd(pa, oa)
        pb, _ = comp(pub_x, g99, shared, level)
        ob, _ = comp(own_x, g116, shared, level)
        m_x = A.tvd(pb, ob)
        raw0 = results[("id_only", level)][0]
        m0 = results[("id_only", level)][1]
        log(f"  {level}: with Rn7sk raw={raw0:.4f} restr={m0:.4f} rel={raw0-m0:.4f}")
        log(f"  {level}: without    raw={raw_x:.4f} restr={m_x:.4f} rel={raw_x-m_x:.4f}")
        rec(f"crossplate_TVD_raw_{level}_no_Rn7sk", "id_only", round(raw_x, 4), "pp",
            "ReadCounts excluding both Rn7sk entries")
        rec(f"release_attributable_pp_{level}_no_Rn7sk", "id_only",
            round(raw_x - m_x, 4), "pp", "ReadCounts excluding both Rn7sk entries",
            "compare with the with-Rn7sk figure to see how much one locus drives it")

    pd.DataFrame(comp_rows).to_csv(f"{OUT}/annotation_release_bracket_composition.tsv",
                                   sep="\t", index=False)
    pd.DataFrame(out).to_csv(f"{OUT}/annotation_release_bracket.tsv", sep="\t", index=False)
    log("")
    log(f"wrote {OUT}/annotation_release_bracket.tsv")
    log(f"wrote {OUT}/annotation_release_bracket_composition.tsv")


if __name__ == "__main__":
    main()
