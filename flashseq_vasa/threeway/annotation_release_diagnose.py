#!/usr/bin/env python3
"""
annotation_release_diagnose.py -- WHY does the gene universe differ between
Ensembl 99 and 116, and why is the restriction test so much more costly for the
own plate (7.0% of reads) than for the published plate (0.94%)?

This exists because two patterns in annotation_release_effect.py's output cannot
be interpreted without knowing their mechanism, and interpreting them wrongly
would mislead every downstream track:

  1. Several sncRNA classes have IDENTICAL gene counts in both releases but only
     partial id overlap (MiscRna 562 vs 562, only 32 shared; rRNA 354 vs 354,
     only 13 shared). Identical counts with disjoint ids is the signature of an
     id REASSIGNMENT of the same loci, not of genes being added or removed. If
     that is what it is, then "Ensembl 99 did not have this gene" is the wrong
     description -- the correct one is "Ensembl 99 had this locus under a
     different id, and an id-keyed comparison cannot see that".
  2. The asymmetry itself. Restricting to the shared set costs the own plate 7x
     more than the published plate. Two candidate causes: (a) the own plate has
     more of the signal that sits on churned genes, or (b) the churn is
     concentrated in classes the own plate happens to detect. These are
     distinguishable and the answer changes the wording of the conclusion.

Read-only. Writes annotation_release_mechanism.tsv.
"""
import os
import sys
from collections import Counter, defaultdict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import annotation_release_effect as A


def log(*a):
    print(*a, flush=True)


def id_series(gid):
    """Ensembl mouse ids come in blocks. 'ENSMUSG00000xxxxxx' is the original
    series; 'ENSMUSG00002xxxxxx' is a later block used for ncRNA models imported
    from external databases. The block prefix is therefore informative about
    WHERE a gene model came from."""
    return gid[:11]


def main():
    out = []

    def rec(metric, scope, value, unit, denominator, note=""):
        out.append({"metric": metric, "scope": scope, "value": value,
                    "unit": unit, "denominator": denominator, "note": note})

    b99 = A.parse_bed(A.BED99, mouse_only=True)
    b116 = A.parse_bed(A.BED116, mouse_only=False)
    g99, g116 = b99["genes"], b116["genes"]
    i99, i116 = set(g99), set(g116)
    shared, o99, o116 = i99 & i116, i99 - i116, i116 - i99

    log("=== id-series composition of each gene set ===")
    for nm, S, G in (("shared", shared, g99), ("only99", o99, g99), ("only116", o116, g116)):
        cc = Counter(id_series(g) for g in S).most_common(5)
        log(f"  {nm:8} n={len(S):>6}  {cc}")
        for pref, n in cc:
            rec("id_series_count", f"{nm} / {pref}", n, "genes", f"{len(S)} genes in {nm}")

    SNC = A.SNCRNA_BIOTYPES
    log("")
    log("=== the sncRNA classes: identical totals, partial overlap ===")
    for bt in sorted(SNC):
        a = {g for g in i99 if g99[g]["biotype"] == bt}
        b = {g for g in i116 if g116[g]["biotype"] == bt}
        if not a and not b:
            continue
        sh = a & b
        # symbol bridge within this biotype
        sa = defaultdict(set)
        sb = defaultdict(set)
        for g in a - i116:
            sa[g99[g]["symbol"]].add(g)
        for g in b - i99:
            sb[g116[g]["symbol"]].add(g)
        br = set(sa) & set(sb)
        bridged = sum(len(sa[s]) for s in br)
        log(f"  {bt:<10} E99={len(a):>5} E116={len(b):>5} shared_id={len(sh):>5} "
            f"only99={len(a - i116):>5} only116={len(b - i99):>5} "
            f"only99_symbol_bridged={bridged:>5}")
        rec("sncRNA_class_id_overlap", bt,
            f"E99={len(a)};E116={len(b)};shared_id={len(sh)};only99={len(a-i116)};"
            f"only116={len(b-i99)};only99_symbol_bridged={bridged}",
            "counts", "genes of this biotype in each release",
            "identical E99/E116 totals with partial id overlap = id reassignment, "
            "not gene gain/loss")

    # pooled sncRNA symbol bridge
    sa = defaultdict(set)
    sb = defaultdict(set)
    for g in o99:
        if g99[g]["biotype"] in SNC:
            sa[g99[g]["symbol"]].add(g)
    for g in o116:
        if g116[g]["biotype"] in SNC:
            sb[g116[g]["symbol"]].add(g)
    br = set(sa) & set(sb)
    n_a = sum(len(v) for v in sa.values())
    bridged = sum(len(sa[s]) for s in br)
    log("")
    log(f"=== pooled sncRNA only-sets: {len(br)} symbols in common, "
        f"{bridged}/{n_a} only99 genes bridged ===")
    log(f"  bridged examples : {sorted(br)[:10]}")
    log(f"  only99 unbridged : {sorted(set(sa) - set(sb))[:10]}")
    log(f"  only116 unbridged: {sorted(set(sb) - set(sa))[:10]}")
    rec("sncRNA_only99_genes_symbol_bridged_to_only116", "pooled sncRNA classes",
        bridged, "genes", f"{n_a} only-E99 sncRNA genes",
        "same symbol present in E116 under a different id")
    rec("sncRNA_only99_genes_symbol_bridged_pct", "pooled sncRNA classes",
        round(100.0 * bridged / n_a, 2) if n_a else 0.0, "%",
        f"{n_a} only-E99 sncRNA genes")

    log("")
    log("=== example only-set members, by class ===")
    for bt in ("MiscRna", "rRNA", "snoRNA", "miRNA", "snRNA", "lncRNA"):
        a = sorted((g, g99[g]["symbol"]) for g in o99 if g99[g]["biotype"] == bt)[:5]
        b = sorted((g, g116[g]["symbol"]) for g in o116 if g116[g]["biotype"] == bt)[:5]
        log(f"  {bt:<9} only99 : {a}")
        log(f"  {bt:<9} only116: {b}")

    # lncRNA: is the 23,236-gene expansion new ids or new series?
    lnc = {g for g in o116 if g116[g]["biotype"] == "lncRNA"}
    cc = Counter(id_series(g) for g in lnc).most_common(4)
    log("")
    log(f"=== only116 lncRNA n={len(lnc)}: id series {cc} ===")
    for pref, n in cc:
        rec("only116_lncRNA_id_series", pref, n, "genes", f"{len(lnc)} only-E116 lncRNA genes")

    # === the asymmetry: where does each plate's signal sit? ==================
    log("")
    log("=== why the restriction costs the two plates differently ===")
    bt99 = {g: g99[g]["biotype"] for g in g99}
    bt116 = {g: g116[g]["biotype"] for g in g116}

    def id_of(lab):
        if "-" in lab or "_" not in lab:
            return None
        return lab.split("_")[0]

    # own plate, ReadCounts, 12 real cells
    hdr = open(A.OWN_READ).readline().rstrip("\n").split("\t")[1:]
    keep = {A.norm_col(cn) for cn in hdr} - A.OWN_BLANKS
    own, wells, _ = A.stream_row_totals(A.OWN_READ, keep_wells=keep)
    # published plate, ReadCounts, mESC wells
    hs, ms = A.stream_percell_species(A.PUB_UFI)
    lab_, _ = A.classify_fig1d(hs, ms)
    mesc = set(lab_.index[lab_ == "mouse"])
    pub, pwells, _ = A.stream_row_totals(A.PUB_READ, keep_wells=mesc)

    rowsx = []
    for plate, tot, btmap, univ in (("own_E116", own, bt116, i116),
                                    ("published_E99", pub, bt99, i99)):
        lost = defaultdict(float)
        kept = defaultdict(float)
        for lab2, v in tot.items():
            if A.is_trna_row(lab2):
                continue
            gid = id_of(lab2)
            if gid is None or not gid.startswith("ENSMUSG") or gid not in btmap:
                continue
            cls = A.coarse_class(btmap[gid])
            if gid in shared:
                kept[cls] += v
            else:
                lost[cls] += v
        tt = sum(kept.values()) + sum(lost.values())
        for cls in sorted(set(kept) | set(lost)):
            rowsx.append({"plate": plate, "coarse_class": cls,
                          "reads_on_shared_genes": kept.get(cls, 0.0),
                          "reads_on_release_only_genes": lost.get(cls, 0.0),
                          "pct_of_plate_lost_by_restriction":
                              round(100.0 * lost.get(cls, 0.0) / tt, 4)})
        tot_lost = 100.0 * sum(lost.values()) / tt
        log(f"  {plate}: restriction removes {tot_lost:.3f}% of single-gene reads")
        for cls in sorted(lost, key=lambda k: -lost[k])[:4]:
            log(f"      {cls:<14} {100.0*lost[cls]/tt:7.3f} pp")
        rec("restriction_read_loss", plate, round(tot_lost, 4), "%",
            "single-gene Ensembl reads on that plate",
            "reads on genes absent from the other release")
    dfx = pd.DataFrame(rowsx)

    OUT = sys.argv[1] if len(sys.argv) > 1 else f"{A.W}/res/threeway"
    dfx.to_csv(f"{OUT}/annotation_release_asymmetry.tsv", sep="\t", index=False)
    pd.DataFrame(out).to_csv(f"{OUT}/annotation_release_mechanism.tsv", sep="\t", index=False)
    log("")
    log(f"wrote {OUT}/annotation_release_mechanism.tsv")
    log(f"wrote {OUT}/annotation_release_asymmetry.tsv")


if __name__ == "__main__":
    main()
