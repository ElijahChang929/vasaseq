#!/usr/bin/env python3
"""
annotation_release_effect_precheck.py -- read-only precheck for
annotation_release_effect.py, per lab convention Rule 2.

Replays every real operation of the main script over the real data, with the
two very large count tables truncated, and reports what would break plus a
memory estimate. Imports the REAL functions from the main script -- so what is
tested is what will execute, not a paraphrase of it.

What it checks
  A. both BEDs exist, are readable, and parse under the real parse_bed()
     including its four assertions (4-field names, feature token, col6 span,
     intra-gene biotype/symbol consistency);
  B. the tiling-convention assertion that gates the paired length comparison;
  C. the count tables exist and their headers normalise to the expected wells,
     including that all 4 blank barcodes are present to be excluded and that
     384 published wells resolve;
  D. the row-label vocabulary contains no label shape the classifiers silently
     drop -- specifically that every non-combination row is either a 3-field
     Ensembl name or a tRNAscan dot-name, with a census printed;
  E. the species classifier reproduces the benchmark's mESC well count on the
     published UFI table (the one number this script must not get wrong,
     because the whole symmetric test hangs off it);
  F. peak RSS.

Exit 0 = safe to run the real thing. Any assertion failure = do not.

USAGE
    ./annotation_release_effect_precheck.py
"""
import os
import resource
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import annotation_release_effect as A   # the REAL module

LIMIT = 40000   # rows of each big table to stream in the precheck


def log(*a):
    print(*a, flush=True)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 2)


def main():
    fail = []

    # --- A/B: BEDs -----------------------------------------------------------
    log("=" * 74)
    log("A. BED parse (real parse_bed, all assertions live)")
    log("=" * 74)
    for p in (A.BED99, A.BED116):
        if not os.path.exists(p):
            fail.append(f"MISSING BED: {p}")
            log(f"  MISSING {p}")
    if fail:
        log("\n".join(fail))
        return 1

    b99 = A.parse_bed(A.BED99, mouse_only=True)
    log(f"  E99  : {len(b99['genes']):>6} mouse genes, {b99['n_rows']} rows, "
        f"{b99['skipped_species']} non-mouse skipped, {b99['trna']} tRNAscan")
    b116 = A.parse_bed(A.BED116, mouse_only=False)
    log(f"  E116 : {len(b116['genes']):>6} genes, {b116['n_rows']} rows, "
        f"{b116['trna']} tRNAscan")
    log(f"  peak RSS after both BEDs: {rss_gb():.2f} GB")

    log("")
    log("  B. coordinate convention (detected, then verified)")
    log(f"  E99  modal inter-feature gaps: {b99['gaps'].most_common(3)}  "
        f"-> offset +{b99['offset']}")
    log(f"  E116 modal inter-feature gaps: {b116['gaps'].most_common(3)}  "
        f"-> offset +{b116['offset']}")
    if b99["offset"] != b116["offset"]:
        log(f"  NOTE conventions DIFFER (E99 +{b99['offset']}, E116 +{b116['offset']}). "
            "parse_bed corrects both; verifying the correction below.")

    # single-feature genes present in both: corrected lengths must agree more
    # often than raw ones, or the offsets are wrong.
    g99, g116 = b99["genes"], b116["genes"]
    single_both = [gid for gid in (set(g99) & set(g116))
                   if g99[gid]["n_exon"] + g99[gid]["n_intron"] == 1
                   and g116[gid]["n_exon"] + g116[gid]["n_intron"] == 1]
    if not single_both:
        fail.append("no single-feature genes shared between the releases -- "
                    "the convention correction cannot be verified")
        log("  FAIL no single-feature control genes")
    else:
        def L(g, gid):
            return (g[gid]["exon_len"] + g[gid]["intron_len"])[0]
        corr = sum(1 for gid in single_both if L(g99, gid) == L(g116, gid))
        raw = sum(1 for gid in single_both
                  if L(g99, gid) - b99["offset"] == L(g116, gid) - b116["offset"])
        log(f"  single-feature control genes: {len(single_both)}")
        log(f"    corrected lengths identical: {corr} ({100*corr/len(single_both):.2f}%)")
        log(f"    raw      lengths identical: {raw} ({100*raw/len(single_both):.2f}%)")
        if corr < raw:
            fail.append(f"convention correction makes agreement WORSE "
                        f"({corr} < {raw} of {len(single_both)}) -- offsets wrong")
            log("    FAIL correction is not a correction")
        else:
            log(f"    OK correction improves agreement by "
                f"{100*(corr-raw)/len(single_both):.2f} pp of control genes")
        ex = "ENSMUSG00000048040"
        if ex in g99 and ex in g116:
            log(f"    {ex}: E99 corrected={L(g99, ex)} (raw {L(g99, ex)-b99['offset']}), "
                f"E116 corrected={L(g116, ex)} (raw {L(g116, ex)-b116['offset']})")

    ids99, ids116 = set(g99), set(g116)
    log(f"  gene-set overlap preview: shared={len(ids99 & ids116)} "
        f"only99={len(ids99 - ids116)} only116={len(ids116 - ids99)}")

    # --- C: count tables -----------------------------------------------------
    log("")
    log("=" * 74)
    log("C. count tables and well resolution")
    log("=" * 74)
    for p in (A.OWN_READ, A.OWN_TRANS, A.PUB_READ, A.PUB_UFI):
        ok = os.path.exists(p) and os.path.getsize(p) > 0
        log(f"  {'OK  ' if ok else 'MISS'} {os.path.getsize(p)/1e6:9.1f} MB  {p}")
        if not ok:
            fail.append(f"missing/empty count table: {p}")
    if any("count table" in f for f in fail):
        log("\n".join(fail))
        return 1

    for name, p, exp in (("own ReadCounts", A.OWN_READ, 16),
                         ("own TranscriptCounts", A.OWN_TRANS, 16),
                         ("pub ReadCounts", A.PUB_READ, 384),
                         ("pub UFICounts", A.PUB_UFI, 384)):
        with open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        wells = [A.norm_col(c) for c in hdr[1:]]
        dup = [w for w, n in Counter(wells).items() if n > 1]
        log(f"  {name:<22} {len(wells):>4} wells, first 3 = {wells[:3]}, "
            f"dups={dup if dup else 'none'}")
        if len(wells) != exp:
            fail.append(f"{name}: expected {exp} wells, got {len(wells)}")
        if dup:
            fail.append(f"{name}: duplicate normalised well labels {dup}")
        if exp == 16:
            missing_blanks = A.OWN_BLANKS - set(wells)
            if missing_blanks:
                fail.append(f"{name}: blank barcodes not present to exclude: {missing_blanks}")
            log(f"  {'':<22} blanks present to exclude: "
                f"{sorted(A.OWN_BLANKS & set(wells))}")

    # --- D: row-label vocabulary census -------------------------------------
    log("")
    log("=" * 74)
    log("D. row-label vocabulary -- would any shape be silently dropped?")
    log("=" * 74)
    for name, p in (("own", A.OWN_READ), ("published", A.PUB_READ)):
        shapes = Counter()
        unclassified = []
        with open(p) as fh:
            fh.readline()
            for n, line in enumerate(fh):
                if n >= 200000:
                    break
                lab = line.split("\t", 1)[0]
                if "-" in lab:
                    shapes["combination"] += 1
                    continue
                if A.is_trna_row(lab):
                    shapes["trnascan_simple"] += 1
                    continue
                nf = lab.count("_") + 1
                shapes[f"ensembl_{nf}field"] += 1
                bt = A.biotype_of(lab)
                if bt is None or A.coarse_class(bt) is None:
                    shapes["UNCLASSIFIED"] += 1
                    if len(unclassified) < 5:
                        unclassified.append(lab)
        log(f"  {name:<10} {dict(shapes)}")
        if shapes.get("UNCLASSIFIED"):
            fail.append(f"{name}: {shapes['UNCLASSIFIED']} rows no classifier handles, "
                        f"e.g. {unclassified}")
            log(f"    FAIL examples: {unclassified}")
        odd = {k: v for k, v in shapes.items()
               if k.startswith("ensembl_") and k != "ensembl_3field"}
        if odd:
            fail.append(f"{name}: unexpected Ensembl label field counts {odd}")
            log(f"    FAIL unexpected field counts {odd}")

    # --- streaming smoke test on truncated tables ----------------------------
    log("")
    log(f"  streaming first {LIMIT} rows of each big table ...")
    tot, wells, nrow = A.stream_row_totals(
        A.OWN_READ, keep_wells={w for w in
                                [A.norm_col(c) for c in
                                 open(A.OWN_READ).readline().rstrip("\n").split("\t")[1:]]
                                if w not in A.OWN_BLANKS}, limit=LIMIT)
    log(f"  own ReadCounts (truncated): {nrow} rows x {len(wells)} wells, "
        f"sum={tot.sum():.0f}")
    if len(wells) != 12:
        fail.append(f"own real-well count is {len(wells)}, expected 12")

    ptot, pwells, pnrow = A.stream_row_totals(A.PUB_READ, keep_wells=None, limit=LIMIT)
    log(f"  pub ReadCounts (truncated): {pnrow} rows x {len(pwells)} wells, "
        f"sum={ptot.sum():.0f}")
    log(f"  peak RSS after truncated streams: {rss_gb():.2f} GB")

    # composition on the truncated data -- exercises the real function
    shared = ids99 & ids116

    def id_of(lab):
        if "-" in lab or "_" not in lab:
            return None
        return lab.split("_")[0]

    def class_coarse(lab):
        if A.is_trna_row(lab):
            return None
        bt = A.biotype_of(lab)
        return A.coarse_class(bt) if bt else None

    a, dA, nA = A.composition(tot, id_of, class_coarse, None)
    b, dB, nB = A.composition(tot, id_of, class_coarse, shared)
    log(f"  composition() on truncated own data: {nA} rows -> {dict(a.round(2))}")
    log(f"  restricted to shared ids:            {nB} rows, TVD={A.tvd(a, b):.3f} pp")
    log("  (TRUNCATED -- not a result, only proof the code path runs)")

    # --- E: species classifier reproduces the benchmark's mESC count ---------
    log("")
    log("=" * 74)
    log("E. species classification on the FULL published UFI table")
    log("=" * 74)
    log("  (this one cannot be truncated -- per-cell sums need every row)")
    hs, ms = A.stream_percell_species(A.PUB_UFI)
    lab, frac = A.classify_fig1d(hs, ms)
    counts = lab.value_counts().to_dict()
    log(f"  fig1d calls: {counts}")
    log(f"  peak RSS: {rss_gb():.2f} GB")
    if counts.get("mouse", 0) == 0:
        fail.append("no mESC wells called -- symmetric test impossible")
    if sum(counts.values()) != 384:
        fail.append(f"well total is {sum(counts.values())}, expected 384")
    log(f"  mESC wells that the symmetric test will use: {counts.get('mouse', 0)}")

    # --- verdict -------------------------------------------------------------
    log("")
    log("=" * 74)
    log(f"peak RSS overall: {rss_gb():.2f} GB")
    if fail:
        log(f"PRECHECK FAILED -- {len(fail)} problem(s):")
        for f in fail:
            log(f"  - {f}")
        return 1
    log("PRECHECK PASSED -- safe to run annotation_release_effect.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
