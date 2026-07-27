#!/usr/bin/env python3
"""
01_compare.py -- cross-check our re-run of the published VASA-seq pipeline
against (a) the deposited count table for the same library and (b) the numbers
actually stated in Salmen & De Jonghe et al. 2022.

Writes, into res/vasaplate/:
    comparison_summary.tsv     one row per named quantity, ours vs reference
    per_cell.tsv               per-well species calls, UFIs, genes
    gene_concordance.tsv       per-gene totals, ours vs published
    biotype_composition.tsv    UFI share by biotype, ours vs published

Run once per pipeline run:  ./01_compare.py rrnav2      (baseline, 0 tRNA rows)
                            ./01_compare.py bedv2       (with tRNA)

The published table is the primary reference. The manuscript is secondary and
mostly silent about this library -- see vp_common.NOT_IN_PAPER.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vp_common as C


def log(*a):
    print(*a, flush=True)


def summarise_run(name, df, tag):
    """Barnyard + per-cell stats for one count table."""
    sp = C.species_vector(df.index)
    h, m = C.per_cell_species(df, sp)
    gh = (df[sp.values == "human"] > 0).sum()
    gm = (df[sp.values == "mouse"] > 0).sum()

    lab_f, frac_u = C.classify_fig1d(h, m)
    lab_m, frac_g = C.classify_methods(gh, gm, h, m)

    per_cell = pd.DataFrame({
        "source": tag,
        "ufi_human": h, "ufi_mouse": m, "ufi_total": h + m,
        "genes_human": gh, "genes_mouse": gm,
        "frac_ufi_human": frac_u, "frac_genes_human": frac_g,
        "call_fig1d": lab_f, "call_methods": lab_m,
    })

    kept = per_cell[per_cell.call_fig1d != "discarded"]
    hum = kept[kept.call_fig1d == "human"]
    mou = kept[kept.call_fig1d == "mouse"]

    stats = {
        f"{name}: barcodes >= {C.MIN_UFI} UFIs": len(kept),
        f"{name}: doublet rate, Fig.1d rule (%)": round(C.doublet_rate(lab_f), 2),
        f"{name}: doublet rate, Methods rule (%)": round(C.doublet_rate(lab_m), 2),
        f"{name}: n human (Fig.1d)": int((lab_f == "human").sum()),
        f"{name}: n mouse (Fig.1d)": int((lab_f == "mouse").sum()),
        f"{name}: n mixed (Fig.1d)": int((lab_f == "mixed").sum()),
        f"{name}: HEK293T median UFIs": int(hum.ufi_total.median()) if len(hum) else np.nan,
        f"{name}: HEK293T median genes": int(hum.genes_human.median()) if len(hum) else np.nan,
        f"{name}: HEK293T mean purity": round(float(hum.frac_ufi_human.mean()), 4) if len(hum) else np.nan,
        f"{name}: mESC median UFIs": int(mou.ufi_total.median()) if len(mou) else np.nan,
        f"{name}: mESC median genes": int(mou.genes_mouse.median()) if len(mou) else np.nan,
        f"{name}: mESC mean purity": round(float(1 - mou.frac_ufi_human.mean()), 4) if len(mou) else np.nan,
    }
    return per_cell, stats


def biotype_share(df, tag):
    """UFI share per biotype, over simple (non-combination) rows only."""
    bt = pd.Series([C.biotype_of(i) for i in df.index], index=df.index)
    keep = bt.notna()
    tot = df[keep.values].sum(axis=1)
    grp = tot.groupby(bt[keep].values).sum()
    out = (100.0 * grp / grp.sum()).sort_values(ascending=False)
    return out.rename(tag)


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else "rrnav2"
    if not C.have(run):
        sys.exit(f"count table not found for run '{run}': {C.run_table(run)}\n"
                 f"has stage 7 finished? check squeue first.")

    os.makedirs(C.RES, exist_ok=True)
    rows = {}

    log(f"loading published : {C.PUBLISHED}")
    pub = C.normalise_columns(C.load_counts(C.PUBLISHED))
    log(f"  {pub.shape[0]} rows x {pub.shape[1]} cells")

    log(f"loading ours ({run}) : {C.run_table(run)}")
    ours = C.normalise_columns(C.load_counts(C.run_table(run)))
    log(f"  {ours.shape[0]} rows x {ours.shape[1]} cells")

    common_cells = sorted(set(pub.columns) & set(ours.columns))
    log(f"cells in common: {len(common_cells)}")
    assert len(common_cells) == 384, "expected all 384 wells to line up"
    pub = pub[common_cells]
    ours = ours[common_cells]

    # --- row-set overlap ------------------------------------------------------
    simple_p = {i for i in pub.index if "-" not in i}
    simple_o = {i for i in ours.index if "-" not in i}
    trna_p = {i for i in simple_p if "tRNA" in i}
    trna_o = {i for i in simple_o if "tRNA" in i}
    rows.update({
        "rows: published (all)": len(pub.index),
        "rows: ours (all)": len(ours.index),
        "rows: published simple": len(simple_p),
        "rows: ours simple": len(simple_o),
        "rows: simple shared": len(simple_p & simple_o),
        "rows: simple published-only": len(simple_p - simple_o),
        "tRNA rows: published": len(trna_p),
        "tRNA rows: ours": len(trna_o),
        "tRNA rows: shared": len(trna_p & trna_o),
    })

    # --- per-gene concordance on shared simple rows ---------------------------
    shared = sorted(simple_p & simple_o)
    gp = pub.loc[shared].sum(axis=1)
    go = ours.loc[shared].sum(axis=1)
    conc = pd.DataFrame({"published": gp, "ours": go})
    conc["log2_ratio"] = np.log2((conc.ours + 1) / (conc.published + 1))
    conc.to_csv(f"{C.RES}/gene_concordance.tsv", sep="\t")

    both = conc[(conc.published > 0) & (conc.ours > 0)]
    rows.update({
        "concordance: genes compared": len(conc),
        "concordance: exactly equal totals": int((conc.published == conc.ours).sum()),
        "concordance: Spearman r (per-gene totals)": round(
            float(conc.published.corr(conc.ours, method="spearman")), 4),
        "concordance: Pearson r (log10 totals)": round(float(
            np.log10(both.published).corr(np.log10(both.ours))), 4),
        "concordance: median log2(ours/published)": round(float(both.log2_ratio.median()), 4),
    })

    # per-cell correlation on shared rows
    pc = [float(np.corrcoef(np.log1p(pub.loc[shared, c]),
                            np.log1p(ours.loc[shared, c]))[0, 1])
          for c in common_cells]
    rows["concordance: median per-cell Pearson r (log1p)"] = round(float(np.median(pc)), 4)

    # --- barnyard, both tables, both rules ------------------------------------
    pcp, sp_stats = summarise_run("published", pub, "published")
    pco, so_stats = summarise_run("ours", ours, f"ours_{run}")
    rows.update(sp_stats)
    rows.update(so_stats)

    per_cell = pd.concat([pcp, pco])
    per_cell.index.name = "well"
    per_cell.to_csv(f"{C.RES}/per_cell.tsv", sep="\t")

    # --- biotype composition --------------------------------------------------
    bshare = pd.concat([biotype_share(pub, "published"),
                        biotype_share(ours, f"ours_{run}")], axis=1).fillna(0.0)
    bshare.index.name = "biotype"
    bshare.to_csv(f"{C.RES}/biotype_composition.tsv", sep="\t")

    snc = [b for b in bshare.index if b in C.SNCRNA_BIOTYPES]
    rows["sncRNA UFI share, published (%)"] = round(float(bshare.loc[snc, "published"].sum()), 3)
    rows[f"sncRNA UFI share, ours (%)"] = round(float(bshare.loc[snc, f"ours_{run}"].sum()), 3)
    rows["  paper states, VASA-plate (%)"] = C.PAPER["sncrna_frac"][0]

    if trna_o:
        rows["tRNA UFI share, ours (%)"] = round(float(
            bshare.loc["tRNA", f"ours_{run}"]) if "tRNA" in bshare.index else 0.0, 4)
    rows["tRNA UFI share, published (%)"] = round(float(
        bshare.loc["tRNA", "published"]) if "tRNA" in bshare.index else 0.0, 4)

    # --- write ----------------------------------------------------------------
    out = pd.DataFrame({"value": pd.Series(rows)})
    out.index.name = "quantity"
    dest = f"{C.RES}/comparison_summary.tsv"
    if os.path.exists(dest) and run != "rrnav2":
        prev = pd.read_csv(dest, sep="\t", index_col=0)
        out = prev.combine_first(out).reindex(
            list(dict.fromkeys(list(prev.index) + list(out.index))))
    out.to_csv(dest, sep="\t")

    log("")
    log(f"wrote {dest}")
    for k, v in rows.items():
        log(f"  {k:<52} {v}")


if __name__ == "__main__":
    main()
