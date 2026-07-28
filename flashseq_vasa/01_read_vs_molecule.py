#!/usr/bin/env python3
"""read_vs_molecule.py -- quantify the read-vs-molecule asymmetry in the
VASA/FLASH-seq comparison, from tables already on disk.

VASA counts deduplicated molecules (UFIs; TranscriptCounts = UFIs with the
binomial collision correction applied). smartseq_noUMI on FLASH-seq will count
READS, because there is no UMI to deduplicate. So "% protein-coding" is not the
same quantity on the two sides, and any ratio between them inherits that.

This measures VASA's own read:UFI ratio -- the amplification/duplication factor
VASA removes and FLASH-seq cannot -- per cell and per biotype, so the asymmetry
is a stated number.

Denominator note: all three tables are the *total* (spliced+unspliced) analysis
matrices over the same 222,412 entries x 12 cells, so per-cell column sums are
directly comparable. rRNA here is the rRNA *biotype* of reads that SURVIVED
step 3 (i.e. the residual after bwa depletion) -- it is NOT the 21.39 % bwa
figure, which is over all reads before depletion. Two different denominators;
kept separate on purpose.

Reproducibility: pandas/numpy versions printed; no randomness, no seed needed.
"""
import json
import sys

import numpy as np
import pandas as pd

AN = "/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out/analysis"
PFX = f"{AN}/ZHA9292A1_analysis_total"

print(f"pandas={pd.__version__} numpy={np.__version__}", file=sys.stderr)

reads = pd.read_csv(f"{PFX}.ReadCounts.tsv", sep="\t", index_col=0)
ufis = pd.read_csv(f"{PFX}.UFICounts.tsv", sep="\t", index_col=0)
tx = pd.read_csv(f"{PFX}.TranscriptCounts.tsv", sep="\t", index_col=0)
meta = pd.read_csv(f"{AN}/gene_metadata.tsv", sep="\t", index_col=0)

# The three matrices must be the same object, else nothing below is comparable.
assert list(reads.columns) == list(ufis.columns) == list(tx.columns), "column mismatch"
assert reads.index.equals(ufis.index) and reads.index.equals(tx.index), "index mismatch"
assert reads.index.equals(meta.index), "gene_metadata index does not match the matrices"
CELLS = list(reads.columns)
print(f"shape={reads.shape} cells={CELLS}", file=sys.stderr)

# ---------------------------------------------------------------- per cell ---
per_cell = pd.DataFrame({
    "reads": reads.sum(axis=0),
    "UFIs": ufis.sum(axis=0),
    "transcripts": tx.sum(axis=0),
})
per_cell["read_per_UFI"] = per_cell["reads"] / per_cell["UFIs"]
per_cell["read_per_transcript"] = per_cell["reads"] / per_cell["transcripts"]
per_cell["UFI_per_transcript"] = per_cell["UFIs"] / per_cell["transcripts"]
# genes detected on reads vs on UFIs -- for the detection-count decision
per_cell["genes_det_reads"] = (reads > 0).sum(axis=0)
per_cell["genes_det_UFIs"] = (ufis > 0).sum(axis=0)
per_cell["genes_det_diff"] = per_cell["genes_det_reads"] - per_cell["genes_det_UFIs"]
per_cell.index.name = "cell"

# ------------------------------------------------------------- per biotype ---
bt = meta["biotype"]
g_reads = reads.groupby(bt).sum()
g_ufis = ufis.groupby(bt).sum()
g_tx = tx.groupby(bt).sum()

per_bt = pd.DataFrame({
    "reads": g_reads.sum(axis=1),
    "UFIs": g_ufis.sum(axis=1),
    "transcripts": g_tx.sum(axis=1),
    "n_entries": bt.value_counts(),
})
per_bt["read_per_UFI"] = per_bt["reads"] / per_bt["UFIs"]
# share of the library on each currency -- this is the number the comparison
# will actually quote, and the point is that the two columns differ
per_bt["pct_of_reads"] = 100 * per_bt["reads"] / per_bt["reads"].sum()
per_bt["pct_of_UFIs"] = 100 * per_bt["UFIs"] / per_bt["UFIs"].sum()
per_bt["pct_of_transcripts"] = 100 * per_bt["transcripts"] / per_bt["transcripts"].sum()
per_bt["pct_read_minus_pct_UFI"] = per_bt["pct_of_reads"] - per_bt["pct_of_UFIs"]
# ratio of the two shares: >1 means the biotype looks BIGGER on reads than on
# molecules, i.e. quoting VASA UFI% against FLASH-seq read% understates it
per_bt["share_ratio_read_over_UFI"] = per_bt["pct_of_reads"] / per_bt["pct_of_UFIs"]
per_bt = per_bt.sort_values("reads", ascending=False)
per_bt.index.name = "biotype"

# per-biotype x per-cell ratio, so the spread across cells is visible too
bt_cell = (g_reads / g_ufis).replace([np.inf, -np.inf], np.nan)
per_bt["read_per_UFI_cell_min"] = bt_cell.min(axis=1)
per_bt["read_per_UFI_cell_max"] = bt_cell.max(axis=1)

# ------------------------------------------- is duplication expression-dep? --
# per-entry: does a more-expressed entry carry more reads per molecule?
ent = pd.DataFrame({"reads": reads.sum(axis=1), "UFIs": ufis.sum(axis=1)})
ent = ent[ent["UFIs"] >= 10]                     # ratio is meaningless at 1-2 UFIs
ent["ratio"] = ent["reads"] / ent["UFIs"]
qs = pd.qcut(ent["UFIs"], 10, labels=False, duplicates="drop")
decile = ent.groupby(qs).agg(
    n_entries=("ratio", "size"),
    UFI_min=("UFIs", "min"),
    UFI_max=("UFIs", "max"),
    reads=("reads", "sum"),
    UFIs=("UFIs", "sum"),
)
decile["read_per_UFI"] = decile["reads"] / decile["UFIs"]
decile.index.name = "UFI_decile"
# Spearman without scipy (envs/vasa has no scipy, and installing it into a
# pandas-2.x-pinned env for one correlation is not worth the risk): Spearman's
# rho IS the Pearson correlation of the average ranks, ties included.
rho = ent["UFIs"].rank().corr(ent["ratio"].rank(), method="pearson")

# --------------------------------------------------------------- write out ---
long = []
for cell, r in per_cell.iterrows():
    long.append({"level": "cell", "key": cell, "biotype": "ALL",
                 "reads": r["reads"], "UFIs": r["UFIs"], "transcripts": r["transcripts"],
                 "read_per_UFI": r["read_per_UFI"],
                 "read_per_transcript": r["read_per_transcript"],
                 "pct_of_reads": np.nan, "pct_of_UFIs": np.nan,
                 "pct_of_transcripts": np.nan, "share_ratio_read_over_UFI": np.nan,
                 "genes_det_reads": r["genes_det_reads"],
                 "genes_det_UFIs": r["genes_det_UFIs"]})
for b, r in per_bt.iterrows():
    long.append({"level": "biotype", "key": b, "biotype": b,
                 "reads": r["reads"], "UFIs": r["UFIs"], "transcripts": r["transcripts"],
                 "read_per_UFI": r["read_per_UFI"],
                 "read_per_transcript": r["reads"] / r["transcripts"] if r["transcripts"] else np.nan,
                 "pct_of_reads": r["pct_of_reads"], "pct_of_UFIs": r["pct_of_UFIs"],
                 "pct_of_transcripts": r["pct_of_transcripts"],
                 "share_ratio_read_over_UFI": r["share_ratio_read_over_UFI"],
                 "genes_det_reads": np.nan, "genes_det_UFIs": np.nan})
out = pd.DataFrame(long)
out.to_csv("read_vs_molecule.tsv", sep="\t", index=False, float_format="%.6g")
per_bt.to_csv("read_vs_molecule_biotype_full.tsv", sep="\t", float_format="%.6g")
decile.to_csv("read_vs_molecule_deciles.tsv", sep="\t", float_format="%.6g")

# --- the rRNA-biotype residual, and whether one locus still dominates it -----
# Commit a9d0c46 documented ~600x the published share of rRNA-biotype UFIs, all
# on ENSMUSG00000106106 (an Ensembl entry stored antisense to the real
# transcript). "No fix applied yet" at that time. If that locus still dominates,
# the rRNA *biotype* fraction in these tables is not a usable comparator and the
# rRNA comparison must come from the bwa stage instead.
rrna_rows = meta.index[meta["biotype"].isin(["rRNA", "MtRrna", "rRNAPseudogene"])]
rrna_tab = pd.DataFrame({
    "reads": reads.loc[rrna_rows].sum(axis=1),
    "UFIs": ufis.loc[rrna_rows].sum(axis=1),
    "biotype": meta.loc[rrna_rows, "biotype"],
}).sort_values("UFIs", ascending=False)
rrna_tab.to_csv("rrna_biotype_entries.tsv", sep="\t")
susp = [i for i in rrna_rows if "ENSMUSG00000106106" in i]
rrna_ufi_total = float(rrna_tab["UFIs"].sum())
susp_ufi = float(ufis.loc[susp].sum().sum()) if susp else 0.0

# --- detection floor: how many entries sit at exactly 1 read / 1 UFI? --------
# A gene seen on 1 read is "detected" on ReadCounts; VASA's dedup can only ever
# move it to 1 UFI, so detection counts differ far less than totals do. This is
# the evidence for the (iii) column choice.
one_read = int((reads == 1).sum().sum())
one_ufi = int((ufis == 1).sum().sum())

summary = {
    "n_cells": len(CELLS),
    "rrna_biotype_pct_of_reads": float(100 * rrna_tab["reads"].sum() / reads.values.sum()),
    "rrna_biotype_pct_of_UFIs": float(100 * rrna_ufi_total / ufis.values.sum()),
    "rrna_top_entry": str(rrna_tab.index[0]),
    "rrna_top_entry_pct_of_rrna_UFIs": float(100 * rrna_tab["UFIs"].iloc[0] / rrna_ufi_total),
    "ENSMUSG00000106106_UFIs": susp_ufi,
    "ENSMUSG00000106106_pct_of_rrna_UFIs": float(100 * susp_ufi / rrna_ufi_total) if rrna_ufi_total else 0.0,
    "n_cellgene_at_1_read": one_read,
    "n_cellgene_at_1_UFI": one_ufi,
    "shape": list(reads.shape),
    "read_per_UFI_min": float(per_cell["read_per_UFI"].min()),
    "read_per_UFI_max": float(per_cell["read_per_UFI"].max()),
    "read_per_UFI_median": float(per_cell["read_per_UFI"].median()),
    "read_per_UFI_cell_argmin": str(per_cell["read_per_UFI"].idxmin()),
    "read_per_UFI_cell_argmax": str(per_cell["read_per_UFI"].idxmax()),
    "read_per_transcript_median": float(per_cell["read_per_transcript"].median()),
    "UFI_per_transcript_median": float(per_cell["UFI_per_transcript"].median()),
    "library_read_per_UFI": float(per_cell["reads"].sum() / per_cell["UFIs"].sum()),
    "spearman_UFI_vs_ratio": float(rho),
    "n_entries_ratio_test": int(len(ent)),
    "decile_ratio_lowest": float(decile["read_per_UFI"].iloc[0]),
    "decile_ratio_highest": float(decile["read_per_UFI"].iloc[-1]),
    "genes_det_diff_median": float(per_cell["genes_det_diff"].median()),
    "genes_det_diff_max": int(per_cell["genes_det_diff"].max()),
    "pandas": pd.__version__,
    "numpy": np.__version__,
}
with open("read_vs_molecule_summary.json", "w") as fh:
    json.dump(summary, fh, indent=1)

pd.set_option("display.width", 220, "display.max_columns", 40)
print("=== per cell ===")
print(per_cell.round(4).to_string())
print("\n=== per biotype (top 14 by reads) ===")
print(per_bt.head(14)[["reads", "UFIs", "read_per_UFI", "pct_of_reads", "pct_of_UFIs",
                       "pct_of_transcripts", "pct_read_minus_pct_UFI",
                       "share_ratio_read_over_UFI",
                       "read_per_UFI_cell_min", "read_per_UFI_cell_max"]].round(4).to_string())
print("\n=== read:UFI by expression decile ===")
print(decile.round(4).to_string())
print(f"\nspearman(UFIs, read:UFI) over {len(ent)} entries with >=10 UFIs = {rho:.4f}")
print("\n=== summary ===")
print(json.dumps(summary, indent=1))
print("DONE_PART_B")
