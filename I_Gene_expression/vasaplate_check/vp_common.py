#!/usr/bin/env python3
"""
vp_common.py -- shared paths, loaders and species logic for the VASA-plate
species-mixing control cross-check.

The library is SRR14783059 / GSM5369495, `vasaplate-HEK293T-mESC`: the one
VASA-plate library in GSE176588. HEK293T is human, mESC is mouse, so it maps to
a concatenated GRCh38+GRCm38 (Ensembl 99) reference and every gene name carries
its species in its Ensembl id -- ENSG* human, ENSMUSG* mouse. There is no
species prefix on the gene names themselves; that is the authors' convention and
we reproduce it.

WHICH TABLE TO COMPARE
----------------------
GEO ships exactly one count table for this sample:

    GSM5369495_vasaplate_HEK293T-mESC_split_total.UFICounts.tsv.gz

so UFICounts is what the comparison has to use, even though UFIs saturate on a
6 nt UMI (4,096 molecules per gene per cell) and TranscriptCounts is the better
quantity for one's own analysis. Both are 384 columns -- all plate wells,
unfiltered, not 384 QC-passing cells.
"""
import gzip
import os

import numpy as np
import pandas as pd

VASASEQ = "/nemo/lab/turnerj/working/guangxin/vasaseq"
FQDIR = f"{VASASEQ}/data/ref/fastq_vasaplate"
TABLES = f"{VASASEQ}/data/ref/processed/tables"
RES = f"{VASASEQ}/res/vasaplate"
FIGS = f"{RES}/figures"

PUBLISHED = f"{TABLES}/GSM5369495_vasaplate_HEK293T-mESC_split_total.UFICounts.tsv.gz"
BC2CELL = f"{TABLES}/GSM5369495_bc2cellID_vasaplate.tsv.gz"

# The two runs being compared. `rrnav2` is the baseline: correct rRNA reference,
# annotation BED with zero tRNA rows. `bedv2` adds the 1,758 tRNA rows and
# changes nothing else -- same BAMs, same everything upstream of step 5.
RUNS = {
    "rrnav2": f"{FQDIR}/vasaplate_out_rrnav2",
    "bedv2": f"{FQDIR}/vasaplate_out_bedv2",
}


def run_table(run, kind="total", counts="UFICounts"):
    return f"{RUNS[run]}_{kind}.{counts}.tsv"


def have(run, kind="total", counts="UFICounts"):
    p = run_table(run, kind, counts)
    return os.path.exists(p) and os.path.getsize(p) > 0


def load_counts(path, index_col=0):
    """Read a count table. int32 halves the memory of a 588k x 384 table."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    ncol = len(header) - 1
    df = pd.read_csv(path, sep="\t", index_col=index_col,
                     dtype={c: np.int32 for c in header[1:]})
    assert df.shape[1] == ncol, f"{path}: {df.shape[1]} != {ncol}"
    return df


def normalise_columns(df):
    """Column labels to bare zero-padded well ids ('001'..'384').

    Ours arrive as 'vasaplate_out_rrnav2/SRR14783059_001'; the published table
    already uses '001'. Making them agree is what lets the two be compared
    cell by cell rather than only in aggregate.
    """
    out = []
    for c in df.columns:
        c = str(c).rsplit("/", 1)[-1]
        c = c.rsplit("_", 1)[-1]
        out.append(c.zfill(3))
    df = df.copy()
    df.columns = out
    return df


# --- species -----------------------------------------------------------------
# Combination rows ('geneA-geneB') are the un-collapsed multi-annotation output
# and are >80% of the published table's rows. A combination naming genes from
# both species cannot be assigned, and is excluded from species arithmetic
# rather than arbitrarily broken -- the paper's Fig. 1d works on UFIs assigned
# to a species, not on every row.

def species_of(idx):
    """'human' | 'mouse' | 'mixed' | 'trna' | 'other' for one row label."""
    parts = idx.split("-")
    h = any(p.startswith("ENSG") for p in parts)
    m = any(p.startswith("ENSMUSG") for p in parts)
    if h and m:
        return "mixed"
    if h:
        return "human"
    if m:
        return "mouse"
    if "tRNA" in idx:
        return "trna"
    return "other"


def species_vector(index):
    return pd.Series([species_of(i) for i in index], index=index)


def per_cell_species(df, sp=None):
    """(human UFIs, mouse UFIs) per cell, from unambiguous rows only."""
    if sp is None:
        sp = species_vector(df.index)
    return df[sp.values == "human"].sum(), df[sp.values == "mouse"].sum()


# --- the paper's two doublet rules -------------------------------------------
# They disagree, and by a lot, so any number quoted has to say which was used.
#
#   Fig. 1d caption (p. 1781): "Barcodes with more than 25% of detected UFIs
#   belonging to the other species were considered doublets/mixed. Detected
#   barcodes with low UFIs (<7,500) were discarded."
#
#   Methods (p. 18): "barcodes with more that 75% of the genes assigned to only
#   one of either mouse or human were considered singlets. Cells with fewer than
#   7,500 UFIs were filtered out."
#
# The first thresholds UFI fraction, the second gene fraction.
MIN_UFI = 7500


def classify_fig1d(h, m, min_ufi=MIN_UFI):
    tot = h + m
    keep = tot >= min_ufi
    frac_h = h.where(keep) / tot.where(keep)
    lab = pd.Series("discarded", index=h.index, dtype=object)
    lab[keep & (frac_h > 0.75)] = "human"
    lab[keep & (frac_h < 0.25)] = "mouse"
    lab[keep & (frac_h >= 0.25) & (frac_h <= 0.75)] = "mixed"
    return lab, frac_h


def classify_methods(gh, gm, h, m, min_ufi=MIN_UFI):
    """Same 7,500-UFI gate, but purity measured on GENES detected."""
    tot_u = h + m
    keep = tot_u >= min_ufi
    gtot = gh + gm
    frac_h = gh.where(keep) / gtot.where(keep)
    lab = pd.Series("discarded", index=gh.index, dtype=object)
    lab[keep & (frac_h > 0.75)] = "human"
    lab[keep & (frac_h < 0.25)] = "mouse"
    lab[keep & (frac_h >= 0.25) & (frac_h <= 0.75)] = "mixed"
    return lab, frac_h


def doublet_rate(lab):
    called = lab[lab != "discarded"]
    if len(called) == 0:
        return float("nan")
    return 100.0 * (called == "mixed").sum() / len(called)


# --- what the paper actually states -------------------------------------------
# Only the entries with a number are checkable. The rest are recorded so the
# report can say plainly that there is nothing to compare against, rather than
# quietly omitting them.
PAPER = {
    "doublet_rate_vasadrop": (3.08, "%", "p.1781 Fig.1d -- VASA-DROP, not plate"),
    "genes_per_cell_75k": (9480, 1252, "p.1781 Fig.1f, VASA-plate HEK293T, 75k trimmed reads/cell"),
    "genes_per_cell_750k": (15248, 1092, "p.1781 Ext.Data Fig.2e, VASA-plate, 750k trimmed reads/cell"),
    "unspliced_frac": (44.1, 10.1, "p.1783, VASA-plate"),
    "sncrna_frac": (1.4, None, "p.1781, VASA-plate"),
}

NOT_IN_PAPER = [
    "doublet rate for the VASA-PLATE mixing control (GSM5369495) -- the paper's "
    "3.08% is VASA-drop (GSM5369496); this library appears only as a source of "
    "cells for the benchmarking panels",
    "ribosomal read fraction for VASA-seq -- Extended Data Fig. 2f is an "
    "un-numbered bar chart and the Discussion is qualitative",
    "tRNA fraction for VASA-seq -- the paper's single tRNA sentence attributes "
    "tRNA detection to Smart-seq-total, and the Methods name no tRNA database",
]

SNCRNA_BIOTYPES = {"MiscRna", "snoRNA", "snRNA", "scaRNA", "miRNA", "ribozyme",
                   "sRNA", "MtTrna", "rRNA", "MtRrna", "vaultRNA"}


def biotype_of(idx):
    """Biotype token of a simple row; None for combinations, 'tRNA' for tRNA."""
    if "-" in idx:
        return None
    if "_" not in idx:
        return "tRNA" if "tRNA" in idx else None
    return idx.rsplit("_", 1)[-1]
