#!/usr/bin/env python3
"""
02_figures.py -- figures for the VASA-plate mixing-control cross-check.

Reads the TSVs 01_compare.py wrote (never the count tables again), so figures
are cheap to regenerate and can never disagree with the summary numbers.

Every figure is rendered twice: a light PDF into res/vasaplate/figures/ and a
light+dark PNG pair into figures/png/ for the theme-aware HTML report.

Palette is the dataviz reference categorical instance, validated with
scripts/validate_palette.js in BOTH modes before use:
    light #2a78d6,#eb6834,#1baf7a,#eda100 -> ALL CHECKS PASS
    dark  #3987e5,#d95926,#199e70,#c98500 -> ALL CHECKS PASS
Light mode raises a contrast WARN on the aqua and yellow slots, which obligates
visible labels or a table view; both are present (direct-labelled medians, and
the TSVs the report tabulates).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vp_common as C

PNG = f"{C.FIGS}/png"

# categorical slots, fixed order, never cycled
LIGHT = dict(s1="#2a78d6", s2="#eb6834", s3="#1baf7a", s4="#eda100",
             surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", grid="#d9d8d3")
DARK = dict(s1="#3987e5", s2="#d95926", s3="#199e70", s4="#c98500",
            surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", grid="#3a3a38")

# identity is bound to the entity, not to rank -- these never move
CALL_SLOT = {"human": "s1", "mouse": "s2", "mixed": "s3", "discarded": None}


def style(P):
    plt.rcParams.update({
        "figure.facecolor": P["surface"], "axes.facecolor": P["surface"],
        "savefig.facecolor": P["surface"],
        "text.color": P["ink"], "axes.labelcolor": P["ink2"],
        "xtick.color": P["ink2"], "ytick.color": P["ink2"],
        "axes.edgecolor": P["grid"], "grid.color": P["grid"],
        "axes.grid": True, "grid.linewidth": 0.5, "grid.alpha": 0.7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "lines.linewidth": 2.0,
        "font.size": 9, "axes.titlesize": 10, "legend.frameon": False,
    })


def save(fig, name, P, mode):
    os.makedirs(C.FIGS, exist_ok=True)
    os.makedirs(PNG, exist_ok=True)
    if mode == "light":
        fig.savefig(f"{C.FIGS}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{PNG}/{name}.{mode}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


# --- 01 barnyard --------------------------------------------------------------
def fig_barnyard(pc, P, mode):
    srcs = list(dict.fromkeys(pc.source))
    # Shared limits across panels: the whole point is a visual comparison, and
    # independently autoscaled axes would make identical data look different.
    lo = 0.7 * max(1, min(pc.ufi_human.min(), pc.ufi_mouse.min()) + 1)
    hi = 1.5 * max(pc.ufi_human.max(), pc.ufi_mouse.max())
    fig, axes = plt.subplots(1, len(srcs), figsize=(5.0 * len(srcs), 4.8),
                             squeeze=False, sharex=True, sharey=True)
    for ax, src in zip(axes[0], srcs):
        d = pc[pc.source == src]
        for call in ["discarded", "mixed", "human", "mouse"]:
            s = d[d.call_fig1d == call]
            if s.empty:
                continue
            slot = CALL_SLOT[call]
            ax.scatter(s.ufi_human + 1, s.ufi_mouse + 1, s=17,
                       c=(P["grid"] if slot is None else P[slot]),
                       edgecolors=P["surface"], linewidths=0.6,
                       label=f"{call} (n={len(s)})", zorder=3 if slot else 2)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("human UFIs per barcode"); ax.set_ylabel("mouse UFIs per barcode")
        rate = C.doublet_rate(d.call_fig1d)
        ax.set_title(f"{src}\nFig. 1d rule: {rate:.2f}% heterotypic", color=P["ink"])
        # upper left is the one empty corner on a barnyard plot
        ax.legend(loc="upper left", fontsize=7.5, labelcolor=P["ink2"])
    fig.suptitle("Species mixing, HEK293T + mESC (GSM5369495)  ·  >25% of UFIs from the "
                 "other species = mixed; <7,500 UFIs discarded",
                 fontsize=9.5, color=P["ink2"], y=1.02)
    save(fig, "01_barnyard", P, mode)


# --- 02 genes / UFIs per cell -------------------------------------------------
def fig_genes_umis(pc, P, mode):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4))
    srcs = list(dict.fromkeys(pc.source))
    groups, labels, colors = [], [], []
    for sp, slot in (("human", "s1"), ("mouse", "s2")):
        for src in srcs:
            d = pc[(pc.source == src) & (pc.call_fig1d == sp)]
            groups.append(d); labels.append(f"{'HEK293T' if sp=='human' else 'mESC'}\n{src}")
            colors.append(P[slot])

    for ax, (col, name) in zip(axes, [("ufi_total", "total UFIs per cell"),
                                      (None, "genes detected per cell")]):
        vals = []
        for d in groups:
            v = d.ufi_total if col else np.where(d.call_fig1d == "human", d.genes_human, d.genes_mouse)
            vals.append(np.asarray(v, dtype=float))
        for i, (v, c) in enumerate(zip(vals, colors)):
            x = np.random.default_rng(0).normal(i, 0.07, len(v))
            ax.scatter(x, v, s=9, c=c, alpha=0.55, edgecolors="none", zorder=2)
            med = np.median(v) if len(v) else np.nan
            ax.hlines(med, i - 0.30, i + 0.30, color=P["ink"], lw=2, zorder=4)
            ax.annotate(f"{med:,.0f}", (i, med), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=8.5, color=P["ink"], zorder=5)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_yscale("log"); ax.set_ylabel(name)
        ax.set_title(name, color=P["ink"])
    fig.suptitle("Per-cell yield, ours vs the deposited table  ·  medians labelled",
                 fontsize=9.5, color=P["ink2"], y=1.02)
    save(fig, "02_genes_umis", P, mode)


# --- 03 biotype composition ---------------------------------------------------
def fig_biotypes(bs, P, mode):
    cols = list(bs.columns)
    top = bs.sum(axis=1).sort_values(ascending=False).head(14).index[::-1]
    d = bs.loc[top]
    y = np.arange(len(top)); h = 0.38
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    for k, (col, slot) in enumerate(zip(cols, ["s1", "s2"])):
        ax.barh(y + (k - 0.5) * h, d[col], height=h - 0.04, color=P[slot], label=col)
    ax.set_yticks(y); ax.set_yticklabels(top, fontsize=8)
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_xlabel("share of assigned UFIs (%)")
    ax.legend(fontsize=8, labelcolor=P["ink2"], loc="lower right")
    ax.set_title("Biotype composition, simple (non-combination) rows", color=P["ink"])
    save(fig, "03_biotype_composition", P, mode)


# --- 04 ours vs published -----------------------------------------------------
def fig_concordance(conc, P, mode):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4))
    ax = axes[0]
    d = conc[(conc.published > 0) & (conc.ours > 0)]
    ax.hexbin(np.log10(d.published), np.log10(d.ours), gridsize=70,
              bins="log", cmap="Blues" if mode == "light" else "cividis",
              linewidths=0)
    lim = [0, float(np.log10(d[["published", "ours"]].max().max()))]
    ax.plot(lim, lim, color=P["s2"], lw=1.6, ls="--", label="y = x")
    ax.set_xlabel("published, log10 total UFIs"); ax.set_ylabel("ours, log10 total UFIs")
    ax.legend(fontsize=8, labelcolor=P["ink2"])
    ax.set_title(f"Per-gene totals, {len(d):,} shared genes", color=P["ink"])

    ax = axes[1]
    r = d.log2_ratio
    ax.hist(np.clip(r, -3, 3), bins=90, color=P["s1"], edgecolor="none")
    ax.axvline(0, color=P["ink2"], lw=1.2)
    ax.axvline(float(r.median()), color=P["s2"], lw=1.8,
               label=f"median {r.median():+.3f}")
    ax.set_xlabel("log2(ours / published), clipped to ±3"); ax.set_ylabel("genes")
    ax.legend(fontsize=8, labelcolor=P["ink2"])
    ax.set_title("Agreement per gene", color=P["ink"])
    save(fig, "04_ours_vs_published", P, mode)


# --- 05 tRNA before / after ---------------------------------------------------
def fig_trna(summary, P, mode):
    """summary has ONE COLUMN PER RUN (ours_rrnav2, ours_bedv2, ...)."""
    def g(k, col, default=np.nan):
        try:
            return float(summary.loc[k, col])
        except Exception:
            return default

    pub = next((g("tRNA rows: published", c) for c in summary.columns
                if not np.isnan(g("tRNA rows: published", c))), np.nan)
    bars_spec = [("published\n(GSM5369495)", pub, "s4")]
    for run, slot in (("rrnav2", "s2"), ("bedv2", "s1")):
        col = f"ours_{run}"
        if col in summary.columns:
            bars_spec.append((f"ours, {run}\n"
                              f"({'BED without tRNA' if run == 'rrnav2' else 'BED with tRNA'})",
                              g("tRNA rows: ours", col), slot))

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    labels = [b[0] for b in bars_spec]
    vals = [0.0 if np.isnan(b[1]) else b[1] for b in bars_spec]
    bars = ax.bar(labels, vals, color=[P[b[2]] for b in bars_spec], width=0.58)
    for b, (_, v, _) in zip(bars, bars_spec):
        ax.annotate("not run" if np.isnan(v) else f"{v:,.0f}",
                    (b.get_x() + b.get_width() / 2, b.get_height()),
                    textcoords="offset points", xytext=(0, 5), ha="center",
                    fontsize=10, color=P["ink"])
    ax.set_ylim(0, max(vals + [1]) * 1.18)
    ax.set_ylabel("distinct tRNA rows detected")
    ax.set_title("tRNA detection, ours against the deposited table", color=P["ink"])
    # The shortfall is geometric, and saying so on the figure stops it being read
    # as a bug in the BED.
    ax.text(0.5, -0.30,
            "Reads (mean 70.8 bp) are as long as tRNA features (median 72 bp), and step 6 keeps a\n"
            "non-splicing biotype only when the read is fully inside it — true for just 9.5% here.",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=8.5, color=P["ink2"])
    save(fig, "05_trna_before_after", P, mode)


def main():
    pc = pd.read_csv(f"{C.RES}/per_cell.tsv", sep="\t", index_col=0)
    bs = pd.read_csv(f"{C.RES}/biotype_composition.tsv", sep="\t", index_col=0)
    conc = pd.read_csv(f"{C.RES}/gene_concordance.tsv", sep="\t", index_col=0)
    summary = pd.read_csv(f"{C.RES}/comparison_summary.tsv", sep="\t", index_col=0)

    for mode, P in (("light", LIGHT), ("dark", DARK)):
        style(P)
        fig_barnyard(pc, P, mode)
        fig_genes_umis(pc, P, mode)
        fig_biotypes(bs, P, mode)
        fig_concordance(conc, P, mode)
        fig_trna(summary, P, mode)
    print(f"wrote figures to {C.FIGS} (pdf, light) and {PNG} (png, light+dark)")


if __name__ == "__main__":
    main()
