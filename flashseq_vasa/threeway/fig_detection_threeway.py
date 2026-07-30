#!/usr/bin/env python3
"""Three-way depth-matched gene-detection saturation figure.

Reads only the tables written by mk_detection_threeway.py, so the figure cannot
disagree with them. Headline slice: table='uniagg', universe='shared'.

Panel a  saturation curves, five tracks. Per-unit curves are drawn faint behind
         each median. A SOLID median means every unit of that track contributes
         at that depth; a DASHED median with open markers means only the units
         that natively reach it do -- i.e. a depth-selected subset (the size of
         that selection effect is `selection_bias_genes` in the tracks table).
         The strip below panel a shows each track's native in-scope depth range,
         which is the asymmetry that makes this necessary: the published plate's
         mouse cells and the own plate's cells do not overlap in depth at all.
Panel b  the same curves expressed as a difference from the published plate, so
         the reader can see the sign and size of every gap on one axis. The grey
         band is where no published-plate cell is deep enough to compare against.
"""
import os
import sys

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

RES = sys.argv[1] if len(sys.argv) > 1 else \
    "/nemo/lab/turnerj/working/guangxin/vasaseq/res/threeway"
OUTPNG = "%s/detection_threeway_saturation.png" % RES

det = pd.read_csv("%s/detection_threeway.tsv" % RES, sep="\t")
sat = pd.read_csv("%s/detection_threeway_saturation.tsv" % RES, sep="\t")
trk = pd.read_csv("%s/detection_threeway_tracks.tsv" % RES, sep="\t")
for d in (det, sat):
    d["key"] = d.dataset + "|" + d.arm + "|" + d.unit

TRACKS = {
    "VASA published (mouse cells)":
        lambda d: (d.dataset == "VASA_published") & (d.label_fig1d == "mouse"),
    "VASA own plate": lambda d: d.dataset == "VASA_own",
    "FLASH-seq native": lambda d: (d.arm == "native") & (d.qc_verdict != "exclude"),
    "FLASH-seq VASA-trimmed": lambda d: (d.arm == "vasalen") & (d.qc_verdict != "exclude"),
    "FLASH-seq 30 pg (trimmed)": lambda d: (d.arm == "vasalen") & (d.group == "30 pg"),
}
NAME = {k: k.replace(" (mouse cells)", ", mouse cells").replace(" (trimmed)", ", trimmed")
        for k in TRACKS}
# One hue per entity, threaded across both panels. VASA arms are the focal pair
# (saturated blue/red); FLASH-seq arms are lower visual weight. No red/green pair.
COL = {"VASA published (mouse cells)": "#4C72B0", "VASA own plate": "#C44E52",
       "FLASH-seq native": "#8C8C8C", "FLASH-seq VASA-trimmed": "#DD8452",
       "FLASH-seq 30 pg (trimmed)": "#937860"}
LW = {"VASA own plate": 2.4, "VASA published (mouse cells)": 2.4}
# label anchors: (depth, dx_pt, dy_pt, ha) placed where that track is separated
ANCH = {"VASA published (mouse cells)": (2e5, -6, 17, "right"),
        "VASA own plate": (5e6, 2, 12, "right"),
        "FLASH-seq VASA-trimmed": (2e7, 7, 4, "left"),
        "FLASH-seq native": (2e7, 7, -5, "left"),
        "FLASH-seq 30 pg (trimmed)": (2e6, 4, -13, "left")}
GREY = "#6E6E6E"
XL = (8.2e3, 2.55e7)

hl = (det.table == "uniagg") & (det.universe == "shared")
su = sat[(sat.table == "uniagg") & (sat.universe == "shared")]
tk = trk[(trk.table == "uniagg") & (trk.universe == "shared")]
piv = tk.pivot_table(index="depth", columns="track", values="median_genes")
REF = "VASA published (mouse cells)"
xmax_pub = float(piv[REF].dropna().index.max())

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                     "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.8, "figure.dpi": 110,
                     "font.family": "sans-serif"})

fig = plt.figure(figsize=(7.4, 6.9))
# bottom margin and panel-b y-range are deliberately generous: this figure is
# rendered on a host whose default sans-serif is wider than the authoring
# machine's, and a tight layout that passed locally failed the overlap assertion
# there (the 30 pg label collided with the x-axis label). Clearance, not tuning.
gs = fig.add_gridspec(4, 1, height_ratios=[2.5, 0.40, 0.22, 1.15], hspace=0.0,
                      left=0.150, right=0.795, top=0.905, bottom=0.105)
ax, axr, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[3])


def letter(a, s):
    a.text(-0.085, 1.02, s, transform=a.transAxes, fontsize=11, fontweight="bold",
           va="bottom", ha="right")


# --- panel a --------------------------------------------------------------
for t, fn in TRACKS.items():
    keys = set(det[hl & fn(det)].key)
    for _, g in su[su.key.isin(keys)].groupby("key"):
        ax.plot(g.depth, g.genes, color=COL[t], lw=0.35, alpha=0.13, zorder=1,
                solid_capstyle="butt")
    tt = tk[tk.track == t]
    solid = tt[tt.support == "all_units"]
    ax.plot(solid.depth, solid.median_genes, color=COL[t], lw=LW.get(t, 1.7),
            marker="o", ms=3.6, zorder=4)
    sub = tt[tt.depth >= solid.depth.max()]
    ax.plot(sub.depth, sub.median_genes, color=COL[t], lw=LW.get(t, 1.7),
            ls=(0, (3.2, 1.6)), marker="o", ms=3.6, mfc="white", mew=1.0, zorder=4)
    xa, dx, dy, ha = ANCH[t]
    ya = float(tt[tt.depth == int(xa)].median_genes.iloc[0])
    ax.annotate(NAME[t], (xa, ya), xytext=(dx, dy), textcoords="offset points",
                ha=ha, va="center", color=COL[t], fontsize=7.5, zorder=6,
                bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.72))
ax.set_xscale("log")
ax.set_xlim(*XL)
ax.set_ylim(3550, 19300)
ax.set_xticklabels([])
ax.set_yticks([4000, 8000, 12000, 16000])
ax.set_yticklabels(["4k", "8k", "12k", "16k"])
ax.set_ylabel("protein-coding genes detected")
ax.set_title("At every matched depth the own VASA plate detects more genes than the "
             "published\nplate, and the gap widens with depth", loc="left", pad=8)
ax.text(0.014, 0.965, "higher = better", transform=ax.transAxes, color=GREY,
        fontsize=7, va="top")
ax.plot([], [], color=GREY, lw=1.5, label="every unit contributes")
ax.plot([], [], color=GREY, lw=1.5, ls=(0, (3.2, 1.6)), marker="o", ms=3.6,
        mfc="white", mew=1.0, label="only units this deep")
ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.02), frameon=False, fontsize=7,
          handlelength=2.4, borderpad=0.1, labelspacing=0.28)
letter(ax, "a")

# --- native depth strip ---------------------------------------------------
for i, t in enumerate(list(TRACKS)[::-1]):
    n = det[hl & TRACKS[t](det)].total_reads
    axr.plot([n.min(), n.max()], [i, i], color=COL[t], lw=3.0,
             solid_capstyle="butt", alpha=0.9)
    axr.plot([n.median()], [i], marker="|", color="white", ms=5.5, mew=1.2)
axr.set_xscale("log")
axr.set_xlim(*XL)
axr.set_ylim(-1.0, 5.2)
axr.set_yticks([])
axr.set_xticklabels([])
axr.tick_params(axis="x", bottom=False, which="both")
for s in ("left", "right", "top", "bottom"):
    axr.spines[s].set_visible(False)
axr.text(-0.012, 0.5, "native depth", transform=axr.transAxes, ha="right",
         va="center", fontsize=7.4)
axr.text(1.008, 0.5, "each unit's own depth\nrange (tick = median)",
         transform=axr.transAxes, color=GREY, fontsize=6.9, va="center")

# --- panel b --------------------------------------------------------------
for t, dy in [("VASA own plate", 0), ("FLASH-seq VASA-trimmed", 5),
              ("FLASH-seq native", -7), ("FLASH-seq 30 pg (trimmed)", 0)]:
    d = piv[[REF, t]].dropna()
    v = d[t] - d[REF]
    ax2.plot(d.index, v, color=COL[t], lw=1.7, marker="o", ms=3.6)
    ax2.annotate(NAME[t], (d.index[-1], v.iloc[-1]), xytext=(7, dy),
                 textcoords="offset points", va="center", color=COL[t], fontsize=7.3)
ax2.axvspan(xmax_pub, XL[1], color="#F2F2F2", zorder=0)
ax2.text(4.5e6, 800, "no published-plate cell\nis this deep", color=GREY,
         fontsize=6.9, ha="center", va="center")
ax2.axhline(0, color=COL[REF], lw=1.5, zorder=2)
ax2.annotate("published plate = 0", (9.4e3, 0), xytext=(0, -11),
             textcoords="offset points", color=COL[REF], fontsize=7.2)
ax2.set_xscale("log")
ax2.set_xlim(*XL)
ax2.set_ylim(-1750, 3450)
ax2.set_xticks([1e4, 1e5, 1e6, 1e7])
ax2.set_xticklabels(["10k", "100k", "1M", "10M"])
ax2.set_yticks([-1000, 0, 1000, 2000, 3000])
ax2.set_yticklabels(["\u22121k", "0", "+1k", "+2k", "+3k"])
ax2.set_xlabel("matched depth (reads on in-scope entries, log scale)")
ax2.set_ylabel("genes minus\npublished plate", labelpad=8)
ax2.set_title("Both VASA plates sit above the 30 pg rung; only the FLASH-seq arms "
              "cross each other,\nand only beyond 13M reads", loc="left", pad=5)
letter(ax2, "b")

# --- claim-titles must be true of every track DRAWN IN THAT PANEL ---------
# A reviewer caught two titles that were true of the three-way comparison but
# false of the 30 pg and native arms drawn beside it. Both titles are therefore
# asserted here against the plotted data, so a title can no longer outrun the
# panel it sits on.
cross = pd.read_csv("%s/detection_threeway_crossings.tsv" % RES, sep="\t")
cs = cross[(cross.table == "uniagg") & (cross.universe == "shared")]

# panel a: "own plate detects more genes than the published plate at every
#           matched depth, and the gap widens with depth"
gap = (piv["VASA own plate"] - piv[REF]).dropna()
assert (gap > 0).all(), "panel-a title false: own plate not above published at %s" % \
    list(gap[gap <= 0].index)
assert (gap.sort_index().diff().dropna() > 0).all(), \
    "panel-a title false: gap does not widen monotonically:\n%s" % gap
# panel b: "both VASA plates sit above the 30 pg rung; only the FLASH-seq arms
#           cross each other, and only beyond 13M reads"
p30 = "FLASH-seq 30 pg (trimmed)"
for v in ("VASA own plate", REF):
    d = piv[[v, p30]].dropna()
    assert (d[v] > d[p30]).all(), "panel-b title false: %s not above 30 pg" % v
xr = cs[cs.crosses == "yes"]
assert set(xr.track_a) <= {"FLASH-seq native", "FLASH-seq VASA-trimmed"} and \
    set(xr.track_b) == {p30}, \
    "panel-b title false: a non-FLASH-seq pair crosses:\n%s" % xr[["track_a", "track_b"]]
assert float(xr.crossing_depth.min()) > 13e6, \
    "panel-b title false: earliest crossing is %.0f, not beyond 13M" % xr.crossing_depth.min()
print("title claims verified: own-vs-published gap positive and widening at all %d "
      "rungs; %d crossing(s), all FLASH-seq-internal, earliest at %s reads"
      % (len(gap), len(xr), format(int(xr.crossing_depth.min()), ",")))

fig.canvas.draw()   # realise text extents; the save happens after the nudge loop

# --- render-then-verify: geometric overlap check --------------------------
# A rotated y-axis label sitting on its own tick labels is not a finding, so those
# two are exempt; everything else must be clear. Series annotations are nudged
# upward and re-measured rather than hard-coded to one machine's font metrics --
# the figure must come out clean on whatever sans-serif the host resolves.
AXIS_LABELS = {"protein-coding genes detected", "genes minus\npublished plate"}


def collisions():
    r = fig.canvas.get_renderer()
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    return [(a, b) for i, (a, ba) in enumerate(tx) for b, bb in tx[i + 1:]
            if ba.overlaps(bb) and a.get_text() not in AXIS_LABELS
            and b.get_text() not in AXIS_LABELS]


for attempt in range(6):
    bad = collisions()
    if not bad:
        break
    for a, b in bad:                       # lift whichever is an annotation
        for t in (a, b):
            if isinstance(t, mpl.text.Annotation):
                dx, dy = t.get_position() if not hasattr(t, "xyann") else t.xyann
                t.set_position((dx, dy + 3.0))
                break
    fig.canvas.draw()
bad = collisions()
assert not bad, "figure still has overlapping labels after nudging: %s" % \
    [(a.get_text()[:26], b.get_text()[:26]) for a, b in bad]

fig.savefig(OUTPNG, dpi=300)
print("wrote %s (%.0f kB); overlap check clean after %d nudge pass(es)"
      % (OUTPNG, os.path.getsize(OUTPNG) / 1e3, attempt))
