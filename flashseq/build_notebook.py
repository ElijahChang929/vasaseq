#!/usr/bin/env python3
"""Generate flashseq_qc.ipynb.

The notebook is generated rather than hand-edited so that its source stays
reviewable as plain text -- a committed .ipynb diffs badly, and this repository
cares about being able to read what changed. Edit THIS file, then:

    source code/flashseq/config.sh
    fs_python code/flashseq/build_notebook.py code/flashseq/flashseq_qc.ipynb

Any execution outputs are discarded by regenerating, which is intended: the
rendered copy with outputs belongs in res/flashseq/flashseq_qc.html, not in git.
"""
import json, sys, pathlib

md = lambda s: {"cell_type": "markdown", "metadata": {}, "source": s.rstrip("\n").split("\n")}
def code(s):
    lines = s.strip("\n").split("\n")
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": [l + "\n" for l in lines[:-1]] + [lines[-1]]}

cells = []

cells.append(md("""
# FLASH-seq RN26038 — data quality

**Run:** `20260325_LH00442_0237_B23GT7GLT3` lane 7, NovaSeq X, 10 individually-indexed
libraries `ZHA8833A1..A10`, paired-end 151 bp, mouse.

**Pipeline:** nf-core/rnaseq 3.22.2, `--aligner star_rsem`, GRCm39 + Ensembl release-116,
STAR index built at `--sjdbOverhang 150`. Driver: `data/flashseq/nfcore_rnaseq_all.sh`.

**Design — this is the thing to hold onto.** The run is not ten comparable samples. It is an
**input-amount titration in duplicate**, spanning 30 ng down to 30 pg: a 1000-fold range.
A mammalian cell carries roughly 10–30 pg of total RNA, so the 30 pg rung is the
single-cell-equivalent one, and it is the rung that matters for comparing against VASA-seq.

The mapping comes from the RN26038 LIMS sheet and is recorded in
`code/flashseq/sample_metadata.tsv`; nothing inside the run directory carries it.
"""))

cells.append(code("""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

ROOT = Path(os.environ.get("FS_ROOT", "/nemo/lab/turnerj/working/guangxin/vasaseq"))
CODE = ROOT / "code/flashseq"
RES  = ROOT / "res/flashseq"
FIG  = RES / "figures"; FIG.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 9, "axes.titlesize": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
    "grid.linewidth": 0.5, "legend.frameon": False,
})

meta = pd.read_csv(CODE / "sample_metadata.tsv", sep="\\t", comment="#")
# Everything is ordered by input amount, high to low, throughout this notebook.
meta = meta.sort_values("input_pg", ascending=False).reset_index(drop=True)
ORDER  = meta["library"].tolist()
LEVELS = meta.drop_duplicates("input_amount")["input_amount"].tolist()

# One colour per input level, dark (high input) to light (low input).
CMAP = plt.get_cmap("viridis")
LEVEL_COLOR = {lv: CMAP(i / max(len(LEVELS) - 1, 1) * 0.85) for i, lv in enumerate(LEVELS)}
lib2level = dict(zip(meta["library"], meta["input_amount"]))
def colors(libs):
    return [LEVEL_COLOR[lib2level[l]] for l in libs]

meta
"""))

cells.append(md("""
## 1. Sequencing and alignment

Everything below comes from `res/flashseq/qc_summary.tsv`, written by
`code/flashseq/00_collect_qc.py` out of the MultiQC tables.
"""))

cells.append(code("""
qc = pd.read_csv(RES / "qc_summary.tsv", sep="\\t")
qc = qc.set_index("library").loc[ORDER].reset_index()

show = qc[["library", "input_amount", "replicate", "well", "raw_pairs",
           "star_uniquely_mapped_pct", "star_multimapped_pct",
           "star_unmapped_tooshort_pct", "rsem_alignable_pct",
           "picard_dup_pct", "qualimap_intergenic_pct", "qualimap_5_3_bias",
           "pct_MT", "insert_size_avg"]].copy()
show["raw_pairs"] = (show["raw_pairs"] / 1e6).round(1)
show = show.rename(columns={"raw_pairs": "raw_pairs_M"})
show.style.format(precision=2).background_gradient(
    subset=["star_uniquely_mapped_pct", "rsem_alignable_pct"], cmap="RdYlGn")
"""))

cells.append(code("""
fig, axes = plt.subplots(1, 3, figsize=(11, 3.1))
panels = [
    ("star_uniquely_mapped_pct", "STAR uniquely mapped (%)"),
    ("star_unmapped_tooshort_pct", "STAR unmapped, too short (%)"),
    ("picard_dup_pct", "duplicate reads (%)"),
]
for ax, (col, title) in zip(axes, panels):
    ax.bar(range(len(qc)), qc[col], color=colors(qc["library"]))
    ax.set_xticks(range(len(qc)))
    ax.set_xticklabels([l.replace("ZHA8833", "") for l in qc["library"]], rotation=0)
    ax.set_title(title)
    ax.set_xlabel("library")
handles = [plt.Rectangle((0, 0), 1, 1, color=LEVEL_COLOR[lv]) for lv in LEVELS]
axes[-1].legend(handles, LEVELS, title="input", loc="upper left", fontsize=7,
                title_fontsize=7)
fig.suptitle("Alignment quality falls off only at the two lowest input rungs", y=1.04)
fig.savefig(FIG / "01_alignment.pdf"); plt.show()
"""))

cells.append(md("""
Mechanically the run is sound: quality scores are clean, no reads flagged as poor quality,
strandedness auto-inferred `unstranded` on all ten (correct for FLASH-seq), mitochondrial
fraction 1–3 %, intergenic 1–2 %, 5′/3′ bias 1.14–1.31.

The two 60 pg libraries and the 30 pg pair carry visibly more unmapped "too short" reads.
Sections 4 and 5 show that this has **two separate causes** that happen to land on
neighbouring rungs — a contaminant confined to two wells, and poly-G that scales with input.

## 2. Sensitivity — genes detected across the titration

The obvious worry with a gene-detection drop is that it is really a depth difference. It is
not: `03_gene_detection.py` also reports the expected number of genes surviving binomial
thinning to the shallowest library's total count, and the split is unchanged.
"""))

cells.append(code("""
det = pd.read_csv(RES / "gene_detection.tsv", sep="\\t")
det = det.set_index("library").loc[ORDER].reset_index()
det[["library", "input_amount", "replicate", "total_counts",
     "genes_ge1", "genes_rarefied", "genes_tpm_gt1"]]
"""))

cells.append(code("""
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))

ax = axes[0]
for lv in LEVELS:
    sub = det[det["input_amount"] == lv]
    ax.scatter(sub["input_pg"], sub["genes_ge1"], s=55, color=LEVEL_COLOR[lv],
               label=lv, zorder=3, edgecolor="white", linewidth=0.8)
    ax.scatter(sub["input_pg"], sub["genes_rarefied"], s=55, facecolor="none",
               edgecolor=LEVEL_COLOR[lv], linewidth=1.4, zorder=3)
ax.set_xscale("log")
ax.set_xlabel("input (pg, log scale)"); ax.set_ylabel("genes detected")
ax.set_title("filled = as sequenced;  open = rarefied to common depth")
ax.legend(title="input", fontsize=7, title_fontsize=7, loc="lower right")

ax = axes[1]
w = 0.38
x = np.arange(len(det))
ax.bar(x - w/2, det["genes_ge1"], w, color=colors(det["library"]), label="as sequenced")
ax.bar(x + w/2, det["genes_rarefied"], w, color=colors(det["library"]),
       alpha=0.45, hatch="//", label="rarefied")
ax.set_xticks(x); ax.set_xticklabels([l.replace("ZHA8833", "") for l in det["library"]])
ax.set_ylabel("genes with >=1 count"); ax.set_xlabel("library")
ax.legend(fontsize=7)
ax.set_title("depth is not the explanation")

fig.suptitle("Gene detection holds flat to 1.5 ng, then halves", y=1.03)
fig.savefig(FIG / "02_sensitivity.pdf"); plt.show()
"""))

cells.append(md("""
## 3. Replicate concordance — where the titration actually breaks

Pearson *r* of log2(TPM+1) between the two replicates of each rung, over genes expressed in
both (TPM > 1).
"""))

cells.append(code("""
conc = pd.read_csv(RES / "replicate_concordance.tsv", sep="\\t")
conc
"""))

cells.append(code("""
tpm = pd.read_csv(ROOT / "data/flashseq/results/star_rsem/rsem.merged.gene_tpm.tsv",
                  sep="\\t", index_col=0).drop(columns=["transcript_id(s)"])

n = len(conc)
fig, axes = plt.subplots(1, n, figsize=(2.5 * n, 2.8), sharex=True, sharey=True)
for ax, (_, r) in zip(np.atleast_1d(axes), conc.iterrows()):
    a, b = np.log2(tpm[r["rep1"]] + 1), np.log2(tpm[r["rep2"]] + 1)
    keep = (tpm[r["rep1"]] > 1) & (tpm[r["rep2"]] > 1)
    ax.scatter(a[keep], b[keep], s=1.2, alpha=0.16,
               color=LEVEL_COLOR[r["input_amount"]], rasterized=True)
    lim = (0, float(max(a.max(), b.max())) * 1.02)
    ax.plot(lim, lim, color="0.35", lw=0.7, ls="--")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_title(f"{r['input_amount']}\\nr = {r['pearson_log2tpm']:.4f}")
    ax.set_xlabel("rep 1, log2(TPM+1)")
np.atleast_1d(axes)[0].set_ylabel("rep 2, log2(TPM+1)")
fig.suptitle("Reproducibility is flat to 1.5 ng and collapses by 60 pg", y=1.10)
fig.savefig(FIG / "03_replicates.pdf"); plt.show()
"""))

cells.append(md("""
30 ng, 3 ng and 1.5 ng are **indistinguishable from one another** — *r* = 0.987 / 0.983 / 0.980
and ~14,380 genes at TPM > 1 in every one of the six libraries. Performance then collapses at
60 pg (*r* = 0.82).

**The design's blind spot:** there is no rung between 1.5 ng and 60 pg — a 25-fold gap — and
that gap is exactly where the method stops working. If the point of the titration was to find
the sensitivity floor, the floor is somewhere inside an interval the experiment did not
sample. Filling it (say 500 pg / 150 pg) is the obvious follow-up.

## 4. rRNA — the pipeline's number is wrong, and wrong in a knowable direction

nf-core reports rRNA as the read fraction over genes annotated `gene_biotype "rRNA"` in the
Ensembl GTF. On GRCm39 that annotation holds 354 rRNA genes which are **essentially all 5S**
(`n-R5s*`) plus one `Rn18s-rs5` relic — there is **no `Rn45s`, `Rn28s` or `Rn5-8s` gene at
all**, because the rDNA array is collapsed out of the primary assembly. So the reported
0.7–1.4 % is a measurement of 5S, not of rRNA.

This is the same defect already found and fixed on the VASA side of this repository
(see `CLAUDE.md`, "Reference provenance"). `01_rrna_kmer_screen.py` re-measures against
`unique_rRNA_mouse.v2.fa`, which carries the NCBI 47S pre-rRNA unit `BK000964.3:1-13403`,
and reports the 47S and Ensembl contributions separately.
"""))

cells.append(code("""
rr = pd.read_csv(RES / "rrna_kmer.tsv", sep="\\t").set_index("library").loc[ORDER].reset_index()
rr = rr.merge(qc[["library", "nfcore_biotype_rRNA_pct", "input_amount"]], on="library")
rr[["library", "input_amount", "nfcore_biotype_rRNA_pct",
    "pct_rRNA_47S", "pct_rRNA_ensembl_only", "pct_rRNA_total"]].round(2)
"""))

cells.append(code("""
fig, ax = plt.subplots(figsize=(7.2, 3.3))
x = np.arange(len(rr)); w = 0.36
ax.bar(x - w/2, rr["nfcore_biotype_rRNA_pct"], w, color="0.68",
       label="nf-core: Ensembl gene_biotype rRNA")
ax.bar(x + w/2, rr["pct_rRNA_ensembl_only"], w, color="#7f9fbf",
       label="k-mer: Ensembl records only")
ax.bar(x + w/2, rr["pct_rRNA_47S"], w, bottom=rr["pct_rRNA_ensembl_only"],
       color="#b3402f", label="k-mer: 47S unit (BK000964.3)")
ax.set_xticks(x); ax.set_xticklabels([l.replace("ZHA8833", "") for l in rr["library"]])
ax.set_ylabel("% of reads"); ax.set_xlabel("library")
ax.set_title("What the annotation-based rRNA number structurally cannot see")
ax.legend(fontsize=7, loc="upper left")
fig.savefig(FIG / "04_rrna.pdf"); plt.show()

share = 100 * rr["pct_rRNA_47S"].sum() / rr["pct_rRNA_total"].sum()
ratio = rr["pct_rRNA_total"].sum() / rr["nfcore_biotype_rRNA_pct"].sum()
print(f"{share:.1f}% of the rRNA signal comes from the 47S unit Ensembl does not annotate.")
print(f"The honest number is {ratio:.1f}x the one nf-core reports.")
print("\\nThis is a LOWER bound: exact 31-mer matching, sampled every 10 bases, no mismatches.")
print("For a number directly comparable with VASA, re-measure with own_version/ribo-bwamem.sh")
print("(bwa aln + bwa mem, both) so the two sides use the same method.")
"""))

cells.append(md("""
## 5. Contamination

FastQC flagged a sequence making up **15.9 % of library A8**. FastQC has no reference to check
against, so it reports only "No Hit". `02_contaminant_check.py` resolves every flagged
sequence against mouse GRCm39, the rRNA reference, ERCC92 and human GRCh38, then **re-counts
each one in every library** — the second step matters because FastQC only reports a sequence
once it exceeds ~0.1 % of a file, so it cannot distinguish "absent" from "below the floor".
"""))

cells.append(code("""
ov = pd.read_csv(RES / "overrepresented.tsv", sep="\\t")

# pct_in_library is max(pct_R1, pct_R2): a fragment appears on one mate, not
# both, so the larger of the two is that fragment's rate in the library.
peak = (ov.groupby(["sequence", "source", "gene", "locus"], dropna=False)
          ["pct_in_library"].max().reset_index()
          .sort_values("pct_in_library", ascending=False))
print(f"{len(peak)} distinct sequences; classification:")
display(peak["source"].value_counts().rename("n_sequences").to_frame())
peak.head(12)
"""))

cells.append(code("""
# _variant picks up the sequencing-error copy of the dominant fragment.
human = ov[ov["source"].str.startswith("human_GRCh38", na=False)].copy()
if len(human):
    gene = ", ".join(sorted(set(human["gene"].dropna()) - {""}))

    # One row per human sequence, one column per library. Do NOT sum down the
    # column: the R1 and R2 sequences of a single fragment are separate rows,
    # and adding them would count that fragment twice.
    wide = (human.pivot_table(index="sequence", columns="library",
                              values="pct_in_library")
                 .reindex(columns=ORDER))
    wide.index = [f"{s[:24]}..." for s in wide.index]

    # Total human content per library WITHOUT double-counting mates: sum pct_R1
    # only. Each fragment sits on exactly one mate, so its R2-mate row
    # contributes ~0 to pct_R1 and the sum comes out right.
    total = (human.groupby("library")["pct_R1"].sum()
                  .reindex(ORDER).rename("human_pct_total").reset_index())

    # The single most abundant human fragment, per library -- this is the one
    # FastQC reported at 15.9% in A8.
    top = wide.max(axis=0).rename("pct").reset_index()
    top = top.merge(meta[["library", "input_amount", "input_pg", "well"]], on="library")

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.4))
    ax = axes[0]
    ax.bar(range(len(top)), top["pct"], color=colors(top["library"]))
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels([f"{l.replace('ZHA8833','')}\\n{w}"
                        for l, w in zip(top["library"], top["well"])], fontsize=7)
    ax.set_ylabel("% of reads on its mate"); ax.set_xlabel("library / plate well")
    ax.set_title(f"human {gene}: adjacent wells G:1, H:1")

    ax = axes[1]
    ax.scatter(top["input_pg"], top["pct"], s=70, color=colors(top["library"]),
               edgecolor="white", linewidth=0.8, zorder=3)
    for _, r in top.iterrows():
        ax.annotate(r["library"].replace("ZHA8833", ""),
                    (r["input_pg"], r["pct"]), fontsize=6.5,
                    xytext=(4, 3), textcoords="offset points")
    ax.set_xscale("log"); ax.set_xlabel("input (pg, log scale)")
    ax.set_ylabel("% of reads"); ax.set_title("no trend with input amount")

    fig.suptitle("The contamination tracks wells, not input", y=1.05)
    fig.savefig(FIG / "05_contamination.pdf"); plt.show()

    print("total human content per library (sum of pct_R1, no mate double-counting):")
    display(total.merge(meta[["library", "input_amount", "well"]], on="library").round(3))
    print("\\nper sequence (pct_in_library = max of the two mates):")
    display(wide.round(3))
else:
    print("no human-assigned overrepresented sequences found")
"""))

cells.append(md("""
The 30 pg pair is **cleaner** than the 30 ng pair, so this is not "low input amplifies
background". A7 and A8 sit in plate wells **G:1 and H:1 — adjacent**; A9 and A10 are A:2 and
B:2. The pattern is localised contamination introduced during prep.

Independent corroboration from the STP's own `fastq_screen`, which never saw our pipeline:
"""))

cells.append(code("""
sc = [c for c in qc.columns if c.startswith("screen_")]
if sc:
    keep = [c for c in ["screen_MOUSE", "screen_HUMAN", "screen_RAT", "screen_VECTOR"] if c in sc]
    display(qc[["library", "input_amount", "well"] + keep].round(1))
else:
    print("fastq_screen output not found -- check FS_SCREEN in config.sh")
"""))

cells.append(md("""
### Technical artefacts: adapter read-through and poly-G

Neither is in any genome, so both are labelled by inspection rather than searched for --
otherwise they would be reported as mystery contaminants.

The FastQC-derived sequences are **not comparable across libraries**: adapter read-through
carries the library's own index, so every library produces a different 50-mer and FastQC only
flags those clearing its threshold. `02_contaminant_check.py` therefore also counts
**index-independent probes** — the shared adapter cores and a 30 nt G-run — in every library,
and those are what is plotted here.
"""))

cells.append(code("""
pr = ov[ov["source"].str.startswith("probe_", na=False)].copy()
if len(pr):
    pr["kind"] = pr["source"].str.replace("probe_", "", regex=False)
    tab = (pr.pivot_table(index="library", columns="kind", values="pct_in_library")
             .reindex(ORDER).fillna(0).reset_index())
    kinds = [c for c in tab.columns if c != "library"]
    kinds = [c for c in kinds if tab[c].max() > 0.05]   # drop probes that never fire

    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    x = np.arange(len(tab)); w = 0.8 / max(len(kinds), 1)
    for k, c in enumerate(kinds):
        ax.bar(x + k * w - 0.4 + w / 2, tab[c], w, label=c)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{l.replace('ZHA8833','')}\\n{a}"
                        for l, a in zip(tab["library"],
                                        meta.set_index("library").loc[ORDER, "input_amount"])],
                       fontsize=7)
    ax.set_ylabel("% of reads (max of the two mates)")
    ax.set_xlabel("library (high input -> low input)")
    ax.set_title("Technical artefacts, measured index-independently")
    ax.legend(fontsize=7)
    fig.savefig(FIG / "06_artefacts.pdf"); plt.show()
    display(tab.round(3))
else:
    print("no probe rows found -- re-run 02_contaminant_check.py")
"""))

cells.append(md("""
Adapter read-through is by far the larger effect: 19–36 % of R1 reads, strictly R1, rising as
input falls (shorter inserts from less material). This is a Nextera library — TrimGalore
auto-detected `CTGTCTCTTATA` and cutadapt reports adapter in 63 % of R1 reads — so
read-through is expected. It is invisible in the MultiQC report because **FastQC only
inspects the first 50 bp of a read** while read-through sits wherever the insert ends.

The question that matters is not how much artefact the raw data holds, but how much survives
trimming into STAR.
"""))

cells.append(code("""
te = pd.read_csv(RES / "trim_effect.tsv", sep="\\t")
seen = te[te["pct_raw"] > 0.01]

summary = (seen.groupby("kind")[["pct_raw", "pct_trimmed"]].mean()
               .assign(n=seen.groupby("kind").size()))
summary["pct_removed_of_total"] = (
    100 * (summary["pct_raw"] - summary["pct_trimmed"]) / summary["pct_raw"])
display(summary.round(3))

fig, ax = plt.subplots(figsize=(6.6, 3.1))
kinds = list(summary.index)
x = np.arange(len(kinds)); w = 0.36
ax.bar(x - w/2, summary["pct_raw"], w, label="raw", color="0.62")
ax.bar(x + w/2, summary["pct_trimmed"], w, label="after TrimGalore", color="#b3402f")
ax.set_xticks(x); ax.set_xticklabels(kinds)
ax.set_ylabel("mean % of reads in files where FastQC reported it")
ax.set_title("What TrimGalore actually removed")
ax.legend(fontsize=7)
fig.savefig(FIG / "07_trim_effect.pdf"); plt.show()
"""))

cells.append(code("""
# Why the adapter sequences survive: they begin MID-mosaic-end, so the pattern
# cutadapt anchors on is not present in them. This is checkable, not inferred.
PATTERN = "CTGTCTCTTATA"          # what TrimGalore auto-detected
adp = sorted(set(te.loc[te["kind"] == "adapter", "sequence"]))
print(f"cutadapt anchors on {PATTERN!r}; do the surviving read-through sequences contain it?")
for s in adp:
    print(f"  {PATTERN in s!s:<6} {s}")
print("\\nNone of them. Reads made entirely of read-through start mid-mosaic-end and cannot")
print("be matched, so they pass through untrimmed; reads with a partial 3' adapter (the 63%)")
print("are trimmed normally.")
"""))

cells.append(md("""
So of the three artefact classes, **none is meaningfully removed**: adapter read-through
survives when the read begins mid-mosaic-end, poly-G survives because TrimGalore's adapter
pass has no concept of it, and the human CALB1 fragments survive because they are real
sequence rather than artefact. All three reach STAR and feed its "unmapped, too short"
fraction.

Fixes, if a re-run is judged worthwhile: `--nextseq-trim 20` / `--2colour 20` for the poly-G,
and a second adapter pattern without the `CTGTCTCT` prefix for the read-through.

A note on reading `overrepresented.tsv` yourself: a fragment's R1 and R2 sequences are
separate rows, and `pct_in_library` is the larger of the two mates, so **summing that column
double-counts every pair** (it makes A8 look 33 % human rather than 16.6 %). Sum `pct_R1`
instead — the R2-mate rows contribute ~0 to it.
"""))

cells.append(md("""
Poly-G is the NovaSeq X two-colour dark-cycle artefact: with no signal, the base caller emits
G. TrimGalore's adapter pass does not remove it, so it survives into STAR and is part of what
gets discarded as "too short". Adding `--nextseq-trim 20` (cutadapt) or `--2colour 20`
(TrimGalore) to the trim step would recover some of those reads — worth quantifying on one
library before deciding whether a re-run is justified.

## 6. What this means for the VASA-seq comparison

VASA-seq operates at single-cell input, ~10–30 pg. The comparable FLASH-seq rungs are
therefore the two lowest — and they are also the two most compromised:

| rung | libraries | usable for the comparison? |
|---|---|---|
| 30 pg | A9, A10 | **yes** — the cleanest low-input pair; no detectable contamination |
| 60 pg | A7, A8 | **A7 with a caveat** (contaminant present); **A8 no** — a sixth of the library is not mouse |
| 1.5 ng and above | A1–A6 | not input-comparable to VASA; use as the ceiling reference |

So the honest single-cell-equivalent comparison uses **A9 and A10**, with A1–A6 showing what
FLASH-seq achieves when RNA is not limiting. Quote FLASH-seq's rRNA fraction from the
k-mer/bwa measurement in section 4, never the nf-core biotype number — the VASA figures were
computed against the 47S-containing reference, and mixing the two methods would compare
different quantities.

## 7. Caveats about the run itself

Two things to know before anything reuses this run's saved outputs.

**The saved GTF is truncated to chromosome 1.**
`results/genome/Mus_musculus.GRCm39.116.filtered.gtf` contains only chr1 (4,742 genes), and
its size is an exact multiple of 64 KiB — the signature of a partial copy. The sibling
`.filtered.bed` and the RSEM reference are both complete and all chromosomes, so **the
quantification in this run is unaffected**. But that GTF would silently produce a chr1-only
analysis if picked up by a later job.

**The STAR index was not saved** despite `--save_reference`: `genome/index/` holds only
`rsem/` and `salmon/`. The work directory `/nemo/lab/turnerj/scratch/zhangg/flashseq` has
already been deleted, so a re-run rebuilds the index (22 min, 65 GB peak, per
`pipeline_info/execution_trace_*.txt`).
"""))

cells.append(code("""
gtf = ROOT / "data/flashseq/results/genome/Mus_musculus.GRCm39.116.filtered.gtf"
bed = ROOT / "data/flashseq/results/genome/Mus_musculus.GRCm39.116.filtered.bed"

def contigs(path, col=0, limit=None):
    seen = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            if limit and i > limit:
                break
            if not line.startswith("#"):
                seen.add(line.split("\\t")[col])
    return seen

g, b = contigs(gtf), contigs(bed)
print(f"filtered.gtf contigs ({len(g)}): {sorted(g)}")
print(f"filtered.bed contigs ({len(b)}): {sorted(b)}")
print(f"\\ngtf size {gtf.stat().st_size:,} B; "
      f"exact multiple of 64 KiB: {gtf.stat().st_size % 65536 == 0}")
assert len(g) == 1, "the truncated-GTF caveat no longer holds -- update section 7"
print("\\nCAVEAT CONFIRMED: the saved GTF is chromosome 1 only.")

star_index = ROOT / "data/flashseq/results/genome/index/star"
print(f"saved STAR index present: {star_index.exists()}")
"""))

for i, c in enumerate(cells):          # nbformat >=4.5 requires cell ids
    c["id"] = f"c{i:02d}"

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
    },
    "nbformat": 4, "nbformat_minor": 5,
}
out = pathlib.Path(sys.argv[1])
out.write_text(json.dumps(nb, indent=1))
print("wrote", out, len(cells), "cells")
