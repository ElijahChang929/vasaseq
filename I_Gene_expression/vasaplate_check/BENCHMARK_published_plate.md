# Gold-standard benchmark — the published VASA-plate through this workflow

Purpose: establish whether the workflow used for the FLASH-seq ↔ VASA comparison
reproduces the **published** VASA-plate library, so that a three-way comparison
(published plate / own plate / FLASH-seq) rests on a validated anchor rather than
on an assumption.

Status 2026-07-29. **The published plate had already been re-run four times in
this repo before this session**, and run 4 already passed against the deposited
count table. This session's contribution is run 5: the one large remaining
discrepancy, traced but never fixed.

---

## What was already established (runs 1–4, `vasaplate_check/`)

Not re-derived here — read from `vasaplate_check/README.md` and
`res/vasaplate/comparison_summary.tsv`.

**The primary reference is the deposited count table, not the manuscript.** The
paper reports no barnyard number, no rRNA percentage and no tRNA percentage for
this library (`GSM5369495`); it appears mainly as a source of cells for the
benchmarking panels. So concordance is measured against
`GSM5369495_vasaplate_HEK293T-mESC_split_total.UFICounts.tsv.gz`.

| quantity | published | run 4 (`vasaplate_out_bedv2`) |
|---|---|---|
| simple gene rows | 77,207 | 80,234 |
| shared rows | — | 72,613 |
| per-gene Spearman r | — | **0.974** |
| per-gene Pearson r (log10) | — | **0.982** |
| median per-cell Pearson r | — | **0.982** |
| median log2(ours/published) | — | **0.000** |
| barcodes ≥ 7,500 UFIs | 353 | **353** |
| HEK293T median UFIs / genes | 190,110 / 18,316 | 185,227 / 18,757 |
| mESC median UFIs / genes | 75,233 / 11,852 | 72,893 / 11,930 |
| sncRNA share | 1.421% | 1.604% (paper states 1.4%) |
| **rRNA biotype** | — | **~600× higher — the open defect** |
| tRNA rows detected | 1,130 | 224 (148 shared) |

A median log2 ratio of exactly 0.000 with r = 0.974–0.982 across 72,613 genes is
a reproduction, not a resemblance. **The workflow is sound on everything except
rRNA.**

### Two discrepancies that are NOT workflow faults

- **tRNA, 224 vs 1,130.** Geometric, not annotational: tRNA features have median
  length 72 bp against a 70.8 bp mean read length, and stage 6 keeps a
  non-splicing biotype only if the read falls *entirely inside* the feature
  (`jS:IN`). 37.5% of reads are longer than the feature and can never be
  contained; 9.5% survive. The authors' tRNA annotation was never deposited and
  their Methods name no tRNA source, so why they recovered 5× more is a
  hypothesis (a precursor annotation with 5′/3′ flanks would fit reads inside),
  not a finding. **Do not present 224 as a reproduction of 1,130.**
- **The two doublet rules disagree ~6× on the same data.** Fig. 1d thresholds the
  UFI fraction, Methods p. 18 the gene fraction: 0.85% vs 4.82% on the published
  table itself. And the authors' own VASA-drop table gives 4.57%/5.42% against
  their printed 3.08%. A barnyard mismatch here is not automatically our fault.

---

## The open defect, and this session's run 5

The `rRNA` biotype is **~600× higher** than published, and it is **one locus**:
`ENSMUSG00000106106_CT010467.1_rRNA`, **95,823 UFIs against 72**. Diagnosed as
two defects that are each harmless alone:

**A.** `riboread-selection.py` emits `r0.seq`, and pysam returns the sequence *as
stored in the BAM* — for a reverse alignment, the reverse complement. So every
read kept by the stranded reprieve enters STAR flipped, inverting sense and
antisense for a stranded protocol.

**B.** One reference entry was stored antisense. `bedtools getfasta -s` gives the
sense strand of the *annotated gene*, which equals transcript sense only where
Ensembl put the gene on the transcript's strand. In Ensembl 99,
`ENSMUSG00000106106` is annotated opposite — so its entry was the 18S region
backwards. **1 of 915.**

Together: a genuine 18S read aligns reverse to that entry, is spared, is written
flipped by A, and STAR then maps it *sense* to the minus-strand gene, where it is
counted as rRNA.

The fix (`build_rrna_reference_mixed.sh` STEP 4, `ORIENT_TO_UNITS=yes`) aligns
every Ensembl entry to the NCBI 47S units, reverse-complements any returning flag
16, then re-aligns and requires zero remaining. It fixes the **class**, not the
instance. Output is `unique_rRNA_human_mouse.v3.fa` — **built in a previous
session and never run.**

### Precheck (Rule 2), job `1103ee07` — exit 0

Before committing 6–10 h over 384 cells, `vasaplate_check/05_v3_precheck.py`
asserted, read-only:

| check | result |
|---|---|
| v3 vs v2 is exactly the documented change | 921 names identical and in order; **920 sequences byte-identical, exactly 1 reverse-complemented, 0 changed any other way** |
| the flipped entry is the blamed locus | `mouse_ENSMUSG00000106106_rRNA(-)` ✓ |
| no antisense entries remain in v3 | 6 unit entries identified; **0** entries match the units better reversed than forward |
| stage 1–2 inputs resolve | 384 `_cbc.fastq.gz`, **0 broken symlinks**, 0.4–67.6 MB |
| bwa index complete and not stale | all five files present, none older than the fasta |

### Expected effect — the two metrics move in OPPOSITE directions

This is the prediction the run tests, recorded before it finished:

| metric | run 4 (v2) | predicted with v3 |
|---|---|---|
| step-3 depletion, 384 cells | 5.06% (9,176,750 / 181,287,059) | ~5.34% (**up** ~0.28 pt) |
| rRNA biotype share in the table | 0.201% (95,823 / 47,732,501 UFI) | ~0.0003% (**down** ~600×) |

Predicted from 8 cells spanning a 20× depth range: 0.278% of reads move from
"kept" to "ribosomal", stable at 0.27–0.51% per cell.

### The run

`START=3 RIBOREF=…v3.fa REFBED=…IntronExonTrna.v2.bed` on
`a_Mapping/submit_vasaplate_map_array.sh`, **unmodified** — overrides via
environment only (Rule 1: fork, never patch). Output `vasaplate_out_v3/`.

Stage chain `51029413` (ribo) → `51029414` (gmap) → `51029415`/`51029416`
(b2bs/b2bm) → `51029417` (cout) → `51029418` (pick). Read from the driver's own
output.

Stages 1–2 are **symlinked, not recomputed** — they are upstream of the rRNA
reference and already concordant at r = 0.974, so re-running them would test
nothing and cost days.

### Three wrapper failures on the way, all mine, none in the pipeline

Recorded because two are new instances of a trap class the conventions already
name — a non-zero exit from a command nobody thinks of as fallible.

1. **Attempt 1** (`51029201`–`51029206`, cancelled after ~1 min): I symlinked
   only the stage-1 `_cbc.fastq.gz`. Stage 3 consumes the stage-2
   `_cbc_trimmed_homoATCG.fq.gz`, so all 384 ribo tasks failed with `fail to open
   file`. **No pipeline product was written** — verified before cancelling. The
   trimmed fastqs live in `vasaplate_out/` and `vasaplate_out_rrnav2/`, not in
   `vasaplate_out_bedv2/`.
2. **Attempt 2** (exit 141): my own readability check `zcat f | head -4` made
   `zcat` die of SIGPIPE when `head` closed the pipe; `pipefail` promoted 141 and
   `set -e` killed the wrapper **before the driver ran**, so nothing was
   submitted. Now reads a bounded byte prefix of the file instead, with
   `pipefail` disabled around it.
3. The first attempt's wrapper also reported exit 1 on an empty `grep` for
   "Submitted batch job" — the driver prints `STAGE N … : <jobid>` instead. The
   driver had actually succeeded.

### Run 5 completed

| stage | job | outcome |
|---|---|---|
| 3 ribo | `51029413` | COMPLETED, 384/384 |
| 4 gmap | `51029414` | COMPLETED, 384/384 |
| 5a/5b | `51029415`/`51029416` | COMPLETED, 384/384 each |
| 6 cout | `51029417` | COMPLETED, **8h48m** (runs 3–4: 8h54m, 8h39m) |
| 7 pick | `51029418` | COMPLETED, **4h52m** (runs 3–4: 4h38m, 4h51m) |
| comparison | `51039884` | COMPLETED, 2m25s |

20 tables written; `mapStats.log` is **21 lines** = complete (trap 2).

> Stage 6 writes its pickle to the **parent** directory with the run-directory
> name as prefix (`vasaplate_out_v3.pickle.gz`, 257 MB; `…dict.pickle`, 1.92 GB),
> not inside the run directory. Looking inside it and finding nothing is not a
> failure — I made that mistake first.

---

## Both predictions: direction CONFIRMED, magnitudes differ and are explained

### Prediction 1 — step-3 depletion UP

| | reads | % of 181,287,059 |
|---|---|---|
| ribosomal, v2 (`rrnav2`) | 9,176,750 | **5.0620%** |
| ribosomal, v3 (fixed) | 9,461,806 | **5.2192%** |
| change | +285,056 | **+0.1572 pt** (predicted +0.28) |

**The v2 baseline reproduces the documented figure exactly** — `5.0620% =
9,176,750 / 181,287,059` against the README's 5.06%. That is itself evidence the
accounting is right.

**Direction is unambiguous: 382 of 384 cells rose, 0 fell, 2 unchanged** (those
two have 1,555 and 3,523 reads).

**Why the magnitude is half.** The shift scales with **rRNA content**, not depth:
Spearman(pct_v2, shift) = **+0.939**, against −0.176 for depth. That is the right
mechanism — v3 recovers reads aligning antisense to one 18S entry, so a cell with
more rRNA has more of them.

The 8-cell pilot that predicted +0.28 sampled **rRNA-rich** cells. The 121 cells
whose shift lands inside the predicted 0.27–0.51 band have median v2-rRNA
**11.16%**, against a plate median of **4.46%**. Taking the 200 cells at or above
that rRNA level gives a read-weighted shift of **+0.354 pt** — reproducing the
pilot. **So the pilot over-estimated by ~1.8× through sampling; +0.157 pt is the
plate-wide truth.**

### Prediction 2 — rRNA biotype share DOWN

| | rRNA UFIs (384 cells) | share of all UFIs |
|---|---|---|
| published | 163 | 0.00030% |
| v2 (run 4) | 95,873 | 0.18417% — **623× published** |
| **v3 (run 5)** | **494** | **0.00095% — 3.2× published** |

**194× fewer rRNA UFIs.** Predicted ~600× down; achieved 194×, because the fix
addresses only one of the two compounding defects.

**It was one locus, and only that locus moved.** From
`rrna_per_locus_v2_v3.tsv` (60 loci):

| | published | v2 | v3 |
|---|---|---|---|
| `ENSMUSG00000106106_CT010467.1_rRNA` | 72 | 95,823 | **444** (216× fewer) |
| all 59 other rRNA loci, summed | 91 | 50 | **50** |

**0 of the 59 other loci changed between v2 and v3.** A one-entry
reverse-complement should move exactly one locus, and it did.

**Why 444 and not 72.** The diagnosis named **two** compounding defects: (A)
`riboread-selection.py` writes reverse-strand reads reverse-complemented, and (B)
one reference entry stored antisense. **v3 fixes B only** — A is upstream code,
left alone under Rule 1. So a residual is expected: A still flips reverse-strand
reads, it just no longer has a backwards reference entry to compound with. 444
UFIs over 384 cells is **1.2 per cell**, against 250 per cell before.

Two ratios, not to be conflated: the **class share** is 3.2× published; the
**locus count** is 6.2× (444/72).

### A third movement, unpredicted and in the right direction

**sncRNA share: 1.604% → 1.408%**, against published **1.421%** and the paper's
stated 1.4%. The error against the deposited table fell from 0.183 pt to
**0.013 pt — 14× closer.** Not predicted, and it follows mechanically: the
mis-assigned reads were inflating a non-protein-coding class.

### Nothing else moved

| | v2 (run 4) | v3 (run 5) |
|---|---|---|
| per-gene Spearman r | 0.9732 | **0.9733** |
| per-gene Pearson r (log10) | 0.9814 | **0.9814** |
| median per-cell Pearson r | 0.9820 | **0.9824** |
| median log2(ours/published) | 0.000 | **0.000** |
| barcodes ≥ 7,500 UFIs | 353 | **353** |
| doublet rate, Fig. 1d rule | 0.57% | **0.57%** |
| shared simple gene rows | 72,760 | 72,753 |

---

## Scrutiny: two things that looked wrong, checked

**`MtRrna` is 0.42× published — a deficit, not an excess.** Checked against run 4:
bedv2 gives **0.43×**, v3 gives **0.44×**. **Pre-existing, not introduced by the
fix.** Explained: the rRNA reference carries **4 mitochondrial rRNA entries**
(12S/16S, human and mouse), so mito-rRNA reads are removed at stage 3 by design
rather than surviving to be annotated `MtRrna`. Whatever the authors used
evidently did not deplete them as thoroughly. Not a workflow fault; worth knowing
before anyone quotes a mito-rRNA fraction from this pipeline.

**`tRNA` is 0.15× published — unchanged from run 4** and already diagnosed as
geometric (`jS:IN` containment vs 72 bp features against 70.8 bp reads), not
annotational. The fix does not touch it and was not expected to.

---

## Verdict: the workflow is validated as a benchmark anchor

Against the **deposited** table, on 72,753 shared genes:
per-gene Spearman **0.9733**, per-cell Pearson **0.9824**, median log2 ratio
**0.000**, barcode count **353 = 353**, sncRNA **1.408% vs 1.421%**.

The one large discrepancy that remained after run 4 is closed: **rRNA from 623×
to 3.2× of published.** The residual is understood and attributable to an
upstream defect deliberately left unpatched.

**Use `vasaplate_out_v3` as the published-plate anchor.** Do not use `bedv2` for
anything rRNA-related.

### Correction made during this analysis

I first wrote "38 other rRNA loci" — that came from the union of a **top-12
listing** in a job log, not from the full table. The correct count is **59**, and
re-deriving it into a committed TSV strengthened the claim rather than weakening
it: 0 of 59 changed, verified over all 60 loci. Every number in the figure now
reads from `rrna_per_locus_v2_v3.tsv`, not from a log.

