# Three-way transcript-body coverage: the own plate's 3' rise is this library's, not VASA's

Every number below is re-derived and asserted by
`verify_coverage_threeway.py` (59 checks, 0 failures →
`res/threeway/verify_coverage_threeway.txt`). Nothing here is transcribed from a
log.

## The question

The own VASA plate shows a sharp coverage rise over the last ~10% of the
transcript (3'/5' = 1.22 in the two-way comparison). The other two axes of this
three-way work both turned out **library-specific rather than protocol-level**:
the own plate's rRNA excess is 1.9x the published plate on the same protocol and
the same pipeline, and 84.3% of its structural-RNA excess is annotation release.
So: is the 3' rise a VASA-protocol property, or a third property of this
particular library?

The discriminating experiment is the **published VASA plate** (SRR14783059 /
GSM5369495) — same protocol, different hands, different cells, different
pipeline run.

## The answer

**The rise is this library's.** The published VASA plate does not merely fail to
show it — it goes the *other* way.

| group | n | 3' rise | 3'/5' | aligned p50 |
|---|---|---|---|---|
| VASA-seq, published plate | 12 cells | **0.649** | 0.790 | 74 nt |
| VASA-seq, own plate | 6 cells | **1.392** | 1.188 | 127 nt |
| FLASH-seq, native | 4 libraries | 0.994 | 0.882 | 150 nt |
| FLASH-seq, VASA-trimmed | 4 libraries | 0.911 | 0.734 | 107 nt |

`rise` = mean(bins 90–99) / mean(bins 40–59) of the midpoint profile; 1.0 = flat.

The two VASA arms fall on **opposite sides of flat**: published is 3'-*depleted*
(0.649, every one of 12 cells below 1.0), own is 3'-*enriched* (1.392, every one
of 6 cells above 1.0). The per-unit ranges are **disjoint** — published max
0.7626, own min 1.3442 — so no published cell approaches any own cell. Gap
0.7424 against a pooled within-group sd of 0.0706, i.e. **Cohen's d = 10.5**.

The exact permutation test gives p = 5.39e-05, which is **1 of 18,564 splits**.
At n=12 vs n=6 that is the *smallest attainable* p: the test is at its
resolution floor, not reporting a small p-value. The separation is what carries
the result, not the p.

The own plate's rise peaks in bin 98 at 1.95x its own body mean. FLASH-seq
native sits at 0.994 — flat — which is worth noting on its own: the own VASA
plate is the only one of the four arms that is 3'-weighted at all.

## Why this is not read length

Read length is the first objection: the published library is a ~75 nt library
(STAR index `sjdbOverhang 73`) and 100.0% of its aligned reads are ≤80 nt,
against the own plate's p50 of 127 nt.

Restricting **every** arm to reads ≤80 nt — the length regime all four share —
does not close the gap. It **widens** it: the own plate's rise goes 1.392 →
1.431 while published stays at 0.649, and the ranges remain disjoint. Longer
reads were suppressing the own plate's rise slightly, not creating it.

## Why this is not the annotation release

The confound is real and is **not removed**: the published arm is quantified on
Ensembl 99/GRCm38 models and the other two on Ensembl 116/GRCm39. On the 4,000
shared genes the models differ substantially — median longest-transcript length
ratio E116/E99 = 1.1129, only 13.8% identical, 37.6% within ±5%, 4 strand
disagreements. Gene order is asserted identical across the two model files, so
only the model differs, but the model difference is not small.

Two things bound it:

1. **The mechanism that could fake this is falsified.** E99's models are the
   *shorter* ones. If that truncates the published arm's annotated 3' end, reads
   at the true transcript end would fall beyond the model and be lost rather
   than binned — artefactually depleting published's last bins. That predicts
   published `lost_txend` > own. Measured: published **0.212%** of aligned
   bases, own **0.297%**. Published is *lower*, so the prediction is falsified
   in the direction that matters.
2. **The gap survives on reads that lose nothing.** Restricted to reads every
   base of which landed inside the model (`mid_fullexon`, 93.9% of published and
   92.8% of own placed reads), published 0.6517 vs own 1.3848, still disjoint.
   No model-edge effect of any kind can reach these reads.

The claim is therefore stated on the **direction** of the rise — published below
flat, own above — which no rescaling of a transcript axis can invert. It is not
stated on the size of the difference, which the release does affect.

## Depth is not the explanation either

The published cells are ~10x shallower (2.46M placed reads over 12 cells vs
14.4M over 6). Recomputing over transcripts with ≥20 placed reads: published
0.616 (n_tx 1,957), own 1.398 (n_tx 3,894), still disjoint. The published arm's
shape is not an artefact of thinly covered transcripts.

## Aligned read length, all three on the record

Pooled over each group's units, from all primary alignments streamed:

| group | p05 | p25 | p50 | p75 | mean | ≤80 nt |
|---|---|---|---|---|---|---|
| VASA-seq, published plate | 46 | 71 | 74 | 75 | 69.8 | 100.0% |
| VASA-seq, own plate | 26 | 100 | 127 | 129 | 109.7 | 15.3% |
| FLASH-seq, native | 51 | 107 | 150 | 151 | 127.1 | 15.3% |
| FLASH-seq, VASA-trimmed | 28 | 71 | 107 | 129 | 96.7 | 30.8% |

→ `res/threeway/coverage_threeway_alnlen.tsv`.

This matters because read length is the input to the **untested edge-clipping
hypothesis** for the base-vs-mid disagreement. The published plate's ~74 nt is a
genuinely different regime from the own plate's ~127 nt and FLASH-seq's ~150 nt,
and it is now measured rather than assumed.

## The base-vs-mid disagreement: still unexplained, and now more interesting

`mid` remains the reported metric and `base` diagnostic only — that decision is
inherited (nf-core's qualimap on the same FLASH-seq BAMs gives 3'/5' = 0.93,
matching `mid` 0.882, not `base` 0.455) and is not re-litigated.

The cause of the 2x disagreement is **still unknown**. Edge clipping remains the
**leading untested hypothesis**; the earlier refutation rested on a wrong VASA
read length and is withdrawn, and is not re-derived here.

What this run adds is one new observation, and it argues for caring:

| group | mid 3'/5' | base 3'/5' |
|---|---|---|
| VASA published | 0.790 | 0.537 |
| VASA own | 1.188 | 0.560 |

Under `base`, the two VASA plates are nearly **equal** (0.537 vs 0.560). The
metric that disagrees with the external tool is also the metric that **erases
the library difference**. That is not a proof that `base` is wrong, but it is a
second, independent reason to distrust it.

The loss decomposition is arithmetically exact (`bases_binned` + every loss
class == `bases_aligned`, exactly, on all 26 units; `lost_modelerr` = 0). Reads
overhanging the annotated 3' end carry 0.370% of aligned bases pooled over all
units — too small to be the whole base-vs-mid gap. That **bounds** the loss; it
does not test the hypothesis.

## A large limitation of the reported metric, stated plainly

`mid` votes at the read's **genomic-span** midpoint. For a junction-spanning
read that point can land in an intron, where no model exon is found and the read
casts **no vote at all** — while still counting in `reads_placed`, because its
bases did place. Measured (`coverage_threeway_middrop.tsv`):

| group | reads placed | cast no midpoint vote |
|---|---|---|
| VASA-seq, published plate | 2,463,966 | **30.554%** |
| VASA-seq, own plate | 14,402,031 | **41.896%** |
| FLASH-seq, native | 12,000,000 | 47.798% |
| FLASH-seq, VASA-trimmed | 12,000,000 | 40.591% |

Per-unit dropout tracks p50 aligned length at **r = 0.982**. So `mid` is not
read-length-neutral in the way the word "midpoint" implies — it is only
*vote-count*-neutral, and the dropout is coupled to exactly the variable the
metric was introduced to neutralise. A transcript-coordinate midpoint (from the
in-transcript positions `profile` already walks) would remove this. **That fix is
not applied here**, so every `mid` number above carries this dropout.

What it can and cannot do: a dropped read contributes nothing to *any* bin, so it
can only reshape a profile if dropped reads are non-uniform along the transcript.
`base` has no midpoint dropout, and under `base` the own plate is still
3'-heavier than the published plate — own/published rise ratio 1.850 by `base`,
2.144 by `mid`. **The direction of the headline claim does not rest on the
dropout**; its magnitude does. This should be fixed before the magnitude is
quoted anywhere.

## Scope and honest limits

- **n=6 own vs n=12 published cells, one plate each.** This distinguishes *this
  library* from *the protocol*; it cannot say which of chemistry, input RNA
  quality, or handling produced the rise. A second own plate would be the next
  discriminating experiment.
- **The published cell set**: 384 barcodes; mouse n=173 by the paper's Fig.1d
  rule (UFI fraction), n=141 by its Methods rule (gene fraction). The two rules
  disagree on 32 barcodes. **The 12 cells used here are mouse under BOTH rules**
  and are the 12 deepest such by mouse UFI. Cells were *not* picked by BAM file
  size — the largest published BAMs are HEK293T (human) cells.
- **Depth caps**: the own-plate and FLASH-seq units hit the 3M-placed-read cap;
  the published cells did not (they are one library across 384 barcodes), so
  `reads_placed` must be read alongside any published number.
- The annotation-release difference is **bounded, not removed** (see above).
- Rule 6: profiles here are unstranded. Sense-only midpoint ratios are reported
  in `coverage_threeway.tsv` (`mid_sense`) and shift all four arms in the same
  direction (published 0.790 → 0.771, own 1.188 → 1.072), so the conclusion does
  not turn on the strand flag.

## Where the numbers live

| file | what |
|---|---|
| `coverage_threeway.tsv` | per-unit, per-metric: profiles' summary stats, depth, loss decomposition, aligned lengths |
| `coverage_threeway_profile.tsv` | the 100-bin profiles, per group and metric |
| `coverage_threeway_alnlen.tsv` | pooled aligned read-length distribution per group |
| `coverage_threeway_robust.tsv` | profiles recomputed over transcripts with ≥20 placed reads |
| `coverage_threeway_pubcells.tsv` | all 384 published barcodes, both rules, which 12 were selected |
| `coverage_threeway_geneset.genes.tsv` | the 4,000 genes with both releases' transcript models |
| `coverage_threeway_release_bound.tsv` | E99-vs-E116 model-agreement summary |
| `coverage_threeway_middrop.tsv` | per-unit midpoint-vote dropout (the limitation above) |
| `coverage_threeway.png` | the figure |
| `verify_coverage_threeway.txt` | the 59 assertions |
