# Where the rRNA numbers come from

Written 2026-07-28 in answer to a direct question: *is the 20%+ rRNA figure my
data or the paper's, how was it counted, and how did the paper count it — n or y?*

Short answer: **it is your data, counted with `stranded=y`, which is exactly what
the original VASA code specifies. The paper publishes no rRNA percentage at all.**

---

## 1. The number is yours, not published

| | value | source |
|---|---|---|
| VASA rRNA fraction | **21.39%** = 19,282,729 / 90,137,383 reads | your run |
| sample | `ZHA9292A1` | |
| run | PM26037 | |
| flowcell | `20260720_LH00442_0273_B23TM55LT4` | symlinked from your STP delivery |
| measured by | `pipeline.sh step3` -> `ribo-bwamem.sh` -> `riboread-selection.py` | job `50788552`, 2026-07-26 |

Nothing in this comparison quotes a published VASA figure. The FLASH-seq side is
likewise your own run `RN26038` (`ZHA8833A1..A10`), measured by
`code/flashseq/05_rrna_bwa.sh`, job `50855065`.

## 2. What your data actually shows

21.39% is a whole-library average over 16 barcodes with a real spread:

| | rRNA % (`stranded=y`) |
|---|---|
| 12 real cells | **17.99 – 25.09%** |
| 4 blanks (001/014/015/016) | 5.86 – 10.49% |
| whole library | **21.39%** |

The blanks being ~3x lower is expected and is a useful control: with little input
RNA there is little rRNA to carry through, so the depletion step has less to
remove.

Subunit composition is stable across cells — 28S is 51.7-55.5% of ribosomal
reads, 5'ETS 15.2-17.8% — which is what makes the FLASH-seq comparison
interpretable at all: a difference in the *total* can be attributed to a
difference in a *subunit* rather than to noise.

Also measured, for reference (this figure had never been computed before):

| flag | VASA rRNA % | what it counts |
|---|---|---|
| `y` (used) | **21.39%** | forward-strand ribosomal hits only |
| `n` | 28.10% | ribosomal hits in either orientation |

## 3. How you counted it

`own_version/config.sh` line 104:

```bash
# Stranded protocol? VASA is stranded -> y
STRANDED="${STRANDED:-y}"
```

That flag is passed through `ribo-bwamem.sh` into `riboread-selection.py`, whose
`stranded == 'y'` branch keeps a read as ribosomal only when it is **mapped and
on the forward strand**. A read whose only ribosomal alignment is reverse-strand
is *not* counted as rRNA and goes on to the non-ribosomal arm.

## 4. How the original method counts it — `y`

This is explicit in the code you were given, not an inference from behaviour:

- `a_Mapping/run_mapping_stepwise.sh` line 41: `STRANDED=y  # VASA is stranded -> y`
- `a_Mapping/README.md`: the parameter "should be set to \"y\" when dealing with
  VASA data", and again for the mapper-assignment step, "where `${stranded}` is
  set to \"y\" in VASA data".
- `riboread-selection.py` is **unmodified from as-received** — `git diff` against
  its first commit (`8bee6f6`) shows only a file-mode change.

So your run used the published logic with the published flag. The paper's own
justification for strand-awareness is consistent: it notes that the reads retain
strand specificity, which improves quantification of overlapping transcripts.

## 5. What the paper does and does not report

Searched the full text of *High-throughput total RNA sequencing in single cells
using VASA-seq* (Nat Biotechnol, `10.1038/s41587-022-01361-8`) for every mention
of rRNA and ribosomal, and every percentage within 400 characters of one.

**There is no numeric rRNA percentage in the paper.** Not one sentence
containing "rRNA" or "ribosomal" anywhere in the text also contains a `%` sign
(checked mechanically over the whole PDF: 0 such sentences). The rRNA statements
are qualitative — they concern the workflow (aRNA is depleted of rRNA) and a
comparison in which both VASA-seq workflows outperformed Smart-seq-total. The
percentages near those passages are about **gene detection and unspliced reads,
not rRNA**.

Two consequences:

1. **There is no published rRNA number to agree or disagree with.** Your 21.39%
   cannot be "wrong relative to the paper"; the paper does not make that claim.
2. **Your unspliced fraction can be compared, and this is the right comparison
   to make.** The paper reports unspliced reads as **44.1 +/- 10.1% for
   VASA-plate** and 56.5 +/- 3.1% for VASA-drop, against Smart-seq3 14.8 +/- 2.5%,
   10x Chromium 17.7 +/- 12.8% and Smart-seq-total 38.1 +/- 4.6%. Your format is
   **VASA-plate**, so 44.1 +/- 10.1% is the figure to use.

   Your 12 real cells give **28.69-33.14%, median 30.9%** — about **1.3 SD below**
   the paper's plate mean. Say it that way; do not say it matches. Three reasons
   not to read it as a defect: the paper's figure is HEK293T (a cultured line)
   against your mouse embryo cells; the paper counted reads aligning to introns
   *or* to exon-intron junctions as unspliced, and whether this pipeline's
   `exon-intron` combination labels are treated the same way needs checking
   before the two are put side by side; and your spread (4.5 points across 12
   cells) is far tighter than the paper's +/- 10.1, which is itself informative
   about library consistency.

   **Do not write "matches the paper".** Same order, same direction relative to
   poly-A methods, not a matched measurement.

## 6. Why FLASH-seq needs the other flag

`stranded=y` is right for VASA and wrong for FLASH-seq, and this is measured, not
assumed:

| | forward strand share of ribosomal reads |
|---|---|
| VASA (`ZHA9292A1`) | **76.1%** (57.7-83.8% across the 12 real cells) |
| FLASH-seq (`ZHA8833*`) | **49.1-50.5%** |

FLASH-seq is unstranded — its ribosomal reads are split evenly, as an unstranded
library's should be. Applying `y` there would discard ~50% of its genuine
ribosomal reads and report half the true rRNA content, penalising it for its
chemistry rather than measuring its rRNA.

**Decision (user, 2026-07-28): each protocol keeps its own correct flag** — VASA
`y`, FLASH-seq `n` — giving **4.5x**. That ratio must always be quoted as "each
under its own correct flag", never as a bare division. Sensitivity, if it is ever
asked for: both `n` gives 5.9x, both `y` gives 9.0x.


---

# Addendum — the published VASA-plate compared directly (2026-07-28)

Prompted by a direct challenge: *are you sure the published plate has the same
rRNA proportion as mine? It must be wrong according to your report.*

The challenge was worth making. **I had said the paper publishes no rRNA
percentage — that is still true — but I had not checked whether the published
plate's own DATA is in this repo. It is**, and a prior session had already run
the rRNA stage on it. So the comparison the paper does not provide can be made
directly from reads.

## What was measured

Published VASA-plate, SRA **`SRR14783059`**, 8 cells, at
`data/ref/fastq_vasaplate/rrna_validation/unique_rRNA_human_mouse.v2/`.
Measured by re-counting the two predicates `riboread-selection.py` itself
uses, over its own `.nsorted.all-ribo.bam` files — the same script, the same
`ribo-bwamem.sh`, the same v2 reference family as your run.

| | published plate (8 cells) | your plate (12 real cells) |
|---|---|---|
| per-cell rRNA % (`stranded=y`) | **1.33 – 26.42%** | **17.99 – 25.09%** |
| median | 9.90% | 21.15% |
| pooled | **9.19%** | **21.39%** |
| max/min across cells | **19.8x** | **1.4x** |
| IQR | 12.59 points | **1.73 points** |

## Two conclusions, and they point in opposite directions

**1. Your rRNA fraction is about 2.3x the published plate's pooled figure.**
9.19% vs 21.39%. So the concern behind the challenge is real: on this evidence
your depletion is leaving more rRNA behind than the published plate did. That is
a genuine finding and it should be treated as one.

**2. But your library is far more CONSISTENT than the published one.** The
published plate's 8 cells span 1.33% to 26.42% — a 19.8-fold range, IQR 12.6
points. Yours span 17.99% to 25.09% — 1.4-fold, IQR 1.7 points. Your worst cell
is better than the published plate's worst cell, and your spread is an order of
magnitude tighter.

Both of those matter. A high mean with a tight spread is a *systematic* offset —
something in the depletion step operating the same way in every well, which is
the kind of thing that can be tuned. A low mean with a 20-fold spread is
*variable* depletion. For a method still being established, the tight spread is
arguably the better starting point, but the offset is real and worth chasing.

## Caveats that must travel with these numbers

- **8 cells, not the whole published plate.** These are the cells a prior session
  pulled for rRNA-reference validation, not a random sample of the published
  experiment, and one of them (005) carries 1.26 M reads while another (002) has
  64 k. The pooled 9.19% is therefore dominated by a few deep cells; the median
  (9.90%) and the per-cell range are the more robust summaries.
- **Different reference file.** The published plate was run against
  `unique_rRNA_human_mouse.v2.fa` (the published experiment is a human/mouse mix);
  yours against `unique_rRNA_mouse.v2.fa`. Same builder, same version, different
  species content. `vasaplate_check/README.md` documents a known orientation
  defect in the *mixed* reference that inflates one antisense rRNA locus ~600x in
  the count tables, fixed in a v3 that has been built but not yet used. That
  defect is in gene *assignment*, not in this ribosomal-read count, and
  `own_version/` is explicitly recorded as unaffected — but it is a reason to
  treat the mixed-reference numbers as provisional.
- **Different biology.** The published plate is HEK293T/mouse mix; yours is mouse
  embryo. rRNA content is not a constant of nature.

## A separate finding, and possibly the more interesting one

**Forward-strand share of ribosomal reads**, now measured per cell on both sides
(job `6689a020`, `own_strand_percell.tsv`; the same script reproduces step 3's
`pct_y` for all 16 barcodes to within 0.005 points, which is the self-check):

| | forward share |
|---|---|
| published plate (8 cells) | **92.7 – 97.4%**, median 95.4% |
| your plate, 12 real cells | **57.7 – 83.8%**, median 74.8% |
| **your plate, 4 blanks** | **86.4 – 93.6%** |

The published plate behaves as a strand-specific protocol should — almost all
ribosomal hits forward. **Your real cells do not**, with roughly a quarter of
ribosomal reads antisense.

**But your blanks do.** They sit at 86-94%, with the published plate rather than
with your own real cells. That is the most informative number here, and it
constrains the cause:

- It is **not the reference, the aligner, or the counting predicate** — the blanks
  went through the identical path and come out strand-specific.
- It is **not a plate-wide reagent or handling failure** — that would move the
  blanks too.
- It **scales with input RNA**: the wells with real RNA lose strand specificity,
  the wells with almost none retain it. A mechanism that acts on abundant
  template — second-strand synthesis, IVT run-on, or template switching during
  amplification — fits; a mechanism that acts on the well does not.

Cell 002 is the extreme (57.7% forward, i.e. nearly symmetric) and also has the
lowest rRNA of the real cells at 17.99%. Worth checking whether antisense share
and rRNA fraction anti-correlate across cells — if they do, some of the
"stranded=y" rRNA differences between cells are strand artefacts rather than
depletion differences.

**A consequence for your own pipeline, not just for this comparison:** with 24%
of ribosomal reads antisense, `stranded=y` sends them to the *non*-ribosomal arm,
where they go on to be gene-assigned. `vasaplate_check/README.md` already
documents exactly this failure mode on the published mixed reference — one
antisense rRNA locus inflated ~600x in the count tables. Your mouse-only
reference is recorded as unaffected, but that was checked when the antisense
fraction was assumed small. **It is worth re-checking `own_version/` count tables
for antisense rRNA leakage given 24%, not the assumed few percent.**

This is worth following up independently of the rRNA total, because strand
specificity is a property of the library chemistry (the UFI/strand-specific
tailing step), not of the depletion. It also means `stranded=y` discards more of
your ribosomal signal than it does of the published plate's — which is exactly
why the flag pairing had to be settled explicitly.

**Neither of these two findings changes the FLASH-seq comparison**, which is
between your two libraries under each protocol's own correct flag. They change
how your VASA library should be described relative to the published method.
