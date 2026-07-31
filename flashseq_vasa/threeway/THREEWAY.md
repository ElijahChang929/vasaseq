# Three-way comparison — published VASA-plate / own VASA-plate / FLASH-seq

Status 2026-07-30. The published plate enters as a **validated benchmark anchor**,
not as an assumption: run 5 (`vasaplate_out_v3`) reproduces the deposited table
`GSM5369495` at per-gene Spearman **0.9733**, median per-cell Pearson **0.9824**,
median log2(ours/published) **0.000** over 72,753 shared genes, and 353 = 353
barcodes above 7,500 UFIs. Details in
`vasaplate_check/BENCHMARK_published_plate.md`.

---

## The three datasets, and what is NOT comparable between them

| | published VASA-plate | own VASA-plate | FLASH-seq |
|---|---|---|---|
| id | SRR14783059 | ZHA9292A1 (PM26037) | RN26038 |
| units | 384 barcodes (HEK293T + mESC) | 16 barcodes (12 real, 4 blank) | 10 libraries (input titration) |
| protocol | VASA, 6 nt UMI | VASA, 6 nt UMI | FLASH-seq, no UMI |
| genome | **GRCm38 + hg38** | GRCm39 | GRCm39 |
| annotation | **Ensembl 99, mixed** | Ensembl 116 | Ensembl 116 |
| rRNA reference | `unique_rRNA_human_mouse.v3.fa` (921) | `unique_rRNA_mouse.v2.fa` (356) | same as own plate |
| strand flag | `y` | `y` | `n` (unstranded) |

**Three axes of difference, not one.** Any gap between the published plate and the
other two mixes:

1. **protocol** — the thing we want to measure,
2. **pipeline** — controlled: all three ran the same VASA scripts, and FLASH-seq
   went through them via the `smartseq_noUMI` branch specifically so the
   comparison is method-matched,
3. **annotation release** — Ensembl 99/GRCm38 vs 116/GRCm39, which changes the
   gene universe, biotype labels and gene models.

The published-vs-own comparison is **protocol-free by construction** (same
protocol, same scripts), which makes it the clean place to size the release
effect. That is what the Annotation-control track exists to do, and every other
track cites its number before attributing a gap to protocol.

### Units: reads, on every side

FLASH-seq has no UMI, so `UFICounts` degenerates to a detection mask and
`TranscriptCounts == UFICounts` elementwise on that branch. **Reads are the only
unit all three protocols measure**, so every cross-dataset table here uses
`ReadCounts`. Molecule-based VASA figures are reported alongside where relevant as
VASA's own biology, never as a cross-protocol quantity.

### What is NOT measured, stated up front

**Subunit composition (18S / 28S / 5.8S / 5S) exists for FLASH-seq only.** In
`rrna_comparison.tsv` all eight subunit columns are empty for **both** VASA
plates — 0 of 8 published rows and 0 of 16 own-plate rows carry any value, while
10 of 10 FLASH-seq rows do. Those cells are *not measured*, not zero. Any
three-way statement about where residual rRNA sits (5′ETS vs 18S vs 28S …) would
require computing it for the VASA sides first.

Related: the published-plate rows in that file are **8 cells**, measured against
the mixed human/mouse reference. The plate-wide 384-barcode figure comes from
`depletion_v2_vs_v3.tsv` (v3 columns), not from those 8.

### "Non-rRNA denominator" means non-**cytoplasmic**-rRNA

rRNA is removed by **alignment to the 47S unit**, not by annotation — the Ensembl
`rRNA` biotype is not an rRNA measurement on GRCm39, where the rDNA array is
collapsed out of the primary assembly.

**Mitochondrial rRNA is NOT removed and stays in the denominator on every side.**
`MtRrna` is a distinct biotype from the cytoplasmic 47S species and the
alignment-based filter does not touch it. `threeway_denominators.tsv` declares
its size per dataset: **published 0.0050%, own plate 0.0155%, FLASH-seq 0.216%**.
Largest in FLASH-seq, and still two orders of magnitude below the structural-RNA
differences discussed below — so it changes no conclusion here. But a reader
comparing `MtRrna` rows across these tables should know they were never filtered.

### Filtering: unfiltered tables, deliberately

Composition uses the **unfiltered `uniaggGenes_*.ReadCounts.tsv`**, not the
filtered `analysis/` matrices. The UMI-ceiling filter drops the 8 most abundant
small RNAs — correct for molecule counting, wrong here, where it shifts
ProteinCoding by +19.19 pp.

---

## Results

**Five of six tracks landed. The coverage leg did NOT complete** — the nemo login
node became unreachable mid-run (ssh exit 255, verified twice from the main
session) and stayed down, so it was stopped rather than left holding the session
open on a connectivity outage. It produced **no artifacts and no partial
results**; nothing below is derived from it, and no coverage number is reported
here.

**What the coverage leg was to answer, still open:** the own plate's transcript
coverage is more 3′-weighted (3′/5′ = 1.22) than FLASH-seq's (0.89), and the
question was whether that 3′ rise is a **VASA-protocol property** or specific to
`ZHA9292A1`. Given that both other axes turned out to be library-specific rather
than protocol-level, this is now a more pointed question than when it was
dispatched — if the published plate is also 3′-weighted, the rise is the
protocol; if it is not, it joins the rRNA and structural-RNA findings as another
property of this particular library. **Re-run this track when the cluster
returns.** It needs the published-plate BAMs at
`vasaplate_out_v3/*_E99_Aligned.out.bam`, mouse cells selected via
`threeway_published_cellcalls.tsv`, and transcript models built from the
**Ensembl 99** BED (not the 116 models — those BAMs are against the mixed
GRCm38/hg38 index).

### The headline: two prior conclusions are overturned

Both of the FLASH-seq↔own-plate findings that looked like **protocol** differences
turn out to be **library** differences, and only the published plate could reveal
that — because published-vs-own shares protocol *and* pipeline.

#### 1. rRNA — the own plate is high, but less than first reported

**Species-matched.** The own plate and FLASH-seq are mouse-only, so the published
plate must be restricted to its **mouse** barcodes (Fig. 1d rule) to be a valid
comparator:

| dataset | n | median residual rRNA |
|---|---|---|
| published VASA-plate, **mouse** | 172 | **11.12%** |
| own VASA-plate | 16 | **20.68%** — **1.9× the published plate** |
| FLASH-seq | 10 | **4.72%** |

For reference, the published plate's **human** barcodes sit at 1.80% and the
pooled all-barcode median is 4.67%.

> **Correction.** An earlier version of this document reported the published
> plate at **4.67%** and concluded that "VASA-the-protocol is indistinguishable
> from FLASH-seq" with the own plate 4.4× above it. That 4.67% was the **pooled
> human+mouse** median — the human cells are far cleaner (1.80%) and drag it down.
> On the species-matched comparison the conclusion inverts in part: **mouse VASA
> runs 2.4× FLASH-seq**, so VASA is *not* indistinguishable from FLASH-seq on
> rRNA, and the own plate's excess over the published plate is **1.9×, not 4.4×**.

What survives: **the own plate still carries roughly twice the residual rRNA of a
published VASA library run through the same pipeline.** That is a real, if
smaller, target — and 11.12% is the achievable figure, not 4.67%.

#### 1b. Subunit composition — direction survives, magnitudes do not

`rrna_threeway.tsv` carries a `subunit_comparable` flag, and **all 384
per-barcode published rows are flagged `no` ("human-record competition")**. The
one row the table designates comparable is `AGGREGATE_mouse_realigned`, realigned
to the mouse-only reference:

| | 5′ETS | 18S | ITS1 | 5.8S | ITS2 | 28S |
|---|---|---|---|---|---|---|
| published (realigned, designated comparable) | **71.1%** | 4.1% | 5.9% | 1.6% | 1.6% | **15.5%** |
| own VASA-plate | **18.2%** | 7.0% | 4.8% | 1.8% | 8.0% | **59.7%** |
| FLASH-seq | 47.1% | 9.3% | 4.4% | 2.2% | 6.7% | 29.7% |

> **Correction.** An earlier version quoted the published plate at 5′ETS 47.8% /
> 28S 30.8% and argued that "the two comparators agree closely with each other".
> Those figures are the median of the rows the table flags **non-comparable**. On
> the designated comparable row the two comparators do **not** agree — they differ
> by 24.1 pt on 5′ETS and 14.2 pt on 28S — so the two-independent-references
> argument does not stand.

What survives, and is now larger: **the own plate is displaced toward mature
large-subunit rRNA and away from the 5′ external spacer** — 28S **3.9×** the
published plate (15.5% → 59.7%), 5′ETS **0.26×** (71.1% → 18.2%). Since 5′ETS is
degraded within minutes of transcription while 28S is the stable mature product,
this is the signature of **retaining mature rRNA** rather than of capturing more
nascent pre-rRNA — which points at the depletion step, not at input quality.

#### 2. Structural RNA — the ~90× is not a protocol effect

| dataset | structural RNA, % of non-rRNA |
|---|---|
| published VASA | **2.25%** |
| own VASA | **20.59%** |
| FLASH-seq | **0.26%** |

The two-way work concluded "structural RNA differs ~79× between protocols". But
**the two VASA plates differ by 9.2× from each other**, on the same protocol and
the same scripts. The published-vs-FLASH-seq contrast is **8.7×**, not 79×.

**The structural gap mostly SURVIVES a matched annotation — the shared-universe
control was biased.** The user asked whether their plate could be re-quantified
against the *same old reference* the published plate used. It was: the own plate
was re-mapped to GRCm38 and re-assigned against the Ensembl 99 BED — using
`combined_genome.fa` and `combined.gtf`, **the published index's own build
inputs**, so the gene models are byte-identical, not merely the same release.

| | structural RNA, % of non-rRNA |
|---|---|
| published plate, E99 | **2.25** |
| own plate, E116 | **20.59** |
| own plate, **E99 (matched)** | **17.78** |

| gap estimate | value | release explains |
|---|---|---|
| raw | 18.34 pp | — |
| shared-gene-universe control | 2.89 pp | 84.3% ← **biased, do not use** |
| **fully matched annotation** | **15.54 pp** | **15.3%** |

> **Correction, and it runs the other way from the last one.** This document
> previously reported that **84.3%** of the structural gap was annotation release,
> leaving a 2.89 pp residual. Under a fully matched annotation the residual is
> **15.54 pp** and release explains only **15.3%**. The shared-gene-universe
> control was not merely imprecise, it was **biased**: it retained only **13.10%**
> of the own plate's 8,390,129 structural reads, and the signal sits precisely in
> the multi-gene combination rows and release-specific rows it discarded. It
> reported **18.6%** of the true gap and missed **81.4%**. An earlier version of
> this document also reported **43.8%** from the annotation scale — all three
> numbers were measured on different input spaces, and only the matched-annotation
> one answers the question.

**Per class, after matching the annotation:**

| class | raw gap | matched gap | release |
|---|---|---|---|
| snRNA | +6.87 pp | **+6.97 pp** | −1.5% |
| MiscRna | +8.57 pp | **+5.34 pp** | 37.7% |
| snoRNA | +2.18 pp | **+2.38 pp** | −9.2% |
| ribozyme | +0.68 pp | **+0.80 pp** | −16.8% |
| scaRNA | +0.04 pp | **+0.05 pp** | −10.4% |
| *lncRNA* | +11.47 pp | *+0.88 pp* | **92.4%** |

**snRNA, snoRNA, ribozyme and scaRNA are unchanged or LARGER under the matched
annotation** — these are real, protocol-driven differences, and snRNA is the
largest surviving effect at 44× the published plate. **MiscRna is mixed** (37.7%
release). **lncRNA collapses to +0.88 pp — 92.4% annotation**, and any lncRNA
claim from the earlier three-way reading should be withdrawn.

**Read length is not the explanation:** measured across the five structural
classes it totals **−0.0007 pp**, four orders of magnitude below 15.54 pp. And
96.48% of snoRNA exon features are shorter than a single 151 nt read, so longer
reads are structurally *worse* on short features — the confound **opposes** the
observed difference rather than creating it.

**One asymmetry, quantified.** The species filter differs between the own-plate
E99 and E116 routes (it does not affect the published-vs-own comparison, where
both arms use the same filter and reference), so the matched gap is reported as a
range: **[15.54, 18.42] pp** — both ends **5.4–6.4×** the shared-universe
control's 2.89 pp.

#### 2b. Transcript coverage — the 3′ rise is this library's too

The own plate's 3′ rise does not merely fail to appear in the published plate —
**it goes the other way.**

| group | n | 3′ rise | 3′/5′ | aligned p50 |
|---|---|---|---|---|
| published VASA-plate | 12 cells | **0.657** | 0.798 | 74 nt |
| own VASA-plate | 6 cells | **1.387** | 1.132 | 127 nt |
| FLASH-seq native | 4 libraries | 0.984 | 0.788 | 150 nt |
| FLASH-seq VASA-trimmed | 4 libraries | 0.895 | 0.659 | 107 nt |

The two VASA arms sit on **opposite sides of flat** — published 3′-*depleted*
(all 12 cells below 1.0, range 0.543–0.763), own 3′-*enriched* (all 6 above 1.0,
range 1.344–1.432). The per-unit ranges are **disjoint**, Cohen's d = **10.5**.
The exact permutation p = 5.4 × 10⁻⁵ is the resolution floor at n = 12 vs 6, so
**the separation carries this, not the p-value**.

Three alternative explanations were tested and each rejected:

- **Not read length.** The published library is ~75 nt, 100% of aligned reads
  ≤80 nt. Restricting every arm to ≤80 nt *widens* the gap (own 1.392 → 1.431).
- **Not the annotation release.** The confound is real and unremoved, but the one
  mechanism that could fake it is falsified: Ensembl 99's models are the
  *shorter* ones, so 3′-end truncation predicts published `lost_txend` > own;
  measured, published is **0.212%** vs own **0.297%** — the wrong way round. The
  gap also survives on reads that lose no bases at all (0.652 vs 1.385).
- **Not depth.** Published cells are ~10× shallower; restricting to transcripts
  with ≥20 placed reads gives 0.616 vs 1.398, still disjoint.

> **Magnitude caveat, from the track's own audit of its metric.** `mid` votes at
> the *genomic-span* midpoint, so a junction-spanning read whose midpoint lands in
> an intron casts no vote while still counting as placed — 30.6% of placed reads
> for the published plate, 41.9% own, 47.8% FLASH-seq native, and per-unit dropout
> tracks aligned length at r = 0.982, i.e. it is coupled to the very variable
> `mid` exists to neutralise. **The direction survives** (own is 3′-heavier under
> `base` too, which has no such dropout), but the *magnitude* depends on it. A
> transcript-coordinate midpoint would fix it and is **not applied**. Do not put a
> rise magnitude in a figure legend before that is done.

**Scope.** One plate each, n = 6 own vs 12 published. This distinguishes *this
library* from *the protocol*; it cannot say whether chemistry, input RNA quality
or handling produced the rise. **A second own plate is the discriminating
experiment.**

#### 3. Detection — ordering is stable, published plate lowest

Shared gene universe, mouse cells, medians:

| depth | own VASA | FLASH-seq (trimmed) | published VASA |
|---|---|---|---|
| 100 k | 10,736 | 9,987 | 8,298 |
| 500 k | 13,690 | 12,785 | 10,835 |
| 1 M | 14,818 | 13,769 | — |
| 5 M | 16,814 | 15,646 | — |

The own plate leads FLASH-seq at every matched depth, reproducing the two-way
result. **Verified against `detection_threeway_crossings.tsv`: 0 of the
cross-dataset pairs cross** — the ordering is stable, not an artefact of the depth
chosen.

Two caveats that the crossings table makes explicit and the medians hide:

- **FLASH-seq's own 30 pg rung falls BELOW the published plate.** The published
  plate is lowest of the three *datasets*, but not lowest of every track: the
  one-cell-equivalent FLASH-seq rung is lower still. `VASA published (mouse
  cells)` leads `FLASH-seq 30 pg (trimmed)` throughout.
- **The only crossings anywhere are within FLASH-seq**, between its 30 pg rung
  and its other arms (at 13.4 M and 19.5 M reads on the shared universe). That is
  a within-protocol depth effect, not a between-protocol one.

**The published plate's cell-selection rule matters, and depth-dependently.** The
paper states two rules that disagree — Fig. 1d thresholds the UFI fraction,
Methods p. 18 the gene fraction — and they call different numbers of mESC cells
(**173 vs 141** of 384). On detection the two rules differ by **16.6 genes at the
shallowest rung and 94.5 at the deepest** (uniagg/all; 16.6→94.0 on the shared
universe). A single "at most ~17 genes" figure would be the shallow end only.
This is small against the 2,000–4,000-gene gaps between datasets, so it does not
disturb the ordering — but it is not depth-independent and should not be quoted
as one number. All figures here use the **Fig. 1d** rule.

---

## What this means for the project

### All three axes point the same way

Every axis tested so far separates **this library** from the VASA protocol, not
VASA from FLASH-seq:

| axis | own plate | published VASA plate | verdict |
|---|---|---|---|
| residual rRNA | 20.68% | 11.12% | **1.9× excess, library-specific** |
| structural RNA | 20.59% | 2.25% | 9.2× raw; **15.54 pp survives a matched annotation** |
| 3′ coverage rise | 1.387 | 0.657 | **opposite sides of flat**, disjoint, d = 10.5 |

Three independent measurements, one plate each, all pointing at `ZHA9292A1`
rather than at the protocol. That is the single most useful result of this
comparison — and it means **a second own plate is the discriminating experiment**
for all three at once.

**The user's library carries roughly twice the residual rRNA of a published VASA
library run through the same pipeline.** Species-matched, the published plate's
mouse cells sit at **11.12%** against `ZHA9292A1`'s **20.68%** — a **1.9×**
excess, on the same protocol, the same scripts and the same benchmark-validated
route. **11.12% is the achievable target.**

Note this is *not* a claim that VASA is as clean as FLASH-seq: mouse VASA runs
**2.4×** FLASH-seq's 4.72%, so some of the gap to FLASH-seq is the protocol. The
1.9× is the part that is specific to this library, and it is the part worth
chasing.

**The 28S / 5′ETS signature is the diagnostic lead.** Against the designated
comparable published row, the own plate carries **3.9×** the mature 28S and
**0.26×** the 5′ETS. Since 5′ETS is degraded within minutes of transcription
while 28S is the stable mature product, this is retention of *mature* rRNA rather
than capture of more nascent pre-rRNA — which points at the depletion step, not
at input RNA quality.

**Any structural-RNA claim must state the annotation release — and a partial
control is not enough.** Under a **fully matched** Ensembl 99 annotation, release
explains only **15.3%** of the apparent VASA-vs-VASA structural difference
(18.34 pp raw → **15.54 pp** surviving). The per-class picture: **snRNA, snoRNA,
ribozyme and scaRNA are unchanged or larger** under matching (real,
protocol-driven), **MiscRna is 37.7% release**, and **lncRNA is 92.4% release**
and should be withdrawn as a finding. No single correction factor applies, and
the shared-gene-universe shortcut is **biased** on these classes — it kept 13.10%
of the structural reads and missed 81.4% of the gap.

---

## Carried-forward findings the three-way must respect

From the completed FLASH-seq ↔ own-plate work (`COMPARISON.md`):

- ~~**Structural RNA differs ~90× between protocols**~~ — **SUPERSEDED by this
  document.** The two-way work reported VASA 20.59% vs FLASH-seq 0.23% summed over
  five classes and read it as a protocol difference, because it had no second VASA
  library to compare against. It does not survive the three-way: the published
  VASA plate sits at **2.25%**, so the two VASA plates differ **9.2×** from each
  other, and the published-vs-FLASH-seq contrast is **8.7×** not ~90×. Under a
  fully matched annotation **15.54 pp of the VASA-vs-VASA gap survives** (release
  explains only 15.3%), so the gap is real — it just separates the two libraries
  rather than the two protocols. What *does* survive from the
  two-way work is the read-length control: **read length moves those same classes
  by ≤0.01 pp**, so the gap — whatever its size — is not a read-length artefact.
- **Read length IS a confound for short-RNA detection**: length-matching recovers
  1.7–2.6× more short-RNA entries, so "FLASH-seq detects fewer short species" is
  not attributable to chemistry without the control.
- **VASA detects more genes at matched depth** — 796–1,226 more than the
  VASA-trimmed FLASH-seq arm at every depth from 100 k to 5 M reads.
- **FLASH-seq is not 3′-biased** (3′/5′ = 0.89); the own VASA plate is the more
  3′-weighted (1.22). Reported on the midpoint metric, which nf-core's qualimap
  independently validates; the per-base metric disagrees by 2× for reasons that
  remain **unknown** — an edge-clipping hypothesis was proposed and then
  **withdrawn** when it turned out to rest on a wrong read length.
- Two low-count classes (*ribozyme* 7 entries, *scaRNA* 18) are unstable across
  library-count changes and are order-of-magnitude only.
