# The three-way comparison in the published figures' own terms

Published VASA-plate / own VASA-plate / FLASH-seq, rendered so they can be read
against the paper directly. Extends
`code/flashseq_vasa/paperfig/paperfig_compare.py` (FLASH-seq vs own plate) by
adding the published plate as a third group.

Every number below is asserted against the source TSV by
`verify_threeway_paperform.py` — **89 checks, 0 failures**. Nothing here is
transcribed.

---

## What is the authors' code and what is my reimplementation

This is the first question to answer about any of these panels, so it is answered
first and it is answered per line, not per figure.

| element | provenance |
|---|---|
| `paper_ubiotype()` — Fig 2b biotype grouping | **AUTHORS'.** The `add_metadata()` one-liners from `b_Analysis/02_scanpy_QCxBiotype.py:109-116`, character for character |
| `paper_id()` | **AUTHORS'.** `add_metadata()` line 107 |
| class list ProteinCoding / lncRNA / smallRNA / tRNA | **AUTHORS'.** `biotype_split()` line 132, verbatim |
| smallRNA membership = snRNA, snoRNA, MiscRna, scaRNA, ribozyme, miRNA | **AUTHORS'.** Line 132; their choice, not mine |
| TF / cofactor test | **AUTHORS'.** `add_metadata()` line 121-122, keyed on their `Mus_musculus_TF.txt` (**1,623** symbols) and `Mus_musculus_TF_cofactors.txt` (**970**) |
| Fig 2b denominator = full per-unit matrix sum, rRNA **not** excluded | **AUTHORS'.** `add_metadata()`'s `n_counts`, line 112 |
| cell-calling rules (Fig 1d UFI-fraction; Methods gene-fraction; 7,500-UFI gate) | **AUTHORS' RULES**, as stated on p.1781 and p.18. The *implementation* is this project's `vasaplate_check/vp_common.py` |
| Fig 1f saturation estimator | **MINE.** The code for this panel was never deposited (grepped every `.py`/`.R` in the repo). Axes, depth grid (5/10/15/20/25/50/75k) and units are the paper's; the estimator is a binomial-thinning expectation, `E[genes] = Σ_g [1 − (1−p)^{c_g}]`, exact in expectation and deterministic |
| mirrored-density rendering of Fig 2b | **MINE**, as a reading of the published panel's form |
| bin counts, axis limits, panel 3 of each figure | **MINE**, and each is a deliberate departure — see below |

`b_Analysis/` and `a_Mapping/` were asserted unmodified **before and after** every
job, inside the job script (`git status --porcelain`, non-empty ⇒ exit 1). Both
jobs exited 0. Nothing in this work imports from, writes to, or edits them.

## Two provenance defects in the two-way version, found by reading it

Neither is silently corrected. Both variants are computed and both are in the TSV,
per Rule 3.

**Defect 1 — the two-way figure compared two different table families.** It
loaded the own plate's `uniaggGenes` table (222,421 rows) against FLASH-seq's
pre-aggregation `_total` table (270,217 rows); its own log records
`VASA 222420 entries` vs `FLASH-seq native 270217 entries`. Gene aggregation is
exactly what collapses multimapper combination entries, so mixing families biases
the combination-entry rate — the statistic the Fig 1f scope flip rests on.
Measured:

| group | combination entries, `uniagg` | pre-aggregation | 
|---|---|---|
| own VASA-plate | 76.43% | 82.87% |
| FLASH-seq native | 79.16% | 83.67% |
| published VASA-plate | 85.13% | 86.32% |

All three groups here use `uniaggGenes`. `table_family` is a column in the TSV and
the pre-aggregation variant is computed alongside.

**Defect 2 — the TF class silently included combination entries.** The authors'
tables are `shortGeneNames_*`, so their `id` field holds gene **symbols**
(verified: their row labels read `0610005C13Rik_lncRNA`). For a combination entry
their `id` is `Sym1-Sym2`, which cannot match a single symbol in the TF list —
**so the authors' TF class contains single-gene entries only.** The fork's
`paper_symbol()` takes `parts[1]` of the full `ENSMUSG_Sym_Biotype` label, which
for a combination returns the *first* gene's symbol and so TF-assigns entries the
authors' code scored `-`. Both are reported: `TF`/`Cofactor` keep the fork's
behaviour so the two-way numbers stay reproducible, `TF_authors`/`Cofactor_authors`
restore the authors' semantics. The fork inflates the TF share by 1.02–1.08×. The
figure plots `TF_authors`.

---

## 1. Fig 2b, three groups — `paperfig2b_threeway_mirrored.png`

Mirror axis is protocol, as in the paper: VASA left, FLASH-seq right, second arm
on each side as a dashed outline. Medians of the plotted (included) units,
**ReadCounts on all sides** (Rule 4 — reads are the only unit both protocols
measure):

| panel (authors' class) | published, mESC, n=173 | own plate, n=12 | FLASH-seq native, n=9 |
|---|---|---|---|
| Protein coding | **0.9278** | **0.6347** | 0.8918 |
| lncRNA | 0.0242 | 0.0577 | 0.0429 |
| Transcription factors | 0.0687 | 0.0388 | 0.0408 |
| sncRNA | 0.0110 | **0.1869** | 0.0026 |

**tRNA is 0.0000 in all three groups, with 0 tRNA entries in every
`uniaggGenes` table.** That is structural, not chemistry: the aggregation stage
does not carry tRNA rows forward, and separately (trap 8) `jS:IN` containment
against ~72 bp tRNA features cannot retain a longer read. It is not evidence about
any protocol's tRNA capture.

**Bin counts deliberately differ between groups, and this is the departure the
larger n makes possible.** The published plate gives n=173 for the first time, so
it gets Freedman-Diaconis per panel (**30 / 20 / 22 / 30** bins for protein
coding / lncRNA / TF / sncRNA); the own plate (n=12) and FLASH-seq (n=9) keep the
two-way figure's fixed 7. FD is *not* used at small n, and that was measured
rather than assumed: on these same panels FD returns 4/2/2/2 bins for the own
plate but 30/4/15/10 for FLASH-seq native, because at n≈10 the IQR collapses
relative to the range. A rule that returns 30 bins for 10 points is not a rule.
Every unit is drawn as a dot regardless, so n is never hidden by the density.

Densities are normalised to their own maxima — with n of 173 against 9 a shared
count scale would make the small groups invisible.

**ZHA8833A8 is excluded from every density and median** (`qc_verdict=exclude`,
18.3% human CALB1) but still drawn, as an open ring: its contamination is a
finding, and the verdict filters interpretation, not QC. FLASH-seq is therefore
n=9 here.

### What the panel shows, and what it cannot attribute

The own plate carries **0.635** protein-coding read fraction against the published
plate's **0.928** — a 0.293 gap — and **17.1×** the sncRNA share (0.1869 vs
0.0110). Those are large. **They cannot be attributed to protocol**, because the
published plate is Ensembl 99 on GRCm38 and mESC, while the own plate is Ensembl
116 on GRCm39 and mouse embryo. Three causes are superimposed: protocol, pipeline,
and annotation release, plus a genuine biological difference in cell type.

What *can* be said: both are the same protocol (`vasa`, 6 nt UMI) through the same
pipeline, so the gap is **not** a protocol difference — which makes it a caution
about reading the own plate's composition against the paper's numbers, not a
result about VASA-seq.

A partial measurement of the release effect (in
`paperform_threeway_release_probe.tsv`): the E99 biotype vocabulary has 32 tokens
against E116's 26, and **all 6 differing tokens are E99-only**
(`IgPseudogene`, `PolymorphicPseudogene`, `RrnaPseudogene`, `TrJGene`,
`TranslatedProcessedPseudogene`, `scRNA`) — E116 introduces none the E99 run
lacks. None of the six is in the authors' four classes, so the release does not
move class *membership*; it moves the gene models and the denominator. Bounding
that properly needs the same cells quantified on both releases, which this
comparison does not have.

**Denominators (Rule 5).** Median assigned reads per unit: published plate
191,455; own plate 3,673,317; FLASH-seq native 23,072,058. Median trimmed reads:
276,377 / 6,387,200 / 26,686,349. The published plate's units are ~19× shallower
than the own plate's and ~120× shallower than a FLASH-seq library.

**The four authors' classes do not cover all reads** — they sum to a median 0.967
(published), 0.891 (own), 0.937 (FLASH-seq native), the remainder being the
authors' own `mixed` and `other` (cross-biotype combinations, pseudogenes,
MtTrna/MtRrna/rRNA, Ig/Tr segments). Writing the self-check as `== 1` is what
surfaced this; it is now asserted as `<= 1` with the residual quantified.

---

## 2. Fig 1f, three groups — `paperfig1f_threeway.png`

Both scopes are kept, because the answer depends on whether multimapper
combination entries count as genes and **the published panel does not say which it
used.** Mean entries detected at 75k trimmed reads, mouse cells:

| scope | published, n=173 | own plate, n=12 | FLASH-seq native pooled, n=9 |
|---|---|---|---|
| all entries | 8,299 | 9,583 | 11,415 |
| single-gene only | 7,699 | 8,785 | 9,558 |

The driver is the combination-entry share **among detected entries**, and the
three-way version corrects the two-way figure's numbers for it (which were
computed across mismatched table families):

| group | median % of detected entries that are combinations |
|---|---|
| published VASA-plate | **10.08%** |
| own VASA-plate | **38.93%** |
| FLASH-seq native | **57.73%** |
| FLASH-seq VASA-trimmed | 68.52% |

### The scope flip is real, but the two-way version located it wrongly

Pooling the 9 qc-ok FLASH-seq libraries averages over a **1,000× input titration**
(30 ng down to 30 pg). Detection rises steeply with input — 11,953 entries at
30 ng against 9,808 at 30 pg — so the pooled mean is dominated by the ng-scale
libraries and is **not a single-cell measurement**. Pooled, FLASH-seq leads on
*both* scopes (+1,832 all entries, +772 single-gene), so the flip does not appear
at all.

Restricting FLASH-seq to the **30 pg rung** — the single-cell-equivalent input per
`sample_metadata.tsv` — recovers it:

| scope | own VASA-plate, n=12 | FLASH-seq 30 pg, n=2 | margin |
|---|---|---|---|
| all entries | 9,583 | 9,808 | **+225 FLASH-seq** |
| single-gene only | 8,785 | 8,056 | **−729 VASA** |

So the qualitative claim survives, at a much smaller magnitude than the pooled
comparison implied, and **at n=2 on the FLASH-seq side it is underpowered** — two
libraries, one input rung, no statistical test attempted. Report it as a direction,
not an effect size.

---

## 3. The published-curve check — `paperfig1f_publishedcheck.png`

**Does our re-run of THEIR data reproduce their published curve?** Only answerable
now the published plate is in the comparison.

The paper's value is **HEK293T — human, not mouse.** p.1781: the HEK293T datasets
were downsampled to assess sensitivity "for all annotated genes", giving
9,480 ± 1,252 genes per cell for VASA-plate at 75,000 trimmed reads; the Fig 1
caption adds that only cells sequenced to at least 75,000 reads were used, n=174
for VASA-plate. Running this on mESC would not be a test of the published number,
so the check uses **HEK293T cells (Fig 1d rule) and human entries**. Our cohort is
n=168 against their n=174.

| depth | scope | ours | published | ratio |
|---|---|---|---|---|
| 75k | all entries | 10,577 ± 729 | 9,480 ± 1,252 | 1.116 |
| 75k | **single-gene** | **9,681 ± 663** | **9,480 ± 1,252** | **1.021** |
| 750k | all entries | 19,379 ± 922 | 15,248 ± 1,092 | 1.271 |
| 750k | **single-gene** | **15,582 ± 719** | **15,248 ± 1,092** | **1.022** |

**On the single-gene scope the reimplementation reproduces the published curve to
1.02× at both depths** — 2.1% at 75k and 2.2% at 750k, each well inside the
published s.d. Two independent depths agreeing to the same 2% is a reproduction,
not a coincidence, and the second depth (750k, Ext. Data Fig 2e) was not used to
tune anything.

**This also answers the scope question the published panel left open.** Counting
combination entries as genes inflates the estimate by 1.093× at 75k and 1.244× at
750k, which would put the 750k figure 27% above the printed value. The published
numbers are therefore consistent with **single-gene counting** and inconsistent
with counting combination entries — so where Fig 1f is concerned, the
single-gene panel is the one to read against the paper.

Caveats that remain, stated rather than buried: the estimator is mine, so a
residual 2% could be the thinning expectation rather than the data; the cohort
differs by 6 cells; and this tests the *published plate through our pipeline*, not
the own plate, so it validates the workflow and the reimplementation, not the
own-plate chemistry.

---

## Files

| file | what |
|---|---|
| `res/threeway/paperform_threeway.tsv` | Fig 2b per-unit class fractions, all groups, both table families, both TF rules |
| `res/threeway/paperform_threeway_fig1f.tsv` | Fig 1f saturation, both scopes, all cohorts incl. the HEK293T check |
| `res/threeway/paperform_threeway_units.tsv` | per-unit reads, detection, combination share, species call |
| `res/threeway/paperform_threeway_release_probe.tsv` | biotype vocabulary and entry counts per group — the E99/E116 measurement |
| `res/threeway/threeway_published_cellcalls.tsv` | 384 barcodes under both of the paper's rules |
| `res/threeway/paperfig2b_threeway_mirrored.png` | Fig 2b, three groups |
| `res/threeway/paperfig1f_threeway.png` | Fig 1f, three groups, both scopes, input-matched panel |
| `res/threeway/paperfig1f_publishedcheck.png` | the published-curve check |
| `code/flashseq_vasa/threeway/threeway_paperform.py` | the analysis; provenance marked per function |
| `code/flashseq_vasa/threeway/threeway_figs.py` | the figures |
| `code/flashseq_vasa/threeway/verify_threeway_paperform.py` | 89 assertions over the TSVs |

## Reproducing

```bash
W=/nemo/lab/turnerj/working/guangxin/vasaseq
PY=/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/python   # pandas 2.3.3, numpy 2.2.6
cd $W
$PY code/flashseq_vasa/threeway/threeway_paperform.py --stage precheck --out res/threeway
$PY code/flashseq_vasa/threeway/threeway_paperform.py --stage main --out res/threeway --also-raw
$PY code/flashseq_vasa/threeway/threeway_figs.py res/threeway
$PY code/flashseq_vasa/threeway/verify_threeway_paperform.py res/threeway
```

Precheck job `a825012a` then `eb305bd1` (exit 0); main pass `4a35ee40` (exit 0),
4 cores on `ncpu`. Job ids read from the jobs' own output. `envs/vasa` must stay
on pandas 2.x.
