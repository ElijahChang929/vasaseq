# Three-way biotype composition on the non-rRNA denominator

**Question this track exists to answer.** The FLASH-seq-vs-own-plate work found
structural RNA classes differing ~90x between protocols while read length moved
them <=0.01 pp. Does the published VASA-plate sit **with** the own plate — making
that a protocol difference — or somewhere else?

**Answer: somewhere else.** The published plate is 8.67x FLASH-seq and the own
plate is 79.44x, so the published plate sits **between** the two, 9.16x below the
own plate. The ~90x is therefore not a clean protocol signature: two libraries
of the *same* protocol, through the *same* pipeline, differ by 9.16x on the same
quantity. About half of that residual (49.1% of the sncRNA gap) is the Ensembl
99-vs-116 confound; a **5.5-5.9x residual survives every control I could apply**
and is not explained.

Every number below is asserted against its source TSV by
`code/flashseq_vasa/threeway/verify_threeway.py` — **72/72 checks pass**. Nothing
is transcribed.

---

## How the table is built (identical on all sides)

`res/threeway/composition_threeway.tsv`, from
`code/flashseq_vasa/threeway/threeway_composition.py`. Construction copied from
`res/flashseq_vasa/mk_vasa_composition.py` so the published side is comparable
with the committed two-way table:

- **Unfiltered `uniaggGenes_total.ReadCounts.tsv`**, never `analysis/`. The
  UMI-ceiling filter drops the 8 most abundant small RNAs (+19.19 pp on
  ProteinCoding) — correct for molecule counting, wrong here.
- **`ReadCounts` on all four sides.** Reads are the only unit VASA
  (`protocol='vasa'`) and FLASH-seq (`smartseq_noUMI`) both measure (Rule 4).
- **rRNA out of the denominator**, measured by alignment at stage 3, never by the
  GRCm39 `rRNA` biotype (trap 9).
- Published plate: **mouse entries only, mouse cells only** (Fig. 1d rule,
  173 cells). The species filter removes 1,700,818 reads = 3.887% of that
  plate's rows.

**Verified reproduction:** this script's own-plate ProteinCoding (64.1192%) and
snRNA (7.0327%) match the committed `composition_flashseq_vs_vasa.tsv` to 4 dp.

**One correction to the two-way table.** Its FLASH-seq side was built from
`_total`, not `uniaggGenes_total` — identified by matching its published
ProteinCoding of 84.319941% rather than by reading the code. On `uniaggGenes`
FLASH-seq ProteinCoding is 89.6426%, **+5.32 pp**. This track uses
`uniaggGenes_total` everywhere, so no cross-side table mismatch remains.

### Denominators (Rule 5)

| dataset | protocol | genome / annotation | units | non-rRNA reads | reads in combination keys |
|---|---|---|---|---|---|
| published VASA-plate | vasa | GRCm38 / Ensembl 99 | 173 mouse cells | 42,050,713 | 3.70% |
| own VASA-plate | vasa | GRCm39 / Ensembl 116 | 12 real cells | 52,823,542 | 20.57% |
| FLASH-seq native | smartseq_noUMI | GRCm39 / Ensembl 116 | 10 libraries | 238,184,504 | 8.32% |
| FLASH-seq VASA-length | smartseq_noUMI | GRCm39 / Ensembl 116 | 10 libraries | 237,986,306 | 10.39% |

---

## The headline result

Structural RNA = MiscRna + snRNA + snoRNA + scaRNA + ribozyme, % of non-rRNA reads:

| dataset | structural | vs FLASH-seq native |
|---|---|---|
| FLASH-seq native | 0.2592% | 1.00x |
| FLASH-seq VASA-length | 0.2585% | 0.997x |
| **published VASA-plate** | **2.2474%** | **8.67x** |
| **own VASA-plate** | **20.5906%** | **79.44x** |

Read length still does almost nothing: the two FLASH-seq arms differ by
**-0.0007 pp**, reproducing the two-way finding.

**The three datasets form three separate bands, with no overlap between the two
VASA plates.** Per cell: the own plate spans 13.83-25.21% (12 cells), the
published plate 0.076-11.63% (173 cells). **0 of 173 published cells reach the
own plate's minimum.** The 9.16x is a property of the plates, not of pooling.

### What does *not* differ

ProteinCoding is 92.22% (published) / 89.64% (FLASH-seq native) / 64.12% (own),
and lncRNA 2.44 / 6.71 / 13.91%. The published plate is closest to **FLASH-seq**
on protein-coding — so the separation is specific to structural classes, not a
global compositional shift.

---

## Which gaps survive, and which do not

Five candidate artefacts were tested, all expressed as the own÷published fold on
structural RNA so they are directly comparable:

| control applied | own ÷ published | verdict |
|---|---|---|
| raw alignments, before any counting decision | **9.38x** | gap already present upstream |
| none (all reads, all rows) | 9.16x | — |
| combination keys re-allocated (unanimous, or split equally) | 9.16x | **no effect** |
| molecules instead of reads (removes PCR duplication) | 8.60x | removes 1.54 of 18.34 pp (8.4%) |
| simple rows, genes present in both releases | **5.94x** | removes ~half |
| + published side relabelled to its E116 biotype | **5.94x** | **no further effect** |

**Allocation is not the explanation.** The biotype rule both upstream scripts use
credits the *last* member of a multi-gene key, and the own plate has 20.57% of
reads in such keys — a real worry. But all three allocation rules give
*identical* structural shares (max spread 0.0), because **0 of 81,333
biotype-discordant own-plate keys and 0 of 16,627 published keys contain a
structural biotype at all.** Structural entries are annotated cleanly or with
same-biotype partners.

**Duplication is not the explanation.** The own plate carries 18x more reads per
unit (4.40M vs 0.24M). Deduplicating both plates to `TranscriptCounts` moves the
gap from 18.34 to 16.80 pp — **8.4%**, leaving 8.60x.

**Geometry is not the explanation, and runs the wrong way.** Trap 8 says read
length gates short features via the `jS:IN` containment requirement. Measured
from the per-cell annotated BEDs using **stage 5's own `jS:IN` tag** (not a
reimplementation), with a deterministic stride: own-plate median aligned span is
**130 nt vs the published plate's 74 nt** — 1.76x *longer*, therefore *harder* to
contain in a short feature. Realised jS:IN on structural features is essentially
equal (**82.97% own vs 82.16% published**). Geometry opposes the observed gap.

**The gap exists before any counting decision.** Structural features' share of
annotated BED rows is **4.835% (own) vs 0.515% (published) = 9.38x** — measured
on raw alignments, with no denominator, allocation rule, annotation release or
filter involved. Only which features the reads landed in.

### The release confound, citing the Annotation control track

The Annotation control track reports (`res/threeway/annotation_crossplate_gap.tsv`):

- sncRNA, all genes: published **0.9666%** vs own **6.3704%**, gap **+5.4037 pp**
  (**6.59x**).
- restricted to genes present in both releases: published **0.6095%** vs own
  **3.3611%**, gap **+2.7515 pp** (**5.51x**).
- **the restriction closes 2.6522 pp = 49.08% of the raw gap.**

My independent restriction on the structural-5 set agrees: 9.16x → **5.94x**.
**Two tracks, different class definitions, same conclusion — about half the
published-vs-own gap is the Ensembl 99-vs-116 confound, and a 5.5-5.9x residual
is not.**

The mechanism is **gene-set membership, not relabelling**. Relabelling the
published side to its E116 biotypes moves structural composition by
**0.000000 pp**, and the churn ledger shows all five structural biotypes have
**0 shared genes lost to or gained from another class** between releases. The
release effect is instead asymmetric in *which genes exist*: 26,583 mouse gene
ids are E116-only against 3,706 E99-only, and by the annotation track's own
figures the own plate loses **3.245%** of its sncRNA reads to release-only genes
against the published plate's **0.363%**.

---

## What this means, stated at the confidence it earns

1. **The two-way "~90x protocol difference" does not survive as stated.** Two
   VASA-plate libraries of the same protocol through the same pipeline differ
   9.16x. The protocol-vs-FLASH-seq contrast is real and large (8.67x on the
   published plate, the conservative anchor), but 79.44x is specific to the own
   plate, not a VASA-seq property.
2. **~Half the published-vs-own gap is annotation release**, cited from the
   Annotation control track (49.08%) and reproduced independently here
   (9.16x → 5.94x).
3. **A 5.5-5.9x residual survives** allocation, deduplication, read-length
   geometry, pooling, biotype relabelling and shared-gene restriction, and is
   visible in the raw alignments. It is **not explained by this analysis.**
4. **Two candidate causes remain untested here, and both are confounded with the
   plates.** (a) *Biology*: the published plate is HEK293T + mESC, the own plate
   mouse embryo cells — different cell types genuinely differ in snRNA/snoRNA
   content. (b) *Chemistry*: the own plate is the 20x-volume manual protocol
   under development. **These cannot be separated with the data in hand**, since
   there is no own-protocol mESC library and no published-protocol embryo
   library. Do not attribute the residual to chemistry without one.

### Underpowered and unfalsifiable claims, flagged

- n=12 own-plate cells against n=173 published. The **non-overlap** is robust
  (0/173), but any *effect size* from 12 cells is weakly determined.
- The published plate is mESC+HEK293T at 0.24M reads/cell; the own plate is
  embryo at 4.40M. Depth, cell type, release and chemistry all differ together.
  Only depth (via deduplication) and release (via shared genes) were separable.
- **tRNA cannot be compared across these datasets.** It never enters the gene
  tables (separate `_tRNA.ReadCounts.tsv`; 0 tRNA rows in all three gene
  tables), and per trap 8 the `jS:IN` containment rule against ~72 bp features
  makes any cross-read-length tRNA claim unfalsifiable. Reported for the record
  only: 0.0017% (published), 0.0256% (own), 0.0003% / 0.0007% (FLASH-seq arms) of
  each non-rRNA denominator.
- **`MtRrna` remains inside the denominator**, reproducing the two-way scripts'
  behaviour (their rRNA set is `{'rRNA','Mt_rRNA'}` but the pipeline writes
  `MtRrna`). Quantified rather than silently corrected (Rule 3): removing it
  changes structural shares by <=0.0032 pp on every dataset. Kept for
  comparability with the committed two-way table.
- **Cell-calling rule matters and is stated.** Fig. 1d (UFI fraction) gives 173
  mouse cells; Methods p. 18 (gene fraction) gives 144, and doublet rate 0.57%
  vs 8.78%. Fig. 1d is primary here; the Methods rule shifts structural
  composition by **-0.88 pp**, far too little to affect any conclusion.

---

## Files

| file | contents |
|---|---|
| `composition_threeway.tsv` | **the deliverable** — biotype x 4 datasets, fractions on the non-rRNA denominator, plus gaps, folds and raw read counts |
| `composition_threeway.png` | publication figure, all three datasets |
| `threeway_denominators.tsv` | per-dataset denominator, units, genome/annotation, combination-key and species-filter shares |
| `threeway_structural.tsv` | structural-5 summary and folds |
| `threeway_structural_per_unit.tsv` | per-cell / per-library structural share (the non-overlap) |
| `threeway_allocation_rules.tsv`, `threeway_structural_provenance.tsv`, `threeway_discordant_combos.tsv` | allocation-rule invariance and where structural reads sit |
| `threeway_reads_vs_molecules.tsv` | reads-vs-molecules on both VASA plates |
| `threeway_readlength_geometry.tsv` | per-cell aligned span and realised `jS:IN`, from stage 5's own tag |
| `threeway_release_control.tsv`, `threeway_biotype_relabel_E99_E116.tsv` | the release control and the full E99/E116 relabel map |
| `threeway_sensitivity.tsv` | table choice, cell-calling rule, MtRrna denominator, tRNA |
| `threeway_published_cellcalls.tsv` | per-barcode species calls under both of the paper's rules |
| `threeway_provenance.tsv` | sha256 of every input table and BED |

Scripts in `code/flashseq_vasa/threeway/`: `threeway_composition.py`,
`threeway_mechanism.py`, `threeway_geometry.py`, `verify_threeway.py`.
`a_Mapping/` asserted clean before and after every job (Rule 1).
