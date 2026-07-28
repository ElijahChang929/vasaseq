# Which VASA table these FLASH-seq tables must be compared against

Decision record. The FLASH-seq tables built by `pipeline_fs.sh` are only half of
a comparison, and pairing them with the wrong VASA table would produce a
difference that is an artefact of two different filter sets rather than of the
two protocols. This file names the comparator and says why.

## Use the UNFILTERED uniagg VASA tables, not `data/PM26037/out/analysis/`

**Comparator column: `V2_vasa_uniagg_spliced`** — uniagg, exon-only, **no
UMI-ceiling filter**. All five variants are tabulated in
`res/flashseq_vasa/pipeline_vs_protocol.tsv` so the choice stays explicit.

VASA's analysis set dropped 8 entries that had saturated the 4096-UFI ceiling:
*Rmrp*, *Rn7sk*, *Rn7s1*/*Rn7s2*, *Rnu1a1*, *Rnu1b6*, *Rnu2.10*, the *Snord3b*
cluster, and a *Cmss1* pair. Measured by the nf-core cross-check track
(commit `e2906ee`): those 8 entries map to **1,372 uniagg rows carrying
12,207,080 reads = 23.02% of the uniagg total**, and dropping them moves the
biotype composition by:

| class | change |
|---|---|
| ProteinCoding | 64.13% -> 83.31% (**+19.19 pp**) |
| misc_RNA | -8.86 pp |
| snRNA | -4.07 pp |
| lncRNA | -4.69 pp |

That filter is correct for molecule counting, which is what it was written for:
those entries were selected *for* saturating the UMI space, i.e. for being the
most abundant small RNAs, and a clamped `bc2trans` value is not a transcript
count. But it is **wrong for a read-based composition**, and the classes it
removes are precisely the ones a total-RNA protocol is supposed to capture — so
applying it would delete the effect the comparison exists to measure, and would
flatter FLASH-seq by ~19 pp of protein-coding share.

## The symmetry, which is the point

The FLASH-seq side has **no ceiling filter to apply, and cannot have one**:
`bc2trans` clamps at UFI >= K, K = 4**len('A') = 4, and UFI on this path is at
most 1. `build_analysis_tables.py` asserts that precondition rather than the
conclusion and records the result in `manifest.json`
(`filters_not_applicable.umi_ceiling`).

So both sides are unfiltered for the same quantity, and the remaining filters
are the same two on both sides (unnamed row; plus all-zero rows here, which
VASA's 16-column frame could not produce). Pairing an unfiltered FLASH-seq table
with VASA's *filtered* analysis set would be the one combination that is not a
comparison.

## Two further things not to charge to the protocol

- **28.76% of VASA's uniagg reads are unspliced** (15,249,326 / 53,022,536), and
  lncRNA is 82.2% unspliced. `featureCounts -t exon` and RSEM's transcriptome
  cannot count those at all. This does not affect a VASA-route-vs-VASA-route
  comparison — both arms here go through the same step 5/6/7 — but if any
  nf-core number is quoted alongside these tables, that difference is an
  annotation-model difference and must not be attributed to the chemistry.
- **nf-core's "rRNA" biotype figure is one 18S relic locus** (*Rn18s-rs5*,
  carrying 99.93-99.999% of the class), not 5S as
  `code/flashseq/README.md` states. Label it "annotation-route 18S relic",
  never "rRNA". The rRNA figure for this comparison comes from the bwa route in
  `res/flashseq/rrna_bwa.tsv`; see `NOUMI_PATH.md` 5c for why the annotation
  route reads low on both sides.

## The spliced+unspliced hazard is asymmetric, and that is now measured

`NOUMI_PATH.md` 3h flagged that the no-UMI `countExonReads`/`countIntronReads`
branches test **substring** containment (`'exon' in k`) where the vasa branch
tests **exact membership**, so a combination label `exon-intron` can be counted
in both and `spliced + unspliced` can exceed `total`.

On the **vasa** branch this does not happen at all: `max(spliced + unspliced -
total) = 0` over all 222,420 uniagg entries x 12 cells (measured, same track).
So a non-zero excess on the FLASH-seq side is **this branch difference**, not a
general property of the tables, and it must not be read as a data problem.
`reconcile.py` reports the per-library excess either way.
