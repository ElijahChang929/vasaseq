# vasaplate — the same plot, on the published VASA-plate library

The figures from `..` (PM26037 / ZHA9292A1) redrawn for the reference dataset:
**`SRR14783059` / `GSM5369495`, `vasaplate-HEK293T-mESC`** — the only VASA-plate
library in GSE176588, 384 CEL-seq2 barcodes.

## Run it

```bash
cd code/I_Gene_expression/own_version/demo_analyze/vasaplate
sbatch -c 16 --mem=8G -t 120 --wrap="./count_demux_reads.sh"   # once; ~30 s
```

Everything after that is driven from the parent:

```bash
cd .. && ./run.sh
```

**There is no plot script here.** All figures for both datasets come from
`../plot_all.R` and land in `../figures/` as `plate_*`. That is on purpose: when
each dataset had its own plot script, restyling one left the other behind.

This folder holds only the counter, `build_tables.R`, and its TSVs. To rebuild
just the tables (~2 s, login node is fine):

```bash
/nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/Rscript build_tables.R
```

`count_demux_reads.sh` measures, `build_tables.R` makes tables. Nothing here
reads or writes anything in `..`; the dependency runs the other way, with
`../plot_all.R` reading this folder's `reads_per_cell.tsv`. So build this folder
before plotting — `../run.sh` does both in the right order.

## What replaces "genotype"

**This library has no genotype design** — it is a species-mixing control, HEK293T
(human) and mESC (mouse) on one plate. The grouping that plays the role genotype
played in PM26037 is the **species call**, read from
`res/vasaplate/per_cell.tsv`, `source = ours_v3`:

| group | n | what it is |
|---|---:|---|
| Human (HEK293T) | 178 | |
| Mouse (mESC) | 173 | |
| Below UFI cutoff | 31 | under 7,500 UFIs — empty/near-empty wells, the nearest thing to a blank control |

**`Mixed (doublet)` is excluded** — 2 wells, 334,101 reads, 0.17% of the library.
Two wells cannot support a comparison, and at that size every bar label rounds to
"0.0M", which reads as an error rather than a number. `DROP` in `build_tables.R`
does it, and the run prints what it removed.

⚠️ **The denominator changed with it.** Everything in this folder is now a share
of **382 wells / 201.4M reads**, not 384 / 201.8M. The barnyard numbers in
`vasaplate_check/` keep all 384 — the doublet rate is the whole point there — so
never carry a percentage between the two.

`ours_v3` is the anchor run — see
`vasaplate_check/BENCHMARK_published_plate.md`, "Use `vasaplate_out_v3` as the
published-plate anchor". Do **not** switch this to `bedv2` or to run 1.

### The doublet rule matters

The paper states **two** doublet rules that disagree ~6× on this data: Fig. 1d
thresholds the UFI fraction, Methods p. 18 thresholds the gene fraction (0.85%
vs 4.82%). These figures use **`call_fig1d`** (UFI fraction). Switching rules is
one line — `RULE` in `build_tables.R` — but any barnyard number quoted from here
must name which was used.

⚠️ **The figures no longer say which rule made them.** The subtitle that carried
that was dropped when the styling was matched to the parent folder. The rule is
recorded here and in `build_tables.R`, but a PNG on its own does not carry it —
so name the rule whenever one of these leaves the folder.

## Numbers

**201,423,886 reads over 382 barcodes** (384 counted, 2 dropped as mixed),
from `vasaplate_out/*_cbc.fastq.gz`. Before the drop it is 201,757,987; step 1's
own log reports 201,775,925 kept, 0.009% higher than that — the same small
bookkeeping gap seen in the parent folder, not a missing file.

Stages 1–4 are symlinked between the run directories, so `vasaplate_out/` is the
right source for read counts no matter which downstream run you are looking at.

| group | total | mean/barcode | median | % of library |
|---|---:|---:|---:|---:|
| Human (HEK293T) | 135,549,069 | 761,512 | 657,699 | 67.3% |
| Mouse (mESC) | 65,611,308 | 379,256 | 310,040 | 32.6% |
| Below UFI cutoff | 263,509 | 8,500 | 3,841 | 0.13% |

**HEK293T barcodes carry ~2× the reads of mESC barcodes** (761k vs 379k mean),
so the plate is not balanced by read depth even though it is nearly balanced by
barcode count (178 vs 173). Expected for a mixing control — HEK293T is a much
larger cell — but it means depth must be handled before any cross-species
comparison.

The empty wells are cleanly separated: 8.5k mean reads against 380–760k, i.e.
~45–90× lower. Compare with the PM26037 blanks at 5–30× lower.

## Trimming — and why only two of the four figures transfer

This library was trimmed by **`a_Mapping/trim.sh`, the published pipeline**, not
by `own_version/`. Two differences decide what can be drawn:

- **There is no pass 0.** Upstream never anchors on the cell barcode, so the own
  library's split into "read through into its own barcode" and "did not" has no
  counterpart. It would also be empty: measured on well 005, only **0.27–0.51%**
  of these reads reach their own barcode (either barcode order), against ~25% in
  the own library. Reads here are **75 nt**; the own library's are 130 nt, so
  read-through is an own-library problem, not a VASA problem.
- **pass 2 is four 6-mer homopolymers** (`polyG1`/`polyC1`/`polyT1`/`polyA1`) at
  the default `-n 1`, with no `--poly-a` and no 5' adapter. Confirmed against the
  run's own cutadapt log, not assumed.

So the plate gets the two figures that do transfer, and not the five-class one.

### Where reads are lost

| group | demultiplexed | pass 1 | pass 2 | kept |
|---|---:|---:|---:|---:|
| Human (HEK293T) | 135.5M | 0.2M (0.2%) | 13.6M (10.0%) | 121.7M (89.8%) |
| Mouse (mESC) | 65.6M | 0.1M (0.2%) | 6.4M (9.8%) | 59.1M (90.0%) |
| Below UFI cutoff | 0.3M | 806 (0.3%) | 40k (15.3%) | 0.2M (84.4%) |

**Upstream keeps ~90% of this library; `own_version` keeps ~50% of PM26037.** The
gap is not a pipeline defect — it is read length and insert length. At 75 nt with
long inserts, trimming rarely takes a read under 15 nt; at 130 nt with a median
insert of ~10 nt before the poly-A, it usually does.

Note pass 1 removes almost nothing here (0.2%) against 4–11% on the own library.

### What cut the reads pass 2 dropped

| group | polyA1 | polyT1 | polyG1 | polyC1 |
|---|---:|---:|---:|---:|
| Human (HEK293T) | 59.8% | 24.1% | 13.0% | ~3% |
| Mouse (mESC) | 61.3% | 19.1% | 15.0% | ~5% |
| Below UFI cutoff | 67.4% | 22.0% | 8.8% | ~2% |

Poly-A dominates here too, but far less lopsidedly than in the own library
(81%): upstream's `AA{5}` fires on any 6 A's, so it also catches ordinary A-rich
sequence, and `polyT1`/`polyG1` take a much larger share. **`polyC1` is not
negligible here (2–5%)** — the one adapter `own_version` dropped as having no
mechanism. On a 6-mer rule it fires anyway, which is consistent with it being an
artefact of the rule rather than of the chemistry.

Unlike the own library, the empty wells are **not** dominated by poly-T
(22.0% vs the own blanks' 59.1%): with no 5' poly-T adapter and no read-through,
upstream simply does not see that failure mode.

## Notes on the figures

Drawn by `../plot_all.R` into `../figures/` as `plate_reads_per_cell` and
`plate_reads_per_group`, in the shared `theme_demo()` style — plain ggplot2,
`theme_classic()`, default colours, centred title, no subtitles. The group
figure carries no title at all.

`plate_reads_per_cell` differs from the own-library one in two ways, both forced
by 384 bars instead of 16: bars are unlabelled (per-well numbers are in
`reads_per_cell.tsv`) and sorted by reads within each group. Facets are **equal
width, not `space = "free_x"`** — the n = 2 and n = 31 groups otherwise collapse to slivers
with clipped strip labels. Equal widths make small n obvious instead: two fat
bars read as "two barcodes".
