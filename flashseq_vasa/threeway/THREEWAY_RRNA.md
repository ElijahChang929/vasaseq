# Residual rRNA across all three datasets, on one measurement method

Written 2026-07-30. Every number here is re-derived and asserted by
`03_threeway_table.py`; `res/threeway/threeway_report.txt` is that script's
output. Nothing in this file is transcribed from an earlier report.

## What was measured, and on what

One method for all three: bwa (`aln` + `mem`) against an rRNA reference
containing the mouse 47S unit `BK000964.3`, classified by the published
`riboread-selection.py` predicate, each dataset under its own correct strand
flag. No dataset's rRNA calls were recomputed here — this reads what those runs
produced and puts them on one axis.

| | published VASA-plate | own VASA-plate | FLASH-seq |
|---|---|---|---|
| id | SRR14783059 | ZHA9292A1 (PM26037) | RN26038 |
| units | 384 barcodes | 16 barcodes (12 real + 4 blank) | 10 libraries |
| strand flag | `y` | `y` | `n` |
| rRNA reference | `unique_rRNA_human_mouse.v3.fa` (921 entries) | `unique_rRNA_mouse.v2.fa` (357) | `unique_rRNA_mouse.v2.fa` (357) |
| genome / annotation | GRCm38 + hg38, Ensembl 99 | GRCm39, Ensembl 116 | GRCm39, Ensembl 116 |
| read length | 74 nt | 151 nt | 151 nt |
| sampling | all reads | all reads | stride 1/64, R1 only |

## Denominator — the same quantity on all three sides

**All three are reads entering stage 3**, i.e. the `Number of reads:` line of
that unit's own `.ribo-map.log`. This was checked per unit rather than assumed:

- published: `depletion_v2_vs_v3.tsv` `input_reads` == log `nreads`, 384/384 cells
- own plate: `rrna_comparison.tsv` `reads_in` == log `nreads`, 16/16 barcodes
- FLASH-seq: `rrna_comparison.tsv` `reads_in` == log `nreads`, 10/10 libraries

So no recomputation was needed and both existing figures stand as they are.
Two scale properties travel with the FLASH-seq denominator and are carried as
columns in the table, not only as prose: it is a deterministic 1/64 stride
sample, and it is R1-only single-end (which is what makes it comparable to VASA,
where there is one biological read per fragment).

One residual, documented rather than smoothed: `ribo_v3` exceeds the log's
summed mapped count by **exactly +1 read per cell** (384 reads plate-wide,
0.00406% of the plate's ribosomal total). `riboread-selection.py` never flushes
the final read group of a file, so `nreads = nunmapped + mapped + 1` by
construction, and the depletion table assigned that read to `ribo`. It moves the
plate figure from 5.219028% to 5.219240%. The table uses the depletion table's
value, which is the brief's anchor and which reproduces exactly:
**5.2192% = 9,461,806 / 181,287,059**.

## Headline numbers

Read-weighted within group, each side under its own correct flag:

| group | n | rRNA % (read-wtd) | median | range |
|---|---|---|---|---|
| published, all 384 wells | 384 | **5.2192** | 4.67 | 0.36–49.06 |
| published, mouse (Fig. 1d rule) | 172 | **11.7728** | 11.12 | 0.60–32.29 |
| published, human (Fig. 1d rule) | 178 | 2.0065 | 1.80 | 0.36–8.19 |
| published, mouse (Methods rule) | 158 | 11.4942 | 11.03 | 0.60–32.29 |
| own plate, real cells | 12 | **21.5733** | 21.14 | 17.99–25.09 |
| own plate, all 16 barcodes | 16 | 21.3926 | 20.68 | 5.86–25.09 |
| FLASH-seq, 9 libraries (A8 excluded) | 9 | **4.7481** | 4.70 | 3.50–6.44 |

Ratios:

- own real cells / published **mouse** cells = **1.83x** read-weighted, 1.90x on medians
- own real cells / FLASH-seq = **4.54x**
- published mouse / FLASH-seq = **2.48x**
- published **mouse / human = 5.87x** — a species effect inside one library, on
  one reference, in one pipeline run

That last row is the most important control in the table. It is larger than
every cross-dataset difference reported here, which sets the scale against
which those differences should be judged.

## The species split (brief item 1)

Rule used for the headline: **the paper's Fig. 1d rule**, as implemented in
`vasaplate_check/vp_common.py::classify_fig1d` — a >=7,500-UFI gate, then species
by UFI purity with a 25/75% doublet band. Both rules are carried in the table
(`species_fig1d`, `species_methods`) because they disagree:

| rule | mouse | human | mixed | discarded |
|---|---|---|---|---|
| Fig. 1d (UFI purity) | 172 | 178 | 3 | 31 |
| Methods (gene purity) | 158 | 178 | 17 | 31 |

The 14 disputed wells are all fig1d-mouse → methods-mixed. Using the Methods
rule instead moves the published mouse figure from 11.7728% to 11.4942% — a
0.28 pp change, so the headline is not sensitive to the choice even though the
doublet counts differ ~6x.

## Where the residual sits — and one correction to my own method

`rrna_comparison.tsv` carries subunit composition for **FLASH-seq only**; its
eight `pct_*` columns are empty for both VASA plates. Subunit composition on the
VASA sides is new work here (`01_vasa_subunits.py`), computed by the same
bucketing `05_rrna_bwa_report.py::composition()` and
`own_version/step3_report.py::parse_bam()` use, with FLASH-seq re-run through the
same function so all three come from one code path.

**I initially asserted that composition was comparable across all three because
the 47S record is byte-identical in both reference files. The data falsified
that.** The record is identical; the competition for reads is not. The published
plate's reference also holds 5 human 45S records, and 18S/28S are strongly
conserved between human and mouse while the 5'ETS spacer is not — so the human
records compete for mouse 18S/28S reads and not for mouse 5'ETS reads. Measured:
**8.47%** of mouse-well ribosomal reads land on human 45S records (a pure mESC
well has no human rRNA), and the raw profile reads 18S **1.1%**, which is not a
biological result.

`02_published_realign.py` therefore realigned those same reads — numerator held
fixed, so this is a composition question and produces no new rRNA percentage —
to the mouse-only reference the other two used. 99.53% of them map there.

Composition as % of 47S-derived reads:

| | 5'ETS | 18S | ITS1 | 5.8S | ITS2 | 28S | 3'ETS |
|---|---|---|---|---|---|---|---|
| published mouse, **as run** (mixed ref) | 79.7 | 1.1 | 4.1 | 1.7 | 1.8 | 11.5 | 0.1 |
| published mouse, **realigned** (mouse ref) | **71.1** | **4.1** | 5.9 | 1.6 | 1.6 | **15.5** | 0.1 |
| own plate, real cells | 18.5 | 6.4 | 5.0 | 2.1 | 8.2 | **58.8** | 1.0 |
| own plate, blanks | 18.5 | 7.9 | 3.7 | 1.4 | 5.5 | 62.4 | 0.6 |
| FLASH-seq, 9 libraries | 44.6 | 8.7 | 4.1 | 2.0 | 6.3 | 34.2 | 0.0 |

Only the realigned row is comparable with the lower three. The as-run row is
kept so the size of the reference effect is visible: realignment moved 18S from
1.1 to 4.1% and 28S from 11.5 to 15.5%, so the effect is real — but it did **not**
overturn the 5'ETS dominance (79.7 → 71.1%). **Reference competition is not what
makes the published plate look different.**

The three-way pattern is a genuine ordering: the published plate's residual is
5'ETS-dominated, the own plate's is 28S-dominated, FLASH-seq sits between them.

### The 5'ETS excess is not the known poly-T artefact

This had to be checked, because a 71% 5'ETS share is exactly what the previously
documented poly-T artefact produced (it inflated 5'ETS to 31–62% of a blank
cell's ribosomal calls). That incident was caught because 88.9% of the hits sat
in **one** 200 nt window. Peak-bin share here:

| | median | range |
|---|---|---|
| published mouse | 27.6% | 17.97–37.35 |
| own plate, real | 11.2% | 9.89–14.03 |
| own plate, blanks | 20.5% | 11.67–32.95 |
| FLASH-seq | 13.8% | 12.31–14.98 |

The published plate is elevated, but its 5'ETS signal is **broad, not spiked**:
its top bins (600/800/1000 nt) are contiguous, 3 bins hold 50% of its 5'ETS
reads, and its positional profile correlates **r = 0.745** with the own plate's.
A single-window artefact would show none of that. So the 5'ETS dominance is a
real property of those reads.

### 5S — measured, and essentially absent

5S is Pol III-transcribed and is **not part of the 47S unit**, so it can only be
seen through the dispersed reference entries. The mouse-only reference names 74
such entries; the mixed reference's mouse entries carry no gene symbol, so 5S
was recovered there by exact sequence identity (98 entries matched). Result:

| | 5S reads | % of ribosomal |
|---|---|---|
| published 384 | 347 | 0.0037% |
| own 16 | 522 | 0.0027% |
| FLASH-seq 10 | 9 | 0.0044% |

Under 0.01% on all three. This is a measurement, not a missing value — but note
it is a **bwa-alignment** measurement against ~120 nt reference entries, and it
is not evidence about 5S in the count tables, where the separate feature-length
containment trap applies.

## The reference asymmetry is larger than "human entries added"

Flagging this as the brief asks, with the size measured rather than asserted.
The two references are not the same reference plus human sequence. Matching the
357 mouse-prefixed entries of the mixed file against the 357 entries of the
mouse-only file **by Ensembl accession**:

- **16 of 357 accessions are shared**; 1 of those 16 differs in sequence
- the 47S record `mouse_rDNA_47S_BK000964.3_1-13403` is byte-identical (13,403 nt)

The dispersed gene sets are therefore almost entirely disjoint — the mixed
reference draws them from Ensembl 99 / GRCm38, the mouse-only one from Ensembl
116 / GRCm39. Consequences:

- the 47S-derived subunit shares are comparable (identical reference record,
  and after realignment, identical competition)
- `pct_other`, `pct_5S` and `pct_human45S` are **not** comparable across the
  reference boundary, and are kept as separate columns rather than folded into
  the comparable set
- the non-47S fraction differs a lot (published 34.1%, own 7.6%, FLASH-seq
  10.7%), and most of the published plate's is the human 45S bucket (25.5% of
  its ribosomal reads) — an artefact of reference content, not of chemistry

## Confounds that are not resolved by this analysis

Stated plainly because the three-way gap mixes them:

1. **Annotation release.** Published is Ensembl 99 / GRCm38; the other two are
   Ensembl 116 / GRCm39. For the 47S-derived numbers this does not bite (the
   47S record is identical and is not part of either Ensembl release), which is
   why the subunit comparison is the more defensible one here. It does bite the
   dispersed buckets.
2. **Read length.** 74 nt vs 151 nt. A shorter read maps to a conserved subunit
   more promiscuously, so read length and reference competition push in the same
   direction and this design cannot separate them. Not resolved.
3. **Biology.** Published is HEK293T + mESC (cultured lines); own is mouse
   embryo. rRNA content is not a constant of nature.
4. **Strand flag.** VASA `y`, FLASH-seq `n`, each correct for its chemistry.
   Any VASA-vs-FLASH-seq ratio must be quoted as "each under its own correct
   flag"; the earlier sensitivity analysis (both `n` = 6.08x, both `y` = 8.95x
   for the own-plate comparison) still applies.

## What I would treat as the finding

The own plate's residual rRNA is **1.83x** the published plate's mouse cells —
smaller than the 2.3x figure from the earlier 8-cell comparison, and now
measured on all 384 published barcodes with the corrected v3 reference and split
by species. But the published plate's own human-vs-mouse spread is **5.87x** on
one reference in one run, so a 1.83x cross-dataset difference is well inside the
range that cell type alone moves this number.

The subunit result is the more robust of the two, because it survives the
reference correction and does not depend on the annotation release: **the own
plate's residual is 28S-dominated (58.8%) while the published plate's is
5'ETS-dominated (71.1%)**. Those are different failure modes of depletion, not
different amounts of the same one, and that is a testable difference in what the
depletion step is leaving behind.

## Files

- `res/threeway/rrna_threeway.tsv` — 411 rows, one per unit, 26 columns
- `res/threeway/subunits_percell.tsv` — per-unit subunit counts, all three datasets
- `res/threeway/ets_profile.tsv` — 200 nt-binned 5'ETS positional profiles
- `res/threeway/published_realigned_47S.tsv` — realigned published-mouse composition
- `res/threeway/rrna_threeway.png` — the figure
- `res/threeway/{threeway,subunits,realign}_report.txt` — the three scripts' own logs
- `code/flashseq_vasa/threeway/` — `01_vasa_subunits.py`, `02_published_realign.py`,
  `03_threeway_table.py`, and this file
