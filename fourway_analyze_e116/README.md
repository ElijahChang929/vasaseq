# fourway_analyze_e116 — the same comparison, one annotation

Products only. The code is `../fourway_analyze/scripts/`, run with two env vars
pointing at the re-mapped plate and at this directory. There is deliberately no
second copy of any script here.

**Read `../fourway_analyze/README.md` first.** Everything it says about the four
datasets, the standardisation, and what is *not* comparable still applies. This
file records only what changed and what the annotation-matched numbers are.

## Why this exists

Three of the four datasets were quantified on **Ensembl 116 / GRCm39**; the
published plate was on **Ensembl 99 / GRCm38** (human+mouse). That single
difference sat on every biotype axis, and it was not small: E116 carries
**32,889** mouse lncRNA genes against E99's **9,959** — 3.30×.

The measurement that forced the issue: applying the paper's own cell filter
(`b_Analysis/filterParams.py`, `lncRNA` window `[0.01, 0.03]`), the plate at
2.66% landed inside it and our libraries at 8.3–8.4% were three times outside.
Either our chemistry was very different, or the annotation was. It was both, and
until the plate moved to E116 there was no way to say in what proportion.

So the plate's 173 mouse wells were re-mapped: `scripts/00b_plate_e116.sh`.
Stages 1–3 are reused — rRNA depletion does not depend on the annotation — and
only stages 4–7 re-run, against `own_version/config.sh`'s own defaults, so the
plate ends on **literally the same STAR index and BED** as the other three.

## Mouse-only, and what it costs — measured, not assumed

The plate is a HEK293T/mESC mixing control, so a mouse-only index gives human
reads nowhere correct to go. Rather than argue about it, two human-called wells
were mapped against this exact index:

| well | call | input reads | uniquely mapped |
|---|---|---|---|
| 225 | human (97.3% human by UFI) | 2,824,850 | **16.88%** |
| 121 | human (97.8%) | 2,198,198 | **17.24%** |
| 002 | mouse | 46,559 | 76.28% |
| 004 | mouse | 491,999 | 81.08% |

Human wells still hold ~2.5% mouse, which at 78% accounts for ~1.9 pp, so
roughly **15% of genuinely human reads mismap onto mouse**. Folded into each
mouse well's own off-species share, the spurious fraction of its uniquely mapped
pool is **1.08% at the median**, 1.79% at the 90th percentile, 6.56% at worst —
against biotype differences of 5–15 pp.

Kept rather than filtered. The alternative — a mixed GRCh38+GRCm39 index —
would put the plate on a *different* index from the other three and trade a
measured 1% for an unmeasured asymmetry in multimapping competition. All 173
mouse wells are carried, matching `datasets.sh`'s plate definition; the ten
wells above 10% off-species can be dropped at table level.

## Two settings that are required, not optional

Both were found by getting them wrong, and both failed in ways that did not say
what was wrong. They are now in `run.sh`'s E116 recipe.

**`MOUSECONTIG=''`** — script 04 restricts plate rows to `^GRCm38_`. The E116
plate's contigs are plain `1`, `2`, …, so the filter matched nothing and
discarded all 173 wells. Every table was still written, with `bed_reads=0`, and
the run got 17 minutes further before dying inside an awk summary with
`division by zero`. 04 now checks for a dataset×class where STAR placed reads
and none reached a gene, and exits naming the variable.

**`PLATE_BED` + `PLATE_PREFIX=''`** — worse, because it did not look broken.
Run against the E99 mixed BED, the E116-mapped plate still produced a full set
of curves — from GRCm38 exon coordinates describing GRCm39 reads. `genes_covered`
fell **14,936 → 1,505**, and since the figure uses genes shared by all four
datasets, the shared set fell **11,377 → 827**: every curve, not just the
plate's, built from a seventh of the genes. Exit code 0. 07 now refuses to write
when any dataset covers fewer than a quarter of the median.

Same root cause both times: a reference-specific constant that degrades to
silence, not to an error, when the reference changes.

## Results

### The comparison this folder was built for

At **75,000 assigned reads per cell** — the paper's quoted depth:

| dataset | genes/cell | n |
|---|---|---|
| **VASA own, 75 nt** | **12,054** | 12 |
| VASA own, 130 nt | 11,823 | 12 |
| FLASH-seq | 10,100 | 9 |
| VASA published | 9,058 | 142 |

Our library detects **33% more genes than the published plate at matched depth**
— a comparison that could not be made before, on annotations differing 3.30× in
lncRNA gene count.

That plate number is also the strongest end-to-end check here. The paper quotes
VASA-plate at **9,480 ± 1,252** genes per cell at 75,000 **trimmed** reads; we
measure **9,058** at 75,000 **assigned** reads, which is the stricter depth
since assignment loses reads. Independent pipeline, independent annotation,
independent implementation of the downsampling — inside their s.d.

x is *assigned* rather than *trimmed* reads deliberately: `saturation_qc.tsv`
carries assigned/trimmed per dataset so the two can be converted, and using it
in the axis would push a per-dataset ratio into the comparison.

### Biotype composition, now on one annotation

Per-cell UFI fractions, `uniaggGenes_total`, the paper's own `ubiotype` rule
(mixed-biotype rows counted in the denominator, in no numerator):

| | ProteinCoding | lncRNA | smallRNA |
|---|---|---|---|
| plate, E116 (173 wells) | 87.50% | **5.49%** | 1.78% |
| plate, E99 (for reference) | 91.97% | **2.66%** | 1.15% |
| own130, real cells (12) | 78.13% | **8.29%** | 6.36% |
| own75, real cells (12) | 75.50% | **8.41%** | 6.77% |

**The annotation explains about half the lncRNA gap and no more.** Moving the
plate from E99 to E116 doubles its lncRNA fraction, 2.66% → 5.49% — but we are
still 1.5× above it. ProteinCoding behaves the same way: the deficit narrows
from 14–16 pp to 9.4–12 pp and does not close.

This agrees with `flashseq_vasa/threeway/E99_MATCHED.md`, which attacked the
same confounder from the opposite direction (our 16 cells re-mapped to E99) and
concluded release explains **15.3%**, leaving 15.54 pp real. Two independent
routes, same verdict: the annotation is a major contributor and the residual
difference is real.

**Do not transplant the paper's absolute thresholds.** Under E116 our real cells
sit at 8.3–8.4% lncRNA against a published window of 1–3%. Use the paper's
*axes* and recalibrate: on this data lncRNA < ~12% and ProteinCoding > 70%
separate real cells from blanks cleanly (blanks: 20.4–21.7% lncRNA, 65–66% PC).
The `smallRNA` window `[0.05, 0.15]` does **not** transfer at all — the plate
measures 1.78% here and the paper's own text states 1.4% for VASA-plate, both
far below its floor. That window is specific to the VASA-drop embryo data.

### Gene body coverage — the finding survived intact

| dataset | 3'/5' | intronic, excluded |
|---|---|---|
| VASA own, 130 nt | 1.01 | 25.42% |
| VASA own, 75 nt | 0.90 | 26.50% |
| VASA published | 0.87 | **41.48%** |
| FLASH-seq | 0.90 | **8.48%** |

All four flat — no 3' bias anywhere, including FLASH-seq across its whole input
titration down to 30 pg. The plate's intronic fraction is 41.48% under E116
against 43.69% under E99, so **the polyA-vs-total-RNA result was never an
annotation artefact**. The other three datasets are numerically identical to the
E99 run, which is the control that says only the plate changed.

Shared gene set **11,586** (E99 run: 11,377).

### Unchanged by the remap

rRNA and probe-scoped depletion are stage-3 measurements and the stage-3
products are shared by symlink, so these tables are byte-identical to the E99
run — which is itself a check that the plumbing is right:

| dataset | rRNA % | probe-target residual |
|---|---|---|
| VASA own, 130 nt | 26.09 | **16.68%** |
| VASA own, 75 nt | 26.44 | 16.78% |
| VASA published | 11.79 | **1.38%** |
| FLASH-seq | 4.76 (unstranded) | 1.95% |

Our RNase H reaction left about **12× more probe-target rRNA** than the
published one. That remains the most actionable number in this folder.

Mapping did change for the plate, as it must:

| dataset | multi % | unmapped % |
|---|---|---|
| VASA own, 130 nt | 29.04 | 12.43 |
| VASA own, 75 nt | 33.63 | 10.06 |
| VASA published | 13.92 | 6.73 |
| FLASH-seq | 12.23 | 8.82 |

(The plate was 14.22% multi on the mixed E99 index; 13.92% here. Dropping the
human genome removes competing loci, so the small fall is expected.)

## Layout and how to re-run

`figures/<step>/` 16 PNG+PDF, `tables/cross/` the TSVs, `logs/` the SLURM logs.
Figure group numbers match `../fourway_analyze/`.

```bash
cd ../fourway_analyze
V=/nemo/lab/turnerj/working/guangxin/vasaseq
export PLATE_CELLS=$V/data/ref/fastq_vasaplate/plate_e116/cells
export FOURWAY_OUT=$V/code/fourway_analyze_e116
export MOUSECONTIG='' PLATE_PREFIX=''
export PLATE_BED=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed

scripts/00b_plate_e116.sh setup && scripts/00b_plate_e116.sh submit   # stages 4-7
# the figures need only stages 4-5, NOT step 7 -- about 30 min, not 10 hours
scripts/01_insilico_depletion.sh
sbatch -c 8  --mem=8G  -t 4:00:00 --wrap="scripts/02_probe_qc.sh"
sbatch -c 16 --mem=8G  -t 120    --wrap="scripts/03_mapped_length_dist.sh"
sbatch -c 16 --mem=64G -t 240    --wrap="scripts/04_step5_biotype.sh"
sbatch -c 4  --mem=8G  -t 240    --wrap="$PY scripts/07_genebody_coverage.py"
sbatch -c 4  --mem=16G -t 240    --wrap="$PY scripts/09_saturation.py"
./run.sh
```

`$PY` is `/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/python3`, named
explicitly: the compute nodes' default `python3` is 3.6.8 and the shebang
resolves to whatever the job inherited, which killed both scanners once.

Measured cost: stage 4 **20m36s**, stage 5 **6m46s**, stage 6 **2h59m**, stage 7
**1h25m** (173 wells; the full 384-well plate was 8h54 + 4h38). Scanners: 03
9m44s, 04 19m25s, 07 31m13s, 09 10m11s, figures 3m04s.
