# demo_analyze — reads and UFIs per barcode, two datasets

Two libraries, one set of scripts, one style:

- **own library** — PM26037 / ZHA9292A1, 16 barcodes, grouped by the
  experimental design (4 genotypes × 4 replicates)
- **published VASA-plate** — SRR14783059, 384 barcodes, grouped by species call
  (see `vasaplate/README.md`)

Tables are built per dataset; **all figures come from one script**,
`plot_all.R`, writing individual files into `figures/`. That is deliberate — the
two datasets used to have their own plot scripts, and restyling one silently
left the other behind.

## Run it

```bash
cd code/I_Gene_expression/own_version/demo_analyze

./run.sh            # everything, ~20 s
./run.sh plots      # just plot_all.R, ~10 s  <- the usual edit-rerun loop
./run.sh tables     # just the TSVs
```

**You do not need an interactive node for any of this.** Measured on the login
node: every script is 2–10 s and peaks under 500 MB.

| script | time | peak RAM |
|---|---:|---:|
| `build_tables.R` | 6.5 s | 455 MB |
| `vasaplate/build_tables.R` | 2.0 s | small |
| `plot_all.R` | 11 s | 240 MB |

The heavy steps are the ones that read FASTQ/BAM/BED. They are `sbatch`, they
are **not** part of `run.sh`, and you only run each once (same list as
`run.sh`'s header):

```bash
sbatch -c 8  --mem=8G  -t 60      --wrap="scripts/01_count_demux_reads.sh"        # 12 GB, ~1.5 min
sbatch -c 16 --mem=8G  -t 120     --wrap="scripts/02_count_demux_reads_plate.sh"  #  9 GB, ~30 s
sbatch -c 4  --mem=8G  -t 0:40:00 --wrap="scripts/05_probe_qc.sh"
sbatch -c 16 --mem=8G  -t 60      --wrap="scripts/10_mapped_length_dist.sh"       # 22 GB, ~2.5 min
sbatch -c 16 --mem=16G -t 90      --wrap="scripts/12_step5_biotype.sh"            # 7.5 GB, ~6 min
```

`r4.3` is the env to use — `reanalysis_R` has no `readxl`, which `build_tables.R`
needs to read the design xlsx. `run.sh` has the path baked in, so nothing needs
activating.

### Faster loop when you are tweaking one figure

About 2 s of every run is just loading packages. If you are iterating on a plot,
keep one R session open and re-`source()` instead:

```bash
/nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/R
```
```r
source("plot_all.R")   # edit the file, press up-arrow, run again
```

Only if you ever want a compute node anyway:

```bash
srun -c 4 --mem=8G -t 2:00:00 --pty bash
```

## 目录结构

```
demo_analyze/
  run.sh                     唯一入口 —— ./run.sh [all|tables|plots]
  scripts/                   代码，文件名按运行顺序编号
  tables/                    表格，每个数据集一个目录
    own130/ own75/ plate/ cross/
  figures/                   图，按流程步骤分组
    01_reads/ 02_trim/ 03_rrna/ 04_mapping/ 05_assign/ cross/
  probe_reference/           探针靶区参考（区间文件 + 制作脚本）
```

`scripts/paths.R` 是**唯一定义路径的地方**，从脚本自身位置推出 `ROOT`/`TAB`/`FIG`，
所以整棵树可以整体移动，脚本不用改。

## 运行步骤

| 步骤 | 脚本 | 说明 |
|---|---|---|
| **1** | `scripts/01_count_demux_reads.sh` | 每条码读段数（解压 12 GB）。**sbatch，一次性** |
| **2** | `scripts/02_count_demux_reads_plate.sh` | 同上，plate（9 GB）。**sbatch，一次性** |
| **3** | `scripts/03_build_tables.R` | own 库的全部表格；`DEMO_RUN`/`DEMO_OUT` 切换 130 nt / 75 nt |
| **4** | `scripts/04_build_tables_plate.R` | plate 的表格 |
| **5** | `scripts/05_probe_qc.sh` | 探针靶区残留（扫 Ribo.bam）。**sbatch，一次性** |
| **6** | `scripts/06_plot_all.R` | 读段 / trim / 五分类 / 长度 / 跨库对比 |
| **7** | `scripts/07_plot_insilico.R` | step 3，in-silico rRNA 去除比例 |
| **8** | `scripts/08_plot_probe_qc.R` | 探针靶区残留 |
| **9** | `scripts/09_plot_step4.R` | step 4，STAR 比对去向 |
| **10** | `scripts/10_mapped_length_dist.sh` | 按比对结果分类的读长（扫 22 GB BAM）。**sbatch，一次性** |
| **11** | `scripts/11_plot_step4_length.R` | step 4，读长 × 比对结果；多重比对率 vs 读长 |
| **12** | `scripts/12_step5_biotype.sh` | step 5，每条读段落到的 biotype 与读长（扫 7.5 GB BED）。**sbatch，一次性** |
| **13** | `scripts/13_plot_step5.R` | step 5，落到基因的比例 / biotype / 落到几个基因 / biotype × 读长 |

`./run.sh` 跑 **3–9、11、13**（几十秒，登录节点即可）。**1、2、5、10、12** 要读
FASTQ/BAM/BED，是 sbatch 一次性任务，命令写在 `run.sh` 头部注释里。

## How a design well becomes a cell id

Nothing is hard-coded to a cell number. `build_tables.R` reads the ordered-primer
design out of `../vasa cbc 6bp.xlsx`, slices the 6 nt cell barcode out of each
sequence (45 nt 5′ handle + `NNNNNN` UMI + **CBC** + polyT → positions 52–57),
and joins that on barcode to `../bc_PM26037_6nt.tsv`, which is what the pipeline
actually detected in the data and how it numbered the cells.

The barcode file writes cell ids 2-wide (`01`), every pipeline output writes them
3-wide (`001`); the script pads before joining, and hard-errors if any cell ends
up without a count — a silent join failure here just looks like an empty plot.

**The design order is not the cell order.** Sorting by design well gives a
scrambled cell sequence, e.g. XO is cells 007, 008, 009, **006** and EpiLCs is
**013**, 010, 012, 011. Read the genotype off `sample_sheet.tsv`, never off the
cell number.

### Genotype blocks

The 16 detected wells, in ascending `seq` order, are 4 consecutive replicates of
each genotype (`GENOTYPE_ORDER` in `build_tables.R`):

| genotype | design wells | cells |
|---|---|---|
| XY | seq1–4 | 002, 003, 004, 005 |
| XO | seq5, 6, 7, 9 | 007, 008, 009, 006 |
| EpiLCs | seq10, 23, 25, 26 | 013, 010, 012, 011 |
| Blank control | seq31, 46, 47, 48 | 001, 014, 015, 016 |

This assignment is confirmed independently of the design: the blank-control
block lands **exactly** on cells 001, 014, 015, 016 — the four cells the pipeline
had already flagged as 5–30× lower than the rest and blank-behaving, before any
of this was known. If that stops being true, the block assignment is wrong.

### The 17th well

The xlsx has **17** wells; `seq8` (`CACTAG`) was never detected in the data and
is dropped, with a message on stderr. The other 16 account for every detected
barcode — the script errors if a detected barcode has no design well.

## Numbers

Reads counted straight from the step-1 demultiplexed FASTQs. Total
**193,993,061** across the 16 barcodes, which agrees cell-for-cell with the run's
own `logs/step2_report.txt`; step 1 separately reports 194,316,570 "reads with
proper barcodes", 0.17% higher, so the two are counting slightly different things
and the FASTQs are the ones used here.

| genotype | total | mean/cell | sd | % of library |
|---|---:|---:|---:|---:|
| XY | 31,073,873 | 7,768,468 | 2,142,587 | 16.0% |
| XO | 58,676,637 | 14,669,159 | 3,941,033 | 30.3% |
| EpiLCs | 99,193,935 | 24,798,484 | 3,260,828 | 51.1% |
| Blank control | 5,048,616 | 1,262,154 | 366,233 | 2.6% |

Depth is very uneven — EpiLCs got **20× the per-cell depth of the blanks** and
3.2× that of XY, and half the flow cell went to EpiLCs alone. Any comparison
across genotypes has to be depth-normalised or downsampled; these read counts are
a library-balance QC, not a biological result.

The blanks are behaving as blanks on the second axis too: 5.7% of their reads
reach a gene, against 25–30% for the real cells.

## Read classes

`../diagnostics/classify_reads.py` puts every read in exactly one of five boxes. It is a
**cross of two axes**, not a flat list — collapsing them leaves reads with
nowhere to go, and the biggest single group here ("never read through *and*
died at pass 2") is exactly such a read:

- **A. did it run through into its own barcode?** recomputed with
  `trim_bc_anchor.find_anchor`, the same function pass 0 ran
- **B. where did it end up?** kept / lost at pass 1 / lost at pass 2, from
  matching read names across the three FASTQs step 2 already wrote

Reads that ran through are then split by whether pass 0 alone already left them
under the floor — separating "the insert was too short, nothing could have
helped" from "trimming inside pass 1 shortened it".

| genotype | ① 读穿保留 | ⑤ 未读穿保留 | ② 插入过短 | ③ pass 1 致短 | ④ pass 2 致短 |
|---|---:|---:|---:|---:|---:|
| XY | 2.2M (6.9%) | 13.0M (41.9%) | 0.6M (2.0%) | 0.6M (1.9%) | 14.7M (47.2%) |
| XO | 3.4M (5.9%) | 25.7M (43.8%) | 1.5M (2.6%) | 1.2M (2.1%) | 26.8M (45.6%) |
| EpiLCs | 6.8M (6.9%) | 45.8M (46.1%) | 2.1M (2.1%) | 2.0M (2.1%) | 42.4M (42.8%) |
| Blank control | 0.2M (3.0%) | 1.2M (23.4%) | **0.5M (9.1%)** | 0.1M (1.9%) | **3.2M (62.5%)** |

**Validation:** classes ① + ⑤ = 98,274,772, which is `step2_report.txt`'s `kept`
to the read. `build_tables.R` warns if that ever stops holding.

What it says: **most reads never read through at all** (⑤ alone is 42–46% in
real cells) — the insert is usually long enough that sequencing never reaches
the barcode. And the blanks are structurally different, not just shallower:
**9.1% die of a too-short insert against 2.0–2.6%**, and 62.5% die at pass 2
against 43–47%.

### What it deliberately does not split

- **③ is not "bad quality".** pass 1 is TrimGalore = adapter removal *and*
  quality trimming in one call, so a read it drops cannot be attributed to one
  or the other from the outputs. The class is named after what is known.
- **④ *is* now attributed** — see below. ③ still is not.

Nothing is re-run to produce this: it reads the three FASTQs step 2 already
wrote and derives the rest. 16 cells in parallel, ~6 min for 194M reads.

```bash
export TRIM_MINLEN=15
ls $CELLS/*_cbc.fastq.gz | sed 's/.*_\([0-9]*\)_cbc.fastq.gz/\1/' \
  | xargs -P 16 -I{} ../diagnostics/classify_reads.py --cell {} $CELLS out_{}.tsv
```

### What cut the reads class ④ lost

`../diagnostics/attribute_pass2_loss.py` credits each dropped read to the adapter that
removed the most bases from it (3' adapters take the match and everything after,
the 5' `polyT5` takes the match and everything before). Totals **87,023,863** —
class ④ to the read.

| genotype | polyA | polyT5 | rt | polyG | no_adapter |
|---|---:|---:|---:|---:|---:|
| XY | 76.5% | 20.5% | 2.1% | 0.4% | 0.6% |
| XO | 81.5% | 15.1% | 2.4% | 0.4% | 0.6% |
| EpiLCs | 85.0% | 11.6% | 2.3% | 0.5% | 0.5% |
| **Blank control** | **39.3%** | **59.1%** | 0.6% | 0.3% | 0.7% |
| whole library | 81.0% | 15.7% | 2.3% | 0.4% | 0.6% |

**The blanks fail for a different reason than the real cells.** In real cells the
poly-A tail is what leaves too little insert behind (76–85%, rising with depth:
XY 76.5% → XO 81.5% → EpiLCs 85.0%). In the blanks the dominant cause flips to
**`polyT5` at 59.1%** — reverse-orientation reads that are mostly poly-T, i.e.
empty-well artefact rather than a short transcript.

`rt` and `polyG` are negligible everywhere (≤2.4% and ≤0.5%), which is worth
knowing: the measured read-through adapter is doing its job at pass 0 and rarely
has to be caught by shape at pass 2.

It streams: the info-file would be ~90 GB over the library, so cutadapt writes
into a FIFO and nothing is stored. 16 cells in parallel, ~11 min.

## Read loss inside step 2

Where trimming loses reads, from step 2's own per-cell outputs
(`*_cbc_bcanchor.log` + `*_cbc_cutadapt.json`) — so these build whenever step 2
has run, whatever state steps 3–7 are in.

**Pass 0 is absent from the figures on purpose.** The barcode anchor has no
length filter, so it removes bases and never a whole read; the read count is flat
across it. Its base-level cost is in `trim_loss_per_cell.tsv` as `anchor_bases`
(e.g. cell 011: 641 Mbp cut). Every read step 2 loses goes at pass 1
(TrimGalore `--length 20` + Phred 20) or pass 2 (cutadapt `-m 20`).

| genotype | 解复用 | pass 1 丢失 | pass 2 丢失 | 保留 |
|---|---:|---:|---:|---:|
| XY | 31.1M | 1.2M (4.0%) | 14.7M (47.2%) | 15.2M (48.8%) |
| XO | 58.7M | 2.8M (4.7%) | 26.8M (45.6%) | 29.2M (49.7%) |
| EpiLCs | 99.2M | 4.2M (4.2%) | 42.4M (42.8%) | 52.6M (53.0%) |
| Blank control | 5.0M | **0.6M (11.0%)** | 3.2M (62.5%) | **1.3M (26.5%)** |

*(at `TRIM_MINLEN=15`; at 20 the same table read 6.7–7.6% / 44–49% / 44–49% for
the real cells, so the floor moves pass 1 much more than pass 2.)*

Two things stand out. **Pass 2 is where the library is spent** — 43–63% of reads
die on the `-m` floor, because the median insert before the poly-A is short
enough that removing the tail leaves almost nothing. And **the blanks have a
different loss profile, not just less of everything**: 11.0% fail at pass 1
against 4.0–4.7% in real cells, and only 26.5% survive against 49–53%. An empty
well is mostly short-insert product, so more reads reach their own barcode and
more fall under the length floor once cleaned — exactly what `step2_report.py`'s
header predicts.

The genotype figure is **100% stacked**, not absolute: the genotypes differ ~20×
in depth (5.0M vs 99.2M), so on a shared absolute axis the blank bar is a sliver
and its three labels collide. Filling to 100% gives every segment room for its
count and percentage; the absolute total sits above each bar.

## Why so much multimapping (step 4)

`step4_mapping.png` says the own library multimaps far more than the published
plate. Counting STAR's "too many loci" as what it is — a multimapper — the gap
is **29.0% / 33.6% against 14.2%** (`step4_mapping_from_bam.tsv`).

`10_mapped_length_dist.sh` classifies every read straight out of the STAR BAMs
(one primary record per read, `NH:i:` and `uT:A:`) and tallies its length. It
re-derives STAR's own totals exactly, on all three datasets, or it refuses to
write the table.

**Multimappers are shorter than uniquely mapped reads, in every library:**

| dataset | class | median | IQR | at full length |
|---|---|---:|---:|---:|
| own, 130 nt | unique | 129 | 110–130 | 41.5% |
| own, 130 nt | **multi** | **110** | **68–130** | 27.3% |
| own, 130 nt | unmapped | 91 | 64–128 | 15.1% |
| own, 75 nt | unique | 75 | 73–75 | 52.8% |
| own, 75 nt | **multi** | **74** | **64–75** | 42.9% |
| own, 75 nt | unmapped | 54 | 40–74 | 22.7% |
| published | unique | 75 | 73–75 | 52.0% |
| published | **multi** | **74** | **72–75** | 45.9% |
| published | unmapped | 74 | 67–75 | 35.6% |

So length is real — but it is **not** what makes this library different. The
answer is in `step4_multi_rate_by_length.png`, which asks the question at
matched length: of the reads STAR *placed* at a given length, what share went to
more than one locus?

| % of placed reads multimapped | 25 nt | 40 nt | 55 nt | 70 nt | 75 nt |
|---|---:|---:|---:|---:|---:|
| own library, 130 nt | 60.8 | 35.0 | 41.5 | 31.7 | 37.8 |
| own library, 75 nt | 60.5 | 31.3 | 33.6 | 29.3 | 32.7 |
| published, mouse wells | 41.2 | 24.4 | 15.8 | 13.8 | 13.5 |

Two things, and they are separate:

1. **Within any library, short reads multimap.** The rate falls from ~60% at
   25 nt to ~13–38% at 75 nt. That is why the 75 nt re-run multimaps *more*
   than the 130 nt one (33.6% vs 29.0%) off the same molecules — cutting the
   read shorter moves reads into the ambiguous zone.
2. **At the same length the own library still multimaps ~2.4× the plate**
   (32.7% vs 13.5% at 75 nt), and the gap widens as reads get longer. **Length
   does not explain the difference between the libraries.** If anything the
   comparison is stacked the other way: the plate is mapped against a *doubled*
   human+mouse MIXED index, which should raise its multimapping, not lower it.

The leading candidate looked like residual rDNA: this library is far more
rRNA-heavy — in-silico depletion removes **26.1%** of its reads against the
plate's 11.8%, and probe-target residual *after* depletion runs **16.68%
against 1.37%** (`insilico_depletion.tsv`, `probe_qc.tsv`, ~12×) — and a repeat
array multimaps by construction. **That hypothesis has now been tested against
step 5's own output and it is wrong**; what the multimappers actually are is
the next section.

### Reading the figures

`step4_length_by_fate.png` is **cumulative**, not a per-length curve. Each of
these libraries puts 41–53% of its reads on the single full-length bar, so on a
raw histogram that one spike owns the y axis and every class looks identical.
The cumulative form keeps the body (the rise) and the spike (the final vertical
jump) in the same picture; "which class is shorter" is just which curve sits
left. The reasoning is in the script header, not only here.

`step4_multi_rate_by_length.png` divides by **placed** reads (unique + multi),
not by all reads at that length — including the unmapped would fold mappability
into a curve meant to be about ambiguity alone. Lengths carrying fewer than
1,000 placed reads are dropped so the tail is not three-read noise.

### One denominator caveat

The plate side here is **173 mouse wells** by our own `ours_v3` calls, matching
`tables/plate/` and `04_build_tables_plate.R`. The older, hand-made
`tables/cross/step4_mapping.tsv` used the **paper's** calls instead — 172 wells,
52,062,837 reads — as does `05_probe_qc.sh`. The two sets differ by one well and
0.09% of the reads, so nothing above turns on it, but the two plate rows are not
literally the same reads. `PSOURCE=published` switches it.

Worth recording separately: **`tables/cross/step4_mapping.tsv` has no generator
script anywhere in the repo.** `step4_mapping_from_bam.tsv` is its reproducible
replacement, reconciled against STAR's logs; `09_plot_step4.R` still reads the
old file, deliberately, so the existing figure does not move without being
asked.

## What the multimappers are (step 5)

`12_step5_biotype.sh` reads the biotype out of the annotation name in the BEDs
step 5 already wrote (`*_singlemappers_genes.bed.gz`,
`*_nsorted.multimappers_genes.bed.gz`) for all three runs. Nothing is re-run.
410 BEDs / 7.5 GB, 6 min.

Everything is counted **per read**, not per row: one row is read × locus ×
overlapping feature, so a read appears 1.2–5.1 times. And a multimapper's loci
usually **disagree** — 82% of our multimapper reads touch more than one biotype
against 66% of the plate's — so there is no single "the biotype of this read".
Each read is therefore counted **once per distinct biotype it touches**, and the
columns below sum past 100%. Stacking them would assert a partition that does
not exist.

### The rDNA hypothesis is wrong

| % of that library's multimapper reads touching… | own 130 nt | own 75 nt | published |
|---|---:|---:|---:|
| **rRNA** | **0.91%** | 0.81% | **0.0018%** |
| MtRrna | 0.004% | 0.005% | 0.008% |

The direction is right — we touch annotated rRNA ~490× more often than the
plate does — but the size is not. Under 1% of multimappers touch an rRNA gene
at all, and **1% cannot move a multimapping rate from 14% to 29%.** Residual
rDNA is not the explanation.

This does not contradict `probe_qc.tsv`'s 16.68% vs 1.37%. That is measured
*before* the genome, on `Ribo.bam`; this is what is still there after step 3
removed it and STAR placed what was left.

### What does separate the two libraries

`step5_biotype_by_class.png`, and the numbers behind it in
`step5_biotype_summary.tsv`:

| biotype touched, % of multimapper reads | own 130 nt | own 75 nt | published | genes in the two BEDs |
|---|---:|---:|---:|---|
| ProteinCoding | 65.4% | 77.2% | 88.5% | 21,818 vs 21,933 |
| lncRNA | 72.2% | 71.5% | 8.0% | **32,889 vs 9,959** |
| MiscRna | 23.4% | 22.1% | 16.6% | 562 vs 562 |
| **snRNA** | **18.2%** | 16.2% | **1.2%** | 1,381 vs 1,385 |
| **miRNA** | **12.9%** | 9.8% | **3.9%** | 2,206 vs 2,207 |
| snoRNA | 4.0% | 3.5% | 1.8% | 1,507 vs 1,507 |
| **ProcessedPseudogene** | **7.2%** | 10.8% | **36.2%** | 9,312 vs 10,003 |

**The two libraries multimap for different reasons.** Ours multimaps on small
RNA — snRNA 16×, miRNA 3.3×, snoRNA 2.2×, MiscRna 1.4× the plate's rate. The
plate's multimapping is the ordinary paralog/pseudogene kind: processed
pseudogenes at 36.2% against our 7.2%. Small-RNA genes are short and come in
multi-copy families, so a read landing in one multimaps by construction —
the same logic the rDNA hypothesis had, with a different repeat family.

**Read these four ratios only alongside the length control two sections down.**
Three of them survive it and one — MiscRna — does not.

**The lncRNA row is not evidence, and its gene-count cell is bolded to say so.**
The two runs are annotated against different Ensembl releases (ours GRCm39
E116, the plate the mixed GRCh38+GRCm38 E99 BED), and those releases agree on
gene counts for almost every biotype — the right-hand column — with one large
exception:
**lncRNA, 32,889 genes against 9,959, a 3.3× difference.** A 9× gap over a 3.3×
annotation gap says nothing clean. Every other row in the table sits on biotypes
the two annotations count within 8% of each other, which is why they can be read
as library composition. This is measured, not assumed; the counts are in
`12_step5_biotype.sh`'s header.

The plate is also scored on **mouse annotation only**. It is mapped against the
MIXED reference, and 12.6–15.8% of a mouse well's multimapper BED rows are on
`GRCh38_` contigs; leaving them in would have made the human annotation part of
the answer. What that discards is written out rather than hidden —
`reads_offspecies` in `step5_assign_totals.tsv`, 3.6% of multimapper and 2.2% of
uniquely-mapped reads. `MOUSECONTIG=''` scores both species.

### Two more things the same scan shows

**Reads reaching a gene run the opposite way in the two libraries**
(`step5_gene_assignment.png`). Of the reads STAR placed, step 5 puts on a gene:

| | uniquely mapped | multi-mapped |
|---|---:|---:|
| own, 130 nt | 83.2% | **93.9%** |
| own, 75 nt | 82.7% | 93.7% |
| published | **89.6%** | 78.7% |

Our multimappers are *more* likely to land on annotation than our unique reads;
the plate's are less. That is what a small-RNA-driven multimapper population
looks like — those reads are in annotated features by definition. The reads lost
between STAR and the BED are the ones with no overlapping feature, antisense to
it (`stranded=y`), or spanning past both its ends (`jS:OUT`); `too many loci`
reads are in neither file and in no denominator here.

**Ours are more ambiguous once there.** Median distinct genes per multimapper
read is **3** (130 nt) / 4 (75 nt) against the plate's **2**, and only 9.9% of
ours resolve to a single gene against 21.2% of the plate's. The distribution
(`step5_genes_per_read.png`) is genuinely jagged — 3, 5, 7, 9 sit above their
even neighbours and there is a spike at 12 — which is what a gene-family size
distribution looks like, not noise.

Exon/intron goes the same way: **74.1%** of our uniquely-mapped assigned reads
are exon-only against the plate's 58.2%, but only **34.4%** of our multimappers
are, against 53.0%. Our multimapper population is both more intronic and more
small-RNA than the plate's.

### Does the small-RNA gap survive length-matching?

It has an innocent explanation available. snRNA/miRNA/snoRNA genes are short, a
short read fits inside one and a long read does not (step 5 drops a read
spanning past both ends of a feature, `jS:OUT`), and our reads are shorter. On
that account "more small RNA" would just be "more short reads" restated — which
is exactly the trap `step4_multi_rate_by_length.png` was built to avoid, so the
same control applies here.

Read length is recovered from the CIGAR the pipeline already parked in col 7 of
the BED (`CG:109M21S` → 130 nt), so this is still a read-only scan. It is checked
against `10_mapped_length_dist.sh`'s independent `length(SEQ)` tally out of the
BAMs: at every length the BED must be a subset of the BAM, and the script
refuses to write the table otherwise. `step5_biotype_by_length.png` is the
per-length curve; every plotted length carries ≥5,000 reads, so its structure is
real and not a thin tail.

Reading five points off a curve is not an answer, so the table below is
**directly standardised** — each own library's per-length rates reweighted by
the **plate's** read-length distribution: *what our number would be if our reads
were as long as theirs.* Window 15–75 nt, where both libraries have reads.
(The plate standardised to its own distribution reproduces its crude rate to the
decimal, which is the formula checking itself.)

| % of multimapper reads, standardised | own 130 nt | own 75 nt | published | own 75 nt ÷ published |
|---|---:|---:|---:|---:|
| **snRNA** | 5.2% | **18.0%** | 1.2% | **15×** |
| **miRNA** | 14.2% | 10.3% | 3.9% | **2.6×** |
| snoRNA | 2.4% | 3.8% | 1.8% | 2.1× |
| ~~MiscRna~~ | 14.8% | 24.6% | 16.6% | 1.5× |
| ProcessedPseudogene | 7.9% | 11.7% | **36.2%** | 0.3× |

**snRNA and miRNA survive, and are the result.** 15× and 2.6× after the length
distributions are made identical — those are not read length, they are the
library.

**MiscRna does not survive, and the earlier 1.4× should not be quoted.**
Standardised, own130 lands at 14.8% against the plate's 16.6% — *below* it. The
crude gap was read length, and `step5_biotype_by_length.png` shows why: the
plate's MiscRna curve is near zero up to ~65 nt and then spikes above 50% at
68–70 nt, a length band the plate has few reads in and we have many.

**The plate's pseudogene excess survives too** — 36.2% against 7.9%/11.7%,
unchanged by standardisation. The two libraries really do multimap on different
things.

**Use the 75 nt row, not the 130 nt one, for the headline.** own130's crude
snRNA is 18.2% over all lengths but only 3.9% inside the 15–75 nt window: its
snRNA multimappers live mostly *above* 75 nt, where the plate has no reads and
no comparison exists. `own75` is the same molecules truncated to the plate's
read length, which is what makes it the honest side-by-side.

### What is still open

Why the library is small-RNA-heavy in the first place. Length is now excluded;
what is left is the protocol and the input material, and nothing measured here
separates those two. Note also that snRNA rising with read length in our library
(the red curve climbs from ~5% at 40 nt to >30% above 100 nt) is the *opposite*
of the "short reads fit in short genes" story, and is unexplained.

## UFIs (panel C)

Counted identically on both sides: column sums of that run's
`*_total.UFICounts.tsv` over rows carrying an Ensembl gene id. On the plate that
is human-gene rows + mouse-gene rows (`vasaplate_check/vp_common.py`'s
`species_of`), dropping rows that span both species and tRNA rows; on PM26037,
mouse-only, it drops tRNA rows — 323 of 305,399.

This is **not** the same as `out/analysis/per_cell_after_filter.tsv`, which is
post gene-filtering and runs ~8–9% lower (cell 002: 491,970 there, 535,428 here).

| group | n | median UFIs |
|---|---:|---:|
| XY | 4 | 679,940 |
| XO | 4 | 1,234,698 |
| EpiLCs | 4 | 2,161,366 |
| Blank control | 4 | 54,239 |
| Human (HEK293T) | 178 | 185,176 |
| Mouse (mESC) | 173 | 72,498 |
| Mixed (doublet) | 2 | 22,044 |
| Below UFI cutoff | 31 | 534 |

y is log10 — the two datasets differ ~12× in median UFIs per barcode and the
blank/empty groups differ from real cells by another 1–2 orders. On a linear
axis everything but EpiLCs flattens into the baseline.

Note the two "empty" groups are **not** equivalent: PM26037's blanks still carry
~54k UFIs (real library, just far less of it), while the plate's below-cutoff
wells carry ~534. The plate wells are genuinely empty; the PM26037 blanks are
better described as very low-input.

## Notes

Figures are plain ggplot2 — `theme_classic()`, default colours, one plot per
file, **centred title, no subtitles**. That lives in exactly one place,
`theme_demo()` in `plot_all.R`, so a style change lands on every figure at once.
Change it there, not per plot.

**Composition figures** (trim loss, read classes, pass-2 cause) put the
**category on the fill** with a white border between segments, and read genotype
off the x axis or the facet strip. They share `stacked_abs()` / `stacked_pct()`,
so all three move together.

That replaced encoding the category on *alpha* over the genotype hue. Alpha works
at three levels and fails at five — the bottom of the ramp is not separable, and
that is exactly where the small classes (② and ③, ~2%) live. Hatching was
considered and rejected: `ggpattern` is not installed, and a texture inside a
2%-tall segment reads as noise at print size.

Genotype keeps its own palette in every figure where genotype **is** the
variable (the reads and UFI figures); it is only displaced where a second
categorical dimension has to be shown at the same time.

`reads_by_stage_genotype.tsv` is **not** a monotone funnel: "Assigned to genes"
is larger than "Uniquely mapped" because the pipeline keeps multimappers
(`--outFilterMultimapNmax 20`) and rescues them at assignment, so a read can
reach a gene without ever mapping uniquely.
