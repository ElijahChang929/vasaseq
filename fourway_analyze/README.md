# fourway_analyze — the same figures, four datasets

Everything `demo_analyze/` measures on three datasets, measured again with
FLASH-seq added as a fourth. **Nothing here writes into `demo_analyze/`.** That
folder's tables and figures are the three-way result and stay as they are; this
folder is a parallel tree with its own `tables/` and `figures/`.

> **Two product trees, one set of scripts.** This one has the plate on the
> published Ensembl 99 reference. **`../fourway_analyze_e116/`** has the same
> figures with the plate re-mapped to Ensembl 116, so all four datasets share an
> annotation — that is the tree to quote biotype and gene-detection numbers
> from. Code root and product root are separate (`FOURWAY_OUT`, `PLATE_CELLS`);
> there is no second copy of any script.

## The four

| key | label | what a "unit" is | n | read len | reference |
|---|---|---|---|---|---|
| `own130` | VASA own, 130 nt | barcoded cell | 16 | 130 nt | GRCm39 **E116** |
| `own75` | VASA own, 75 nt | barcoded cell (same 16 barcodes) | 16 | 75 nt | GRCm39 **E116** |
| `plate` | VASA published, mouse | plate well, mouse call only | 173 | 74 nt | mixed GRCh38+GRCm38 **E99** |
| `fs` | FLASH-seq | **library**, not a cell | 10 | 151 nt | GRCm39 **E116** |

`fs` is the `native` arm only. The `vasalen` arm — FLASH-seq reads truncated to
VASA's length — is deliberately **not** a fifth dataset: length is already
controlled properly by direct standardisation (see below), and a truncated arm
would put a real protocol next to a synthetic one on the same axis.

### "Unit", not "cell"

VASA's units are barcoded cells inside one library. FLASH-seq's are ten
separate libraries in an input titration (30 ng → 30 pg, two replicates each),
each of which is **bulk RNA, not a cell**. Calling them all "cells" would make
the reads-per-unit figure read as a depth-per-cell comparison, which it is not.
They are units, and every figure says which is which.

## Four things that are not comparable, and what is done about each

This is the important section. A four-way bar chart makes everything look like
one axis; these are the places where it is not.

### 1. FLASH-seq had no rRNA-depletion stage. Fixed, by adding one.

`flashseq_vasa/pipeline_fs.sh` goes prep → map → assign. There is no step 3 in
it, so FLASH-seq's STAR input contained its rRNA and VASA's did not — VASA's
step 3 removes ~26% of its reads before STAR ever sees them.

Harmless for a count table (rRNA has no gene to be assigned to) and **fatal for
everything this folder measures**: rRNA is transcribed off a repeat array, so
rRNA reads multimap by construction. Leaving them in one library's STAR input
and not the other's puts a protocol difference and a pipeline difference on the
same axis — and multimapping rate and biotype composition are exactly the
demo_analyze conclusions that would corrupt.

`scripts/00_fs_ribo.sh` runs **VASA's own `ribo-bwamem.sh`** on FLASH-seq's
trimmed reads, same script, same rRNA reference, same two aligners. Map and
assign then run off the depleted fastq.

### 2. The strandedness flag is not the same on both sides

FLASH-seq is genuinely unstranded (49.1–50.5% of its ribosomal reads on the
forward strand, all ten libraries), so `stranded=n` is correct for it. VASA ran
`y`, and VASA is measurably not perfectly stranded either (76.1% forward). So
`y` discards ~24% of VASA's ribosomal reads where the same flag would discard
~50% of FLASH-seq's.

**The two in-silico depletion percentages are therefore not on one footing, and
this README says so wherever the number appears.** The flag is not split the
other way because the fastq that feeds STAR can only be made once. The
`.nsorted.all-ribo.bam` is kept so the `y` count stays recoverable from the same
alignment without re-running bwa.

### 3. Annotation release: E116 vs E99 — **now measured, in `../fourway_analyze_e116/`**

`plate` here is mapped to the published mixed E99 reference. E116 annotates
**32,889 lncRNA genes against E99's 9,959 (3.3×)**, while every other biotype
agrees within 8%, so in *this* folder a lncRNA comparison against `plate` mixes
annotation with biology, and is labelled as such rather than dropped.

This folder keeps E99 deliberately — it is what the paper cross-check needs.
But the confounder itself has since been measured, from both directions:

- **`scripts/00b_plate_e116.sh` re-maps the plate's 173 mouse wells to E116**,
  putting all four datasets on literally the same index and BED. Products in
  `../fourway_analyze_e116/` — read its README before quoting any biotype
  number.
- `flashseq_vasa/threeway/E99_MATCHED.md` does the opposite: our 16 cells
  re-mapped to E99.

They agree. Release explains roughly **half** the lncRNA gap — the plate moves
2.66% → 5.49% on E116, still 1.5× below our 8.3–8.4% — and about **15%** of the
structural gap. **The residual is real.** Neither direction alone shows that:
one of them was needed to prove the other was not the whole story.

### 4. FLASH-seq has no UMI

SmartSeq-like, `smartseq_noUMI`. Every UFI/UMI figure in `demo_analyze/` has no
FLASH-seq counterpart at all, and is **not** redrawn here (see below).

## Read length is controlled by standardisation, not by truncation

The four differ in read length by design (151 / 130 / 75 / 74 nt), and length
drives both multimapping rate and biotype composition. Each length-sensitive
figure therefore carries a **directly standardised** companion: one library's
per-length rates reweighted by another's read-length distribution — "what our
number would be if our reads were as long as theirs". The plate's distribution
is the reference weight, over the window where all four have reads.

A rate that survives standardisation is real. One that does not is length, and
`demo_analyze/README.md` already has a worked example of each.

## Three figures deliberately NOT reproduced

Not oversights:

- **read classes** — the class definitions come from VASA's three-pass trim;
  FLASH-seq trims in one cutadapt call, so there is no shared axis.
- **pass-2 loss attribution** — same reason; FLASH-seq has no pass 2.
- **UFI per unit** — FLASH-seq has no UMI.

The trim funnel is reported **per protocol**, never pooled, for the same reason.

## Depth is measured at step 3's input, on every side

`reads_in` = reads entering rRNA depletion = post-trim reads. That is the first
point where all four have run identically: before it, VASA has been
demultiplexed off a barcode read and FLASH-seq has not, so "demultiplexed reads"
has no FLASH-seq counterpart. Raw and post-trim counts are carried alongside
where they exist.

## A8 is kept

`FS_native_library_metadata.tsv` marks `ZHA8833A8` "exclude" for 18.3% human
CALB1 contamination and `A7` "caveat" at 3.6%. Both are processed and reported
identically to the rest and flagged in the `qc` column — the same rule this
project already applies to VASA's four blank cells. Dropping a unit because it
is odd is how an artefact stops being visible.

## Where the inputs come from — the part this folder does NOT do

**Read this before assuming `./run.sh` rebuilds anything.** Every script in
`scripts/` is a *reader*. They scan BAMs, BEDs and logs that four separate
upstream pipelines have already produced, and none of them can regenerate those
inputs. Rerunning the comparison from raw FASTQ means rerunning four different
entry points, in three different working directories, none of which live here.

That is a deliberate split — the four datasets are mapped by the pipelines that
own them, and duplicating those launchers here would be a second implementation
that drifts. But it does mean the chain has to be written down, and until
2026-08-06 it was not. It is below.

| dataset | raw input | entry point | output the scanners read |
|---|---|---|---|
| `own130` | `data/PM26037/*.fastq.gz` | `own_version/pipeline.sh` | `data/PM26037/out/cells/` |
| `own75` | `own130`'s step-1 output | `own_version/make75.sh`, then `pipeline.sh` | `data/PM26037/out75/cells/` |
| `plate` | `data/ref/fastq_vasaplate/SRR14783059_R{1,2}.fastq.gz` | `a_Mapping/submit_vasaplate_map_array.sh` | `data/ref/fastq_vasaplate/vasaplate_out_v3/` |
| `fs` | `data/flashseq_vasa/` raw | `flashseq_vasa/pipeline_fs.sh` + `scripts/00_fs_ribo.sh` | `data/flashseq_vasa/run/nonribo/cells/` |

### own130 — the user's own library

```bash
cd code/I_Gene_expression/own_version
./pipeline.sh check                 # verifies every path, tool and reference FIRST
sbatch -c 16 --mem=120G -t 24:00:00 --wrap="$PWD/pipeline.sh all"
```

`config.sh` is the only file meant to be edited. Per-step sbatch sizings are in
`pipeline.sh`'s own step headers (step 6 wants `-c 8 --mem=200G`). Steps 1–5 are
per cell; 6–7 run once over all of them. Full detail: `own_version/README.md`.

### own75 — the same 16 cells, truncated to the plate's read length

```bash
cd code/I_Gene_expression/own_version
./make75.sh                                        # 130 nt -> 75 nt, step 1 NOT re-run
OUTDIR=…/data/PM26037/out75 ./pipeline.sh step2    # …through step7
```

`make75.sh` truncates the **biological** read after step 1 has already stripped
the 21 nt technical prefix, so `-l 75` leaves 75 nt of insert — directly
comparable to the plate. Cutting the raw 151 nt R2 to 75 instead would leave 54
nt after the skip and would not be comparable. Step 1 is not re-run because
barcode extraction reads R1 and is unaffected by R2's length.

> `make75.sh` and `analyse75.sh` lived in `data/PM26037/out75/` — **outside the
> git repo**, since the repo root is `code/`. They were the only record of how
> one of the four datasets was built, sitting in the directory they write into.
> Copied into `own_version/` on 2026-08-06; the copies under `data/` are now
> historical duplicates and should be deleted so they cannot drift. This is the
> same defect that cost this project `build_mouse_reference.sh` — see
> "Reference provenance" in the root `CLAUDE.md`.

### plate — the published VASA-plate mixing control

The scanners read `vasaplate_out_v3`, which is a **fifth** run and is not in
`vasaplate_check/README.md`'s run table (that table stops at run 4,
`vasaplate_out_bedv2`, and only *plans* a v3 under a different name,
`vasaplate_out_rrnav3`). Its invocation was never written down. Reconstructed
2026-08-06 from Slurm's retained `SubmitLine` and the STAR command line in
`*_E99_Log.txt`:

```bash
cd data/ref/fastq_vasaplate
START=3 \
RIBOREF=$VASA_REFS/mixed/unique_rRNA_human_mouse.v3.fa \
REFBED=$VASA_REFS/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed \
${p2s}/submit_vasaplate_map_array.sh SRR14783059 MIXED 74 vasaplate_out_v3 vasaplate_out_v3 f
```

`START=3` — stages 1–2 were not rerun; the trimmed FASTQs come from an earlier
run. STAR index `mixed/star_index_74`, `stranded=y`, `--genomeLoad LoadAndKeep`.
Jobs 51029201 (ribo) → 51029414 (gmap) → 51029415/51029416 (b2bs/b2bm) →
51029417 (cout) → 51029418 (pick), 2026-07-29/30.

**Recover a lost invocation this way, not by guessing:**
```bash
sacct -j <jobid> --format=SubmitLine%3000 -X -n -P
```
This cluster retains `SubmitLine`, and STAR writes its full command line into
`Log.txt`. Both survived where the human record did not.

### fs — FLASH-seq

The one chain this folder does own a piece of, because `pipeline_fs.sh` was
missing an rRNA-depletion stage entirely (see §1 above). Commands below.

## Run it

Steps 0–4 read FASTQ/BAM/BED and are sbatch one-offs. Steps 5–6 are seconds on
the login node and are what `run.sh` drives.

```bash
FS=…/code/flashseq_vasa ; R=…/data/flashseq_vasa/run
E="FSV_OUTDIR=$R/nonribo FSV_ARM=native FSV_SAMPLE=FS_nonribo"

# 1. trim (native arm, unchanged). Writes to the lab share -- scratch gets purged.
sbatch --chdir=$FS -c 32 --mem=32G -t 240 \
  --wrap="FSV_SCRATCH=$R FSV_ARM=native FSV_NCORES=8 ./pipeline_fs.sh prep"

# 2. THE MISSING STAGE -- rRNA depletion. Resumable; FORCE=1 redoes all.
sbatch -c 32 --mem=48G -t 8:00:00 --wrap="NPAR=4 scripts/00_fs_ribo.sh"

# 3. map + assign off the DEPLETED reads.
#    FSV_SAMPLE is set explicitly: its default FS_$FSV_ARM would overwrite the
#    existing FS_native_* tables.
sbatch --chdir=$FS -c 16 --mem=64G -t 8:00:00 --wrap="$E ./pipeline_fs.sh map"
sbatch --chdir=$FS -c 10 --mem=32G -t 8:00:00 --wrap="$E FSV_NCORES=5 ./pipeline_fs.sh assign"

# 4. count tables. pickle1 is per library -- step 6's cost is linear in BED size.
for L in ZHA8833A1 … ; do sbatch --chdir=$FS --wrap="$E FSV_LIBS=$L ./pipeline_fs.sh pickle1"; done
sbatch --chdir=$FS --wrap="$E ./pipeline_fs.sh pickle_merge tables recon"

# 5. the four-way scans, which read all four datasets at once
scripts/01_insilico_depletion.sh                          # seconds, login node
sbatch -c 8  --mem=8G  -t 4:00:00 --wrap="scripts/02_probe_qc.sh"
sbatch -c 16 --mem=8G  -t 120    --wrap="scripts/03_mapped_length_dist.sh"
sbatch -c 16 --mem=64G -t 240    --wrap="scripts/04_step5_biotype.sh"

# 6. tables + every figure
./run.sh            # or: ./run.sh units | ./run.sh plots
```

`scripts/datasets.sh` is the manifest every scanner sources — one place that
knows what the four datasets are and where their files live. Its `ds_units`
**errors rather than returning nothing**: a glob that silently matches zero
files is the failure mode here, and it looks exactly like an empty facet.

## Layout

```
scripts/    numbered by run order; 00-04 shell (sbatch), 05-06 R
tables/     TSVs; tables/cross/ for the multi-dataset scans
figures/    01_reads/ 03_rrna/ 04_mapping/ 05_assign/ 06_coverage/   (02 = demo_analyze-only, see paths.R)
```

## Results

All four scans completed 2026-08-06 across all four datasets. Every number below
is from `tables/cross/`; read `insilico_depletion.tsv` and
`step5_biotype_length_standardised.tsv` before quoting any of it.

### rRNA — three different questions, three different answers

| dataset | in-silico depleted | flag | probe-target residual |
|---|---|---|---|
| VASA own, 130 nt | 26.09% | `y` | **16.68%** |
| VASA own, 75 nt | 26.44% | `y` | 16.78% |
| VASA published, mouse | 11.79% | `y` | **1.38%** |
| FLASH-seq | 4.76% | `n` | 1.95% |

The depletion column is **not** a like-for-like comparison — see §2, the flags
differ and `y` understates VASA. The probe column is, for the two VASA columns:
same protocol, same flag, same 47S contig, and **our own reaction left 12× more
probe-target rRNA than the published one** (16.68% vs 1.38%). That is a wet-lab
result about our RNase H step, not a pipeline artefact.

FLASH-seq's 1.95% is the no-reaction baseline, and it sits *below* our reacted
library. It is not evidence that skipping the probes is better: FLASH-seq is
polyA-primed, so rRNA is suppressed by a different mechanism, before any probe
exists. The right reading is that a working probe reaction (plate) and polyA
priming (FLASH-seq) land in the same place, and ours did not.

### Step 4 — multimapping

| dataset | to STAR | multi | unmapped |
|---|---|---|---|
| VASA own, 130 nt | 72,632,400 | **29.04%** | 12.43% |
| VASA own, 75 nt | 71,252,083 | **33.63%** | 10.06% |
| VASA published, mouse | 52,108,881 | 14.22% | 5.23% |
| FLASH-seq | 264,213,160 | 12.23% | 8.82% |

FLASH-seq multimaps at 12.23% — essentially the plate's rate, on 151 nt reads,
and less than half ours. Our excess is not a VASA-protocol property: the
published VASA plate sits with FLASH-seq, not with us.

### Step 5 — what the multimappers are, standardised to the plate's read length

`vs_plate` = standardised rate ÷ the plate's, so 1.0 means "same as published".

| biotype | own130 | own75 | plate | FLASH-seq |
|---|---|---|---|---|
| snRNA | 4.3× | **15.0×** | 1.0 | **0.0** |
| miRNA | 3.6× | 2.6× | 1.0 | 0.6× |
| snoRNA | 1.3× | 2.1× | 1.0 | **0.0** |
| MiscRna | 0.9× | 1.5× | 1.0 | 0.1× |
| ProcessedPseudogene | 0.2× | 0.3× | 1.0 | **1.3×** |

**FLASH-seq has essentially no snRNA and no snoRNA multimappers at all** — 0.0
against the plate's 1.0 and our 4–15×. That is the total-RNA/polyA distinction
appearing exactly where the biology says it must: snRNA and snoRNA are not
polyadenylated, so polyA priming never captures them. It is the cleanest
confirmation in this whole comparison that the four datasets are measuring what
they are supposed to measure.

With small RNA absent, FLASH-seq's multimapper pool is dominated by processed
pseudogenes (46.6% standardised, the highest of the four) — the same class that
dominates the plate. Our library's snRNA excess is what displaces it.

### Gene body coverage — nothing is 3'-biased, and that is the result

`figures/06_coverage/`, 11,377 protein-coding genes shared by all four.

| dataset | 5' 10% | middle | 3' 10% | 3'/5' | intronic (excluded) |
|---|---|---|---|---|---|
| VASA own, 130 nt | 0.0094 | 0.0104 | 0.0094 | **1.00** | 25.42% |
| VASA own, 75 nt | 0.0103 | 0.0102 | 0.0092 | 0.89 | 26.50% |
| VASA published, mouse | 0.0085 | 0.0105 | 0.0077 | 0.90 | **43.69%** |
| FLASH-seq | 0.0086 | 0.0108 | 0.0078 | 0.90 | **8.48%** |

All four sit within 0.89–1.00 on 3'/5' and hug the uniform line (0.0100) across
the body. **No protocol here shows the 3' pile-up that oligo-dT priming produces
on degraded input** — including FLASH-seq, which is the one that could have. The
input titration goes down to 30 pg and the curve does not tilt, so the RNA was
intact at every input level.

The one shape difference is at the 5' end: both VASA-own curves carry a bump
around 5–15% of transcript length that the plate and FLASH-seq do not. It is our
library's, not VASA's — the published plate does not have it.

**The intronic column is the real protocol signal, and it is not a caveat about
the curve — it is the finding.** FLASH-seq assigns 8.48% of its gene reads to
introns; VASA runs 25–44%. That is polyA priming versus total-RNA capture, and
it is the same distinction that shows up as FLASH-seq's zero snRNA/snoRNA
multimappers. Note the published plate (43.69%) captures *more* intronic signal
than our own library (25.42%) — worth understanding, since unspliced capture is
a headline VASA claim.

### Caveats that apply to the tables above

- Every `plate` comparison crosses an annotation release (E116 vs E99). lncRNA
  is uninterpretable for that reason (§3) and is omitted from the table above.
- `own75` is the honest headline row for length-sensitive claims: `own130`'s
  snRNA multimappers live above 75 nt, where the plate has no reads to weight.
- FLASH-seq units are libraries, not cells (§"Unit, not cell"), so the depth
  figures are not depth-per-cell comparisons.

## Where this stands (2026-08-06)

Done: prep (job 51252438, 5 min — the adapter rate reproduces TrimGalore's to
0.00 pt on all ten libraries), and rRNA depletion for 9 of 10 libraries.

In flight: `fsribo` retry for `ZHA8833A1` → `fsmap` → `fsassign`, chained on
`afterok`. Then steps 6–7, then the four scans, then `./run.sh`.

### The 2026-08-05 incident, and the two bugs it exposed

Three jobs on three nodes died within 100 s of each other at ~22:20. Ruled out
by measurement: memory (this cluster runs `jobacct_gather/cgroup`, so `MaxRSS`
is the **whole job** — peaks were 36/200 GB and 15/48 GB), walltime, disk, node
reboot, and `scancel`. Neither job reached its own error handler. The trigger is
not identifiable from an unprivileged account; both reruns went straight past
the failure points, so it was transient and external.

Two genuine bugs were hiding behind it, both in `00_fs_ribo.sh`, both fixed:

1. **`xargs -P` abandons all remaining input the moment one child dies by a
   signal.** Reproduced: 10 inputs at `-P2`, one child SIGTERMed → 3 ran, 7 were
   never launched, exit 125. A `|| true` on that line turned 125 into 0, so
   SLURM recorded **`COMPLETED` for a job that did 4 of 10 libraries**. Without
   the verify gate, four-tenths of FLASH-seq would have entered the comparison
   looking like all of it. Replaced with a job-control fan-out.

2. **Incomplete cleanup between attempts.** Three tools in this stage refuse to
   overwrite and fail instead — `gzip`, `samtools merge` ("Please apply '-f'"),
   and `samtools sort` (its `.tmp.NNNN.bam` shards are named deterministically).
   Only some were cleared, so leftovers from the interrupted run made every
   later attempt on `ZHA8833A1` fail on files it had itself written.

The lesson that generalises, and the reason the gate exists: **success is what
is on disk, never an exit code.**
