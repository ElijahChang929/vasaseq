# own_version — minimal VASA pipeline for my data

A stripped-down, run-one-step-at-a-time version of the VASA-plate mapping
pipeline, for a dataset that carries a **21 nt uninformative prefix at the 5′
end of both mates**.

Currently configured for **PM26037 / ZHA9292A1** — see "This dataset" below.

Files here:

| File | What it is |
|---|---|
| `config.sh` | **The only file you normally edit.** All paths and settings. |
| `pipeline.sh` | The seven steps. Run one at a time, or `all`. |
| `concatenator.py` | Forked demultiplexer — the barcode geometry change. |
| `trim.sh` | Forked step 2. See "Step 2" below — this one is worth reading. |
| `bc_PM26037_6nt.tsv` | Cell-barcode whitelist for this library (16 × 6 nt). |
| `trimtest/` | The step-2 trimming benchmark and its results. |
| `README.md` | This file. |

⚠️ `build_mouse_reference.sh` is referenced below but **is not in this
directory and is not anywhere on the filesystem** — it was lost, not
committed, or never written. The three reference files it was supposed to
build do all exist and `./pipeline.sh check` passes, so nothing is blocked;
but if the reference ever has to be rebuilt, that script has to be written
again first.

Everything else (rRNA removal, gene assignment, count tables) is called **in
place** from `../a_Mapping/`. Those scripts are unchanged and already work, so
copying them here would only create a second thing to keep in sync.

---

## This dataset (PM26037, run 20260720_LH00442_0273_B23TM55LT4)

One library, `ZHA9292A1`, delivered as
`ZHA9292A1_S181_L007_R{1,2}_001.fastq.gz` (15 GB each, 2 × 151 nt). Everything
below was **measured from the reads**, not assumed.

| Property | Value | How it was established |
|---|---|---|
| 5′ prefix | **21 nt on both mates** | Per-position base composition is 98–99 % fixed for positions 1–21 and goes random at 22. R1 = `GAGTTCTACAGTCCGACGATC` (RA5 3′ end), R2 = `CCTTGGCACCCGAGAATTCCA` (revcomp of RA3). |
| UMI | **6 nt**, R1 pos 22–27 | All 4096 possible 6-mers observed; flat composition. |
| Cell barcode | **6 nt**, R1 pos 28–33 | Strongly structured; 16 sequences dominate. The 8-mer at 28–35 is always `<barcode>TT`, i.e. barcode + start of polyT — that is what rules out an 8 nt barcode. |
| polyT | from R1 pos 34 | 93 % T, decaying over ~24 nt. |
| Biological read | **R2, 130 nt** after the skip | 151 − 21. → STAR `sjdbOverhang = 129`. |
| Species | **mouse** | `fastq_screen` on R2: 7.8 % unique to MOUSE vs 0.6 % HUMAN. |

Note the barcode is **6 nt, not the CEL-seq2 8 nt** — `../a_Mapping/bc_celseq2.tsv`
(384 × 8 nt) does not apply to this library.

### Barcodes

16 barcodes, in `bc_PM26037_6nt.tsv`, numbered **alphabetically** so the ids are
stable and independent of which read subset you count. That is *not* abundance
order, so here is the map. Counts are from a 400 k read-pair slice:

| cell | barcode | reads | | cell | barcode | reads |
|---|---|---|---|---|---|---|
| 001 | ACTCGA | 2 169 | | 009 | CATGTC | 28 388 |
| 002 | AGACTC | 9 798 | | 010 | GTCTAG | 48 466 |
| 003 | AGCTAG | 19 376 | | 011 | GTGACA | 51 597 |
| 004 | AGCTCA | 15 559 | | 012 | GTTGCA | 43 829 |
| 005 | AGCTTC | 12 667 | | 013 | TCACAG | 38 058 |
| 006 | CAGATC | 19 953 | | 014 | TGCAGA | 2 845 |
| 007 | CATGAG | 36 566 | | 015 | TGTCAC | 3 185 |
| 008 | CATGCA | 23 455 | | 016 | TGTCGA | 1 532 |

**Cells 001, 014, 015 and 016 are 5–30× lower than the rest.** They are real
barcodes, not artefacts — each carries the same polyT anchor and none is a
1-nt-shifted neighbour of an abundant barcode — but check them against your
plate map before trusting them, and be ready to drop them at analysis time.
Ranks 17 and below in the raw data are 3× lower again and *are* shift
artefacts (`GTCTAG`→`TCTAGT`, `GTGACA`→`TGACAT`, …), so the whitelist stops
at 16.

`BC_HAMMING` is 0 and should stay there: several pairs differ only in their
last two bases (`AGCTAG`/`AGCTCA`, `CATGAG`/`CATGCA`, `TGTCAC`/`TGTCGA`), so
1-mismatch variants would be ambiguous.

### Reference

Mouse-only, **GRCm39 + Ensembl 116** — the same build as the nf-core rnaseq
run, so the two analyses stay comparable. Built once by
`sbatch build_mouse_reference.sh` into
`/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/`:

- `star_index_130/` — `sjdbOverhang 129`
- `unique_rRNA_mouse.fa` (+ bwa index) — 356 rRNA / Mt_rRNA gene seqs
- `Mus_musculus.GRCm39.116.homemade_IntronExonTrna.bed` — 718 272 rows

The human+mouse `mixed/` reference built for the published species-mixing
control is **the wrong reference here** — wrong species set *and* wrong read
length (`sjdbOverhang 73`).

---

## Quick start

```bash
# 0. once only: build the mouse reference (~1 h, mostly the STAR index)
sbatch build_mouse_reference.sh

# 1. config.sh is already set for ZHA9292A1 -- edit it only for a new library
# 2. verify everything exists BEFORE running anything
./pipeline.sh check

# 3. run step by step, looking at the output as you go
./pipeline.sh step1
./pipeline.sh step2
...
./pipeline.sh step7

# or all in one go
./pipeline.sh all

# how far along am I?
./pipeline.sh status
```

Every step writes to both the screen and `$OUTDIR/logs/stepN.log`, and every
step can be re-run — it overwrites its own outputs.

---

## The seven steps

Each step's output is the next step's input. That is the whole design.

| Step | Does | Produces |
|---|---|---|
| `step1` extract | strips the 21 nt prefix, moves cell barcode + UMI onto the **read name**, splits into one fastq per cell | `cells/<sample>_<cell>_cbc.fastq.gz` |
| `step2` trim | TrimGalore (adapters, Q<20) then cutadapt: 3' read-through, tail located by the cell barcode, poly-A/G/T — see "Step 2" | `..._cbc_trimmed_homoATCG.fq.gz` |
| `step3` ribo | removes rRNA using **both** `bwa aln` and `bwa mem` | `....nonRibo.fastq.gz` |
| `step4` map | STAR to the genome, multimappers kept | `..._E99_Aligned.out.bam` |
| `step5` assign | intersects reads with the annotation BED, unique + multi separately | `*_genes.bed.gz` |
| `step6` pickle | collapses **all** cells into one UMI-aware structure | `<sample>.pickle.gz` |
| `step7` tables | final spliced / unspliced / total / tRNA tables | `<sample>_*.tsv` |

Steps 1–5 are per cell. Steps 6 and 7 run **once over everything**, so steps
2–5 must have finished for every cell before you start step 6.

⚠️ The `_E99_` in the step 4 BAM name is an upstream filename token, kept so
the names match `../a_Mapping/`. It does **not** track the annotation release —
this reference is Ensembl **116**. Nothing reads it; don't infer from it.

---

## Debugging

**Work on a handful of cells first.** This is the single most useful trick —
a full run is hours, six cells is much less:

```bash
MAXCELLS=6 ./pipeline.sh step2
MAXCELLS=6 ./pipeline.sh step3
MAXCELLS=6 ./pipeline.sh step4
```

`MAXCELLS` affects steps 2–5. Step 1 always demultiplexes everything (it has
to read the whole fastq anyway). Note `MAXCELLS` takes the *first* cells in the
manifest, which here means the alphabetical barcodes — cell 001 is one of the
four low-abundance ones, so a `MAXCELLS=6` run is not representative of volume.

**Override any setting without editing `config.sh`:**

```bash
SKIP5=0 NCORES=4 OUTDIR=/somewhere/else ./pipeline.sh step1
```

**Make step 1 single-process** so it is easy to follow — set `NCORES=1` and it
takes a plain one-call path instead of the shard-and-merge one:

```bash
NCORES=1 ./pipeline.sh step1
```

**Run it on the cluster** rather than the login node. Steps 4 and 6 are heavy
(step 4 needs enough memory to hold the STAR index, step 6 is memory-hungry):

```bash
sbatch -p ncpu -c 16 --mem=120G -t 24:00:00 --wrap="$PWD/pipeline.sh all"
```

Sized for this library: 16 cores matches both `NCORES` and the 16 cells, and
the GRCm39 STAR index needs ~30 GB resident on top of step 6's appetite.

**Watch the disk.** Step 1 shards both fastqs to `$OUTDIR/.step1_work` before
merging, so it transiently needs roughly 3× the 30 GB of input on top of the
final per-cell output. Budget ~120 GB in `$OUTDIR`; it is cleaned up when
step 1 finishes.

---

## The one real thing to get right: read geometry

`./pipeline.sh check` prints this, and it is worth reading carefully:

```
first read: R1=151 nt, R2=151 nt
5' skip   : R1=21 nt, R2=21 nt
geometry OK: after skipping 21, the barcode block still fits in R1
-> build/verify your STAR index with sjdbOverhang = 129
```

After the 21 nt prefix is removed, R1 must still hold the barcode block:

```
R1 = [ 21 nt junk ][ UMI 6 ][ CBC 6 ][ polyT ... ]   needs >= 33 nt
R2 = [ 21 nt junk ][ real cDNA, 130 nt ]
```

If `check` says R1 is too short, then **the barcode is probably not behind the
prefix** and you want `SKIP5R1=0` with the prefix stripped from R2 only. You
can set that per mate:

```bash
SKIP5R1=0 SKIP5R2=21 ./pipeline.sh step1     # only R2 has the prefix
```

(`concatenator.py` takes `--skip5`, `--skip5r1`, `--skip5r2` directly if you
want to test it by hand.)

**The STAR index must match your read length.** `check` computes the right
`sjdbOverhang` for you — it is `(R2 length − SKIP5 − 1)`. An index built for a
different read length will still run and will still produce a BAM; it will just
map worse, silently. This is the easiest thing to get wrong.

### Sanity numbers

After step 1, look at `$OUTDIR/logs/step1_extract.summary`:

```
total sequenced reads: 40000
reads with proper barcodes: 36237 (0.9059)
```

**~0.9 is healthy.** If that fraction collapses toward zero, the barcode is not
where the config says it is — check `SKIP5`, `LEN_UMI`, `LEN_CBC`, `UMI_FIRST`
and that `BC_WHITELIST` is the right barcode list for your protocol.

For this library the settings were validated on a 400 k read-pair slice before
any full run: **357 443 / 400 000 = 0.8936** assigned, all 16 barcodes present,
and the emitted reads exactly 130 nt. That is in line with the published
VASA-plate library (0.8983), so the geometry is right.

---

## Step 2: what is trimmed and why

Changed 2026-07-26. The benchmark that produced every number below is in
[`trimtest/`](trimtest/) — scripts versioned here, inputs and outputs under
`../../../data/PM26037/trimtest/`. It is re-runnable with `sbatch`.

### The problem

Upstream's second trimming pass is

```
cutadapt -m 15 --trim-n -a polyG1=GG{5} -a polyC1=CC{5} -a polyT1=TT{5} -a polyA1=AA{5}
```

"Cut at a run of 6 identical bases and discard everything after it." As a 3'
adapter a 6-mer matches the **first** time it occurs, so a good 130 nt read
with an internal `AAAAAA` at position 10 is cut to 10 nt and then dropped by
`-m 15`. Measured on this library, restricted to reads carrying no poly-A tail
at all: **49% of them are truncated, median loss 120 of 130 bases**, and 30%
are cut below 20 nt.

Worse, it never removed what was actually in the way. The 3' read-through here
is

```
[insert][poly-A][12 nt = revcomp(CBC+UMI)][revcomp(R1 5' prefix) + Nextera]
```

TrimGalore auto-detects only the Nextera half, which starts 21 nt too late, so
every short-insert read went into STAR still carrying ~40 nt of adapter and
failed `outFilterMatchNminOverLread`. Naming that adapter is where nearly all
of the gain comes from.

The adapter was **measured, not assumed**: anchor reads on the 16 nt that is
unambiguously revcomp(R1's 5' prefix), take a per-position consensus of what
follows (68,250 reads, ≥97% agreement to position 55). A first guess of
`...GAACTCTGAAC` — extending the RA5 adapter — had five wrong bases; what
actually follows is the Nextera mosaic end. Re-derive this for a new library.

### The result

Per 300,000 input reads from cell 011 (a real cell). *poly-A-only* = uniquely
mapped reads whose **aligned** block (soft clips excluded) is ≥80% A or T, i.e.
the poly-A tail stuck to a genomic A-tract; *genuine* = unique − poly-A-only.

| setting | kept | unique | poly-A-only | **genuine** |
|---|---|---|---|---|
| `TRIM_MODE=legacy` (upstream) | 154,152 | 74,779 | 4.8% | 71,163 |
| adapter + poly-A (round 5) | 164,798 | 84,860 | 3.7% | 81,702 |
| **adopted (round 7/8)** | 135,539 | 74,808 | — | — |
| adapter only, no poly-A trim | 292,736 | 108,553 | 16.6% | 90,509 |

**Do not read that table as the answer** — uniquely mapped reads is the wrong
final score, and the adopted setting deliberately scores lower on it. See
"Judging by yield alone is wrong" below. The honest scoreboard is:

| setting | unique | in annotation | % | **protein-coding exonic** | % |
|---|---|---|---|---|---|
| upstream | 74,779 | 63,911 | 85.5% | 39,252 | 52.5% |
| adapter + poly-A | 84,860 | 69,750 | 82.2% | 40,968 | 48.3% |
| + barcode anchor | 78,171 | 65,746 | 84.1% | 40,602 | 51.9% |
| **adopted, + 5' poly-T** | 74,808 | 63,896 | 85.4% | **40,678** | **54.4%** |

+3.6% protein-coding exonic reads over upstream at the highest purity of any
setting tried. Cell 007 tracks cell 011 throughout. Splice junctions 18,700 →
19,400.

The third row is the trap. Dropping `TRIM_POLYA` looks like +45% unique reads,
but the extra alignments are poly-A landing on a handful of genomic A-tracts —
18,044 such reads over 3,637 100-kb bins with 1,054 in the top one. Keep the
poly-A trim.

### Judging by yield alone is wrong

Rounds 1–6 scored variants on uniquely mapped reads, on the grounds that a
variant cannot fake a bigger absolute count. That was too naive, and round 7
found the hole. Consider a read like

```
TCACATTCGA AAAAAAAA TGTCAC TTTGTA
|--insert--|-poly-A-|-rc(CBC)-rc(UMI)-|
```

10 nt of real sequence and 19 nt of junk. The round-5 setting kept it — the
poly-A run is only 8 nt so `A{20}` cannot fire — and STAR aligns the 29 nt
happily. It is not A-rich, so `aligned_composition.py` passes it too. The count
goes up and the data gets worse.

`annot_fraction.sh` is the test that sees it: junk has no reason to land inside
a gene. Comparing the round-5 setting against the adopted one, the 10,052 extra
uniquely mapped reads it brings contain **290 extra protein-coding exonic
reads** — 2.9%, against a 54% baseline. They are ~20× depleted in exons, i.e.
they are overwhelmingly intronic and intergenic. For VASA that is the worst
possible place to add noise, because the intronic signal *is* the measurement.

So the adopted setting scores **lower** on uniquely mapped reads than round 5
and is still the better one. Use `annot_fraction.sh`, not read counts.

### What is actually left in a trimmed read

Measured on the adopted setting, cell 011, 164,798 output reads:

| residue | reads |
|---|---|
| adapter still present | 0.26% |
| contains a poly-A run ≥20 nt | **0.00%** |
| ends in A≥6 | 1.88% |
| revcomp(cell barcode) in the last 25 nt | **13.58%** |

So the adapter and every long poly-A are gone, but the 12 nt
revcomp(CBC+UMI) remnant survives in 13.6% of reads. Understand *why*, because
it explains the whole design: `A{20}` is a **3' adapter**, so removing the
poly-A takes everything downstream — the 12 nt included — with it. That only
happens when there is a ≥20 nt run to match. Where the poly-A is shorter, there
is no match and both the short run and the 12 nt stay.

The barcode and UMI are already on the read name from step 1, so the remnant is
pure junk. Two ways of removing it by *shape* were tried and both cost more
than they gained (see the bullets below). What works is removing it by
**identity**:

### The barcode anchor

The tail does not have to be recognised by its shape — we know what follows it.
`concatenator.py` put the cell barcode and UMI on the read name at step 1, and
within one cell's fastq the barcode is **constant** (verified: a per-cell file
has exactly one distinct `CB` tag). So per cell the whole tail is a fixed
28 nt pattern with only the UMI unknown:

```
revcomp(CBC)  N x LEN_UMI  <first 16 nt of the read-through adapter>
 6 specific     6 any            16 specific        = 22 specific of 28
```

`trim.sh` reads the barcode off the first read's `CB:` tag — not the filename,
not the whitelist order — so it cannot silently pair a cell with the wrong
barcode. `min_overlap` is set to the full 28, which forbids partial 3'-end
matching; that is the difference between this and the round-2 wildcard disaster
below. With full-length matching required, 22 specific bases make a chance hit
impossible.

The point is what it unlocks. Once the anchor is gone the poly-A is at the very
3' **end** of the read, which is the only place `--poly-a` looks — so tails
*shorter* than 20 nt, which `A{20}` structurally cannot catch, get cleaned too.
Barcode remnant in the output: **13.6% → 2.2%**.

It does not increase yield — protein-coding exonic reads go 40,968 → 40,602,
i.e. flat. What it buys is that the reads which disappear are the ones that were
mostly barcode. That is the whole argument, and it is why the scoreboard above
had to change from "unique reads" to "exonic reads".

### 5' poly-T, and the reverse orientation

A read in the opposite orientation reads the poly-A tail as poly-T at its
*start*, so it is a 5' problem: `-g "T{20}"` removes the match and everything
before it. On a real cell it deletes ~3,400 uniquely mapped reads of which
~98% were non-exonic, leaving the protein-coding exonic count flat — it removes
junk and essentially nothing else. On the blank barcode 016 it is the
difference between 47,284 and 12,396 uniquely mapped reads: **the blank finally
looks like a blank.**

The mirror-image barcode anchor for reverse-orientation reads —
`revcomp(adapter) N{6} CBC` as a 5' adapter — was built and measured
(`trimtest/bench_trim8.sh`) and **fired 13 times in 300,000 reads**. This
library has no meaningful reverse population; all of the gain is the poly-T.
It is not in `trim.sh`; the benchmark is kept so nobody rebuilds it.

### What did *not* work, so you don't retry it

- **Just raising the homopolymer run 6 → 20.** Barely moves anything (73,867
  unique); the adapter is still there.
- **A wildcard adapter** `A{8}NNNNNNNNNNNN<adapter>` to describe the whole
  construct in one pattern. It *lost* ~8,000 unique reads. cutadapt lets a 3'
  adapter match partially at the read end and `N` matches any base, so the last
  ~20 bases of almost any read ending in a couple of A's get eaten for free.
  Wildcards and 3'-partial matching do not mix.
- **`--nextseq-trim=20`**: costs reads, gains nothing on this flow cell.
- **`--poly-a` instead of the strict `A{20}` adapter.** It yields *more* reads
  (95,439 unique, 86,276 genuine, against 84,878/81,717) but nearly triples the
  poly-A-only rate, 3.7% → 9.6%. The reason is mechanical: `--poly-a` trims a
  poly-A tail at the **very 3' end**, and after the read-through adapter is
  removed the 3' end is the 12 nt barcode remnant, not poly-A. It reaches past
  that only sometimes, so **40.5% of its output still contains a run of ≥20 A**
  against 0.00% for `A{20}`. Adding `--poly-a` *on top of* `A{20}` changes
  nothing (84,309 vs 84,860) — there is nothing left for it to do.
- **Chasing the 12 nt barcode remnant** with `A{10}` + 12 wildcards, at
  `min_overlap` equal to the full 22 nt so partial 3'-end matching is
  impossible. It does clean it (13.58% → 3.23% of reads) and it is *still* a
  net loss: 81,702 → 73,414 genuine reads, and the poly-A-only rate goes **up**,
  3.7% → 5.9%. Requiring a full match stops the wildcards eating read ends, but
  the pattern then fires on any internal `AAAAAAAAAA` with 12 bases after it and
  truncates real mRNA there — upstream's bug again with a longer run. You
  cannot tell a poly-A tail from an internal A-run once the anchor that
  distinguished them (the adapter) has been cut off.
- **Upgrading cutadapt on its own.** 1.18 and 5.1 produce **byte-identical**
  output for the upstream parameters (md5-checked, two cells). 5.1 is needed
  for per-adapter `;min_overlap=`, not for different behaviour.

### Reading the numbers, and what they say about the library

Trimming is not why the yield is what it is. **61.5% of reads carry a poly-A
run of ≥15 nt, and the median insert ahead of it is 10 nt** — only ~21% of all
reads have ≥15 nt of insert before the poly-A. Something like 40% of this
library is short-insert or empty product and no trimming setting recovers it.

The four low-count barcodes (001, 014, 015, 016) are worse than that, and the
trimmed data says so plainly rather than by read count alone. For cell 016 the
median **aligned** length is 39 nt against 123 nt for cell 011, and 53.5% of
its uniquely mapped reads are poly-A-only against 3.7%. Useful QC columns:
`aligned_composition.py` prints both.

### Knobs

All in `config.sh`; `TRIM_MODE=legacy` restores upstream byte-for-byte
(verified by md5 against the benchmark's `v0`).

```bash
TRIM_MODE=legacy   ./pipeline.sh step2    # upstream behaviour
TRIM_ANCHOR_BC=no  ./pipeline.sh step2    # no barcode anchor
TRIM_POLYT5=       ./pipeline.sh step2    # no 5' poly-T trim
TRIM_POLYA=        ./pipeline.sh step2    # disable the poly-A trim (don't)
```

Judge any change with `trimtest/annot_fraction.sh`, not with read counts.

Step 2's second pass needs **cutadapt ≥ 4.4** for per-adapter `;min_overlap=`.
The module tree stops at 1.18 (2018), so cutadapt 5.1 was pip-installed into
the `vasa` conda env and is called by absolute path (`TRIM_CUTADAPT`).
TrimGalore's own pass still drives the module's 1.18.

---

## What was actually changed vs. the published pipeline

1. **`concatenator.py`** — strips `--skip5` nt from both mates before parsing
   anything, then reads the barcode from the new start of R1. Adds an up-front
   geometry check that fails immediately (exit 1) on impossible settings, and
   records the skip lengths in the log. `--skip5 0` reproduces the original
   byte for byte (verified on 384 cells).
2. **`trim.sh`** — step 2's second pass trims the measured 3' read-through
   adapter, the tail located by the **cell barcode** rather than by shape, and
   real poly-A/poly-G/5'-poly-T runs, instead of cutting at the first run of 6
   identical bases. +3.6% protein-coding exonic reads at the highest purity of
   any setting tried; see "Step 2" above. `TRIM_MODE=legacy` restores upstream
   byte for byte (md5-verified on two cells).
3. **STAR is run with the genome loaded once** into shared memory for the whole
   run instead of reloading it per cell, because the index is far larger than
   any single cell's reads. A trap frees it even if you Ctrl-C.
4. **`pipeline.sh check` measures read length with `sed -n '2{p;q}'`, not
   `sed -n 2p`.** Without the `q` sed drains the whole stream, so `check`
   decompressed all 15 GB of each fastq just to read one line — minutes instead
   of milliseconds. Same output, and it matters at this file size.
5. Everything else is the published pipeline, called unchanged.
