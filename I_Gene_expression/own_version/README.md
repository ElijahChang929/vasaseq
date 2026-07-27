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
| `trim.sh` | Forked step 2 — three passes. See "Step 2" below. |
| `trim_bc_anchor.py` | Step 2 pass 0: cuts the 3' tail by finding the read's own barcode. |
| `step2_report.py` | Per-cell trimming table, written by step 2 itself. |
| `step3_report.py` | Per-cell rRNA table, written by step 3 itself. |
| `build_rrna_reference.sh` | Builds the rRNA fasta. See "Step 3" and "Reference". |
| `bc_PM26037_6nt.tsv` | Cell-barcode whitelist for this library (16 × 6 nt). |
| `trimtest/` | The step-2 trimming benchmark and its results. |
| `README.md` | This file. |

⚠️ `build_mouse_reference.sh` — which built the STAR index and the annotation
BED — **was deleted and is nowhere on the filesystem.** Both outputs still
exist and `./pipeline.sh check` passes, so nothing is blocked, but neither has
any provenance and neither can be rebuilt without writing the script again.
Its third output, the rRNA fasta, turned out to be **wrong** as well as
unreproducible; that one has been rebuilt and now has a tracked builder
(`build_rrna_reference.sh`). See "Step 3".

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
run, so the two analyses stay comparable. Three files, in two places:

| What | Where | Built by |
|---|---|---|
| STAR index | `reference/genomes/mus_musculus/GRCm39/star_index_151_r116/` | the nf-core rnaseq run (shared, not a second copy) |
| rRNA fasta | `reference/vasaseq/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa` | `build_rrna_reference.sh` |
| annotation BED | `reference/vasaseq/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.bed` — 718 272 rows | *(builder lost, see below)* |

`config.sh` points at all three; `./pipeline.sh check` verifies them.

The STAR index has `sjdbOverhang 150`, not 129. It is the FLASH-seq run's index,
reused rather than duplicated — same `genome.fa`, same Ensembl 116 GTF, same
562 855 junctions. An overhang **larger** than `readLength-1` is harmless (STAR
stores flanking sequence it never uses); one too small silently costs
junction-spanning sensitivity, which is why `pipeline.sh check` enforces `>=`
and not `==`.

> ⚠️ **`build_mouse_reference.sh` no longer exists.** It was deleted once its
> outputs were on disk, which left the STAR index and the BED with no record of
> how they were made — their provenance had to be reverse-engineered from a
> SLURM log. If either goes missing it has to be rebuilt by hand. **Keep
> reference builders in this repo.** `build_rrna_reference.sh` is the one that
> is done right: tracked here, idempotent, and it writes a
> `rrna_v2.provenance.txt` next to its output.

The human+mouse `mixed/` reference built for the published species-mixing
control is **the wrong reference here** — wrong species set *and* wrong read
length (`sjdbOverhang 73`).

---

## Quick start

```bash
# 0. once only: build the rRNA reference (~2 min). The STAR index and the BED
#    already exist -- see "Reference" above.
./build_rrna_reference.sh

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
| `step3` ribo | removes rRNA using **both** `bwa aln` and `bwa mem` — see "Step 3" | `....nonRibo.fastq.gz`, `logs/step3_report.txt` |
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

Step 4's throwaway files go to `$SCRATCH` (lab scratch by default), not
`$OUTDIR` — they are only the logs STAR insists on writing for its
`--genomeLoad LoadAndExit` / `Remove` calls. Nothing downstream reads
`$SCRATCH`; delete it whenever you like.

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
| + barcode anchor, in cutadapt | 78,171 | 65,746 | 84.1% | 40,602 | 51.9% |
| + 5' poly-T | 75,105 | 63,902 | 85.1% | 40,678 | 54.2% |
| **adopted: anchor in pass 0** | 73,369 | 63,926 | **87.1%** | **40,684** | **55.5%** |

+3.7% protein-coding exonic reads over upstream at the highest purity of any
setting tried, from the simplest of the settings tried. Cell 007 tracks cell 011 throughout. Splice junctions 18,700 →
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

### The barcode anchor — `trim_bc_anchor.py`

The tail never had to be recognised by its shape. Step 1 put the cell barcode
and the UMI on the read name, so for any given read the 12 nt after the poly-A
are a **known literal string**, not a pattern:

```
[insert][poly-A][ revcomp(CBC) revcomp(UMI) ][adapter ...]
                 \____ 12 nt, known per read ____/
```

Pass 0 finds it, drops everything from there to the 3' end, and walks back over
the poly-A. That is the whole algorithm. No adapter needed, no `min_overlap`,
no threshold to tune. Chance-hit probability is 4⁻¹², about 2 expected false
hits per 300,000 reads.

It fires on **29.6%** of reads, of which 79% match exactly and the rest within
one mismatch or as a partial anchor at the read end.

**This replaced a much more elaborate version, and the history is worth
keeping.** Rounds 7–9 did the anchor inside cutadapt, which cannot take a
per-read pattern. The UMI therefore became 6 wildcards; wildcards had to run at
`min_overlap` = the whole pattern length or they ate the end of every read; and
that forced 21 nt of adapter into the pattern, so it only fired when the read
had sequenced that far. Every one of those constraints was a tool limitation,
not biology. Head to head on cell 011:

| | anchor fires | unique | in annotation | **exonic (PC)** | purity |
|---|---|---|---|---|---|
| anchor inside cutadapt | 26.9% | 75,105 | 85.1% | 40,678 | 54.2% |
| **`trim_bc_anchor.py`** | **29.6%** | 73,369 | **87.1%** | **40,684** | **55.5%** |
| both together | — | 72,822 | 87.3% | 40,651 | 55.8% |

Equal exonic yield, higher purity, and far simpler. Running both adds nothing,
so the cutadapt pattern was deleted along with its `TRIM_ANCHOR_ADLEN`
parameter and the sweep that justified it.

Note what pass 0 reaches that the pattern could not: a read that **ends partway
through the barcode**. cutadapt needed the whole pattern inside the read; a
literal string can be matched as a prefix running off the 3' end.

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

### How short a read is still worth keeping (TRIM_MINLEN)

The floor matters more here than it would in most libraries, because the median
insert before the poly-A is 10 nt. Swept with `trimtest/bench_trim11.sh`, cell
011, protein-coding exonic reads per 300k:

| floor | unique | in annotation | **exonic (PC)** | purity |
|---|---|---|---|---|
| 12 | 76,258 | 87.1% | 41,015 | 53.8% |
| 15 | 75,909 | 87.1% | 40,938 | 53.9% |
| 18 | 74,504 | 87.1% | 40,784 | 54.7% |
| **20** | 73,369 | 87.1% | **40,684** | **55.5%** |
| 25 | 71,267 | 87.3% | 40,403 | 56.7% |
| 30 | 69,567 | 87.6% | 40,165 | 57.7% |

**Going below 20 is not worth it, and going above it is not either.** Dropping
the floor from 20 to 12 recovers 2,889 uniquely mapped reads and only 331 more
exonic ones — +0.8% yield for −1.7 points of purity. Raising it to 30 costs
1.3% of exonic reads to buy 2.2 points of purity. The curve is flat: 20 sits in
the middle of a shallow optimum, not on a cliff.

What the marginal reads are actually made of:

| band | extra unique | extra exonic | exonic share |
|---|---|---|---|
| 12–14 | 349 | 77 | 22.1% |
| 15–17 | 1,405 | 154 | 11.0% |
| 18–19 | 1,135 | 100 | 8.8% |
| 20–24 | 2,102 | 281 | 13.4% |
| 25–29 | 1,700 | 238 | 14.0% |

against a **55.5%** baseline for the library as a whole. Every band in this
range is 4–6× depleted in exons — these are not reads that were being unfairly
discarded, they are the same junk population the whole of step 2 is about.

Some arithmetic for why that is unsurprising: a random k-mer occurs ~2.7e9 / 4^k
times in a mouse genome — 0.63 times at k=16, 0.04 at k=18, 0.0025 at k=20. So
20 is roughly where uniqueness begins, before accounting for repeats, which
make it worse. STAR does not protect you either: `outFilterMatchNminOverLread`
is **relative** (0.66), so a 16 nt read needs only 11 matching bases to be
reported as an alignment.

⚠️ **`TRIM_MINLEN` did nothing below 20 until 2026-07-26.** `trim_galore`
hardcodes `$length_cutoff = 20` when `--length` is not passed, and `trim.sh`
did not pass it — so pass 1 had already deleted every short read before pass 2's
`-m` saw anything. `trim.sh` now passes `--length "$TRIM_MINLEN"` explicitly;
at the default 20 this is a verified no-op (md5-identical output). The same
applies to **upstream**, whose `cutadapt -m 15` is therefore dead code: its
effective floor has always been 20. `TRIM_MODE=legacy` reproduces that
faithfully, dead code included.

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
TRIM_ANCHOR_BC=no  ./pipeline.sh step2    # skip pass 0 (the barcode anchor)
TRIM_POLYT5=       ./pipeline.sh step2    # no 5' poly-T trim
TRIM_POLYA=        ./pipeline.sh step2    # disable the poly-A trim (don't)
```

Judge any change with `trimtest/annot_fraction.sh`, not with read counts.

Step 2's second pass needs **cutadapt ≥ 4.4** for per-adapter `;min_overlap=`.
The module tree stops at 1.18 (2018), so cutadapt 5.1 was pip-installed into
the `vasa` conda env and is called by absolute path (`TRIM_CUTADAPT`).
TrimGalore's own pass still drives the module's 1.18.

---

## Step 3: rRNA depletion, and why the reference was rebuilt

Step 3 maps every trimmed read against a small rRNA fasta with **both**
`bwa aln` and `bwa mem`, and keeps only the reads that neither aligner placed.
`riboread-selection.py` makes the call, and two of its properties are worth
knowing before you tune anything:

- **MAPQ is never consulted.** One hit of any quality is fatal. This is
  deliberate — rRNA-derived mapping artefacts otherwise pollute small-ncRNA
  quantification — but it means the *reference contents* decide everything.
- **Stranded (`STRANDED=y`) gives one reprieve:** a read whose only hits are on
  the reverse strand is kept, since VASA is stranded and a real rRNA read must
  be sense.

There is also a harmless off-by-one: a read group is only flushed when the
*next* read name appears, so the final group of each file is counted in the log
total but written nowhere. `in = ribo + kept + 1`, one read per cell.

### The reference was wrong, and it was wrong in a specific way

v1 (`unique_rRNA_mouse.fa`) was built by grepping the Ensembl 116 GTF for
`gene_biotype` in `{rRNA, Mt_rRNA}`. That cannot work, because **the rDNA repeat
array is collapsed in the GRCm39 primary assembly**. Ensembl 116 mouse has no
`Rn28s`, no `Rn45s` and no `Rn5-8s` gene at all — only dispersed 5S/5.8S copies
and fragments, median length 116 nt. The single 18S entry (`Rn18s-rs5`) is a
dispersed related-sequence copy, not the rDNA locus.

The paper did something different. Methods says the rRNA sequences came from
**NCBI**, and `../a_Mapping/README.md` names them: mouse `Rn45s, Rn6s, 12s, 16s,
47s`. That is a handful of *full-length pre-rRNA transcripts*, not an annotation
dump.

Measured, by tiling 130 nt reads across the true subunits and running them
through the real `ribo-bwamem.sh`:

| subunit | caught by v1 | caught by v2 |
|---|---|---|
| 18S | 27 / 27 | 27 / 27 |
| 5.8S | 1 / 1 | 1 / 1 |
| **28S** | **0 / 71** | **71 / 71** |
| **5'ETS + ITS1** | **1 / 60** | **60 / 60** |

28S is 4 730 of the 13 400 nt transcript and the most abundant rRNA by mass.

### v2 = v1 + the NCBI 47S unit

`build_rrna_reference.sh` (tracked here, idempotent) keeps all 356 Ensembl
sequences unchanged and adds one: `BK000964.3:1-13403`, the **transcribed**
portion of the mouse rDNA repeating unit. Two deliberate choices:

- **One sequence, not seven.** Reads straddling a subunit boundary (18S/ITS1,
  ITS2/28S) still align cleanly. QC decomposition is recovered by *position*
  instead, via `rrna_intervals.tsv`.
- **IGS excluded.** `BK000964.3` is 45 306 bp; everything past 13 403 is
  intergenic spacer, never transcribed, dense with SINE/LINE repeats. Including
  it would let any repeat-containing mRNA be deleted — and remember MAPQ is not
  checked.

**Adding 13.4 kb costs nothing in specificity.** 20 000 simulated
protein-coding-exon reads give **0 false positives under both v1 and v2**, at
130 nt *and* at 50 nt.

### Both aligners are load-bearing — don't "simplify" to one

`bwa mem` will not report an alignment scoring below 30 (`-T` default, match
score 1), so **a read under ~30 nt is invisible to it** however well it matches.
**20.2% of this library's trimmed reads are under 30 nt.** The aln-only group is
therefore large and short-read-dominated — on cell 001, 18 478 reads averaging
24.4 nt, 87% of them under 30 nt, which is 49% of that cell's entire rRNA
detection. Delete `bwa aln` and all of it leaks silently into step 4.

### Reading the output

`step3_report.py` runs automatically at the end of step 3 and writes
`logs/step3_report.txt`; re-run it by hand any time with
`step3_report.py $CELLDIR`. Table 1 is depletion plus the aligner split, table 2
is where the ribosomal reads landed within the 47S unit.

**Blank wells are processed and reported exactly like real ones** — same
reference, same thresholds, no special-casing. Their numbers *are* the control.

A high 5'ETS + ITS share is not automatically an error: those sequences exist
only in unprocessed pre-rRNA, and a total-RNA protocol is supposed to see them.
But check *where* in the 5'ETS the reads land before believing it — see the
poly-T leak below.

### Results, ZHA9292A1 (job 50788552, 2026-07-27 — post poly-T fix)

90 137 383 reads in, **19 282 729 ribosomal (21.39%)**, 70 854 638 kept. Full
tables in `logs/step3_report.txt`; this is the summary.

| cell | in | ribo | ribo% | | cell | in | ribo | ribo% |
|---|---|---|---|---|---|---|---|---|
| 001 ° | 284 597 | 29 868 | 10.49% | | 009 | 7 017 728 | 1 458 166 | 20.78% |
| 002 | 2 623 273 | 471 928 | 17.99% | | 010 | 13 279 354 | 2 855 786 | 21.51% |
| 003 | 4 565 233 | 1 145 347 | 25.09% | | 011 | 13 634 453 | 3 066 372 | 22.49% |
| 004 | 3 542 790 | 732 589 | 20.68% | | 012 | 11 082 329 | 2 487 380 | 22.44% |
| 005 | 3 045 052 | 629 745 | 20.68% | | 013 | 10 454 276 | 2 291 476 | 21.92% |
| 006 | 4 748 979 | 915 911 | 19.29% | | 014 ° | 335 982 | 19 683 | 5.86% |
| 007 | 9 192 416 | 2 019 038 | 21.96% | | 015 ° | 393 457 | 33 811 | 8.59% |
| 008 | 5 756 671 | 1 114 134 | 19.35% | | 016 ° | 180 793 | 11 495 | 6.36% |

° the four low-count wells.

The twelve real cells sit in a tight **18.0–25.1%** band, which is the reassuring
part — a per-cell rRNA fraction that varied wildly would mean the reference or
the barcodes were wrong. Composition is equally consistent across them: 28S
51.7–55.5%, 5'ETS 15.2–17.8%, mito 1.8–3.4%.

### The blanks used to look different, and it was a poly-T artefact — now fixed

**This is resolved; the section is kept because the diagnosis is the useful
part.** Before the fix, the low-count wells reported a *lower* ribo% but a
wildly different composition: 5'ETS **31–62%** against the real cells' 17–20%.
That was not extra pre-rRNA. Binning the 5'ETS hits by position:

- cell 016 (blank): **88.9% of them in one 200 nt window** (400–599)
- cell 011 (real): spread across the whole 4 kb, busiest bin only 10.5%

The reads in that window were **pure poly-T**, aligning to a T-rich stretch at
~414–424. Poly-T as a share of each cell's ribo calls, before the fix:

| | 016 | 014 | 015 | 001 | | 011 | 007 |
|---|---|---|---|---|---|---|---|
| ≥90% T | **52.0%** | 38.2% | 25.8% | 19.6% | | 1.1% | 1.0% |

Absolute counts were similar everywhere (7k–34k per cell), i.e. a constant
background that only *looked* large in a near-empty well.

**Root cause was in step 2, not step 3.** `-g "polyT5=T{20}"` is a fixed 20-mer,
not a variable-length run, and `-n 3` caps cutadapt at three passes — so at most
**3 × 20 = 60 nt** of poly-T could ever be removed. Verified on synthetic pure-T
reads: ≤70 nt are consumed and dropped, 100 nt leaves 40, 130 nt leaves 70. The
longest poly-T surviving in the real data was **exactly 70 nt** = 130 − 60, and
130 nt is this library's biological read length. 14.28% of cell 016's trimmed
reads were ≥90% T; poly-A was 0.00%, as it should be.

**Fix, measured against 2 000 real reads from cell 011:**

| | poly-T surviving | real reads kept | real bases kept |
|---|---|---|---|
| `-n 3` (old) | 2 of 9 | 1 997 | 198 535 |
| **`-n 10`** | **0** | 1 997 | **198 535** (identical) |
| `T{130}`, `-n 3` | 0 | 1 997 | 198 523 (12 bp over-trimmed) |

`-n 10` removes all of it and touches nothing else, so that is what `trim.sh`
now uses (commit `5c6d879`).

**Confirmed on the rerun** (job 50788552, the table above). The prediction was
that the real cells would barely move and the blanks' 5'ETS share would stop
dominating, and that is what happened:

| | before (50787728) | after (50788552) |
|---|---|---|
| blanks' ribo% | 8.7–12.1% | **5.9–10.5%** |
| blanks' 5'ETS share | **31–62%** | **13.5–20.8%** |
| real cells' 5'ETS share | 17–20% | 15.2–17.8% |
| whole-library ribo% | 21.53% | 21.39% |

The blanks now sit inside the real cells' 5'ETS band. Cell 001's aligner split
corroborates the mechanism independently: aln-only detection falls 18 478 →
10 682, because the reads that disappeared were exactly the short poly-T
population that only `bwa aln` could see.

---

## Step 5: assigning reads to genes, and two bugs found there

Step 5 turns "this read is at chr1:3,276,400" into "this read belongs to
Xkr4, exonic". `deal_with_singlemappers.sh` and `deal_with_multimappers.sh` run
off the **same** step-4 BAM and split it by how many places the read mapped —
unique reads get counted directly, multimappers go through step 6's rescue
hierarchy. Both are forked here; see "Where the fork differs" below.

### What each script does, in four stages

1. **Select reads, and move the CIGAR + mismatch count onto the read name.**
   `$1 = $1 ";CG:" $6 ";nM:" nm`. This is necessary because the next stage
   (`bamtobed`) flattens the SAM to 6 BED columns and the CIGAR would be lost.
   The read name is the only channel that survives bedtools — the same trick
   step 1 uses for the barcode and UMI.
2. **`bamtobed`, then `bedtools sort`.** Note there is **no `-split`**: a read
   spanning an intron becomes one interval covering the whole span, not two.
   That is why the CIGAR had to be preserved in the name.
3. **`bedtools intersect -wa -wb` against the annotation BED**, then a long awk.
   Output is one row per *(read, overlapping feature)* pair — a read touching
   three features gives three rows, deliberately; step 6 resolves them. The awk
   keeps only same-strand overlaps (`STRANDED=y`), tags each row with `jS:`
   (`IN` = read inside the feature, `OUT` = read spans past both ends →
   **discarded**, `5`/`3` = hangs off one end), and computes `x`, the read's
   normalised position along the *gene*.
4. **gzip.** The multi script additionally `sort -k4` (by read name) first, so
   step 6 sees each read's candidate loci contiguously. The single script does
   not need it; upstream's commented-out line says `# unnecessary`.

### The 9 output columns, and which ones step 6 actually reads

| # | name | read by step 6? |
|---|---|---|
| 1–3 | Chromosome / Start / End | no |
| 4 | Name (read name up to `;CG:`) | **yes** — UMI and cell id live here |
| 5 | Strand | **no** |
| 6 | Gene (`ENSMUSG..._Xkr4_ProteinCoding_exon`) | **yes** — split into Gene / Biotype / Label |
| 7 | Info (`CG:...;nM:...;jS:...`) | **yes** — `nM` and `jS` extracted |
| 8 | Length | **no** |
| 9 | Cov (`x`) | **no** |

Worth knowing before you trust a column: `Strand`, `Length` and `Cov` appear in
`countTables_2pickle_cellsSpliced.py` **only** in the `read_csv(names=[...])`
list and are never referenced again. `gene_assignment` uses just `nMs`, `jSs`,
`Biotype`, `Gene`, `Label` and `Name`. And `jS` is only ever tested as
`jSs == 'IN'` — `jS:5` and `jS:3` are never distinguished.

### Where the fork differs

Both bugs are upstream. `a_Mapping/` is left untouched (published code gets
comments, not logic changes); the fixed copies live here and `config.sh` points
at them via `ASSIGN_SINGLE_SH` / `ASSIGN_MULTI_SH`. Point those back at
`${VASA_SCRIPTS}` for upstream behaviour.

**1. `NH:i:10`–`NH:i:19` were silently dropped by both scripts.** Selection was
done by matching raw SAM text:

```awk
$0 ~ /NH:i:1\tHI:i:1\t/     # single
$0 ~ /NH:i:[2-9]/           # multi
```

`[2-9]` is a *single character*, so only the first digit after `NH:i:` is ever
inspected. `NH:i:10` fails the single test (`1` is followed by `0`, not a TAB)
and fails the multi test (`1` is not in `[2-9]`). `NH:i:20` passed the multi
test only by accident, `NH:i:2` being a prefix of it. STAR here runs
`--outFilterMultimapNmax 20`, so 10–19 is a legal range.

Measured on cell 011, first 3,000,000 alignments:

| NH | alignments | fate |
|---|---|---|
| 0 (unmapped, `--outSAMunmapped Within`) | 135,463 | correctly excluded |
| 1 | 899,120 | single ✓ |
| 2–9 | 1,680,905 | multi ✓ |
| **10–19** | **273,092** | **dropped by both** |
| 20 | 11,420 | multi, by accident |

≈21,400 reads, **~1.5% of all reads and ~5% of multimapping reads**. The fork
parses `NH` as a number (`nh==1`, `nh>=2`) instead of matching it as text.

**2. `if (readstrand="+")` is an assignment, not a comparison** — 4 places per
script. A single `=` in awk assigns and evaluates to the assigned value, so the
condition is always true, the `else` is dead code, and `readstrand` itself is
overwritten with `"+"`. The line exists because "which end does the read hang
off" has to be translated from genome terms (left/right) into gene terms
(5′/3′), and that translation depends on the strand:

```
                 |------------ exon -----------|
Lypla1 (+)       5' ══════════════════════════> 3'     left = 5' end
Xkr4   (−)       3' <══════════════════════════ 5'     left = 3' end
```

So a read hanging off the left edge should be `jS:5` on Lypla1 but `jS:3` on
Xkr4. With the `=` it is `jS:5` in both cases. Two knock-on effects: the
correctly-written `==` on the *next* line then computes `x` from the wrong end
of the gene, and the printed Strand column becomes `+`.

**This one does not change any count.** All three things it corrupts —
`jS:5`/`jS:3`, `Strand`, `Cov` — are among the columns step 6 never reads (see
the table above). It is fixed because it is wrong, and because anyone reaching
for `Cov` or `Strand` as a QC column would be misled. Reads tagged `jS:IN` or
`jS:OUT` never enter those branches at all.

### Verified head-to-head on real data

Both versions run against cell 001's step-4 BAM (job `50803918`), upstream vs
fork, same reference, same arguments:

**singlemappers — 50,752 rows both ways, and every column step 6 reads is
byte-identical.** The `nh==1` test is therefore exactly equivalent to the old
`NH:i:1\tHI:i:1\t` text match, as expected: fix 1 only ever concerned
multimappers. Differences are confined to the three dead columns, on 18.7% of
rows — those are the reads hanging off a feature edge on a minus-strand gene:

| column | rows differing |
|---|---|
| 5 `Strand` | 9,504 (18.7%) |
| 7 `Info` (the `jS:` part only) | 9,504 (18.7%) |
| 9 `Cov` | 9,500 (18.7%) |

**multimappers — 30,153 → 34,550 reads (+4,397), 83,658 → 110,290 rows
(+31.8%), and no read lost.** Every one of the 4,397 recovered reads has an
`NH` between 10 and 19 in the BAM — nothing outside that range appeared, which
is the point:

| NH | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 |
|---|---|---|---|---|---|---|---|---|---|---|
| reads | 481 | 345 | 330 | 340 | 1583 | 492 | 319 | 173 | 119 | 215 |

Reproduce with `$SCRATCH/step5cmp/cmp_step5.sh` (it copies the BAM into two
directories and runs each version in its own, since both write their outputs
next to the input).

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
   any single cell's reads. A trap frees it even if you Ctrl-C. The two
   bookkeeping calls (`LoadAndExit`, `Remove`) still write logs nobody reads;
   those go to `$SCRATCH`, not `$OUTDIR`. Note the trap's scope is the subshell
   `pipeline.sh` runs each step in, so `./pipeline.sh step4` releases the
   segment immediately while `./pipeline.sh all` holds it until step 7 exits.
4. **`pipeline.sh check` measures read length with `sed -n '2{p;q}'`, not
   `sed -n 2p`.** Without the `q` sed drains the whole stream, so `check`
   decompressed all 15 GB of each fastq just to read one line — minutes instead
   of milliseconds. Same output, and it matters at this file size.
5. **The rRNA reference was rebuilt** (`unique_rRNA_mouse.v2.fa`). The Ensembl
   biotype dump used first contained no 28S and no ETS/ITS, because the rDNA
   array is collapsed in GRCm39; the paper used full-length NCBI pre-rRNA
   instead. v2 restores that by adding `BK000964.3:1-13403`, at zero cost in
   specificity. See "Step 3" above. The two upstream scripts it feeds
   (`ribo-bwamem.sh`, `riboread-selection.py`) are unchanged — except for a
   pre-existing path bug in the former, fixed earlier (`${fq##*/}`, not
   `${fq#*/}`, which broke any absolute path). `riboread-selection.py`'s bare
   `gzip` is a real re-run hazard but is **deliberately left alone**; it is
   handled from this side instead. See "Counting files is not proof" below.
6. **`step2_report.py` / `step3_report.py`** — per-cell tables generated by the
   run itself, so no number in this README has to be reconstructed by hand.
7. **`deal_with_singlemappers.sh` / `deal_with_multimappers.sh`** — forked for
   two upstream bugs, both measured; see "Step 5" above. `NH:i:10`–`NH:i:19`
   matched neither script's text pattern and were dropped by both (~5% of this
   library's multimapping reads; the fork recovers 4,397 reads on cell 001,
   all of them in that NH range). And `if (readstrand="+")` is an assignment,
   not a comparison — fixed, though it changes no count, because the three
   columns it corrupts are ones step 6 never reads. `ASSIGN_SINGLE_SH` /
   `ASSIGN_MULTI_SH` in `config.sh` point back at `${VASA_SCRIPTS}` for upstream
   behaviour.
8. Everything else is the published pipeline, called unchanged.

---

## Where this stands (2026-07-27)

The 2026-07-26 chain (`50788551` step2 → `50788552` step3 → `50788553` step4)
**all reported COMPLETED, exit 0:0.** Steps 2 and 3 are good and their numbers
are in this README. **Step 4 was not** — see below — and was rerun.

| job | step | result |
|---|---|---|
| `50788551` | step2 trim | ✅ 21m38s, 193.99 M in → 90.14 M kept (46.5%) |
| `50788552` | step3 ribo | ✅ 21m46s, 21.39% ribosomal, poly-T fix confirmed |
| `50788553` | step4 map | ❌ **mapped stale input** — superseded |
| `50803001` | gzip repair | ✅ 4m47s, recompressed step 3's real output over the stale `.gz` |
| `50803002` | step4 map | ✅ 10m08s, MaxRSS 29.0 GB — **verified correct input** |

Step 4 is now good: for **all 16 cells** STAR's `Number of input reads` equals
`step3_report.txt`'s `kept` column exactly. Real cells map 54.6–62.5% uniquely
with 24.7–30.9% multi-mapping; the four blanks map 17.1–23.0% uniquely with
56.8–64.1% "too short", which is what a blank should look like.

`# sbatch -c 16 --mem=64G -N 1 -t 6:00:00 pipeline.sh step4` — measured on
`50803002`: 10m08s elapsed, MaxRSS 29.0 GB. Do **not** cut the request toward
30 GB on that number: MaxRSS does not fully account for the shared-memory
genome segment, which is charged to the cgroup.

### Counting files is not proof — the step-4 staleness incident

Step 4 read the **previous** run's `.nonRibo.fastq.gz`, not `50788552`'s. It was
caught by comparing STAR's `Number of input reads` against `step3_report.txt`:

| cell 001 | reads |
|---|---|
| new step-3 kept (`50788552`) | 254 728 |
| old step-3 kept (`50787728`) | 272 515 |
| **STAR input reads** | **272 514** ← the old run, minus the known 1-read flush |

All 16 `.nonRibo.fastq.gz` were timestamped 20:52–21:06 (the *old* run), while
step 3 ran 23:48–00:09 and left its real output sitting **uncompressed** as
`.nonRibo.fastq`, 21 GB of it.

**Cause:** `riboread-selection.py` ended with
`os.system('gzip '+output+'.nonRibo.fastq')`. Without `-f`, gzip refuses to
replace an existing `.gz` and exits 2; the return value was never checked, and
`do_ribo` sends the whole stage's stderr to `/dev/null`. Nothing appeared in any
log, `pipeline.sh status` said 16/16, and step 3's own "done: 16 nonRibo files"
line said 16 — because all three **counted files rather than checking
freshness**. That is the convention's "confirm a stage has really finished"
rule failing one layer deeper than it is usually stated: the job *had* finished,
and its output was still stale.

**This is a re-run hazard, not an error in the published analysis.** In the
upstream workflow every stage runs once into a clean directory, where `gzip` and
`gzip -f` are identical and nothing is affected. It only bites when a stage is
re-run over its own previous output — which `own_version` does constantly and
upstream did not. So the published numbers are not in question; the published
scripts are simply not re-run-safe.

That is why **`a_Mapping/` was left untouched**, per this repo's rule that
published scripts get comments and not logic changes. `riboread-selection.py`,
`concatenator.py`, `deal_with_*mappers.sh` and
`countTables_2pickle_cellsSpliced.py` all carry the same bare `gzip`; all four
stages are protected from this side instead, in `pipeline.sh`:

1. **Detection** — `step3_ribo()` no longer counts files. It checks, per cell,
   that the `.nonRibo.fastq.gz` is **newer than its trimmed input** and that no
   uncompressed `.nonRibo.fastq` was left beside it, and `return 1`s so an
   `afterok` chain stops instead of feeding step 4.
2. **Prevention** — `rm_stale` deletes a step's own previous outputs before it
   runs, so the bare `gzip` never meets an existing `.gz`. Wired into all four
   affected stages: step 1 (`concatenator.py`), step 3 (`riboread-selection.py`),
   step 5 (`deal_with_{single,multi}mappers.sh`) and step 6
   (`countTables_2pickle_cellsSpliced.py`). It removes **per cell**, not by
   wiping the directory, so `MAXCELLS` runs and partial failures only touch the
   cells actually being processed.

   Step 5 is the case that mattered most: step 6 *globs*
   `*.singlemappers_genes.bed.gz`, so a stale BED there would have been folded
   into the count tables with nothing to show for it.

   Verified against the real failure mode, not asserted: with a stand-in for the
   upstream scripts (write data, then bare `gzip`), a second run leaves the
   *first* run's data in the `.gz` and the second run's uncompressed beside it —
   reproducing the bug exactly — while the same sequence with `rm_stale` in
   front yields the new data and no leftover.

Step 3 did not have to be recomputed: all 16 uncompressed `.nonRibo.fastq`
files were verified complete (read counts matched `step3_report.txt` exactly,
no truncation), so the repair was just `gzip -f` over the stale `.gz`
(job `50803001`), followed by the step-4 rerun.

**Generalise it:** whenever a stage is re-run over a directory that already
holds its previous output, an output file existing proves nothing. Compare it
to the *input's* mtime, or to a count the stage itself reported.

### Next: step 5 → 6 → 7

Steps 1–4 are done and verified. Steps 6 and 7 need *every* cell finished, so
run step 5 for all 16 before starting step 6.

**The standing check after any mapping stage** — cheap, and it is the one that
caught the stale-input run:

```bash
# STAR's input reads must equal step3_report.txt's `kept` column, per cell
grep -H "Number of input reads" $CELLDIR/*_E99_Log.final.txt
cat $OUTDIR/logs/step3_report.txt
```

Do the equivalent after step 5 (row counts per cell BED against the BAM) and
after step 6 (cells in the pickle == cells on disk). `rm_stale` should make a
stale read impossible now, but the check costs seconds and does not rely on
believing that.

Sizing: step 5 has no `# sbatch` line yet — measure it on this run and fill it
in, the way steps 1–4 now are. Step 6 is the memory-hungry one (upstream's own
note asks for ~160 GB on a full plate; 16 cells should need far less, but
measure rather than assume).

### Known-open, not urgent

- The STAR index and the annotation BED still have **no build script** (see
  "Reference"). Only the rRNA fasta does.
- `../a_Mapping/riboread-selection.py` never flushes its last read group, so one
  read per cell is counted but never written. Harmless; documented in "Step 3"
  so nobody rediscovers it as a bug.
- The published `mixed/` rRNA reference still has the missing-28S defect this
  library's reference was fixed for. It was left alone deliberately — that run
  is a species-mixing control and does not depend on depletion depth.
