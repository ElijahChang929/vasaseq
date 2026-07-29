# FLASH-seq vs VASA-seq — rRNA content and RNA species detected

Status 2026-07-29. **Two of four comparison legs are complete and verified; two are
blocked on compute that was still running when the cluster became unreachable.**
Everything below was measured on the user's own data. Nothing is quoted from a
publication. Where a number is provisional, it says so.

Datasets, both the user's own:

| | VASA-seq | FLASH-seq |
|---|---|---|
| sample | `ZHA9292A1`, run PM26037 | `ZHA8833A1–A10`, run RN26038 |
| unit | 16 plate barcodes = 12 real cells + 4 confirmed blanks | 10 bulk libraries, an input titration in duplicate |
| reads | 90.1 M (0.3–13.3 M per barcode) | ~17–23 M per library |
| chemistry | total RNA, fragmented + poly-A tailed, 6 nt UFI, stranded | poly-A primed, no UMI, unstranded, 151 nt PE |
| reference | GRCm39 + Ensembl 116, v2 annotation BED | same |

**The comparison point that matters is the 30 pg rung (A9, A10)** — roughly one
cell-equivalent of input RNA. A1–A6 (1.5–30 ng) are the not-RNA-limited ceiling.
**A8 is `qc_verdict=exclude`** (18.3% human *CALB1*, a well effect at H:1) but is
processed and reported everywhere; the verdict filters interpretation, not QC.

---

## 1. rRNA — COMPLETE

**FLASH-seq carries 4.5× less rRNA than this VASA plate.**

| | median | range | n |
|---|---|---|---|
| VASA own, real cells (`stranded=y`) | **21.14%** | 17.99–25.09% | 12 |
| FLASH-seq, A8 excluded (`stranded=n`) | **4.70%** | 3.50–6.44% | 9 |
| FLASH-seq 30 pg rung | 4.64% | 4.02–5.27% | 2 |
| VASA published `SRR14783059` (`stranded=y`) | 9.90% | 1.33–26.42% | 8 |

### Why this is one measurement and not two

`code/flashseq/05_rrna_bwa.sh` runs **VASA's own** `ribo-bwamem.sh` and
`riboread-selection.py`, unmodified, over FLASH-seq reads against the same
`unique_rRNA_mouse.v2.fa` (carrying the NCBI 47S unit `BK000964.3`). Both
aligners, `bwa aln` and `bwa mem`, as upstream runs them. That was deliberate —
`config.sh` warns that reimplementing it would silently create a second method.

### The strand flag, stated because it moves the answer twofold

**Each protocol is measured under its own correct flag: VASA `y`, FLASH-seq `n`.**

`stranded=y` counts a read as ribosomal only if it maps forward. That is correct
for VASA — the published code specifies it (`run_mapping_stepwise.sh:41`, and
`a_Mapping/README.md` twice), and `riboread-selection.py` is unmodified from
as-received. It is wrong for FLASH-seq, which is unstranded: its forward strand
carries **49.1–50.5%** of ribosomal reads, so `y` would discard half of its
genuine rRNA and penalise it for its chemistry.

**So 4.5× is a ratio of two correctly-flagged measurements, not a bare
division. Never quote it as the latter.** For sensitivity: both `n` gives **6.08×**,
both `y` gives **8.95×**. The asymmetry is real and measured — VASA's forward share is
**76.1%**, not the ~50% an unstranded library shows, so the same flag removes
23.9% of one library's ribosomal reads and 50.6% of the other's.

### Where the residual rRNA is

Median % of FLASH-seq ribosomal reads, 9 libraries:

| segment | median | range |
|---|---|---|
| **5'ETS** (pre-rRNA) | **41.1%** | 32.6–46.9% |
| 28S | 25.7% | 24.3–44.9% |
| 18S | 7.9% | 5.7–9.5% |
| other | 6.5% | 4.4–9.3% |
| ITS2 | 5.9% | 4.4–6.5% |
| mito | 4.0% | 1.9–6.6% |
| ITS1 | 3.8% | 2.9–4.4% |
| 5.8S | 1.9% | 1.5–2.2% |

**Most of what survives poly-A selection in FLASH-seq is pre-rRNA, not mature
rRNA.** That is mechanistically sensible: the 47S precursor is transient and
carries non-templated ends, so it is more poly-A-accessible than mature 28S. The
practical consequence is that the gap between the protocols is largest in mature
28S — where VASA's in-silico depletion has most to remove.

### Do not use the annotation route for rRNA

Ensembl's `rRNA` biotype reads **0.47–1.40%** where the alignment route reads
3.50–6.44%. On GRCm39 the rDNA array is collapsed out of the primary assembly,
so there is no `Rn45s`, `Rn28s` or `Rn5-8s` gene, and **99.93–99.999% of that
class's signal comes from a single 18S relic locus, `Rn18s-rs5`**
(`ENSMUSG00000119584`). Mechanism: RSEM's effective length — implied mean
fragment 195.1 nt, and 99.7% of rRNA features are shorter, so the one 1,849 bp
locus absorbs the whole class. **Label it "annotation-route 18S relic", never
"rRNA" and never "5S".** (This corrects `code/flashseq/README.md`, which called
it a 5S measurement.)

### Verification

`verify_rrna_crossplate.py` (job `0f559630`, exit 0) re-derives every cross-plate
number from the BAMs and asserts it — **15 of 15 pass**. Its load-bearing check is
that the `stranded=y` predicate reimplemented in it reproduces `step3.log`'s own
`ribo%` for all 16 barcodes, worst deviation **0.0048 percentage points**. As an
independent check, the read-weighted all-16 figure in
`rrna_comparison_summary.tsv` recomputes to **21.393%** against step 3's 21.39%.

---

## 2. RNA species — the read-length confound, and why it had to be measured first

**PARTIAL. The mechanism is established; the magnitude is provisional on 4
libraries.**

The naive result was that FLASH-seq detects **zero tRNA**. That is not a finding
about chemistry. VASA's assignment rule `jS:IN` requires a read to be *contained*
within the feature, and at 151 nt a FLASH-seq read cannot sit inside most small
RNAs: **98.5% of tRNA features are shorter than one read** (rRNA 99.7%, miRNA
99.2%, snoRNA 96.5%, snRNA 84.6%). In the A9 dry run six reads overlapped a tRNA
and **all six were discarded by containment**, while VASA cell 005 had 229 tRNA
rows. Read length alone produces "FLASH-seq detects fewer short species".

So both arms were quantified — identical in every setting **except read length**:

- `native` — adapter-trimmed R1 at natural length
- `vasalen` — additionally hard-trimmed to VASA's aligned-length distribution

### Result: length-matching recovers much of the gap

**6 libraries complete in both arms** (A1, A3, A4, A5, A8, A9 — including A9, the
30 pg rung). Entries detected:

| biotype | native 151 nt | VASA-trimmed | fold | fold at n=4 |
|---|---|---|---|---|
| ribozyme | 5 | 11 | **2.20×** | 2.50× |
| snRNA | 106 | 198 | **1.87×** | 2.00× |
| miRNA | 359 | 669 | **1.86×** | 1.92× |
| snoRNA | 246 | 417 | **1.70×** | 1.64× |
| MiscRna | 151 | 247 | 1.64× | 1.63× |
| scaRNA | 16 | 24 | 1.50× | 1.38× |

**The conclusion is stable under the n=4 → n=6 upgrade**: every fold moved by
≤0.30 and the ordering is essentially unchanged, which is the check that matters
for a provisional result. At n=4 the tRNA table went 27 → 46 rows (361 → 754
reads).

Note `MiscRna` recovers **entries** (1.64×) but not **reads** (0.96×): trimming
finds more distinct loci while splitting the same read mass across them. So
"recovery" means detection breadth, not extra signal, for the abundant classes.

Mechanism confirmed directly on the BED rows that fed step 6 (A9):

| | native | vasalen | VASA cell 005 |
|---|---|---|---|
| reads ≥140 nt | 82.4% | **55.4%** | 44.5% |
| contained in feature (`jS:IN`) | 45.7% | **52.0%** | 55.7% |

Trimming moves both metrics toward VASA, in the predicted direction.

### What this does and does not license

**Established:** the mechanism operates, and it is large. A 151 nt library
under-reports short RNA classes for reasons that have nothing to do with the
protocol's chemistry.

**Not established:** that read length explains the *whole* gap. Recovery is
partial — the trimmed arm still sits short of VASA on both metrics, and tRNA
remains absent from the biotype scan even though the dedicated tRNA table nearly
doubled. **The short-species axis is therefore confounded, not resolved.**

**Consequence for how to state the result:** "FLASH-seq detects fewer short
non-poly-A RNAs than VASA" is **not yet a claim this data supports**. What it
supports is: *at matched read length, FLASH-seq still detects fewer, but the
native-length difference overstates it by roughly 1.4–2.5× depending on class.*

---

## 3. Composition on the non-rRNA denominator — PARTIAL (6 of 10 libraries)

Measured on the 6 matched libraries, as % of non-rRNA reads (`rRNA` and `Mt_rRNA`
rows removed from the denominator, per the user's decision):

| biotype | native | VASA-trimmed | Δ pp |
|---|---|---|---|
| ProteinCoding | **84.15%** | 81.53% | −2.61 |
| lncRNA | 10.59% | 11.61% | +1.02 |
| ProcessedPseudogene | 3.35% | 4.81% | +1.45 |
| UnprocessedPseudogene | 0.88% | 0.88% | +0.01 |
| TranscribedProcessedPseudogene | 0.48% | 0.57% | +0.10 |
| MiscRna | 0.26% | 0.25% | −0.01 |
| MtRrna | 0.22% | 0.22% | −0.00 |
| miRNA | 0.042% | 0.080% | +0.04 |

**84% protein-coding is what a poly-A library should look like**, and the
comparison against VASA is the piece still missing. Note that read length alone
moves ProteinCoding by −2.6 pp and the pseudogene classes up: shorter reads are
less uniquely assignable, so some protein-coding mass redistributes to paralogous
pseudogenes. **Any FLASH-seq↔VASA composition difference smaller than ~3 pp is
within the read-length effect and must not be attributed to protocol.**

The remaining four libraries (A2, A6, A7, A10) still need step-6 pickles. **10 of 20 step-6 pickles are built** (all
`rc=0`, 39–59 GB, 7h40m–10h05m each); the remainder were still running when the
nemo login node became unreachable, and the sub-agent driving them had already
died on a network error (its jobs survived).

Two things are settled for when it resumes:

**Denominator (user decision):** rRNA as % of all reads; every other biotype as %
of the non-rRNA remainder, on both sides.

**Comparator (from the nf-core cross-check):** use the **unfiltered `uniagg`**
VASA tables, *not* `data/PM26037/out/analysis/`. The UMI-ceiling filter drops 8
entries — *Rmrp*, *Rn7sk*, *Rn7s1*/*Rn7s2*, *Rnu1a1*, *Rnu1b6*, *Rnu2.10*, the
*Snord3b* cluster, a *Cmss1* pair — that carry **23.02% of uniagg reads** and
whose removal shifts ProteinCoding by **+19.19 pp**. They were selected precisely
*for* saturating the UMI ceiling, i.e. for being the most abundant small RNAs,
which is what a total-RNA protocol exists to capture. Correct for molecule
counting; wrong here.

---

## 4. Depth-matched detection — BLOCKED

Same dependency. The methodological prerequisite (the read-length control) is
done; the downsampling pass was not run.

---

## The read-vs-molecule asymmetry — intrinsic, must always be stated

VASA counts **deduplicated molecules**; FLASH-seq via `smartseq_noUMI` counts
**reads**. This cannot be corrected away.

Measured: VASA read:UFI = **2.13–2.96 per cell**, median 2.57, pooled 2.69 — and
it is **expression-dependent** (Spearman +0.498 over 35,758 entries; 1.69 → 2.91
across expression deciles). Because the bias runs in **both directions depending
on the class**, no single factor fixes it.

**Therefore: compare composition on VASA `ReadCounts`** (reads vs reads), and
report VASA's molecule-based figures alongside as its own biology.

On the FLASH-seq side, **only `ReadCounts` is meaningful at all.** With one
literal UMI key `'A'`, `UFICounts` degenerates to a 0/1 detection mask (verified:
distinct values `[1]`), and `bc2trans` is the identity on a mask, so
`TranscriptCounts == UFICounts` elementwise. Nine TSVs of each are still written
and **carry no abundance information**. The UMI ceiling does not apply here
either — `K` is inferred from UMI length, giving `K=4`, not 4096.

---

## Reproducing this

```bash
cd /nemo/lab/turnerj/working/guangxin/vasaseq
# rRNA leg (both sides, one method)
code/flashseq_vasa/verify_rrna_crossplate.py \
    data/PM26037/out/cells \
    data/ref/fastq_vasaplate/rrna_validation/unique_rRNA_human_mouse.v2 \
    data/PM26037/out/logs/step3.log \
    res/flashseq_vasa
# FLASH-seq quantification, one arm at a time
FSV_ARM=native  code/flashseq_vasa/pipeline_fs.sh check prep map assign pickle1 pickle_merge tables recon
FSV_ARM=vasalen code/flashseq_vasa/pipeline_fs.sh check prep map assign pickle1 pickle_merge tables recon
```

`config.sh` is the only file meant to be edited.

**NOT REPRODUCIBLE NOW, CAUSE STILL UNKNOWN — do not treat as diagnosed.** The
`check` stage reported `MISS` for all 10 libraries in both arms on 2026-07-28
(`check_native.txt` 21:03:06, `check_vasalen.txt` 21:03:07).

What is actually established:

- **It does not reproduce.** Re-running `pipeline_fs.sh check` for both arms today
  gives `OK` for all 10 libraries in each, and reproducing `find_r1` step by step
  returns `n=1`, rc=0.
- **The `MISS` message names the same path `$FSV_FASTQ` resolves to now**
  (`.../RN26038/20260325_LH00442_0237_B23GT7GLT3/fastq`), so the checker was not
  looking somewhere else.
- **`find -L` handles the symlink correctly** — `$FSV_FASTQ` points into
  `/nemo/stp/sequencing/`, and `config.sh:87` documents the `-L` for exactly that
  reason.
- **No result depends on it**: the ten libraries produced trimmed FASTQs, BAMs,
  BEDs and pickles whose read counts reconcile, so the pipeline was not starved of
  data.

What is **not** established: *why* it failed. I asserted in an earlier version of
this file that the delivery filesystem was "briefly unreadable at 21:03" — **that
was inference presented as diagnosis, and I never tested it.** No probe in this
work measured filesystem availability at that time; the mount's `stat` shows only
its March mtime, and the `cutadapt_*.log` / `prep_*.tsv` files I first cited as
corroboration **do not exist at that path**. A transient mount failure remains the
most plausible explanation — the same window in which the driving agent lost its
connection — but it is a hypothesis.

**Three errors on the way to this, each worse than the last.** (1) I blamed the
glob for missing the lane field, when `ZHA8833A9_S*_R1_001.fastq.gz` matches
`ZHA8833A9_S108_L007_R1_001.fastq.gz` uniquely — and the `ls` I cited used a
*broader* pattern and never ran the checker's own glob, so my evidence contradicted
my claim. (2) I then called a 20/20 validator failure cosmetic. (3) I then named an
untested cause and marked the item resolved.

**Method notes:** to test a glob, run that glob. When a validator fails on every
input, suspect the environment before the validator — but *suspecting* is not
*diagnosing*. "Not reproducible" and "cause known" are different claims.

## Outputs

| file | what |
|---|---|
| `res/flashseq_vasa/rrna_comparison.tsv` | 34 rows: 10 FLASH-seq libraries, 16 VASA barcodes, 8 published cells, with the flag and its rationale on every row |
| `res/flashseq_vasa/rrna_comparison_summary.tsv` | subset medians, ranges, read-weighted figures |
| `res/flashseq_vasa/verify_rrna_crossplate.txt` | the 15 assertions and the predicate self-check |
| `shortbiotype_readlength.tsv` | per-biotype recovery, native vs length-matched |
| `res/flashseq_vasa/provenance.tsv` | component-by-component: same fasta, same scripts, which STAR index |
| `res/flashseq_vasa/read_vs_molecule*.tsv` | the asymmetry, per cell and per biotype |
| `res/flashseq_vasa/nfcore_composition.tsv` + 11 more | the nf-core cross-check |
| `code/flashseq_vasa/NOUMI_PATH.md` | what the `smartseq_noUMI` branch does, with line numbers |
| `code/flashseq_vasa/STEPS_2026-07-28_flashseq_vasa.md` | every step, and whether it made sense |
