# FLASH-seq vs VASA-seq — rRNA content and RNA species detected

Status 2026-07-29. **Three of four comparison legs are complete and verified.** All 20
step-6 pickles finished, so the composition leg is now answered on all ten libraries;
only depth-matched detection remains unrun.
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

| biotype | native 151 nt | VASA-trimmed | fold (n=10) | n=6 | n=4 |
|---|---|---|---|---|---|
| ribozyme | 7 | 18 | **2.57×** | 2.20× | 2.50× |
| snRNA | — | — | **1.90×** | 1.87× | 2.00× |
| miRNA | — | — | **1.77×** | 1.86× | 1.92× |
| MiscRna | — | — | 1.71× | 1.64× | 1.63× |
| snoRNA | — | — | 1.69× | 1.70× | 1.64× |
| scaRNA | — | — | 1.67× | 1.50× | 1.38× |

**Correction to a stability claim — twice over.** At n=6 I wrote "every fold moved by
≤0.30"; at n=10 that is false, *ribozyme* moved **+0.37**. I then wrote that "the five
classes with double- or triple-digit counts all moved by ≤0.09" — **also false**:
*scaRNA* moved **+0.17** on 18→30 entries.

Recomputed from the two tables, the pattern is cleaner than either framing and tracks
**entry count**, not digit count:

| biotype | fold n=6 → n=10 | Δ | entries detected (native, n=10) |
|---|---|---|---|
| ribozyme | 2.20 → 2.57 | **+0.37** | 7 |
| scaRNA | 1.50 → 1.67 | **+0.17** | 18 |
| miRNA | 1.86 → 1.77 | −0.09 | 443 |
| MiscRna | 1.64 → 1.71 | +0.07 | 175 |
| snRNA | 1.87 → 1.90 | +0.03 | 138 |
| snoRNA | 1.70 → 1.69 | −0.01 | 287 |

**The two classes with fewer than ~20 detected entries moved; the four with more than
100 moved by ≤0.09.** So: quote the range as **1.7–2.6×**, rest the claim on *miRNA*,
*snoRNA*, *snRNA* and *MiscRna*, and treat **both** *ribozyme* and *scaRNA* as
order-of-magnitude only.

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

## 3. Composition on the non-rRNA denominator — COMPLETE

All 20 step-6 pickles finished 2026-07-29, so this leg is now answered on **all 10
FLASH-seq libraries against 12 VASA cells**, reads vs reads.

| biotype | FLASH-seq (trimmed) | VASA | gap (pp) | read-length effect (pp) |
|---|---|---|---|---|
| ProteinCoding | 81.77% | 64.12% | **+17.65** | −2.55 |
| lncRNA | 11.48% | 13.92% | −2.44 | +0.97 |
| **MiscRna** | 0.22% | **9.82%** | **−9.60** | −0.01 |
| **snRNA** | 0.0006% | **7.03%** | **−7.03** | 0.00 |
| snoRNA | 0.007% | 2.60% | −2.59 | 0.00 |
| ribozyme | 0.0009% | 1.08% | −1.08 | 0.00 |
| ProcessedPseudogene | 4.76% | 0.97% | **+3.79** | +1.44 |
| UnprocessedPseudogene | 0.88% | 0.09% | +0.79 | +0.01 |
| scaRNA | 0.0004% | 0.06% | −0.06 | 0.00 |

**Four gaps clear the ±3 pp read-length floor**, and they are the answer to the
question:

- **ProteinCoding +17.65 pp**, **ProcessedPseudogene +3.79 pp** — poly-A selection
  concentrates FLASH-seq on mRNA; shorter reads spread some of it onto paralogous
  pseudogenes.
- **MiscRna −9.60 pp**, **snRNA −7.03 pp** — VASA captures what FLASH-seq does not.
  The five structural-RNA classes together are **20.59% of VASA's non-rRNA reads
  against 0.23% of FLASH-seq's: 90×**, where read length moves those classes by
  ≤0.01 pp.

**This is the finding that survives the read-length confound.** The *short-species
detection* axis is confounded (§2 — length-matching recovers 1.7–2.6× more entries).
The *structural-RNA composition* gap is 90× and read length cannot touch it. Measured
on the same reference, annotation, counting code and denominator.

`snRNA`, `ribozyme` and `scaRNA` sit at **0.0004–0.0009%** of FLASH-seq's non-rRNA
reads — trace signal four orders of magnitude below VASA, **not zero**. Reported as
measured; a handful of reads is a different claim from none.

**VASA side uses the unfiltered `uniagg` tables**, not `analysis/`: the UMI-ceiling
filter drops 8 entries carrying **23.02% of uniagg reads** and shifts ProteinCoding by
**+19.19 pp**, and those entries were selected *for* saturating the ceiling — i.e. for
being the most abundant small RNAs. Correct for molecule counting, wrong here.

## 4. Depth-matched detection — still not run

The methodological prerequisite (the read-length control) is complete; the
downsampling pass to VASA's median depth was not run.


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
