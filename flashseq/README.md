# FLASH-seq RN26038 — data quality

QC of the FLASH-seq input titration, for comparison against the VASA-seq library
(`PM26037`/`ZHA9292A1`). Code lives here; everything it writes lands in `res/flashseq/`,
which is **not** version-controlled.

## What the run is

`20260325_LH00442_0237_B23GT7GLT3` lane 7, NovaSeq X, 10 individually-indexed libraries
`ZHA8833A1..A10`, paired-end 151 bp, mouse. Mapped by nf-core/rnaseq 3.22.2 in `star_rsem`
mode against GRCm39 + Ensembl release-116 (`data/flashseq/nfcore_rnaseq_all.sh`), STAR index
built at `--sjdbOverhang 150`.

**It is an input-amount titration in duplicate, not ten comparable samples.** This mapping
comes from the RN26038 LIMS sheet — nothing in `data/flashseq/` carries it, and without it
the analysis is uninterpretable:

| libraries | input | LIMS names | wells (plate `eASF_A07136`) |
|---|---|---|---|
| A1, A2 | 30 ng | `30ng_1`, `30ng_2` | A:1, B:1 |
| A3, A4 | 3 ng | `3ng_1`, `3ng_2` | C:1, D:1 |
| A5, A6 | 1.5 ng | `1pt5ng_1`, `1pt5ng_2` | E:1, F:1 |
| A7, A8 | 60 pg | `60pg_1`, `60pg_2` | G:1, H:1 |
| A9, A10 | 30 pg | `30pg_1`, `30pg_2` | A:2, B:2 |

It is recorded in `sample_metadata.tsv`, which every script joins against. A mammalian cell
carries ~10–30 pg of total RNA, so **30 pg is the single-cell-equivalent rung** and is the one
that matters for the VASA comparison.

## Findings

### The titration is flat to 1.5 ng, then falls off a cliff

Replicate Pearson *r* of log2(TPM+1) over genes expressed in both, and genes at TPM > 1:

| input | *r* | genes TPM>1 | genes ≥1 count | genes ≥1, rarefied |
|---|---|---|---|---|
| 30 ng | 0.9866 | 14355 / 14380 | 25851 / 25646 | 25289 / 25545 |
| 3 ng | 0.9832 | 14386 / 14350 | 23744 / 23881 | 23713 / 23491 |
| 1.5 ng | 0.9799 | 14375 / 14387 | 22786 / 22661 | 22528 / 22737 |
| 60 pg | **0.8218** | 11797 / 9830 | 14739 / 14374 | 14842 / 14356 |
| 30 pg | **0.8435** | 10263 / 9903 | 13558 / 12860 | 13308 / 12913 |

The six ng-scale libraries are indistinguishable from each other. The gene-detection drop at
pg scale is **not** a depth artefact: the rarefied column is the expected detection count
after binomial thinning to the shallowest library's total (16.5 M), and the split survives it.

**The design's blind spot:** there is no rung between 1.5 ng and 60 pg — a 25-fold gap — and
that is exactly where the method stops working. If the aim was to locate the sensitivity
floor, the floor lies inside an interval the experiment never sampled. Adding ~500 pg and
~150 pg would close it.

### rRNA is ~5× higher than the pipeline reports

nf-core reports rRNA as the read fraction over Ensembl `gene_biotype "rRNA"` genes. On GRCm39
that annotation holds **354 genes which are essentially all 5S** (`n-R5s*`) plus one
`Rn18s-rs5` relic — there is **no `Rn45s`, `Rn28s` or `Rn5-8s` gene at all**, because the rDNA
array is collapsed out of the primary assembly. The reported figure measures 5S, not rRNA.
Check it yourself:

```bash
awk -F'\t' '$3=="gene" && $9~/gene_biotype "rRNA"/' \
    /nemo/lab/turnerj/working/guangxin/reference/genomes/mus_musculus/GRCm39/annotation/release-116/gtf/Mus_musculus.GRCm39.116.gtf \
  | grep -o 'gene_name "[^"]*"' | sort | uniq -c | sort -rn | head
```

This is the same defect found and fixed on the VASA side (see `CLAUDE.md`, "Reference
provenance"). Re-measured by exact 31-mer containment against `unique_rRNA_mouse.v2.fa`,
which carries the NCBI 47S unit `BK000964.3:1-13403`, on every 64th R1 read of each library
(362 k–605 k reads):

| | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 |
|---|---|---|---|---|---|---|---|---|---|---|
| nf-core biotype | 0.81 | 0.70 | 0.85 | 0.83 | 0.85 | 0.73 | 0.73 | 1.40 | 1.07 | 0.73 |
| **47S unit** | **6.17** | **4.90** | **4.33** | **4.43** | **3.97** | **3.49** | **3.16** | **4.58** | **4.58** | **3.79** |
| Ensembl records only | 0.15 | 0.19 | 0.24 | 0.19 | 0.21 | 0.24 | 0.23 | 0.30 | 0.12 | 0.08 |

**95.7 %** of the rRNA signal comes from the 47S unit Ensembl does not annotate. These are
exact matches sampled every 10 bases with no mismatch tolerance, so they are a **lower
bound**; `05_rrna_bwa.sh` closes that gap by aligning instead — see below.

### FLASH-seq carries 4.3× less rRNA than VASA, and the gap is in mature 28S

The number above is a lower bound from one method, and the VASA figure it needs to be
compared with came from another, so `05_rrna_bwa.sh` runs **the VASA pipeline's own rRNA
stage** over the FLASH-seq reads: `a_Mapping/ribo-bwamem.sh` and `riboread-selection.py`
unmodified, `bwa aln` **and** `bwa mem`, the same `unique_rRNA_mouse.v2.fa`.
`05_rrna_bwa_report.py` even imports `step3_report.py`'s own `parse_log`, so the arithmetic
is shared rather than re-derived. Trimmed arm, every 64th read:

| | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 |
|---|---|---|---|---|---|---|---|---|---|---|
| input | 30ng | 30ng | 3ng | 3ng | 1.5ng | 1.5ng | 60pg | 60pg | 30pg | 30pg |
| **bwa rRNA %** | **6.44** | **5.23** | **4.70** | **4.75** | **4.31** | **3.84** | **3.50** | **5.02** | **5.26** | **4.02** |
| k-mer lower bound | 6.33 | 5.09 | 4.57 | 4.63 | 4.18 | 3.73 | 3.40 | 4.88 | 4.69 | 3.87 |
| `stranded=y` would say | 3.23 | 2.66 | 2.36 | 2.37 | 2.13 | 1.95 | 1.73 | 2.50 | 2.62 | 1.98 |

**FLASH-seq 3.50–6.44 % against VASA-seq's 21.39 %** (`ZHA9292A1`, whole library;
18.0–25.1 % across its 12 real cells) — **4.3×**, measured by one script against one
reference. Five things worth keeping from it:

1. **The gap is in mature rRNA, not pre-rRNA.** As a share of *all* reads — which is the
   comparison that means anything, since a subunit can be a larger slice of a much smaller
   cake — **28S is 11.64 % of VASA reads and 1.43 % of FLASH-seq (8.1×)**, while 5'ETS is
   3.66 % vs 1.86 % (only 2.0×). That is what a poly-A-primed protocol against a total-RNA
   protocol predicts: mature rRNA is what oligo-dT priming misses.
2. **`stranded=n` was the right flag, and the cost of getting it wrong is a factor of two.**
   Re-selecting the same BAM with `y` — VASA's setting — halves every figure. The forward
   strand carries **49.1–50.5 %** of ribosomal reads in all ten libraries, which is the
   measurement that says the library is unstranded rather than the assumption.
3. **Trimming is irrelevant here: +0.01 points on average.** The worry that ~30 % adapter
   read-through was deflating the denominator was reasonable and turns out to be worth
   nothing, because those reads cannot map to rRNA either way.
4. **The k-mer screen was already within 3.5 % of the truth** (+0.16 points). Both screened
   the *same* reads — deterministic stride, identical counts, checked by the report — so
   that is a pure difference of method, and it is smaller than expected for an exact-match
   lower bound.
5. **`bwa aln` contributes essentially nothing at 151 nt** — 0 aln-only reads on the raw arm,
   42–92 on the trimmed. That is the exact mirror of VASA, where reads are short and aln-only
   was ~49 % of detection on cell 001. Both aligners are still run; the asymmetry *is* the
   reason the stage uses two.

Two cautions before quoting this. Composition varies more between FLASH-seq libraries than
within VASA (28S 24.3–44.9 % of ribosomal reads, against VASA's 51.7–55.5 % per cell), and
the two replicates of a rung differ more than the rungs do — A9 5.26 % vs A10 4.02 %, both
30 pg. And the 5'ETS share is high (32–47 %), which on the VASA side once turned out to be a
poly-T artefact; it is not one here — the report's `ets5pk` column shows the 5'ETS hits
spread across the whole 4 kb, busiest 200 nt window 12.3–16.0 %, against 88.9 % in the case
that *was* an artefact.

### The 60 pg pair is contaminated with human CALB1 — and it is a well effect

FastQC flagged one sequence as 15.9 % of A8 and could only say "No Hit". Tracked down, and
**A8 is 18.3 % human in total** once every fragment is counted across the whole file:

- R1 `TCTAGCCTGTGAGGGAACACGTTGAAGAAAAACTAAGCAGCAGGTAATGG` (5,060,926 reads) /
  R2 `TCCTCAGTTTCTATGAAGCCACTGTGGTCAGTATCATATTTTCTCCATGT` (5,029,903 reads); a second
  fragment adds 2.3 %.
- Absent from GRCm39 (full scan, both strands), from `unique_rRNA_mouse.v2.fa`, and from
  ERCC92.
- NCBI megablast: 100 % identity over 50 bp to human chr8q21 genomic clones (AF049895.3,
  AF068862.1, AC123779.5; *Pongo abelii* 98 %), **no hit in refseq_rna**.
- Located in GRCh38: all four sequences fall inside `ENSG00000104327` / ***CALB1***,
  chr8:90,058,608–90,095,475 (−). Fragment 1 pairs at 90,069,005 (exon) / 90,069,296
  (intron) — contiguous genomic, i.e. **unspliced**, which is why refseq_rna misses it.
  Fragment 2 pairs at 90,063,405 / 90,082,641, both exonic and 19 kb apart — a **spliced**
  transcript.

Counted in every library rather than trusting FastQC's ~0.1 % reporting floor:

| | A1 | A2 | A3 | A4 | A5 | A6 | **A7** | **A8** | A9 | A10 |
|---|---|---|---|---|---|---|---|---|---|---|
| input | 30ng | 30ng | 3ng | 3ng | 1.5ng | 1.5ng | **60pg** | **60pg** | 30pg | 30pg |
| well | A:1 | B:1 | C:1 | D:1 | E:1 | F:1 | **G:1** | **H:1** | A:2 | B:2 |
| CALB1 unspliced % | 0.48 | 0.35 | 0.32 | 0.06 | 0.64 | 0.24 | **3.54** | **15.90** | 0.05 | 0.00 |
| CALB1 spliced % | 0.01 | 0.01 | 0.01 | 0.00 | 0.00 | 0.01 | **0.06** | **2.30** | 0.00 | 0.00 |
| **all human, sum `pct_R1`** | 0.49 | 0.36 | 0.34 | 0.06 | 0.64 | 0.25 | **3.62** | **18.32** | 0.05 | 0.00 |

The 30 pg pair is *cleaner* than the 30 ng pair, so this is **not** "low input amplifies
background". A7 and A8 are wells **G:1 and H:1 — adjacent**; A9/A10 are A:2 and B:2. This is
localised contamination introduced during prep. The STP's own `fastq_screen`, which never saw
this pipeline, agrees: A8 is 79.9 % mouse / 42.4 % human against ~98 % / ~26 % for A1–A6.

### Adapter read-through is large, and trimming does not catch all of it

This is a **Nextera** library — TrimGalore auto-detected it
(`CTGTCTCTTATA`, the transposase mosaic end) and cutadapt reports adapter in **63 % of R1
reads**, removing 17 % of bases. So inserts are short and read-through is expected.

MultiQC understates how much: **FastQC only inspects the first 50 bp of a read**, while
read-through sits wherever the insert ends. Measured on raw reads against the
index-independent junction core `CTCGTGGGCTCGGAGATGTGTATAAG` (matched on either strand —
the reads carry its reverse complement), `probe_adapter_P7_readthrough` in
`res/flashseq/overrepresented.tsv`:

| | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 |
|---|---|---|---|---|---|---|---|---|---|---|
| R1 % | 19.6 | 15.7 | 24.0 | 28.5 | 28.1 | 28.3 | 24.3 | 23.2 | **31.9** | **32.1** |
| R2 % | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Strictly R1, and it rises as input falls (16–20 % at 30 ng, ~32 % at 30 pg) — shorter inserts
from less material.

**A subset survives trimming.** Comparing FastQC's overrepresented tables before and after
trimming — like-for-like, same method both sides (`04_trim_effect.py`):

| class | n read files | mean % raw | mean % trimmed | removed |
|---|---|---|---|---|
| adapter | 8 | 0.149 | 0.150 | ~0 % |
| polyG | 10 | 1.767 | 1.705 | 3.5 % |
| other (human CALB1) | 17 | 2.796 | 2.802 | ~0 % |

The reason the adapter ones survive is checkable, not guesswork: every one of those eight
sequences begins **mid**-mosaic-end (`CTTATACACATCT…`) and so does **not contain**
`CTGTCTCTTATA`, the pattern cutadapt anchors on. Reads that consist entirely of read-through
therefore pass through untrimmed, while reads with a partial 3′ adapter (the 63 %) are
trimmed normally. They are, however, a much smaller population than FastQC's table implies —
1.14 % of A10's R1, measured below.

**Poly-G is not removed at all.** MultiQC shows the exact 50×G read at 1.0–2.8 % of R2;
counting reads containing *any* 30 nt G-run gives R1: A1 3.3, A2 3.0, A3 4.8, A4 3.8, A5 5.3,
A6 5.0, A7 4.2, A8 3.5, A9 5.7, **A10 7.1 %**. After trimming it is still 2.71 % of A10 R2
(from 2.79 % raw). This is the NovaSeq X two-colour dark-cycle artefact — no signal is called
G — and TrimGalore's adapter pass has no concept of it.

Both therefore reach STAR and contribute to its "unmapped, too short" fraction.

#### …and the two artefacts turn out to be the same reads

Quantified by `06_trim_options.sh` on **A10** — the worst poly-G library, and a 30 pg rung —
over 90,404 read pairs. Three schemes, identical except for the flags under test:

| scheme | pairs kept | bases written | R1 still carrying mosaic end at its 3′ |
|---|---|---|---|
| current | 90,130 | 22,532,184 | 0.10 % |
| `+ --nextseq-trim=20` | 87,689 | 21,755,261 | 0.11 % |
| `+ also -a CTTATACACATCT` | 87,680 | 22,134,623 | **37.18 %** |

Three results, and two of them overturn what this section used to recommend:

1. **`--nextseq-trim=20` removes 2,441 pairs (2.7 %)** that currently reach STAR — matching
   the 2.80 % of R2 reads that are ≥ 80 % G almost exactly, so it takes the poly-G population
   and little else.
2. **The read-through and the poly-G are one artefact.** Only **1.14 %** of R1 actually begins
   mid-mosaic-end (1,032 reads — far less than FastQC's table implies, because FastQC reads
   only the first 50 bp), and **1,019 of those 1,032 have a poly-G R2 mate**. They are a
   no-insert pair seen from its two ends. `--nextseq-trim=20` already removes them; the extra
   adapter adds 9 pairs.
3. **The second adapter pattern is actively harmful — that proposal is withdrawn.**
   `CTTATACACATCT` is 13 nt against `CTGTCTCTTATA`'s 12, so on the ordinary read-through that
   the current settings trim correctly it wins cutadapt's best-match contest and cuts 7 nt
   *later*, leaving mosaic-end sequence on 37 % of R1 reads instead of 0.1 %. It is also why
   that row writes *more* bases than the row above it, which otherwise reads like a recovery.

So the fix is **`--nextseq-trim 20` / `--2colour 20` alone**. What it buys is a cleaner
denominator, not recovered data — those pairs have no insert to recover. Whether that
justifies a re-run is a judgement call, and putting a number on the mapping-rate gain needs
the STAR index this run did not save (see the caveats below).

**Reading `overrepresented.tsv`:** a fragment's R1 and R2 sequences are separate rows, and
`pct_in_library` is the larger of the two mates — so **summing that column double-counts
every pair** (it makes A8 look 36 % human instead of 18.3 %). Sum `pct_R1` instead; the
R2-mate rows contribute ~0 to it.

### The first N reads of a fastq are not a sample of it

Found 2026-07-27, while validating that `05_rrna_bwa.sh` reproduced TrimGalore. It does —
but the check failed anyway, and the reason was the sampling, not the trimming. One identical
cutadapt call on `ZHA8833A1`, reading the adapter-containing rate off the head of the file:

| first 20 k | first 400 k | first 2 M | every 64th | whole library |
|---|---|---|---|---|
| 58.5 % | 57.7 % | 56.4 % | **55.2 %** | **55.1 %** |

A fastq is ordered by flowcell position, and adapter content — hence insert length, hence
anything downstream of it — drifts along it. Stride sampling also reproduces the whole-library
quality-trimmed rate (0.7 %) and bases-written rate (89.2 %); head sampling reproduces
neither. The whole-library column is the saved TrimGalore report, so it is an independent
witness rather than another run of our own code.

`01_rrna_kmer_screen.py` and `02_contaminant_check.py` both took the first 400 k reads and
have been switched to `FS_STRIDE`. On the rRNA numbers the correction is small and one-signed
— **−2.1 % relative on average**, 9 of 10 libraries down:

| | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 |
|---|---|---|---|---|---|---|---|---|---|---|
| head 400 k | 6.43 | 5.24 | 4.74 | 4.77 | 4.25 | 3.84 | 3.50 | 4.97 | 4.73 | 3.87 |
| every 64th | 6.33 | 5.09 | 4.57 | 4.63 | 4.18 | 3.73 | 3.40 | 4.88 | 4.69 | 3.87 |

Nothing in the conclusions above moved — the CALB1 finding is a 25× difference between
libraries and could not care about 2 % — but rates are exactly the quantity this distorts, so
the tables now carry the stride numbers.

The stride is deterministic, so scripts sharing `FS_STRIDE` screen **the same reads** without
depending on one another. That is what makes the k-mer-vs-bwa gap below a difference of
method alone. The cost is one full decompress per file (~57 s), which is why `01` gained an
sbatch wrapper.

**Not fixable the same way: FastQC.** Its overrepresented-sequence table comes from a
head-of-file sample by design, and it only inspects the first 50 bp. That is why
`04_trim_effect.py` uses FastQC on *both* sides of trimming and reports only the difference —
the absolute levels there are not comparable with `02`'s whole-file rates.

### Otherwise the run is mechanically sound

Quality scores clean, zero reads flagged poor quality, strandedness auto-inferred
`unstranded` on all ten (correct for FLASH-seq), exonic 87.5–89.8 %, intronic 6.9–8.8 %,
intergenic 2.8–3.9 %, mitochondrial 1.1–3.2 %, 5′/3′ bias 1.14–1.31, duplication 43–59 %.

## Two caveats about the run's saved outputs

**`results/genome/Mus_musculus.GRCm39.116.filtered.gtf` is truncated to chromosome 1.**
393,637 lines, 4,742 gene ids, and a size of 167,837,696 B — an exact multiple of 64 KiB, the
signature of a partial copy. The sibling `.filtered.bed` (481,956 rows, all chromosomes) and
the RSEM reference (481,956 transcripts) are complete, so **this run's quantification is
unaffected** — but that GTF would silently produce a chr1-only analysis if reused. The
notebook asserts on this so it cannot be forgotten.

**The STAR index was not saved** despite `--save_reference`: `genome/index/` holds only
`rsem/` and `salmon/`. The work dir `/nemo/lab/turnerj/scratch/zhangg/flashseq` is already
deleted, so a re-run rebuilds it (22 min, 65 GB peak, per
`pipeline_info/execution_trace_2026-07-23_14-10-50.txt`, where `--sjdbOverhang 150` is
confirmed to have applied).

## What this means for the VASA comparison

**A8 is excluded — decided 2026-07-27.** It is recorded as `qc_verdict = exclude` in
`sample_metadata.tsv`, which every script joins against, and carried into
`qc_summary.tsv`, so the decision travels with the data rather than living only here.

| rung | libraries | verdict |
|---|---|---|
| 30 pg | A9, A10 | **`ok`** — the cleanest low-input pair, no detectable contamination |
| 60 pg | A7 | **`caveat`** — 3.6 % human CALB1; the reads are exactly identifiable and can be filtered |
| 60 pg | A8 | **`exclude`** — 18.3 % of the library is not mouse |
| ≥1.5 ng | A1–A6 | **`ok`**, but not input-comparable to VASA; use as the ceiling reference |

Two consequences of dropping A8, and they pull in opposite directions:

- **The 60 pg rung loses its replicate.** `replicate_concordance.tsv`'s 60 pg row
  (*r* = 0.8218) is a correlation *against an excluded library* and must not be quoted.
  There is no other 60 pg pair.
- **The titration conclusion does not depend on it.** The cliff between 1.5 ng and pg scale
  is shown independently by the **clean** 30 pg pair, *r* = 0.8435 against 0.980–0.987 for
  the ng rungs. So "flat to 1.5 ng, then falls off a cliff" survives A8's removal intact;
  what is lost is the ability to say anything about 60 pg specifically.

**A8 is still processed and reported everywhere, identically to the others.** The verdict is
a filter for interpretation, not permission to drop it from QC — its contamination is a
finding, and removing it would destroy the well effect (G:1/H:1 adjacent) that identified the
cause as prep rather than input amount.

**Quote FLASH-seq's rRNA fraction from `rrna_bwa.tsv`** — never the nf-core biotype number,
which measures 5S, and in preference to the k-mer lower bound. Only the bwa figures were
produced by the same script, the same two aligners and the same reference as the VASA ones;
everything else compares different quantities. On the 30 pg rung that is **A9 5.26 % / A10
4.02 %** against VASA's 21.39 %.

## Running it

```bash
cd /nemo/lab/turnerj/working/guangxin/vasaseq
source code/flashseq/config.sh
fs_check_python                              # confirms the interpreter

# on the login node -- these read nf-core's own output, not the FASTQs
fs_python code/flashseq/00_collect_qc.py     # seconds  -> qc_summary.tsv
fs_python code/flashseq/03_gene_detection.py # <1 min   -> gene_detection.tsv,
                                             #             replicate_concordance.tsv
fs_python code/flashseq/04_trim_effect.py    # seconds  -> trim_effect.tsv

# these three stream the FASTQs end to end -- submit them, do not run them here
sbatch code/flashseq/01_rrna_kmer_screen.sbatch    # -> rrna_kmer.tsv
sbatch code/flashseq/02_contaminant_check.sbatch   # -> overrepresented.tsv
sbatch code/flashseq/05_rrna_bwa.sbatch            # -> res/flashseq/rrna_bwa/
fs_python code/flashseq/05_rrna_bwa_report.py      # -> rrna_bwa.tsv

fs_python -m nbconvert --execute --to html \
    --output-dir res/flashseq code/flashseq/flashseq_qc.ipynb
```

`config.sh` is the only file meant to be edited — every path, `FS_STRIDE`, and the tool
locations live there. Every script is re-runnable and overwrites its own outputs; `05` clears
each library's directory before rebuilding it, because `riboread-selection.py` ends in a bare
`gzip` that refuses to clobber and would otherwise leave a rerun looking successful while
writing nothing.

**`05_rrna_bwa_report.py` refuses to build a partial table.** If any library's
`ribo-map.log` is missing it exits rather than tabulating what happens to be there — partial
output in that directory is indistinguishable from complete output otherwise.

**The notebook is generated, not hand-edited.** `flashseq_qc.ipynb` is committed without
outputs (a committed `.ipynb` full of base64 images diffs badly); its source of truth is
`build_notebook.py`. To change the report, edit that and regenerate:

```bash
fs_python code/flashseq/build_notebook.py code/flashseq/flashseq_qc.ipynb
```

The rendered copy *with* outputs is `res/flashseq/flashseq_qc.html`, outside git.

**Python:** the cluster's Anaconda3 base interpreter, used directly — it already has numpy,
pandas, matplotlib, scipy and the full nbconvert stack, so nothing is built or installed.
`fs_python` clears `PYTHONPATH` before calling it, which is **required**: if an EasyBuild
module is loaded in your shell, its packages shadow Anaconda's and an old `typing_extensions`
breaks `jupyter_client` on import. `config.sh` records the three environments tried before
this one and why each was rejected — in particular `envs/vasa` was deliberately left alone,
since it is the live VASA pipeline's environment.

**`05_rrna_bwa.sh` is the one exception, and it only borrows.** It runs `bwa`/`samtools` from
the EasyBuild modules and `riboread-selection.py` under `envs/vasa`, because that script needs
pysam and Anaconda base has none — which is exactly how the VASA side runs the same stage.
Nothing is installed into `envs/vasa`; it is activated and read. Its report script,
`05_rrna_bwa_report.py`, is back on `fs_python` like everything else here.

## Where this stands (end of 2026-07-27)

**Everything in this directory has been run end to end and every number above was produced
by the code next to it.** Nothing is queued, nothing is half-finished, no SLURM job is
outstanding. Outputs in `res/flashseq/`:

| file | written by | notes |
|---|---|---|
| `qc_summary.tsv` | `00_collect_qc.py` | 10 rows × 56 cols |
| `gene_detection.tsv`, `replicate_concordance.tsv` | `03_gene_detection.py` | |
| `rrna_kmer.tsv` | `01_rrna_kmer_screen.py` | job `50855262`, stride 64 |
| `overrepresented.tsv` | `02_contaminant_check.py` | job `50855331`, stride 64 |
| `trim_effect.tsv` | `04_trim_effect.py` | seconds |
| `rrna_bwa.tsv` + `rrna_bwa/` | `05_rrna_bwa.sh` → `05_rrna_bwa_report.py` | job `50855065`, 46 min, MaxRSS 778 MB |
| `trim_options/` | `06_trim_options.sh` | ~3 min, login node |
| `flashseq_qc.html` + `figures/*.pdf` | the notebook | 8 figures, 0 cell errors |

**`01` and `02` were re-run on 2026-07-27 after the head-of-file sampling defect was found**
(see the finding above); every number in this README is from those re-runs. `05` and `06` are
new in the same pass. **A8 was excluded the same day** — recorded as `qc_verdict` in
`sample_metadata.tsv` and carried into `qc_summary.tsv`, not only stated in prose. Re-running from scratch is three login-node commands plus three
`sbatch` — see "Running it" above.

### Three unit/method traps that were found the hard way

Recorded because each produced a plausible-looking wrong number first:

1. **MultiQC's `qualimap_*` columns are read counts in millions, not percentages.** For A1,
   `reads_aligned_exonic` = 39.46 against `reads_aligned` = 27.23 — impossible as a percent.
   Reading them as percentages gives "intergenic 1.7 %" instead of the correct 3.8 %.
2. **`multiqc_picard_dups.txt`'s `PERCENT_DUPLICATION` is a fraction** (0.458) while the
   general-stats copy is already scaled (45.8).
3. **Sequence searches must check both strands.** The P7 read-through core appears in these
   reads only as its reverse complement; a forward-only count scored it **0.00 % in all ten
   libraries** when the true figure is 19–36 %. `02_contaminant_check.py` now matches both
   orientations everywhere.

And two comparison traps:

**FastQC only inspects the first 50 bp of a read**, and draws its overrepresented-sequence
table from a head-of-file sample, so its percentages are not comparable with a whole-read,
whole-file search. Comparing the two is what briefly made trimming look effective on adapter
when it is not. Compare like with like — `04_trim_effect.py` uses FastQC on both sides.

**A longer adapter pattern is not a stricter one.** cutadapt picks the best-scoring adapter,
so adding a 13 nt pattern alongside a 12 nt one lets the longer pattern win on reads the
shorter one was trimming correctly — and cut 7 nt later. Adding what looked like a strictly
extra safeguard left mosaic-end sequence on 37 % of R1 instead of 0.1 % (`06_trim_options.sh`).
Measure a trimming change on reads it was *not* aimed at, not only on the ones it was.

## Next steps

Roughly in priority order.

1. **Start the actual FLASH-seq ↔ VASA comparison** once the VASA count tables are final
   (that is the other agent's `own_version` steps 2–7; step 6 was still running at
   2026-07-27 18:30, job `50836432`, expecting ~3 h from 16:01). The FLASH-seq side is ready:
   per library QC keyed by input amount, with A9/A10 (30 pg) as the single-cell-equivalent
   comparison point and A1–A6 as the not-RNA-limited ceiling; **A8 is excluded**. The rRNA
   leg of that comparison is already done and like-for-like — see the bwa section above.
2. **Decide whether a re-run with `--2colour 20` is worth it.** `06_trim_options.sh` has
   quantified the gain: 2.7 % fewer junk read pairs reaching STAR on A10, and nothing
   recovered, because those pairs have no insert. Putting a number on the mapping-rate change
   needs the STAR index this run did not save (22 min to rebuild). Do **not** add the second
   adapter pattern — measured harmful, see above.
3. **Fill the 1.5 ng → 60 pg gap** (~500 pg, ~150 pg) if locating the sensitivity floor is
   the actual experimental question — the current design brackets it but does not sample it.
   This is now the *only* way to get a 60 pg measurement at all, since excluding A8 leaves
   that rung with a single unreplicated library.
4. **Consider re-running nf-core with `--save_reference` fixed**, or at least regenerate the
   truncated `filtered.gtf`, if this run's saved reference will ever be reused. Rebuilding
   the STAR index is also the prerequisite for putting a number on item 2.

Not started, deliberately: the comparison itself (item 1), which is blocked on the VASA count
tables, and anything that would modify `data/flashseq/` — every pass so far only reads that
directory.
