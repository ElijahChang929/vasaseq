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
which carries the NCBI 47S unit `BK000964.3:1-13403`, on 400 k R1 reads per library:

| | A1 | A2 | A3 | A4 | A5 | A6 | A7 | A8 | A9 | A10 |
|---|---|---|---|---|---|---|---|---|---|---|
| nf-core biotype | 0.81 | 0.70 | 0.85 | 0.83 | 0.85 | 0.73 | 0.73 | 1.40 | 1.07 | 0.73 |
| **47S unit** | **6.26** | **5.03** | **4.49** | **4.57** | **4.04** | **3.60** | **3.25** | **4.66** | **4.62** | **3.77** |
| Ensembl records only | 0.17 | 0.21 | 0.24 | 0.20 | 0.20 | 0.24 | 0.25 | 0.31 | 0.11 | 0.09 |

**95.6 %** of the rRNA signal comes from the 47S unit Ensembl does not annotate. These are
exact matches sampled every 10 bases with no mismatch tolerance, so they are a **lower
bound** — see "Open" below for making them directly comparable to the VASA numbers.

### The 60 pg pair is contaminated with human CALB1 — and it is a well effect

FastQC flagged one sequence as 15.9 % of A8 and could only say "No Hit". Tracked down:

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
| CALB1 unspliced % | 0.43 | 0.29 | 0.30 | 0.06 | 0.56 | 0.22 | **3.25** | **14.67** | 0.04 | 0.00 |
| CALB1 spliced % | 0.01 | 0.01 | 0.01 | 0.00 | 0.00 | 0.01 | **0.07** | **1.84** | 0.00 | 0.00 |

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
| R1 % | 23.0 | 19.5 | 28.0 | 32.7 | 32.5 | 31.8 | 28.3 | 26.4 | **35.6** | **35.8** |
| R2 % | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Strictly R1, and it rises as input falls (19–23 % at 30 ng, ~36 % at 30 pg) — shorter inserts
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
trimmed normally.

**Poly-G is not removed at all.** MultiQC shows the exact 50×G read at 1.0–2.8 % of R2;
counting reads containing *any* 30 nt G-run on 400 k reads gives R1: A1 4.0, A2 3.7, A3 5.8,
A4 4.5, A5 6.2, A6 5.9, A7 5.2, A8 4.1, A9 6.5, **A10 8.2 %**. After trimming it is still
2.71 % of A10 R2 (from 2.79 % raw). This is the NovaSeq X two-colour dark-cycle artefact —
no signal is called G — and TrimGalore's adapter pass has no concept of it.

Both therefore reach STAR and contribute to its "unmapped, too short" fraction.
`--nextseq-trim 20` / `--2colour 20` addresses the poly-G; the untrimmed read-through needs a
second adapter pattern without the `CTGTCTCT` prefix. Quantify on one library before
deciding whether a re-run is justified.

**Reading `overrepresented.tsv`:** a fragment's R1 and R2 sequences are separate rows, and
`pct_in_library` is the larger of the two mates — so **summing that column double-counts
every pair** (it makes A8 look 33 % human instead of 16.6 %). Sum `pct_R1` instead; the
R2-mate rows contribute ~0 to it.

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

| rung | libraries | usable? |
|---|---|---|
| 30 pg | A9, A10 | **yes** — the cleanest low-input pair, no detectable contamination |
| 60 pg | A7, A8 | A7 with a caveat; **A8 no** — a sixth of the library is not mouse |
| ≥1.5 ng | A1–A6 | not input-comparable to VASA; use as the ceiling reference |

Quote FLASH-seq's rRNA fraction from the k-mer (or bwa) measurement, never the nf-core
biotype number — the VASA figures were computed against the 47S-containing reference, and
mixing methods compares different quantities.

## Running it

```bash
cd /nemo/lab/turnerj/working/guangxin/vasaseq
source code/flashseq/config.sh
fs_check_python                              # confirms the interpreter

fs_python code/flashseq/00_collect_qc.py     # seconds  -> qc_summary.tsv
fs_python code/flashseq/03_gene_detection.py # <1 min   -> gene_detection.tsv,
                                             #             replicate_concordance.tsv
fs_python code/flashseq/01_rrna_kmer_screen.py   # ~5 min -> rrna_kmer.tsv
sbatch    code/flashseq/02_contaminant_check.sbatch  # ~18 min -> overrepresented.tsv
fs_python code/flashseq/04_trim_effect.py    # seconds -> trim_effect.tsv

fs_python -m nbconvert --execute --to html \
    --output-dir res/flashseq code/flashseq/flashseq_qc.ipynb
```

`config.sh` is the only file meant to be edited — every path and the subsample size live
there. All four scripts are re-runnable and overwrite their own outputs.

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

## Where this stands (end of 2026-07-27)

**Everything in this directory has been run end to end and every number above was produced
by the code next to it.** Nothing is queued, nothing is half-finished, no SLURM job is
outstanding. Outputs in `res/flashseq/`:

| file | written by | notes |
|---|---|---|
| `qc_summary.tsv` | `00_collect_qc.py` | 10 rows × 56 cols |
| `gene_detection.tsv`, `replicate_concordance.tsv` | `03_gene_detection.py` | |
| `rrna_kmer.tsv` | `01_rrna_kmer_screen.py` | 400 k R1 reads/library |
| `overrepresented.tsv` | `02_contaminant_check.py` | job `50837920`, ~18 min |
| `trim_effect.tsv` | `04_trim_effect.py` | seconds |
| `flashseq_qc.html` + `figures/*.pdf` | the notebook | 7 figures, 0 cell errors |

Re-running from scratch is four commands plus one `sbatch` — see "Running it" above.

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

And one comparison trap: **FastQC only inspects the first 50 bp of a read**, so its
percentages are not comparable with a whole-read search. Comparing the two is what briefly
made trimming look effective on adapter when it is not. Compare like with like —
`04_trim_effect.py` uses FastQC on both sides.

## Next steps

Roughly in priority order.

1. **Re-measure rRNA with `bwa` so FLASH-seq and VASA are directly comparable.** Do not
   extend `01_rrna_kmer_screen.py` — run `code/I_Gene_expression/own_version/ribo-bwamem.sh`
   (bwa `aln` **and** `mem`, both, for the reason in `CLAUDE.md`) against these FASTQs
   against `unique_rRNA_mouse.v2.fa`. The k-mer figures are exact-match lower bounds and will
   move up. This is the one blocking item for a like-for-like FLASH-seq vs VASA rRNA claim.
2. **Decide A7/A8's fate.** A8 is not a usable 60 pg data point (16.6 % human CALB1). A7 is
   at 3.3 %; the reads are exactly identifiable, so filtering them is straightforward if you
   want to keep it. This is a judgement call, not a technical one.
3. **Start the actual FLASH-seq ↔ VASA comparison** once the VASA count tables are final
   (that is the other agent's `own_version` steps 2–7). The FLASH-seq side is ready: per
   library QC keyed by input amount, with A9/A10 (30 pg) as the single-cell-equivalent
   comparison point and A1–A6 as the not-RNA-limited ceiling.
4. **Quantify what better trimming would recover**, on one library, before deciding whether a
   re-run is justified: `--nextseq-trim 20` / `--2colour 20` for poly-G, plus an adapter
   pattern without the `CTGTCTCT` prefix for the read-through that currently survives.
   Both artefacts feed STAR's "unmapped, too short" (23–24 % in A8/A9).
5. **Fill the 1.5 ng → 60 pg gap** (~500 pg, ~150 pg) if locating the sensitivity floor is
   the actual experimental question — the current design brackets it but does not sample it.
6. **Consider re-running nf-core with `--save_reference` fixed**, or at least regenerate the
   truncated `filtered.gtf`, if this run's saved reference will ever be reused. Not needed
   for anything above.

Not started, deliberately: the comparison itself (item 3), and anything that would modify
`data/flashseq/` — this pass only reads that directory.
