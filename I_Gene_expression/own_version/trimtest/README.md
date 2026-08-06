# trimtest — the benchmark behind step 2's settings

Evidence for the step-2 change described in `../README.md` ("Step 2: what is
trimmed and why"). Re-runnable; each script is `sbatch`-able and skips work
that already exists, so a re-run after a partial failure is cheap.

**Scripts live here (versioned); inputs and outputs live in the data tree**, at
`/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/trimtest/`. Every
script `cd`s there itself, so run them from anywhere.

## Method

Every variant starts from the **same 300,000 reads** taken from the head of four
demultiplexed cells — 011 and 007 (high count, real cells) and 001 and 016 (two
of the four low-count barcodes) — and from the **same TrimGalore pass**. Only
the second, cutadapt pass varies. The score is therefore the **absolute number
of uniquely mapped reads**, not the mapping rate: a variant that discards more
reads can always show a better rate, but it cannot fake a bigger count.

That reasoning is not enough, and round 7 broke it. A read that is 10 nt of
insert followed by 19 nt of poly-A and barcode aligns as a 29 nt unit and
inflates the count with an alignment that is mostly junk;
`aligned_composition.py` does not catch it either, because such a read is not
A-rich. **`annot_fraction.sh` is the metric that decides**: junk has no reason
to land inside a gene, so the score is uniquely mapped reads overlapping a
protein-coding exon. The adopted setting scores *lower* on raw unique reads
than round 5 and is still the better one.

## Files

| file | what it does |
|---|---|
| `bench_trim.sh` | round 1: current setting, cutadapt 1.18 vs 5.1, longer homopolymer, a first "model the construct" attempt, and a no-homopolymer control |
| `bench_trim2.sh` | round 2: build up from the control one modifier at a time |
| `bench_trim3.sh` | round 3: literal adapters with a high `min_overlap`, after wildcards proved harmful |
| `bench_trim4.sh` | round 4: same, with the read-through adapter **measured** from the reads |
| `bench_trim5.sh` | round 5: is poly-G free? (yes). **v17 = the adopted setting** |
| `bench_trim6.sh` | round 6: can the leftover 12 nt barcode remnant be removed *by shape*? (no) |
| `bench_trim7.sh` | round 7: locate the tail by the **cell barcode** instead. **v23** |
| `bench_trim8.sh` | round 8: reverse orientation (absent) and 5' poly-T (worth it). **v25** |
| `bench_trim9.sh` | round 9: sweep how much adapter the cutadapt anchor carries |
| `bench_trim10.sh` | round 10: do the anchor in python instead. **v26 = adopted** |
| `bench_trim11.sh` | round 11: sweep the minimum-length floor (20 is a shallow optimum) |
| `bench_trim12.sh` | round 12: `--nextseq-trim` vs the poly-G adapter, against **current** production. **Rejected** |
| `bench_trim13.sh` | round 13: does pass 0 need its own poly-A trimmer? **No — remove `strip_polya`** |
| `annot_fraction.sh` | where alignments land in the annotation — **the metric that decides** |
| `bench_bam.sh` | re-maps selected variants keeping the BAM |
| `aligned_composition.py` | poly-A-only / short / soft-clip breakdown from those BAMs |
| `summarise_trim.sh` | the comparison table; `VARS="v1 v13 v17" ./summarise_trim.sh` |

## Variants

| | second pass |
|---|---|
| v0/v1 | upstream: `-m 15 -a XX{5}` on A/C/G/T (v0 = cutadapt 1.18, v1 = 5.1 — **byte-identical**) |
| v2 | as v0 but the run raised 6 → 20 |
| v3 | wildcard construct + `--poly-a` + `--nextseq-trim` |
| v4 | control: TrimGalore only |
| v5–v8 | v4 plus, one at a time: wildcard adapter, `--poly-a`, `--nextseq-trim`, poly-G |
| v9–v12 | literal adapters, `min_overlap=8`; adapter guessed |
| v13–v16 | as v9–v12 with the **measured** adapter |
| **v17** | **v14 + poly-G — what `trim.sh` runs in `TRIM_MODE=vasa`** |
| v18–v20 | v17 plus attempts to remove the 12 nt remnant by shape; all net losses |
| v21–v23 | barcode-anchored tail removal; **v23** = anchor + v17 + `--poly-a` |
| v24 | v23 + reverse-orientation anchor + 5' poly-T |
| v25 | v24 minus the reverse anchor (it fired 13 times) |
| L8…L26 | v25 with the anchor carrying 8/12/16/21/26 nt of adapter |
| L21 | = v25 at the principled boundary; best the cutadapt anchor can do |
| **v26** | **`../trim_bc_anchor.py` pass 0 + cutadapt without the anchor pattern — what `trim.sh` runs** |
| M12…M30 | v26 at minimum-length floors of 12/15/18/20/25/30 (**M20 = adopted**) |
| v27 | v26 keeping the cutadapt anchor as well — adds nothing, hence deleted |
| G0–G3 | round 12: poly-G by `--nextseq-trim` instead of the `G{20}` adapter — **all rejected**, see below |

## Round 12 — `--nextseq-trim` is the documented tool, and it is still wrong here

cutadapt documents `--nextseq-trim` as the way to handle two-colour chemistry,
where "basecalls without any signal are called as high-quality G bases", and
this **is** a two-colour run (flowcell LH00442 = NovaSeq X Plus). The artefact
is real and `-q 20` provably cannot reach it: of the reads in cell 011 ending in
a run of ≥10 G, **75% have mean Phred ≥ 25**.

It was still rejected, because the reads it removes are mostly not that artefact.

| variant | second pass |
|---|---|
| G0 | current production, unchanged — the baseline |
| G1 | G0 + `--nextseq-trim=20` |
| G2 | G1 with `-a polyG` removed |
| G3 | pass 1 `trim_galore --2colour 20` instead of `-q 20`, pass 2 unchanged |

cell 011 (real) / cell 016 (blank), uniquely mapped → in annotation → protein-coding exonic:

| | uniq 011 | exonPC 011 | uniq 016 | exonPC 016 |
|---|---:|---:|---:|---:|
| **G0** | **73,304** | **40,680 (55.5%)** | **13,930** | **3,190 (22.9%)** |
| G1 | 73,094 | 40,632 (55.6%) | 13,497 | 3,167 (23.5%) |
| G2 | 73,097 | 40,634 (55.6%) | 13,503 | 3,169 (23.5%) |
| G3 | 73,227 | 40,655 (55.5%) | 13,791 | 3,160 (22.9%) |

**Every variant loses protein-coding exonic reads** (−48, −46, −25 on cell 011)
and the fraction does not move to pay for it: in-annotation is 87.2% for G0, G1
and G2 alike, exonic 55.5% → 55.6%. By the round-7 rule — a variant whose extra
reads are junk *dilutes* the fraction — G0's extra reads are not junk, so
`--nextseq-trim` is deleting legitimate data.

Supporting measurement, on 200,000 reads of cell 011: `--nextseq-trim=20` drops
~1,970 reads that G0 keeps, and those reads are **not** poly-G — 4% are ≥50% G,
none ≥80%. It is behaving as general quality trimming, not as poly-G removal.

Two things worth keeping from round 12 even though the answer was no:

- Poly-G is genuinely marginal here — 0.217% of reads end in ≥10 G, and the
  `G{20}` adapter accounts for only 506 reads per 200,000 (0.25%).
- Once `--nextseq-trim` is on, `-a polyG` is redundant (G1 vs G2 differ by 3
  uniquely mapped reads). That matters only if a future library makes
  `--nextseq-trim` worth adopting — then drop the adapter rather than keep both.

## Round 13 — pass 0's `strip_polya` is redundant *and* subtly wrong: remove it

`trim_bc_anchor.py` cuts at the barcode anchor and then walks back over the
poly-A itself, via `strip_polya()`. Pass 2's `--poly-a` does the same job. The
overlap is real — and pass 0's copy is the weaker implementation.

**The defect is one missing rule.** cutadapt's algorithm is (a) score each
suffix +1 per A, −2 per non-A and take the max, *and* (b) "exclude all suffixes
from consideration that have more than 20% non-A". `strip_polya()` implements
(a) only. Scoring alone permits up to 33% non-A (a suffix of *a* A's and *b*
non-A's scores `a−2b`, positive iff `b/(a+b) < 1/3`), so between 20% and 33% the
two disagree and pass 0 eats sequence cutadapt would keep. On 400,000 reads of
cell 011 (118,163 anchor hits): **99.910% identical, 106 disagreements, and all
106 are pass 0 over-trimming**, mean 24 nt — worst cases 89→3, 126→49, 125→59 nt.

Not dead code either: `strip_polya` fires on ~97% of anchor hits, 24.5% of all
reads library-wide.

| variant | pass 0 |
|---|---|
| S0 | current, `strip_polya` as written — the baseline |
| S1 | `strip_polya` **removed**; pass 2's `--poly-a` does it |
| S2 | `strip_polya` **fixed** (20% non-A guard added), kept |

| | uniq 011 | exonPC 011 | uniq 016 | exonPC 016 |
|---|---:|---:|---:|---:|
| S0 | 73,304 | 40,680 (55.5%) | 13,930 | 3,190 (22.9%) |
| **S1** | **73,450** | **40,704 (55.4%)** | 13,908 | **3,201 (23.0%)** |
| S2 | 73,319 | 40,682 (55.5%) | 13,931 | 3,190 (22.9%) |

**S1 wins on the decision metric in both cells** (+24 and +11 protein-coding
exonic reads) while in-annotation holds exactly — 87.2% and 75.3% for every
variant. Cell 016 shows the ideal signature: uniquely mapped goes *down* 22
while exonic goes *up* 11, i.e. junk alignments traded for real ones.

**S2 recovers almost nothing** (+2 exonic on 011, 0 on 016; it changes 36 reads
per 300,000). So the guard is not the real story — the function is simply
redundant, and pass 2 additionally gets to decide *after* TrimGalore's quality
trim, which is the better place to decide.

The effect is small (+0.06% exonic). The case for removing it is correctness and
one less reimplementation of a reference algorithm, not yield.

Safe to remove: nothing parses the `poly-A also stripped` log line —
`step2_report.py` reads only `reads`, `anchor found` and `of which exact`.

## Rebuilding the inputs

```bash
cd /nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037
mkdir -p trimtest/in
for c in 011 007 016 001; do
  zcat out/cells/ZHA9292A1_${c}_cbc.fastq.gz | head -n 1200000 \
    | gzip -1 > trimtest/in/${c}.fastq.gz
done
```

Then `sbatch bench_trim.sh` (it does the shared TrimGalore pass) before any of
the later rounds.

## Deriving the read-through adapter for a new library

Do not reuse `TRIM_ADAPTER3` blindly — it is specific to this library's R1
primer. Anchor on the 16 nt that is revcomp of R1's 5' prefix and take a
per-position consensus of what follows:

```python
import gzip, collections
ANCH = "GATCGTCGGACTGTAG"          # revcomp(R1 5' prefix), first 16 nt
cnt = [collections.Counter() for _ in range(80)]
with gzip.open("in/011.fastq.gz", "rt") as f:
    for i, l in enumerate(f):
        if i % 4 != 1: continue
        p = l.strip().find(ANCH)
        if p < 0: continue
        for j, ch in enumerate(l.strip()[p:][:80]): cnt[j][ch] += 1
print("".join(c.most_common(1)[0][0] for c in cnt if c))
```

68,250 of 300,000 reads anchor, ≥97% per-position agreement to position 55.
