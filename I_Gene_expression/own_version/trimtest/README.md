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

That number alone is not enough, because this library's poly-A tail can align
to genomic A-tracts and inflate it. `aligned_composition.py` reads the BAM and
splits uniquely mapped reads into genuine versus **poly-A-only** — aligned
block (soft clips excluded) ≥80% A or T. The recommendation maximises genuine.

## Files

| file | what it does |
|---|---|
| `bench_trim.sh` | round 1: current setting, cutadapt 1.18 vs 5.1, longer homopolymer, a first "model the construct" attempt, and a no-homopolymer control |
| `bench_trim2.sh` | round 2: build up from the control one modifier at a time |
| `bench_trim3.sh` | round 3: literal adapters with a high `min_overlap`, after wildcards proved harmful |
| `bench_trim4.sh` | round 4: same, with the read-through adapter **measured** from the reads |
| `bench_trim5.sh` | round 5: is poly-G free? (yes). **v17 = the adopted setting** |
| `bench_trim6.sh` | round 6: can the leftover 12 nt barcode remnant be removed? (not worth it) |
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
| v18–v20 | v17 plus attempts to remove the leftover 12 nt barcode remnant; all net losses |

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
