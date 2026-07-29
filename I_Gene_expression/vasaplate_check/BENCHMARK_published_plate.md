# Gold-standard benchmark — the published VASA-plate through this workflow

Purpose: establish whether the workflow used for the FLASH-seq ↔ VASA comparison
reproduces the **published** VASA-plate library, so that a three-way comparison
(published plate / own plate / FLASH-seq) rests on a validated anchor rather than
on an assumption.

Status 2026-07-29. **The published plate had already been re-run four times in
this repo before this session**, and run 4 already passed against the deposited
count table. This session's contribution is run 5: the one large remaining
discrepancy, traced but never fixed.

---

## What was already established (runs 1–4, `vasaplate_check/`)

Not re-derived here — read from `vasaplate_check/README.md` and
`res/vasaplate/comparison_summary.tsv`.

**The primary reference is the deposited count table, not the manuscript.** The
paper reports no barnyard number, no rRNA percentage and no tRNA percentage for
this library (`GSM5369495`); it appears mainly as a source of cells for the
benchmarking panels. So concordance is measured against
`GSM5369495_vasaplate_HEK293T-mESC_split_total.UFICounts.tsv.gz`.

| quantity | published | run 4 (`vasaplate_out_bedv2`) |
|---|---|---|
| simple gene rows | 77,207 | 80,234 |
| shared rows | — | 72,613 |
| per-gene Spearman r | — | **0.974** |
| per-gene Pearson r (log10) | — | **0.982** |
| median per-cell Pearson r | — | **0.982** |
| median log2(ours/published) | — | **0.000** |
| barcodes ≥ 7,500 UFIs | 353 | **353** |
| HEK293T median UFIs / genes | 190,110 / 18,316 | 185,227 / 18,757 |
| mESC median UFIs / genes | 75,233 / 11,852 | 72,893 / 11,930 |
| sncRNA share | 1.421% | 1.604% (paper states 1.4%) |
| **rRNA biotype** | — | **~600× higher — the open defect** |
| tRNA rows detected | 1,130 | 224 (148 shared) |

A median log2 ratio of exactly 0.000 with r = 0.974–0.982 across 72,613 genes is
a reproduction, not a resemblance. **The workflow is sound on everything except
rRNA.**

### Two discrepancies that are NOT workflow faults

- **tRNA, 224 vs 1,130.** Geometric, not annotational: tRNA features have median
  length 72 bp against a 70.8 bp mean read length, and stage 6 keeps a
  non-splicing biotype only if the read falls *entirely inside* the feature
  (`jS:IN`). 37.5% of reads are longer than the feature and can never be
  contained; 9.5% survive. The authors' tRNA annotation was never deposited and
  their Methods name no tRNA source, so why they recovered 5× more is a
  hypothesis (a precursor annotation with 5′/3′ flanks would fit reads inside),
  not a finding. **Do not present 224 as a reproduction of 1,130.**
- **The two doublet rules disagree ~6× on the same data.** Fig. 1d thresholds the
  UFI fraction, Methods p. 18 the gene fraction: 0.85% vs 4.82% on the published
  table itself. And the authors' own VASA-drop table gives 4.57%/5.42% against
  their printed 3.08%. A barnyard mismatch here is not automatically our fault.

---

## The open defect, and this session's run 5

The `rRNA` biotype is **~600× higher** than published, and it is **one locus**:
`ENSMUSG00000106106_CT010467.1_rRNA`, **95,823 UFIs against 72**. Diagnosed as
two defects that are each harmless alone:

**A.** `riboread-selection.py` emits `r0.seq`, and pysam returns the sequence *as
stored in the BAM* — for a reverse alignment, the reverse complement. So every
read kept by the stranded reprieve enters STAR flipped, inverting sense and
antisense for a stranded protocol.

**B.** One reference entry was stored antisense. `bedtools getfasta -s` gives the
sense strand of the *annotated gene*, which equals transcript sense only where
Ensembl put the gene on the transcript's strand. In Ensembl 99,
`ENSMUSG00000106106` is annotated opposite — so its entry was the 18S region
backwards. **1 of 915.**

Together: a genuine 18S read aligns reverse to that entry, is spared, is written
flipped by A, and STAR then maps it *sense* to the minus-strand gene, where it is
counted as rRNA.

The fix (`build_rrna_reference_mixed.sh` STEP 4, `ORIENT_TO_UNITS=yes`) aligns
every Ensembl entry to the NCBI 47S units, reverse-complements any returning flag
16, then re-aligns and requires zero remaining. It fixes the **class**, not the
instance. Output is `unique_rRNA_human_mouse.v3.fa` — **built in a previous
session and never run.**

### Precheck (Rule 2), job `1103ee07` — exit 0

Before committing 6–10 h over 384 cells, `vasaplate_check/05_v3_precheck.py`
asserted, read-only:

| check | result |
|---|---|
| v3 vs v2 is exactly the documented change | 921 names identical and in order; **920 sequences byte-identical, exactly 1 reverse-complemented, 0 changed any other way** |
| the flipped entry is the blamed locus | `mouse_ENSMUSG00000106106_rRNA(-)` ✓ |
| no antisense entries remain in v3 | 6 unit entries identified; **0** entries match the units better reversed than forward |
| stage 1–2 inputs resolve | 384 `_cbc.fastq.gz`, **0 broken symlinks**, 0.4–67.6 MB |
| bwa index complete and not stale | all five files present, none older than the fasta |

### Expected effect — the two metrics move in OPPOSITE directions

This is the prediction the run tests, recorded before it finished:

| metric | run 4 (v2) | predicted with v3 |
|---|---|---|
| step-3 depletion, 384 cells | 5.06% (9,176,750 / 181,287,059) | ~5.34% (**up** ~0.28 pt) |
| rRNA biotype share in the table | 0.201% (95,823 / 47,732,501 UFI) | ~0.0003% (**down** ~600×) |

Predicted from 8 cells spanning a 20× depth range: 0.278% of reads move from
"kept" to "ribosomal", stable at 0.27–0.51% per cell.

### The run

`START=3 RIBOREF=…v3.fa REFBED=…IntronExonTrna.v2.bed` on
`a_Mapping/submit_vasaplate_map_array.sh`, **unmodified** — overrides via
environment only (Rule 1: fork, never patch). Output `vasaplate_out_v3/`.

Stage chain `51029413` (ribo) → `51029414` (gmap) → `51029415`/`51029416`
(b2bs/b2bm) → `51029417` (cout) → `51029418` (pick). Read from the driver's own
output.

Stages 1–2 are **symlinked, not recomputed** — they are upstream of the rRNA
reference and already concordant at r = 0.974, so re-running them would test
nothing and cost days.

### Three wrapper failures on the way, all mine, none in the pipeline

Recorded because two are new instances of a trap class the conventions already
name — a non-zero exit from a command nobody thinks of as fallible.

1. **Attempt 1** (`51029201`–`51029206`, cancelled after ~1 min): I symlinked
   only the stage-1 `_cbc.fastq.gz`. Stage 3 consumes the stage-2
   `_cbc_trimmed_homoATCG.fq.gz`, so all 384 ribo tasks failed with `fail to open
   file`. **No pipeline product was written** — verified before cancelling. The
   trimmed fastqs live in `vasaplate_out/` and `vasaplate_out_rrnav2/`, not in
   `vasaplate_out_bedv2/`.
2. **Attempt 2** (exit 141): my own readability check `zcat f | head -4` made
   `zcat` die of SIGPIPE when `head` closed the pipe; `pipefail` promoted 141 and
   `set -e` killed the wrapper **before the driver ran**, so nothing was
   submitted. Now reads a bounded byte prefix of the file instead, with
   `pipefail` disabled around it.
3. The first attempt's wrapper also reported exit 1 on an empty `grep` for
   "Submitted batch job" — the driver prints `STAGE N … : <jobid>` instead. The
   driver had actually succeeded.

TO BE FILLED — run 5 results, and whether the two predicted movements occurred.
