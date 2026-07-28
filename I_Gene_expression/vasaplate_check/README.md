# vasaplate_check — reproducing the published VASA-seq pipeline, and checking it against the paper

Cross-check of our re-run of the *published* mapping pipeline against
(a) the count table the authors deposited for the same library and
(b) the numbers actually stated in Salmen & De Jonghe et al. 2022.

Library: **`SRR14783059` / `GSM5369495`, `vasaplate-HEK293T-mESC`** — the only
VASA-plate library in GSE176588. HEK293T is human, mESC is mouse, so it maps to a
concatenated GRCh38+GRCm38 (Ensembl 99) reference.

Nothing here touches `own_version/` (the PM26037 work) or `flashseq/`.

---

## Read this first: what the paper does and does not let you check

The manuscript is nearly silent about this library, and getting this wrong wastes
a day.

- **The paper reports no barnyard number for it.** The 3.08% heterotypic doublet
  rate (p. 1781, Fig. 1d) is **VASA-drop** (`GSM5369496`). GSM5369495 appears only
  as a source of cells for the benchmarking panels.
- **No rRNA percentage for VASA-seq anywhere** — Extended Data Fig. 2f is an
  un-numbered bar chart, the Discussion is qualitative.
- **No tRNA percentage.** Its one tRNA sentence attributes tRNA detection to
  Smart-seq-total, and the Methods name no tRNA database.

So the **deposited table is the primary reference**, not the manuscript.

**The paper states two doublet rules that disagree ~6× on this data.** Fig. 1d
thresholds the *UFI* fraction (>25% from the other species = mixed); Methods p. 18
thresholds the *gene* fraction (>75% from one species = singlet). On the published
table: **0.85% vs 4.82%**. Any barnyard number must name its rule.

**Their own deposited table does not reproduce their own headline.** The VASA-drop
table gives 4.57% / 5.42% against the printed 3.08%. A barnyard mismatch is not
automatically a fault in this pipeline.

Manuscript numbers that *are* checkable here: sncRNA share 1.4%, unspliced
fraction 44.1 ± 10.1%, genes/cell 9,480 ± 1,252 at 75k and 15,248 ± 1,092 at 750k
trimmed reads per cell (the last two need downsampling to be comparable at all).

---

## The four runs, and which one is authoritative

All in `data/ref/fastq_vasaplate/`. Stage 1–4 products are **symlinked** between
run directories, never recomputed.

| run | output prefix | what changed | jobs | state |
|---|---|---|---|---|
| 1 | `vasaplate_out` | first full run | `50542435`–`50542441` | ❌ **stage 7 FAILED** |
| 2 | `rerun_fixednames/` | gene-name fix, stages 6–7 | `50606988`–`50606989` | ✅ superseded |
| 3 | `vasaplate_out_rrnav2` | rRNA reference v2, stages 3–7 | `50788060`–`50788065` | ✅ baseline |
| 4 | `vasaplate_out_bedv2` | annotation BED **with tRNA**, stages 5–7 | `50861367`–`50861370` | ✅ **current** |

> ⚠️ **Run 1's partial tables are still on disk and must not be used.** Its stage 7
> died with an `IndexError` in `reduceGeneName`. The tell: `vasaplate_out_mapStats.log`
> is 8 lines where a complete run writes 21, and every `uniaggGenes_*` /
> `shortGeneNames_*` table is missing. **`status.sh` still defaults to run 1's job
> ids**, so running it bare shows you the failed run.

Runs 3 and 4 were launched with overrides on `a_Mapping/submit_vasaplate_map_array.sh`:

```bash
cd data/ref/fastq_vasaplate
R=/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mixed
p2s=/nemo/lab/turnerj/working/guangxin/vasaseq/code/I_Gene_expression/a_Mapping
START=5 \
REFBED=$R/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed \
RIBOREF=$R/unique_rRNA_human_mouse.v2.fa \
${p2s}/submit_vasaplate_map_array.sh SRR14783059 MIXED 74 vasaplate_out_bedv2 vasaplate_out_bedv2 f
```

`REFBED` was added for run 4 (`START` and `RIBOREF` already existed). Only stages
5–7 read the BED, so `START=5` is the right pairing.

**Setting up a new run directory:** symlink `*_cbc.fastq.gz` (the driver rebuilds
`.cells.manifest` from that glob, so it must be present even at `START=5`) and
`*_E99_Aligned.out.bam` from the previous run's directory. Copying the manifest
alone does not work — it gets overwritten.

---

## References

`$VASA_REFS = /nemo/lab/turnerj/working/guangxin/reference/vasaseq`, all under `mixed/`.

| file | rows / seqs | built by | status |
|---|---|---|---|
| `Human_Mouse_ensembl99.homemade_IntronExonTrna.bed` | 1,165,949, **0 tRNA** | `mixed/build/gtf_to_homemade_bed.py` (untracked) | superseded |
| `Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed` | 1,167,707 = 1,165,949 + **1,758 tRNA** | `own_version/build_annotation_bed_mixed.sh` | **in use, run 4** |
| `unique_rRNA_human_mouse.fa` (v1) | 915 | — | broken (0.27% depletion) |
| `unique_rRNA_human_mouse.v2.fa` | 921 | `own_version/build_rrna_reference_mixed.sh` | used by runs 3–4 |
| `unique_rRNA_human_mouse.v3.fa` | 921 | same, `ORIENT_TO_UNITS=yes` | **built, not yet used** |

### The BED with tRNA

`build_annotation_bed_mixed.sh` reuses `gtf2bed_vasa.py` unmodified by running it
on the already-species-prefixed `combined.gtf`. It proves itself: run in the old
BED's configuration it reproduces all **1,165,949** rows, sorted-set identical —
which is what establishes the rules of the untracked builder it replaces.

Two decisions were settled **by measuring against the published table**, not by
preference:

- **tRNA names are bare, no species tag.** Published rows carry no species marker,
  and the two GtRNAdb sets (619 hg38 + 1,139 mm10) collide on exactly **one** name,
  `16.tRNA1-LysNNN`. Only the *contig* column is species-prefixed, because that has
  to match the STAR index.
- **`BED_COORD=asis`.** If the authors had used true 0-based while we use `asis`,
  short non-splicing biotypes would sit ~8% *below* protein-coding in a per-gene
  ours/published ratio. Measured on 72,613 shared genes: **+1.12%** — wrong sign,
  wrong size. So `asis` is this script's default, the **inverse** of the mouse
  `build_annotation_bed.sh` default, deliberately.
- **`TRNA_SET=all`**, because the published table has 1,130 distinct simple tRNA
  rows and `high` offers at most 561+465 = 1,026 loci.

**Their tRNA locus ids cannot be reproduced — do not try.** Only 503 of the 1,130
published names exist in current GtRNAdb; 110 published rows use `Undet`, the
retired tRNAscan-SE spelling, where current GtRNAdb writes `Und`, and the loci have
been renumbered. Compare at **isotype** level, where 59 of 62 classes (95%) match —
and isotype is exactly what step 7 collapses tRNA to.

> `_tRNA.*Counts.tsv` **merges the two species.** `countTables_fromPickle.py` groups
> tRNA rows by `rsplit('.')[-1]`, so human and mouse `ValAAC` land in one row. True
> of the published tables too. Species-resolved per-locus tRNA survives only in the
> `total` tables, which keep the full row name.

---

## Results (run 4, with tRNA; run 3 alongside)

From `res/vasaplate/comparison_summary.tsv`. **Never hand-transcribe these — read
the TSV.**

| quantity | published | ours |
|---|---|---|
| simple gene rows | 77,207 | 80,234 |
| shared | — | **72,613** |
| per-gene Spearman r | — | **0.974** |
| per-gene Pearson r (log10) | — | **0.982** |
| median per-cell Pearson r | — | **0.982** |
| median log2(ours/published) | — | **0.000** |
| barcodes ≥ 7,500 UFIs | 353 | 353 |
| doublet rate, Fig. 1d rule | 0.85% | 0.57% |
| doublet rate, Methods rule | 4.82% | 8.78% |
| HEK293T median UFIs / genes | 190,110 / 18,316 | 185,227 / 18,757 |
| mESC median UFIs / genes | 75,233 / 11,852 | 72,893 / 11,930 |
| sncRNA share | 1.421% | **1.604%** (paper states **1.4%**) |
| **tRNA rows detected** | **1,130** | **224** (148 shared) |
| **tRNA UFI share** | **0.0241%** | **0.0037%** |

Adding tRNA did not disturb anything else: per-gene Spearman 0.9739 → **0.9732**,
median per-cell Pearson 0.9821 → **0.9820**, shared rows 72,613 → **72,760**.

### The tRNA shortfall is geometric, not an annotation defect

| | |
|---|---|
| tRNA feature length, median | **72 bp** (1,758 loci) |
| read length, mean | **70.8 bp** |
| reads LONGER than the feature | 37.5% — can never be contained |
| reads equal in length | 8.3% |
| reads actually kept (`jS:IN`) | **9.5%** (33 of 349, five cells) |

Step 6 keeps a non-splicing biotype only when the read falls **entirely inside**
the feature. A tRNA gene is about the length of a read in this library, so
containment is close to geometrically impossible.

**Why the authors got 5× more is unresolved.** Same FASTQ, so the same read
lengths — which points at their tRNA features being longer than GtRNAdb's
mature-tRNA bounds (a precursor annotation with 5′ leader and 3′ trailer would fit
reads inside). Their annotation was never deposited and the Methods name no tRNA
source, so this is a hypothesis. **Do not present the 224 as a reproduction of the
1,130.**

> A retracted claim, kept so it is not re-derived: an earlier pass measured "62.8%
> of tRNA reads overhang by exactly 1 bp, so `BED_COORD=fix` would recover them".
> That was computed on the wrong columns — in the stage-5 output, column 7 is the
> `CG:…;nM:…;jS:…` info string and column 8 is gene length, **not** gene
> start/end. The 9-column layout is documented in `own_version/README.md`
> ("The 9 output columns"). `fix` would not recover these reads.

The Methods-rule rates diverge more than the Fig. 1d rates because that rule counts
*genes*, so it is dominated by genes detected at one or two counts, and we detect
more rows overall. The UFI-based rule — the one Fig. 1d is drawn with — agrees.

---

## The residual rRNA, diagnosed

Every biotype agrees with the published table within 2× **except `rRNA`**, ~600×
higher in ours — and it is **one locus**: `ENSMUSG00000106106_CT010467.1_rRNA`,
**95,823 UFIs against 72**. Two defects, each harmless alone.

**A. `riboread-selection.py` writes reverse-strand reads reverse-complemented.**
It emits `r0.seq`, and pysam's `.seq` returns the sequence *as stored in the BAM* —
for a reverse alignment, the reverse complement. `.qual` is reversed too. Every
read kept by the stranded reprieve enters STAR flipped, which for a stranded
protocol inverts sense and antisense. Cell 005: 2,222 of 45,880 reads with an rRNA
hit (4.8%).

**B. One reference entry was stored antisense.** `bedtools getfasta -s` gives the
sense strand of the *annotated gene*, which equals transcript sense only where
Ensembl put the gene on the transcript's strand. In Ensembl 99, `ENSMUSG00000106106`
is annotated opposite — so its entry was the 18S region backwards. **1 of 915.**

Together: a genuine 18S read aligns reverse to that entry, is spared, is written
flipped by A, and STAR then maps it *sense* to the minus-strand gene, where it is
counted as rRNA.

Why the backwards entry wins: it and the 47S unit capture the **same** reads (30/30
on the test set), but it scores **higher** (mean AS 67.3 vs 64.7), being the exact
GRCm38 locus the reads came from where BK000964.3 is a curated consensus. `bwa mem`
reports one primary alignment, so the backwards entry takes it and returns REVERSE.

**Fixed in `build_rrna_reference_mixed.sh` STEP 4** (`ORIENT_TO_UNITS=yes`, default):
align every Ensembl entry to the NCBI units, reverse-complement any returning flag
16, re-align and require zero remaining. Fixes the **class**, not the instance.
Output is `v3`; verified v3 vs v2 = same 921 names, 920 sequences byte-identical,
exactly 1 reverse-complemented, 0 other changes.

Expected effect of switching to v3 — note the two move in **opposite** directions:

| metric | now (v2) | with v3 |
|---|---|---|
| step-3 depletion, 384 cells | 5.06% (9,176,750 / 181,287,059) | ~5.34% (**up** ~0.28 pt) |
| rRNA biotype share in the table | 0.201% (95,823 / 47,732,501 UFI) | ~0.0003% (**down** ~600×) |

Measured on 8 cells spanning a 20× depth range: 0.278% of reads move from "kept" to
"ribosomal", stable at 0.27–0.51% per cell.

> **`own_version/` is NOT affected.** `unique_rRNA_mouse.v2.fa` (GRCm39 / Ensembl 116)
> has **0 antisense entries** among its 356 — checked. PM26037's 21.39% stands. The
> fix was deliberately not applied there.

**What we cannot say:** that the authors "used only NCBI sequences and no Ensembl
dump". Their FASTA was never deposited, the Methods name no accessions, and
`a_Mapping/README.md`'s list is prefixed "e.g." — an example, not a specification.
The published table argues against it too: on every rRNA gene *except* this one
their residual is equal to or higher than ours (8 vs 2 on Gm25908; 13 and 11 on
human RNA5.8SN2/N4 where we have none). The supportable statement is narrower:
whatever their reference was, it did not carry this locus in an orientation that
produces the artefact.

---

## The scripts

| file | what it does |
|---|---|
| `vp_common.py` | paths, loaders, species logic, both doublet rules, the paper's stated numbers |
| `01_compare.py <run>` | the analysis. Writes 4 TSVs into `res/vasaplate/` |
| `02_figures.py` | 5 figures, PDF (light) + PNG (light & dark) for the report |
| `03_report.py` | assembles the self-contained HTML |
| `04_diagnose_residual_rrna.sh [cell]` | reproduces all 5 pieces of the rRNA diagnosis from files on disk |

```bash
PY=/camp/apps/eb/software/Anaconda3/2024.10-1/bin/python   # PYTHONPATH must be cleared
env -u PYTHONPATH $PY 01_compare.py bedv2   # ~50 GB RAM, ~1 min; sbatch it
env -u PYTHONPATH $PY 02_figures.py
env -u PYTHONPATH $PY 03_report.py
```

Outputs land in `res/vasaplate/`: `vasaplate_paper_check.html`, `figures/`,
`comparison_summary.tsv`, `per_cell.tsv`, `gene_concordance.tsv`,
`biotype_composition.tsv`.

Figure palette is the dataviz reference categorical instance, validated in **both**
light and dark before use. Light mode raises a contrast WARN on two slots, which
obligates visible labels or a table view — both are present.

---

## Where this stands (2026-07-27, end of session)

**Run 4 is mid-flight.** Stage 5 finished 384/384 with no errors; tRNA confirmed
flowing (12 sampled cells: 906 tRNA rows vs 0 under the old BED).

| stage | job | state |
|---|---|---|
| 5 b2bs / b2bm | `50861367` / `50861368` | ✅ COMPLETED 384/384 |
| 6 cout | `50861369` | 🔄 RUNNING, ~9 h |
| 7 pick | `50861370` | ⏳ PENDING (dependency), ~5 h |

Check with `squeue -u $USER`; cancel with `scancel 50861369 50861370`.
**Confirm the job is gone from `squeue` before reading any output** — per-cell files
appear while a stage is still running, so counting files is not proof.

### To do next, in order

1. **When run 4 finishes**, verify first:
   - `vasaplate_out_bedv2_mapStats.log` shows **non-zero** `Total reads assigned to tRNA`
     (run 3 shows 0);
   - `vasaplate_out_bedv2_tRNA.ReadCounts.tsv` has rows (run 3's is header-only);
   - the log is the full **22 lines**, not 8 (that is how run 1's failure looked).
2. ~~Regenerate the report against run 4~~ — **done**. Outputs in `res/vasaplate/`.
3. **Decide what to do about the tRNA shortfall.** The geometric cause is measured;
   the open question is whether to test a padded/precursor tRNA annotation. That is
   inventing annotation unless a source can be named, so it needs a decision, not a
   default.
4. **Then, as a separate run**, switch to `unique_rRNA_human_mouse.v3.fa` and re-run
   stages 3–7 into `vasaplate_out_rrnav3` (`START=3`, `RIBOREF=…v3.fa`,
   `REFBED=…v2.bed`). Keep it separate from the tRNA change so the two variables stay
   separable. Success = the 95,823 UFIs on `ENSMUSG00000106106` collapse to
   published-level, and depletion rises 5.06% → ~5.34%.

### Known-open

- **Defect A is not fixed.** `riboread-selection.py` still writes reverse-strand
  reads flipped. It is an upstream published script, so under this repo's conventions
  a genuine bug there is fixed in place and noted — but it changes every run's step 3,
  and the authors ran the same code, so fixing it moves *away* from reproducing the
  paper. Left as a documented decision, not an oversight.
- The mixed rRNA FASTA omits the paper's human `45S9` / `45S1-17` entries.
- The mixed BED still lacks the paper's small-ncRNA overlap-reduction hierarchy.
- `submit_vasaplate_map*.sh` still hardcode the **v1** rRNA fasta and the **v1** BED;
  v2/v3 are reachable only via `RIBOREF=` / `REFBED=`.
- `mixed/build/build_mixed_reference.sh` and `mixed/build/gtf_to_homemade_bed.py` are
  untracked builders sitting next to their outputs — the exact arrangement that lost
  `build_mouse_reference.sh`.
- **Do not delete `mixed/build/`**: `flashseq/config.sh` reads `human.genome.fa` and
  `Homo_sapiens.GRCh38.99.gtf.gz` from it, and `build_annotation_bed_mixed.sh` reads
  `combined.gtf`.
