# What VASA's `smartseq_noUMI` branch actually does

Scope: the branch exists in upstream `code/I_Gene_expression/a_Mapping/` and had
never been run in this repo before this note. Everything below was read out of
the source at the commit this repo has.

**Provenance of every number here, stated explicitly.** Claims are labelled:

- **[code]** — read off the source. No execution needed.
- **[precheck]** — computed by `noumi_precheck.py` on real FLASH-seq read names
  and the real v2 BED (job `7d2eead1`, report `noumi_precheck_A9.txt`).
- **[dryrun]** — computed by the end-to-end slice and
  `noumi_dryrun_report.py` (report `noumi_dryrun_A9.txt`).

Anything not so labelled is a **prediction from [code] that the dry run is
designed to test**, and is marked *predicted*. Line numbers refer to
`a_Mapping/countTables_2pickle_cellsSpliced.py` (130 lines, "2pickle") and
`a_Mapping/countTables_fromPickle.py` (287 lines, "fromPickle").

Companion files in this directory:

| file | what it is |
|---|---|
| `noumi_precheck.py` | read-only validator; replays the branch's string ops over real read names and the whole v2 BED before anything is submitted |
| `dryrun_a9_noumi.sh` | end-to-end slice on ZHA8833A9, paired-end **and** single-end arms |
| `noumi_dryrun_report.py` | sanity report over the tables step 7 writes |

---

## 1. The read-name contract

`protocol` reaches the pickle builder as `sys.argv[3]` (2pickle line 16) and is
passed through `gene_assignment` into `get_UMI` (lines 47-49). The branch is
2pickle lines 42-44:

```python
elif protocol == 'smartseq_noUMI':
    name_info = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in name.rsplit(';')[1:]}
    umi = 'A'
```

`name_info` is **built and then discarded** — the returned UMI is the literal
`'A'` either way. But the dict comprehension still executes, so the read name
must not make it raise. It raises `IndexError` iff some `;`-separated field
*after the first* contains no `:`. Note the asymmetry with `protocol='ramda'`
(line 40-41), which is `umi = 'A'` and nothing else: `ramda` imposes no contract
at all, `smartseq_noUMI` imposes this vestigial one.

**Requirement, stated exactly:** a read name is acceptable iff every
`;`-separated field after the first contains at least one `:`. A name with **no
`;` at all** yields `rsplit(';')[1:] == []`, an empty comprehension, and cannot
raise.

**FLASH-seq satisfies this untransformed.** Read names in
`ZHA8833A9_S108_L007_R1_001.fastq.gz` are plain Illumina:

```
@LH00442:237:23GT7GLT3:7:1101:32183:1048 1:N:0:TTAAGTGC+NCCTGTCG
```

**[precheck]** Over 18,902 names sampled at stride 2048 across all 38,710,566
read records, in all three forms the name can reach `get_UMI` in — raw header,
STAR-truncated (the BAM QNAME, first whitespace token), and SM-tagged —
**0 raise, 18,902 pass**, and no name contains a `;` of its own. The same names
through `protocol='vasa'` raise `KeyError: 'RX'` 2000/2000, which is what the
branch buys.

**So no read-name transformation is needed for the UMI field.** VASA's step 1
(`a_Mapping/concatenator.py`, called via `extractBC.sh`) is what normally
creates the tagged shape — line 175:

```python
name = ';'.join([n1] + [':'.join(x) for x in zip(['SS','CB','QT','RX','RQ','SM'], ...)])
```

producing `…:13708;SS:AGCTTC;CB:AGCTTC;QT:jjjjjj;RX:AGTTAC;RQ:jjjjjj;SM:005`.
For FLASH-seq that step is not needed and must not be run: there is no cell
barcode to extract and no UMI to record.

### The one place a name *is* needed: the column name

Step 6 derives the column (cell) id two ways (2pickle lines 89-92):

```python
if cellidFROMfilename == 'f':
    cellid = cellfile[:cellfile.index('_cbc')]
else:
    cellid = {x.rsplit(':')[0]: x.rsplit(':')[1] for x in df.iloc[0]['Name'].rsplit(';') if ':' in x}['SM']
```

**[precheck]**, replayed both ways over real strings:

| mode | input | result |
|---|---|---|
| `f` | `ZHA8833A9_cbc_noumi_…_genes.bed.gz` | `'ZHA8833A9'` |
| `f` | `ZHA8833A9_E99_…_genes.bed.gz` (no `_cbc`) | `ValueError: substring not found` |
| `r` | `LH00442:237:23GT7GLT3:7:1101:32183:1048` | `KeyError: 'SM'` |
| `r` | same `+ ';SM:ZHA8833A9'` | `'ZHA8833A9'` |

An untagged Illumina name *has* colons, so the dict builds fine; it simply has
no `SM` key. So the honest choice is:

- **`CELLID_FROM=f`** (what the dry run uses) — costs one filename convention:
  the per-library stem must contain the literal `_cbc`. No read is touched.
  Caveat: step 6 keeps the folder string it globbed, so the column name comes
  out as `cells/<stem>` — strip it downstream.
- **`CELLID_FROM=r`** — needs `;SM:<lib>` appended to every read name, i.e. a
  rewriting pass over 29 M reads per library.

`f` is the minimal, honest transformation, and it transforms a filename rather
than data.

---

## 2. What the branch does differently from `protocol='vasa'`

**There is no deduplication anywhere, on either path.** VASA does not
deduplicate by UMI either — it *counts distinct UMIs* as a separate table
(`countTotalUMI`, fromPickle 55-56). The no-UMI branch changes how reads are
retrieved from the per-gene dict, not whether anything is collapsed.

Every read is filed under the single literal key `'A'`, so a cell entry is
`{'A': Counter({'exon': 7, 'intron': 3})}` instead of VASA's
`{'AGTTAC': Counter({'exon': 1}), 'GCGGAC': Counter({'intron': 1}), …}`.

Three counters branch on `protocol`, and three do not:

| counter | lines | branches? | no-UMI behaviour |
|---|---|---|---|
| `countTotalReads` | 47-54 | yes | `sum(x['A'].values())` instead of summing over all UMIs |
| `countExonReads` | 59-68 | yes | `'exon' in k` (**substring**) instead of `'intron' not in [...]` (**exact membership**) |
| `countIntronReads` | 73-82 | yes | `'intron' in k` (**substring**) instead of `'intron' in [...]` |
| `countTotalUMI` | 55-56 | **no** | `len(x)` — always 1 for a detected gene |
| `countExonUMI` | 69-71 | **no** | uses VASA's exact-membership test on the single key |
| `countIntronUMI` | 83-85 | **no** | ditto |

Consequences, in order of how much they matter:

1. **`UFICounts` degenerates to a detection mask.** `countTotalUMI` has no
   `protocol` argument, so on a no-UMI pickle it returns `len({'A': …}) == 1`
   for every detected gene and 0 otherwise. **[code]**; *predicted* for the real
   tables, and `noumi_dryrun_report.py` §2 checks it by listing the distinct
   values the written `_total.UFICounts.tsv` actually holds.
2. **`TranscriptCounts` is that same mask, and is not a quantification.**
   fromPickle line 87 reads the UMI length off the *first UMI it finds in the
   pickle*, then line 88 sets `K = 4**len(umi)`. With `umi == 'A'`,
   **K = 4**, not 4096. `bc2trans` (lines 89-96) is then applied to the UFI
   table at lines 102, 239 and 247-249. `bc2trans(0) = 0` and
   `bc2trans(1) = 1.0` exactly, at K = 4 and at K = 4096 alike — **[precheck]**,
   which tabulates `bc2trans` at both K. Since UFI ∈ {0,1} on this path, every
   `TranscriptCounts` cell should equal its `UFICounts` cell (*predicted*;
   `noumi_dryrun_report.py` §2 tests it elementwise on the written tables).

   **This is the critical point.** `bc2trans` is the `bc2trans` collision
   correction and it *is* still computed and written — nine TSVs of it. It is
   not applied to read counts (it is applied to the UFI table, which is the
   mask), so it is not the "garbage transform on read counts" failure mode; it
   is the *identity on a mask*. Either way the output carries no abundance
   information. **Do not use any `*.TranscriptCounts.tsv` from the no-UMI path
   as a quantification. Use `ReadCounts`.**
3. **The exon/intron branches are not equivalent.** VASA tests exact list
   membership (`'intron' not in ['-'.join(set(k.rsplit('-')))…]`), the no-UMI
   branch tests substring containment (`'exon' in k`). A combination label
   `'exon-intron'` — one read spanning two genes, one exonically and one
   intronically — is not the literal string `'intron'`, so VASA counts it as
   **exon** only; the no-UMI branch matches both substrings and counts it in
   **exon and intron**. So `spliced + unspliced` can exceed `total` on the
   no-UMI path. **[code]**, demonstrated on a synthetic entry by **[precheck]**
   §4 (`vasa` total=5 exon=5 intron=0; `smartseq_noUMI` total=5 exon=5
   intron=5). How much it bites on real data is **[dryrun]** §3 — step 7 writes
   the spliced/unspliced tables only for uni-genes plus single-label
   multi-genes, so the excess is an empirical question, not a derivable one.
4. **`countTotalReads` on the wrong pickle.** The no-UMI branch would raise
   `KeyError: 'A'` on a VASA pickle; the `vasa` branch on a no-UMI pickle gives
   the *same* total (one UMI to sum over). Only the exon/intron split differs.
5. **`gene_assignment` (2pickle 47-77) is protocol-independent.** Best `nM`,
   `jS:IN` priority, non-spliceable-biotype priority, `'-'.join` of gene names
   on a tie — identical on both paths. The read→gene assignment is the same
   measurement on both sides of the comparison.

### The assignment shell scripts assume neither barcode nor UMI

`deal_with_singlemappers.sh` / `deal_with_multimappers.sh` never parse the read
name. They append to it (`$1=$1";CG:"$6";nM:"nm`) and later split it on the
first `;CG:` (`sx=match(readname, /;CG:/)`) to separate name from tags. Both
work on any QNAME. Two things they *do* care about:

- **`stranded`** (`$3`). FLASH-seq is unstranded — the forward strand carries
  49.1-50.5% of ribosomal reads (measured, `code/flashseq/README.md`) — so it
  must be called with `n`. VASA runs with `y`. Passing `y` here would halve
  every biotype figure.
- **`NH` parsing.** Upstream's `/NH:i:1\tHI:i:1\t/` and `/NH:i:[2-9]/` both miss
  `NH:i:10..19`; `own_version/deal_with_singlemappers.sh` parses NH as a number.
  The dry run uses the `own_version` fork for singles and upstream for multis,
  exactly as `own_version/pipeline.sh` does.

---

## 3. The UMI ceiling does not apply

On the VASA side, 8 genes saturated the 6 nt UMI space (K = 4096) and had to be
dropped because `bc2trans` clamps at `x >= K` (line 90) to the constant
62356.12277. Without UMIs:

- UFI per gene per column is **at most 1** (one literal key `'A'`). **[code]**
- The clamp branch needs UFI ≥ K. At K = 4 that is unreachable from UFI ≤ 1.
  (For the record the clamp value at K = 4 is 28.8306 — **[precheck]** computed
  it — but no cell can reach it.)

**So the ceiling cannot arise, and the FLASH-seq analysis-set filter must
differ from VASA's.** VASA's filter removed genes whose `TranscriptCounts` were
clamped; that criterion is undefined here. FLASH-seq must be filtered on
`ReadCounts`, and only `ReadCounts`.

The other filter to set deliberately rather than inherit is step 7's
`filt_unigenes` (argv[4]). At `y` it computes
`ncells = max(5, round(0.01*len(cntdf.columns)))`, written for a 384-cell plate
(5/384 = 1.3%). **[precheck]** evaluates it at the counts that matter here:

| columns | ncells | as % of columns |
|---|---|---|
| 1 | 5 | 500% |
| 2 | 5 | 250% |
| 10 | 5 | 50% |
| 12 | 5 | 41.7% |
| 16 | 5 | 31.2% |
| 384 | 5 | 1.3% |

At 10 FLASH-seq libraries the threshold is 50% of them, and at 1 no gene can
pass at all — `uni_genes_filt` is empty and `reduceGeneName` silently loses its
"exactly one component is a known unigene" rule. **Pass `n`** and apply a
filter chosen for this comparison downstream.

---

## 4. Paired-end

**Answer, measured: paired-end input does NOT work. See §5a for the mechanism
and the rates.** The path is single-end only.

The reason is not where it looks like it should be. VASA reads are short
single-end after barcode extraction; FLASH-seq is 151 nt PE. Nothing in step 6
or the assignment scripts inspects the SAM FLAG or the mate fields, so PE input
is not rejected at any of the obvious places, and the read→gene logic would in
fact have handled it:

- STAR emits one alignment record per mate, both sharing one QNAME.
- Both mates carry the same `NH`, so both survive the singlemapper filter and
  both become BED rows.
- Step 6 groups the BED by the `Name` column, which is the QNAME
  (2pickle line 101: `gdf = df.groupby('Name')`), so the two mates of a
  fragment fall into **one** group and are assigned **once** — a fragment
  count, which is what you want.
- But STAR's `Log.final.txt` "Number of input reads" counts *pairs*, and the
  two mates can land on different features, in which case
  `gene_assignment`'s tie-breaking treats them as a multi-feature read.

All four points above are **[code]**, and all four hold. What kills PE is none of
them: it is `bedtools bamtobed` appending the `/1`,`/2` mate suffix into the
middle of the tag string that step 6 must `int()` — **[dryrun]**, §5a, 100% of
1,033,528 PE rows.

That is why the dry run maps both arms instead of reasoning about one: the four
plausible failure modes above were all fine, and the actual one was in a line
nobody would have inspected. `noumi_precheck.py` §5 now tests exactly this
contract on real step-5 BEDs, so it is a five-second check rather than an
hour-long job.

---

## 5. What the dry run actually did — [dryrun]

ZHA8833A9, stride-64 subsample = **604,853 read pairs** out of 38,710,566 read
records. Untrimmed, deliberately (see the driver header, note 7). Jobs:
`7d2eead1` (precheck), `5b727925` (the PE failure), `3b69a369` (the run that
completed, exit 0). Report: `noumi_dryrun_A9.txt`.

**Verdict: the branch works, single-end only.** Step 7 wrote all 20 tables and
every internal consistency check passed. Three predictions confirmed exactly:

| prediction | result |
|---|---|
| `UFICounts` ∈ {0,1}, a detection mask | distinct values = `[1]`; `sum(UFI) − count(reads>0) = 0` |
| `TranscriptCounts` == `UFICounts` elementwise | `True` |
| protein-coding dominates a poly-A library | **86.47%** of the non-rRNA remainder |

Scale, SE arm: 604,853 input reads → 313,679 uniquely mapped → 534,978
singlemapper BED rows → **336,369 reads in the tables across 15,781 genes**
(0.556 of STAR input; the shortfall is the 42.17% unmapped-too-short plus reads
overlapping no annotated feature). Unspliced fraction 0.0822. Top genes are
*Ftl1*, *Hsp90ab1*, *Eef1a1*, *Hspa8*, *Dppa5a*, ribosomal proteins —
plausible for an embryonic poly-A library.

`A9dry_mapStats.log` is **21 lines**, i.e. complete (the repo's `CLAUDE.md` says
22; it is wrong — 21 unconditional `fout.write` calls). It reconciles: 336,369
total mapped reads, 15,781 gene entries collapsing to 12,866 after aggregation,
290,422 reads on 10,734 uni-genes before aggregation and 316,621 on 10,891
after, 19,748 on 1,975 multi-genes. **0 reads assigned to tRNA** (see §5b).

**The exon/intron double-count did not materialise here:** spliced + unspliced =
336,369 = total, excess **0**. The mechanism in §2.3 is real, but step 7 writes
those tables only for uni-genes plus single-label multi-genes, and mapStats
confirms why it could not bite in this slice: *"multi-genes after aggregation
that have multiple labels (exon-intron etc) = 0"*. Not a general guarantee — the
check stays in the report script, and on a library where that line is non-zero
the excess will be real.

### 5a. Paired-end does not work, and this is upstream's contract, not a miscall

Step 6 line 97 is:

```python
df['nMs'] = df.apply(lambda x: int(x['Info'].rsplit(';nM:')[1].rsplit(';jS:')[0]), axis=1)
```

`bedtools bamtobed` appends `/1` or `/2` to the read name of a paired-end mate.
`deal_with_singlemappers.sh` appends `;CG:<cigar>;nM:<n>` to the QNAME **before**
bamtobed runs, so the mate suffix lands on the **end of the nM value**:

```
INFO=CG:116M1D35M;nM:0/2;jS:IN     <- PE
INFO=CG:62M1D89M;nM:0;jS:IN        <- SE
```

`int('0/2')` raises. Measured: **1,033,528 of 1,033,528 PE rows (100%)
unparseable, 0 of 534,978 SE rows (0%)**. Step 6 dies on the first PE row with
`ValueError: invalid literal for int() with base 10: '0/2'`. It is 100%, so
there is no partial-usability question. `noumi_precheck.py` §5 now checks this
contract against real BEDs; the PE outputs are quarantined in
`cells_pe_unusable/`, not deleted.

The rest of the PE arm is sound — 653,581 mate1 and 653,581 mate2 records,
1,033,528 BED rows collapsing to 292,614 distinct QNAMEs, i.e. the QNAME
grouping *would* have given a fragment count. Only the nM field blocks it.

### 5b. The tRNA tables came out EMPTY, and the reason threatens the comparison

`A9dry_tRNA.ReadCounts.tsv` has a header and **no rows**. Not a crash — a
structural consequence, and the most important thing this dry run found.

Chain of evidence:

1. Only **6** rows of the SE step-5 BED overlap a tRNA feature (VASA cell 005,
   for comparison: **229**).
2. Of those 6, **0 have `jS:IN`** — 5 are `jS:5`, 1 is `jS:3`.
3. Step 6 line 100 drops any row whose biotype is not in `biotypeWsplicing`
   unless `jSs == 'IN'`. tRNA is not in that list, so all 6 are dropped.
4. `jS:IN` requires the read to be **contained**: `readstart >= refstart &&
   readend <= refend` (`deal_with_singlemappers.sh`). A 151 nt read cannot be
   contained in a feature shorter than 151 nt.
5. In the v2 BED, the fraction of features shorter than 151 nt is:

   | biotype | features | under 151 nt |
   |---|---|---|
   | rRNA | 354 | **99.7%** |
   | miRNA | 2,206 | **99.2%** |
   | tRNA | 1,137 | **98.5%** |
   | snoRNA | 1,507 | **96.5%** |
   | snRNA | 1,381 | **84.6%** |
   | ProteinCoding | 455,326 | 29.5% |
   | lncRNA | 224,327 | 22.8% |

6. And the read populations differ accordingly: **92.5%** of FLASH-seq SE BED
   rows span ≥140 nt versus **43.2%** of VASA cell 005's; `jS:IN` rate 39.7%
   (FLASH-seq) versus 56.7% (VASA).

**Why this matters:** the comparison's second axis is *which RNA species are
detected*. The short non-poly-A species — tRNA, snoRNA, snRNA, miRNA — are
exactly where VASA is expected to win, and read length alone suppresses them on
the FLASH-seq side independently of any biology. VASA detects 90 tRNA isotypes
in its deposited table; FLASH-seq detects 0 here.

**Do not read that as "FLASH-seq detects no tRNA".** Two effects are
superimposed and this dry run cannot separate them: poly-A priming genuinely
depletes non-poly-A RNA, *and* the containment rule structurally discards long
reads over short features. What is proven is that the structural mechanism
operates — every tRNA-overlapping read present was discarded by containment, not
missing — not its magnitude relative to the biological depletion. Resolving that
needs a read-length-matched control. Flagged, quantified, not silently filtered.

### 5c. Do not use the annotation route for rRNA

The gene tables put rRNA at **0.80%** of reads. The rRNA leg of this comparison
is already measured by the bwa route at **3.50-6.44%** for FLASH-seq
(`res/flashseq/rrna_bwa.tsv`). The annotation route is low for two independent
reasons: 99.7% of rRNA features are under 151 nt (so §5b applies in full), and
the Ensembl annotation lacks the 47S unit at all — the finding already recorded
as "rRNA is not 0.8%" in `code/flashseq/README.md`. This dry run reproduces that
0.8% figure from the other direction. **Use `res/flashseq/rrna_bwa.tsv`.**

---

## 6. Changes needed before the ten-library run

In `code/flashseq_vasa/`, as forks — nothing in `a_Mapping/` is to be edited.

1. **Single-end input, or a fork that strips the mate suffix.** PE is a hard
   stop (§5a). Either map R1 only, or fork `deal_with_singlemappers.sh` to run
   `bamtobed` before the tag append, or to `sed 's|/[12];|;|'` the name column.
   Mapping R1 only is the smaller change and the one the dry run validated;
   it halves the usable bases.
2. **Trim before mapping.** 42-46% of reads were unmapped-too-short untrimmed,
   consistent with the 55.1% adapter read-through measured on ZHA8833A1. Reuse
   the cutadapt call in `code/flashseq/05_rrna_bwa.sh`, which already exists
   precisely so this is one recipe rather than two.
3. **Decide the read-length control for §5b, and state it in the report.** The
   options, in the repo's own idiom of not creating a second method: hard-trim
   FLASH-seq R1 to VASA's aligned-length distribution before mapping (matches
   what `05_rrna_bwa.sh` did for the rRNA leg — one measurement, both sides), or
   report the short-biotype axis as structurally non-comparable at 151 nt and
   restrict the species claim to long biotypes. **This is a decision for the
   user, not for the code.**
4. **`filt_unigenes='n'`, and choose the detection filter explicitly** (§3).
   Filter on `ReadCounts`; there is no UMI ceiling to filter on.
5. **`stranded='n'`** on every assignment call. `y` would halve every figure.
6. **Filename stems must contain `_cbc`** for `CELLID_FROM=f`, and the
   `cells/` prefix must be stripped from the column names downstream.
7. **Do not carry any `*.TranscriptCounts.tsv` forward** (§2.2).

Sizing from this slice: 604,853 reads took ~1 min of STAR and ~20 s of
assignment per arm. A full library is ~29 M reads, so ~48x — budget roughly an
hour per library for steps 4-5, and note that step 6's cost is superlinear in
BED size (the VASA side: 69 MB → 28 min, 254 MB → 1 h 47 m).

---

## 7. How to call it

```bash
# step 6
countTables_2pickle_cellsSpliced.py cells <SAMPLE> smartseq_noUMI f
# step 7
countTables_fromPickle.py <SAMPLE>.pickle.gz <SAMPLE> smartseq_noUMI n
```

with `deal_with_singlemappers.sh <bam> <v2.bed> n <samtools> <bedtools>` before
them (note the `n`), and per-library BAM stems containing the literal `_cbc`.

**Tables to use / not use, from the no-UMI path:**

| table | use it? |
|---|---|
| `*.ReadCounts.tsv` | **yes** — this is the quantification |
| `*.UFICounts.tsv` | only as a detection mask (0/1) |
| `*.TranscriptCounts.tsv` | **no** — identical to the mask; no abundance information |
| `*_spliced/_unspliced.ReadCounts.tsv` | yes, with the double-count caveat in §2.3 quantified |
