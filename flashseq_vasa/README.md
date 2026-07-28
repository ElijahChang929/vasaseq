# flashseq_vasa/ -- is the comparison actually a comparison?

This directory holds the comparability audit run BEFORE any VASA-vs-FLASH-seq
comparison is built, plus the scripts that produced every number in it. Nothing
here reimplements a measurement that already exists; where a number was already
measured (`res/flashseq/rrna_bwa.tsv`, `data/PM26037/out/logs/step3_report.txt`,
`data/PM26037/out/analysis/`) it is read, not recomputed.

    00_vasa_ribo_strand.sh   the one measurement that was MISSING (see below)
    01_read_vs_molecule.py   VASA's own read:UFI ratio, per cell and per biotype
    provenance.tsv           one row per shared component, and whether it matches
    read_vs_molecule.tsv     the per-cell / per-biotype ratios
    vasa_strand_check.tsv    stranded=y vs n over VASA's 16 surviving all-ribo BAMs

## Verdict: PARTIAL. One blocker, and it is not the one that was expected.

Twelve of eighteen components are byte-identical between the two runs. The rRNA
leg really is one measurement: the same `unique_rRNA_mouse.v2.fa`
(md5 `177913504bea6d1d0b992c018f84383b`, one file, one inode -- not two copies),
the same bwa index sidecars, the same `riboread-selection.py`
(md5 `f09e4ac4844ed181bf6adfee628c2183`), the same `ribo-bwamem.sh`, the same bwa
0.7.17-r1188 and samtools 1.11, the same conda prefix. Both runs postdate the
2026-07-26 20:34 reference rebuild, so neither used v1.

`riboread-selection.py` has an mtime (2026-07-27 09:57) that falls BETWEEN the
VASA run and the FLASH-seq run, which looks alarming and is not: its git blob is
`1ca5ff1c0fc4e6fcfddea887b6d018bdedeeaf89` in the working tree, in HEAD, and in
8bee6f6 (2026-07-16), and `git diff HEAD` is empty. No commit touches either
ribo script between the two runs. The mtime is cosmetic.

**The STAR index is ONE index, not two.** This was the expected blocker and it
is not one. `star_index_151_r116/genomeParameters.txt` is nf-core's own
`genomeGenerate` command line (`/opt/conda/bin/STAR-avx2 ... --sjdbGTFfile
Mus_musculus.GRCm39.116.filtered.gtf --sjdbOverhang 150`), and VASA's own STAR
`Log.txt` names that directory as `genomeDir`. nf-core built it, it was moved
into the reference tree, and VASA then used it. 78,348 genes, 481,956
transcripts, 61 contigs, 562,855 GTF junctions, versionGenome 2.7.4a. Biotype
assignment cannot diverge from a difference that does not exist.

## The blocker: stranded=y/n is NOT the symmetric factor of two the README implies

`05_rrna_bwa.sh` was careful here -- it ran FLASH-seq with `n` and re-selected
the SAME merged BAM with `y`, so its side is measured. The forward strand carries
**49.1-50.5 %** of FLASH-seq ribosomal reads and `y/n = 0.491-0.509` across all
ten libraries: genuinely unstranded, and `n` is the right flag. That part holds.

What had never been computed is the mirror image. VASA only ever ran `y`, so
"what would VASA say under `n`" was an assumption. All 16 `.nsorted.all-ribo.bam`
survive, so it is measurable. `00_vasa_ribo_strand.sh` counts the two predicates
`riboread-selection.py` itself uses over the same name-sorted BAM, and
reproduces its logged counts exactly (`ribo_stranded_y == log_nmapped_sum` for
all 16 barcodes, and `n_groups == log_nreads`; that equality is the script's own
self-check).

    VASA, stranded=y  21.39 %  (19,282,729 / 90,137,383)   <- the published figure
    VASA, stranded=n  28.10 %  (25,328,057 / 90,137,383)   <- never computed before

The forward strand carries **76.1 %** of VASA's ribosomal reads (57.7-83.8 %
across the 12 real cells), not ~50 %. So `y` costs VASA 23.9 % of its ribosomal
reads, where the same flag costs FLASH-seq 50.6 %. **The two published
percentages are not on the same footing**, and the ratio between them depends on
a flag choice by roughly a factor of two:

| both sides on | VASA | FLASH-seq (pooled, trimmed) | ratio |
|---|---|---|---|
| as published (VASA `y`, FS `n`) | 21.39 % | 4.78 % | **4.5x** |
| both `n` | 28.10 % | 4.78 % | **5.9x** |
| both `y` | 21.39 % | 2.39 % | **9.0x** |

**The published 4.3x does not reproduce, and I could not find a derivation that
gives it.** Pooled-against-pooled is the right denominator against VASA's own
`ALL` row -- sum(ribo)/sum(reads_in) over the trimmed arm is 4.780 %, so
21.39/4.780 = **4.48x**. Every other candidate lands in the same band: trimmed
mean of per-library rates 4.54x, trimmed median 4.53x, raw pooled 4.48x, raw mean
4.55x, the 30 pg pair alone 4.46x, excluding A8 4.51x. Against VASA's 12-real-cell
21.57 % the band is 4.49-4.59x. To land on 4.30x the FLASH-seq rate would have to
be 4.974 %, which is not any single library (range 3.496-6.439 %) nor any central
measure of them.

So the gap is ~4 % relative and it is unexplained. It changes no conclusion --
FLASH-seq carries several-fold less rRNA than VASA under every derivation -- but
**quote 4.5x (pooled, and state the flag pairing), not 4.3x**, until whoever wrote
the 4.3x can say where it came from. Two things it is NOT: it is not the
`stranded` pairing (that moves the ratio to 5.9x or 9.0x, not to 4.3x), and it is
not the arithmetic of the pooling, which reproduces exactly.

**Recommendation, not a decision:** report `both n` (5.9x). `n` counts every
read that aligns to rRNA regardless of orientation, which is the same question
on both protocols, and it does not require VASA to be perfectly stranded --
which, at 76.1 % forward, it measurably is not. The 76.1 % is itself worth a
look: for a stranded protocol it is lower than expected, and commit a9d0c46
already documented one mechanism (reverse-strand reads written flipped). Whether
that is the same 23.9 % is not established here. Quote whichever pairing you
choose, but quote the flag with the number, and do not mix them.

## Read-vs-molecule: 2.6x, and it is expression-dependent

VASA counts molecules; `smartseq_noUMI` will count reads. Measured on the 12 real
cells (222,412 x 12 total tables):

* **read:UFI = 2.13-2.96 per cell, median 2.57**, pooled 2.69. Cell 002 lowest,
  011 highest. read:transcript median 2.42; UFI:transcript 0.942 (the collision
  correction adds ~6 %).
* **Not uniform across biotypes.** Spearman(UFIs, read:UFI) = **+0.498** over
  35,758 entries with >=10 UFIs; by expression decile the ratio climbs
  **1.69 -> 2.91 (1.73x)**. Highly-expressed entries are more duplicated, as
  expected.
* **Direction of the bias -- set by per-class duplication, NOT by abundance.**
  The sign of each class's shift is decided by its own read:UFI against the
  library's pooled 2.689: `share_ratio = (class read:UFI) / 2.689` exactly
  (verified to 1e-3 on every row). A class more duplicated than the library
  average takes a LARGER share of reads than of molecules; a less duplicated one
  takes a smaller share. Abundance does not set the sign -- the small-RNA classes
  are rare AND heavily duplicated, so they move the same way as the abundant ones.

  | class | % of reads | % of UFIs | read:UFI | share ratio | effect of using UFI-shares |
  |---|---|---|---|---|---|
  | MiscRna | 0.960 | 0.450 | 5.73 | 2.131 | **understates by 113 %** |
  | rRNA (residual) | 0.486 | 0.273 | 4.79 | 1.780 | understates by 78 % |
  | snRNA | 2.972 | 2.229 | 3.58 | 1.333 | understates by 33 % |
  | ribozyme | 0.157 | 0.131 | 3.24 | 1.205 | understates by 20 % |
  | ProteinCoding | 82.646 | 81.734 | 2.72 | 1.011 | understates by 1 % |
  | snoRNA | 1.525 | 1.672 | 2.45 | 0.912 | overstates by 10 % |
  | lncRNA | 6.970 | 7.910 | 2.37 | 0.881 | overstates by 13 % |
  | ProcessedPseudogene | 0.699 | 0.939 | 2.00 | 0.745 | overstates by 34 % |
  | miRNA | 0.076 | 0.137 | 1.50 | 0.557 | **overstates by 80 %** |

  So quoting VASA UFI-shares against FLASH-seq read-shares **understates the
  heavily-duplicated classes (small nuclear/misc RNA, residual rRNA) by 20-113 %
  relative, and overstates the lightly-duplicated ones (miRNA, pseudogenes,
  lncRNA) by 10-80 %**. Protein-coding is almost exactly unbiased at 1 %, because
  it dominates the pooled average it is being compared against. The classes a
  total-RNA-vs-poly-A comparison is actually about are the ones where this is
  worst, and it runs in BOTH directions depending on the class -- so it cannot be
  absorbed into a single correction factor.

  Note that the three worst-affected small-RNA classes are the RESIDUAL after
  step 7 already dropped 8 UMI-ceiling genes (Rmrp, Rnu1a1/1b6/2.10, Rn7sk,
  Rn7s1/2, Snord3b1-4 -- see `analysis/manifest.json`). Their read:UFI is high
  because short, highly-expressed loci saturate a 6 nt UMI space; the ceiling
  filter removed the saturated cases, not the tendency.
* **Gene detection is the exception: it is immune.** Genes detected on reads and
  on UFIs are **identical in all 12 cells** (difference 0, not merely small):
  a gene seen at all is seen on >=1 read and >=1 UFI. 225,716 cell-gene entries
  sit at exactly 1 read, 274,258 at exactly 1 UFI.

## Which VASA column each comparison must use

| comparison | column | why |
|---|---|---|
| (i) rRNA fraction | **neither** -- use the bwa stage | The rRNA % lives in step 3, upstream of any counting. Do NOT use the rRNA *biotype* rows: they are the residual that survived depletion (0.506 % of reads / 0.304 % of UFIs), a different denominator, and 87.1 % of it is one entry (`ENSMUSG00000119584_Rn18s.rs5_rRNA`). Report with the stranded flag stated. |
| (ii) biotype composition | **ReadCounts** (headline) + **TranscriptCounts** (VASA's biology) | Like against like: both sides reads. Carrying both is not hedging -- the gap between them IS the measured bias, which runs from 2.13x (MiscRna) down to 0.56x (miRNA) and so cannot be corrected with one factor. |
| (iii) gene detection | **UFICounts** or ReadCounts -- measurably identical | Difference is 0 in all 12 cells, so the choice is free. Use UFICounts for consistency with VASA's own tables. Depth-match instead: that is the real confound (VASA spans 6.2x of depth and 32,591-65,137 genes). |

TranscriptCounts is the right column for VASA's own biology (it is UFIs with the
binomial collision correction, +5.8 % median) but it is the WRONG comparator
against FLASH-seq: it corrects for a saturation FLASH-seq has no equivalent of.

## Two further cautions carried forward

The annotation is the same Ensembl 116 release on both sides (source GTF
md5 `0f9ab91d5ed1be2c7538589d6950f3af`) but not the same OBJECT: VASA uses
`Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed` (719,409 rows, one entry
per gene, merged-exon union + explicit introns + 1,137 GtRNAdb tRNA rows), while
nf-core used a per-transcript filtered GTF. Biotype labels agree because both
read `gene_biotype` from the same GTF. This only bites if nf-core/RSEM is used as
the comparator, which decision 1 already rules out.

`results/genome/Mus_musculus.GRCm39.116.filtered.gtf` as published is a
**truncated copy** -- exactly 167,837,696 bytes (2561 x 64 KiB), last line cut
mid-attribute, chromosome 1 only, 4,742 genes. The run itself was NOT affected:
RSEM's merged gene matrix has all 78,348 genes, `genome.transcripts.fa` and the
filtered BED both have 481,956 transcripts across 38 contigs, and the STAR index
has 78,348 genes. The truncation is in the copy on disk, not in what was
computed. Worth fixing so nobody rebuilds from it.
