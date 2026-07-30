# Annotation-release effect: Ensembl 99/GRCm38 vs Ensembl 116/GRCm39

**Track**: control track for the three-way comparison.
**Question**: how much of any published-plate-vs-own-plate difference is caused by
the annotation release rather than by protocol or pipeline?
**Answer in one line**: **36-58% of the biotype-composition gap** and
**essentially all of the raw gene-detection gap** could be release rather than
protocol. Both bounds are measured, not argued; the interval is the honest answer
and the reason it is an interval is explained below.

Every number here is asserted against its source table by
`verify_annotation_release_effect.py`. Nothing is transcribed.

---

## The two headline sentences, for citation by other tracks

> **Biotype composition.** The raw composition gap between the published plate
> (Ensembl 99) and the own plate (Ensembl 116), measured as total variation
> distance on ReadCounts over single-gene mouse Ensembl rows, is **9.97 pp**.
> Forcing the two plates onto a common gene universe removes between **3.57 pp
> (35.8%)** and **5.83 pp (58.5%)** of it, depending on whether genes re-issued
> under a new Ensembl ID are treated as the same gene. So **at least a third and
> possibly nearly three-fifths of the composition gap is annotation release, not
> protocol** — and the residual **4.14-6.40 pp** is the most that protocol,
> biology and depth can jointly explain.

> **Gene detection.** The own plate detects 51,868 single-gene entries and the
> published plate 32,430 mouse ones, but **34.81%** of the own plate's detected
> genes have IDs that do not exist in Ensembl 99 at all, against **2.30%** in the
> other direction. A raw "own plate detects more genes" comparison is therefore
> **uninterpretable as a chemistry statement**: the newer annotation supplies
> 1.41x as many mouse gene IDs (78,348 vs 55,471), and 3.30x as many lncRNA genes
> (32,889 vs 9,959). Detection must be compared on the shared gene set or not at
> all.

---

## 1. What was and was not measured

**Measured.** The gene universes, biotype labels and feature lengths of the two
annotation BEDs *the pipeline actually consumed*, and the biotype composition of
both plates recomputed under three different definitions of "the same gene".

**Not measured, and no track should imply otherwise.**

- **The GRCm38 -> GRCm39 assembly change is not measured here.** Both plates were
  mapped to their own assembly and nothing in this track re-maps a read or lifts
  over a coordinate. A gene counts as "shared" when its ID is in both BEDs, even
  if the underlying sequence moved or changed. Every figure here is therefore a
  **lower bound on the total release effect**, because the assembly component is
  invisible to it.
- **The FLASH-seq arm has no release term.** It shares Ensembl 116/GRCm39 with the
  own plate, so for own-vs-FLASH the release effect is zero by construction and
  this document does not apply.
- **Nothing here is a protocol measurement.** Both plates are VASA. The residual
  after harmonisation bounds the protocol+biology+depth term from above; it does
  not decompose it.

**Definitions used throughout.** ReadCounts on all sides (Rule 4: reads are the
only unit both protocols measure, and this track must be comparable with the
cross-protocol ones). Own plate: 12 real cells, blanks `001/014/015/016`
excluded, UMI-ceiling genes **retained** (Rule 3 — on a read-based composition
they carry real reads). Published plate: the **173 mESC wells** called by
`classify_fig1d` on the v3 UFICounts with MIN_UFI=7500, the rule copied verbatim
from `vp_common.py` (178 HEK293T, 2 mixed, 31 discarded). Composition is over
**single-gene mouse Ensembl rows only**; combination rows and tRNAscan rows are
excluded on both sides, for reasons given in §6.

---

## 2. Gene-universe overlap

| | Ensembl 99 (mouse) | Ensembl 116 | shared | 99-only | 116-only |
|---|---|---|---|---|---|
| all genes | 55,471 | 78,348 | 51,765 | 3,706 | 26,583 |
| % of that release | — | — | 93.32% of E99, 66.07% of E116 | 6.68% | 33.93% |

Jaccard index **0.6309**. Ensembl 116 carries **1.41x** as many mouse gene IDs.

**Protein-coding is the stable part of the annotation:**

| | E99 | E116 | shared, same class | 99-only | 116-only | ID present but reclassified |
|---|---|---|---|---|---|---|
| ProteinCoding | 21,933 | 21,818 | 21,513 | 199 | 158 | 221 (E99→other), 147 (E116→other) |

**98.09%** of Ensembl 99's protein-coding genes are still protein-coding in 116.
A protein-coding-restricted comparison is therefore only mildly release-exposed.

**lncRNA is where the universe explodes:** 9,959 → 32,889 genes (**3.30x**), with
**23,236** lncRNA genes present only in Ensembl 116. Any lncRNA comparison
between these two plates is dominated by the release.

**Vocabulary churn.** Ensembl 99 has 38 biotype tokens, 116 has 37.
`PolymorphicPseudogene` was **retired** — all 88 of its shared-gene members were
relabelled, 100% churn — and no new class appeared. A composition table that
lists `PolymorphicPseudogene` for the published plate and 0 for the own plate is
reporting a **vocabulary change, not a biological absence**.

---

## 3. Biotype-label churn among genes present in both releases

Of the 51,765 shared genes:

- **1,761 (3.40%)** changed biotype label;
- **379 (0.73%)** changed even after lumping into coarse classes
  (ProteinCoding / lncRNA / sncRNA / pseudogene / TEC / IgTr);
- **2,911 (5.62%)** changed gene **symbol** while keeping the same ID — which
  silently breaks any symbol-keyed join between the two datasets.

The churn is overwhelmingly **inside the pseudogene family**, and it is mostly one
transition: **`UnprocessedPseudogene` → `TranscribedUnprocessedPseudogene`, 671
genes**, the largest single transition in the matrix.
`UnprocessedPseudogene` lost 25.68% of its members; `ProcessedPseudogene` lost
687 genes (6.88%). `ProteinCoding` lost only 1.02%.

**Consequence.** A fine-grained biotype comparison between these releases is
contaminated at the **3.4%** level, but a *coarse* comparison is contaminated at
only **0.73%** — because most churn is pseudogene-subclass reshuffling that the
coarse lumping absorbs. **Use coarse classes for any cross-release biotype claim.**

---

## 4. Feature-length churn — and a file-format trap that had to be fixed first

**The trap.** The two BEDs do not share a coordinate convention. The Ensembl 99
mixed BED is **1-based inclusive** (modal inter-feature gap 1) and the Ensembl 116
BED is **0-based half-open** (modal gap 0). Taking `end - start` in both — the
obvious thing to do — makes **every** Ensembl 99 feature read exactly 1 bp shorter
than its Ensembl 116 counterpart.

The correction is verified on the 19,092 genes that are single-feature in both
releases, where an unchanged model must give an identical length:

| | lengths identical |
|---|---|
| raw `end - start` | **0.026%** (5 of 19,092) |
| convention-corrected | **98.04%** (18,718 of 19,092) |

`ENSMUSG00000048040`: E99 raw 1531, E116 raw 1532; corrected, both 1532.
**Uncorrected, this artefact would have been reported as "Ensembl 116 lengthened
every gene model by 1 bp".** It is detected from the data, not assumed, and the
detection is asserted in the analysis script.

**With the convention corrected, the aggregate looks stable — but the aggregate is
misleading, and this is the one place in this track where the headline number hides
the finding.** Across all 51,765 shared genes: **55.64%** have byte-identical
summed exonic length, **71.83%** have an unchanged exon count, and the median
log2(E116/E99) exonic ratio is **0.0**. Read alone, that says "nothing moved".

Split by biotype, two very different regimes appear:

| biotype (shared genes) | n | median exonic E99 | median exonic E116 | median log2 ratio | % byte-identical |
|---|---|---|---|---|---|
| ProteinCoding | 21,734 | 3,532.5 bp | 4,744.0 bp | **+0.4253** | 24.43% |
| lncRNA | 9,584 | 1,045.0 bp | 1,896.5 bp | **+0.8592** | 38.16% |
| snoRNA | 1,065 | 129.0 bp | 129.0 bp | 0.0000 | **100.00%** |
| snRNA | 789 | 107.0 bp | 107.0 bp | 0.0000 | **100.00%** |
| scaRNA | 34 | 136.0 bp | 136.0 bp | 0.0000 | **100.00%** |
| miRNA | 1,124 | 80.5 bp | 80.5 bp | 0.0000 | 99.91% |
| MiscRna | 32 | 238.0 bp | 238.0 bp | 0.0000 | 96.88% |
| rRNA | 13 | 119.0 bp | 119.0 bp | 0.0000 | 92.31% |

**The classes that carry the reads grew; the classes that carry the counts of
short RNAs did not move at all.** A shared protein-coding gene's median exonic
length rose by **34%** (2^0.4253 = 1.34) and a shared lncRNA's by **81%**
(2^0.8592 = 1.81). The overall median of 0.0 is an artefact of the many thousands
of short, frozen ncRNA genes outvoting them.

**Two consequences, in opposite directions.**

1. **Longer models capture more reads.** A gene whose annotated exonic footprint
   grew by a third has more sequence on which a read can be contained, so part of
   the own plate's protein-coding read count reflects a bigger target, not more
   transcript. This is an *additional* release effect on composition beyond the
   gene-universe effect quantified in §5 — and §5 does **not** capture it, because
   §5 harmonises *which genes exist*, not *how large each one is*. The
   composition-gap bracket in §5 is therefore, in this respect too, a **lower
   bound**.
2. **Trap 8 is a constant here, not a release effect.** VASA's `jS:IN` rule
   requires a read to be contained in a feature, so at 151 nt the short classes
   are structurally near-invisible: for shared genes, **96.62%** of snoRNA exon
   features, **80.23%** of snRNA and **100.00%** of miRNA are shorter than 151 nt
   — and those percentages are **identical to three decimal places in both
   releases**. The suppression cannot explain any plate-to-plate difference in
   short-RNA detection. It remains a real confound for any short-RNA claim; it is
   just not a release confound.

---

## 5. The decisive test: composition restricted to a common gene universe

This is the number the task asked for. Both plates' biotype composition is
recomputed on only those genes both releases contain, and the gap that closes is
the release effect.

**Why there are two bounds rather than one number.** An ID-keyed restriction
treats a gene re-issued under a new Ensembl stable ID as two different genes and
discards **both** copies. That happens here, demonstrably:

- **Rn7sk** carries **669,205** reads on the own plate as `ENSMUSG00002076161`
  and **4,199** reads on the published plate as `ENSMUSG00000065037`. Same gene,
  different ID, discarded from both sides.
- **Ptprg** (`ENSMUSG00000121513` vs `ENSMUSG00000021745`), **Gm25360**,
  **Gm22988** and **Rn7s2** do the same.

So an ID-only restriction **overstates** the release effect. Bridging additionally
on gene symbol **understates** it, because a placeholder symbol (`Gm*`, `*Rik`)
can be reused for a genuinely different locus. The truth is bracketed:

| gene universe | restricted TVD | release-attributable | as % of raw gap |
|---|---|---|---|
| raw (no restriction) | 9.97 pp | — | — |
| ID only — **upper bound** | 4.14 pp | **5.83 pp** | **58.47%** |
| ID + informative symbol | 5.66 pp | 4.31 pp | 43.19% |
| ID + any symbol — **lower bound** | 6.40 pp | **3.57 pp** | **35.78%** |

(923 symbol-bridged 1:1 pairs; 200 of them on informative symbols. Coarse and
fine biotype levels agree to within 0.01 pp, so the choice of granularity does not
move this.)

**Per-class, ID-only universe (the numbers other tracks will cite):**

| coarse class | published | own | raw gap | published (shared) | own (shared) | gap (shared) | gap closed |
|---|---|---|---|---|---|---|---|
| ProteinCoding | 95.41% | 85.83% | **−9.58 pp** | 96.01% | 92.18% | **−3.84 pp** | 5.74 pp |
| sncRNA | 0.97% | 6.37% | **+5.40 pp** | 0.61% | 3.36% | **+2.75 pp** | 2.65 pp |
| lncRNA | 2.49% | 7.06% | **+4.57 pp** | 2.28% | 3.67% | **+1.39 pp** | 3.18 pp |
| pseudogene | 1.13% | 0.74% | −0.39 pp | 1.09% | 0.79% | −0.30 pp | 0.09 pp |

**Read this as:** the own plate's apparent 9.58 pp protein-coding deficit shrinks
to 3.84 pp once the gene universe is harmonised — **60% of it was the newer
annotation offering more non-coding places for reads to land**. Likewise the
own plate's apparent sncRNA excess halves and its lncRNA excess falls by ~70%.
The residuals are real and are what the protocol/depth tracks must explain.

**One locus carries a quarter of it.** Dropping `Rn7sk` from both plates lowers the
raw gap from 9.97 to **8.40 pp** and the ID-only release term from 5.83 to
**4.26 pp** — so **1.57 pp of the 5.83 pp, i.e. 27%, is one 7SK RNA locus** that
both releases annotate under different IDs. It alone is **25.62%** of the own
plate's restriction loss. **The release effect on composition is not diffuse; it is
a few-locus effect with a long tail**, and any track quoting the 5.83 pp should say
so.

**Asymmetry, explained.** Restriction costs the own plate **7.02%** of its
single-gene reads but the published plate only **0.94%** — a **7.5x** difference,
and it is not incidental. The own plate's loss sits on **lncRNA (3.64 pp)** and
**sncRNA (3.25 pp)**, i.e. exactly the classes where Ensembl 116 added 23,236
lncRNA genes and re-issued the sncRNA IDs. Protein-coding contributes 0.12 pp.

**Sensitivity to counting unit.** On TranscriptCounts the own-plate composition
shift is larger (TVD 7.72 pp coarse vs 6.40 pp on ReadCounts; ProteinCoding
+7.64 pp vs +6.35 pp), because UMI collision correction upweights the abundant
short non-coding genes that the restriction removes. ReadCounts is the reported
unit per Rule 4; the TranscriptCounts columns are in the table as a sensitivity
check, and they move the release effect **up**, not down.

---

## 6. Why the sncRNA gene sets look bizarre — and what it means

Several sncRNA classes have **identical gene counts in both releases** but nearly
disjoint ID sets:

| class | E99 | E116 | shared IDs | length-multiset match | shuffled control |
|---|---|---|---|---|---|
| snRNA | 1,385 | 1,381 | 789 | **100.00%** | 72.25% |
| snoRNA | 1,507 | 1,507 | 1,065 | **100.00%** | 56.97% |
| scaRNA | 51 | 51 | 34 | **100.00%** | 51.18% |
| miRNA | 2,207 | 2,206 | 1,124 | 99.91% | 77.78% |
| MiscRna | 562 | 562 | 32 | 99.81% | 94.92% |
| rRNA | 354 | 354 | 13 | 99.41% | 97.26% |
| ProteinCoding | 199 | 158 | — | 15.82% | 1.58% |

The "length-multiset match" is the fraction of the smaller only-set whose exonic
lengths are reproduced exactly in the other release's only-set; the control is the
same statistic against 20 size-matched draws from that biotype's whole length pool
(seed 20260730). The ID blocks are disjoint: only-E99 sncRNA IDs span
**64,393-106,670**, only-E116 span **118,674-2,076,992**.

**Interpretation, with its limit stated.** Identical class totals, disjoint and
later ID blocks, and length multisets that match at 100% against controls of
51-72% are together consistent with **the same loci re-issued under new stable
IDs**, not with genes being added or removed. This is **evidence about
length-distribution identity, not proof of locus identity** — no liftover was
performed, and the two releases are on different assemblies, so per-locus identity
was not established. For `MiscRna` and `rRNA` the control is already at 92-97%
(their length distributions are narrow), so those two rows carry much less
information than the snRNA/snoRNA/scaRNA rows.

The contrast with `ProteinCoding` (15.82% match, control 1.58%) is informative in
the other direction: the protein-coding only-sets are **genuinely different
genes**, not re-issued IDs.

**For trap 9.** The `rRNA` biotype behaves exactly as the conventions warn. On the
own plate the second-largest single contributor to restriction loss is
`ENSMUSG00000119584_Rn18s.rs5_rRNA` with 197,472 reads — the 18S relic locus. This
track does not measure rRNA content and nothing here should be read as doing so.

**Excluded row classes, and why.** tRNAscan rows carry no Ensembl ID, so they can
never be "shared" in the ID sense; leaving them in the unrestricted composition
and dropping them from the restricted one would book their entire share as a
release effect. They are also not release-comparable — the E99 BED carries 1,758
tRNAscan rows for human+mouse pooled, the E116 BED 1,137 for mouse alone — so they
are held out of every composition on both sides and reported separately (own plate
0.025% of reads, published 0.0016%). Combination (multi-gene) rows are excluded
too, and they are **large**: 30.10% of own-plate reads, of which **82.85%** name at
least one gene that exists only in Ensembl 116. A combination-row-inclusive
comparison between these releases would be **more** release-contaminated than
anything reported above, not less.

---

## 7. What this licenses and forbids

**Licensed:**

- Quoting **3.57-5.83 pp (35.8-58.5%)** as the annotation-release share of the
  published-vs-own biotype-composition gap, with the bracket, on ReadCounts,
  coarse classes — noting it is a **lower bound**, because it harmonises which
  genes exist but not how large each one is (§4: shared protein-coding models grew
  34%, lncRNA 81%).
- Quoting **4.14-6.40 pp** as the upper limit on what protocol, biology and depth
  can jointly explain in that gap.
- Comparing **protein-coding** biotype composition across releases with only mild
  correction (98.09% class stability), preferably still on the shared set.
- Stating that the own plate's **gene-detection advantage is not established**:
  34.81% of its detected genes are absent from Ensembl 99 by ID.

**Forbidden:**

- Presenting the raw **9.97 pp** composition gap, or the raw ProteinCoding
  **−9.58 pp**, as a protocol difference.
- Any **lncRNA** comparison between these two plates without harmonising the gene
  set (3.30x universe expansion).
- Any **short-RNA** (tRNA / snoRNA / miRNA) comparison between the plates without
  a read-length-matched control — feature lengths are unchanged between releases,
  so this is trap 8, not a release effect.
- Treating the **5.83 pp** figure as diffuse: 27% of it is `Rn7sk` alone.
- Reading anything here as an **rRNA** measurement (trap 9) or as covering the
  **GRCm38→GRCm39 assembly** change (not measured).
- Comparing **per-gene expression** (not just composition) across the two plates
  without per-gene length normalisation, or without restricting to the 24.43% of
  shared protein-coding genes whose exonic length is byte-identical. A shared gene
  is not the same measurement target in the two releases if its annotated footprint
  grew by a third.

---

## 8. Files

Written to `res/threeway/`:

| file | contents |
|---|---|
| `annotation_release_effect.tsv` | headline table, 114 metrics, each with its denominator |
| `annotation_release_bracket.tsv` | the upper/lower bound bracket and Rn7sk leave-one-out |
| `annotation_release_bracket_composition.tsv` | per-class composition under all three universes |
| `annotation_gene_universe_biotype.tsv` / `_coarse.tsv` | universe overlap per class |
| `annotation_biotype_churn_matrix.tsv` | every 99→116 biotype transition with counts |
| `annotation_biotype_churn_ledger.tsv` | per-class gained/lost/kept ledger |
| `annotation_feature_length.tsv` / `_paired.tsv` | length churn, per class and per gene |
| `annotation_composition_shift.tsv` | all-genes vs shared-genes composition, both plates |
| `annotation_crossplate_gap.tsv` | the per-class gap decomposition |
| `annotation_release_mechanism.tsv` / `_asymmetry.tsv` | ID-series and loss-attribution diagnostics |
| `annotation_release_reid.tsv` | length-multiset re-ID test with shuffle controls |
| `annotation_release_loss_top.tsv` | the genes carrying the restriction loss |
| `annotation_release_effect_provenance.json` | exact inputs, wells used, versions |

Scripts in `code/flashseq_vasa/threeway/`: `annotation_release_effect.py` (main),
`annotation_release_effect_precheck.py` (Rule 2 precheck),
`annotation_release_diagnose.py`, `annotation_release_reid_test.py`,
`annotation_release_bridge.py` (the bracket),
`verify_annotation_release_effect.py` (asserts every number in this note).

Peak RSS 0.44 GB, wall 41 s for the main analysis. Environment `envs/vasa`:
python 3.10.20, pandas 2.3.3, numpy 2.2.6.
