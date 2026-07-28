# What nf-core/RSEM says about FLASH-seq composition, on its own terms

Run RN26038, `ZHA8833A1..A10`, nf-core/rnaseq `--aligner star_rsem`, Ensembl 116 /
GRCm39. Scripts `02_nfcore_crosscheck.py`, `04_nfcore_mechanisms.py`,
`05_nfcore_final.py`, `06_rrna_locus_units.py`, `07_vasa_comparator.py`. Tables in
`$W/res/flashseq_vasa/`.

Purpose: establish the nf-core composition **independently of the VASA route**, so
that a FLASH-seq vs VASA difference can be split into a *quantifier* part and a
*chemistry* part. This is a cross-check, not the comparator — decision 1 in
`provenance.tsv` already settled that.

Two premises I was given turned out to be wrong. Both are corrected below on
measured evidence, and both matter for how the numbers are quoted.

---

## Correction 1 — nf-core's rRNA biotype figure is not a 5S measurement. It is one 18S relic locus.

The brief said `gene_biotype "rRNA"` in GRCm39 "holds ~354 genes that are
essentially all `n-R5s*` plus one `Rn18s-rs5` relic". The gene *count* is roughly
that shape; the **signal** is not.

By count, the 354 rRNA-biotype genes are 279 `Gm*` (median 113 bp), 72 `n-R5s*`
(5S, median 119 bp), 2 `Rn5s*`, and 1 `Rn18s-rs5` (1,849 bp). There is indeed no
`Rn45s`/`Rn28s`/`Rn5-8s`.

By signal, `Rn18s-rs5` alone carries **99.979–99.999 %** of the RSEM
rRNA-biotype `expected_count` in all ten libraries, and `n-R5s*` genes carry
**0.0000–0.0061 %**. Only 2–9 of the 354 genes are nonzero in any library.

Confirmed independently against the featureCounts side by indexed region query on
the STAR BAMs (`samtools view -c -F 0x100 -q 255`, i.e. unique primary):

| library | locus | length | unique reads | % of featureCounts `rRNA` row |
|---|---|---|---|---|
| ZHA8833A9  | `Rn18s-rs5` | 1,849 bp | 450,934 | **99.933 %** |
| ZHA8833A9  | 353 rRNA loci < 400 bp, pooled | ~116 bp | 11,055 | 2.450 % |
| ZHA8833A10 | `Rn18s-rs5` | 1,849 bp | 225,319 | **99.947 %** |
| ZHA8833A10 | 353 rRNA loci < 400 bp, pooled | ~116 bp | 7,027 | 3.117 % |

The two columns sum past 100 % because a region query counts any overlap while
featureCounts `-B -C` additionally requires exon overlap, both mates assigned and
no biotype ambiguity — so a region query is an **upper bound**. Reconciliation:
summed locus reads / featureCounts row = **1.0238** (A9 rRNA), **1.0060** (A9
Mt_rRNA), **1.0306** (A10 rRNA), **1.0137** (A10 Mt_rRNA). All slightly above 1,
which is the expected direction.

Mechanism: RSEM's `effective_length`. A transcript much shorter than the fragment
length cannot generate the observed fragments and so receives almost nothing. The
implied mean fragment length is **195.1 nt** (median `length − effective_length`
over 18,779 protein-coding genes > 1 kb, +1). For A9:

| biotype | n genes | median length | median eff_length | eff/length | % genes eff < 10 | % genes shorter than 195 nt |
|---|---|---|---|---|---|---|
| protein_coding | 21,818 | 2,148.9 | 1,955.2 | 0.910 | 0.16 % | 0.3 % |
| rRNA | 354 | 116.0 | 12.3 | 0.106 | 27.40 % | 99.7 % |
| snRNA | 1,381 | 107.0 | 9.7 | 0.090 | 72.19 % | 99.9 % |
| snoRNA | 1,507 | 127.0 | 16.0 | 0.126 | 32.05 % | 97.7 % |
| miRNA | 2,206 | 89.5 | 5.4 | 0.061 | 72.48 % | 99.8 % |
| misc_RNA | 562 | 265.5 | 96.2 | 0.362 | 13.52 % | 31.5 % |
| Mt_tRNA | 22 | 68.5 | 2.0 | 0.029 | 100.00 % | 100 % |

`Rn18s-rs5` is the only rRNA-biotype gene with an adequate effective length
(1,654.9 at 1,849 bp), so it takes essentially the whole class.

**How to quote it.** Not "rRNA", and not "5S". Write *"nf-core annotation-route
rRNA = a single 18S relic locus (`Rn18s-rs5`, chr17:40,157,244–40,159,092)"*. The
bwa route reads 3.50–6.44 % because it sees the whole 47S unit; the annotation
route has no 5'ETS and no 28S in the primary assembly at all, and the bwa route
measures those at 37.3–46.9 % and 24.3–31.2 % of ribosomal reads respectively.
Consistency check: A9's bwa route gives ribo 5.256 % of reads with 18S = 9.3 % of
those, i.e. 0.489 % of all reads from 18S alone, against 0.698 % (RSEM) and
1.066 % (featureCounts). Same order of magnitude — which is as much as this
comparison supports, since the annotated 18S gene also collects reads from the
flanking 5'ETS/ITS1 that the bwa route assigns to those segments.

## Correction 2 — `results/genome/Mus_musculus.GRCm39.116.filtered.gtf` is a truncated file. Do not use it.

It is **167,837,696 bytes = exactly 160.0625 MiB = 2,561 × 64 KiB**, contains
**zero `gene` feature lines**, covers **only contig 1**, and ends mid-record
without a trailing newline. A block-aligned size with a mid-record cut is a
truncated write, not a filter product.

The GTF used at runtime was complete. Evidence: the featureCounts biotype tables
span the full Ensembl vocabulary including `Mt_rRNA` (MT contig) and `Mt_tRNA`;
the `Chr` field of the `rRNA` row lists all 21 assembled contigs
(1–19, X, Y); and `params.gtf` points at
`reference/genomes/mus_musculus/GRCm39/annotation/release-116/gtf/Mus_musculus.GRCm39.116.gtf`
with **md5 `0f9ab91d5ed1be2c7538589d6950f3af`**, identical to the source GTF that
VASA's v2 BED was built from (the md5 `provenance.tsv` records).

Consequence: every biotype join here uses the **source** GTF, and
`results/genome/*.filtered.gtf` should be treated as a damaged artefact of the
`--save_reference` copy. Script `02_` initially reported "0 gene rows" for it and
that looked like a parse bug; it was the file.

---

## Step 1 — nf-core composition

Denominators, verified before use rather than assumed:

- The featureCounts biotype table sums **exactly** to `Assigned` in all ten
  libraries (difference 0 for every one), so its denominator is assigned
  fragments only.
- nf-core's own `percent_rRNA` is reproduced to 1.1e-16 as
  `rRNA / sum(biotypes)` — rRNA only, `Mt_rRNA` **not** included.
- The table is in **reads (mates)**, not fragments: measured
  `reads / fragments = 1.988–2.003` at every rRNA locus. This matters — comparing
  a fragment count against it halves the answer, which is how the 18S locus first
  appeared to be "50 %" of the rRNA row.
- featureCounts discards multireads (`Unassigned_MultiMapping`) and
  biotype-ambiguous fragments (`Unassigned_Ambiguity`). Multireads are
  **31.50–43.48 %** of alignable reads (A8, the excluded library, is the 43.48 %).
- RSEM reports fractional `expected_count` over 78,348 genes; pooled total
  192,583,929. Every RSEM gene got a biotype from the source GTF (0 unannotated).

Composition, per the user's denominator decision (rRNA as % of all reads,
everything else as % of the non-rRNA remainder). Range is across all ten
libraries; the 30 pg rung (A9+A10) is the single-cell-equivalent comparison point:

| biotype | RSEM range | RSEM 30 pg | featureCounts range | featureCounts 30 pg |
|---|---|---|---|---|
| protein_coding | 95.431–96.873 % | 96.497 % | 97.285–97.987 % | 97.942 % |
| lncRNA | 1.631–2.193 % | 1.696 % | 1.218–1.530 % | 1.259 % |
| processed_pseudogene | 0.710–1.059 % | 0.898 % | 0.537–0.929 % | 0.629 % |
| transcribed_processed_pseudogene | 0.280–0.476 % | 0.423 % | 0.104–0.226 % | 0.115 % |
| unprocessed_pseudogene | 0.227–0.501 % | 0.228 % | 0.004–0.009 % | 0.005 % |
| misc_RNA | 0.084–0.276 % | 0.172 % | 0.0003–0.0011 % | 0.0003 % |
| TEC | 0.035–0.061 % | 0.049 % | 0.023–0.053 % | 0.037 % |
| miRNA | 0.0089–0.0165 % | 0.0161 % | 0.0029–0.0132 % | 0.0050 % |
| snoRNA | 0.0016–0.0035 % | 0.0032 % | 0.0023–0.0042 % | 0.0025 % |
| snRNA | 0.0001–0.0004 % | 0.0002 % | 0.0003–0.0016 % | 0.0009 % |
| ribozyme | 0.0002–0.0016 % | 0.0011 % | 0.0000–0.0001 % | 0.0000 % |
| Mt_tRNA | 0.0012–0.0031 % | 0.0019 % | 0.0030–0.0059 % | 0.0035 % |
| **rRNA** (18S relic, % of ALL) | 0.469–0.905 % | 0.594 % | 0.700–1.402 % | 0.898 % |
| **Mt_rRNA** (% of ALL) | 0.095–0.406 % | 0.121 % | 0.118–0.533 % | 0.151 % |

`well`, `input_amount`, `input_pg`, `replicate` and `qc_verdict` are carried on
every row of `nfcore_composition.tsv`. A8 (`exclude`, 18.3 % human CALB1, well
H:1) is the max of both rRNA columns and of the multiread share; it is reported,
not dropped.

## Step 2 — pipeline vs protocol

The protocol term is **not yet computable**. It needs the VASA-route
quantification of the same FLASH-seq reads (`native` arm), which the other track
is still producing; `pipeline_vs_protocol.tsv` carries
`vasa_route_flashseq_native` and `protocol_term_pp` as explicit empty columns
rather than a guess. What *is* measured are three pipeline terms:

**A — quantifier, identical reads.** RSEM vs featureCounts on the same BAMs.
`A_quantifier_pp`, `A_quantifier_log2`.

**B — annotation model, identical reads.** The VASA route counts intronic reads;
featureCounts `-t exon` and RSEM's transcriptome cannot. **28.76 %** of VASA's
uniagg reads are unspliced (15,249,326 / 53,022,536), concentrated in lncRNA
(82.2 % of its reads unspliced) and protein_coding (27.1 %). `B_intron_model_pp`
is the shift from restricting VASA to spliced-only, which is the fair analogue of
what nf-core can see.

**B2 — the analysis-set UMI-ceiling filter.** This is the largest single term and
it is easy to miss. The eight entries the analysis filter drops
(`Rmrp`, `Rn7sk`, `Rn7s1`/`Rn7s2`, `Rnu1a1`, `Rnu1b6`, `Rnu2.10`, the `Snord3b`
cluster, and a `Cmss1`-`ENSMUSG00000127915` pair) map to 1,372 uniagg rows
carrying **12,207,080 reads = 23.02 %** of the uniagg total. They were selected
for saturating a 4,096-UFI ceiling, i.e. precisely for being the most abundant
small RNAs. Dropping them is correct for molecule counting and wrong for a
read-based composition, and it moves protein_coding by **+19.19 pp**:

| biotype | V1 uniagg total | V1 minus ceiling | V3 analysis (pure) | B2 shift |
|---|---|---|---|---|
| ProteinCoding | 64.13 % | 83.31 % | 85.84 % | +19.19 pp |
| lncRNA | 13.92 % | 9.22 % | 7.24 % | −4.69 pp |
| misc_RNA | 9.82 % | 0.96 % | 1.00 % | −8.86 pp |
| snRNA | 7.03 % | 2.96 % | 3.09 % | −4.07 pp |
| snoRNA | 2.60 % | 1.53 % | 1.58 % | −1.06 pp |
| ribozyme | 1.08 % | 0.16 % | 0.16 % | −0.92 pp |

`pipeline_vs_protocol.tsv` therefore carries **all three** VASA columns (`V1`
uniagg total, `V2` uniagg spliced-only, `V3` analysis-set pure) with the ceiling
variants, so the comparator choice is explicit. **Recommendation: use `V2`
(uniagg spliced-only, no ceiling filter) against nf-core.** It matches nf-core's
exon-only annotation model and does not remove the abundant small RNAs that are
the whole point of a total-RNA protocol.

### What does not decompose cleanly

- **Counting currency.** RSEM = fractional EM `expected_count`; featureCounts =
  assigned fragments; VASA route = reads. Every figure above is a *within-route*
  fraction, where the unit cancels. A cross-route absolute count would not be
  meaningful and none is quoted.
- **VASA's combination entries have no nf-core counterpart.** 170,001 / 222,412
  analysis rows (76.44 %) are `'-'`-joined multi-gene entries carrying 7.08 % of
  reads, up to 34 genes deep (read-weighted mean 1.217 genes/entry). featureCounts
  discards such a fragment; RSEM splits it fractionally; the VASA route assigns it
  to a hierarchy-resolved combination label. Three different objects. Both
  allocation rules (`pure`-only, and equal `split`) are in the tables; neither is
  the truth.
- **tRNA.** VASA's BED carries 1,137 GtRNAdb rows that Ensembl has no counterpart
  for, so nf-core has no tRNA class at all. `Mt_tRNA` is the only overlap.
- **TEC** appears on the nf-core side only — VASA's BED build does not emit it.
- Trap (h) from `NOUMI_PATH.md` was checked on the real tables:
  `max(spliced + unspliced − total) = 0` over all 222,420 uniagg entries and 12
  cells. The vasa branch's exact-membership test holds; the double-counting risk
  is specific to the no-UMI branch.

### Where the two quantifiers disagree — and why the prediction is only half right

The prediction was that disagreement concentrates in the multi-locus small-RNA
classes. **In relative terms yes; in absolute terms no; and the mechanism is two
competing effects, not one.**

Largest **relative** disagreements (30 pg rung) are small-RNA and pseudogene
classes: misc_RNA 682×, unprocessed_pseudogene 48×, ribozyme 29×,
transcribed_unprocessed_pseudogene 19×, snRNA 5.9×,
transcribed_processed_pseudogene 3.7×, miRNA 3.0×. Largest **absolute**
disagreements are the abundant classes: protein_coding −1.448 pp, lncRNA
+0.435 pp, rRNA −0.314 pp, transcribed_processed_pseudogene +0.310 pp. Median
|log2 ratio| is 0.936 for the multi-locus small-RNA classes (n = 9) vs 0.469 for
the rest (n = 8) — a 2× separation — but median |Δpp| runs the other way, 0.0016
vs 0.2458 pp.

The class-level correlation does **not** reach significance and should not be
quoted as support: |log2(RSEM/featureCounts)| vs VASA read-weighted
`n_genes_in_entry` gives Spearman **rho = 0.327, permutation p = 0.203**
(n = 17 biotypes, seed 0, 20,000 permutations); against GTF median paralog family
size, rho = 0.388, p = 0.128. With 17 classes there is not enough power, and the
honest statement is that the *gene-level mechanism* is demonstrable while the
*class-level trend* is not established.

At gene level the two mechanisms are clear, and they push in opposite directions
with the crossover at the ~195 nt fragment length:

- **Above the fragment length, EM rescue dominates and RSEM reads higher.**
  misc_RNA is carried by `Rn7s1`/`Rn7s2` (7SL SRP RNA, 300 bp, effective_length
  123.7). These are near-identical paralogs, so featureCounts calls every read a
  multiread and discards it: **A9 featureCounts misc_RNA = 109 reads** against
  **RSEM 47,442 expected_count — a 435× ratio in one library**, driven by two
  genes. On the VASA side the same pair is the top misc_RNA carrier (70.8 % of the
  class) and appears as a `'-'`-joined two-gene entry, i.e. all three quantifiers
  see the same multi-locus problem and resolve it three different ways.
- **Below the fragment length, the effective-length floor dominates and RSEM
  reads lower.** snRNA (median 107 bp, effective_length 9.67, 72.2 % of genes
  below 10) gets 0.0002 % from RSEM against 0.0010 % from featureCounts — RSEM
  **5.9× lower**. Same for Mt_tRNA (1.9× lower, 100 % of genes below eff_length
  10) and for the 5S-class rRNA genes, which is Correction 1.

So a class can be multi-locus *and* short, and then the two effects fight; which
wins is set by length relative to the fragment, not by locus multiplicity alone.
For the FLASH-seq vs VASA comparison the practical consequence is that **neither
nf-core quantifier can be used for the short non-poly-A classes** — RSEM floors
them, featureCounts discards their multireads — and the VASA route is the only one
of the three that reports them at all.

## Files

| file | contents |
|---|---|
| `nfcore_composition.tsv` | 740 rows: library × biotype × route, counts and both denominators, with `well`/`input_amount`/`qc_verdict` |
| `nfcore_biotype_summary.tsv` | 74 rows: per biotype per route, pooled/min/max/median plus the 30 pg rung |
| `pipeline_vs_protocol.tsv` | 37 biotypes × 21 cols: all four nf-core columns, five VASA columns, terms A/B/B2, `decomposes_cleanly` |
| `multiread_disagreement.tsv` | 37 biotypes ranked by quantifier disagreement, with `n_genes_in_entry` and paralog family size |
| `rsem_effective_length.tsv` | per biotype: length, effective_length, fraction below the floor |
| `rrna_biotype_genes.tsv` | all 356 rRNA/Mt_rRNA loci with coordinates and lengths |
| `rsem_rrna_signal_by_gene.tsv` | per library: which gene carries the rRNA signal, and the `n-R5s*` share |
| `rrna_locus_attribution.tsv` | region counts per locus in **both** reads and fragments, with the reconciliation ratio |
| `vasa_intron_axis.tsv` | VASA uniagg spliced/unspliced split per biotype |
| `vasa_class_carriers.tsv` | top VASA carriers per class, flagged for the UMI-ceiling filter |
| `vasa_spliced_percell_spread.tsv` | per-cell min/median/max of the spliced-only composition |
| `nfcore_denominator_check.tsv` | the denominator verification and full featureCounts read fates |

## Pending

`vasa_route_flashseq_native` — the VASA-route quantification of the FLASH-seq
`native` arm. Until it lands, the protocol term is unmeasured and only the
pipeline terms above are quotable. When it arrives, the protocol term is
`V2_vasa_uniagg_spliced − vasa_route_flashseq_native` per biotype, and the
`nfcore_rsem_30pg − vasa_route_flashseq_native` difference gives the quantifier
term on identical reads for a third time as a consistency check.
