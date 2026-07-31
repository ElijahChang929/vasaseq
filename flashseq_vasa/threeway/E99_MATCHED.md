# The own plate under the published plate's own reference (Ensembl 99 / GRCm38)

**Question.** The published VASA plate was quantified on Ensembl 99 / GRCm38
(human+mouse); the own plate on Ensembl 116 / GRCm39 (mouse only). So every
published-vs-own difference mixed protocol biology with annotation release. Does
the earlier three-way conclusion survive when both plates are put under ONE
annotation?

**Answer: partly, and the correction we applied earlier was badly wrong.**
The shared-gene-universe control said annotation release explained **84.3%** of
the structural-RNA gap, leaving 2.89 pp. Under a truly matched annotation the
gap is **15.54 pp** (bracket 15.54–18.42, see *Caveats*), so release explains
only **15.3%**. The partial control did not under-correct slightly: it reported
18.6% of the real gap and missed **81.4%** of it (12.65 of 15.54 pp) — the true
value is **5.4x** what it gave.

Everything below is re-derived and asserted by `verify_e99_matched.py`
(30 checks, 0 failures). No number here was transcribed.

---

## What was run

Route chosen: **re-map**, not liftover. Liftover was rejected because it is lossy
in exactly the class under study (short, repeat-adjacent snRNA/snoRNA/misc-RNA
loci), and it would have introduced an artefact of its own into the measurement.

Cost was much lower than feared, for two reasons found on disk:
- `combined_genome.fa` (5.93 GB) and `combined.gtf` (2.08 GB) — **the exact
  inputs of the published run's own index** — still existed, so the E99 genome
  and gene models are byte-identical inputs, not merely the same release.
- Stages 1–3 were reusable, so only stages 4–7 were re-run on 16 cells.

| step | detail | time |
|---|---|---|
| index build | `star_index_151`, sjdbOverhang **150** | 35 min, 54 GB |
| stage 4 STAR | 16 cells vs E99 mixed index | 12 min |
| stage 5 assign | E99 mixed BED, **upstream** `a_Mapping/` scripts | 12 min |
| stage 6 pickle | 8 workers (`NCORES` in the driver) | 3 h 46 min, 65.8 GB peak |
| stage 7 tables | 23 tables, `mapStats.log` = 21 lines (complete) | 21 min |

Job `51075957` was allocated **16 CPUs / 220 GB** (`sacct`), with observed peak
65.8 GB. Note that stage 6 runs **8** concurrent workers, not 16 — the extra
cores served stages 4–5 only. Prior sizing work on this pipeline found 8 workers
already optimal for stage 6 (going to 16 raises worst-case peak memory without
reducing wall time), so a future re-run can safely request `-c 8 --mem=150G`.

Outputs went to a **new tree** `data/PM26037/out_E99/`; the validated E116 tables
were never written to (still dated 2026-07-28). `a_Mapping/` is untouched
(Rule 1).

> **A trap for whoever reads the E116 run next.** The existing own-plate files are
> named `..._E99_Aligned.out.bam`, but their STAR `Log.txt` records
> `genomeDir=GRCm39/star_index_151_r116`. **That `E99` is a misnomer — those files
> are Ensembl 116.** The genuinely-E99 products carry the tag `_e99mixed_` and
> live in `out_E99/`. The old name was left alone rather than renamed under a
> result already in use.

### Why sjdbOverhang differs between the two arms, and why that is correct

The published index is `sjdbOverhang 73`, right for its 75 nt reads. The own
plate's reads are **130 nt** after the 21 nt prefix skip, needing ≥129. Overhang
is baked into STAR's `Genome`/`SA` files and cannot be overridden at map time.
Reusing the published index would not error — it would **silently drop
junction-spanning alignments**, i.e. manufacture the very artefact this control
exists to remove. So it is sized per library (150 vs 73). Everything else is
held constant: same genome FASTA, same GTF, same BED, same STAR 2.7.7a, same
STAR parameters, `stranded=y` both sides, same step-5 code, same counting code.

---

## 1. The key number

Structural RNA = MiscRna + snRNA + snoRNA + scaRNA + ribozyme, as % of the
**non-rRNA read** denominator (Rule 5), `ReadCounts` on both sides (Rule 4),
mouse entries only, published-plate mouse-called cells only (n=173, Fig. 1d
rule), own plate n=12 real cells.

| arm | structural % |
|---|---|
| published plate, E99 | **2.25** |
| own plate, E116 (as before) | **20.59** |
| own plate, **E99 (matched)** | **17.78** |

| gap estimate | value | release explains |
|---|---|---|
| raw (own E116 − pub E99) | 18.34 pp | — |
| shared-gene-universe control | 2.89 pp | 84.3% *(this was wrong)* |
| **fully matched annotation** | **15.54 pp** | **15.3%** |

## 2. Why the partial control was biased, not merely approximate

Measured, not argued: of the own plate's 8,390,129 structural reads under E99,
the shared-universe rule retains only **13.10%**. It discards 3,061,250 reads on
multi-gene combination rows and 4,229,690 on gene ids absent from one release.

The own plate's structural signal sits **precisely on the rows that rule throws
away**. So restricting to the shared simple-row universe did not remove the
release term — it removed the measurement. The 2.89 pp figure is an artefact of
the control, and any conclusion resting on it needs revising.

## 3. What the annotation cannot explain, per class

| class | raw gap (pp) | matched gap (pp) | release explains |
|---|---|---|---|
| MiscRna | +8.57 | **+5.34** | 37.7% |
| snRNA | +6.87 | **+6.97** | −1.5% |
| snoRNA | +2.18 | **+2.38** | −9.2% |
| scaRNA | +0.04 | **+0.05** | −10.4% |
| ribozyme | +0.68 | **+0.80** | −16.8% |
| ProteinCoding | −28.10 | **−15.25** | 45.7% |
| lncRNA | +11.47 | **+0.88** | 92.4% |

Positive = own plate above published. Reading:

- **snRNA, snoRNA, scaRNA, ribozyme are real and not annotation.** Their gaps are
  *unchanged or slightly larger* under matched annotation (negative "explains").
  snRNA at +6.97 pp — a 44-fold excess over the published plate's 0.16% — is the
  single largest surviving effect.
- **MiscRna is genuinely mixed**: 37.7% release, 5.34 pp real.
- **lncRNA was almost entirely annotation** (92.4%). Ensembl 116 annotates far
  more lncRNA than 99; that leg of the earlier comparison should be dropped.
- **ProteinCoding**: 45.7% release. The own plate still sits 15.25 pp *below* the
  published plate, which is the mirror of the structural excess — reads are being
  re-attributed, not lost (see §5).

## 4. Gene detection is release-free and near-equal

Depth-matched at 10,000 in-scope reads (the deepest rung every unit of both E99
arms reaches natively), single-gene protein-coding mouse entries, the identical
deterministic estimator `E[genes] = Σ_g 1−(1−p)^{c_g}` used in
`mk_detection_threeway.py`:

| arm | n | median genes |
|---|---|---|
| published, E99 | 173 | 4,069 |
| own, E99 | 12 | 4,245 |
| own, E116 | 12 | 4,592 |

Own/published = **1.04x** under one annotation. **n=12 vs 173: descriptive, not a
test.** Note the own plate's own E116→E99 change (4,592→4,245, 1.08x) is larger
than the between-plate difference — most of the apparent detection advantage in
the earlier comparison was the newer annotation.

## 5. Assignment churn (diagnostic)

Own plate, simple rows, E99 vs E116:

- assigned under E99: 42,897,518 reads on 33,289 gene ids
- assigned under E116: 41,955,553 reads on 52,386 gene ids
- on ids **only in E99**: 4,526,008 (10.55% of E99)
- on ids **only in E116**: 5,848,536 (13.94% of E116)
- shared ids: 31,501; net change **+941,965 reads (+2.2%)**

Total assigned reads barely move (E99/E116 all-row ratio 0.9487) while 10–14%
change gene identity: the annotation **re-attributes** reads rather than gaining
or losing them.

## Caveats — what is still not controlled

1. **Species filter (quantified, bracketed).** The own plate is mouse-only but
   the E99 reference is human+mouse. The mouse-only filter removes 6.18% of its
   rows (1.91% human ids, 4.28% mixed). For the key number this cannot bias
   anything — *the identical filter is applied to both arms on the identical
   reference* (it removes 3.89% of the published arm). For the secondary
   own-E99-vs-own-E116 delta it is a genuine asymmetry, since a mouse-only E116
   reference has no human rows to drop. Recomputing with the filter **off** gives
   structural 20.67% and a matched gap of **18.42 pp**. So:
   **matched gap bracket = [15.54, 18.42] pp.** Both ends are 5.4–6.4x the
   2.89 pp shared-universe figure, so the conclusion is insensitive to the choice
   even though the point estimate is not. The own plate's human-id share (1.91%)
   is *lower* than the published plate's mouse cells (2.62%), so the mixed
   reference costs the own plate no more than it cost the published plate.
2. **Read length: 130 nt vs 75 nt, not equalised — but measured, and small.**
   Trap 8 applies: VASA's `jS:IN` rule requires a read to be *contained* in the
   feature, and under Ensembl 99 **96.48%** of mouse snoRNA and **99.23%** of
   miRNA exon features are shorter than one 151 nt read
   (`annotation_feature_length.tsv`, `pct_exon_features_shorter_than_151`), vs
   84.55% of snRNA and 23.49% of MiscRna. A longer read is therefore structurally
   *worse* at short features, so this confound works **against** the observed
   own-plate excess rather than creating it. Its measured size is also tiny: the
   FLASH-seq native-vs-VASA-length contrast (`readlen_effect_pp` in
   `composition_threeway.tsv`) sums to **−0.0007 pp** across the five structural
   classes (largest single class |0.0029| pp, MiscRna). So read length cannot
   account for a 15.54 pp gap — it is four orders of magnitude too small.
3. **Depth: ~4.4M vs ~0.5M reads per cell.** Addressed by depth-matching in §4;
   not addressable in composition.
4. **n=12 vs n=173.** Every own-plate figure here rests on 12 cells.
5. **Step-6/7 code is shared, so not a confound** — but note the own plate's
   combination-key share is much higher under E116 (20.57%) than E99 (9.08%),
   which is itself part of what the release changes.

## Bottom line

The user's instinct was right, and acting on it changed a conclusion. Of the
earlier three-way findings:

- **Survives, strengthened:** the own plate's snRNA/snoRNA/scaRNA/ribozyme excess
  over the published plate is real, protocol-driven, and *not* annotation. The
  headline structural gap is 15.54 pp (bracket to 18.42), not 2.89 pp.
- **Survives, reduced:** MiscRna — 37.7% of it was release.
- **Does not survive:** the lncRNA gap (92.4% annotation) and most of the gene
  detection advantage.
- **Needs revising:** any statement citing "84.3% of the structural gap is
  annotation release". The shared-gene-universe control is biased for this
  quantity because the signal lives on the rows it discards, and it should not be
  used for composition claims about the short-RNA classes.

## Files

Scripts (`code/flashseq_vasa/threeway/`): `build_e99_index_oh150.sh`,
`run_own_e99.sh`, `e99_own_precheck.py`, `e99_matched.py`,
`verify_e99_matched.py`.

Tables (`res/threeway/`): `own_plate_E99.tsv`, `e99_matched_composition.tsv`,
`e99_matched_structural.tsv`, `e99_matched_assignment.tsv`,
`e99_matched_detection.tsv`, `e99_matched_provenance.tsv`,
`e99_matched_report.txt`. Figure: `e99_matched.png`.

Jobs: index `51075072`; stages 4–7 `51075957`. (A first index attempt,
`51074873`, died in my own disk-headroom precheck on a `df -P --output=` option
conflict before STAR started.)
