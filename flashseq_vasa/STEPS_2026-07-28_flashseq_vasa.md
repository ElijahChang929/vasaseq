# Step record — FLASH-seq ↔ VASA comparison (2026-07-28)

Every step, in order, with **what it was**, **what it produced**, and **whether it
made sense** — including the ones that did not. Job ids and commits are given so
each number can be traced to the run that produced it.

Conventions followed throughout: fork rather than patch upstream `a_Mapping/`;
all new code in `code/flashseq_vasa/`; report-and-flag rather than silently
filter; never state a number that was not computed.

---

## Phase 0 — establish comparability before spending compute

### Step 0.1 — Read what already existed *before* proposing anything

**Did:** read `code/flashseq/README.md` (544 lines), `sample_metadata.tsv`,
`config.sh`, the committed result tables, and the git log for `flashseq/`.

**Found:** the FLASH-seq QC was already complete end-to-end, and — critically —
`05_rrna_bwa.sh` already measured FLASH-seq rRNA **by invoking VASA's own
`ribo-bwamem.sh` + `riboread-selection.py`** against the same reference,
deliberately, so the two percentages would be one measurement rather than two
methods. `config.sh` carries a warning that reimplementing it would silently
create a second method.

**Sensible?** Yes, and it was the highest-value step in the session. The user's
request said "do from the beginning again"; taken literally that would have
rebuilt an existing like-for-like measurement as a second, non-comparable one.
Reading first turned a 3-phase plan into a 3-phase plan with one leg already
done.

**Reflection:** the lesson generalises — in a repo with this much prior work, the
first action is always "what does the git log and the README already know", not
"what can I compute".

### Step 0.2 — Ask before planning, not after

**Did:** three `ask_user` calls: quantification route (method-matched vs nf-core
vs both), rRNA denominator, comparison unit.

**Answers:** method-matched through VASA's own step 5–7; dual denominator (rRNA
as % of all reads, everything else as % of non-rRNA); per-library vs per-cell
plus a depth-matched subsample.

**Sensible?** Yes. Each of the three changes what the numbers *mean*, not just
how they are computed. Guessing the quantification route in particular would have
produced a defensible-looking table that confounds protocol with pipeline —
exactly the failure the user asked to guard against by cross-checking nf-core.

### Step 0.3 — Provenance audit (delegated)

**Did:** checksummed the rRNA fasta, `riboread-selection.py` and
`ribo-bwamem.sh` across both runs; compared STAR indices and GTF releases;
measured the read-vs-molecule asymmetry.

**Produced:** `provenance.tsv`, `read_vs_molecule.tsv`,
`read_vs_molecule_biotype_full.tsv`, `vasa_strand_check.tsv`.

**Found two asymmetries, one of which was invisible until measured:**

1. **Strand flags were mismatched.** FLASH-seq was measured `stranded=n`
   (correct — its forward strand carries 49.1–50.5% of ribosomal reads);
   VASA's 21.39% is `stranded=y`. VASA's forward share had **never been
   measured**, and turned out to be **76.1%**, not ~50%. So the same flag discards
   23.9% of VASA's ribosomal reads and 50.6% of FLASH-seq's. The published 4.3×
   was a mixed-flag ratio.
2. **Counting currency.** VASA read:UFI = 2.13–2.96 per cell (median 2.57), and
   it is expression-dependent (Spearman +0.498; 1.69 → 2.91 across expression
   deciles). The bias on biotype fractions runs in **both** directions by class,
   so no single correction factor exists.

**Sensible?** Yes, and this is the step that most changed the result. Without it
the headline number would have been a ratio of two differently-defined
quantities. The self-check inside it (the predicate reproducing step 3's own log
for all 16 barcodes) is what makes the 28.10% trustworthy.

### Step 0.4 — Validate the `smartseq_noUMI` path (delegated)

**Did:** read the branch in `countTables_2pickle_cellsSpliced.py` and
`countTables_fromPickle.py`; wrote `noumi_precheck.py` in the style of
`step7_precheck.py`; ran a stride-64 dry run of A9 end-to-end.

**Produced:** `NOUMI_PATH.md`, `noumi_precheck.py`, `dryrun_a9_noumi.sh`,
`noumi_dryrun_report.py`. Dry run passed for single-end: 20 tables, 21-line
mapStats reconciling, protein-coding 86.47% of non-rRNA, top genes
*Ftl1*/*Hsp90ab1*/*Eef1a1*/*Hspa8*/*Dppa5a*.

**Found four things that would each have wasted the ten-library run:**

1. **Paired-end fails structurally.** `bedtools bamtobed` appends `/1`,`/2` to
   the read name, and `deal_with_*.sh` appends `;CG:…;nM:<n>` *before* bamtobed
   runs, so the mate suffix lands on the end of the `nM` value and step 6's
   `int()` raises. Measured 1,033,528/1,033,528 PE rows unparseable vs
   0/534,978 SE. Upstream contract → `a_Mapping/` untouched, map R1 only.
2. **`UFICounts` degenerates to a 0/1 detection mask** (one literal UMI key
   `'A'`; verified distinct values `[1]`), and `bc2trans` is the identity on a
   mask, so `TranscriptCounts == UFICounts` elementwise. Nine TSVs of it are
   still written and **carry no abundance information**. Use `ReadCounts`.
3. **The UMI ceiling does not apply** — `K` is inferred from the UMI length, so
   `K=4` not 4096. FLASH-seq's filter set therefore differs from VASA's.
4. **The no-UMI branch tests substring containment** where the vasa branch tests
   exact membership, so an `exon-intron` label can be counted in both and
   spliced+unspliced can exceed total. Excess was 0 in the dry run — not
   guaranteed, so it is checked per library.

**Sensible?** Emphatically. Step 7 had already destroyed one run in this repo for
a reason invisible beforehand; the precheck convention exists because of it. Two
failures happened during the dry run and were **correctly distinguished**: one
was the agent's own calling error (`$P2SAMTOOLS` is a directory, not a binary),
the other was upstream behaviour. Only the first was fixed.

### Step 0.5 — The read-length confound (found, not sought)

**Found:** the dry run's tRNA table came back **empty — not because the reads are
absent but because they are too long**. `jS:IN` requires containment, and 98.5%
of tRNA features are shorter than one 151 nt read (rRNA 99.7%, miRNA 99.2%,
snoRNA 96.5%, snRNA 84.6%). Six FLASH-seq reads overlapped a tRNA; all six were
discarded. VASA cell 005 has 229 tRNA rows.

**Why it matters:** read length **alone** produces "FLASH-seq detects fewer short
species than VASA". Without a control that claim is unfalsifiable — and it is
half of what the user asked for.

**Decision (user):** hard-trim FLASH-seq R1 to VASA's aligned-length
distribution, mirroring what `05_rrna_bwa.sh` did for rRNA — one measurement,
both sides.

**Sensible?** This is the most important finding of the session. A pipeline that
silently cannot see a biotype will report its absence as biology.

---

## Interlude — the user's challenge, and what it exposed

### Step I.1 — "是你的数据还是原文?" (provenance of 21.39%)

**Did:** traced the number to `ZHA9292A1`, run PM26037, flowcell
`20260720_LH00442_0273_B23TM55LT4`, job `50788552`; read `STRANDED=y` in
`own_version/config.sh:104`; confirmed `run_mapping_stepwise.sh:41` and
`a_Mapping/README.md` both specify `y` for VASA; confirmed
`riboread-selection.py` unmodified from as-received (`git diff` vs `8bee6f6` =
mode change only). Fetched the paper and searched it.

**Found:** the number is the user's own data, counted with the same flag the
original method specifies. **The paper reports no rRNA percentage at all** — zero
sentences containing "rRNA"/"ribosomal" also contain a `%`.

**My error:** I had earlier framed a choice as if a published 20%+ figure existed
to preserve. It does not. Corrected in the same turn.

**Sensible?** The check was necessary and I should have run it before implying a
published comparator existed. Fetching the paper rather than answering from
memory is the only defensible way to make a claim about what a paper says.

### Step I.2 — "确定原文是 20%+?" → measure the published plate

**Did:** the challenge prompted a search for the published plate's *data*, not
its prose. Found SRA **`SRR14783059`** already in `data/ref/fastq_vasaplate/`
with an rRNA validation run from a prior session. Counted both predicates over
its 8 cells' own `.nsorted.all-ribo.bam`, then per-cell on all 16 of the user's
barcodes (job `6689a020`).

**Found — two things, pointing opposite ways:**

| | published (n=8) | own (n=12 real) |
|---|---|---|
| per-cell rRNA % (`y`) | 1.33–26.42% | 17.99–25.09% |
| median | 9.90% | 21.14% (**2.14×**) |
| pooled | 9.19% | 21.39% (**2.33×**) |
| spread max/min | **19.79×** | **1.39×** |

and **forward-strand share**: published 92.7–97.4% (median 95.4%); own real cells
57.7–83.8% (median 74.8%); **own blanks 86.4–93.6%**.

**The blanks are the informative part.** They took the identical path and came
out strand-specific, which rules out reference, aligner, predicate and
plate-wide handling, and points at a mechanism acting on abundant template.

**Sensible?** Yes — and it is the clearest case in the session of a user
challenge being worth more than the answer I would have given. I was
technically right that no published percentage exists, and would have stopped
there. The data was on disk the whole time.

### Step I.3 — Verify the interactive numbers with a committed script

**Did:** wrote `verify_rrna_crossplate.py` and ran it (job `0f559630`, exit 0).
15 assertions against every reported value, plus the predicate self-check.

**Result: PASS.** Predicate reproduces step 3's log for all 16 barcodes, worst
deviation **0.0048 pp**. All 15 reported values reproduce. Ratios: median
**2.136×**, pooled **2.329×**, spread **14.2×** tighter.

**Sensible?** Necessary. The cross-plate numbers were first measured in a shell
heredoc from a kernel whose state is gone — not a reproducible provenance for a
number now in a report. This step converts "I ran something and got 2.1×" into
"here is a script that gets 2.1× or fails loudly". It also surfaced that
**median and pooled ratios differ (2.14× vs 2.33×)** because the published cells
span 64k–1.26M reads, so pooling weights the deep ones — they are not
interchangeable and the report must say which is meant.

---

## Phase 1 — quantification

### Step 1.1 — nf-core cross-check (delegated, complete)

**Produced:** commit `e2906ee`, `NFCORE_CROSSCHECK.md` + 12 tables.

**Found five things that change the comparison's design:**

1. **nf-core's rRNA biotype is one 18S relic locus, not 5S.** `Rn18s-rs5`
   (`ENSMUSG00000119584`) carries 99.93–99.999% of the class in all 10
   libraries; `n-R5s*` genes carry 0.0000–0.0061%. **This corrects the existing
   `flashseq/README.md`**, which called it a 5S measurement. Mechanism: RSEM's
   effective length — implied mean fragment 195.1 nt, and 99.7% of rRNA features
   are shorter, so the one 1,849 bp locus takes the whole class.
2. **The on-disk filtered GTF is truncated** — exactly 160.0625 MiB, zero `gene`
   lines, contig 1 only, ends mid-record. The *runtime* GTF was complete (md5
   matches VASA's source GTF). Do not use the on-disk copy.
3. **The UMI-ceiling filter is the largest confound for a read-based
   comparison.** The 8 dropped entries carry **23.02% of uniagg reads** and
   shift ProteinCoding by **+19.19 pp**. They were dropped for saturating the
   UMI ceiling — i.e. *selected for being the most abundant small RNAs*, which is
   what a total-RNA protocol is for. Correct for molecule counting, wrong here.
4. **28.76% of VASA's uniagg reads are unspliced** and structurally invisible to
   both nf-core quantifiers (`featureCounts -t exon`, RSEM transcriptome).
   lncRNA is 82.2% unspliced, so comparing VASA-total against nf-core charges an
   annotation-model difference to the protocol.
5. **The two nf-core quantifiers fail in opposite directions**, crossover at the
   195.1 nt fragment length: EM rescue inflates RSEM for long multi-locus
   classes (misc_RNA 682×, from the *Rn7s1*/*Rn7s2* 7SL pair — 109 featureCounts
   reads vs 47,442 RSEM expected_count in A9); the effective-length floor
   deflates it for short ones (snRNA 5.9×). **Neither is usable for the short
   non-poly-A classes** — the VASA route is the only one of the three that
   reports them at all.

**Two things this track did that are worth naming:**

- It **refused to invent the protocol term.** The VASA-route quantification of
  FLASH-seq reads did not exist yet, so `pipeline_vs_protocol.tsv` carries those
  columns as **declared-empty (verified all-NaN)** rather than estimated.
- It **reported its own hypothesis as underpowered.** The predicted multi-locus
  class-level trend gave Spearman ρ=0.327, permutation p=0.203 (n=17, seed 0,
  20k permutations) — reported as not supported, with the gene-level mechanism
  (which *is* demonstrable) given separately.

**Sensible?** Yes. Finding 3 in particular means a filter agreed earlier for
VASA's own biology must **not** be applied to this comparison — a case where
carrying a prior decision forward unexamined would have produced a 19 pp error.

### Step 1.2 — Map and assign (delegated, running at time of writing)

Ten libraries × two arms (`native`, `vasalen`). Status at time of writing: two
650G step-6 pickle builds running (`50951230`, `50951233`), an earlier pair
completed at 1h08m each. Step 6 is superlinear in BED size.

**Not yet complete — nothing from it is reported.**

---

## Open items carried forward

1. **The protocol term is pending** on the mapping track's native arm. Until it
   lands, no FLASH-seq↔VASA composition difference should be attributed to
   chemistry: the three measured *pipeline* terms are large enough to swamp it.
2. **Use `V2_vasa_uniagg_spliced`** for the nf-core comparison (uniagg,
   exon-only, no ceiling filter) — matches nf-core's annotation model and keeps
   the abundant small RNAs. All five VASA column variants are in the tables so
   the choice is explicit and reversible.
3. **tRNA cannot be compared through nf-core at all** — VASA's BED carries 1,137
   GtRNAdb rows with no Ensembl counterpart. `Mt_tRNA` is the only overlap. TEC
   is the mirror case (nf-core only).
4. **Re-check `own_version/` count tables for antisense rRNA leakage.** With 24%
   of VASA's ribosomal reads antisense, `stranded=y` routes them to the
   *non*-ribosomal arm where they are gene-assigned. `vasaplate_check/README.md`
   documents this exact failure on the published mixed reference (one antisense
   locus inflated ~600×). The mouse-only reference was assessed as unaffected —
   but that assessment predates knowing the antisense fraction is 24%.
5. **The strand-specificity finding deserves its own follow-up.** It is a
   property of the library chemistry, not the depletion, and the blank-vs-real
   split localises it to a template-abundance-dependent step.
6. **`flashseq/README.md` needs the 5S → 18S-relic correction.**
7. **The `check` stage of `pipeline_fs.sh` fails for all 10 libraries in both
   arms, and the cause is UNKNOWN.** It reports *"no unique
   `ZHA8833A9_S*_R1_001.fastq.gz`"* for a directory that contains exactly one
   matching file, `ZHA8833A9_S108_L007_R1_001.fastq.gz`, which that glob does
   match uniquely.

   **My first explanation — that the glob misses the lane field — was wrong**, and
   wrong in an instructive way: I "confirmed" it with an `ls` that used a
   *different, broader* pattern (`ZHA8833A9*R1*.fastq.gz`) and never ran the
   checker's own glob. The listing I cited as evidence actually contradicted the
   claim. Then I compounded it by calling the failure cosmetic, which is the one
   conclusion a validator failing 20/20 does not license.

   Untested candidates: `$FSV_FASTQ` resolving elsewhere at run time; the glob not
   expanding in its actual context (quoting, `nullglob`, subshell `cwd`); or the
   uniqueness test miscounting. Diagnose by running the checker's exact command in
   its exact context, not a paraphrase of it.

   Separately established: the run itself resolved the right inputs (all ten
   libraries produced reconciling trimmed FASTQs, BAMs, BEDs and pickles). So no
   result here rests on the broken check — but the check cannot be trusted to
   catch a genuine missing input until it is fixed.

   **Method note for future sessions: to test a glob, run that glob.** A broader
   pattern matching is not evidence the narrower one matches.
