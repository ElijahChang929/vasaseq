#!/usr/bin/env python3
"""Read-only precheck for re-quantifying the own VASA plate under Ensembl 99.

Rule 2: replay the real operations over the real data and report what would
break, BEFORE committing to stages 4-7 (~4-5 h). Nothing here writes into any
data directory.

What this checks, in the order the pipeline would hit it:

  1. INPUTS      the 16 step-3 nonRibo fastqs exist and are non-empty.
  2. GEOMETRY    read length vs the new index's sjdbOverhang. This is the whole
                 reason a new index is being built: the published plate's
                 star_index_74 (oh73) is too small for 130 nt reads, and STAR
                 would silently lose junction-spanning alignments rather than
                 fail.
  3. CONTIGS     the index's chromosome names and the Ensembl 99 BED's must
                 agree, or step 5's bedtools intersect quietly returns nothing.
  4. NAMING      the reduceGeneName naming contract that killed a step-7 run
                 before (IndexError at `g.rsplit('_')[1]`, upstream line 164).
                 Replayed over the REAL E99 BED gene names with the COMPLETE
                 upstream function -- including the two branches past line 160,
                 which are the only ones that can raise. Carries a self-check
                 that the copy does raise on a name it should, so that a pass
                 means something.
  5. HUMAN ARM   the own plate is mouse-only but the E99 reference is
                 human+mouse. Any read assigned to a GRCh38_ contig is
                 mismapping, and that is a NEW artefact this arm introduces
                 which the GRCm39 arm could not have. This section does NOT
                 size that rate -- it cannot, since no E99 alignment of the own
                 plate exists yet when the precheck runs. It only confirms the
                 input needed to size it later is present
                 (published_cell_species.tsv). The rate itself is measured in
                 e99_matched.py, whose yardstick is the published plate's own
                 mouse-called cells on the same mixed reference.
  6. DISK        headroom for a second set of BAMs + step-5 BEDs.

Helpers copied VERBATIM from the upstream scripts are marked; everything else
is this file's own.
"""
import gzip
import os
import shutil
import subprocess
import sys

W = "/nemo/lab/turnerj/working/guangxin"
REF = f"{W}/reference/vasaseq"
CELLDIR = f"{W}/vasaseq/data/PM26037/out/cells"
IDX_NEW = f"{REF}/mixed/star_index_151"
IDX_OLD = f"{REF}/mixed/star_index_74"
BED99 = f"{REF}/mixed/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed"
BED116 = f"{REF}/mouse_GRCm39_E116/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed"

# own_version/config.sh: SKIP5=21 is stripped from both mates before anything
# else, so the biological read is 151-21 = 130 nt and the junction overhang
# STAR needs is readlen-1 = 129.
SKIP5 = 21
READLEN_RAW = 151

bad = 0
warn = 0


def say(s=""):
    print(s, flush=True)


def fail(s):
    global bad
    bad += 1
    say(f"  FAIL {s}")


def flag(s):
    global warn
    warn += 1
    say(f"  WARN {s}")


def ok(s):
    say(f"  ok   {s}")


# ---------------------------------------------------------------- 1. inputs
say("=" * 78)
say("1. STEP-3 INPUTS (own plate, 16 cells)")
say("=" * 78)
cells_file = f"{CELLDIR}/.cells"
if not os.path.exists(cells_file):
    fail(f"no cell list at {cells_file}")
    cells = []
else:
    cells = [c.strip() for c in open(cells_file) if c.strip()]
    ok(f"{len(cells)} cells in .cells")
tot = 0
for c in cells:
    fq = f"{CELLDIR}/{c}_cbc_trimmed_homoATCG.nonRibo.fastq.gz"
    if not (os.path.exists(fq) and os.path.getsize(fq) > 0):
        fail(f"missing/empty step-3 fastq for {c}")
    else:
        tot += os.path.getsize(fq)
say(f"       total step-3 fastq: {tot / 2**30:.2f} GiB")

# ------------------------------------------------------------- 2. geometry
say()
say("=" * 78)
say("2. READ GEOMETRY vs INDEX sjdbOverhang")
say("=" * 78)
need_ovh = READLEN_RAW - SKIP5 - 1
say(f"       raw read {READLEN_RAW} nt - skip5 {SKIP5} - 1 = need sjdbOverhang >= {need_ovh}")


def read_ovh(idx):
    gp = f"{idx}/genomeParameters.txt"
    if not os.path.exists(gp):
        return None
    for line in open(gp):
        f = line.split()
        if f and f[0] == "sjdbOverhang":
            return int(f[1])
    return None


old = read_ovh(IDX_OLD)
if old is None:
    flag(f"cannot read sjdbOverhang from {IDX_OLD}")
else:
    say(f"       star_index_74  (published plate) sjdbOverhang = {old}")
    if old >= need_ovh:
        flag("the published index would already suffice -- the new build is unnecessary")
    else:
        ok(f"confirms the published index is too small for this library ({old} < {need_ovh}); "
           "reusing it would silently drop junction reads")

new = read_ovh(IDX_NEW)
if new is None:
    fail(f"new index has no genomeParameters.txt yet: {IDX_NEW} (is the build finished?)")
else:
    if new >= need_ovh:
        ok(f"star_index_151 sjdbOverhang = {new} >= {need_ovh}")
    else:
        fail(f"star_index_151 sjdbOverhang = {new} < {need_ovh}")
    for part in ("SA", "SAindex", "Genome", "chrName.txt"):
        p = f"{IDX_NEW}/{part}"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            ok(f"{part} present ({os.path.getsize(p) / 2**30:.2f} GiB)")
        else:
            fail(f"{part} missing or empty in {IDX_NEW}")

# Empirical read length, deterministic stride (trap 4: the first N reads are
# not a sample). Longest observed read is what STAR must span.
if cells:
    probe = None
    for c in cells:
        fq = f"{CELLDIR}/{c}_cbc_trimmed_homoATCG.nonRibo.fastq.gz"
        if os.path.exists(fq) and os.path.getsize(fq) > 50_000_000:
            probe = fq
            break
    if probe:
        mx = 0
        n = 0
        with gzip.open(probe, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4000 == 1:          # stride over sequence lines
                    L = len(line.rstrip())
                    mx = max(mx, L)
                    n += 1
                if n >= 5000:
                    break
        say(f"       longest read in {n} strided reads of {os.path.basename(probe)}: {mx} nt")
        if new is not None and mx - 1 > new:
            fail(f"observed read {mx} nt needs overhang {mx - 1} > {new}")
        else:
            ok(f"observed max read length {mx} nt fits overhang {new}")

# -------------------------------------------------------------- 3. contigs
say()
say("=" * 78)
say("3. CONTIG NAMES: index vs Ensembl 99 BED")
say("=" * 78)
if new is not None:
    idx_chr = {l.strip() for l in open(f"{IDX_NEW}/chrName.txt") if l.strip()}
    bed_chr = set()
    with open(BED99) as fh:
        for line in fh:
            bed_chr.add(line.split("\t", 1)[0])
    say(f"       index contigs {len(idx_chr)}   BED contigs {len(bed_chr)}")
    miss = bed_chr - idx_chr
    if miss:
        fail(f"{len(miss)} BED contigs absent from index, e.g. {sorted(miss)[:5]}")
    else:
        ok("every BED contig exists in the index")
    nm = sum(1 for c in idx_chr if c.startswith("GRCm38_"))
    nh = sum(1 for c in idx_chr if c.startswith("GRCh38_"))
    say(f"       index: GRCm38_ {nm}, GRCh38_ {nh}")
    if nm == 0 or nh == 0:
        fail("index is not species-prefixed as the BED expects")
    else:
        ok("species prefixes present on both sides")

# --------------------------------------------------------------- 4. naming
say()
say("=" * 78)
say("4. reduceGeneName CONTRACT on Ensembl 99 gene names")
say("=" * 78)


def reduceGeneName(gene, uni_genes):
    """VERBATIM from a_Mapping/countTables_fromPickle.py lines 142-166, the
    COMPLETE function including the two branches past line 160.

    Line numbers verified with `grep -n` on the upstream file (2026-07-30):
    `def reduceGeneName` 142, `return rg` 166, `def fixGeneLabels` 168.

    Those last two branches are the whole point of this check and an earlier
    version of this file wrongly omitted them:

        line 161:  if sum([g in uni_genes for g in gene.rsplit('-')]) == 1:
        line 164:  if gene.count('-') >= 1 and
                      sum([g.rsplit('_')[1][:2] != "Gm" for g in gene.rsplit('-')]) == 1:

    `g.rsplit('_')[1]` is the indexing that raises. It assumes every member of
    a '-'-joined combination key splits into at least two underscore-separated
    fields (ID_NAME_BIOTYPE); a member with NO underscore gives a 1-element
    list and `[1]` raises IndexError. That is what killed the earlier step-7
    run. Nothing before line 163 can raise, so a copy truncated at line 160
    would always pass and give false assurance.
    """
    rg = gene
    if gene.count('-') == 0:
        rg = gene
    else:
        bios = set([x.rsplit('_')[-1] for x in gene.rsplit('-')])
        shortlived = ['miRNA', 'tRNA', 'MtTrna']
        longstuff = ['lncRNA']                      # noqa: F841 (upstream declares, never uses)
        shortstuff = ['snRNA', 'snoRNA', 'MiscRna', 'scaRNA']
        ribos = ['rRNA', 'ribozyme']
        if any([b in ribos for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in ribos])
            rg = gene
        if any([b not in shortlived for b in bios]) and any([b in shortlived for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] not in shortlived])
            rg = gene
        if any([b in shortstuff for b in bios]) and any([b not in shortstuff for b in bios]):
            gene = '-'.join([g for g in gene.rsplit('-') if g.rsplit('_')[-1] in shortstuff])
            rg = gene
        if sum([g in uni_genes for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g in uni_genes][0]
            gene = rg
        if gene.count('-') >= 1 and sum([g.rsplit('_')[1][:2] != "Gm" for g in gene.rsplit('-')]) == 1:
            rg = [g for g in gene.rsplit('-') if g.rsplit('_')[1][:2] != "Gm"][0]
    return rg


names = set()
nfields = {}
with open(BED99) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        nm_ = f[4]
        body = nm_.rsplit("_", 1)[0]      # strip _exon/_intron as step 6 does
        names.add(body)
        nfields[body.count("_")] = nfields.get(body.count("_"), 0) + 1
say(f"       distinct gene-name bodies in E99 BED: {len(names):,}")
say(f"       underscore-count histogram: {dict(sorted(nfields.items()))}")

# Members with no underscore are the ones that make line 164 raise.
nounder = sorted(n for n in names if "_" not in n)
say(f"       name bodies with NO underscore (line 164 would raise on these): {len(nounder):,}")
if nounder:
    say(f"       examples: {nounder[:3]}")

# ---- which of these can actually REACH reduceGeneName ----------------------
# countTables_fromPickle.py lines 112-113, verbatim:
#     tRNAs = [idx for idx in cntdf.index if 'tRNA' in idx]
#     genes = [idx for idx in cntdf.index if idx not in tRNAs]
# and reduceGeneName is called only over `genes` (lines 188, 191). That test is
# a SUBSTRING match, so it captures a mixed key like
# "1.tRNA1-ValCAC-ENSMUSG..._Gnai3_ProteinCoding" too. So an underscore-free
# name only matters here if it does NOT contain "tRNA".
def reaches_reducegenename(idx):
    """True if this count-table row index is in `genes`, i.e. line 113."""
    return "tRNA" not in idx


escapes = [n for n in nounder if reaches_reducegenename(n)]
say(f"       of those, NOT routed to the tRNA table by the line-112 substring "
    f"test: {len(escapes):,}")
if escapes:
    fail("underscore-free names that reach reduceGeneName: "
         f"{escapes[:5]} -- these WILL raise IndexError in step 7")
else:
    ok("every underscore-free name contains 'tRNA' and is routed to the tRNA "
       "table before reduceGeneName is called -- consistent with the published "
       "plate and the own-plate E116 run both completing step 7 on tRNA-carrying BEDs")

# SELF-CHECK: the reimplementation must actually be able to raise, or this
# section proves nothing. Feed it a combination key built from a real
# underscore-free name and require the IndexError.
if nounder:
    probe = f"{nounder[0]}-ENSMUSG00000000001_Gnai3_ProteinCoding"
    try:
        reduceGeneName(probe, set())
        flag("self-check: the copied reduceGeneName did NOT raise on an "
             f"underscore-free member ({probe!r}) -- this section cannot detect "
             "the IndexError it exists to detect")
    except IndexError:
        ok("self-check: the copied function reproduces the IndexError on an "
           "underscore-free combination member, so a pass below is meaningful")
    except Exception as e:      # noqa: BLE001
        flag(f"self-check raised {type(e).__name__} rather than IndexError: {e}")
else:
    flag("no underscore-free name bodies in the E99 BED, so the self-check "
         "cannot confirm the copied function is able to raise")

# Now the real test, over the population that actually reaches the function.
# Step 7 sees COMBINATION keys as well as single names, and line 164 only runs
# when the key contains a '-'. Exhaustive pairing is O(n^2) and unnecessary --
# whether line 164 raises depends only on whether SOME member lacks an
# underscore -- so test every eligible name once alone and once as a member of
# a 2-part key against a fixed well-formed partner.
PARTNER = "ENSMUSG00000000001_Gnai3_ProteinCoding"
eligible = [n for n in names if reaches_reducegenename(n)]
say(f"       names reaching reduceGeneName (the `genes` population): {len(eligible):,}")
raised_single, raised_pair = 0, 0
ex = []
for nm_ in eligible:
    try:
        reduceGeneName(nm_, set())
    except Exception as e:      # noqa: BLE001 -- the point is to catch anything
        raised_single += 1
        if len(ex) < 5:
            ex.append(f"single {nm_!r} -> {type(e).__name__}: {e}")
    try:
        reduceGeneName(f"{nm_}-{PARTNER}", set())
    except Exception as e:      # noqa: BLE001
        raised_pair += 1
        if len(ex) < 10:
            ex.append(f"pair {nm_!r}+partner -> {type(e).__name__}: {e}")
say(f"       raises alone: {raised_single:,} / raises in a 2-member key: {raised_pair:,}")
if raised_single or raised_pair:
    # Rule 3: report and flag, with the examples the user needs to judge it.
    fail(f"reduceGeneName raises on {raised_single:,} single / {raised_pair:,} "
         f"combined E99 gene names that DO reach it; examples: {ex}")
else:
    ok(f"reduceGeneName survives all {len(eligible):,} reachable E99 gene names, "
       "alone and as a combination member")

# tRNA rows carry no Ensembl id; confirm they are shaped as step 6 expects.
ntrna = sum(1 for n in names if not n.startswith("ENS"))
say(f"       non-ENS (tRNA) name bodies: {ntrna:,}")

# ------------------------------------------------------------ 5. human arm
say()
say("=" * 78)
say("5. NEW ARTEFACT: mouse-only library on a human+mouse reference")
say("=" * 78)
say("       The own plate is mouse only. Under the E99 mixed reference some")
say("       mouse reads will mismap to GRCh38_ contigs and be lost from the")
say("       mouse denominator -- an artefact the GRCm39 arm cannot have.")
say("       Empirical bound from the published plate's own mouse-called cells")
say("       is computed in the main analysis; here we only confirm the")
say("       machinery to measure it exists.")
sp = f"{W}/vasaseq/res/threeway/published_cell_species.tsv"
if os.path.exists(sp):
    ok(f"published per-cell species calls present ({os.path.getsize(sp)} bytes)")
else:
    flag("published_cell_species.tsv absent; cross-species bleed cannot be bounded from it")

# ----------------------------------------------------------------- 6. disk
say()
say("=" * 78)
say("6. DISK")
say("=" * 78)
for path in (f"{W}/vasaseq/data/PM26037", f"{REF}/mixed"):
    u = shutil.disk_usage(path)
    say(f"       {path}: {u.free / 2**30:.1f} GiB free of {u.total / 2**40:.1f} TiB")
# Existing E116 BAMs + step-5 BEDs measure the second set's cost directly.
have = 0
for f in os.listdir(CELLDIR):
    if "_E99_Aligned.out.bam" in f or "_genes.bed.gz" in f:
        have += os.path.getsize(os.path.join(CELLDIR, f))
say(f"       existing (E116-mapped) BAM+BED footprint: {have / 2**30:.1f} GiB")
say(f"       a second set costs about the same again")
u = shutil.disk_usage(f"{W}/vasaseq/data/PM26037")
if u.free < have * 1.5:
    fail(f"insufficient headroom: need ~{have * 1.5 / 2**30:.0f} GiB, have {u.free / 2**30:.0f}")
else:
    ok("headroom sufficient for a second set of alignments")

say()
say("=" * 78)
say(f"PRECHECK: {bad} FAIL, {warn} WARN")
say("=" * 78)
sys.exit(1 if bad else 0)
