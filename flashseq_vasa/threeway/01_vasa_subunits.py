#!/usr/bin/env python3
"""
01_vasa_subunits.py -- measure WHERE the residual rRNA sits, on all three
datasets, by ONE code path.

WHY THIS EXISTS
---------------
`res/flashseq_vasa/rrna_comparison.tsv` carries subunit composition for
FLASH-seq only. Its eight pct_* columns are EMPTY for both VASA plates -- the
published plate (8 cells) and the own plate (16 barcodes). Subunit composition
was therefore a FLASH-seq-only measurement, and any "three-way subunit
comparison" built on that file would have been reporting empty cells as zeros.

This script computes it for the two VASA sides, and RE-computes it for
FLASH-seq through the same function, so all three come from one code path
rather than from two that are merely believed to agree.

THE ONE MEASUREMENT METHOD
--------------------------
`bucket()` below is the same bucketing as:
  - code/flashseq/05_rrna_bwa_report.py :: composition()   (FLASH-seq side)
  - code/I_Gene_expression/own_version/step3_report.py :: parse_bam()  (VASA side)
i.e. an alignment against the 47S unit is assigned to a subunit by its
1-based leftmost position, using rrna_intervals.tsv; anything on an Mt_rRNA
record is mito; everything else is a dispersed rRNA entry. Two buckets are
ADDED here because the mixed reference has records the mouse-only one does not
(see THE REFERENCE ASYMMETRY): `human45S`, and mito split by species.

Input is the .Ribo.bam that riboread-selection.py wrote (VASA sides) or the
.riboloci.tsv that 05_rrna_bwa.sh wrote (FLASH-seq side). Both hold one record
per read already CALLED ribosomal under that dataset's own strand flag, so this
script never re-decides what is ribosomal -- it only asks where those reads
landed. That keeps the strand-flag decision exactly where the brief put it
(VASA y, FLASH-seq n) and out of this script.

THE REFERENCE ASYMMETRY -- MEASURED, NOT ASSUMED
------------------------------------------------
The published plate was run against `unique_rRNA_human_mouse.v3.fa` (921
entries, human+mouse); the own plate and FLASH-seq against
`unique_rRNA_mouse.v2.fa` (357 entries, mouse only). That is not only a
species difference. Matching the 357 mouse-prefixed entries of the mixed
reference to the 357 entries of the mouse-only reference BY ENSEMBL ACCESSION
gives only 16 shared accessions, and 13 of those 16 differ in sequence -- the
mixed reference draws its dispersed rRNA genes from Ensembl 99 / GRCm38 and the
mouse-only one from Ensembl 116 / GRCm39, which are largely disjoint gene sets.

The 47S record is the exception and it is the load-bearing one:
`mouse_rDNA_47S_BK000964.3_1-13403` is byte-identical (13,403 nt) in both
files. So:

  * The seven 47S-derived subunits ARE comparable across all three datasets,
    because they are positions on an identical reference record. This script
    therefore reports composition renormalised to 47S-derived reads
    (`pct47_*`) as the comparable quantity.
  * The dispersed buckets are NOT comparable across the reference boundary.
    They are reported as their own columns (`pct_other`, `pct_5S`,
    `pct_human45S`) with that stated, never folded into the comparable set.

5S is a separate case and is handled honestly rather than approximated: 5S is
Pol III-transcribed and is NOT part of the 47S unit, so it can only come from
the dispersed entries. The mouse-only reference names them (72 entries whose
symbol matches n-R5s/Rn5s), so 5S is measurable for the own plate and for
FLASH-seq -- which share that reference exactly. The mixed reference's mouse
entries carry no gene symbol, so 5S is not nameable there. This script attempts
a sequence-identity rescue (exact match, either orientation, against the
mouse-only 5S set) and REPORTS ITS COVERAGE, so the writeup can say how much of
the published plate's dispersed signal could and could not be resolved instead
of guessing.

THE 5'ETS ARTEFACT CHECK IS NOT OPTIONAL
----------------------------------------
05_rrna_bwa_report.py::composition() returns a 5'ETS peak-bin share for a
specific reason recorded there: a poly-T population that survived trimming
aligned to a T-rich stretch of the 5'ETS and inflated it to 31-62% of a blank
cell's ribosomal calls, and was caught because 88.9% of those hits sat in ONE
200 nt window against 10.5% in a real cell. A high 5'ETS share is legitimate
for a total-RNA protocol, so the share alone cannot tell the two apart -- the
concentration can. This script reproduces that check (`ets_peak_bin_pct`) for
every unit on all three datasets. Any 5'ETS claim without it is uninterpretable.

DENOMINATORS
------------
This script reports composition as a share of that unit's own ribosomal reads,
so it is denominator-free by construction. The rRNA PERCENTAGE denominators are
assembled in 02_threeway_table.py, which is where that argument belongs.

It also validates `depletion_v2_vs_v3.tsv` against the per-cell ribo-map.log
files, because that table is the source for the plate-wide published figure and
nothing had checked it against the logs it was derived from.

USAGE
    01_vasa_subunits.py <outdir>
Environment (all defaulted, so the job script stays short):
    VP_V3_DIR, OWN_CELLS_DIR, FS_RRNA_BWA_DIR, MIXED_FA, MOUSE_FA, INTERVALS
"""
import collections
import csv
import glob
import os
import re
import sys

import pysam

ROOT = "/nemo/lab/turnerj/working/guangxin/vasaseq"
REFS = "/nemo/lab/turnerj/working/guangxin/reference/vasaseq"

V3_DIR = os.environ.get("VP_V3_DIR", f"{ROOT}/data/ref/fastq_vasaplate/vasaplate_out_v3")
OWN_DIR = os.environ.get("OWN_CELLS_DIR", f"{ROOT}/data/PM26037/out/cells")
FS_DIR = os.environ.get("FS_RRNA_BWA_DIR", f"{ROOT}/res/flashseq/rrna_bwa")
MIXED_FA = os.environ.get("MIXED_FA", f"{REFS}/mixed/unique_rRNA_human_mouse.v3.fa")
MOUSE_FA = os.environ.get("MOUSE_FA", f"{REFS}/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa")
INTERVALS = os.environ.get("INTERVALS", f"{REFS}/mouse_GRCm39_E116/rrna_intervals.tsv")
DEPLETION = os.environ.get("DEPLETION", f"{ROOT}/data/ref/fastq_vasaplate/depletion_v2_vs_v3.tsv")

ETS_BIN = 200          # nt; the window the poly-T artefact was caught in
SUB_ORDER = ["5ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3ETS"]
EXTRA_ORDER = ["47S_outside", "human45S", "mito_mouse", "mito_human", "5S", "other"]


# ---------------------------------------------------------------------------
# reference bookkeeping
# ---------------------------------------------------------------------------
def load_fasta(path):
    seqs, h = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                h = line[1:]
                seqs[h] = []
            elif h is not None:
                seqs[h].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_intervals(path):
    """[(name, start, end)], 1-based inclusive. Same loader as step3_report.py."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def build_5s_sets():
    """(named 5S refs in the mouse-only file, 5S refs rescued in the mixed file).

    The rescue is by exact sequence identity in either orientation. It is
    reported, not trusted: main() prints how many dispersed entries and how many
    dispersed READS it resolved, so 'not measurable' can be said with a number.
    """
    mouse = load_fasta(MOUSE_FA)
    mixed = load_fasta(MIXED_FA)
    named = {h for h in mouse if re.search(r"_(n-R5s|Rn5s)", h, re.I)}
    seq5s = set()
    for h in named:
        seq5s.add(mouse[h])
        seq5s.add(revcomp(mouse[h]))
    rescued = {h for h, s in mixed.items() if s in seq5s}
    return named, rescued, mouse, mixed


# ---------------------------------------------------------------------------
# THE bucketing -- one function, all three datasets
# ---------------------------------------------------------------------------
def bucket(records, intervals, five_s):
    """records: iterable of (ref_name, 1-based pos). -> (counts, ets_peak_pct).

    `five_s` is the set of reference names to call 5S. Everything not on the
    47S unit, not an Mt_rRNA record, not a human 45S record and not in `five_s`
    lands in `other`, which is deliberately a residual rather than a claim.
    """
    c = collections.Counter()
    ets_bins = collections.Counter()
    for ref, pos in records:
        if "rDNA_47S" in ref:
            for name, st, en in intervals:
                if st <= pos <= en:
                    c[name] += 1
                    if name == "5ETS":
                        ets_bins[(pos // ETS_BIN) * ETS_BIN] += 1
                    break
            else:
                c["47S_outside"] += 1
        elif "_45S_" in ref and ref.startswith("human"):
            c["human45S"] += 1
        elif "Mt_rRNA" in ref:
            c["mito_human" if ref.startswith("human") else "mito_mouse"] += 1
        elif ref in five_s:
            c["5S"] += 1
        else:
            c["other"] += 1
    n_ets = sum(ets_bins.values())
    peak = 100.0 * max(ets_bins.values()) / n_ets if n_ets else 0.0
    return c, peak


def from_bam(path, intervals, five_s):
    bam = pysam.AlignmentFile(path, "rb", check_sq=False)

    def gen():
        for r in bam.fetch(until_eof=True):
            yield (r.reference_name or ""), r.reference_start + 1

    out = bucket(gen(), intervals, five_s)
    bam.close()
    return out


def from_riboloci(path, intervals, five_s):
    """FLASH-seq's per-read locus table: flag / ref / pos, written by 05_rrna_bwa.sh."""
    def gen():
        with open(path) as fh:
            header = next(fh, None)
            if header is None or not header.startswith("flag"):
                sys.exit(f"FATAL: {path} has no flag/ref/pos header")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3:
                    yield f[1], int(f[2])
    return bucket(gen(), intervals, five_s)


def parse_ribo_log(path):
    """-> (nreads, nunmapped, mapped_sum) from a .ribo-map.log."""
    body = open(path, errors="replace").read()
    total = int(re.search(r"Number of reads:\s*(\d+)", body).group(1))
    unmapped = int(re.search(r"Number of unmapped reads:\s*(\d+)", body).group(1))
    rg = {k: int(v) for k, v in re.findall(r"^\t(\S+):\s*(\d+)", body, re.M)}
    return total, unmapped, sum(rg.values())


# ---------------------------------------------------------------------------
def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(outdir, exist_ok=True)
    log = []

    def emit(s=""):
        print(s, flush=True)
        log.append(s)

    intervals = load_intervals(INTERVALS)
    assert [n for n, _, _ in intervals] == SUB_ORDER, \
        f"interval names {[n for n,_,_ in intervals]} != expected {SUB_ORDER}"

    emit("=" * 78)
    emit("REFERENCE ASYMMETRY (measured here, so the writeup can state it)")
    emit("=" * 78)
    named5s, rescued5s, mouse_ref, mixed_ref = build_5s_sets()
    mm = {h for h in mixed_ref if h.startswith("mouse_")}
    emit(f"  mixed  {os.path.basename(MIXED_FA)}: {len(mixed_ref)} entries "
         f"({len(mm)} mouse-prefixed, {len(mixed_ref)-len(mm)} human)")
    emit(f"  mouse  {os.path.basename(MOUSE_FA)}: {len(mouse_ref)} entries")

    def acc(h):
        h2 = h[len("mouse_"):] if h.startswith("mouse_") else h
        m = re.match(r"(ENSMUSG\d+)", h2)
        st = "(-)" if h2.endswith("(-)") else ("(+)" if h2.endswith("(+)") else "")
        return (m.group(1) + st) if m else h2 + st

    a_mixed = {acc(h): mixed_ref[h] for h in mm}
    a_mouse = {acc(h): mouse_ref[h] for h in mouse_ref}
    shared = set(a_mixed) & set(a_mouse)
    seqdiff = [k for k in shared if a_mixed[k] != a_mouse[k]]
    emit(f"  mouse accessions shared between the two references: {len(shared)} / "
         f"{len(a_mouse)}   of which sequence-differing: {len(seqdiff)}")
    k47 = "mouse_rDNA_47S_BK000964.3_1-13403"
    same47 = mixed_ref.get(k47) == mouse_ref.get(k47) and mixed_ref.get(k47) is not None
    emit(f"  47S record identical in both: {same47} (len {len(mixed_ref.get(k47,''))})")
    assert same47, "47S record differs between references -- pct47_* would NOT be comparable"
    emit(f"  5S entries named in the mouse-only reference: {len(named5s)}")
    emit(f"  5S entries rescued by exact sequence match in the mixed reference: "
         f"{len(rescued5s)}")
    emit("  -> pct47_* (renormalised to 47S-derived reads) is comparable across all")
    emit("     three; pct_other / pct_5S / pct_human45S are NOT comparable across")
    emit("     the reference boundary and are reported as their own columns.")

    rows = []

    # ---------------- published plate, 384 cells ---------------------------
    emit()
    emit("=" * 78)
    emit("PUBLISHED PLATE -- SRR14783059, vasaplate_out_v3, mixed reference")
    emit("=" * 78)
    bams = sorted(glob.glob(os.path.join(V3_DIR, "*_cbc_trimmed_homoATCG.Ribo.bam")))
    emit(f"  {len(bams)} Ribo.bam found")
    dep_check = []
    resolved_reads = unresolved_reads = 0
    for i, b in enumerate(bams, 1):
        m = re.search(r"_(\d{3})_cbc", os.path.basename(b))
        cell = m.group(1) if m else os.path.basename(b)
        c, peak = from_bam(b, intervals, rescued5s)
        rows.append(("published", cell, c, peak))
        resolved_reads += c["5S"]
        unresolved_reads += c["other"]
        lg = b.replace(".Ribo.bam", ".ribo-map.log")
        if os.path.exists(lg):
            nreads, nunmapped, mapped = parse_ribo_log(lg)
            dep_check.append((cell, nreads, nunmapped, mapped, sum(c.values())))
        if i % 96 == 0:
            emit(f"    ... {i}/{len(bams)}")
    emit(f"  5S rescue on the published plate: {resolved_reads:,} dispersed reads "
         f"resolved as 5S, {unresolved_reads:,} left unresolved in `other`")

    # ---------------- own plate, 16 barcodes ------------------------------
    emit()
    emit("=" * 78)
    emit("OWN PLATE -- ZHA9292A1 (PM26037), mouse-only reference")
    emit("=" * 78)
    obams = sorted(glob.glob(os.path.join(OWN_DIR, "*_cbc_trimmed_homoATCG.Ribo.bam")))
    emit(f"  {len(obams)} Ribo.bam found")
    for b in obams:
        m = re.search(r"_(\d{3})_cbc", os.path.basename(b))
        cell = m.group(1) if m else os.path.basename(b)
        c, peak = from_bam(b, intervals, named5s)
        rows.append(("own", cell, c, peak))
        tot = sum(c.values())
        emit(f"    {cell}  ribo_records={tot:9,}  5'ETS={100*c['5ETS']/tot:5.1f}%  "
             f"28S={100*c['28S']/tot:5.1f}%  peak_bin={peak:5.1f}%")

    # ---------------- FLASH-seq, same code path ---------------------------
    emit()
    emit("=" * 78)
    emit("FLASH-seq -- RN26038, trimmed arm, mouse-only reference, SAME bucket()")
    emit("=" * 78)
    loci = sorted(glob.glob(os.path.join(FS_DIR, "ZHA8833A*", "*.trimmed.riboloci.tsv")))
    emit(f"  {len(loci)} riboloci.tsv found under {FS_DIR}")
    for p in loci:
        lib = re.search(r"(ZHA8833A\d+)", os.path.basename(p)).group(1)
        c, peak = from_riboloci(p, intervals, named5s)
        rows.append(("flashseq", lib, c, peak))
        tot = sum(c.values())
        emit(f"    {lib:<11} ribo_records={tot:8,}  5'ETS={100*c['5ETS']/tot:5.1f}%  "
             f"28S={100*c['28S']/tot:5.1f}%  peak_bin={peak:5.1f}%")

    # ---------------- validate the depletion table ------------------------
    emit()
    emit("=" * 78)
    emit("CHECK -- does depletion_v2_vs_v3.tsv agree with the ribo-map.logs?")
    emit("=" * 78)
    emit("  riboread-selection.py never flushes the LAST read group of a file, so")
    emit("  nunmapped + sum(mapped) == nreads - 1 by construction. Whether the")
    emit("  depletion table absorbed that one read into `ribo` changes the")
    emit("  published figure by ~1 read per cell; it is checked, not waved away.")
    if os.path.exists(DEPLETION) and dep_check:
        dep = {}
        with open(DEPLETION) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                dep[str(r["cell"]).zfill(3)] = r
        off = collections.Counter()
        sum_in = sum_ribo_log = sum_ribo_tab = sum_bam = 0
        worst = None
        for cell, nreads, nunmapped, mapped, nbam in dep_check:
            d = dep.get(cell)
            if not d:
                continue
            tab_in, tab_ribo = int(d["input_reads"]), int(d["ribo_v3"])
            off[(tab_in - nreads, tab_ribo - mapped, nbam - mapped)] += 1
            sum_in += nreads
            sum_ribo_log += mapped
            sum_ribo_tab += tab_ribo
            sum_bam += nbam
            if worst is None or abs(tab_ribo - mapped) > worst[1]:
                worst = (cell, abs(tab_ribo - mapped))
        emit(f"  cells compared: {len(dep_check)}")
        for k, n in off.most_common():
            emit(f"    (table_in - log_nreads, table_ribo - log_mapped, "
                 f"bam_records - log_mapped) = {k}  x{n} cells")
        emit(f"  plate totals: log nreads={sum_in:,}  log mapped={sum_ribo_log:,}  "
             f"table ribo_v3={sum_ribo_tab:,}  Ribo.bam records={sum_bam:,}")
        emit(f"  table-vs-log ribo excess = {sum_ribo_tab - sum_ribo_log:,} reads "
             f"({100.0*(sum_ribo_tab-sum_ribo_log)/sum_ribo_tab:.5f}% of the table's ribo)")
        emit(f"  plate-wide pct from table = {100.0*sum_ribo_tab/sum_in:.6f}%")
        emit(f"  plate-wide pct from logs  = {100.0*sum_ribo_log/sum_in:.6f}%")
    else:
        emit("  SKIPPED: depletion table or logs unavailable")

    # ---------------- write ------------------------------------------------
    out_tsv = os.path.join(outdir, "subunits_percell.tsv")
    with open(out_tsv, "w") as fh:
        cols = (["dataset", "unit", "ribo_records", "ets_peak_bin_pct"]
                + [f"n_{k}" for k in SUB_ORDER + EXTRA_ORDER]
                + [f"pct_{k}" for k in SUB_ORDER + EXTRA_ORDER]
                + [f"pct47_{k}" for k in SUB_ORDER])
        fh.write("\t".join(cols) + "\n")
        for ds, unit, c, peak in rows:
            tot = sum(c.values())
            t47 = sum(c[k] for k in SUB_ORDER)
            vals = [ds, unit, str(tot), f"{peak:.2f}"]
            vals += [str(c[k]) for k in SUB_ORDER + EXTRA_ORDER]
            vals += [f"{100.0*c[k]/tot:.4f}" if tot else "" for k in SUB_ORDER + EXTRA_ORDER]
            vals += [f"{100.0*c[k]/t47:.4f}" if t47 else "" for k in SUB_ORDER]
            fh.write("\t".join(vals) + "\n")
    emit()
    emit(f"wrote {out_tsv} ({len(rows)} rows)")

    with open(os.path.join(outdir, "subunits_report.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
