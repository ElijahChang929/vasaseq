#!/usr/bin/env python3
"""
02_published_realign.py -- make the published plate's subunit composition
actually comparable, and test whether its 5'ETS excess is an artefact.

WHY THIS EXISTS
---------------
01_vasa_subunits.py measured subunit composition for all three datasets and I
asserted in it that `pct47_*` is comparable across datasets because the 47S
record `mouse_rDNA_47S_BK000964.3_1-13403` is byte-identical in both reference
files. **That assertion was wrong, and this script exists because the data
falsified it.**

The record is identical; the COMPETITION for reads is not. The published plate
was aligned against a reference that also contains five human 45S records, and
18S/28S are strongly conserved between human and mouse while the 5'ETS spacer is
not. So the human records compete for mouse 18S/28S reads and do not compete for
mouse 5'ETS reads. Measured in 01:

    mouse-called wells: 8.47% of their ribosomal reads land on human 45S records
                        (a pure mESC well has no human rRNA -- these are
                        cross-species mismappings)
    human-called wells: 18.04% of their ribosomal reads land on the MOUSE 47S
                        record (the reverse leak, which contaminates the bucket
                        the comparison reads from)
    published mouse-cell 47S profile: 5'ETS 79.7%, 18S 1.1%, 28S 11.5%
    own-plate real-cell 47S profile:  5'ETS 18.5%, 18S 6.4%, 28S 58.8%

An 18S share of 1.1% is not a biological result; it is a reference-competition
artefact. Reporting the two profiles side by side without correcting for this
would be exactly the "present a reference difference as a protocol difference"
error the brief forbids.

WHAT THIS SCRIPT DOES
---------------------
TRACK A -- realign, to remove the competition.
  Take the reads the published plate ALREADY called ribosomal, and realign them
  against the mouse-only reference the other two datasets used, by the same
  ribo-bwamem.sh + riboread-selection.py route with the same stranded=y flag.
  The numerator is held fixed -- this asks only "where do these same reads sit
  on the same reference the other two were measured on", which is a composition
  question, not a depletion question. It therefore does NOT produce a new rRNA
  percentage and this script deliberately does not report one.

  Restricted to the 172 fig1d mouse-called wells: they are the only ones
  comparable to the other two datasets, and it keeps the job small.

  Reads that fail to map to the mouse-only reference are counted, not dropped
  silently -- that count is the human contribution plus GRCm38/GRCm39 dispersed
  differences, and it is reported.

TRACK B -- the 5'ETS positional histogram, for all three datasets.
  05_rrna_bwa_report.py records that a poly-T population once inflated 5'ETS to
  31-62% of a blank cell's ribosomal calls, caught because 88.9% of those hits
  sat in ONE 200 nt window. 01 found the published plate's 5'ETS peak-bin share
  is median 27.6% (48/172 mouse wells above 30%) against 9.9-14.0% for the own
  plate's real cells. That is the artefact signature, so the position of the
  peak has to be looked at rather than inferred. This track writes the full
  200 nt-binned 5'ETS profile per dataset so the peaks can be compared by
  LOCATION, not just by height.

  Read length is a third difference and is reported alongside: the published
  plate is 74 nt (v3_driver.log), the own plate and FLASH-seq 151 nt. A shorter
  read maps to a conserved subunit more promiscuously, so read length and
  reference competition push the same direction and are not separable by this
  design. Stated, not resolved.

USAGE
    02_published_realign.py <outdir>
"""
import collections
import glob
import os
import re
import subprocess
import sys

import pysam

ROOT = "/nemo/lab/turnerj/working/guangxin/vasaseq"
REFS = "/nemo/lab/turnerj/working/guangxin/reference/vasaseq"
V3_DIR = os.environ.get("VP_V3_DIR", f"{ROOT}/data/ref/fastq_vasaplate/vasaplate_out_v3")
OWN_DIR = os.environ.get("OWN_CELLS_DIR", f"{ROOT}/data/PM26037/out/cells")
FS_DIR = os.environ.get("FS_RRNA_BWA_DIR", f"{ROOT}/res/flashseq/rrna_bwa")
MOUSE_FA = os.environ.get("MOUSE_FA", f"{REFS}/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa")
INTERVALS = os.environ.get("INTERVALS", f"{REFS}/mouse_GRCm39_E116/rrna_intervals.tsv")
MOUSE_CELLS = os.environ.get("MOUSE_CELLS", "")   # comma-separated 3-digit well ids

ETS_BIN = 200
SUB_ORDER = ["5ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "3ETS"]
BWA = os.environ["BWA"]
SAMTOOLS = os.environ["SAMTOOLS"]


def load_intervals(path):
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        out.append((f[0], int(f[1]), int(f[2])))
    return out


def sub_of(pos, intervals):
    for name, st, en in intervals:
        if st <= pos <= en:
            return name
    return "47S_outside"


def run(cmd, **kw):
    """Never let a silent tool failure become a NaN -- check=True, always."""
    return subprocess.run(cmd, shell=True, check=True, text=True,
                          capture_output=True, **kw)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(outdir, exist_ok=True)
    intervals = load_intervals(INTERVALS)
    log = []

    def emit(s=""):
        print(s, flush=True)
        log.append(s)

    wells = [w.strip().zfill(3) for w in MOUSE_CELLS.split(",") if w.strip()]
    emit(f"mouse-called wells supplied: {len(wells)}")
    assert wells, "MOUSE_CELLS is empty -- nothing to realign"

    # ---------------- TRACK A ---------------------------------------------
    emit()
    emit("=" * 78)
    emit("TRACK A -- realign published mouse-well ribosomal reads to the")
    emit("           MOUSE-ONLY reference, removing human-record competition")
    emit("=" * 78)

    fq = os.path.join(outdir, "pub_mouse_ribo.fastq")
    n_in = 0
    with open(fq, "w") as fh:
        for w in wells:
            hits = glob.glob(os.path.join(V3_DIR, f"*_{w}_cbc_trimmed_homoATCG.Ribo.bam"))
            if not hits:
                emit(f"  ! no Ribo.bam for well {w}")
                continue
            bam = pysam.AlignmentFile(hits[0], "rb", check_sq=False)
            for r in bam.fetch(until_eof=True):
                seq = r.query_sequence
                qual = r.qual
                if not seq:
                    continue
                # reads were stored as aligned; restore original orientation
                if r.is_reverse:
                    seq = seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
                    qual = qual[::-1] if qual else qual
                fh.write(f"@{r.qname}_{w}\n{seq}\n+\n{qual or 'I'*len(seq)}\n")
                n_in += 1
            bam.close()
    emit(f"  extracted {n_in:,} already-ribosomal reads from {len(wells)} mouse wells")
    assert n_in > 0, "no reads extracted"

    # same two aligners ribo-bwamem.sh uses, same reference, same reasons:
    # bwa mem cannot report an alignment scoring under its -T 30 default, so
    # short reads are invisible to it and aln is needed to catch them.
    sam_aln = os.path.join(outdir, "aln.sam")
    sam_mem = os.path.join(outdir, "mem.sam")
    run(f'{BWA} aln -t 8 "{MOUSE_FA}" "{fq}" > "{outdir}/aln.sai" 2>"{outdir}/bwa_aln.err"')
    run(f'{BWA} samse "{MOUSE_FA}" "{outdir}/aln.sai" "{fq}" > "{sam_aln}" 2>>"{outdir}/bwa_aln.err"')
    run(f'{BWA} mem -t 8 "{MOUSE_FA}" "{fq}" > "{sam_mem}" 2>"{outdir}/bwa_mem.err"')

    # riboread-selection.py's stranded=y predicate: a read is ribosomal on this
    # reference if it has at least one FORWARD mapped alignment. Group by qname
    # across both aligners, exactly as the pipeline does over its merged BAM.
    best = {}
    unmapped_any = set()
    for sam in (sam_aln, sam_mem):
        af = pysam.AlignmentFile(sam, "r", check_sq=False)
        for r in af.fetch(until_eof=True):
            if r.is_unmapped:
                unmapped_any.add(r.qname)
                continue
            if r.is_reverse:
                continue                      # stranded=y drops reverse hits
            ref = r.reference_name or ""
            pos = r.reference_start + 1
            # prefer the 47S record when a read maps to both it and a dispersed
            # entry, so composition is not decided by aligner tie-breaking
            key = (0 if "rDNA_47S" in ref else 1, -(r.query_alignment_length or 0))
            if r.qname not in best or key < best[r.qname][0]:
                best[r.qname] = (key, ref, pos)
        af.close()

    c = collections.Counter()
    ets_bins = collections.Counter()
    for qname, (_, ref, pos) in best.items():
        if "rDNA_47S" in ref:
            nm = sub_of(pos, intervals)
            c[nm] += 1
            if nm == "5ETS":
                ets_bins[(pos // ETS_BIN) * ETS_BIN] += 1
        elif "Mt_rRNA" in ref:
            c["mito_mouse"] += 1
        else:
            c["other"] += 1
    n_fwd = len(best)
    n_nofwd = n_in - n_fwd
    emit(f"  reads with >=1 FORWARD hit on the mouse-only reference: {n_fwd:,} "
         f"({100.0*n_fwd/n_in:.2f}%)")
    emit(f"  reads with none: {n_nofwd:,} ({100.0*n_nofwd/n_in:.2f}%) -- human rRNA")
    emit(f"      plus GRCm38-vs-GRCm39 dispersed-entry differences and reverse-only hits")
    t47 = sum(c[k] for k in SUB_ORDER)
    emit()
    emit(f"  47S-derived reads: {t47:,}")
    emit("  REALIGNED published mouse-well 47S composition (pct of 47S-derived):")
    for k in SUB_ORDER:
        emit(f"      {k:<6} {c[k]:>9,}  {100.0*c[k]/t47:6.2f}%" if t47 else f"      {k}: n/a")
    emit(f"      [non-47S: mito_mouse={c['mito_mouse']:,}  other={c['other']:,}]")
    pk = 100.0 * max(ets_bins.values()) / sum(ets_bins.values()) if ets_bins else 0.0
    emit(f"  5'ETS peak-bin share after realignment: {pk:.2f}%")

    with open(os.path.join(outdir, "published_realigned_47S.tsv"), "w") as fh:
        fh.write("subunit\tn\tpct47\n")
        for k in SUB_ORDER:
            fh.write(f"{k}\t{c[k]}\t{100.0*c[k]/t47:.4f}\n" if t47 else f"{k}\t0\t\n")
        fh.write(f"__n_input\t{n_in}\t\n__n_fwd_mouse\t{n_fwd}\t\n__n_nofwd\t{n_nofwd}\t\n")
        fh.write(f"__mito_mouse\t{c['mito_mouse']}\t\n__other\t{c['other']}\t\n")

    for f in (sam_aln, sam_mem, fq, f"{outdir}/aln.sai"):
        try:
            os.remove(f)
        except OSError:
            pass

    # ---------------- TRACK B ---------------------------------------------
    emit()
    emit("=" * 78)
    emit("TRACK B -- 5'ETS positional profile, by LOCATION not just height")
    emit("=" * 78)
    prof = {}

    def add_bam(tag, paths, restrict=None):
        b = collections.Counter()
        for p in paths:
            if restrict:
                m = re.search(r"_(\d{3})_cbc", os.path.basename(p))
                if not m or m.group(1) not in restrict:
                    continue
            bam = pysam.AlignmentFile(p, "rb", check_sq=False)
            for r in bam.fetch(until_eof=True):
                ref = r.reference_name or ""
                if "rDNA_47S" in ref:
                    pos = r.reference_start + 1
                    if sub_of(pos, intervals) == "5ETS":
                        b[(pos // ETS_BIN) * ETS_BIN] += 1
            bam.close()
        prof[tag] = b

    add_bam("published_mouse",
            sorted(glob.glob(os.path.join(V3_DIR, "*_cbc_trimmed_homoATCG.Ribo.bam"))),
            restrict=set(wells))
    own = sorted(glob.glob(os.path.join(OWN_DIR, "*_cbc_trimmed_homoATCG.Ribo.bam")))
    add_bam("own_real", [p for p in own
                         if re.search(r"_(\d{3})_cbc", p).group(1) not in
                         {"001", "014", "015", "016"}])
    add_bam("own_blank", [p for p in own
                          if re.search(r"_(\d{3})_cbc", p).group(1) in
                          {"001", "014", "015", "016"}])
    fsb = collections.Counter()
    for p in sorted(glob.glob(os.path.join(FS_DIR, "ZHA8833A*", "*.trimmed.riboloci.tsv"))):
        with open(p) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and "rDNA_47S" in f[1]:
                    pos = int(f[2])
                    if sub_of(pos, intervals) == "5ETS":
                        fsb[(pos // ETS_BIN) * ETS_BIN] += 1
    prof["flashseq"] = fsb
    prof["published_mouse_realigned"] = ets_bins

    with open(os.path.join(outdir, "ets_profile.tsv"), "w") as fh:
        fh.write("dataset\tbin_start\tn\tpct_of_5ETS\n")
        for tag, b in prof.items():
            t = sum(b.values())
            for k in sorted(b):
                fh.write(f"{tag}\t{k}\t{b[k]}\t{100.0*b[k]/t:.4f}\n")
    emit(f"  wrote ets_profile.tsv")
    emit()
    emit(f"  {'dataset':<28} {'n_5ETS':>10} {'peak_bin':>9} {'peak%':>7} {'top3 bins'}")
    for tag, b in prof.items():
        t = sum(b.values())
        if not t:
            emit(f"  {tag:<28} {'0':>10}")
            continue
        top = b.most_common(3)
        emit(f"  {tag:<28} {t:>10,} {top[0][0]:>9} {100.0*top[0][1]/t:6.2f}%  "
             + ", ".join(f"{k}({100.0*v/t:.1f}%)" for k, v in top))
    emit()
    emit("  Same peak LOCATION across datasets => a real feature of the 5'ETS")
    emit("  (or of the reference). A peak unique to one dataset => that dataset's")
    emit("  own artefact. This is what distinguishes the two.")

    with open(os.path.join(outdir, "realign_report.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
