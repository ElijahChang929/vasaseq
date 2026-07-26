#!/usr/bin/env python3
"""
step3_report.py [CELLDIR] [INTERVALS_TSV]

Per-cell summary of what step 3 (in-silico rRNA depletion) did. pipeline.sh runs
this at the end of step3, so the tables land in the step log; run it by hand any
time afterwards to get the same numbers back.

Blank wells are reported exactly like real ones -- same reference, same
thresholds, no special-casing anywhere. Their numbers ARE the control, so
excluding or flagging them in the pipeline would defeat the purpose.

Inputs, both already written by the run, one pair per cell:

  <CELLDIR>/<sample>_<cell>_cbc_trimmed_homoATCG.ribo-map.log
  <CELLDIR>/<sample>_<cell>_cbc_trimmed_homoATCG.Ribo.bam

INTERVALS_TSV defaults to rrna_intervals.tsv next to $RRNA_FASTA. It maps
positions within the 47S unit to subunits; without it table 2 is skipped and
table 1 still prints.

TABLE 1 -- depletion, and which aligner found it
  in         reads out of step 2
  ribo       reads called ribosomal, i.e. records in the .Ribo.bam
  ribo%      ribo / in -- the number you actually care about
  kept       reads written to .nonRibo.fastq.gz, the input to step 4
  aln / mem / both
             which of the two bwa runs caught each ribo read. `bwa mem` refuses
             to report any alignment scoring under 30 (its -T default, match
             score 1), so a read shorter than ~30 nt is invisible to it no
             matter how well it matches -- those show up as aln-only. On this
             library ~20% of trimmed reads are under 30 nt, so aln-only is
             expected to be a large, short-read-dominated slice. That is the
             whole reason ribo-bwamem.sh runs both aligners; do not "simplify"
             it to one.
  aln_len    mean read length of the aln-only group. Expect it well under 30.

  in = ribo + kept + 1, not ribo + kept. riboread-selection.py only flushes a
  read group when it sees the NEXT read name, so the final group of each file is
  counted in its "Number of reads" total but written nowhere. One read per cell.

TABLE 2 -- where the ribosomal reads landed
  Columns are the subunits of the 47S pre-rRNA transcript, plus:
  mito       mt-Rnr1 / mt-Rnr2. Depleted here on purpose, matching the paper.
  other      the remaining Ensembl rRNA entries (dispersed 5S/5.8S/18S copies).

  A high 5'ETS + ITS share is not an error: those are only present in
  unprocessed pre-rRNA, and a total-RNA protocol is supposed to see them. It is
  a useful readout in its own right.
"""
import glob
import os
import re
import sys

SUFFIX = "_cbc_trimmed_homoATCG"


def load_intervals(path):
    """[(name, start, end)] from rrna_intervals.tsv, 1-based inclusive."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def parse_log(path):
    """-> (total, unmapped, {rg: count}) from a .ribo-map.log"""
    body = open(path, errors="replace").read()
    total = int(re.search(r"Number of reads:\s*(\d+)", body).group(1))
    unmapped = int(re.search(r"Number of unmapped reads:\s*(\d+)", body).group(1))
    rg = {k: int(v) for k, v in re.findall(r"^\t(\S+):\s*(\d+)", body, re.M)}
    return total, unmapped, rg


def parse_bam(path, intervals):
    """-> (subunit counts, mito, other, mean length of aln-only reads)"""
    try:
        import pysam
    except ImportError:
        return None, 0, 0, None
    sub = {name: 0 for name, _, _ in intervals}
    mito = other = 0
    aln_n = aln_bp = 0
    bam = pysam.AlignmentFile(path, "rb", check_sq=False)
    for r in bam.fetch(until_eof=True):
        ref = r.reference_name or ""
        if "rDNA_47S" in ref:
            pos = r.reference_start + 1  # pysam is 0-based
            for name, st, en in intervals:
                if st <= pos <= en:
                    sub[name] += 1
                    break
        elif "Mt_rRNA" in ref:
            mito += 1
        else:
            other += 1
        try:
            if r.get_tag("RG") == "aln":
                aln_n += 1
                aln_bp += r.query_length or 0
        except KeyError:
            pass
    bam.close()
    return sub, mito, other, (aln_bp / aln_n if aln_n else None)


def main():
    celldir = sys.argv[1] if len(sys.argv) > 1 else "cells"

    if len(sys.argv) > 2:
        iv_path = sys.argv[2]
    else:
        ref = os.environ.get("RRNA_FASTA", "")
        iv_path = os.path.join(os.path.dirname(ref), "rrna_intervals.tsv") if ref else ""
    intervals = load_intervals(iv_path) if iv_path and os.path.exists(iv_path) else []

    rows = []
    for log in sorted(glob.glob(os.path.join(celldir, f"*{SUFFIX}.ribo-map.log"))):
        base = os.path.basename(log)[: -len(".ribo-map.log")]
        m = re.search(r"_(\d+)" + re.escape(SUFFIX) + r"$", base)
        cell = m.group(1) if m else base
        total, unmapped, rg = parse_log(log)

        bam = os.path.join(celldir, base + ".Ribo.bam")
        sub = mito = other = aln_len = None
        if intervals and os.path.exists(bam):
            sub, mito, other, aln_len = parse_bam(bam, intervals)
        rows.append((cell, total, unmapped, rg, sub, mito, other, aln_len))

    if not rows:
        sys.exit(f"step3_report: no *{SUFFIX}.ribo-map.log under {celldir}")

    # ---- table 1 ------------------------------------------------------------
    hdr = (f"{'cell':>5} {'in':>13} {'ribo':>12} {'ribo%':>7} {'kept':>13} "
           f"{'aln':>11} {'mem':>10} {'both':>11} {'aln_len':>8}")
    print("TABLE 1 -- rRNA depletion per cell")
    print(hdr)
    print("-" * len(hdr))
    T = dict(n=0, ribo=0, kept=0, aln=0, mem=0, both=0)
    for cell, total, unmapped, rg, sub, mito, other, aln_len in rows:
        aln, mem = rg.get("aln", 0), rg.get("mem", 0)
        # mem_mem etc. are supplementary records from the same aligner
        both = sum(v for k, v in rg.items() if k not in ("aln", "mem"))
        ribo = aln + mem + both
        print(f"{cell:>5} {total:13,} {ribo:12,} {100*ribo/total:6.2f}% {unmapped:13,} "
              f"{aln:11,} {mem:10,} {both:11,} "
              f"{(f'{aln_len:8.1f}' if aln_len else '       ?')}")
        T["n"] += total; T["ribo"] += ribo; T["kept"] += unmapped
        T["aln"] += aln; T["mem"] += mem; T["both"] += both
    print("-" * len(hdr))
    print(f"{'ALL':>5} {T['n']:13,} {T['ribo']:12,} {100*T['ribo']/T['n']:6.2f}% "
          f"{T['kept']:13,} {T['aln']:11,} {T['mem']:10,} {T['both']:11,}")

    # ---- table 2 ------------------------------------------------------------
    if not intervals or rows[0][4] is None:
        print("\n(TABLE 2 skipped: no rrna_intervals.tsv, or pysam unavailable)")
        return
    names = [n for n, _, _ in intervals]
    print("\nTABLE 2 -- composition of the ribosomal reads (% of that cell's ribo)")
    hdr2 = f"{'cell':>5} " + " ".join(f"{n:>7}" for n in names) + f" {'mito':>7} {'other':>7}"
    print(hdr2)
    print("-" * len(hdr2))
    tot = {n: 0 for n in names}
    tmito = tother = 0
    for cell, total, unmapped, rg, sub, mito, other, aln_len in rows:
        d = sum(sub.values()) + mito + other
        if not d:
            continue
        print(f"{cell:>5} " + " ".join(f"{100*sub[n]/d:6.1f}%" for n in names)
              + f" {100*mito/d:6.1f}% {100*other/d:6.1f}%")
        for n in names:
            tot[n] += sub[n]
        tmito += mito; tother += other
    d = sum(tot.values()) + tmito + tother
    print("-" * len(hdr2))
    print(f"{'ALL':>5} " + " ".join(f"{100*tot[n]/d:6.1f}%" for n in names)
          + f" {100*tmito/d:6.1f}% {100*tother/d:6.1f}%")


if __name__ == "__main__":
    main()
