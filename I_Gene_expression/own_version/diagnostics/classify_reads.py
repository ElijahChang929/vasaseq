#!/usr/bin/env python3
"""
classify_reads.py [CELLDIR] [OUTFILE]

Every read of step 2, put in exactly one box. Written for the question "where
does the library actually go", which the stage-by-stage counts answer only
half of -- they say how many reads die at each pass, not what kind of read it was.

THE TWO AXES
------------
The classification is a cross of two independent things, because collapsing
them into one list leaves reads with nowhere to go (the biggest single group
here is "never read through AND died at pass 2", which a one-dimensional
scheme has no slot for):

  A. did the read run through into its own barcode?
     Recomputed with trim_bc_anchor.find_anchor on the step-1 read -- the same
     function pass 0 ran, so the answer is identical to what actually happened.

  B. where did it end up?  kept / lost at pass 1 / lost at pass 2
     Read names are matched across the three FASTQs step 2 already wrote. A
     read is in exactly one of the three sets by construction, and the script
     asserts that they sum to the input.

Reads that ran through are then split by whether pass 0 alone already left them
under the length floor. That separates "the insert was too short to survive,
nothing could have helped" from "trimming inside pass 1 shortened it", which
look identical in the stage counts.

WHAT THIS DOES NOT SPLIT
------------------------
* pass 1 is TrimGalore = adapter removal AND quality trimming in one call, so a
  read it drops cannot be attributed to one or the other from the outputs. The
  category is named after what is known ("pass 1 trimmed it short"), not after
  a guess.
* pass 2 losses are not attributed to a particular adapter. Doing that needs
  cutadapt --info-file, i.e. re-running pass 2; deliberately out of scope here.

Nothing is re-run and nothing is written into the run directory: this reads the
three FASTQs step 2 already produced and derives the rest.

  ./classify_reads.py <CELLDIR> read_classes.tsv
"""
import glob
import gzip
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trim_bc_anchor import rc, find_anchor          # noqa: E402

MINLEN = int(os.environ.get("TRIM_MINLEN", "15"))

# name -> (ran_through, fate, short_after_pass0 or None)
CLASSES = [
    ("1 read-through, kept",        (True,  "kept",  None)),
    ("2 insert too short",          (True,  "pass1", True)),
    ("3 pass 1 trimmed it short",   (None,  "pass1", False)),
    ("4 pass 2 trimmed it short",   (None,  "pass2", None)),
    ("5 no read-through, kept",     (False, "kept",  None)),
]


def read_names(path):
    """Read ids only -- the part of the header before the first ';'."""
    s = set()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                s.add(line.split(";", 1)[0][1:])
    return s


def classify_cell(celldir, sample, cell):
    fq_in = f"{celldir}/{sample}_{cell}_cbc.fastq.gz"
    fq_p1 = f"{celldir}/{sample}_{cell}_cbc_trimmed.fq.gz"
    fq_p2 = f"{celldir}/{sample}_{cell}_cbc_trimmed_homoATCG.fq.gz"
    for p in (fq_in, fq_p1, fq_p2):
        if not os.path.exists(p):
            raise SystemExit(f"classify_reads: missing {p}")

    kept = read_names(fq_p2)
    after_p1 = read_names(fq_p1)

    cnt = Counter()
    n = 0
    with gzip.open(fq_in, "rt") as f:
        while True:
            name = f.readline()
            if not name:
                break
            seq = f.readline().rstrip("\n")
            f.readline()
            f.readline()
            n += 1

            rid = name.split(";", 1)[0][1:]
            tags = dict(t.split(":", 1)
                        for t in name.rstrip("\n").split(";")[1:] if ":" in t)
            cb, rx = tags.get("CB"), tags.get("RX")

            p = find_anchor(seq, rc(cb) + rc(rx)) if (cb and rx) else -1
            ran_through = p >= 0
            # length pass 0 leaves behind; an unanchored read is untouched by it
            len_after_p0 = p if ran_through else len(seq)

            if rid in kept:
                fate = "kept"
            elif rid in after_p1:
                fate = "pass2"
            else:
                fate = "pass1"

            if fate == "kept":
                key = "1 read-through, kept" if ran_through else "5 no read-through, kept"
            elif fate == "pass2":
                key = "4 pass 2 trimmed it short"
            elif ran_through and len_after_p0 < MINLEN:
                key = "2 insert too short"
            else:
                key = "3 pass 1 trimmed it short"
            cnt[key] += 1

    assert sum(cnt.values()) == n, "classes do not sum to the input"
    return cell, n, cnt


def main():
    # `--cell NNN` does one cell only, so a caller can run all 16 in parallel
    # (find_anchor over 194M reads is ~2 h serial, ~10 min 16-way).
    args = sys.argv[1:]
    only = None
    if "--cell" in args:
        i = args.index("--cell")
        only = args[i + 1]
        del args[i:i + 2]
    celldir = args[0] if len(args) > 0 else "cells"
    outfile = args[1] if len(args) > 1 else "read_classes.tsv"

    cells = sorted(
        os.path.basename(p)[:-len("_cbc.fastq.gz")].rsplit("_", 1)
        for p in glob.glob(os.path.join(celldir, "*_cbc.fastq.gz"))
    )
    if only is not None:
        cells = [c for c in cells if c[1] == only]
        if not cells:
            raise SystemExit(f"classify_reads: no such cell {only} in {celldir}")
    if not cells:
        raise SystemExit(f"classify_reads: no *_cbc.fastq.gz under {celldir}")

    order = [k for k, _ in CLASSES]
    rows = []
    for sample, cell in cells:
        cell, n, cnt = classify_cell(celldir, sample, cell)
        sys.stderr.write(f"  {cell}: {n:,} reads\n")
        sys.stderr.flush()
        for k in order:
            rows.append((cell, k, cnt[k], 100 * cnt[k] / n if n else 0.0))

    with open(outfile, "w") as fh:
        fh.write("cell\tclass\treads\tpct\n")
        for cell, k, v, pct in rows:
            fh.write(f"{cell}\t{k}\t{v}\t{pct:.4f}\n")
    sys.stderr.write(f"wrote {outfile}  (minlen={MINLEN})\n")


if __name__ == "__main__":
    main()
