#!/usr/bin/env python3
"""
attribute_pass2_loss.py --cell NNN CELLDIR OUTFILE

Which adapter is responsible for each read class 4 loses?

classify_reads.py says 44.9% of the library dies at pass 2. That is the single
biggest thing that happens to this library and "cutadapt -m" is not an
explanation, so this asks cutadapt which adapter did the cutting, per read.

HOW
---
cutadapt --info-file writes one line per adapter match per round (with -n 10 a
read can produce several). The columns are:

  1 read name  2 errors  3 start  4 end  5 before  6 match  7 after
  8 adapter name  9-12 qualities

Bases a match removes depends on the adapter's end:
  3' adapter (rt, polyA, polyG) removes the match and everything AFTER it
  5' adapter (polyT5)           removes the match and everything BEFORE it

The read is credited to whichever adapter removed the most bases from it; a read
where no adapter matched at all is credited to "no adapter" (it was shortened by
--poly-a or --trim-n, or arrived already short).

WHY IT STREAMS
--------------
The info-file is ~490 bytes per read, so the full library would be ~90 GB on
disk. cutadapt writes into a FIFO and this reads it as it goes; nothing is
stored. Measured on cell 016: 730k reads -> 1.56M lines, 355 MB, 25 s.

The set of reads that pass 2 dropped is not recomputed -- it is the difference
between the read names of pass 1's output and pass 2's output, both already on
disk, which is exactly what classify_reads.py calls class 4.

  ./attribute_pass2_loss.py --cell 016 <CELLDIR> out_016.tsv
"""
import gzip
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

CUTADAPT = os.environ.get(
    "TRIM_CUTADAPT", "/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/cutadapt")
MINLEN = os.environ.get("TRIM_MINLEN", "15")
ADAPTER3 = os.environ.get(
    "TRIM_ADAPTER3", "GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT")

# Two pass-2s exist and they share nothing but the -m flag, so the mode has to
# be explicit rather than inferred:
#
#   vasa    own_version/trim.sh -- measured read-through adapter, 20-mer
#           homopolymers with min_overlap, a 5' poly-T, --poly-a, -n 10
#   legacy  a_Mapping/trim.sh (the published pipeline, and what the VASA-plate
#           library was actually run with) -- four 6-mer homopolymers, default
#           -n 1, no --poly-a, no 5' adapter. Verified against the run's own
#           cutadapt log before being written down here.
#
# `ends` says which side each adapter trims from, which is how bases-removed is
# counted: a 3' adapter takes the match and everything after, a 5' adapter takes
# the match and everything before.
MODES = {
    "vasa": dict(
        args=lambda: ["-m", MINLEN, "--trim-n", "-n", "10", "--poly-a",
                      "-a", f"rt={ADAPTER3};min_overlap=8",
                      "-a", "polyA=A{20};min_overlap=10",
                      "-a", "polyG=G{20};min_overlap=10",
                      "-g", "polyT5=T{20};min_overlap=10"],
        ends={"rt": 3, "polyA": 3, "polyG": 3, "polyT5": 5},
    ),
    "legacy": dict(
        args=lambda: ["-m", MINLEN, "--trim-n",
                      "-a", "polyG1=GG{5}", "-a", "polyC1=CC{5}",
                      "-a", "polyT1=TT{5}", "-a", "polyA1=AA{5}"],
        ends={"polyG1": 3, "polyC1": 3, "polyT1": 3, "polyA1": 3},
    ),
}


def read_names(path):
    s = set()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                s.add(line.split(";", 1)[0][1:])
    return s


def main():
    args = sys.argv[1:]
    cell = None
    mode = "vasa"
    if "--cell" in args:
        i = args.index("--cell")
        cell = args[i + 1]
        del args[i:i + 2]
    if "--mode" in args:
        i = args.index("--mode")
        mode = args[i + 1]
        del args[i:i + 2]
    if mode not in MODES:
        raise SystemExit(f"--mode must be one of {sorted(MODES)}")
    if cell is None or len(args) < 2:
        raise SystemExit(__doc__)
    celldir, outfile = args[0], args[1]
    ENDS = MODES[mode]["ends"]

    hits = [p for p in os.listdir(celldir)
            if p.endswith(f"_{cell}_cbc_trimmed.fq.gz")]
    if len(hits) != 1:
        raise SystemExit(f"attribute_pass2_loss: cannot find cell {cell} in {celldir}")
    sample = hits[0][: -len(f"_{cell}_cbc_trimmed.fq.gz")]
    fq_p1 = os.path.join(celldir, f"{sample}_{cell}_cbc_trimmed.fq.gz")
    fq_p2 = os.path.join(celldir, f"{sample}_{cell}_cbc_trimmed_homoATCG.fq.gz")

    dropped = read_names(fq_p1) - read_names(fq_p2)
    sys.stderr.write(f"  {cell}: {len(dropped):,} reads dropped at pass 2\n")
    sys.stderr.flush()

    tmpdir = tempfile.mkdtemp()
    fifo = os.path.join(tmpdir, "info")
    os.mkfifo(fifo)
    cmd = ([CUTADAPT] + MODES[mode]["args"]()
           + ["--info-file", fifo, "-o", os.devnull, fq_p1])
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL, env=env)

    dominant = Counter()      # dropped reads, by the adapter that cut the most
    bases = Counter()         # total bases each adapter removed from them
    cur_name, cur = None, defaultdict(int)
    cur_dropped = False       # only dropped reads are counted -- a kept read has
                              # an empty `cur` too and would otherwise be filed
                              # under "no_adapter", inflating it to every read

    def flush():
        if cur_name is None or not cur_dropped:
            return
        if cur:
            top = max(cur, key=lambda k: cur[k])
            dominant[top] += 1
            for k, v in cur.items():
                bases[k] += v
        else:
            dominant["no_adapter"] += 1

    with open(fifo) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name = f[0].split(";", 1)[0]
            if name != cur_name:
                flush()
                cur_name, cur = name, defaultdict(int)
                cur_dropped = name in dropped
            if not cur_dropped:
                continue
            if f[1] == "-1":                 # no adapter matched this round
                continue
            ad = f[7]
            before, match, after = f[4], f[5], f[6]
            cur[ad] += (len(match) + len(after)) if ENDS.get(ad, 3) == 3 \
                else (len(before) + len(match))
    flush()
    proc.wait()
    if proc.returncode != 0:
        raise SystemExit(f"attribute_pass2_loss: cutadapt failed on cell {cell}")
    os.remove(fifo)
    os.rmdir(tmpdir)

    total = sum(dominant.values())
    if total != len(dropped):
        sys.stderr.write(f"  WARNING {cell}: attributed {total:,} of "
                         f"{len(dropped):,} dropped reads\n")
    with open(outfile, "w") as fh:
        fh.write("cell\tadapter\treads\tpct_of_dropped\tbases_removed\n")
        for ad, n in dominant.most_common():
            fh.write(f"{cell}\t{ad}\t{n}\t{100*n/total if total else 0:.4f}"
                     f"\t{bases.get(ad, 0)}\n")


if __name__ == "__main__":
    main()
