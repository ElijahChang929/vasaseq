#!/usr/bin/env python3
"""
trim_bc_anchor.py  IN.fastq.gz  OUT.fastq.gz  [--log FILE]

Step 2, pass 0. Cut the 3' technical tail off a demultiplexed VASA read by
FINDING it, not by recognising its shape.

The idea
--------
The 3' construct is

    [insert][poly-A][revcomp(CBC)][revcomp(UMI)][adapter ...]

and after step 1 both the cell barcode and the UMI are on the read name
(CB: and RX:). So for THIS read the 12 nt that follow the poly-A are known
exactly -- not a pattern, a literal string. Find it, throw away everything from
there to the 3' end, then walk back over the poly-A that precedes it.

Why not just do this in cutadapt
--------------------------------
cutadapt cannot take a per-read adapter. Doing it there means the UMI becomes
6 wildcards, and wildcards must run at min_overlap = the pattern's full length
or they eat the end of every read -- which then forces the adapter to be
present too, so the pattern only fires when ~21 nt of adapter was sequenced.
Here the anchor is 12 fully specific bases and needs no adapter at all, so it
also fires on reads whose tail is truncated. Chance hit probability is 4^-12,
about 2 expected false hits per 300,000 reads.

What it does NOT do
-------------------
No quality trimming, no adapter trimming, no length filtering. Reads with no
anchor pass through untouched. This runs BEFORE TrimGalore and the regular
cutadapt pass, which still handle everything else -- there is one place that
decides minimum length and it is not here.
"""
import argparse
import gzip
import sys

# xopen arrives with cutadapt and wraps isal/igzip, which is several times
# faster than the stdlib gzip at both ends. Fall back if it is ever absent --
# the output is identical either way, only the wall time changes.
try:
    from xopen import xopen

    def _open(path, mode):
        return xopen(path, mode, threads=2, compresslevel=1)
except ImportError:  # pragma: no cover
    def _open(path, mode):
        return gzip.open(path, mode, compresslevel=1)

RC = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(RC)[::-1]


def find_anchor(seq, anchor, max_mm=1, min_partial=4):
    """Leftmost position where `anchor` starts, or -1.

    Three cases, in order of confidence:
      1. exact match anywhere
      2. <=max_mm mismatches anywhere, anchored on an exact match of the first
         half so the scan stays cheap
      3. the read ENDS partway through the anchor -- a suffix of the read that
         is a prefix of the anchor, >=min_partial long. This is the case the
         cutadapt version structurally cannot catch, because it needs the whole
         pattern to fit inside the read.
    """
    p = seq.find(anchor)
    if p >= 0:
        return p

    n = len(anchor)
    if max_mm:
        half = anchor[: n // 2]
        start = 0
        while True:
            i = seq.find(half, start)
            if i < 0 or i + n > len(seq):
                break
            mm = sum(1 for a, b in zip(seq[i:i + n], anchor) if a != b)
            if mm <= max_mm:
                return i
            start = i + 1

    for k in range(min(n - 1, len(seq)), min_partial - 1, -1):
        if seq.endswith(anchor[:k]):
            return len(seq) - k
    return -1


def strip_polya(seq, min_run=3):
    """Remove a trailing poly-A run, tolerating sequencing errors.

    Scored like cutadapt's poly-A trimmer: walking left from the 3' end, +1 for
    an A and -2 for anything else, cut at the best score. That way a run of 20 A
    with one miscalled base is still removed whole, while a genuine A-ending
    coding sequence is left alone.
    """
    best = score = 0
    cut = len(seq)
    for i in range(len(seq) - 1, -1, -1):
        score += 1 if seq[i] == "A" else -2
        if score > best:
            best, cut = score, i
    return seq[:cut] if len(seq) - cut >= min_run else seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--log", default=None)
    ap.add_argument("--max-mm", type=int, default=1,
                    help="mismatches allowed in the 12 nt anchor (default 1)")
    ap.add_argument("--min-partial", type=int, default=4,
                    help="shortest anchor prefix accepted at the read end")
    a = ap.parse_args()

    n = hit = exact = polya = 0
    lost = 0
    with _open(a.infile, "rt") as fi, _open(a.outfile, "wt") as fo:
        while True:
            name = fi.readline()
            if not name:
                break
            seq = fi.readline().rstrip("\n")
            plus = fi.readline()
            qual = fi.readline().rstrip("\n")
            n += 1

            tags = dict(t.split(":", 1) for t in name.rstrip("\n").split(";")[1:] if ":" in t)
            cb, rx = tags.get("CB"), tags.get("RX")
            if cb and rx:
                anchor = rc(cb) + rc(rx)
                p = find_anchor(seq, anchor, a.max_mm, a.min_partial)
                if p >= 0:
                    hit += 1
                    if seq[p:p + len(anchor)] == anchor:
                        exact += 1
                    lost += len(seq) - p
                    seq, qual = seq[:p], qual[:p]
                    t = strip_polya(seq)
                    if len(t) != len(seq):
                        polya += 1
                        seq, qual = t, qual[:len(t)]
            fo.write(f"{name}{seq}\n{plus}{qual}\n")

    msg = (f"reads {n}\n"
           f"anchor found {hit} ({hit/n:.2%} of reads), of which exact {exact}\n"
           f"poly-A also stripped {polya} ({polya/n:.2%})\n"
           f"bases removed at the anchor {lost}\n")
    sys.stderr.write(msg)
    if a.log:
        with open(a.log, "w") as fh:
            fh.write(msg)


if __name__ == "__main__":
    main()
