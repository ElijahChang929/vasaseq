#!/usr/bin/env python3
"""trimview - IGV-like side-by-side view of what each pipeline stage does to a read.

For one cell, the mapping pipeline writes a chain of FASTQs where every stage
only ever *shortens* a read or *removes* it entirely:

    _cbc.fastq.gz                     step1  extract  (barcode moved to read name)
    _cbc_trimmed.fq.gz                step2a TrimGalore (adapter + quality, 3')
    _cbc_trimmed_homoATCG.fq.gz       step2b cutadapt homopolymers (5' and 3')
    _cbc_trimmed_homoATCG.nonRibo.*   step3  rRNA depletion (read kept or dropped)

This shows the same read across all of them, aligned in the coordinate frame of
the original untrimmed read, so the trimmed-away bases stay visible (dimmed).

Above the sequence runs a feature track marking the technical parts of the read:
the 3' construct is [insert][poly-A][rc(CBC)][rc(UMI)][adapter][Illumina tail],
and both the barcode and the UMI are on the read name, so their position is
known exactly rather than guessed. The anchor is located by calling
trim_bc_anchor.py's own find_anchor(), so what is drawn is what step 2 pass 0
actually searched for. Reverse-orientation reads carry the same construct
revcomp'd at the 5' end and are labelled as such.

The adapter side is matched the way cutadapt matches it -- whole or cut short
by the read end, with mismatches allowed -- and carries on past TRIM_ADAPTER3
through the rest of the P5 arm (Read1 primer, i5 index, P5), so no part of a
read that ran into the adapter is left looking like sequence. Everything that
is adapter is struck through in the sequence rows as well as marked on the
track.

Order: the first three files are in identical read order, so they are streamed
together as a merge join -- no seeking, no index, constant memory, and paging
forward is instant.  The nonRibo file is name-sorted by samtools and therefore
in a *different* order, so presence there is answered from a hash index that is
built once per cell and cached under ~/.cache/trimview/.

Usage (from anywhere; --dir defaults to the current directory):

    ./trimview.py 001                 # interactive pager, 3 reads at a time
    ./trimview.py 001 -n 5            # 5 reads per page
    ./trimview.py 001 --skip 10000    # start further into the file
    ./trimview.py 001 --dropped       # only reads that get dropped somewhere
    ./trimview.py 001 --find 12345:1083
    ./trimview.py 001 --stats 200000  # no alignment view, just the numbers
    ./trimview.py --list              # which cells are available

Interactive keys: Enter/n next page, n<k> k reads, s<k> skip k, f<txt> find,
d toggle dropped-only, t toggle trimmed-only, c<cell> switch cell, q quit.
"""

import argparse
import bisect
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from array import array

# ---------------------------------------------------------------------------
# stage definitions: (short label, filename suffix, stream|indexed)
# ---------------------------------------------------------------------------

STAGES = [
    ("cbc",      "_cbc.fastq.gz",                            "stream"),
    ("trimmed",  "_cbc_trimmed.fq.gz",                       "stream"),
    ("homoATCG", "_cbc_trimmed_homoATCG.fq.gz",              "stream"),
    ("nonRibo",  "_cbc_trimmed_homoATCG.nonRibo.fastq.gz",   "indexed"),
]

# What each stage does, and the only way it can make a read disappear.
STAGE_DOC = {
    "cbc":      ("step1 extract: barcode+UMI moved onto the read name", None),
    "trimmed":  ("step2 pass0 barcode anchor + pass1 TrimGalore",
                 "shorter than TRIM_MINLEN after trimming"),
    "homoATCG": ("step2 pass2 cutadapt: read-through, poly-A/G 3', poly-T 5'",
                 "shorter than TRIM_MINLEN after trimming"),
    "nonRibo":  ("step3 rRNA depletion (bwa aln + bwa mem)",
                 "mapped to the rRNA reference"),
}

# Assignment stage (step5) -- read name -> gene, looked up on demand.
BED_STAGES = [
    ("single", "_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out"
               ".singlemappers_genes.bed.gz"),
    ("multi",  "_cbc_trimmed_homoATCG.nonRibo_E99_Aligned.out.nsorted"
               ".multimappers_genes.bed.gz"),
]

# The stock adapters TrimGalore auto-detects. On this library pass 1 picks
# *Nextera* (`CTGTCTCTTATA`, 30.9% of cell 002's reads -- see its
# *_trimming_report.txt), which is no coincidence: it is nt 22-40 of the
# read-through construct below, so most of these hits are that construct seen
# without its first 21 nt.
ADAPTERS = [
    ("TruSeq/Illumina", "AGATCGGAAGAGC"),
    ("Nextera",         "CTGTCTCTTATACACATCT"),
    ("smallRNA",        "TGGAATTCTCGG"),
]

# The 3' read-through construct this library carries (trim_bc_anchor.py):
#     [insert][poly-A][rc(CBC)][rc(UMI)][adapter ...]
# Default adapter is config.sh's TRIM_ADAPTER3; override with --adapter.
DEFAULT_ADAPTER3 = "GATCGTCGGACTGTAGAACTCCTGTCTCTTATACACATCT"
MIN_ADAPTER_HIT = 8          # shortest adapter match worth marking on its own
ADAPTER_ERR = 0.1            # mismatches tolerated, as cutadapt's -e default

# TRIM_ADAPTER3 stops after 40 nt, but the molecule does not: a read that runs
# into it keeps going through the rest of the P5 arm, which cutadapt never
# names and which would otherwise be drawn as if it were sequence.
#
# Measured over the first 400,000 reads of cell 002 (54,860 of them reach the
# adapter): the 53 nt following TRIM_ADAPTER3 are a fixed tail -- one consensus
# base per column, 82% pure at the worst column and 95% at the median, on
# >=15,457 reads -- and it reads as rc(Read1 primer) + rc(i5 index, 10 nt) +
# rc(P5). The two primer literals are constant for any Illumina library (and
# are matched with mismatches allowed, hence the 82%); the 10 nt between them
# is this library's index, so it is found as the gap rather than matched.
ADAPTER_TAIL = [
    ("rc(Read1 primer)", "GACGCTGCCGACGA"),
    ("rc(P5)",           "GTGTAGATCTCGGTGGTCGCCGTATCATT"),
]
MAX_INDEX_GAP = 16           # widest gap between tail parts still called index

# feature track characters
F_CBC, F_UMI, F_ADAPT, F_POLY = "C", "U", "=", "~"
F_CBC_ALONE, F_UMI_ALONE = "c", "u"
F_ADAPT_PART, F_TECH = ":", "+"

# inline styling of the base rows: barcode and UMI as before, and everything
# that is adapter struck through, so a technical stretch is obvious in the
# sequence itself and not only on the track. (SGR 9 is ignored by a few old
# terminals -- there the feature track still carries it.)
MARK_STYLE = {F_CBC: "7;", F_UMI: "4;", F_CBC_ALONE: "7;", F_UMI_ALONE: "4;",
              F_ADAPT: "9;", F_ADAPT_PART: "9;", F_TECH: "9;"}

# short names for the summary line, so a cut can be described by what it
# covered rather than guessed at from the fragment
MARK_NAME = [(F_POLY, "poly-A"), (F_CBC, "barcode"), (F_UMI, "UMI"),
             (F_CBC_ALONE, "barcode?"), (F_UMI_ALONE, "UMI?"),
             (F_ADAPT, "adapter"), (F_ADAPT_PART, "adapter"),
             (F_TECH, "Illumina tail")]

# trim_bc_anchor.py is the pipeline's own pass 0. Import its functions rather
# than reimplementing them, so the view cannot drift away from what actually
# ran. Falls back to a local copy of rc() if the file ever moves.
# realpath, not abspath: this is normally invoked through the `trimview`
# symlink in $OUTDIR/cells, and abspath would resolve to that directory.
# trim_bc_anchor.py lives one level up, in own_version/, so the parent goes on
# the path too.
_HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(1, os.path.dirname(_HERE))
try:
    from trim_bc_anchor import rc, find_anchor
    HAVE_ANCHOR = True
except ImportError:                                     # pragma: no cover
    HAVE_ANCHOR = False
    _RC = str.maketrans("ACGTN", "TGCAN")

    def rc(s):
        return s.translate(_RC)[::-1]

def _polya_cut(seq, min_run=3):
    """Where cutadapt's --poly-a would cut. Display only -- pass 2 does the real
    trimming.

    Both of cutadapt's rules, unlike the half-implementation this replaced:
      1. score each suffix, +1 per A and -2 per non-A, take the maximum
      2. exclude any suffix that is more than 20% non-A
    Rule 2 is the one that matters: scoring alone permits up to 33% non-A, which
    is how the old version ate 89 nt reads down to 3.
    """
    best = score = 0
    cut = len(seq)
    n_a = n_other = 0
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == "A":
            n_a += 1
            score += 1
        else:
            n_other += 1
            score -= 2
        if score > best and n_other <= 0.2 * (n_a + n_other):
            best, cut = score, i
    return seq[:cut] if len(seq) - cut >= min_run else seq


CACHE_DIR = os.environ.get("TRIMVIEW_CACHE",
                           os.path.expanduser("~/.cache/trimview"))

# ---------------------------------------------------------------------------
# colour
# ---------------------------------------------------------------------------


class Palette(object):
    def __init__(self, enabled):
        self.on = enabled

    def _w(self, code, s):
        if not self.on:
            return s
        return "\033[%sm%s\033[0m" % (code, s)

    def dim(self, s):
        return self._w("2;90", s)

    def bold(self, s):
        return self._w("1", s)

    def head(self, s):
        return self._w("1;36", s)

    def good(self, s):
        return self._w("32", s)

    def bad(self, s):
        return self._w("31", s)

    def warn(self, s):
        return self._w("33", s)

    def label(self, s):
        return self._w("36", s)

    # IGV-ish base colours
    BASE = {"A": "32", "C": "34", "G": "33", "T": "31", "N": "35"}

    def base(self, ch, style=""):
        """One base, coloured by nucleotide; `style` prefixes SGR attributes
        (reverse video for the cell barcode, underline for the UMI)."""
        if not self.on:
            return ch
        return self._w(style + self.BASE.get(ch.upper(), "37"), ch)

    def bases(self, s):
        if not self.on:
            return s
        return "".join(self.base(c) for c in s)

    def dim_base(self, ch, style=""):
        # trimmed-away sequence: keep it readable but clearly secondary
        if not self.on:
            return ch.lower()
        return self._w(style + "2;90", ch.lower())

    def dim_bases(self, s):
        return self.dim(s.lower())


# ---------------------------------------------------------------------------
# fastq streaming
# ---------------------------------------------------------------------------


def open_fastq(path):
    """Return a text line iterator plus a closer, using pigz/zcat when present."""
    prog = shutil.which("pigz") or shutil.which("zcat")
    if prog:
        args = [prog, "-dc", path] if prog.endswith("pigz") else [prog, path]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL)

        def closer():
            try:
                proc.stdout.close()
                proc.kill()
                proc.wait()
            except Exception:
                pass

        return (l.decode("ascii", "replace") for l in proc.stdout), closer

    fh = gzip.open(path, "rt")
    return fh, fh.close


class FastqStream(object):
    """One-record lookahead over a FASTQ, so several can be merge-joined."""

    def __init__(self, path):
        self.path = path
        self.lines, self._close = open_fastq(path)
        self.rec = None          # (name, seq, qual)
        self.consumed = 0
        self.eof = False
        self._read()

    def _read(self):
        try:
            name = next(self.lines)
            seq = next(self.lines)
            next(self.lines)          # '+'
            qual = next(self.lines)
        except StopIteration:
            self.rec = None
            self.eof = True
            return
        self.rec = (name[1:].rstrip("\n"), seq.rstrip("\n"), qual.rstrip("\n"))

    def take(self):
        r = self.rec
        self.consumed += 1
        self._read()
        return r

    def close(self):
        self._close()


# ---------------------------------------------------------------------------
# nonRibo presence index (that file is name-sorted, not in pipeline order)
# ---------------------------------------------------------------------------


def _h64(name):
    return int.from_bytes(hashlib.blake2b(name.encode(), digest_size=8).digest(),
                          "little")


class NameIndex(object):
    """Sorted array of 64-bit name hashes; ~8 bytes/read, cached on disk."""

    def __init__(self, path, quiet=False):
        self.path = path
        self.arr = array("Q")
        st = os.stat(path)
        key = hashlib.blake2b(
            ("%s|%d|%d" % (os.path.realpath(path), st.st_size,
                           int(st.st_mtime))).encode(), digest_size=12).hexdigest()
        self.cache = os.path.join(CACHE_DIR, key + ".idx")
        if os.path.exists(self.cache):
            with open(self.cache, "rb") as fh:
                self.arr.frombytes(fh.read())
        else:
            self._build(quiet)

    def _build(self, quiet):
        if not quiet:
            sys.stderr.write(
                "[trimview] indexing %s (one time, cached)... "
                % os.path.basename(self.path))
            sys.stderr.flush()
        lines, close = open_fastq(self.path)
        n = 0
        try:
            for i, line in enumerate(lines):
                if i & 3:
                    continue
                self.arr.append(_h64(line[1:].rstrip("\n")))
                n += 1
        finally:
            close()
        # array.sort() does not exist; go through a list once
        self.arr = array("Q", sorted(self.arr))
        try:
            os.makedirs(CACHE_DIR, exist_ok=True)
            tmp = self.cache + ".tmp%d" % os.getpid()
            with open(tmp, "wb") as fh:
                fh.write(self.arr.tobytes())
            os.replace(tmp, self.cache)
        except OSError:
            pass
        if not quiet:
            sys.stderr.write("%d reads\n" % n)

    def __contains__(self, name):
        h = _h64(name)
        i = bisect.bisect_left(self.arr, h)
        return i < len(self.arr) and self.arr[i] == h

    def __len__(self):
        return len(self.arr)


# ---------------------------------------------------------------------------
# alignment of a child read onto its parent
# ---------------------------------------------------------------------------


def find_offset(pseq, pqual, cseq, cqual):
    """Offset of the child read within its parent.

    Trimming only removes bases from the ends, so the child is always a
    substring.  Sequence alone is ambiguous on homopolymer reads (a poly-T
    child matches a poly-T parent at many offsets), so the quality string --
    effectively random -- is required to match at the same offset too.
    Returns (offset, ambiguous) or (None, False) if it is not a substring.
    """
    n = len(cseq)
    if n == 0:
        return (len(pseq), False)
    hits = []
    start = 0
    while True:
        i = pseq.find(cseq, start)
        if i < 0:
            break
        if pqual[i:i + n] == cqual:
            hits.append(i)
            if len(hits) > 1:
                break
        start = i + 1
    if hits:
        return (hits[0], len(hits) > 1)
    # quality can be rewritten in principle; fall back to sequence only
    i = pseq.find(cseq)
    return (i, True) if i >= 0 else (None, False)


def classify_cut(frag, end, marks=None):
    """Name the reason a fragment was removed.

    Prefers the feature track over the sequence: if the removed window is
    mostly covered by features that were located from the read name and the
    adapter (barcode, UMI, poly-A, adapter), say so -- that is knowledge, where
    matching the fragment against a list of adapters is inference.
    """
    if not frag:
        return ""
    if marks:
        covered = sum(1 for m in marks if m)
        if covered >= 0.5 * len(marks):
            names = []
            for ch, name in MARK_NAME:
                if ch in marks and name not in names:
                    names.append(name)
            if names:
                return " + ".join(names)
    if end == "3":
        for name, ad in ADAPTERS:
            if frag.startswith(ad[:min(len(ad), len(frag))]) and len(frag) >= 8:
                return name + " adapter"
    body = frag.upper()
    for b in "ACGT":
        if body.count(b) >= 0.9 * len(body) and len(body) >= 6:
            return "poly-%s" % b
    if len(frag) <= 3:
        return "quality"
    return ""


# ---------------------------------------------------------------------------
# where the barcode, UMI and adapter sit inside the read
# ---------------------------------------------------------------------------


def _paint(marks, start, end, ch, only_empty=False):
    """Paint a feature onto the track. `only_empty` leaves existing marks
    alone, which is what the adapter scans use: the barcode and the UMI are
    known exactly from the read name, an adapter hit is only a match."""
    for i in range(max(0, start), min(len(marks), end)):
        if not (only_empty and marks[i]):
            marks[i] = ch


def _mismatches(a, b, limit):
    """Hamming distance of two equal-length strings, or None past `limit`."""
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > limit:
                return None
    return n


def find_adapter(seq, pat, min_overlap=MIN_ADAPTER_HIT, err=ADAPTER_ERR):
    """Leftmost 3' adapter match: (pos, length, mismatches) or None.

    The two shapes cutadapt's -a allows: the whole adapter somewhere in the
    read, or a prefix of it flush with the 3' end because the read ran out
    first. Up to err*length mismatches -- without that, one sequencing error
    inside the adapter ends the match there and the rest of it gets drawn as
    if it were insert (read 21 of cell 002: 34 of 40 nt, the other 6 lost to a
    single miscall). Indels are not modelled; cutadapt allows them, but a
    block drawn shifted would be harder to read than one that stops early.
    """
    best = None
    for i in range(len(seq)):
        k = min(len(pat), len(seq) - i)
        if k < min_overlap:
            break
        mm = _mismatches(seq[i:i + k], pat[:k], int(k * err))
        if mm is None:
            continue
        score = k - 2 * mm
        if best is None or score > best[0]:
            best = (score, i, k, mm)
    return None if best is None else best[1:]


def find_adapter_5p(seq, pat, min_overlap=MIN_ADAPTER_HIT, err=ADAPTER_ERR):
    """Same, for a read in the reverse orientation: the construct arrives
    revcomp'd, so the adapter is rc(pat) -- either whole somewhere in the read,
    or its tail flush with the read's 5' end because the read starts partway
    through it. Returns (pos, length, mismatches) or None."""
    rcp = rc(pat)
    hit = find_adapter(seq, rcp, min_overlap, err)
    if hit is not None and hit[1] == len(rcp):
        return hit
    for j in range(min(len(rcp), len(seq)), min_overlap - 1, -1):
        mm = _mismatches(seq[:j], rcp[-j:], int(j * err))
        if mm is not None:
            return (0, j, mm)
    return None


def _paint_adapters(seq, marks, notes, adapter, anchor_end=None):
    """Everything technical at the 3' end: the configured read-through adapter,
    the fixed Illumina tail behind it, and any stock adapter left over."""
    hit = find_adapter(seq, adapter) if adapter else None
    rev = find_adapter_5p(seq, adapter) if adapter and hit is None else None
    if hit is not None:
        i, k, mm = hit
        _paint(marks, i, i + k, F_ADAPT, only_empty=True)
        notes.append("adapter at %d (%d of %d nt%s)"
                     % (i + 1, k, len(adapter),
                        ", %d mismatch" % mm if mm else ""))
        _paint_tail(seq, marks, notes, i + k)
    elif rev is not None:
        i, k, mm = rev
        _paint(marks, i, i + k, F_ADAPT, only_empty=True)
        notes.append("adapter (rc) at %d (%d of %d nt%s)"
                     % (i + 1, k, len(adapter),
                        ", %d mismatch" % mm if mm else ""))
    elif (adapter and anchor_end is not None
          and 0 < len(seq) - anchor_end < MIN_ADAPTER_HIT):
        # too few bases left to call an adapter on its own -- but they sit
        # immediately behind the anchor, where the adapter is what comes next,
        # so context makes them real where a bare 4-mer on its own would not be
        k = len(seq) - anchor_end
        if seq[anchor_end:] == adapter[:k]:
            _paint(marks, anchor_end, len(seq), F_ADAPT_PART, only_empty=True)
            notes.append("adapter starts %d nt before the read ends" % k)

    # a stock adapter that the construct above did not already explain: this is
    # what TrimGalore matched when pass 1 cut the read
    for name, ad in ADAPTERS:
        got = find_adapter(seq, ad)
        if got is None:
            continue
        i, k, mm = got
        if all(marks[p] for p in range(i, min(len(marks), i + k))):
            continue
        _paint(marks, i, i + k, F_ADAPT, only_empty=True)
        notes.append("%s adapter at %d (%d of %d nt)"
                     % (name, i + 1, k, len(ad)))


def _paint_tail(seq, marks, notes, start):
    """The rest of the P5 arm behind TRIM_ADAPTER3 (see ADAPTER_TAIL)."""
    blocks = []
    for name, pat in ADAPTER_TAIL:
        got = find_adapter(seq[start:], pat, min_overlap=min(6, len(pat)))
        if got is None:
            continue
        i, k, _ = got
        blocks.append((start + i, start + i + k, name))
    if not blocks:
        return
    blocks.sort()
    named = []
    prev_end = None
    for b0, b1, name in blocks:
        # the gap between two known parts is the sample index (10 nt here)
        if prev_end is not None and 0 < b0 - prev_end <= MAX_INDEX_GAP:
            _paint(marks, prev_end, b0, F_TECH, only_empty=True)
            named.append("i5 index %d nt" % (b0 - prev_end))
        _paint(marks, b0, b1, F_TECH, only_empty=True)
        named.append(name)
        prev_end = b1
    notes.append("Illumina tail from %d: %s" % (blocks[0][0] + 1,
                                                " + ".join(named)))


def annotate(seq, tags, adapter):
    """Feature track for one read: (marks list, notes list).

    The cell barcode and UMI are on the read name, so for this read the 12 nt
    that follow the poly-A are a known literal string -- that is exactly what
    trim_bc_anchor.py searches for, and this calls the same function.
    """
    marks = [None] * len(seq)
    notes = []
    cb, umi = tags.get("CB"), tags.get("RX")

    anchored = False
    anchor_end = None
    if cb and umi and HAVE_ANCHOR:
        anchor = rc(cb) + rc(umi)
        p = find_anchor(seq, anchor, 1, 4)
        if p >= 0:
            anchored = True
            obs = seq[p:p + len(anchor)]
            if obs == anchor:
                how = "exact"
            elif len(obs) < len(anchor):
                how = "partial, %d of %d nt (read ends inside it)" % (
                    len(obs), len(anchor))
            else:
                mm = sum(1 for a, b in zip(obs, anchor) if a != b)
                how = "%d mismatch" % mm
            _paint(marks, p, p + len(cb), F_CBC)
            _paint(marks, p + len(cb), p + len(anchor), F_UMI)
            anchor_end = min(len(seq), p + len(anchor))
            notes.append("anchor rc(CBC)+rc(UMI) at %d (%s)" % (p + 1, how))
            # The poly-A in front of the anchor. Pass 0 leaves it alone (it used
            # to strip it; removed 2026-08-03, see trimtest/README.md round 13)
            # -- cutting at the anchor makes it a clean 3' suffix and pass 2's
            # cutadapt --poly-a removes it there. Shown here because it is worth
            # seeing, using cutadapt's algorithm so the picture matches what
            # actually happens downstream.
            head = seq[:p]
            kept = _polya_cut(head)
            if len(kept) != len(head):
                _paint(marks, len(kept), p, F_POLY)
                notes.append("poly-A %d nt before it (removed in pass 2)"
                             % (p - len(kept)))

    # A read in the opposite orientation carries the same construct revcomp'd,
    # so it STARTS with [UMI][CBC][poly-T] instead of ending with the anchor.
    # The read usually begins partway through it, hence the suffix scan.
    fwd_hit = False
    if not anchored and cb and umi:
        fwd = umi + cb
        for k in range(len(fwd), 7, -1):
            if seq.startswith(fwd[-k:]):
                umi_end = max(0, len(umi) - (len(fwd) - k))
                _paint(marks, 0, umi_end, F_UMI)
                _paint(marks, umi_end, k, F_CBC)
                notes.append("5' [UMI][CBC], %d of %d nt -- read is in the "
                             "reverse orientation" % (k, len(fwd)))
                fwd_hit = True
                break

    if not anchored and not fwd_hit and cb:
        # no full anchor: still show a lone barcode hit, which is usually why
        # pass 0 did not fire (UMI miscalled, or a chance 6-mer)
        for pat, ch, what in ((rc(cb), F_CBC_ALONE, "rc(CBC)"),
                              (cb, F_CBC_ALONE, "CBC as-is")):
            i = seq.find(pat)
            if i >= 0:
                _paint(marks, i, i + len(pat), ch)
                notes.append("%s at %d, no anchor" % (what, i + 1))
                break
        if umi and not fwd_hit:
            for pat, what in ((rc(umi), "rc(UMI)"), (umi, "UMI as-is")):
                i = seq.find(pat)
                if i >= 0:
                    _paint(marks, i, i + len(pat), F_UMI_ALONE)
                    notes.append("%s at %d" % (what, i + 1))
                    break

    _paint_adapters(seq, marks, notes, adapter, anchor_end)
    return marks, notes


# ---------------------------------------------------------------------------
# step5 gene assignment, looked up on demand
# ---------------------------------------------------------------------------


def lookup_genes(directory, lib, cell, names):
    """name -> [(single|multi, gene, cigar), ...] for a handful of read names.

    Scans the two *_genes.bed.gz once per page with grep -F rather than holding
    an index in memory: a page is a few reads, and these BEDs run to hundreds of
    MB on a deep cell.  Cost is one pass per page, so it is opt-in (--genes).
    """
    out = dict((n, []) for n in names)
    if not names:
        return out
    prefix = os.path.join(directory, "%s_%s" % (lib, cell))
    pattern = "\n".join(names) + "\n"
    for kind, suffix in BED_STAGES:
        path = prefix + suffix
        if not os.path.exists(path):
            continue
        dec = shutil.which("pigz") or shutil.which("zcat")
        # grep -f wants its patterns in a file and the data on stdin
        with tempfile.NamedTemporaryFile("w", suffix=".pat", delete=False) as pf:
            pf.write(pattern)
            patfile = pf.name
        try:
            cmd = "%s %s %s | grep -F -f %s" % (
                dec, "-dc" if dec.endswith("pigz") else "", path, patfile)
            res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.DEVNULL)
            for line in res.stdout.decode().splitlines():
                f = line.split("\t")
                if len(f) < 7:
                    continue
                if f[3] in out:
                    out[f[3]].append((kind, f[5], f[6]))
        finally:
            os.unlink(patfile)
    return out


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------


def render_read(idx, base_rec, rows, pal, width, show_qual, genes=None,
                adapter=DEFAULT_ADAPTER3):
    """rows: list of (label, seq, qual, offset, present, ambiguous)."""
    name, oseq, oqual = base_rec
    L = len(oseq)
    out = []

    rid, _, tags = name.partition(";")
    tagd = dict(t.split(":", 1) for t in tags.split(";") if ":" in t)
    out.append(pal.head("── read #%d  %s" % (idx, rid)) + "   " +
               pal.dim("CB=%s UMI=%s cell=%s" %
                       (tagd.get("CB", "?"), tagd.get("RX", "?"),
                        tagd.get("SM", "?"))))

    marks, notes = annotate(oseq, tagd, adapter)
    if notes:
        out.append("   " + pal.label("found: ") + pal.dim("; ".join(notes)))

    # ---- per-stage summary table -----------------------------------------
    lab_w = max(len(r[0]) for r in rows) + 1
    prev = None
    alive = True
    for (label, seq, qual, off, present, amb) in rows:
        if not present:
            if alive:
                why = STAGE_DOC.get(label, (None, None))[1]
                out.append("   %-*s %s" %
                           (lab_w, label,
                            pal.bad("DROPPED") +
                            ("  (%s)" % why if why else "")))
                alive = False
            else:
                out.append("   %-*s %s" % (lab_w, label, pal.dim("-")))
            prev = None
            continue
        if prev is None:
            out.append("   %-*s %s" % (lab_w, label,
                                       pal.bold("%3d nt" % len(seq))))
        else:
            plabel, pseq, poff = prev
            if off is None:
                note = pal.warn("not a substring of %s" % plabel)
                out.append("   %-*s %s   %s" %
                           (lab_w, label, pal.bold("%3d nt" % len(seq)), note))
            else:
                c5 = off - poff
                c3 = (poff + len(pseq)) - (off + len(seq))
                bits = []
                if c5:
                    frag = oseq[poff:poff + c5]
                    why = classify_cut(frag, "5", marks[poff:poff + c5])
                    bits.append(pal.warn("-%d 5'" % c5) +
                                (" (%s)" % why if why else ""))
                if c3:
                    frag = oseq[off + len(seq):poff + len(pseq)]
                    why = classify_cut(frag, "3",
                                       marks[off + len(seq):poff + len(pseq)])
                    bits.append(pal.warn("-%d 3'" % c3) +
                                (" (%s)" % why if why else ""))
                if not bits:
                    bits.append(pal.good("unchanged"))
                if amb:
                    bits.append(pal.dim("[ambiguous placement]"))
                out.append("   %-*s %s   %s" %
                           (lab_w, label, pal.bold("%3d nt" % len(seq)),
                            "  ".join(bits)))
        prev = (label, seq, off)

    # ---- step5 gene assignment (only when --genes) -------------------------
    if genes is not None:
        hits = genes.get(name, [])
        if not hits:
            out.append("   %-*s %s" % (lab_w, "genes",
                                       pal.dim("no annotation hit")))
        else:
            kinds = sorted(set(h[0] for h in hits))
            seen = []
            for kind, gene, cig in hits:
                if gene not in seen:
                    seen.append(gene)
            head = "%s, %d %s" % ("+".join(kinds), len(seen),
                                  "locus" if len(seen) == 1 else "loci")
            out.append("   %-*s %s" % (lab_w, "genes", pal.label(head)))
            for gene in seen[:6]:
                cigs = sorted(set(h[2] for h in hits if h[1] == gene))
                out.append("   %-*s   %s %s" %
                           (lab_w, "", gene, pal.dim(" ".join(cigs[:2]))))
            if len(seen) > 6:
                out.append("   %-*s   %s" %
                           (lab_w, "", pal.dim("... %d more" % (len(seen) - 6))))

    # ---- alignment blocks -------------------------------------------------
    gut = lab_w + 4
    avail = max(20, width - gut)
    for bstart in range(0, L, avail):
        bend = min(L, bstart + avail)
        # ruler
        ticks = [" "] * (bend - bstart)
        nums = [" "] * (bend - bstart)
        for p in range(bstart, bend):
            if (p + 1) % 10 == 0:
                ticks[p - bstart] = "|"
                s = str(p + 1)
                st = p - bstart - len(s) + 1
                if st >= 0:
                    nums[st:st + len(s)] = list(s)
        out.append(" " * gut + pal.dim("".join(nums)))
        out.append(" " * gut + pal.dim("".join(ticks)))
        if any(marks):
            track = "".join(marks[p] or " " for p in range(bstart, bend))
            out.append("   %-*s %s" % (lab_w, "", pal.label(track)))

        for (label, seq, qual, off, present, amb) in rows:
            if not present:
                out.append("   %-*s %s" % (lab_w, label,
                                           pal.dim("-" * (bend - bstart))))
                continue
            if off is None:
                off = 0
            line = []
            for p in range(bstart, bend):
                style = MARK_STYLE.get(marks[p] if p < len(marks) else None, "")
                if off <= p < off + len(seq):
                    line.append(pal.base(seq[p - off], style))
                else:
                    line.append(pal.dim_base(oseq[p], style) if p < L else " ")
            out.append("   %-*s %s" % (lab_w, label, "".join(line)))
            if show_qual:
                q = []
                for p in range(bstart, bend):
                    if off <= p < off + len(seq):
                        q.append(qual[p - off])
                    else:
                        q.append(" ")
                out.append("   %-*s %s" % (lab_w, "", pal.dim("".join(q))))
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# the viewer
# ---------------------------------------------------------------------------


class Viewer(object):
    def __init__(self, directory, lib, cell, quiet=False):
        self.dir = directory
        self.lib = lib
        self.cell = cell
        prefix = os.path.join(directory, "%s_%s" % (lib, cell))
        self.streams = []
        self.index = None
        self.labels = []
        missing = []
        for label, suffix, kind in STAGES:
            path = prefix + suffix
            if not os.path.exists(path):
                missing.append(label)
                continue
            self.labels.append(label)
            if kind == "stream":
                self.streams.append((label, FastqStream(path)))
            else:
                self.index = (label, NameIndex(path, quiet=quiet))
        if not self.streams:
            raise SystemExit("no FASTQ stages found for %s_%s in %s"
                             % (lib, cell, directory))
        self.missing = missing
        self.n = 0

    def close(self):
        for _, s in self.streams:
            s.close()

    def next_read(self):
        """Advance one read of the base stage; return (base_rec, rows) or None."""
        base_label, base = self.streams[0]
        if base.rec is None:
            return None
        rec = base.take()
        name, oseq, oqual = rec
        self.n += 1
        rows = [(base_label, oseq, oqual, 0, True, False)]
        prev = (oseq, oqual, 0)
        for label, st in self.streams[1:]:
            if st.rec is not None and st.rec[0] == name:
                _, seq, qual = st.take()
                if prev is None:
                    rows.append((label, seq, qual, None, True, False))
                    prev = (seq, qual, 0)
                else:
                    pseq, pqual, poff = prev
                    rel, amb = find_offset(pseq, pqual, seq, qual)
                    off = None if rel is None else poff + rel
                    rows.append((label, seq, qual, off, True, amb))
                    prev = (seq, qual, off if off is not None else 0)
            else:
                rows.append((label, "", "", None, False, False))
                prev = None
        if self.index is not None:
            label, idx = self.index
            present = name in idx
            if present and prev is not None:
                rows.append((label, prev[0], prev[1], prev[2], True, False))
            elif present:
                rows.append((label, "", "", None, True, False))
            else:
                rows.append((label, "", "", None, False, False))
        return rec, rows


def changed(rows):
    """True if any stage trimmed or dropped this read."""
    lens = [len(r[1]) for r in rows if r[4]]
    return (not all(r[4] for r in rows)) or (len(set(lens)) > 1)


def dropped(rows):
    return not all(r[4] for r in rows)


# ---------------------------------------------------------------------------
# stats mode
# ---------------------------------------------------------------------------


def run_stats(v, n, pal):
    labels = v.labels
    present = dict((l, 0) for l in labels)
    total_len = dict((l, 0) for l in labels)
    dropped_at = dict((l, 0) for l in labels)
    cut5 = dict((l, 0) for l in labels)
    cut3 = dict((l, 0) for l in labels)
    hist = {}
    seen = 0
    while seen < n:
        got = v.next_read()
        if got is None:
            break
        rec, rows = got
        seen += 1
        prev_alive = True
        prev = None
        for (label, seq, qual, off, ok, amb) in rows:
            if ok:
                present[label] += 1
                total_len[label] += len(seq)
                if prev is not None and off is not None:
                    poff, pseq = prev
                    cut5[label] += off - poff
                    cut3[label] += (poff + len(pseq)) - (off + len(seq))
                prev = (off if off is not None else 0, seq)
            else:
                if prev_alive:
                    dropped_at[label] += 1
                prev_alive = False
                prev = None
        # length distribution of reads that make it through every stage --
        # reads dropped on the way have no final length and are excluded
        if rows[-1][4]:
            L = len(rows[-1][1])
            hist[L // 10 * 10] = hist.get(L // 10 * 10, 0) + 1

    print(pal.head("stats over %d reads of cell %s" % (seen, v.cell)))
    print("")
    print("  %-10s %10s %8s %8s %9s %9s" %
          ("stage", "reads", "% of cbc", "mean len", "cut 5'/rd", "cut 3'/rd"))
    base = seen or 1
    for l in labels:
        p = present[l]
        print("  %-10s %10d %7.2f%% %8.1f %9.2f %9.2f" %
              (l, p, 100.0 * p / base,
               (total_len[l] / p) if p else 0.0,
               (cut5[l] / p) if p else 0.0,
               (cut3[l] / p) if p else 0.0))
    print("")
    print("  reads lost at each stage:")
    for l in labels[1:]:
        print("    %-10s %10d  (%.2f%%)" %
              (l, dropped_at[l], 100.0 * dropped_at[l] / base))
    print("")
    print("  length distribution of reads surviving every stage (10 nt bins):")
    for k in sorted(hist):
        bar = "#" * min(60, int(60.0 * hist[k] / max(hist.values())))
        print("    %3d-%3d %8d  %s" % (k, k + 9, hist[k], bar))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def detect_lib(directory):
    libs = set()
    for f in os.listdir(directory):
        m = re.match(r"(.+?)_(\d{3})_cbc\.fastq\.gz$", f)
        if m:
            libs.add(m.group(1))
    if len(libs) == 1:
        return libs.pop()
    if not libs:
        return None
    raise SystemExit("several libraries here (%s); pass --lib"
                     % ", ".join(sorted(libs)))


def list_cells(directory, lib):
    cells = []
    for f in sorted(os.listdir(directory)):
        m = re.match(re.escape(lib) + r"_(\d{3})_cbc\.fastq\.gz$", f)
        if m:
            cells.append(m.group(1))
    return cells


def main():
    ap = argparse.ArgumentParser(
        description="IGV-like per-read view of the trim/ribo stages",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage")[-1])
    ap.add_argument("cell", nargs="?", help="cell id, e.g. 001")
    ap.add_argument("--dir", default=".", help="directory with the per-cell FASTQs")
    ap.add_argument("--lib", help="library prefix (auto-detected)")
    ap.add_argument("-n", "--num", type=int, default=3, help="reads per page")
    ap.add_argument("--skip", type=int, default=0, help="skip this many reads first")
    ap.add_argument("--find", help="jump to first read whose name contains this")
    ap.add_argument("--dropped", action="store_true",
                    help="only show reads dropped at some stage")
    ap.add_argument("--trimmed-only", action="store_true",
                    help="only show reads that actually changed")
    ap.add_argument("--qual", action="store_true", help="show quality strings")
    ap.add_argument("--adapter", default=DEFAULT_ADAPTER3,
                    help="3' read-through adapter to mark "
                         "(default = config.sh's TRIM_ADAPTER3); "
                         "pass '' to mark only the stock Illumina adapters")
    ap.add_argument("--genes", action="store_true",
                    help="also show the step5 gene assignment for each read "
                         "(one scan of the *_genes.bed.gz per page)")
    ap.add_argument("--width", type=int, default=0, help="columns (0 = terminal)")
    ap.add_argument("--no-color", action="store_true")
    ap.add_argument("--stats", nargs="?", type=int, const=100000, default=None,
                    metavar="N", help="print summary over N reads instead")
    ap.add_argument("--no-ribo", action="store_true",
                    help="skip the nonRibo stage (avoids building its index)")
    ap.add_argument("--list", action="store_true", help="list available cells")
    ap.add_argument("--batch", action="store_true",
                    help="print one page and exit (no prompt)")
    args = ap.parse_args()

    directory = os.path.abspath(args.dir)
    lib = args.lib or detect_lib(directory)
    if lib is None:
        raise SystemExit("no *_NNN_cbc.fastq.gz in %s -- wrong --dir?" % directory)

    pal = Palette(not args.no_color and sys.stdout.isatty())
    width = args.width or shutil.get_terminal_size((160, 40)).columns

    cells = list_cells(directory, lib)
    if args.list or not args.cell:
        print("library %s in %s" % (lib, directory))
        print("cells: %s" % " ".join(cells))
        if not args.cell:
            print("\npick one, e.g.:  %s %s" %
                  (os.path.basename(sys.argv[0]), cells[0] if cells else "001"))
        return

    global STAGES
    if args.no_ribo:
        STAGES = [s for s in STAGES if s[2] != "indexed"]

    cell = args.cell.zfill(3)
    v = Viewer(directory, lib, cell)

    if args.stats is not None:
        run_stats(v, args.stats, pal)
        v.close()
        return

    if v.missing:
        print(pal.warn("note: missing stage(s) %s" % ", ".join(v.missing)))

    state = {"num": args.num, "dropped": args.dropped,
             "trimmed": args.trimmed_only}

    def skip(k):
        for _ in range(k):
            if v.next_read() is None:
                return False
        return True

    def page(k, find=None):
        batch = []
        while len(batch) < k:
            got = v.next_read()
            if got is None:
                break
            rec, rows = got
            if find and find not in rec[0]:
                continue
            if state["dropped"] and not dropped(rows):
                continue
            if state["trimmed"] and not changed(rows):
                continue
            batch.append((v.n, rec, rows))
            find = None
        genes = None
        if args.genes and batch:
            genes = lookup_genes(directory, lib, v.cell,
                                 [b[1][0] for b in batch])
        for n, rec, rows in batch:
            print(render_read(n, rec, rows, pal, width, args.qual, genes,
                              args.adapter))
        if len(batch) < k:
            print(pal.dim("-- end of file --"))
            return False
        return True

    if args.skip:
        print(pal.dim("skipping %d reads..." % args.skip))
        skip(args.skip)

    print(pal.head("cell %s of %s   (%s)" % (cell, lib, directory)))
    for l in v.labels:
        print(pal.dim("  %-9s %s" % (l, STAGE_DOC.get(l, ("", ""))[0])))
    print(pal.dim("  bases still present are coloured; "
                  "trimmed-away bases stay in place, dim lower-case"))
    print(pal.dim("  feature track: %s cell barcode (reverse video)  "
                  "%s UMI (underlined)  %s poly-A before the anchor"
                  % (F_CBC, F_UMI, F_POLY)))
    print(pal.dim("  %s adapter  %s adapter cut short by the read end  "
                  "%s Illumina tail behind it (index, P5) "
                  "-- all struck through in the sequence"
                  % (F_ADAPT, F_ADAPT_PART, F_TECH)))
    print(pal.dim("  lower-case %s/%s = a lone barcode/UMI hit with no anchor "
                  "-- 6 nt can also match by chance" % (F_CBC_ALONE,
                                                        F_UMI_ALONE)))
    if not HAVE_ANCHOR:
        print(pal.warn("  trim_bc_anchor.py not importable -- "
                       "anchor track disabled"))
    print("")

    page(state["num"], find=args.find)
    if args.batch or not sys.stdin.isatty():
        v.close()
        return

    HELP = ("  Enter / n   next page        n5   next 5 reads\n"
            "  s1000       skip 1000        f<txt> find read name containing txt\n"
            "  d           dropped-only toggle    t  changed-only toggle\n"
            "  c003        switch to cell 003     q  quit")
    while True:
        try:
            cmd = input(pal.dim("[read %d] > " % v.n)).strip()
        except (EOFError, KeyboardInterrupt):
            print("")
            break
        if cmd in ("q", "quit", "exit"):
            break
        if cmd in ("?", "h", "help"):
            print(HELP)
            continue
        if cmd == "" or cmd == "n":
            page(state["num"])
            continue
        if cmd == "d":
            state["dropped"] = not state["dropped"]
            print(pal.dim("dropped-only: %s" % state["dropped"]))
            continue
        if cmd == "t":
            state["trimmed"] = not state["trimmed"]
            print(pal.dim("changed-only: %s" % state["trimmed"]))
            continue
        m = re.match(r"^n\s*(\d+)$", cmd)
        if m:
            state["num"] = int(m.group(1))
            page(state["num"])
            continue
        m = re.match(r"^s\s*(\d+)$", cmd)
        if m:
            skip(int(m.group(1)))
            print(pal.dim("at read %d" % v.n))
            continue
        m = re.match(r"^f\s*(\S+)$", cmd)
        if m:
            print(pal.dim("scanning forward..."))
            if not page(1, find=m.group(1)):
                print(pal.warn("'%s' not found -- f only scans forward from "
                               "here; use c%s to restart this cell"
                               % (m.group(1), v.cell)))
            continue
        m = re.match(r"^c\s*(\d+)$", cmd)
        if m:
            newcell = m.group(1).zfill(3)
            v.close()
            v = Viewer(directory, lib, newcell)
            print(pal.head("cell %s" % newcell))
            page(state["num"])
            continue
        print(pal.dim("unknown command; ? for help"))
    v.close()


if __name__ == "__main__":
    main()
