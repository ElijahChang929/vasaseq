#!/usr/bin/env python3
"""
aligned_composition.py BAM [BAM ...]

For every UNIQUELY mapped read (NH==1), look only at the bases STAR actually
aligned (soft clips excluded) and report:

  A-rich    aligned block is >=80% A or >=80% T -- i.e. the "alignment" is
            really the poly-A tail (or its complement) stuck to a genomic
            A-tract. These are the false positives to worry about.
  short     aligned block < 30 nt
  clip      how much of the read STAR had to throw away to make the alignment

A trimming setting that inflates the read count by inflating A-rich is not
gaining anything; one that inflates it while A-rich stays flat is real.
"""
import sys, collections, pysam

def frac(s, b):
    return s.count(b) / len(s) if s else 0.0

print(f"{'bam':28s} {'uniq':>9s} {'A-rich':>9s} {'%':>6s} "
      f"{'short<30':>9s} {'%':>6s} {'medAln':>7s} {'medClip':>8s}")
for fn in sys.argv[1:]:
    n = arich = short = 0
    alen = []
    clip = []
    with pysam.AlignmentFile(fn, "rb") as bam:
        for r in bam:
            if r.is_unmapped or r.get_tag("NH") != 1:
                continue
            n += 1
            q = r.query_alignment_sequence or ""
            alen.append(len(q))
            clip.append(r.query_length - len(q))
            if len(q) < 30:
                short += 1
            if frac(q, "A") >= 0.8 or frac(q, "T") >= 0.8:
                arich += 1
    if not n:
        print(f"{fn:28s} {'-- empty --'}")
        continue
    alen.sort(); clip.sort()
    print(f"{fn:28s} {n:9d} {arich:9d} {100*arich/n:5.1f}% "
          f"{short:9d} {100*short/n:5.1f}% {alen[len(alen)//2]:7d} {clip[len(clip)//2]:8d}")
