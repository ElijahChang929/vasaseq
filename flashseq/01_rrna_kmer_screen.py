#!/usr/bin/env python3
"""Measure the real rRNA content of each FLASH-seq library.

Why this script exists
----------------------
nf-core/rnaseq reports rRNA as the fraction of reads over genes annotated
gene_biotype "rRNA" in the Ensembl GTF. On GRCm39 / Ensembl 116 that annotation
contains 354 rRNA genes which are essentially all 5S (n-R5s*), plus a single
Rn18s-rs5 relic -- there is NO Rn45s, Rn28s or Rn5-8s gene at all, because the
rDNA array is collapsed out of the primary assembly. The reported 0.7-1.4% is
therefore a measurement of 5S, not of rRNA. Verify it yourself:

    awk -F'\\t' '$3=="gene" && $9~/gene_biotype "rRNA"/' $FS_MOUSE_GTF \\
      | grep -o 'gene_name "[^"]*"' | sort | uniq -c

This is the same defect that was found and fixed on the VASA side; see the
"Reference provenance" section of the repository CLAUDE.md and
own_version/build_rrna_reference.sh. The fix there was to add the NCBI 47S
pre-rRNA unit BK000964.3:1-13403 to the rRNA reference, which is what
unique_rRNA_mouse.v2.fa carries.

What this measures
------------------
Exact 31-mer containment against unique_rRNA_mouse.v2.fa, split into two sets
that are reported separately:

  * the 47S unit (BK000964.3)  -- what Ensembl does not have
  * everything else            -- the Ensembl rRNA-biotype genes, i.e. 5S

The split IS the point: it shows how much of the signal the annotation-based
number structurally cannot see.

Being exact-match and strand-aware but mismatch-intolerant, and sampling every
FS_KMER_STRIDE bases rather than every position, this is a LOWER BOUND. A real
alignment finds more. That comparison has since been made: 05_rrna_bwa.sh runs
the VASA pipeline's own ribo-bwamem.sh (bwa aln + bwa mem) over the same reads,
and res/flashseq/rrna_bwa.tsv reports both side by side. Do not extend this
script to try to close the gap -- a second half-alignment method is worth less
than one number from the same script the VASA side used.

Which reads get screened
------------------------
Every FS_STRIDE-th read pair across the whole file, NOT the first N. Head
sampling is measurably biased -- see the FS_STRIDE comment in config.sh, where
the adapter rate reads 58.5 / 57.7 / 56.4% off the first 20k / 400k / 2M reads
of ZHA8833A1 against 55.1% for the whole library, because a fastq is ordered by
flowcell position. Stride sampling reproduces the whole-library figure.

The stride is deterministic, so this screens EXACTLY the reads 05_rrna_bwa.sh
aligns, as long as both run with the same FS_STRIDE. That is what makes the
k-mer-vs-bwa gap a difference of method rather than of input, and it is why
neither script has to depend on the other to get it.

Output
------
    res/flashseq/rrna_kmer.tsv

Runtime: ~4 min per library -- the whole file is decompressed, since reaching
the tail is the entire point. Use the sbatch wrapper, 01_rrna_kmer_screen.sbatch.
"""
from __future__ import annotations

import csv
import gzip
import os
import re
import sys
from pathlib import Path

ROOT = Path(os.environ.get("FS_ROOT", "/nemo/lab/turnerj/working/guangxin/vasaseq"))
CODE = Path(os.environ.get("FS_CODE", ROOT / "code/flashseq"))
OUT = Path(os.environ.get("FS_OUT", ROOT / "res/flashseq"))
FASTQ = Path(os.environ.get(
    "FS_FASTQ",
    "/nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/RN26038/"
    "20260325_LH00442_0237_B23GT7GLT3/fastq",
))
RRNA_FA = Path(os.environ.get(
    "FS_RRNA_FA",
    "/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa",
))
K = int(os.environ.get("FS_K", "31"))
# Two different strides, and confusing them would be easy:
#   STRIDE      -- bases between sampled k-mers WITHIN a read
#   READ_STRIDE -- reads between sampled reads WITHIN a file
STRIDE = int(os.environ.get("FS_KMER_STRIDE", "10"))
READ_STRIDE = int(os.environ.get("FS_STRIDE", "64"))

# The record name build_rrna_reference.sh gives the NCBI 47S unit. If this stops
# matching, the 47S/Ensembl split silently collapses into one bucket, so fail
# loudly rather than report a wrong split.
R47S_PREFIX = "mouse_rDNA_47S"

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def kmers(seq: str) -> set[str]:
    """Every K-mer of seq, both strands, skipping any containing N."""
    out = set()
    for s in (seq, revcomp(seq)):
        for i in range(len(s) - K + 1):
            km = s[i:i + K]
            if "N" not in km:
                out.add(km)
    return out


def load_reference() -> tuple[set[str], set[str]]:
    """Return (47S k-mers, Ensembl-only k-mers) with the 47S set taking priority.

    A k-mer shared between the 47S unit and an Ensembl record is attributed to
    the 47S set, so the two columns stay additive and the Ensembl column answers
    'what could the annotation-based method have seen that the 47S does not
    already explain'.
    """
    if not RRNA_FA.exists():
        sys.exit(f"FATAL: missing rRNA reference {RRNA_FA}")
    records: dict[str, list[str]] = {}
    name = None
    with open(RRNA_FA) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
                records[name] = []
            else:
                records[name].append(line.strip().upper())

    r47 = [n for n in records if n.startswith(R47S_PREFIX)]
    if not r47:
        sys.exit(f"FATAL: no record starting '{R47S_PREFIX}' in {RRNA_FA} -- this is "
                 f"unique_rRNA_mouse.fa (v1), not v2. v1 has no 47S unit and would "
                 f"reproduce exactly the bug this script exists to expose.")

    k47: set[str] = set()
    for n in r47:
        k47 |= kmers("".join(records[n]))
    kens: set[str] = set()
    for n, parts in records.items():
        if n in r47:
            continue
        kens |= kmers("".join(parts))
    kens -= k47
    return k47, kens


def hits(read: str, kset: set[str]) -> bool:
    """Does any sampled K-mer of read appear in kset?"""
    return any(read[i:i + K] in kset for i in range(0, len(read) - K + 1, STRIDE))


def libraries() -> list[tuple[str, Path]]:
    """(library, R1 path), ordered A1..A10 rather than lexically."""
    found = []
    for p in FASTQ.glob("ZHA8833A*_R1_001.fastq.gz"):
        m = re.match(r"(ZHA8833A(\d+))_S", p.name)
        if m:
            found.append((m.group(1), int(m.group(2)), p))
    if not found:
        sys.exit(f"FATAL: no R1 fastqs under {FASTQ}")
    found.sort(key=lambda t: t[1])
    return [(lib, p) for lib, _, p in found]


def main() -> None:
    k47, kens = load_reference()
    print(f"reference {RRNA_FA.name}: {len(k47):,} 47S k-mers, "
          f"{len(kens):,} Ensembl-only k-mers (K={K}, stride={STRIDE})", flush=True)

    rows = []
    for lib, r1 in libraries():
        n47 = nens = total = seen = 0
        with gzip.open(r1, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                seen += 1
                # every READ_STRIDE-th read, across the whole file. Deliberately
                # no `break`: stopping early is exactly the head-of-file bias
                # this replaced, and the tail of the flowcell is where the
                # sample stops being representative.
                if (seen - 1) % READ_STRIDE:
                    continue
                total += 1
                read = line.strip()
                if hits(read, k47):
                    n47 += 1
                elif hits(read, kens):
                    nens += 1
        rows.append({
            "library": lib,
            "reads_in_file": seen,
            "read_stride": READ_STRIDE,
            "reads_screened": total,
            "pct_rRNA_47S": 100 * n47 / total,
            "pct_rRNA_ensembl_only": 100 * nens / total,
            "pct_rRNA_total": 100 * (n47 + nens) / total,
        })
        print(f"{lib:<12} n={total:>8,}  47S={rows[-1]['pct_rRNA_47S']:6.2f}%"
              f"  ensembl_only={rows[-1]['pct_rRNA_ensembl_only']:5.2f}%"
              f"  total={rows[-1]['pct_rRNA_total']:6.2f}%", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "rrna_kmer.tsv"
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {dest}")

    frac47 = sum(r["pct_rRNA_47S"] for r in rows)
    fracall = sum(r["pct_rRNA_total"] for r in rows)
    print(f"across all libraries, {100 * frac47 / fracall:.1f}% of the rRNA signal "
          f"comes from the 47S unit Ensembl does not annotate.")


if __name__ == "__main__":
    main()
