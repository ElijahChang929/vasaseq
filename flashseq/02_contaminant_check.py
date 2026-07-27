#!/usr/bin/env python3
"""Identify and quantify overrepresented sequences across all FLASH-seq libraries.

Why this script exists
----------------------
FastQC flagged one sequence as 15.9% of library ZHA8833A8. FastQC itself says
only "No Hit" -- it has no reference to check against. Tracking it down by hand
established that it is human CALB1, not mouse, and that it is a WELL effect
(G:1, H:1) rather than a low-input effect. This script does that work for every
library so the next run does not need the manual chase.

Two stages:

  1. CLASSIFY. Take every overrepresented sequence FastQC reported, in any
     library, and locate it by exact streaming search (both strands) against
     mouse GRCm39, the rRNA reference, ERCC92, and human GRCh38. A hit is then
     resolved to a gene via the matching GTF. Poly-base runs are labelled
     without searching -- a 50-mer of G is a sequencer artefact, and searching
     for it in a genome is meaningless.

  2. RECOUNT. Count every classified sequence in every library. This is the step
     that matters: FastQC only reports a sequence once it exceeds ~0.1% of a
     file, so a contaminant present at 0.04% is invisible in the per-library
     FastQC tables. Without this stage you cannot tell "absent" from "below the
     reporting floor", and that distinction is exactly what showed the CALB1
     contamination is confined to two adjacent wells.

Output
------
    res/flashseq/overrepresented.tsv
        one row per (sequence, library), with source / gene / locus and the
        measured percentage in that library

Which reads stage 2 counts
--------------------------
Every FS_STRIDE-th read of every file, across the whole file -- not the first N.
Head sampling is measurably biased (see config.sh), and while that bias never
came close to threatening the CALB1 result, which is a 25x difference between
libraries, the adapter read-through and poly-G rates this script also reports
are precisely the quantity it distorts.

Runtime: stage 1 is ~4 min per genome pass (mouse ~2.7 GB, human ~3.1 GB
gzipped). Stage 2 now decompresses all 20 fastqs end to end to reach their
tails, which is the bulk of the job. Submit it, do not run it on a login node:

    sbatch code/flashseq/02_contaminant_check.sbatch
"""
from __future__ import annotations

import csv
import gzip
import os
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(os.environ.get("FS_ROOT", "/nemo/lab/turnerj/working/guangxin/vasaseq"))
RESULTS = Path(os.environ.get("FS_RESULTS", ROOT / "data/flashseq/results"))
OUT = Path(os.environ.get("FS_OUT", ROOT / "res/flashseq"))
FASTQ = Path(os.environ.get(
    "FS_FASTQ",
    "/nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/RN26038/"
    "20260325_LH00442_0237_B23GT7GLT3/fastq",
))
# Reads between sampled reads. Every READ_STRIDE-th read of every file is
# screened, rather than the first N -- see config.sh's FS_STRIDE comment for the
# measurement that made that change necessary. 01_rrna_kmer_screen.py and
# 05_rrna_bwa.sh use the same variable, so with the default all three screen the
# same reads.
READ_STRIDE = int(os.environ.get("FS_STRIDE", "64"))

REFS = [
    # (label, fasta path, gtf path or None)
    ("mouse_GRCm39", os.environ.get("FS_MOUSE_FA"), os.environ.get("FS_MOUSE_GTF")),
    ("rRNA_47S_v2", os.environ.get("FS_RRNA_FA"), None),
    ("ERCC92", os.environ.get("FS_ERCC_FA"), None),
    ("human_GRCh38", os.environ.get("FS_HUMAN_FA"), os.environ.get("FS_HUMAN_GTF")),
]

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def opener(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# --- stage 0: gather the sequences FastQC flagged ----------------------------

def fastqc_overrepresented() -> dict[str, list[tuple[str, float]]]:
    """{sequence: [(library_readfile, pct), ...]} from the raw FastQC zips."""
    zips = sorted((RESULTS / "fastqc/raw").glob("*_fastqc.zip"))
    if not zips:
        sys.exit(f"FATAL: no FastQC zips under {RESULTS/'fastqc/raw'}")
    seqs: dict[str, list[tuple[str, float]]] = {}
    for z in zips:
        stem = z.name[:-len("_fastqc.zip")]
        with zipfile.ZipFile(z) as zf:
            with zf.open(f"{stem}_fastqc/fastqc_data.txt") as fh:
                inblock = False
                for raw in fh:
                    line = raw.decode()
                    if line.startswith(">>Overrepresented sequences"):
                        inblock = True
                        continue
                    if inblock and line.startswith(">>END_MODULE"):
                        break
                    if not inblock or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 3:
                        seqs.setdefault(f[0], []).append((stem, float(f[2])))
    print(f"FastQC flagged {len(seqs)} distinct sequences across {len(zips)} read files")
    return seqs


def is_homopolymer(seq: str) -> str | None:
    """Label obvious sequencer artefacts so they are not searched for.

    Poly-G on a NovaSeq X is not a sequence at all: the two-colour chemistry
    calls G when there is no signal. Searching a genome for it is meaningless.
    """
    for base in "ACGT":
        if seq.count(base) >= 0.9 * len(seq):
            return f"poly{base}"
    return None


# Library-construction oligos. Adapter read-through is in no genome, so without
# this table it would be reported as 'unidentified' and look like a mystery
# contaminant. Matched as a substring, either strand.
#
# Keep these CORES SHORT. A read that runs into the adapter usually shows only
# part of it, so a full-length oligo will not match: the flagged sequences here
# end at '...GGAGATGTGTATAAG' and the full 34 nt Nextera R2 oligo missed all
# eight of them. Every pattern below was checked against the actual flagged
# sequences before being committed.
ADAPTERS = {
    "P7_index_read_through": "CTCGTGGGCTCGGAGATGTGTATAAG",
    "P5": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "TruSeq_read_through": "AGATCGGAAGAGC",
    "Nextera_R1": "TCGTCGGCAGCGTC",
    "Nextera_ME": "AGATGTGTATAAGAGACAG",
    "polyA_tail": "AAAAAAAAAAAAAAAAAAAA",
}


def is_adapter(seq: str) -> str | None:
    """Name the library-construction oligo this sequence contains, if any."""
    rc = revcomp(seq)
    for name, oligo in ADAPTERS.items():
        if oligo in seq or oligo in rc:
            return f"adapter_{name}"
    return None


# Index-independent probes, counted in EVERY library whether or not FastQC
# flagged anything matching them.
#
# Necessary because the FastQC-derived sequences are not comparable across
# libraries: adapter read-through carries the library's own index, so each
# library produces a different 50-mer and FastQC only flags the ones that clear
# its threshold. Summing those per library gives 14.3% for A3 and exactly 0.00%
# for A1 -- and the 0.00% is a reporting-floor artefact, not an absence of
# adapter. These probes are the shared, index-free cores, so their rates mean
# the same thing in every library.
PROBES = {
    "probe_adapter_P7_readthrough": "CTCGTGGGCTCGGAGATGTGTATAAG",
    "probe_adapter_TruSeq": "AGATCGGAAGAGC",
    "probe_polyG_30": "G" * 30,
    "probe_polyA_30": "A" * 30,
}


def mismatches(a: str, b: str, cap: int) -> int:
    """Hamming distance between equal-length strings, giving up past cap."""
    if len(a) != len(b):
        return cap + 1
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > cap:
                return n
    return n


def nearest_classified(seq: str, classified: dict[str, dict], cap: int = 2):
    """Find an already-identified sequence within `cap` mismatches of seq.

    A high-abundance fragment produces sequencing-error variants of itself that
    also clear FastQC's reporting threshold -- this run has one, a single G->C
    substitution in the dominant fragment. Searching a genome for such a variant
    fails (exact matching), and calling it 'unidentified' would invent a second
    mystery sequence out of a base-call error. Inheriting the neighbour's label
    is both cheaper and closer to the truth.
    """
    for other, meta in classified.items():
        if other == seq or meta["source"] in ("unknown", "unidentified"):
            continue
        if mismatches(seq, other, cap) <= cap:
            return other, meta
    return None, None


# --- stage 1: locate each sequence in a reference ----------------------------

def scan_reference(label: str, fasta: Path, queries: dict[str, str]) -> dict[str, tuple]:
    """Exact streaming search for every query (and its revcomp) in one fasta.

    Streams in overlapping 4 Mb blocks so memory stays flat regardless of genome
    size; the 200 bp overlap is longer than any query, so a match spanning a
    block boundary is not missed. No index is needed or built.

    Returns {query_name: (contig, 1-based_start, strand)}.
    """
    pats = {}
    for name, q in queries.items():
        pats[(name, "+")] = q
        pats[(name, "-")] = revcomp(q)
    maxlen = max(len(p) for p in pats.values())
    overlap = maxlen + 10

    found: dict[str, tuple] = {}
    buf, contig, offset = "", None, 0

    def flush(buf: str, contig: str | None, offset: int) -> None:
        if not contig:
            return
        for (name, strand), pat in pats.items():
            if name in found:
                continue
            i = buf.find(pat)
            if i >= 0:
                found[name] = (contig, offset + i + 1, strand)

    with opener(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(buf, contig, offset)
                contig, buf, offset = line[1:].split()[0], "", 0
                continue
            buf += line.strip().upper()
            if len(buf) > 4_000_000:
                flush(buf, contig, offset)
                offset += len(buf) - overlap
                buf = buf[-overlap:]
        flush(buf, contig, offset)

    print(f"  {label}: located {len(found)}/{len(queries)}", flush=True)
    return found


def gene_at(gtf: Path, contig: str, pos: int) -> tuple[str, str]:
    """(gene_name, gene_id) of the gene spanning contig:pos, '' if intergenic."""
    if not gtf or not Path(gtf).exists():
        return "", ""
    with opener(Path(gtf)) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "gene" or f[0] != contig:
                continue
            if int(f[3]) <= pos <= int(f[4]):
                gid = re.search(r'gene_id "([^"]+)"', f[8])
                gnm = re.search(r'gene_name "([^"]+)"', f[8])
                return (gnm.group(1) if gnm else ""), (gid.group(1) if gid else "")
    return "", ""


# --- stage 2: count every classified sequence in every library ---------------

def libraries() -> list[tuple[str, Path, Path]]:
    """(library, R1, R2), ordered A1..A10 rather than lexically."""
    found = []
    for p in FASTQ.glob("ZHA8833A*_R1_001.fastq.gz"):
        m = re.match(r"(ZHA8833A(\d+))_S", p.name)
        if not m:
            continue
        r2 = Path(str(p).replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
        if not r2.exists():
            sys.exit(f"FATAL: {p.name} has no R2 mate")
        found.append((m.group(1), int(m.group(2)), p, r2))
    if not found:
        sys.exit(f"FATAL: no R1 fastqs under {FASTQ}")
    found.sort(key=lambda t: t[1])
    return [(lib, r1, r2) for lib, _, r1, r2 in found]


def count_in_libraries(seqs: list[str]) -> dict[str, dict[str, dict[str, float]]]:
    """{library: {sequence: {'R1': pct, 'R2': pct}}} by exact substring match.

    Two things this stage does that FastQC's own tables cannot:

      * It beats the ~0.1% reporting floor, so 'absent' is distinguishable from
        'present but below the floor'. That distinction is what showed the CALB1
        contamination is confined to two adjacent wells rather than being a
        general low-input effect.
      * It scans BOTH mates. A fragment's R1 and R2 sequences are different, so
        a sequence FastQC flagged in the _2 file occurs essentially never in R1
        -- counting it against R1 alone would report 0% and read as 'absent'
        when it is in fact abundant on the other mate.

    Matching is done in BOTH orientations. An adapter core is written in one
    strand's convention but the reads carry whichever orientation sequencing
    produced -- the P7 read-through core appears in these reads only as its
    reverse complement, so a forward-only search scored it 0.00% in all ten
    libraries when it is in fact up to 14%. `is_adapter` already checked both
    strands, so forward-only counting here was also internally inconsistent.

    Both-strand matching is safe for the mate-pair bookkeeping: R1 and R2 come
    from opposite ends of a ~1300 bp fragment (measured insert size), so they do
    not overlap and neither mate can contain the other's reverse complement.
    """
    # Precompute both orientations once rather than per read.
    both = {s: (s, revcomp(s)) for s in seqs}

    out: dict[str, dict[str, dict[str, float]]] = {}
    for lib, r1, r2 in libraries():
        per_mate: dict[str, dict[str, float]] = {}
        for mate, path in (("R1", r1), ("R2", r2)):
            counts = dict.fromkeys(seqs, 0)
            total = seen = 0
            with gzip.open(path, "rt") as fh:
                for i, line in enumerate(fh):
                    if i % 4 != 1:
                        continue
                    seen += 1
                    # Every READ_STRIDE-th read across the whole file. This used
                    # to take the first N and stop, which is biased: a fastq is
                    # ordered by flowcell position and adapter content drifts
                    # along it (config.sh has the measurement). The bias is only
                    # a few percent relative -- it never threatened the CALB1
                    # finding, which is a 25x effect -- but the adapter
                    # read-through and poly-G rates below ARE the quantity it
                    # distorts, so they are worth getting right.
                    if (seen - 1) % READ_STRIDE:
                        continue
                    total += 1
                    for s, (fwd, rev) in both.items():
                        if fwd in line or rev in line:
                            counts[s] += 1
            per_mate[mate] = {s: 100 * c / total for s, c in counts.items()}
        out[lib] = {s: {"R1": per_mate["R1"][s], "R2": per_mate["R2"][s]} for s in seqs}
        print(f"  {lib}: counted {len(seqs)} sequences in {total:,} reads x 2 mates",
              flush=True)
    return out


def main() -> None:
    flagged = fastqc_overrepresented()

    # Classify: homopolymers by inspection, the rest by search.
    classification: dict[str, dict] = {}
    to_search: dict[str, str] = {}
    for i, seq in enumerate(sorted(flagged)):
        name = f"seq{i:03d}"
        # Homopolymers and library oligos are named by inspection; searching a
        # genome for either is meaningless and would leave them 'unidentified'.
        label = is_homopolymer(seq) or is_adapter(seq)
        classification[seq] = {"name": name, "source": label or "unknown",
                               "gene": "", "locus": ""}
        if not label:
            to_search[name] = seq
    by_name = {v["name"]: k for k, v in classification.items()}
    print(f"{len(flagged) - len(to_search)} labelled by inspection "
          f"(homopolymer/adapter), {len(to_search)} to search against references")

    # Search references in order; first hit wins, so mouse beats human for any
    # sequence conserved between them, and a real contaminant is only called
    # human once mouse has been ruled out.
    remaining = dict(to_search)
    for label, fa, gtf in REFS:
        if not remaining:
            break
        if not fa or not Path(fa).exists():
            print(f"  {label}: SKIPPED, {fa} not readable")
            continue
        hits = scan_reference(label, Path(fa), remaining)
        for name, (contig, pos, strand) in hits.items():
            seq = by_name[name]
            gname, gid = gene_at(Path(gtf), contig, pos) if gtf else ("", "")
            classification[seq].update({
                "source": label,
                "gene": gname or gid,
                "locus": f"{contig}:{pos}({strand})",
            })
            remaining.pop(name, None)
    # Anything still unplaced: check whether it is a sequencing-error variant of
    # something already identified before declaring it a mystery.
    for name in list(remaining):
        seq = by_name[name]
        near, meta = nearest_classified(seq, classification)
        if meta:
            classification[seq].update({
                "source": f"{meta['source']}_variant",
                "gene": meta["gene"],
                "locus": meta["locus"],
            })
            print(f"  {name}: {mismatches(seq, near, 2)} mismatch(es) from a "
                  f"{meta['source']} sequence -- labelled as its variant")
        else:
            classification[seq]["source"] = "unidentified"

    # The index-independent probes join the count list as first-class entries,
    # so the output carries one comparable number per library per artefact type.
    for label, probe in PROBES.items():
        classification.setdefault(probe, {"name": label, "source": label,
                                          "gene": "", "locus": ""})

    # Recount all of them everywhere.
    print("\ncounting every classified sequence in every library:")
    counts = count_in_libraries(sorted(classification))

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "overrepresented.tsv"
    fields = ["library", "sequence", "pct_R1", "pct_R2", "pct_in_library",
              "source", "gene", "locus", "fastqc_flagged_in"]

    def peak_of(lib: str, seq: str) -> float:
        """A fragment shows up on one mate, not both -- the larger is its rate."""
        return max(counts[lib][seq]["R1"], counts[lib][seq]["R2"])

    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for lib in counts:
            for seq, meta in sorted(classification.items(),
                                    key=lambda kv: -peak_of(lib, kv[0])):
                w.writerow({
                    "library": lib,
                    "sequence": seq,
                    "pct_R1": f"{counts[lib][seq]['R1']:.4f}",
                    "pct_R2": f"{counts[lib][seq]['R2']:.4f}",
                    "pct_in_library": f"{peak_of(lib, seq):.4f}",
                    "source": meta["source"],
                    "gene": meta["gene"],
                    "locus": meta["locus"],
                    # probes were not flagged by FastQC; they are ours
                    "fastqc_flagged_in": ";".join(f for f, _ in flagged.get(seq, [])),
                })
    print(f"\nwrote {dest}")
    print("NOTE on adding these up: a fragment's R1 and R2 sequences are separate rows,")
    print("  so summing pct_in_library (= max of the two mates) DOUBLE-COUNTS every pair.")
    print("  Sum pct_R1 alone -- the R2-mate rows contribute ~0 to it, so the total is right.")

    print("\nsequences reaching >1% in any library:")
    for seq, meta in classification.items():
        peak = max(peak_of(l, seq) for l in counts)
        if peak > 1:
            worst = max(counts, key=lambda l: peak_of(l, seq))
            print(f"  {meta['source']:<22} {meta['gene']:<8} {meta['locus']:<20} "
                  f"peak {peak:5.2f}% in {worst}")


if __name__ == "__main__":
    main()
