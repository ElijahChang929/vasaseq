#!/usr/bin/env python3
"""What TrimGalore actually removed, per overrepresented sequence.

Why this is a separate step
---------------------------
Knowing that a library is 30% adapter read-through says nothing on its own --
what matters is whether the trimmer took it out. Two artefacts in this run look
similar in the raw FASTQ and behave completely differently:

  * Adapter read-through IS removed. TrimGalore auto-detected Nextera
    ('CTGTCTCTTATA', the transposase mosaic end) and trimmed 63% of reads.
  * Poly-G is NOT removed. TrimGalore's adapter pass has no concept of it, so
    it survives essentially untouched into STAR, where it becomes part of the
    'unmapped, too short' fraction.

Comparing the two FastQC passes nf-core already ran is enough to tell them
apart, and needs no FASTQ streaming -- this runs in seconds.

Two limits inherited from FastQC, neither fixable here
------------------------------------------------------
FastQC inspects only the first 50 bp of each read, and (its own documented
behaviour) builds the overrepresented-sequence table from a sample taken at the
START of the file rather than the whole of it. So these percentages are not on
the same footing as 02_contaminant_check.py's, which counts every
FS_STRIDE-th read across the entire file and is the number to quote for a rate.
What FastQC is good for is exactly what it is used for here: the SAME
measurement on both sides of trimming, so the before/after difference is
meaningful even where the absolute level is not.

The comparison is per read FILE, not per library, because poly-G is strongly
mate-specific here (it is much heavier on R2) and averaging the mates would
hide that.

Output
------
    res/flashseq/trim_effect.tsv

Usage:
    source code/flashseq/config.sh
    fs_python code/flashseq/04_trim_effect.py
"""
from __future__ import annotations

import csv
import os
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(os.environ.get("FS_ROOT", "/nemo/lab/turnerj/working/guangxin/vasaseq"))
RESULTS = Path(os.environ.get("FS_RESULTS", ROOT / "data/flashseq/results"))
CODE = Path(os.environ.get("FS_CODE", ROOT / "code/flashseq"))
OUT = Path(os.environ.get("FS_OUT", ROOT / "res/flashseq"))


def overrepresented(zip_path: Path) -> dict[str, float]:
    """{sequence: percent} from one FastQC zip."""
    stem = zip_path.name[: -len("_fastqc.zip")]
    out: dict[str, float] = {}
    with zipfile.ZipFile(zip_path) as zf:
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
                    out[f[0]] = float(f[2])
    return out


def key_of(zip_path: Path) -> str | None:
    """'ZHA8833A10_1' from either a raw or a trimmed FastQC zip name.

    Raw:     ZHA8833A10_raw_1_fastqc.zip
    Trimmed: ZHA8833A10_trimmed_1_val_1_fastqc.zip
    The mate number is the one attached to 'raw'/'trimmed', not the trailing
    '_val_N' that TrimGalore appends.
    """
    m = re.match(r"(ZHA8833A\d+)_(?:raw|trimmed)_([12])", zip_path.name)
    return f"{m.group(1)}_{m.group(2)}" if m else None


def classify(seq: str) -> str:
    """Coarse label, enough to group the rows meaningfully."""
    for base in "ACGT":
        if seq.count(base) >= 0.9 * len(seq):
            return f"poly{base}"
    rc = seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    for probe in ("CTCGTGGGCTCGGAGATGTGTATAAG", "CTGTCTCTTATACACATCT",
                  "AGATCGGAAGAGC"):
        if probe in seq or probe in rc:
            return "adapter"
    return "other"


def main() -> None:
    raw_dir, trim_dir = RESULTS / "fastqc/raw", RESULTS / "fastqc/trim"
    for d in (raw_dir, trim_dir):
        if not d.is_dir():
            sys.exit(f"FATAL: missing {d}")

    raw: dict[str, dict[str, float]] = {}
    trim: dict[str, dict[str, float]] = {}
    for d, store in ((raw_dir, raw), (trim_dir, trim)):
        for z in sorted(d.glob("*_fastqc.zip")):
            k = key_of(z)
            if k:
                store[k] = overrepresented(z)
    if not raw or not trim:
        sys.exit("FATAL: could not read both FastQC passes")
    print(f"read {len(raw)} raw and {len(trim)} trimmed FastQC reports")

    with open(CODE / "sample_metadata.tsv") as fh:
        meta = {r["library"]: r for r in csv.DictReader(
            [ln for ln in fh if not ln.startswith("#")], delimiter="\t")}

    rows = []
    for k in sorted(raw, key=lambda s: (int(re.search(r"A(\d+)", s).group(1)), s)):
        lib, mate = k.rsplit("_", 1)
        # Union of both passes: a sequence can appear only after trimming, when
        # removing adapter from other reads lifts it over the reporting floor.
        for seq in sorted(set(raw[k]) | set(trim.get(k, {}))):
            before = raw[k].get(seq, 0.0)
            after = trim.get(k, {}).get(seq, 0.0)
            rows.append({
                "library": lib,
                "mate": f"R{mate}",
                "input_amount": meta.get(lib, {}).get("input_amount", ""),
                "sequence": seq,
                "kind": classify(seq),
                "pct_raw": f"{before:.4f}",
                "pct_trimmed": f"{after:.4f}",
                "pct_removed": f"{before - after:.4f}",
                # Guard against dividing by a floor-grazing 'before'.
                "fraction_removed": f"{(before - after) / before:.3f}" if before > 0.01 else "",
            })

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "trim_effect.tsv"
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {dest} ({len(rows)} rows)\n")

    print("what survives trimming, by artefact class "
          "(mean over read files where FastQC reported it):")
    print(f"{'kind':<10} {'n':>3}  {'mean pct raw':>12}  {'mean pct trimmed':>16}  {'removed':>8}")
    for kind in ("adapter", "polyG", "other"):
        sel = [r for r in rows if r["kind"] == kind and float(r["pct_raw"]) > 0.01]
        if not sel:
            continue
        b = sum(float(r["pct_raw"]) for r in sel) / len(sel)
        a = sum(float(r["pct_trimmed"]) for r in sel) / len(sel)
        print(f"{kind:<10} {len(sel):>3}  {b:>12.3f}  {a:>16.3f}  "
              f"{100 * (b - a) / b:>7.1f}%")

    print("\nThe contrast is the point: adapter is removed, poly-G is not, and")
    print("'other' (the human CALB1 fragments) is untouched because it is real")
    print("sequence rather than an artefact.")


if __name__ == "__main__":
    main()
