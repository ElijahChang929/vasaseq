#!/usr/bin/env python3
"""Tabulate what 05_rrna_bwa.sh measured, and put it beside the VASA number.

This script does no detection of its own. 05_rrna_bwa.sh ran the VASA
pipeline's own rRNA stage over the FLASH-seq reads; this reads what it wrote.

Why the arithmetic is imported rather than rewritten
----------------------------------------------------
`parse_log` and `load_intervals` come from own_version/step3_report.py -- the
same functions that produced the VASA percentages this table is compared
against. Reimplementing them here would create a second method that merely
looks like the first, which is the exact failure this whole exercise exists to
correct on the reference side. If step3_report.py ever moves, this fails loudly
rather than falling back to a copy.

`parse_bam` is NOT imported: it needs pysam, which the Anaconda interpreter the
rest of code/flashseq/ runs on does not have. 05_rrna_bwa.sh therefore dumps
flag/ref/pos per ribosomal read with `samtools view` into a .riboloci.tsv, and
the binning below works off that. The binning rule -- first interval containing
the 1-based leftmost position -- is step3_report.parse_bam's, and the interval
file is literally the same rrna_intervals.tsv.

What the columns mean
---------------------
    reads_in            reads the arm handed to bwa
    ribo / ribo_pct     called ribosomal by bwa aln + bwa mem together
    ribo_pct_stranded_y the same BAM re-selected with stranded=y, i.e. what you
                        would get by copying VASA's flag across unexamined. For
                        an unstranded library this should be about half of
                        ribo_pct; that it is IS the evidence that `n` is the
                        right flag here, and it is measured, not assumed.
    aln_only/mem_only/both
                        which aligner caught each read. At 151 nt aln_only is
                        expected to be tiny -- the mirror image of VASA, where
                        reads are short and aln_only was ~49% of detection on
                        cell 001. Both aligners are still run: dropping one
                        because it looks redundant on THIS library is how the
                        upstream short-read leak happened.
                        `both` is "everything that is not aln-only or mem-only",
                        which is step3_report.py's convention and so includes
                        the handful of `mem_mem` groups -- supplementary records
                        from mem alone, ~1% of the column. Kept as-is rather
                        than corrected, because diverging from the VASA side's
                        arithmetic would defeat the point of the comparison.
    fwd_pct             share of ribosomal reads on the forward strand
    5ETS..3ETS, mito, other
                        composition, as % of that arm's ribosomal reads
    kmer_pct            01_rrna_kmer_screen.py's exact-31-mer lower bound

Output
------
    res/flashseq/rrna_bwa.tsv        one row per library x arm
    stdout                           the two tables worth reading

Usage:
    source code/flashseq/config.sh
    fs_python code/flashseq/05_rrna_bwa_report.py
"""
from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path

ROOT = Path(os.environ.get("FS_ROOT", "/nemo/lab/turnerj/working/guangxin/vasaseq"))
CODE = Path(os.environ.get("FS_CODE", ROOT / "code/flashseq"))
OUT = Path(os.environ.get("FS_OUT", ROOT / "res/flashseq"))
WORK = Path(os.environ.get("FS_RRNA_BWA_DIR", OUT / "rrna_bwa"))
OWN = Path(os.environ.get("FS_VASA_OWN", ROOT / "code/I_Gene_expression/own_version"))
INTERVALS = Path(os.environ.get(
    "FS_RRNA_INTERVALS",
    "/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/rrna_intervals.tsv",
))

# The VASA library this is compared against: data/PM26037/out/logs/
# step3_report.txt, job 50788552, the ALL rows of its tables 1 and 2. Trimmed
# reads, bwa aln + mem, stranded=y, unique_rRNA_mouse.v2.fa -- i.e. the same
# scripts and the same reference 05_rrna_bwa.sh just ran over FLASH-seq.
#
# Transcribed by hand, so it can drift if that run is ever redone. Re-read it
# with:  step3_report.py data/PM26037/out/cells
VASA = {
    "label": "ZHA9292A1 (VASA-seq, 16 cells)",
    "ribo_pct": 21.39,
    "real_cell_lo": 18.0, "real_cell_hi": 25.1,
    "composition": {"5ETS": 17.1, "18S": 5.9, "ITS1": 4.6, "5.8S": 1.9,
                    "ITS2": 7.6, "28S": 54.4, "3ETS": 0.9,
                    "mito": 2.7, "other": 4.9},
}

sys.path.insert(0, str(OWN))
try:
    from step3_report import load_intervals, parse_log
except ImportError as exc:  # pragma: no cover -- a moved file, not a bad input
    sys.exit(f"FATAL: cannot import step3_report from {OWN} ({exc}).\n"
             f"       This script deliberately shares that module's arithmetic "
             f"with the VASA side rather than copying it. Point FS_VASA_OWN at "
             f"own_version/ rather than reimplementing it here.")

ARMS = ("raw", "trimmed")


def libraries() -> list[str]:
    found = []
    for p in sorted(WORK.glob("ZHA8833A*")):
        m = re.fullmatch(r"ZHA8833A(\d+)", p.name)
        if m and p.is_dir():
            found.append((int(m.group(1)), p.name))
    if not found:
        sys.exit(f"FATAL: no library directories under {WORK} -- run "
                 f"05_rrna_bwa.sh (or sbatch 05_rrna_bwa.sbatch) first.")
    return [name for _, name in sorted(found)]


def read_metadata() -> dict[str, dict[str, str]]:
    path = CODE / "sample_metadata.tsv"
    rows: dict[str, dict[str, str]] = {}
    with open(path) as fh:
        reader = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in reader:
            rows[r["library"]] = r
    return rows


def read_kmer() -> tuple[dict[str, float], dict[str, int]]:
    """01's lower bound and the read count it screened, if it has been run.

    The read count is the point: stride sampling is deterministic, so if 01 and
    05 used the same FS_STRIDE they screened the SAME reads, and the k-mer/bwa
    difference is a pure difference of method. main() checks that rather than
    hedging about it -- an unverified 'these are comparable' is worth nothing.
    Absent is not an error; a stale file from a head-sampled run is caught by
    the count not matching.
    """
    path = OUT / "rrna_kmer.tsv"
    if not path.exists():
        return {}, {}
    pct, n = {}, {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pct[r["library"]] = float(r["pct_rRNA_total"])
            if "reads_screened" in r:
                n[r["library"]] = int(r["reads_screened"])
    return pct, n


ETS_BIN = 200  # nt, the window size the VASA poly-T artefact was caught in


def composition(loci: Path, intervals):
    """(subunit counts, mito, other, forward, total, 5'ETS peak-bin share).

    The last value exists because of a specific incident on the VASA side: a
    poly-T population that survived trimming aligned to a T-rich stretch of the
    5'ETS and inflated it to 31-62% of a blank cell's ribosomal calls. It was
    caught by binning 5'ETS hits by position -- 88.9% of them sat in one 200 nt
    window, against a busiest bin of 10.5% in a real cell. A high 5'ETS share is
    legitimate for a total-RNA protocol, so the share alone cannot tell the two
    apart; the concentration can. Anything much above ~30% here means look at
    the reads before believing the number.
    """
    sub = {name: 0 for name, _, _ in intervals}
    ets_end = next((en for name, _, en in intervals if name == "5ETS"), 0)
    ets_bins: dict[int, int] = {}
    mito = other = fwd = total = 0
    with open(loci) as fh:
        header = next(fh, None)
        if header is None or not header.startswith("flag"):
            sys.exit(f"FATAL: {loci} has no header -- rerun 05_rrna_bwa.sh")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            flag, ref, pos = int(f[0]), f[1], int(f[2])
            total += 1
            if not flag & 0x10:
                fwd += 1
            if "rDNA_47S" in ref:
                for name, st, en in intervals:
                    if st <= pos <= en:
                        sub[name] += 1
                        if name == "5ETS":
                            b = (pos // ETS_BIN) * ETS_BIN
                            ets_bins[b] = ets_bins.get(b, 0) + 1
                        break
            elif "Mt_rRNA" in ref:
                mito += 1
            else:
                other += 1
    n_ets = sum(ets_bins.values())
    peak = 100 * max(ets_bins.values()) / n_ets if n_ets else 0.0
    return sub, mito, other, fwd, total, peak


def read_len_stats(path: Path) -> tuple[float, int]:
    """(mean, median) from a two-column length histogram."""
    lens: list[tuple[int, int]] = []
    with open(path) as fh:
        for line in fh:
            a, b = line.split()
            lens.append((int(a), int(b)))
    n = sum(c for _, c in lens)
    if not n:
        return 0.0, 0
    mean = sum(l * c for l, c in lens) / n
    seen, median = 0, lens[-1][0]
    for l, c in sorted(lens):
        seen += c
        if seen >= n / 2:
            median = l
            break
    return mean, median


def main() -> None:
    intervals = load_intervals(INTERVALS) if INTERVALS.exists() else []
    if not intervals:
        sys.exit(f"FATAL: no intervals at {INTERVALS} -- the composition columns "
                 f"would silently be empty.")
    names = [n for n, _, _ in intervals]
    meta = read_metadata()
    kmer, kmer_n = read_kmer()

    rows = []
    for lib in libraries():
        d = WORK / lib
        for arm in ARMS:
            # Every per-arm file is checked up front. 05_rrna_bwa.sh writes them
            # in sequence, so a run interrupted mid-library leaves a directory
            # that is indistinguishable from a complete one until you look for
            # the file that is missing.
            need = {k: d / f"{lib}.{arm}.{k}" for k in
                    ("ribo-map.log", "stranded_y.ribo-map.log",
                     "riboloci.tsv", "readlen.tsv")}
            for path in need.values():
                if not path.exists():
                    sys.exit(f"FATAL: {path} missing -- 05_rrna_bwa.sh did not finish "
                             f"{lib}/{arm}. Check the job log before reading anything "
                             f"in {WORK}; partial output here looks exactly like "
                             f"complete output.")

            total, unmapped, rg = parse_log(need["ribo-map.log"])
            aln, mem = rg.get("aln", 0), rg.get("mem", 0)
            both = sum(v for k, v in rg.items() if k not in ("aln", "mem"))
            ribo = aln + mem + both

            ytot, yunmapped, yrg = parse_log(need["stranded_y.ribo-map.log"])
            yribo = sum(yrg.values())

            sub, mito, other, fwd, nloci, ets_peak = composition(
                need["riboloci.tsv"], intervals)
            mean_len, med_len = read_len_stats(need["readlen.tsv"])

            m = meta.get(lib, {})
            row = {
                "library": lib,
                "input_amount": m.get("input_amount", ""),
                "input_pg": m.get("input_pg", ""),
                "replicate": m.get("replicate", ""),
                "well": m.get("well", ""),
                "arm": arm,
                "reads_in": total,
                "read_len_mean": round(mean_len, 1),
                "read_len_median": med_len,
                "ribo": ribo,
                "ribo_pct": round(100 * ribo / total, 3) if total else 0.0,
                "ribo_pct_stranded_y": round(100 * yribo / ytot, 3) if ytot else 0.0,
                "kept": unmapped,
                "aln_only": aln,
                "mem_only": mem,
                "both": both,
                "fwd_pct": round(100 * fwd / nloci, 1) if nloci else 0.0,
                "kmer_pct": kmer.get(lib, ""),
            }
            denom = sum(sub.values()) + mito + other
            for n in names:
                row[f"pct_{n}"] = round(100 * sub[n] / denom, 1) if denom else 0.0
            row["pct_mito"] = round(100 * mito / denom, 1) if denom else 0.0
            row["pct_other"] = round(100 * other / denom, 1) if denom else 0.0
            row["ets5_peak_bin_pct"] = round(ets_peak, 1)
            # ...and the same thing as a share of ALL reads. Composition
            # percentages alone mislead when the totals differ 3x: a subunit can
            # be a larger slice of a much smaller cake. abs_* is what actually
            # compares against VASA subunit for subunit.
            for n in names:
                row[f"abs_{n}"] = round(row["ribo_pct"] * row[f"pct_{n}"] / 100, 3)
            rows.append(row)

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "rrna_bwa.tsv"
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ---- table 1: the headline ---------------------------------------------
    print("TABLE 1 -- rRNA by bwa aln + bwa mem, unique_rRNA_mouse.v2.fa\n"
          "           (the VASA pipeline's own ribo stage, run over FLASH-seq reads)\n")
    hdr = (f"{'library':<11} {'input':>7} {'well':>5} {'arm':>8} {'reads':>9} "
           f"{'len':>4} {'ribo%':>7} {'y-flag%':>8} {'kmer%':>7} "
           f"{'aln':>6} {'mem':>7} {'both':>7} {'fwd%':>6}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        km = f"{r['kmer_pct']:6.2f}%" if r["kmer_pct"] != "" else "      -"
        print(f"{r['library']:<11} {r['input_amount']:>7} {r['well']:>5} {r['arm']:>8} "
              f"{r['reads_in']:>9,} {r['read_len_median']:>4} {r['ribo_pct']:>6.2f}% "
              f"{r['ribo_pct_stranded_y']:>7.2f}% {km} "
              f"{r['aln_only']:>6,} {r['mem_only']:>7,} {r['both']:>7,} "
              f"{r['fwd_pct']:>5.1f}%")

    # ---- table 2: composition ----------------------------------------------
    print("\nTABLE 2 -- where the ribosomal reads landed, trimmed arm "
          "(% of that library's ribosomal reads)\n"
          "           ets5pk = share of the 5'ETS hits in their busiest 200 nt "
          "window; see composition()\n")
    cols = names + ["mito", "other"]
    hdr2 = (f"{'library':<11} {'input':>7} " + " ".join(f"{c:>7}" for c in cols)
            + f" {'ets5pk':>7}")
    print(hdr2)
    print("-" * len(hdr2))
    tr = [r for r in rows if r["arm"] == "trimmed"]
    for r in tr:
        print(f"{r['library']:<11} {r['input_amount']:>7} "
              + " ".join(f"{r[f'pct_{c}']:6.1f}%" for c in cols)
              + f" {r['ets5_peak_bin_pct']:6.1f}%")
    print("-" * len(hdr2))
    print(f"{'VASA ALL':<11} {'':>7} "
          + " ".join(f"{VASA['composition'][c]:6.1f}%" for c in cols))

    # ---- table 3: the same thing as a share of all reads --------------------
    # This is the table that actually compares. A subunit's share OF THE
    # RIBOSOMAL READS can rise simply because there are fewer of them.
    print("\nTABLE 3 -- rRNA subunits as % of ALL reads, trimmed arm "
          "(ribo% x composition)\n")
    hdr3 = f"{'library':<11} {'input':>7} {'ribo%':>7} " + " ".join(f"{c:>7}" for c in names)
    print(hdr3)
    print("-" * len(hdr3))
    for r in tr:
        print(f"{r['library']:<11} {r['input_amount']:>7} {r['ribo_pct']:6.2f}% "
              + " ".join(f"{r[f'abs_{c}']:6.2f}%" for c in names))
    print("-" * len(hdr3))
    print(f"{'VASA ALL':<11} {'':>7} {VASA['ribo_pct']:6.2f}% "
          + " ".join(f"{VASA['ribo_pct'] * VASA['composition'][c] / 100:6.2f}%"
                     for c in names))

    # ---- the comparison this whole script exists for ------------------------
    lo = min(r["ribo_pct"] for r in tr)
    hi = max(r["ribo_pct"] for r in tr)
    raw = [r for r in rows if r["arm"] == "raw"]
    print(f"\nFLASH-seq, trimmed arm: {lo:.2f}-{hi:.2f}% ribosomal "
          f"across {len(tr)} libraries")
    print(f"VASA-seq {VASA['label']}: {VASA['ribo_pct']:.2f}% whole-library, "
          f"{VASA['real_cell_lo']:.1f}-{VASA['real_cell_hi']:.1f}% across the 12 real cells")
    print(f"  -> VASA carries {VASA['ribo_pct'] / ((lo + hi) / 2):.1f}x the rRNA, "
          f"same aligners, same reference, same script.")
    # and where that difference actually sits
    f28 = sum(r["abs_28S"] for r in tr) / len(tr)
    v28 = VASA["ribo_pct"] * VASA["composition"]["28S"] / 100
    fets = sum(r["abs_5ETS"] for r in tr) / len(tr)
    vets = VASA["ribo_pct"] * VASA["composition"]["5ETS"] / 100
    print(f"  -> mature 28S: {v28:.2f}% of VASA reads vs {f28:.2f}% of FLASH-seq "
          f"({v28 / f28:.1f}x). 5'ETS: {vets:.2f}% vs {fets:.2f}% "
          f"({vets / fets:.1f}x). The gap is in mature rRNA, not pre-rRNA.")
    d_trim = sum(r["ribo_pct"] for r in tr) / len(tr) - \
        sum(r["ribo_pct"] for r in raw) / len(raw)
    print(f"\ntrimming moved the FLASH-seq figure by {d_trim:+.2f} points on average "
          f"(raw -> trimmed)")
    if kmer:
        gaps = [r["ribo_pct"] - r["kmer_pct"] for r in raw if r["kmer_pct"] != ""]
        # Same reads or not? Stride sampling is deterministic, so equal screened
        # counts mean 01 and 05 saw the identical reads and the gap below is a
        # difference of method alone. Unequal means one of them sampled
        # differently -- most likely a stale rrna_kmer.tsv from the head-sampled
        # era -- and the comparison is not clean.
        same = all(kmer_n.get(r["library"]) == r["reads_in"] for r in raw) and bool(kmer_n)
        if gaps:
            print(f"\nbwa finds {sum(gaps) / len(gaps):+.2f} points more than the k-mer "
                  f"lower bound on the raw arm (n={len(gaps)}), i.e. the exact-match "
                  f"screen was already within "
                  f"{100 * sum(gaps) / sum(r['ribo_pct'] for r in raw):.1f}% of it.")
            if same:
                print("  Both screened the SAME reads (identical stride, identical "
                      "counts), so that is a difference of method alone.")
            else:
                print("  WARNING: rrna_kmer.tsv screened different read counts from "
                      "this run, so the two are not on the same reads. Re-run "
                      "01_rrna_kmer_screen.sbatch with the same FS_STRIDE before "
                      "quoting the gap.")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
