#!/usr/bin/env python3
"""Collect the nf-core/rnaseq QC numbers into one per-library table.

MultiQC already computed everything here; the problem is that it is spread over
a dozen files, keyed by library id with no idea that the run is a titration, and
that two of the files need special handling (see the parsers below). This script
joins it all onto sample_metadata.tsv and writes one tidy table.

    res/flashseq/qc_summary.tsv   one row per library, input amount attached

Reads only. Runs in seconds; no need for sbatch.

Usage:
    source code/flashseq/config.sh && fs_load_python
    $FS_PY code/flashseq/00_collect_qc.py
"""
from __future__ import annotations

import ast
import csv
import os
import re
import sys
from pathlib import Path

# --- paths, all overridable from config.sh ----------------------------------
ROOT = Path(os.environ.get("FS_ROOT", "/nemo/lab/turnerj/working/guangxin/vasaseq"))
MQC = Path(os.environ.get("FS_MQC", ROOT / "data/flashseq/results/multiqc/star_rsem/multiqc_report_data"))
SCREEN = Path(os.environ.get(
    "FS_SCREEN",
    "/nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/RN26038/"
    "20260325_LH00442_0237_B23GT7GLT3/fastq_screen",
))
CODE = Path(os.environ.get("FS_CODE", ROOT / "code/flashseq"))
OUT = Path(os.environ.get("FS_OUT", ROOT / "res/flashseq"))


def read_metadata() -> list[dict]:
    """sample_metadata.tsv, minus its '#' header comment."""
    rows = []
    with open(CODE / "sample_metadata.tsv") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        rows.append(r)
    if len(rows) != 10:
        sys.exit(f"FATAL: expected 10 libraries in sample_metadata.tsv, got {len(rows)}")
    return rows


def read_mqc(name: str) -> dict[str, dict[str, str]]:
    """A MultiQC *.txt table, keyed by its first column."""
    path = MQC / name
    if not path.exists():
        sys.exit(f"FATAL: missing {path}")
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        key = rdr.fieldnames[0]
        return {r[key]: r for r in rdr}


def num(row: dict | None, field: str):
    """MultiQC leaves blanks for metrics a row does not carry."""
    if row is None:
        return ""
    v = row.get(field, "")
    if v in ("", None):
        return ""
    return float(v)


def idxstats_fractions(row: dict) -> dict[str, float]:
    """Per-contig read fractions from multiqc_samtools_idxstats.txt.

    Every cell in that file is a Python literal '[mapped_reads, contig_length]',
    not a number -- csv gives it back as the string "[3467045, 195154279]", and
    float() on it raises. ast.literal_eval is the parser that works.
    """
    counts = {}
    for contig, cell in row.items():
        if contig == "Sample" or not cell:
            continue
        try:
            counts[contig] = ast.literal_eval(cell)[0]
        except (ValueError, SyntaxError, IndexError, TypeError):
            continue
    total = sum(counts.values())
    if not total:
        return {}
    primary = {str(i) for i in range(1, 20)} | {"X", "Y", "MT"}
    unplaced = sum(v for k, v in counts.items() if k not in primary)
    return {
        "pct_MT": 100 * counts.get("MT", 0) / total,
        "pct_X": 100 * counts.get("X", 0) / total,
        "pct_Y": 100 * counts.get("Y", 0) / total,
        "pct_unplaced": 100 * unplaced / total,
    }


def fastq_screen(library: str) -> dict[str, float]:
    """The STP's own fastq_screen result for R1, as percent-mapped per genome.

    Independent of our pipeline entirely, which is what makes it useful as a
    cross-check on the mouse/human split. The file reports %Unmapped, so mapped
    is 100 - that. Filenames carry the S-number, so glob rather than construct:
    'ZHA8833A1_S100_...' must not be matched by a pattern for A10.
    """
    hits = sorted(SCREEN.glob(f"{library}_S*_R1_001_screen.txt"))
    if not hits:
        return {}
    out = {}
    with open(hits[0]) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 4 or f[0] in ("Genome", "#Fastq_screen"):
                continue
            try:
                out[f"screen_{f[0]}"] = 100 - float(f[3])
            except ValueError:
                continue
    return out


def main() -> None:
    meta = read_metadata()

    gen = read_mqc("multiqc_general_stats.txt")
    star = read_mqc("multiqc_star.txt")
    rsem = read_mqc("multiqc_rsem.txt")
    idx = read_mqc("multiqc_samtools_idxstats.txt")
    dups = read_mqc("multiqc_picard_dups.txt")
    fqc = read_mqc("multiqc_fastqc_fastqc_raw.txt")
    cut = read_mqc("multiqc_cutadapt.txt")

    rows = []
    for m in meta:
        lib = m["library"]
        g, s, r = gen.get(lib), star.get(lib), rsem.get(lib)
        # FastQC and cutadapt are per read file, so R1 is '<lib>_1'.
        f1, c1 = fqc.get(f"{lib}_1"), cut.get(f"{lib}_1")

        row = {
            "library": lib,
            "lims_name": m["lims_name"],
            "input_amount": m["input_amount"],
            "input_pg": int(m["input_pg"]),
            "replicate": m["replicate"],
            "well": m["well"],
            # sequencing
            "raw_pairs": num(f1, "Total Sequences"),
            "pct_gc": num(f1, "%GC"),
            "pct_trimmed_bp": num(c1, "percent_trimmed"),
            # alignment
            "star_uniquely_mapped_pct": num(s, "uniquely_mapped_percent"),
            "star_multimapped_pct": num(s, "multimapped_percent"),
            "star_unmapped_tooshort_pct": num(s, "unmapped_tooshort_percent"),
            "star_mismatch_rate": num(s, "mismatch_rate"),
            "star_num_splices": num(s, "num_splices"),
            "rsem_alignable_pct": num(r, "alignable_percent"),
            "rsem_unique": num(r, "Unique"),
            "rsem_multi": num(r, "Multi"),
            # duplication / library complexity
            # PERCENT_DUPLICATION in multiqc_picard_dups.txt is a FRACTION
            # (0.458); the general-stats copy is already scaled to a percent.
            # Take the scaled one so the column name is not a lie.
            "picard_dup_pct": num(g, "picard_mark_duplicates-PERCENT_DUPLICATION"),
            "dupradar_intercept": num(g, "custom_content_dupradar-dupRadar_intercept"),
            # composition
            # MultiQC's qualimap columns are read COUNTS IN MILLIONS, not
            # percentages -- for A1 'reads_aligned_exonic' is 39.46 against
            # 'reads_aligned' 27.23, which cannot be a percentage. They are
            # carried through as counts and converted below.
            "qualimap_exonic_M": num(g, "qualimap_rnaseq-reads_aligned_exonic"),
            "qualimap_intronic_M": num(g, "qualimap_rnaseq-reads_aligned_intronic"),
            "qualimap_intergenic_M": num(g, "qualimap_rnaseq-reads_aligned_intergenic"),
            "qualimap_5_3_bias": num(g, "qualimap_rnaseq-5_3_bias"),
            "insert_size_avg": num(g, "samtools_stats-insert_size_average"),
            "properly_paired_pct": num(g, "samtools_stats-reads_properly_paired_percent"),
            # rRNA AS NF-CORE REPORTS IT. This is Ensembl gene_biotype "rRNA"
            # only, which on GRCm39 is essentially 5S -- it is a ~5x
            # underestimate. 01_rrna_kmer_screen.py produces the real number;
            # the column is kept so the notebook can show them side by side.
            "nfcore_biotype_rRNA_pct": num(g, "custom_content_biotype_counts-percent_rRNA"),
        }
        # Genomic-origin percentages, over the three qualimap buckets only
        # (they partition the aligned reads; 'overlapping_exon' is a subset of
        # exonic and would double-count).
        origin = (row["qualimap_exonic_M"] + row["qualimap_intronic_M"]
                  + row["qualimap_intergenic_M"])
        for part in ("exonic", "intronic", "intergenic"):
            row[f"qualimap_{part}_pct"] = (
                100 * row[f"qualimap_{part}_M"] / origin if origin else "")

        row.update(idxstats_fractions(idx.get(lib, {})))
        row.update(fastq_screen(lib))
        rows.append(row)

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "qc_summary.tsv"
    fields = list(rows[0].keys())
    for r in rows:  # fastq_screen genome list is stable, but do not assume it
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"wrote {dest} ({len(rows)} libraries, {len(fields)} columns)")
    print("\ninput      library      raw_pairs   STARuniq%  RSEM%   dup%   nfcore_rRNA%")
    for r in rows:
        print(f"{r['input_amount']:<10} {r['library']:<12} {r['raw_pairs']:>10,.0f}"
              f"   {r['star_uniquely_mapped_pct']:>8.1f}"
              f"  {r['rsem_alignable_pct']:>6.1f}"
              f"  {r['picard_dup_pct']:>5.1f}"
              f"  {r['nfcore_biotype_rRNA_pct']:>12.2f}")


if __name__ == "__main__":
    main()
