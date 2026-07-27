#!/usr/bin/env python3
"""Sensitivity and replicate concordance across the input titration.

Two questions this answers, both of which the MultiQC report cannot:

  1. Is the gene-detection drop at low input real, or just sequencing depth?
     Answered by rarefying every library to the shallowest one's total count.
     Rather than actually resampling, the expected number of genes surviving
     binomial thinning at rate p is  sum over genes of  1 - (1-p)^count, which
     is exact in expectation and deterministic -- no seed to argue about.

  2. How reproducible is each input level? Answered by the Pearson correlation
     of log2(TPM+1) between the two replicates of each rung, over genes
     expressed in both. This is the number that shows where the titration
     breaks down.

Both outputs are keyed by input amount, which requires sample_metadata.tsv --
the count matrices themselves know nothing about the design.

Output
------
    res/flashseq/gene_detection.tsv     per library
    res/flashseq/replicate_concordance.tsv   per input level

Runs in under a minute on the merged RSEM tables. No sbatch needed.
"""
from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("FS_ROOT", "/nemo/lab/turnerj/working/guangxin/vasaseq"))
RESULTS = Path(os.environ.get("FS_RESULTS", ROOT / "data/flashseq/results"))
CODE = Path(os.environ.get("FS_CODE", ROOT / "code/flashseq"))
OUT = Path(os.environ.get("FS_OUT", ROOT / "res/flashseq"))

COUNTS = RESULTS / "star_rsem/rsem.merged.gene_counts.tsv"
TPM = RESULTS / "star_rsem/rsem.merged.gene_tpm.tsv"


def read_metadata() -> list[dict]:
    with open(CODE / "sample_metadata.tsv") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def read_matrix(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    """RSEM merged table -> (gene_ids, sample_names, values[gene][sample]).

    Column 0 is gene_id and column 1 is the comma-joined transcript list, so the
    data starts at column 2.
    """
    if not path.exists():
        sys.exit(f"FATAL: missing {path}")
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        samples = hdr[2:]
        genes, values = [], []
        for row in rdr:
            genes.append(row[0])
            values.append([float(x) for x in row[2:]])
    return genes, samples, values


def pearson(a: list[float], b: list[float]) -> float:
    n = len(a)
    if n < 2:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    return num / (da * db) if da and db else float("nan")


def main() -> None:
    meta = read_metadata()
    order = [m["library"] for m in meta]
    by_lib = {m["library"]: m for m in meta}

    genes, samples, counts = read_matrix(COUNTS)
    col = {s: i for i, s in enumerate(samples)}
    missing = [l for l in order if l not in col]
    if missing:
        sys.exit(f"FATAL: libraries in metadata but not in {COUNTS.name}: {missing}")

    totals = {l: sum(row[col[l]] for row in counts) for l in order}
    target = min(totals.values())
    print(f"{len(genes):,} genes; rarefying to {target:,.0f} counts "
          f"(the total of {min(totals, key=totals.get)})")

    rows = []
    for lib in order:
        i = col[lib]
        p = target / totals[lib]
        rows.append({
            "library": lib,
            "input_amount": by_lib[lib]["input_amount"],
            "input_pg": int(by_lib[lib]["input_pg"]),
            "replicate": by_lib[lib]["replicate"],
            "well": by_lib[lib]["well"],
            "total_counts": round(totals[lib]),
            "genes_ge1": sum(1 for r in counts if r[i] >= 1),
            "genes_ge5": sum(1 for r in counts if r[i] >= 5),
            # expected detections after binomial thinning to the common depth
            "genes_rarefied": round(sum(1 - (1 - p) ** r[i] for r in counts if r[i] > 0)),
            "rarefy_fraction": round(p, 4),
        })

    _, tpm_samples, tpm = read_matrix(TPM)
    tcol = {s: i for i, s in enumerate(tpm_samples)}
    for row in rows:
        i = tcol[row["library"]]
        row["genes_tpm_gt1"] = sum(1 for r in tpm if r[i] > 1)

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "gene_detection.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT/'gene_detection.tsv'}")

    # replicate concordance, one row per input level
    levels: dict[str, list[str]] = {}
    for m in meta:
        levels.setdefault(m["input_amount"], []).append(m["library"])

    conc = []
    for amount, libs in sorted(levels.items(), key=lambda kv: -int(by_lib[kv[1][0]]["input_pg"])):
        if len(libs) != 2:
            print(f"  skipping {amount}: {len(libs)} replicate(s), need 2")
            continue
        i, j = tcol[libs[0]], tcol[libs[1]]
        both = [(math.log2(r[i] + 1), math.log2(r[j] + 1))
                for r in tpm if r[i] > 1 and r[j] > 1]
        either = sum(1 for r in tpm if r[i] > 1 or r[j] > 1)
        conc.append({
            "input_amount": amount,
            "input_pg": int(by_lib[libs[0]]["input_pg"]),
            "rep1": libs[0],
            "rep2": libs[1],
            "pearson_log2tpm": round(pearson([x for x, _ in both], [y for _, y in both]), 4),
            "genes_in_both": len(both),
            "genes_in_either": either,
            "jaccard": round(len(both) / either, 4) if either else float("nan"),
        })

    with open(OUT / "replicate_concordance.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(conc[0]), delimiter="\t")
        w.writeheader()
        w.writerows(conc)
    print(f"wrote {OUT/'replicate_concordance.tsv'}\n")

    print("input      genes>=1  rarefied  TPM>1     r(reps)")
    for c in conc:
        libs = [r for r in rows if r["input_amount"] == c["input_amount"]]
        print(f"{c['input_amount']:<10} "
              f"{'/'.join(str(l['genes_ge1']) for l in libs):>13} "
              f"{'/'.join(str(l['genes_rarefied']) for l in libs):>13} "
              f"{'/'.join(str(l['genes_tpm_gt1']) for l in libs):>13} "
              f"{c['pearson_log2tpm']:>9.4f}")


if __name__ == "__main__":
    main()
