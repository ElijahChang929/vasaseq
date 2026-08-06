#!/usr/bin/env python3
"""
step2_report.py [CELLDIR] [LOGDIR]

Per-cell summary of what step 2 did. pipeline.sh runs this at the end of step2,
so the table lands in the step log; run it by hand any time afterwards to get
the same numbers back.

It reads these per cell, all already written by the run:

  <CELLDIR>/<sample>_<cell>_cbc_bcanchor.log    trim_bc_anchor.py --log (pass 0)
  <CELLDIR>/<sample>_<cell>_cbc_cutadapt.json   pass 2's --json report
  <LOGDIR>/step2_<sample>_<cell>.log            trim.sh's stdout (passes 1 and 2)

The JSON is preferred for pass 2's in/kept because it is keyed rather than
scraped; the log is the fallback, so cells trimmed before --json was added
(2026-08-03) still report identically.

Columns, and why there are two filters rather than one:

  in        reads out of step 1
  anchor    reads where pass 0 found this read's own barcode and cut there
  exact%    of those, how many matched the 12 nt with no mismatch
  ->tg      reads surviving TrimGalore. NOTE this is a real filter: TrimGalore
            applies its own --length 20 AFTER trimming, and on this library that
            alone removes ~20% of reads. It is easy to miss because the cutadapt
            summary TrimGalore prints says "100%" -- that is cutadapt's output,
            before TrimGalore's length cutoff is applied.
  kept      reads surviving pass 2 (cutadapt), i.e. the input to step 3
  end%      kept as a fraction of `in` -- the only end-to-end number

Expect the four low-count barcodes to stand out with a HIGH anchor rate and a
LOW end%: an empty well is mostly short-insert product, so more reads reach
their own barcode and more of them fall under the length floor once cleaned.
"""
import glob
import os
import json
import re
import sys


def main():
    celldir = sys.argv[1] if len(sys.argv) > 1 else "cells"
    logdir = sys.argv[2] if len(sys.argv) > 2 else "logs"

    rows = []
    for anchor_log in sorted(glob.glob(os.path.join(celldir, "*_cbc_bcanchor.log"))):
        base = os.path.basename(anchor_log)[: -len("_cbc_bcanchor.log")]
        m = re.search(r"_(\d+)$", base)
        cell = m.group(1) if m else base

        t = open(anchor_log, errors="replace").read()
        n = int(re.search(r"^reads (\d+)", t, re.M).group(1))
        anc = int(re.search(r"anchor found (\d+)", t).group(1))
        exact = int(re.search(r"of which exact (\d+)", t).group(1))

        # Pass 2's own numbers. Prefer the JSON report trim.sh now asks cutadapt
        # for: it is keyed, so it does not care how either tool words its output.
        # Fall back to scraping the text log for cells trimmed before --json was
        # added, so old runs still report.
        tg_in = kept = None
        cut_json = os.path.join(celldir, f"{base}_cbc_cutadapt.json")
        if os.path.exists(cut_json):
            try:
                with open(cut_json) as fh:
                    j = json.load(fh)
                tg_in = j["read_counts"]["input"]
                kept = j["read_counts"]["output"]
            except (ValueError, KeyError) as e:
                print(f"# warning: {cut_json} unreadable ({e}), falling back to log",
                      file=sys.stderr)

        step_log = os.path.join(logdir, f"step2_{base}.log")
        if tg_in is None and os.path.exists(step_log):
            body = open(step_log, errors="replace").read()
            # "Total reads processed" appears once per cutadapt invocation:
            # TrimGalore's (1.18) first, then pass 2's (5.1). Pass 2's input is
            # what survived TrimGalore's length cutoff.
            totals = re.findall(r"Total reads processed:\s+([\d,]+)", body)
            if len(totals) >= 2:
                tg_in = int(totals[-1].replace(",", ""))
            written = re.findall(r"Reads written \(passing filters\):\s+([\d,]+)", body)
            if written:
                kept = int(written[-1].replace(",", ""))
        rows.append((cell, n, anc, exact, tg_in, kept))

    if not rows:
        sys.exit(f"step2_report: no *_cbc_bcanchor.log under {celldir}")

    hdr = f"{'cell':>5} {'in':>13} {'anchor':>12} {'anc%':>6} {'exact%':>7} {'->tg':>13} {'kept':>13} {'end%':>6}"
    print(hdr)
    print("-" * len(hdr))
    tn = ta = tg = tk = 0
    for cell, n, anc, exact, tgi, kept in rows:
        f = lambda v: f"{v:,}" if v is not None else "?"
        pct = lambda a, b: f"{100*a/b:5.1f}%" if (a is not None and b) else "     ?"
        print(f"{cell:>5} {n:13,} {anc:12,} {pct(anc,n)} {pct(exact,anc)} "
              f"{f(tgi):>13} {f(kept):>13} {pct(kept,n)}")
        tn += n
        ta += anc
        tg += tgi or 0
        tk += kept or 0
    print("-" * len(hdr))
    print(f"{'ALL':>5} {tn:13,} {ta:12,} {100*ta/tn:5.1f}% {'':>7} "
          f"{tg:13,} {tk:13,} {100*tk/tn:5.1f}%")


if __name__ == "__main__":
    main()
