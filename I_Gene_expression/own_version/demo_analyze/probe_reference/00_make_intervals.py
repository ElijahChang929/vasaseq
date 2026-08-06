#!/usr/bin/env python3
"""00_make_intervals.py -- HOW probe_target_intervals.mouse.tsv was produced.

This is the provenance script, kept so the interval file is reproducible rather
than a magic input. Running it regenerates probe_target_intervals.mouse.tsv.

    python3 00_make_intervals.py            # needs: biopython, network (NCBI)

WHAT IT DOES, AND WHY
---------------------
The wet-lab depletion uses the Adiconis et al. 2013 RNase H oligo set: 195
50-mers tiling HUMAN rRNA (Nat Methods 10:623, doi:10.1038/nmeth.2483). Our
sample is MOUSE. A human 50-mer only drives RNase H cleavage where it still
base-pairs, so the probe-reachable footprint on mouse is not "all of mouse rRNA":
it must be computed.

Method: fetch the human rRNA templates and the mouse reference contigs, align
each mouse contig to its best-matching human template, and keep every 50 nt
window retaining >= 90 % identity (<= 5 mismatches). Those windows are the
scoreable footprint. Everything else is excluded from the metric because the
probes cannot reach it.

Nuclear subunits stay largely addressable; mitochondrial 12S/16S are much less
conserved and are mostly out of reach. That is a probe-design limit, not a
failure of the reaction, which is exactly why the metric is windowed.

NOTE: the probe SEQUENCES themselves were never obtainable -- the supplementary
table of Adiconis 2013 could not be downloaded from PMC, EuropePMC, or the
publisher (all returned HTML or empty responses, checked 2026-08-02). So the
footprint is computed from the human rRNA TEMPLATE SPACE the 195 50-mers tile,
not from the individual oligos. This is an approximation: it assumes the oligo
set tiles its targets without large gaps. Consequence: the footprint is an upper
bound on true probe reach, so the reported residual is, if anything, slightly
optimistic. Treat cross-plate comparisons (which share the assumption) as sound
and the absolute value as approximate.
"""
import numpy as np, pandas as pd, urllib.parse, urllib.request, os
from collections import Counter
from Bio import SeqIO, Align
from io import StringIO

THR, W = 0.90, 50          # identity threshold, probe length
HERE = os.path.dirname(os.path.abspath(__file__))
MOUSE_FA = "/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa"

# human rRNA templates the Adiconis oligo set tiles, + mouse 47S for coordinates
# accessions verified at NCBI 2026-07-31 and re-verified 2026-08-04.
# Mito rRNAs are subranges of the human mitochondrial genome, NOT separate RefSeq
# records: NC_012920.1:648-1601 is gene RNR1 / product "s-rRNA" (12S, 954 nt) and
# 1671-3229 is RNR2 / "l-rRNA" (16S, 1559 nt) -- checked against the GenBank
# feature table so the two are not transposed.
ACC = {"18S":("NR_003286.4", None, None),
       "5.8S":("NR_003285.3", None, None),
       "28S":("NR_003287.4", None, None),
       "5S":("V00589.1", None, None),
       "mt12S":("NC_012920.1", 648, 1601),
       "mt16S":("NC_012920.1", 1671, 3229)}
M47S_ACC, M47S_RANGE = "BK000964.3", (1, 13403)
SUB47 = {"18S":(4008,5877), "5.8S":(6878,7034), "28S":(8123,12852)}   # in 47S coords

def efetch(acc, start=None, stop=None):
    p = dict(db="nuccore", id=acc, rettype="fasta", retmode="text")
    if start: p.update(seq_start=start, seq_stop=stop)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(p)
    with urllib.request.urlopen(url, timeout=120) as f:
        return str(next(SeqIO.parse(StringIO(f.read().decode()), "fasta")).seq).upper()

print("fetching human templates from NCBI ...")
HUMAN = {k: efetch(a, ss, se) for k, (a, ss, se) in ACC.items()}
M47S  = efetch(M47S_ACC, *M47S_RANGE)
refseqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(MOUSE_FA, "fasta")}
c47 = [k for k in refseqs if "rDNA_47S" in k][0]
assert refseqs[c47] == M47S, "47S contig differs from NCBI BK000964.3:1-13403"
print(f"  {len(refseqs)} mouse contigs; 47S verified against {M47S_ACC}")

def mk(mode):
    a = Align.PairwiseAligner(); a.mode = mode
    a.match_score, a.mismatch_score = 2, -1
    a.open_gap_score, a.extend_gap_score = -5, -0.5
    return a
glob, loc = mk("global"), mk("local")

def best_windows(match):
    """Per-base best 50 nt window identity."""
    n = len(match)
    if n < W: return np.full(n, match.mean())
    cs = np.concatenate([[0], np.cumsum(match)])
    st = np.arange(n - W + 1)
    wid = (cs[st + W] - cs[st]) / W
    best = np.zeros(n)
    for s, v in zip(st, wid): best[s:s+W] = np.maximum(best[s:s+W], v)
    return best

def project(template, seq, aligner):
    """Align seq to template; return per-seq-base match bool."""
    aln = aligner.align(template, seq)[0]
    t_al, s_al = str(aln[0]), str(aln[1])
    m, i = np.zeros(len(seq), bool), 0
    for tc, sc in zip(t_al, s_al):
        if sc != "-":
            m[i] = (tc == sc); i += 1
    return m, aln.score

def runs(mask, offset=0, min_len=20):
    out, s = [], None
    for i, v in enumerate(list(mask) + [False]):
        if v and s is None: s = i
        elif not v and s is not None:
            if i - s >= min_len: out.append((s + offset + 1, i + offset))
            s = None
    return out

bed = []
# 47S nuclear subunits: global projection, in 47S coordinates
for k in ("18S", "5.8S", "28S"):
    lo, hi = SUB47[k]
    m, _ = project(HUMAN[k], M47S[lo-1:hi], glob)
    for s, e in runs(best_windows(m) >= THR, offset=lo-1):
        bed.append((c47, s, e, f"probe_{k}"))

# every other contig: assign to its best-scoring human template by local alignment
print("assigning remaining contigs to templates ...")
n_by = Counter()
for cid, seq in refseqs.items():
    if cid == c47: continue
    scored = [(project(HUMAN[t], seq, loc), t) for t in HUMAN]
    (m, sc), t = max(scored, key=lambda x: x[0][1])
    n_by[t] += 1
    for s, e in runs(best_windows(m) >= THR):
        bed.append((cid, s, e, f"probe_{t}"))

df = pd.DataFrame(bed, columns=["contig","start","end","probe_target"])
df["length"] = df["end"] - df["start"] + 1
df = df.sort_values(["contig","start"])
df.to_csv(f"{HERE}/probe_target_intervals.mouse.tsv", sep="\t", index=False)
print(f"\nwrote probe_target_intervals.mouse.tsv: {len(df)} intervals, "
      f"{df['length'].sum():,} bp over {df['contig'].nunique()} contigs")
print(df.groupby("probe_target")["length"].agg(["size","sum"]).to_string())
