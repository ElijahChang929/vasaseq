#!/usr/bin/env python3
"""
02_nfcore_crosscheck.py -- what nf-core/RSEM says about FLASH-seq composition,
on its own terms, and what part of a FLASH-seq vs VASA difference is the
QUANTIFIER rather than the CHEMISTRY.

Three routes are read here. They are NOT interchangeable and the whole point of
the script is to keep them apart:

  A. nf-core RSEM        transcriptome alignment, EM-distributed multireads,
                         fractional expected_count per gene.
  B. nf-core featureCounts  genomic BAM, -g gene_biotype -t exon -s 0 -p -B -C,
                         multireads DISCARDED (Unassigned_MultiMapping), a read
                         spanning two biotypes DISCARDED (Unassigned_Ambiguity).
  C. VASA route          genomic BAM -> BED overlap with a per-GENE merged-exon
                         + explicit-intron + GtRNAdb-tRNA BED, biotype-hierarchy
                         multimapper rescue, '-'-joined combination entries.

A and B are the SAME READS under two quantifiers -> quantifier (pipeline) axis.
C on VASA reads vs C on FLASH-seq reads would be the chemistry (protocol) axis;
the FLASH-seq arm of C is produced by the other track and is NOT read here.

Denominators follow the user's decision: rRNA is quoted as % of ALL reads,
everything else as % of the NON-rRNA remainder.

Seeds: none (no stochastic step).
"""
import gzip, json, os, re, sys, hashlib
import pandas as pd
import numpy as np

print("pandas", pd.__version__, "numpy", np.__version__, file=sys.stderr)
assert hasattr(pd.DataFrame, "applymap"), "pandas 3.x would break the VASA pipeline env"

W    = "/nemo/lab/turnerj/working/guangxin/vasaseq"
FSR  = f"{W}/data/flashseq/results"
FC   = f"{FSR}/star_rsem/featurecounts"
VA   = f"{W}/data/PM26037/out/analysis"
SRC_GTF = "/nemo/lab/turnerj/working/guangxin/reference/gtf/Mus_musculus.GRCm39.116.gtf.gz"
NF_GTF  = f"{FSR}/genome/Mus_musculus.GRCm39.116.filtered.gtf"
OUT  = os.environ.get("OUTDIR", ".")
os.makedirs(OUT, exist_ok=True)

LIBS = [f"ZHA8833A{i}" for i in range(1, 11)]
# Annotation classes that are rRNA in the annotation sense. Named explicitly so
# the denominator is auditable rather than implied by a regex.
RRNA_CANON = {"rrna", "mtrrna", "rrnapseudogene"}

def canon(s):
    """Fold Ensembl snake_case and VASA CamelCase onto one key.
    '_' is dropped; '-' is NOT (it is VASA's combination separator)."""
    return s.replace("_", "").lower()

def md5(path, blocks=None):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()

log = []
def note(msg):
    log.append(msg)
    print(msg, file=sys.stderr)

# ---------------------------------------------------------------- 1. annotation
# gene_id -> gene_biotype from the SOURCE Ensembl 116 GTF, i.e. the same GTF
# VASA's v2 BED was built from (provenance.tsv: md5 0f9ab91d..., 78348 gene rows).
gid_re = re.compile(r'gene_id "([^"]+)"')
gbt_re = re.compile(r'gene_biotype "([^"]+)"')
gnm_re = re.compile(r'gene_name "([^"]+)"')

def gtf_gene_biotypes(path):
    op = gzip.open if path.endswith(".gz") else open
    d, names = {}, {}
    with op(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            g = gid_re.search(f[8]); b = gbt_re.search(f[8])
            if g and b:
                d[g.group(1)] = b.group(1)
                n = gnm_re.search(f[8])
                names[g.group(1)] = n.group(1) if n else g.group(1)
    return d, names

src_bt, src_nm = gtf_gene_biotypes(SRC_GTF)
nf_bt,  _      = gtf_gene_biotypes(NF_GTF)
note(f"source GTF {os.path.basename(SRC_GTF)}: {len(src_bt)} gene rows")
note(f"nf-core filtered GTF: {len(nf_bt)} gene rows (md5 {md5(NF_GTF)})")
shared = set(src_bt) & set(nf_bt)
disagree = [g for g in shared if src_bt[g] != nf_bt[g]]
note(f"gene_biotype agreement source vs nf-core-filtered: {len(shared)-len(disagree)}/{len(shared)} "
     f"shared genes agree; {len(disagree)} disagree; "
     f"{len(set(src_bt))-len(shared)} source-only, {len(set(nf_bt))-len(shared)} nfcore-only")

# 5S-only check for the rRNA biotype, re-derived rather than quoted
rrna_genes = sorted(g for g, b in src_bt.items() if b == "rRNA")
rrna_names = [src_nm.get(g, g) for g in rrna_genes]
n5s   = sum(1 for n in rrna_names if n.lower().startswith("n-r5s"))
n45s  = [n for n in rrna_names if re.match(r"(?i)rn(45s|28s|18s|5-?8s)", n)]
note(f"gene_biotype rRNA in Ensembl 116/GRCm39: {len(rrna_genes)} genes; "
     f"{n5s} named n-R5s*; non-5S names present: {sorted(set(n45s))}")

# ------------------------------------------------------- 2. sample metadata
meta = pd.read_csv(f"{W}/code/flashseq/sample_metadata.tsv", sep="\t", comment="#")
meta = meta.set_index("library")
assert set(LIBS) <= set(meta.index), "sample_metadata.tsv does not cover all 10 libraries"

# --------------------------------------- 3. nf-core featureCounts biotype tables
fc_rows, fc_summary = [], []
for lib in LIBS:
    p = f"{FC}/{lib}.biotype_counts_mqc.tsv"
    d = {}
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            k, v = line.rstrip("\n").split("\t")
            d[k] = float(v)
    fc_rows.append(pd.Series(d, name=lib))
    s = pd.read_csv(f"{FC}/{lib}.featureCounts.tsv.summary", sep="\t", index_col=0).iloc[:, 0]
    s.name = lib
    fc_summary.append(s)
fc = pd.DataFrame(fc_rows).fillna(0.0)                 # libraries x biotype, READS
fcs = pd.DataFrame(fc_summary).T                       # status x libraries
note(f"featureCounts biotype table: {fc.shape[0]} libraries x {fc.shape[1]} biotypes")

# Does the biotype table sum to Assigned? (tests the denominator before trusting it)
chk = pd.DataFrame({"biotype_sum": fc.sum(axis=1), "Assigned": fcs.loc["Assigned"]})
chk["diff"] = chk["biotype_sum"] - chk["Assigned"]
note("biotype_sum - Assigned per library: " +
     ", ".join(f"{k}:{v:+.0f}" for k, v in chk['diff'].items()))

# Reproduce nf-core's own percent_rRNA to learn its denominator
nf_pct = {}
for lib in LIBS:
    with open(f"{FC}/{lib}.biotype_counts_rrna_mqc.tsv") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Sample"):
                continue
            nf_pct[lib] = float(line.rstrip().split("\t")[1])
nf_pct = pd.Series(nf_pct)
cand = pd.DataFrame({
    "reported":        nf_pct,
    "rRNA/biotype_sum": 100 * fc["rRNA"] / fc.sum(axis=1),
    "rRNA/Assigned":    100 * fc["rRNA"] / fcs.loc["Assigned"],
    "(rRNA+Mt_rRNA)/biotype_sum": 100 * (fc["rRNA"] + fc["Mt_rRNA"]) / fc.sum(axis=1),
})
cand["resid_biotype_sum"] = (cand["reported"] - cand["rRNA/biotype_sum"]).abs()
note("nf-core percent_rRNA denominator: max |reported - rRNA/biotype_sum| = "
     f"{cand['resid_biotype_sum'].max():.3e}  (rRNA only, sum-of-biotypes denominator)")

# ------------------------------------------------------- 4. nf-core RSEM counts
rs = pd.read_csv(f"{FSR}/star_rsem/rsem.merged.gene_counts.tsv", sep="\t")
rs = rs.set_index("gene_id").drop(columns=["transcript_id(s)"])[LIBS]
note(f"RSEM merged gene_counts: {rs.shape[0]} genes x {rs.shape[1]} libraries; "
     f"pooled expected_count total {rs.values.sum():,.0f}")
rs_bt = pd.Series({g: src_bt.get(g) for g in rs.index})
missing = rs_bt.isna().sum()
note(f"RSEM genes with no gene_biotype in the source GTF: {missing} "
     f"(carrying {rs[rs_bt.isna().values].values.sum():,.1f} expected_count)")
rsem_by_bt = rs.groupby(rs_bt.fillna("UNANNOTATED").values).sum().T   # libs x biotype

# ------------------------------------------------------------ 5. VASA route
vr = pd.read_csv(f"{VA}/ZHA9292A1_analysis_total.ReadCounts.tsv", sep="\t", index_col=0)
gm = pd.read_csv(f"{VA}/gene_metadata.tsv", sep="\t", index_col=0)
assert vr.index.equals(gm.index), "ReadCounts and gene_metadata rows are not aligned"
CELLS = list(vr.columns)
note(f"VASA analysis ReadCounts: {vr.shape[0]} entries x {vr.shape[1]} cells "
     f"({CELLS[0]}..{CELLS[-1]}); pooled reads {vr.values.sum():,.0f}")

vr_tot = vr.sum(axis=1)                                  # pooled reads per entry
bt_lab = gm["biotype"].astype(str)
ncomp  = gm["n_genes_in_entry"].astype(int)
iscomb = gm["is_combination"].astype(str).eq("True")

# Two allocation rules for '-'-joined combination entries. Both are stated, both
# are reported; neither is presented as the truth.
#   PURE : keep only single-biotype entries  -> analogue of featureCounts, which
#          discards a read that spans two biotypes (Unassigned_Ambiguity)
#   SPLIT: divide an entry's reads equally among its distinct biotypes -> a crude,
#          prior-free analogue of RSEM's EM spreading a multiread over loci
pure_mask = ~bt_lab.str.contains("-")
vasa_pure = vr_tot[pure_mask].groupby(bt_lab[pure_mask]).sum()

split = {}
for lab, v in vr_tot.groupby(bt_lab).sum().items():
    parts = lab.split("-")
    for p in parts:
        split[p] = split.get(p, 0.0) + v / len(parts)
vasa_split = pd.Series(split)

# how much of each biotype's reads arrive via combination entries
comb_reads = {}
for lab, v in vr_tot[iscomb].groupby(bt_lab[iscomb]).sum().items():
    for p in lab.split("-"):
        comb_reads[p] = comb_reads.get(p, 0.0) + v / len(lab.split("-"))
comb_reads = pd.Series(comb_reads)

vasa_reads_total = float(vr_tot.sum())
note(f"VASA combination entries: {int(iscomb.sum())}/{len(iscomb)} rows "
     f"({100*iscomb.sum()/len(iscomb):.2f}%) carrying "
     f"{100*vr_tot[iscomb].sum()/vasa_reads_total:.2f}% of pooled reads; "
     f"n_genes_in_entry max {int(ncomp.max())}, "
     f"read-weighted mean {float((vr_tot*ncomp).sum()/vasa_reads_total):.3f}")

# per-cell VASA fractions (for a range, not just the pool)
vasa_percell_pure = vr[pure_mask.values].groupby(bt_lab[pure_mask].values).sum()

# --------------------------------------------------- 6. bwa-route rRNA, for context
bwa = None
try:
    bwa = pd.read_csv(f"{W}/res/flashseq/rrna_bwa.tsv", sep="\t")
    note(f"rrna_bwa.tsv columns: {list(bwa.columns)}")
except Exception as e:
    note(f"rrna_bwa.tsv not read: {e}")

# ============================================================ OUTPUT 1
# nfcore_composition.tsv -- long form, one row per library x biotype x route
def to_long(df, route, unit):
    tot = df.sum(axis=1)
    rr  = df[[c for c in df.columns if canon(c) in RRNA_CANON]].sum(axis=1)
    nonrr = tot - rr
    out = []
    for lib in df.index:
        for bt in df.columns:
            v = float(df.loc[lib, bt])
            out.append(dict(
                route=route, unit=unit, library=lib, biotype=bt,
                biotype_canon=canon(bt), count=v,
                pct_of_all=100 * v / tot[lib] if tot[lib] else np.nan,
                pct_of_nonrRNA=(np.nan if canon(bt) in RRNA_CANON
                                else (100 * v / nonrr[lib] if nonrr[lib] else np.nan)),
                total_all=float(tot[lib]), total_nonrRNA=float(nonrr[lib]),
            ))
    return pd.DataFrame(out)

comp = pd.concat([
    to_long(rsem_by_bt, "nfcore_rsem", "expected_count(EM reads)"),
    to_long(fc,         "nfcore_featurecounts", "assigned read pairs"),
], ignore_index=True)
for col in ["input_amount", "input_pg", "replicate", "well", "qc_verdict"]:
    comp[col] = comp["library"].map(meta[col])
comp["annotation_rRNA_is_5S_only"] = comp["biotype_canon"].eq("rrna")
comp = comp.sort_values(["route", "input_pg", "library", "count"],
                        ascending=[True, False, True, False])
comp.to_csv(f"{OUT}/nfcore_composition.tsv", sep="\t", index=False)
note(f"wrote nfcore_composition.tsv  {comp.shape[0]} rows x {comp.shape[1]} cols")

# ============================================================ OUTPUT 2
# nfcore_biotype_summary.tsv -- wide, per biotype: range across the 10 libraries
def summarise(df, route):
    tot = df.sum(axis=1)
    rr  = df[[c for c in df.columns if canon(c) in RRNA_CANON]].sum(axis=1)
    nonrr = tot - rr
    frac_all = 100 * df.div(tot, axis=0)
    frac_non = 100 * df.div(nonrr, axis=0)
    sc = ["ZHA8833A9", "ZHA8833A10"]
    rows = []
    for bt in df.columns:
        isrr = canon(bt) in RRNA_CANON
        f = frac_all[bt] if isrr else frac_non[bt]
        rows.append(dict(
            route=route, biotype=bt, biotype_canon=canon(bt),
            denominator="all reads" if isrr else "non-rRNA reads",
            pooled_pct=100 * df[bt].sum() / (tot.sum() if isrr else nonrr.sum()),
            min_pct=f.min(), max_pct=f.max(), median_pct=f.median(),
            pct_30pg_A9=f["ZHA8833A9"], pct_30pg_A10=f["ZHA8833A10"],
            mean_pct_30pg=f[sc].mean(),
            pooled_count=float(df[bt].sum()),
        ))
    return pd.DataFrame(rows)

summ = pd.concat([summarise(rsem_by_bt, "nfcore_rsem"),
                  summarise(fc, "nfcore_featurecounts")], ignore_index=True)
summ = summ.sort_values(["route", "pooled_pct"], ascending=[True, False])
summ.to_csv(f"{OUT}/nfcore_biotype_summary.tsv", sep="\t", index=False)
note(f"wrote nfcore_biotype_summary.tsv  {summ.shape[0]} rows")

# ============================================================ OUTPUT 3
# pipeline_vs_protocol.tsv
def fracs(series, name):
    s = series.copy()
    s.index = [canon(i) for i in s.index]
    s = s.groupby(level=0).sum()
    tot = s.sum()
    rr = s[[i for i in s.index if i in RRNA_CANON]].sum()
    non = tot - rr
    out = {}
    for i, v in s.items():
        out[i] = 100 * v / tot if i in RRNA_CANON else 100 * v / non
    return pd.Series(out, name=name)

sc = ["ZHA8833A9", "ZHA8833A10"]
f_rsem_all = fracs(rsem_by_bt.sum(axis=0), "nfcore_rsem_all10")
f_rsem_sc  = fracs(rsem_by_bt.loc[sc].sum(axis=0), "nfcore_rsem_30pg")
f_fc_all   = fracs(fc.sum(axis=0), "nfcore_featurecounts_all10")
f_fc_sc    = fracs(fc.loc[sc].sum(axis=0), "nfcore_featurecounts_30pg")
f_v_pure   = fracs(vasa_pure, "vasa_route_vasa_pure")
f_v_split  = fracs(vasa_split, "vasa_route_vasa_split")

pvp = pd.concat([f_rsem_all, f_rsem_sc, f_fc_all, f_fc_sc, f_v_pure, f_v_split], axis=1)
pvp.index.name = "biotype_canon"
pvp["vasa_route_flashseq_native"] = np.nan     # other track; not invented here
_disp = {}
for b in set(src_bt.values()) | set(bt_lab.unique()):
    for p in b.split("-"):
        _disp.setdefault(canon(p), p)
pvp["display_biotype"] = [_disp.get(i, i) for i in pvp.index]

# quantifier disagreement, measured on the SAME READS (RSEM vs featureCounts)
pvp["quantifier_delta_pp_30pg"]  = pvp["nfcore_rsem_30pg"] - pvp["nfcore_featurecounts_30pg"]
pvp["quantifier_ratio_30pg"]     = pvp["nfcore_rsem_30pg"] / pvp["nfcore_featurecounts_30pg"]
# total observed FLASH-seq(nf-core) vs VASA(VASA-route) gap: pipeline + protocol, confounded
pvp["total_gap_pp_rsem_vs_vasa"] = pvp["nfcore_rsem_30pg"] - pvp["vasa_route_vasa_pure"]
# VASA's own allocation-rule sensitivity: how much of its number is the '-' entries
pvp["vasa_allocation_delta_pp"]  = pvp["vasa_route_vasa_split"] - pvp["vasa_route_vasa_pure"]
pvp["vasa_pct_reads_via_combination_entries"] = pd.Series(
    {canon(k): 100 * v / vasa_reads_total for k, v in comb_reads.items()})
pvp["vasa_pooled_reads"] = pd.Series({canon(k): v for k, v in
                                      vr_tot.groupby(bt_lab).sum().items()})
pvp["decomposable"] = np.where(
    pvp[["nfcore_rsem_30pg", "vasa_route_vasa_pure"]].notna().all(axis=1),
    "pipeline term PENDING (vasa_route_flashseq_native not yet produced); "
    "only the confounded total and the within-nf-core quantifier delta are measured",
    "one side absent in this class")
pvp = pvp.sort_values("nfcore_rsem_30pg", ascending=False)
pvp.to_csv(f"{OUT}/pipeline_vs_protocol.tsv", sep="\t")
note(f"wrote pipeline_vs_protocol.tsv  {pvp.shape[0]} rows x {pvp.shape[1]} cols")

# ============================================================ OUTPUT 4
# multiread_disagreement.tsv -- where the two quantifiers differ, and whether it
# is the multi-locus small-RNA classes as predicted
mr = pd.DataFrame({
    "rsem_pct_30pg": pvp["nfcore_rsem_30pg"],
    "featurecounts_pct_30pg": pvp["nfcore_featurecounts_30pg"],
    "delta_pp": pvp["quantifier_delta_pp_30pg"],
    "log2_ratio": np.log2(pvp["quantifier_ratio_30pg"].replace(0, np.nan)),
})
mr["rsem_pooled_count"] = rsem_by_bt.sum(axis=0).groupby(
    [canon(c) for c in rsem_by_bt.columns]).sum()
mr["featurecounts_pooled_count"] = fc.sum(axis=0).groupby(
    [canon(c) for c in fc.columns]).sum()
mr["vasa_pct_reads_via_combination_entries"] = pvp["vasa_pct_reads_via_combination_entries"]
# mean n_genes_in_entry, read-weighted, for each biotype on the VASA side
w = {}
for lab, sub in vr_tot.groupby(bt_lab):
    n = ncomp.loc[sub.index]
    for p in lab.split("-"):
        w.setdefault(p, [0.0, 0.0])
        w[p][0] += float((sub * n).sum()) / len(lab.split("-"))
        w[p][1] += float(sub.sum()) / len(lab.split("-"))
mr["vasa_read_weighted_n_genes_in_entry"] = pd.Series(
    {canon(k): (a / b if b else np.nan) for k, (a, b) in w.items()})
SMALL = {"snorna", "snrna", "mirna", "scarna", "srna", "scrna", "miscrna",
         "rrna", "mtrrna", "rrnapseudogene", "trna", "mttrna", "ribozyme"}
mr["class"] = np.where(mr.index.isin(SMALL), "multi-locus small RNA", "other")
mr = mr.sort_values("delta_pp", key=abs, ascending=False)
mr.to_csv(f"{OUT}/multiread_disagreement.tsv", sep="\t")
note(f"wrote multiread_disagreement.tsv  {mr.shape[0]} rows")

# ============================================================ OUTPUT 5 diagnostics
diag = pd.concat([chk, cand, fcs.T], axis=1)
diag.index.name = "library"
for col in ["input_amount", "well", "qc_verdict"]:
    diag[col] = diag.index.map(meta[col])
diag.to_csv(f"{OUT}/nfcore_denominator_check.tsv", sep="\t")
note(f"wrote nfcore_denominator_check.tsv  {diag.shape[0]} rows x {diag.shape[1]} cols")

with open(f"{OUT}/crosscheck_log.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")

# ------------------------------------------------------------------ stdout facts
print("### rRNA biotype (5S-only) evidence")
print(f"rRNA genes in Ensembl 116: {len(rrna_genes)}; n-R5s* named: {n5s}; "
      f"non-5S names: {sorted(set(n45s))}")
print("\n### nf-core percent_rRNA reproduced")
print(cand.round(4).to_string())
print("\n### RSEM biotype fractions, 30 pg rung (A9+A10 pooled), top 14")
print(f_rsem_sc.sort_values(ascending=False).head(14).round(4).to_string())
print("\n### featureCounts biotype fractions, 30 pg rung, top 14")
print(f_fc_sc.sort_values(ascending=False).head(14).round(4).to_string())
print("\n### VASA route (pure-entry rule), pooled 12 cells, top 14")
print(f_v_pure.sort_values(ascending=False).head(14).round(4).to_string())
print("\n### largest quantifier disagreements (RSEM - featureCounts, pp, 30 pg)")
print(mr[["rsem_pct_30pg", "featurecounts_pct_30pg", "delta_pp", "log2_ratio",
          "vasa_read_weighted_n_genes_in_entry", "class"]].head(16).round(4).to_string())
print("\n### featureCounts fate of reads, per library")
print(fcs.loc[["Assigned", "Unassigned_MultiMapping", "Unassigned_Ambiguity",
               "Unassigned_NoFeatures", "Unassigned_Unmapped"]].T.to_string())
print("\n### multiread share of alignable reads, per library")
mm = fcs.loc["Unassigned_MultiMapping"]
al = fcs.loc["Assigned"] + fcs.loc["Unassigned_MultiMapping"] + \
     fcs.loc["Unassigned_Ambiguity"] + fcs.loc["Unassigned_NoFeatures"]
print((100 * mm / al).round(3).to_string())
print("\n### VASA combination-entry read share, top classes")
print((100 * comb_reads / vasa_reads_total).sort_values(ascending=False).head(12).round(4).to_string())
if bwa is not None:
    print("\n### bwa-route rRNA (the actual rRNA measurement), head")
    print(bwa.head(24).to_string())
