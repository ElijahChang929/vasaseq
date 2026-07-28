#!/usr/bin/env python3
"""
07_vasa_comparator.py -- which VASA table is the fair comparator for nf-core,
and what each choice costs. Then write the final pipeline_vs_protocol.tsv.

The VASA side has THREE defensible read-based compositions and they differ by
more than the FLASH-seq/VASA effect being measured, so the choice cannot be made
silently:

  V1 uniagg total      per-GENE, uniquely-assigned, exon+intron, no filters.
  V2 uniagg spliced    same, exon only -> the fair analogue of featureCounts
                       -t exon and of RSEM's transcriptome, neither of which can
                       count an intronic read.
  V3 analysis pure     the analysis set: UMI-ceiling genes dropped, 4 blank cells
                       dropped, single-biotype entries only.

V3 is the table the rest of this project uses, but its UMI-ceiling filter removes
8 entries chosen because they SATURATE a 4096-UFI ceiling -- and those 8 are
Rmrp, Rn7sk, Rn7s1/Rn7s2, three Rnu1/Rnu2 snRNAs and the Snord3b cluster, i.e.
the most abundant small RNAs in the library. That filter is correct for molecule
counting and wrong for a read-based composition comparison. This script measures
how much it moves each class instead of asserting it.

Also identifies the genes behind the largest quantifier disagreement, to test the
multi-locus prediction at gene rather than class level.

Seeds: none.
"""
import gzip, json, os, re, sys
import pandas as pd, numpy as np

W   = "/nemo/lab/turnerj/working/guangxin/vasaseq"
FSR = f"{W}/data/flashseq/results"
RAW = f"{W}/data/PM26037/out"
VA  = f"{RAW}/analysis"
SRC_GTF = "/nemo/lab/turnerj/working/guangxin/reference/gtf/Mus_musculus.GRCm39.116.gtf.gz"
OUT = os.environ.get("OUTDIR", ".")
LIBS = [f"ZHA8833A{i}" for i in range(1, 11)]
SC   = ["ZHA8833A9", "ZHA8833A10"]
CELLS = ["%03d" % i for i in range(2, 14)]
RRNA_CANON = {"rrna", "mtrrna", "rrnapseudogene"}
log = []
def note(m): log.append(str(m)); print(m, file=sys.stderr, flush=True)
def canon(s): return s.replace("_", "").lower()

gid = re.compile(r'gene_id "([^"]+)"'); gbt = re.compile(r'gene_biotype "([^"]+)"')
gnm = re.compile(r'gene_name "([^"]+)"')
src_bt, src_nm = {}, {}
with gzip.open(SRC_GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene": continue
        g = gid.search(f[8]).group(1)
        src_bt[g] = gbt.search(f[8]).group(1)
        n = gnm.search(f[8]); src_nm[g] = n.group(1) if n else g

def rd(p):
    d = pd.read_csv(p, sep="\t", index_col=0)
    d = d[~pd.isna(d.index)]; d.index = d.index.astype(str)
    return d.reindex(columns=CELLS).fillna(0.0)

u_tot = rd(f"{RAW}/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv")
u_spl = rd(f"{RAW}/ZHA9292A1_uniaggGenes_spliced.ReadCounts.tsv").reindex(u_tot.index).fillna(0.0)
ubt = pd.Series({e: e.split("_")[-1] for e in u_tot.index})

CEIL = json.load(open(f"{VA}/manifest.json"))["filters"]["umi_ceiling_genes_dropped"]
note(f"UMI-ceiling entries dropped by the analysis filter ({len(CEIL)}): {CEIL}")
# map each ceiling entry back onto uniagg per-gene rows (they are '-'-joined multiagg keys)
ceil_genes = set()
for e in CEIL:
    for part in e.split("-"):
        m = re.match(r"(ENSMUSG\d+)_", part)
        if m: ceil_genes.add(m.group(1))
ceil_rows = [i for i in u_tot.index if i.split("_")[0] in ceil_genes]
note(f"uniagg rows matching those genes: {len(ceil_rows)}; "
     f"they carry {u_tot.loc[ceil_rows].values.sum():,.0f} reads = "
     f"{100*u_tot.loc[ceil_rows].values.sum()/u_tot.values.sum():.2f}% of uniagg total")
note("per-entry: " + ", ".join(
    f"{i.split('_')[1]}:{u_tot.loc[i].sum():,.0f}"
    for i in sorted(ceil_rows, key=lambda x: -u_tot.loc[x].sum())[:12]))

def compose(s, bt, name):
    g = s.groupby(bt.reindex(s.index)).sum()
    g.index = [canon(i) for i in g.index]; g = g.groupby(level=0).sum()
    tot = g.sum(); rr = g[[i for i in g.index if i in RRNA_CANON]].sum(); non = tot - rr
    return pd.Series({i: (100*v/tot if i in RRNA_CANON else 100*v/non) for i, v in g.items()},
                     name=name)

keep = [i for i in u_tot.index if i not in set(ceil_rows)]
V1  = compose(u_tot.sum(axis=1), ubt, "V1_vasa_uniagg_total")
V2  = compose(u_spl.sum(axis=1), ubt, "V2_vasa_uniagg_spliced")
V1c = compose(u_tot.loc[keep].sum(axis=1), ubt, "V1_minus_umi_ceiling")
V2c = compose(u_spl.loc[keep].sum(axis=1), ubt, "V2_minus_umi_ceiling")

a_tot = pd.read_csv(f"{VA}/ZHA9292A1_analysis_total.ReadCounts.tsv", sep="\t", index_col=0)
gm = pd.read_csv(f"{VA}/gene_metadata.tsv", sep="\t", index_col=0)
abt = gm["biotype"].astype(str); pure = ~abt.str.contains("-")
V3 = compose(a_tot.sum(axis=1)[pure], abt[pure], "V3_vasa_analysis_pure")

# top carriers per small-RNA class on the VASA side (gene-level, so the reader can check)
carriers = []
for cls in ["MiscRna", "snRNA", "snoRNA", "ribozyme", "rRNA", "scaRNA", "miRNA", "lncRNA"]:
    sub = u_tot[ubt.reindex(u_tot.index).eq(cls).values].sum(axis=1).sort_values(ascending=False)
    for e in sub.head(4).index:
        carriers.append(dict(biotype=cls, entry=e, symbol=e.split("_")[1],
                             reads=float(sub[e]),
                             pct_of_class=100*sub[e]/sub.sum() if sub.sum() else np.nan,
                             dropped_by_umi_ceiling=e in set(ceil_rows)))
car = pd.DataFrame(carriers)
car.to_csv(f"{OUT}/vasa_class_carriers.tsv", sep="\t", index=False)
note("VASA top carriers per class:\n" + car.round(3).to_string(index=False))

# ---- nf-core sides
fc = pd.DataFrame([{k: float(v) for k, v in
    (l.rstrip("\n").split("\t") for l in
     open(f"{FSR}/star_rsem/featurecounts/{lib}.biotype_counts_mqc.tsv") if not l.startswith("#"))}
    for lib in LIBS], index=LIBS).fillna(0.0)
rs = pd.read_csv(f"{FSR}/star_rsem/rsem.merged.gene_counts.tsv", sep="\t").set_index("gene_id")
rs = rs.drop(columns=["transcript_id(s)"])[LIBS]
rsb = rs.groupby(pd.Series({g: src_bt.get(g, "UNANNOTATED") for g in rs.index}).values).sum().T

def cw(df, name, libs=None):
    s = (df.loc[libs] if libs else df).sum(axis=0)
    s.index = [canon(i) for i in s.index]; s = s.groupby(level=0).sum()
    tot = s.sum(); rr = s[[i for i in s.index if i in RRNA_CANON]].sum(); non = tot - rr
    return pd.Series({i: (100*v/tot if i in RRNA_CANON else 100*v/non) for i, v in s.items()},
                     name=name)
R_sc, R_al = cw(rsb, "nfcore_rsem_30pg", SC), cw(rsb, "nfcore_rsem_all10")
F_sc, F_al = cw(fc, "nfcore_featurecounts_30pg", SC), cw(fc, "nfcore_featurecounts_all10")

# gene-level test of the multiread prediction on the worst class
note("\n=== gene-level: who carries misc_RNA, the largest relative disagreement ===")
g9 = pd.read_csv(f"{FSR}/star_rsem/ZHA8833A9.genes.results", sep="\t").set_index("gene_id")
g9["biotype"] = pd.Series(src_bt); g9["gene_name"] = pd.Series(src_nm)
mm = g9[g9.biotype.eq("misc_RNA")].sort_values("expected_count", ascending=False).head(6)
note(mm[["gene_name", "length", "effective_length", "expected_count"]].round(2).to_string())
note(f"featureCounts A9 misc_RNA total reads: {fc.loc['ZHA8833A9','misc_RNA']:.0f}  "
     f"vs RSEM {g9[g9.biotype.eq('misc_RNA')].expected_count.sum():,.0f} expected_count "
     f"-- ratio {g9[g9.biotype.eq('misc_RNA')].expected_count.sum()/max(fc.loc['ZHA8833A9','misc_RNA'],1):,.0f}x")

# ============================================================ final deliverable
pvp = pd.concat([R_al, R_sc, F_al, F_sc, V1, V2, V1c, V2c, V3], axis=1)
pvp.index.name = "biotype_canon"
disp = {}
for b in set(src_bt.values()) | set(abt.unique()) | set(ubt.unique()):
    for p in str(b).split("-"): disp.setdefault(canon(p), p)
pvp.insert(0, "display_biotype", [disp.get(i, i) for i in pvp.index])
pvp["denominator"] = np.where(pvp.index.isin(RRNA_CANON), "pct of ALL reads",
                              "pct of NON-rRNA reads")
pvp["A_quantifier_pp"] = pvp.nfcore_rsem_30pg - pvp.nfcore_featurecounts_30pg
pvp["A_quantifier_log2"] = np.log2((pvp.nfcore_rsem_30pg / pvp.nfcore_featurecounts_30pg
                                    .replace(0, np.nan)).replace(0, np.nan))
pvp["B_intron_model_pp"] = pvp.V1_vasa_uniagg_total - pvp.V2_vasa_uniagg_spliced
pvp["B2_umi_ceiling_filter_pp"] = pvp.V1_minus_umi_ceiling - pvp.V1_vasa_uniagg_total
pvp["C_total_gap_pp"] = pvp.nfcore_rsem_30pg - pvp.V1_vasa_uniagg_total
pvp["C_exon_matched_gap_pp"] = pvp.nfcore_rsem_30pg - pvp.V2_vasa_uniagg_spliced
pvp["vasa_route_flashseq_native"] = np.nan
pvp["protocol_term_pp"] = np.nan            # = V2 - vasa_route_flashseq_native, pending
def cls(i):
    if i in RRNA_CANON: return "rRNA (annotation route; see rrna_biotype_note)"
    if i in {"snorna","snrna","mirna","scarna","srna","scrna","miscrna","trna","mttrna","ribozyme"}:
        return "short non-poly-A"
    if "pseudogene" in i: return "pseudogene"
    if i in {"proteincoding","lncrna","tec"}: return "long poly-A / lncRNA"
    return "other"
pvp["class"] = [cls(i) for i in pvp.index]
pvp["decomposes_cleanly"] = np.where(
    pvp[["nfcore_rsem_30pg", "V2_vasa_uniagg_spliced"]].notna().all(axis=1),
    "PARTIAL: A (quantifier) and B (annotation model) measured; protocol term needs "
    "vasa_route_flashseq_native, not yet produced",
    "NO: class present on only one side")
pvp = pvp.sort_values("nfcore_rsem_30pg", ascending=False)
pvp.to_csv(f"{OUT}/pipeline_vs_protocol.tsv", sep="\t")

KEY = ["proteincoding","lncrna","miscrna","snrna","snorna","ribozyme","rrna","mtrrna",
       "processedpseudogene","transcribedprocessedpseudogene","unprocessedpseudogene",
       "mirna","scarna","mttrna","tec"]
show = [k for k in KEY if k in pvp.index]
note("\n=== FINAL pipeline_vs_protocol, key classes ===\n" + pvp.loc[show, [
    "display_biotype","nfcore_rsem_30pg","nfcore_featurecounts_30pg",
    "V1_vasa_uniagg_total","V2_vasa_uniagg_spliced","V3_vasa_analysis_pure",
    "A_quantifier_pp","B_intron_model_pp","B2_umi_ceiling_filter_pp",
    "C_exon_matched_gap_pp"]].round(4).to_string())
note(f"\nwrote pipeline_vs_protocol.tsv {pvp.shape}")
with open(f"{OUT}/comparator_log.txt", "w") as fh: fh.write("\n".join(log) + "\n")
