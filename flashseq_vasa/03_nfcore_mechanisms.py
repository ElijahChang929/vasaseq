#!/usr/bin/env python3
"""
03_nfcore_mechanisms.py -- WHY nf-core/RSEM and the VASA route disagree.

Four mechanisms, each measured rather than asserted:
  (a) what gene_biotype "rRNA" actually contains in Ensembl 116/GRCm39, and
      which of those genes carry nf-core's rRNA signal  -> is the 5S label right?
  (b) what nf-core's FILTERED GTF contains (the 0-gene-rows result from script 02
      was a parse artefact: the filter emits no `gene` feature lines)
  (c) RSEM's effective_length floor: a transcript shorter than the fragment
      length gets effective_length ~0 and CANNOT receive reads. This is a
      quantifier property, not chemistry, and it hits exactly the short
      non-poly-A species.
  (d) the exon-vs-exon+intron axis: the VASA route counts intronic reads,
      nf-core (featureCounts -t exon; RSEM transcriptome) cannot.
Plus the quantitative test of the multi-locus prediction.
"""
import gzip, os, re, sys
import pandas as pd, numpy as np

SEED = 0   # only stochastic element: the permutation null for the Spearman p-value

def _rank(a):
    a = np.asarray(a, float); o = a.argsort(kind="mergesort")
    r = np.empty(len(a), float); r[o] = np.arange(1, len(a) + 1)
    # average ties, so the rho below is the tie-corrected Spearman
    u, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    for i in np.where(cnt > 1)[0]:
        m = inv == i
        r[m] = r[m].mean()
    return r

def spearmanr(x, y, n_perm=20000, seed=SEED):
    """Tie-corrected Spearman rho with a permutation p-value (no scipy in envs/vasa)."""
    rx, ry = _rank(x), _rank(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        if abs(np.corrcoef(rng.permutation(rx), ry)[0, 1]) >= abs(rho) - 1e-12:
            ge += 1
    return rho, (ge + 1) / (n_perm + 1)

W   = "/nemo/lab/turnerj/working/guangxin/vasaseq"
FSR = f"{W}/data/flashseq/results"
RAW = f"{W}/data/PM26037/out"
VA  = f"{RAW}/analysis"
SRC_GTF = "/nemo/lab/turnerj/working/guangxin/reference/gtf/Mus_musculus.GRCm39.116.gtf.gz"
NF_GTF  = f"{FSR}/genome/Mus_musculus.GRCm39.116.filtered.gtf"
OUT = os.environ.get("OUTDIR", ".")
LIBS = [f"ZHA8833A{i}" for i in range(1, 11)]
log = []
def note(m):
    log.append(m); print(m, file=sys.stderr)

def canon(s): return s.replace("_", "").lower()

# ---------------------------------------------------------------- (a) rRNA genes
gid = re.compile(r'gene_id "([^"]+)"'); gbt = re.compile(r'gene_biotype "([^"]+)"')
gnm = re.compile(r'gene_name "([^"]+)"'); tid = re.compile(r'transcript_id "([^"]+)"')
src_bt, src_nm = {}, {}
with gzip.open(SRC_GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "gene": continue
        g = gid.search(f[8]); b = gbt.search(f[8]); n = gnm.search(f[8])
        src_bt[g.group(1)] = b.group(1)
        src_nm[g.group(1)] = n.group(1) if n else g.group(1)

rr = pd.DataFrame({"gene_id": [g for g, b in src_bt.items() if b == "rRNA"]})
rr["gene_name"] = rr.gene_id.map(src_nm)
def fam(n):
    nl = n.lower()
    if nl.startswith("n-r5s"): return "n-R5s* (5S)"
    if re.match(r"^rn18s", nl): return "Rn18s* (18S relic)"
    if re.match(r"^rn(45s|28s|5-?8s|5s)", nl): return "Rn45s/28S/5.8S/5S"
    if nl.startswith("gm"):  return "Gm* (unnamed)"
    if nl.startswith("ensmusg"): return "no gene_name"
    return "other"
rr["family"] = rr.gene_name.map(fam)
note("gene_biotype=rRNA name families:\n" + rr.family.value_counts().to_string())

# nf-core's own rRNA signal, per gene: which rRNA genes actually get counts?
gr = pd.read_csv(f"{FSR}/star_rsem/ZHA8833A9.genes.results", sep="\t")
gr = gr.set_index("gene_id")
gr["biotype"] = pd.Series(src_bt)
gr["gene_name"] = pd.Series(src_nm)
rrg = gr[gr.biotype.eq("rRNA")].sort_values("expected_count", ascending=False)
note(f"A9 RSEM rRNA-biotype expected_count total {rrg.expected_count.sum():,.1f}; "
     f"{int((rrg.expected_count>0).sum())}/{len(rrg)} genes nonzero")
note("A9 top 12 rRNA-biotype genes by RSEM expected_count:\n" +
     rrg.head(12)[["gene_name", "length", "effective_length", "expected_count"]].to_string())
top5s = rrg.head(20).gene_name.str.lower().str.startswith("n-r5s").sum()
share5s = rrg[rrg.gene_name.str.lower().str.startswith("n-r5s")].expected_count.sum() / \
          max(rrg.expected_count.sum(), 1e-12)
note(f"share of A9 RSEM rRNA-biotype signal on n-R5s* genes: {100*share5s:.2f}%")

# ------------------------------------------------------- (b) filtered GTF content
ftypes, nf_tx, nf_genes = {}, set(), set()
with open(NF_GTF) as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.split("\t", 9)
        if len(f) < 9: continue
        ftypes[f[2]] = ftypes.get(f[2], 0) + 1
        g = gid.search(f[8]); t = tid.search(f[8])
        if g: nf_genes.add(g.group(1))
        if t: nf_tx.add(t.group(1))
note(f"nf-core filtered GTF feature types: {ftypes}")
note(f"filtered GTF: {len(nf_genes)} distinct gene_id, {len(nf_tx)} distinct transcript_id")
note(f"source GTF genes not in filtered GTF: {len(set(src_bt)-nf_genes)}")
lost = pd.Series({g: src_bt[g] for g in set(src_bt) - nf_genes})
if len(lost):
    note("biotypes of genes dropped by nf-core's GTF filter:\n" +
         lost.value_counts().head(12).to_string())

# ------------------------------------------- (c) RSEM effective_length floor
gr["eff_zero"] = gr.effective_length <= 1.0
byb = gr.groupby("biotype").agg(
    n_genes=("length", "size"),
    median_length=("length", "median"),
    median_eff_length=("effective_length", "median"),
    n_eff_le1=("eff_zero", "sum"),
    expected_count=("expected_count", "sum"))
byb["pct_genes_eff_le1"] = 100 * byb.n_eff_le1 / byb.n_genes
byb = byb.sort_values("expected_count", ascending=False)
byb.to_csv(f"{OUT}/rsem_effective_length.tsv", sep="\t")
SMALL = ["rRNA", "snRNA", "snoRNA", "miRNA", "scaRNA", "misc_RNA", "Mt_tRNA",
         "Mt_rRNA", "ribozyme", "scRNA", "sRNA", "protein_coding", "lncRNA"]
note("RSEM effective_length by biotype (A9), key classes:\n" +
     byb.loc[[b for b in SMALL if b in byb.index]].round(2).to_string())

# ---------------------------------------------- (d) exon vs exon+intron on VASA
tot = pd.read_csv(f"{RAW}/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv", sep="\t", index_col=0)
spl = pd.read_csv(f"{RAW}/ZHA9292A1_uniaggGenes_spliced.ReadCounts.tsv", sep="\t", index_col=0)
uns = pd.read_csv(f"{RAW}/ZHA9292A1_uniaggGenes_unspliced.ReadCounts.tsv", sep="\t", index_col=0)
CELLS = ["%03d" % i for i in range(2, 14)]
tot, spl, uns = [d.reindex(index=tot.index)[CELLS].fillna(0) for d in (tot, spl, uns)]
excess = (spl + uns - tot)
note(f"uniagg spliced+unspliced-total: max excess {float(excess.values.max()):.1f}, "
     f"n cells with any excess {int((excess.values>0).any(axis=0).sum())} "
     f"(vasa branch uses exact membership; nonzero would be trap (h))")
vb = pd.Series({e: e.split("_")[-1] for e in tot.index})
agg = pd.DataFrame({"total": tot.sum(axis=1), "spliced": spl.sum(axis=1),
                    "unspliced": uns.sum(axis=1), "biotype": vb}).groupby("biotype").sum(
                    numeric_only=True)
agg["pct_unspliced"] = 100 * agg.unspliced / agg.total.replace(0, np.nan)
agg["pct_of_all_reads"] = 100 * agg.total / agg.total.sum()
agg = agg.sort_values("total", ascending=False)
agg.to_csv(f"{OUT}/vasa_intron_axis.tsv", sep="\t")
note("VASA uniagg reads, intron share by biotype (top 14):\n" +
     agg.head(14).round(3).to_string())
note(f"VASA uniagg overall intron (unspliced) share of reads: "
     f"{100*agg.unspliced.sum()/agg.total.sum():.2f}%")

# ------------------------------------- multi-locus prediction, quantitative test
mr = pd.read_csv(f"{OUT}/multiread_disagreement.tsv", sep="\t", index_col=0)
# second, VASA-independent multi-locus measure: paralog family size from the GTF
base = pd.Series({g: re.sub(r"[-_.]?\d+$", "", src_nm[g].lower()) for g in src_bt})
famsize = base.value_counts()
gene_fam = pd.DataFrame({"biotype": pd.Series(src_bt), "fam": base})
gene_fam["famsize"] = gene_fam.fam.map(famsize)
fs = gene_fam.groupby("biotype").famsize.median()
mr["gtf_median_paralog_family_size"] = pd.Series({canon(k): v for k, v in fs.items()})
m = mr.dropna(subset=["log2_ratio", "vasa_read_weighted_n_genes_in_entry"])
m = m[(m.rsem_pooled_count > 0) & (m.featurecounts_pooled_count > 0)]
rho1, p1 = spearmanr(m.log2_ratio.abs(), m.vasa_read_weighted_n_genes_in_entry)
rho2, p2 = spearmanr(m.log2_ratio.abs(), m.gtf_median_paralog_family_size)
note(f"|log2(RSEM/featureCounts)| vs VASA read-weighted n_genes_in_entry: "
     f"Spearman rho={rho1:.3f} perm-p={p1:.4g} (seed {SEED}, 20k perms) over n={len(m)} biotypes")
note(f"|log2(RSEM/featureCounts)| vs GTF median paralog family size: "
     f"rho={rho2:.3f} perm-p={p2:.4g}")
for cl, sub in m.groupby("class"):
    note(f"  class={cl}: n={len(sub)} median |log2 ratio|={sub.log2_ratio.abs().median():.3f} "
         f"median |delta_pp|={sub.delta_pp.abs().median():.4f}")
mr.to_csv(f"{OUT}/multiread_disagreement.tsv", sep="\t")

with open(f"{OUT}/mechanisms_log.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")
print("\n".join(log))
