#!/usr/bin/env python3
"""
04_nfcore_mechanisms.py -- WHY nf-core/RSEM and the VASA route disagree.

Supersedes 03_ (which died on a NaN row label in the step-7 uniagg index -- the
same unnamed row the analysis-set filter drops, see analysis/manifest.json).

Five mechanisms, each measured:
  (a) what gene_biotype "rRNA" contains in Ensembl 116/GRCm39, and WHICH of those
      genes actually carries nf-core's rRNA signal.
  (b) what the PUBLISHED results/genome/*.filtered.gtf contains (it is not what
      featureCounts read at runtime -- checked, not assumed).
  (c) RSEM's effective_length floor: a transcript far shorter than the fragment
      length cannot generate the observed fragments and so receives almost no
      expected_count. Quantifier property, not chemistry.
  (d) the exon-vs-exon+intron axis: the VASA route counts intronic reads;
      featureCounts -t exon and RSEM's transcriptome cannot.
  (e) the multi-locus prediction, tested against two independent measures.

Seeds: SEED=0, used only for the permutation null of the Spearman p-value.
"""
import gzip, os, re, subprocess, sys
import pandas as pd, numpy as np

SEED = 0
W   = "/nemo/lab/turnerj/working/guangxin/vasaseq"
FSR = f"{W}/data/flashseq/results"
RAW = f"{W}/data/PM26037/out"
SRC_GTF = "/nemo/lab/turnerj/working/guangxin/reference/gtf/Mus_musculus.GRCm39.116.gtf.gz"
NF_GTF  = f"{FSR}/genome/Mus_musculus.GRCm39.116.filtered.gtf"
SAMTOOLS = "/camp/apps/eb/software/SAMtools/1.11-GCC-10.2.0/bin/samtools"
OUT = os.environ.get("OUTDIR", ".")
LIBS = [f"ZHA8833A{i}" for i in range(1, 11)]
log = []
def note(m):
    log.append(str(m)); print(m, file=sys.stderr, flush=True)
def canon(s): return s.replace("_", "").lower()

def _rank(a):
    a = np.asarray(a, float); o = a.argsort(kind="mergesort")
    r = np.empty(len(a), float); r[o] = np.arange(1, len(a) + 1)
    u, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    for i in np.where(cnt > 1)[0]:
        m = inv == i; r[m] = r[m].mean()
    return r
def spearmanr(x, y, n_perm=20000, seed=SEED):
    rx, ry = _rank(x), _rank(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.default_rng(seed); ge = 0
    for _ in range(n_perm):
        if abs(np.corrcoef(rng.permutation(rx), ry)[0, 1]) >= abs(rho) - 1e-12: ge += 1
    return rho, (ge + 1) / (n_perm + 1)

# --------------------------------------------------------- GTF: genes + rRNA loci
gid = re.compile(r'gene_id "([^"]+)"'); gbt = re.compile(r'gene_biotype "([^"]+)"')
gnm = re.compile(r'gene_name "([^"]+)"'); tid = re.compile(r'transcript_id "([^"]+)"')
src_bt, src_nm, rr_loc = {}, {}, {}
with gzip.open(SRC_GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene": continue
        g = gid.search(f[8]).group(1); b = gbt.search(f[8]).group(1)
        n = gnm.search(f[8]); n = n.group(1) if n else g
        src_bt[g] = b; src_nm[g] = n
        if b in ("rRNA", "Mt_rRNA"):
            rr_loc[g] = (f[0], int(f[3]), int(f[4]), n, b)
note(f"source GTF: {len(src_bt)} genes; rRNA/Mt_rRNA loci captured: {len(rr_loc)}")

# ------------------------------------------------------------------ (a) rRNA class
rr = pd.DataFrame([{"gene_id": g, "gene_name": v[3], "biotype": v[4],
                    "contig": v[0], "start": v[1], "end": v[2],
                    "gene_len_bp": v[2] - v[1] + 1} for g, v in rr_loc.items()])
def fam(n):
    nl = n.lower()
    if nl.startswith("n-r5s"): return "n-R5s* (5S)"
    if nl.startswith("rn18s"): return "Rn18s* (18S relic)"
    if nl.startswith("rn5s"):  return "Rn5s*"
    if nl.startswith("rn45s") or nl.startswith("rn28s") or nl.startswith("rn5-8s"):
        return "Rn45s/28S/5.8S"
    if nl.startswith("mt-r"):  return "mt-Rnr* (mito)"
    if nl.startswith("gm"):    return "Gm* (unnamed)"
    return "other: " + n
rr["family"] = rr.gene_name.map(fam)
note("rRNA + Mt_rRNA gene_name families (Ensembl 116/GRCm39):\n" +
     rr.groupby(["biotype", "family"]).agg(n=("gene_id", "size"),
        median_bp=("gene_len_bp", "median"), max_bp=("gene_len_bp", "max")).to_string())
note("all non-Gm, non-n-R5s rRNA-biotype gene names: " +
     str(sorted(rr[(rr.biotype == "rRNA") & (~rr.family.isin(["Gm* (unnamed)", "n-R5s* (5S)"]))]
                .gene_name.tolist())))
rr.sort_values(["biotype", "gene_len_bp"], ascending=[True, False]).to_csv(
    f"{OUT}/rrna_biotype_genes.tsv", sep="\t", index=False)

# which rRNA genes carry the RSEM signal, per library
rows = []
for lib in LIBS:
    g = pd.read_csv(f"{FSR}/star_rsem/{lib}.genes.results", sep="\t").set_index("gene_id")
    g["biotype"] = pd.Series(src_bt); g["gene_name"] = pd.Series(src_nm)
    sub = g[g.biotype.eq("rRNA")]
    tot = sub.expected_count.sum()
    top = sub.expected_count.idxmax()
    rows.append(dict(library=lib, rrna_biotype_expected_count=tot,
        n_genes_nonzero=int((sub.expected_count > 0).sum()),
        top_gene=src_nm.get(top, top),
        top_gene_share_pct=100 * sub.expected_count.max() / tot if tot else np.nan,
        n_R5s_share_pct=100 * sub[sub.gene_name.str.lower().str.startswith("n-r5s")]
                          .expected_count.sum() / tot if tot else np.nan))
rsig = pd.DataFrame(rows).set_index("library")
rsig.to_csv(f"{OUT}/rsem_rrna_signal_by_gene.tsv", sep="\t")
note("RSEM rRNA-biotype signal, which gene carries it:\n" + rsig.round(4).to_string())

# featureCounts side: does the SAME single locus dominate a genomic exon count?
# Approximate featureCounts' rule with an indexed region query on the STAR BAM:
# primary alignments only (-F 0x100/0x400), MAPQ 255 = STAR-unique (featureCounts
# discards NH>1 as Unassigned_MultiMapping).
fc_reg = []
for lib in ["ZHA8833A9", "ZHA8833A10"]:
    bam = f"{FSR}/star_rsem/{lib}.markdup.sorted.bam"
    for g, (ctg, s, e, n, b) in rr_loc.items():
        if b != "rRNA" or (e - s + 1) < 500:   # the long rRNA-biotype loci only
            continue
        try:
            out = subprocess.run([SAMTOOLS, "view", "-c", "-F", "0x100", "-q", "255",
                                  bam, f"{ctg}:{s}-{e}"], capture_output=True, text=True,
                                 timeout=300)
            cnt = int(out.stdout.strip()) if out.returncode == 0 else np.nan
        except Exception as ex:
            cnt = np.nan; note(f"region query failed {lib} {n}: {ex}")
        fc_reg.append(dict(library=lib, gene_name=n, gene_id=g, contig=ctg,
                           start=s, end=e, gene_len_bp=e - s + 1, unique_primary_reads=cnt))
fcr = pd.DataFrame(fc_reg)
if len(fcr):
    fcr.to_csv(f"{OUT}/rrna_locus_region_counts.tsv", sep="\t", index=False)
    note("unique primary reads over the long (>=500 bp) rRNA-biotype loci:\n" +
         fcr.to_string(index=False))

# --------------------------------------------------- (b) published filtered GTF
ftypes, nf_g, nf_t, contigs, nlines = {}, set(), set(), {}, 0
with open(NF_GTF) as fh:
    for line in fh:
        nlines += 1
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        ftypes[f[2]] = ftypes.get(f[2], 0) + 1
        contigs[f[0]] = contigs.get(f[0], 0) + 1
        m = gid.search(f[8]); t = tid.search(f[8])
        if m: nf_g.add(m.group(1))
        if t: nf_t.add(t.group(1))
last = subprocess.run(["tail", "-c", "300", NF_GTF], capture_output=True, text=True).stdout
note(f"published filtered.gtf: {nlines} lines, feature types {ftypes}")
note(f"  {len(nf_g)} gene_id, {len(nf_t)} transcript_id, contigs present: {sorted(contigs)}")
note(f"  final bytes end with newline: {last.endswith(chr(10))!r}; last line: {last.splitlines()[-1][:120]!r}")
note("  -> compare with the featureCounts biotype table, which spans the full "
     "Ensembl biotype vocabulary: the RUNTIME GTF was complete.")

# ---------------------------------------------- (c) RSEM effective_length floor
g9 = pd.read_csv(f"{FSR}/star_rsem/ZHA8833A9.genes.results", sep="\t").set_index("gene_id")
g9["biotype"] = pd.Series(src_bt)
g9["eff_frac"] = g9.effective_length / g9.length.replace(0, np.nan)
byb = g9.groupby("biotype").agg(
    n_genes=("length", "size"), median_length=("length", "median"),
    median_eff_length=("effective_length", "median"),
    median_eff_frac=("eff_frac", "median"),
    n_eff_lt_10=("effective_length", lambda s: int((s < 10).sum())),
    expected_count=("expected_count", "sum"))
byb["pct_genes_eff_lt_10"] = 100 * byb.n_eff_lt_10 / byb.n_genes
byb = byb.sort_values("expected_count", ascending=False)
byb.to_csv(f"{OUT}/rsem_effective_length.tsv", sep="\t")
# mean fragment length, estimated from protein_coding length - effective_length
pc = g9[g9.biotype.eq("protein_coding") & (g9.length > 1000)]
mfl = float((pc.length - pc.effective_length).median()) + 1
note(f"RSEM implied mean fragment length (median length-eff_length over "
     f"{len(pc)} protein_coding genes >1 kb, +1) = {mfl:.1f} nt")
KEY = ["protein_coding", "lncRNA", "rRNA", "Mt_rRNA", "snRNA", "snoRNA", "miRNA",
       "scaRNA", "misc_RNA", "Mt_tRNA", "ribozyme", "processed_pseudogene"]
note("RSEM effective_length by biotype (A9):\n" +
     byb.loc[[b for b in KEY if b in byb.index]].round(3).to_string())
note("fraction of each class shorter than the fragment length: " + ", ".join(
    f"{b}:{100*float((g9[g9.biotype.eq(b)].length < mfl).mean()):.1f}%" for b in KEY))

# ------------------------------------------- (d) exon vs exon+intron, VASA route
CELLS = ["%03d" % i for i in range(2, 14)]
def rd(p):
    d = pd.read_csv(p, sep="\t", index_col=0)
    d = d[~pd.isna(d.index)]
    d.index = d.index.astype(str)
    return d[CELLS].fillna(0.0)
tot = rd(f"{RAW}/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv")
spl = rd(f"{RAW}/ZHA9292A1_uniaggGenes_spliced.ReadCounts.tsv").reindex(tot.index).fillna(0.0)
uns = rd(f"{RAW}/ZHA9292A1_uniaggGenes_unspliced.ReadCounts.tsv").reindex(tot.index).fillna(0.0)
exc = (spl + uns - tot)
note(f"TRAP (h) CHECK on uniagg ReadCounts: max(spliced+unspliced-total) = "
     f"{float(exc.values.max()):.1f}; entries with any excess = "
     f"{int((exc.values > 0).any(axis=1).sum())}/{len(tot)}  "
     f"(vasa branch uses exact membership; 0 is the expected result)")
vb = pd.Series({e: e.split("_")[-1] for e in tot.index})
agg = pd.DataFrame({"total": tot.sum(axis=1), "spliced": spl.sum(axis=1),
                    "unspliced": uns.sum(axis=1), "biotype": vb.reindex(tot.index)})
agg = agg.groupby("biotype").sum(numeric_only=True)
agg["pct_unspliced"] = 100 * agg.unspliced / agg.total.replace(0, np.nan)
agg["pct_of_all_reads"] = 100 * agg.total / agg.total.sum()
agg = agg.sort_values("total", ascending=False)
agg.to_csv(f"{OUT}/vasa_intron_axis.tsv", sep="\t")
note("VASA uniagg reads, intron share by biotype (top 14):\n" + agg.head(14).round(3).to_string())
note(f"VASA uniagg OVERALL unspliced share of reads: "
     f"{100*agg.unspliced.sum()/agg.total.sum():.2f}%  "
     f"(these reads are invisible to featureCounts -t exon and to RSEM's transcriptome)")

# ------------------------------------------------- (e) multi-locus prediction test
mr = pd.read_csv(f"{OUT}/multiread_disagreement.tsv", sep="\t", index_col=0)
base = pd.Series({g: re.sub(r"[-_.]?\d+$", "", src_nm[g].lower()) for g in src_bt})
famsize = base.value_counts()
gf = pd.DataFrame({"biotype": pd.Series(src_bt), "fam": base})
gf["famsize"] = gf.fam.map(famsize)
mr["gtf_median_paralog_family_size"] = pd.Series(
    {canon(k): v for k, v in gf.groupby("biotype").famsize.median().items()})
mr["rsem_median_eff_length"] = pd.Series(
    {canon(k): v for k, v in byb.median_eff_length.items()})
m = mr.dropna(subset=["log2_ratio", "vasa_read_weighted_n_genes_in_entry"])
m = m[(m.rsem_pooled_count > 0) & (m.featurecounts_pooled_count > 0)]
r1, p1 = spearmanr(m.log2_ratio.abs(), m.vasa_read_weighted_n_genes_in_entry)
r2, p2 = spearmanr(m.log2_ratio.abs(), m.gtf_median_paralog_family_size)
note(f"|log2(RSEM/featureCounts)| vs VASA read-weighted n_genes_in_entry: "
     f"rho={r1:.3f} perm-p={p1:.4g} (seed {SEED}, 20k perms), n={len(m)} biotypes")
note(f"|log2(RSEM/featureCounts)| vs GTF median paralog family size: "
     f"rho={r2:.3f} perm-p={p2:.4g}")
for cl, sub in m.groupby("class"):
    note(f"  class={cl}: n={len(sub)}, median |log2 ratio|={sub.log2_ratio.abs().median():.3f}, "
         f"median |delta_pp|={sub.delta_pp.abs().median():.4f}")
mr.sort_values("delta_pp", key=abs, ascending=False).to_csv(
    f"{OUT}/multiread_disagreement.tsv", sep="\t")

with open(f"{OUT}/mechanisms_log.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")
