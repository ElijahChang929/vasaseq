#!/usr/bin/env python3
"""
05_nfcore_final.py -- assemble the nf-core cross-check deliverables.

Adds to scripts 02/04 the two things those runs showed were needed:

  1. FRAGMENT-level region counts over the rRNA-biotype loci. Script 04's
     samtools call failed silently (bare binary, libbz2 not on the loader path
     without `module load`), so the locus attribution was NaN. Redone here with
     -f 0x40 (first-in-pair) so the unit is fragments, matching featureCounts -p.

  2. An EXON-RESTRICTED VASA composition. featureCounts -t exon and RSEM's
     transcriptome cannot count intronic reads at all; the VASA route can and
     does. Comparing VASA-total against nf-core therefore charges an annotation-
     model difference to the protocol. The spliced-only VASA column is the fair
     analogue, and the gap between the two columns is that pipeline term,
     measured rather than assumed.

Unit note carried through every table: featureCounts counts assigned FRAGMENTS
(-p -B -C, multireads and biotype-ambiguous fragments discarded); RSEM reports
EM-distributed fractional expected_count per gene; the VASA route counts READS
(step 7 ReadCounts). Fragment vs read is a factor ~2 on a paired library and is
NOT corrected anywhere below -- every figure is a WITHIN-route fraction, where
the factor cancels.

Seeds: none used here (the permutation test lives in 04_).
"""
import gzip, os, re, subprocess, sys
import pandas as pd, numpy as np

W   = "/nemo/lab/turnerj/working/guangxin/vasaseq"
FSR = f"{W}/data/flashseq/results"
RAW = f"{W}/data/PM26037/out"
VA  = f"{RAW}/analysis"
SRC_GTF = "/nemo/lab/turnerj/working/guangxin/reference/gtf/Mus_musculus.GRCm39.116.gtf.gz"
OUT = os.environ.get("OUTDIR", ".")
LIBS = [f"ZHA8833A{i}" for i in range(1, 11)]
SC   = ["ZHA8833A9", "ZHA8833A10"]           # the 30 pg single-cell-equivalent rung
CELLS = ["%03d" % i for i in range(2, 14)]   # 12 real cells
RRNA_CANON = {"rrna", "mtrrna", "rrnapseudogene"}
log = []
def note(m): log.append(str(m)); print(m, file=sys.stderr, flush=True)
def canon(s): return s.replace("_", "").lower()

# ---------------------------------------------------------------- annotation
gid = re.compile(r'gene_id "([^"]+)"'); gbt = re.compile(r'gene_biotype "([^"]+)"')
gnm = re.compile(r'gene_name "([^"]+)"')
src_bt, src_nm, rr_loc = {}, {}, []
with gzip.open(SRC_GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene": continue
        g = gid.search(f[8]).group(1); b = gbt.search(f[8]).group(1)
        n = gnm.search(f[8]); n = n.group(1) if n else g
        src_bt[g] = b; src_nm[g] = n
        if b in ("rRNA", "Mt_rRNA"):
            rr_loc.append((g, n, b, f[0], int(f[3]), int(f[4])))
note(f"source GTF: {len(src_bt)} genes, {len(rr_loc)} rRNA/Mt_rRNA loci")

# ------------------------------- 1. fragment-level counts over the rRNA loci
ST = os.environ.get("SAMTOOLS", "samtools")
v = subprocess.run([ST, "--version"], capture_output=True, text=True)
assert v.returncode == 0 and v.stdout.strip(), f"samtools not usable: {v.stderr[:200]}"
note(f"samtools: {v.stdout.splitlines()[0]}")

rows = []
for lib in SC:
    bam = f"{FSR}/star_rsem/{lib}.markdup.sorted.bam"
    for g, n, b, ctg, s, e in rr_loc:
        L = e - s + 1
        if L < 400:              # the short 5S-class loci individually; pooled below
            continue
        q = lambda extra: int(subprocess.run(
            [ST, "view", "-c", "-f", "0x40", "-F", "0x100"] + extra + [bam, f"{ctg}:{s}-{e}"],
            capture_output=True, text=True, check=True).stdout.strip())
        rows.append(dict(library=lib, gene_name=n, gene_id=g, biotype=b, contig=ctg,
                         start=s, end=e, gene_len_bp=L,
                         unique_fragments=q(["-q", "255"]), all_primary_fragments=q([])))
    # the short (<400 bp) rRNA-biotype loci, pooled, so nothing is hidden
    short = [(ctg, s, e) for g, n, b, ctg, s, e in rr_loc
             if b == "rRNA" and (e - s + 1) < 400]
    regs = [f"{ctg}:{s}-{e}" for ctg, s, e in short]
    tot_u = 0
    for i in range(0, len(regs), 300):
        tot_u += int(subprocess.run([ST, "view", "-c", "-f", "0x40", "-F", "0x100",
                                     "-q", "255", bam] + regs[i:i + 300],
                                    capture_output=True, text=True, check=True).stdout.strip())
    rows.append(dict(library=lib, gene_name=f"[{len(short)} rRNA loci <400 bp, pooled]",
                     gene_id="-", biotype="rRNA", contig="-", start=-1, end=-1,
                     gene_len_bp=int(np.median([e - s + 1 for _, s, e in short])),
                     unique_fragments=tot_u, all_primary_fragments=np.nan))
loc = pd.DataFrame(rows)

# featureCounts' own rRNA total, and RSEM's, for the same libraries
fcr, rsr = {}, {}
for lib in SC:
    d = {}
    for line in open(f"{FSR}/star_rsem/featurecounts/{lib}.biotype_counts_mqc.tsv"):
        if line.startswith("#"): continue
        k, val = line.rstrip("\n").split("\t"); d[k] = float(val)
    fcr[lib] = d["rRNA"]
    g = pd.read_csv(f"{FSR}/star_rsem/{lib}.genes.results", sep="\t").set_index("gene_id")
    g["biotype"] = pd.Series(src_bt)
    rsr[lib] = float(g[g.biotype.eq("rRNA")].expected_count.sum())
loc["featurecounts_rRNA_biotype_total"] = loc.library.map(fcr)
loc["pct_of_featurecounts_rRNA_total"] = 100 * loc.unique_fragments / loc.library.map(fcr)
loc["rsem_rRNA_biotype_total"] = loc.library.map(rsr)
loc.to_csv(f"{OUT}/rrna_locus_attribution.tsv", sep="\t", index=False)
note("rRNA-biotype locus attribution (fragments, first-in-pair, primary):\n" +
     loc[["library", "gene_name", "gene_len_bp", "unique_fragments",
          "featurecounts_rRNA_biotype_total", "pct_of_featurecounts_rRNA_total"]]
     .to_string(index=False))

# -------------------- 2. VASA composition: total vs spliced-only (exon analogue)
def rd(p):
    d = pd.read_csv(p, sep="\t", index_col=0)
    d = d[~pd.isna(d.index)]; d.index = d.index.astype(str)
    return d.reindex(columns=CELLS).fillna(0.0)
u_tot = rd(f"{RAW}/ZHA9292A1_uniaggGenes_total.ReadCounts.tsv")
u_spl = rd(f"{RAW}/ZHA9292A1_uniaggGenes_spliced.ReadCounts.tsv").reindex(u_tot.index).fillna(0.0)
u_uns = rd(f"{RAW}/ZHA9292A1_uniaggGenes_unspliced.ReadCounts.tsv").reindex(u_tot.index).fillna(0.0)
ubt = pd.Series({e: e.split("_")[-1] for e in u_tot.index})

def compose(series_by_entry, biotypes, name):
    s = series_by_entry.groupby(biotypes.reindex(series_by_entry.index)).sum()
    s.index = [canon(i) for i in s.index]; s = s.groupby(level=0).sum()
    tot = s.sum(); rr = s[[i for i in s.index if i in RRNA_CANON]].sum(); non = tot - rr
    return pd.Series({i: (100 * v / tot if i in RRNA_CANON else 100 * v / non)
                      for i, v in s.items()}, name=name)

v_tot  = compose(u_tot.sum(axis=1), ubt, "vasa_uniagg_total_reads")
v_spl  = compose(u_spl.sum(axis=1), ubt, "vasa_uniagg_spliced_only")
v_uns  = compose(u_uns.sum(axis=1), ubt, "vasa_uniagg_unspliced_only")
note(f"uniagg pooled reads: total {u_tot.values.sum():,.0f}, "
     f"spliced {u_spl.values.sum():,.0f}, unspliced {u_uns.values.sum():,.0f} "
     f"(unspliced share {100*u_uns.values.sum()/u_tot.values.sum():.2f}%)")

# per-cell spread of the spliced-only composition, for a range not just a pool
percell = {}
for cell in CELLS:
    percell[cell] = compose(u_spl[cell], ubt, cell)
v_spl_pc = pd.DataFrame(percell)

# the analysis-set (multiagg, combination-entry) table, for the allocation-rule pair
a_tot = pd.read_csv(f"{VA}/ZHA9292A1_analysis_total.ReadCounts.tsv", sep="\t", index_col=0)
gm = pd.read_csv(f"{VA}/gene_metadata.tsv", sep="\t", index_col=0)
abt = gm["biotype"].astype(str)
a_sum = a_tot.sum(axis=1)
pure = ~abt.str.contains("-")
v_pure = compose(a_sum[pure], abt[pure], "vasa_analysis_pure_entries")

# ------------------------------------------------- nf-core sides, recomputed
fc = pd.DataFrame([{k: float(v) for k, v in
    (l.rstrip("\n").split("\t") for l in open(f"{FSR}/star_rsem/featurecounts/{lib}.biotype_counts_mqc.tsv")
     if not l.startswith("#"))} for lib in LIBS], index=LIBS).fillna(0.0)
rs = pd.read_csv(f"{FSR}/star_rsem/rsem.merged.gene_counts.tsv", sep="\t").set_index("gene_id")
rs = rs.drop(columns=["transcript_id(s)"])[LIBS]
rs_bt = pd.Series({g: src_bt.get(g, "UNANNOTATED") for g in rs.index})
rsb = rs.groupby(rs_bt.values).sum().T

def compose_wide(df, name, libs=None):
    d = df.loc[libs] if libs else df
    s = d.sum(axis=0); s.index = [canon(i) for i in s.index]; s = s.groupby(level=0).sum()
    tot = s.sum(); rr = s[[i for i in s.index if i in RRNA_CANON]].sum(); non = tot - rr
    return pd.Series({i: (100 * v / tot if i in RRNA_CANON else 100 * v / non)
                      for i, v in s.items()}, name=name)

f_rsem_sc = compose_wide(rsb, "nfcore_rsem_30pg", SC)
f_fc_sc   = compose_wide(fc,  "nfcore_featurecounts_30pg", SC)
f_rsem_al = compose_wide(rsb, "nfcore_rsem_all10")
f_fc_al   = compose_wide(fc,  "nfcore_featurecounts_all10")

# ============================================================== deliverable
pvp = pd.concat([f_rsem_al, f_rsem_sc, f_fc_al, f_fc_sc,
                 v_tot, v_spl, v_uns, v_pure], axis=1)
pvp.index.name = "biotype_canon"
disp = {}
for b in set(src_bt.values()) | set(abt.unique()) | set(ubt.unique()):
    for p in str(b).split("-"): disp.setdefault(canon(p), p)
pvp.insert(0, "display_biotype", [disp.get(i, i) for i in pvp.index])
pvp["denominator"] = np.where(pvp.index.isin(RRNA_CANON), "all reads", "non-rRNA reads")

# --- the three measured terms -------------------------------------------------
# A. QUANTIFIER, identical reads: RSEM vs featureCounts on the same FLASH-seq BAMs
pvp["A_quantifier_pp"]    = pvp.nfcore_rsem_30pg - pvp.nfcore_featurecounts_30pg
pvp["A_quantifier_log2"]  = np.log2(pvp.nfcore_rsem_30pg / pvp.nfcore_featurecounts_30pg
                                    .replace(0, np.nan))
# B. ANNOTATION MODEL, identical reads (VASA side): what the intron rows add
pvp["B_intron_model_pp"]  = pvp.vasa_uniagg_total_reads - pvp.vasa_uniagg_spliced_only
# C. the confounded FLASH-seq(nf-core) vs VASA(VASA-route) total gap
pvp["C_total_gap_pp"]     = pvp.nfcore_rsem_30pg - pvp.vasa_uniagg_total_reads
# residual after removing the two measurable pipeline terms
pvp["C_minus_A_minus_B_pp"] = pvp.C_total_gap_pp - pvp.A_quantifier_pp.abs() - pvp.B_intron_model_pp
pvp["exon_matched_gap_pp"] = pvp.nfcore_rsem_30pg - pvp.vasa_uniagg_spliced_only
pvp["vasa_route_flashseq_native"] = np.nan   # other track; deliberately not invented

def cls(i):
    if i in RRNA_CANON: return "rRNA (annotation route -- see rrna_biotype_note)"
    if i in {"snorna","snrna","mirna","scarna","srna","scrna","miscrna","trna",
             "mttrna","ribozyme"}: return "short non-poly-A"
    if "pseudogene" in i: return "pseudogene"
    if i in {"proteincoding","lncrna","tec"}: return "long poly-A / lncRNA"
    return "other"
pvp["class"] = [cls(i) for i in pvp.index]
pvp["decomposes_cleanly"] = np.where(
    pvp[["nfcore_rsem_30pg", "vasa_uniagg_spliced_only"]].notna().all(axis=1),
    "partial: A and B measured, protocol term needs vasa_route_flashseq_native (pending)",
    "no: class absent on one side")
pvp = pvp.sort_values("nfcore_rsem_30pg", ascending=False)
pvp.to_csv(f"{OUT}/pipeline_vs_protocol.tsv", sep="\t")
note(f"wrote pipeline_vs_protocol.tsv {pvp.shape}")

KEY = ["proteincoding","lncrna","rrna","mtrrna","processedpseudogene",
       "transcribedprocessedpseudogene","unprocessedpseudogene","miscrna",
       "snrna","snorna","mirna","scarna","ribozyme","mttrna","tec"]
show = [k for k in KEY if k in pvp.index]
note("\n=== pipeline_vs_protocol, key classes ===\n" + pvp.loc[show, [
    "display_biotype","denominator","nfcore_rsem_30pg","nfcore_featurecounts_30pg",
    "vasa_uniagg_total_reads","vasa_uniagg_spliced_only","A_quantifier_pp",
    "B_intron_model_pp","C_total_gap_pp","exon_matched_gap_pp"]].round(4).to_string())

# per-cell spread on the VASA spliced-only side
spread = pd.DataFrame({"min": v_spl_pc.min(axis=1), "median": v_spl_pc.median(axis=1),
                       "max": v_spl_pc.max(axis=1), "pooled": v_spl})
spread.to_csv(f"{OUT}/vasa_spliced_percell_spread.tsv", sep="\t")
note("VASA spliced-only per-cell spread (12 cells), key classes:\n" +
     spread.reindex(show).round(4).to_string())

with open(f"{OUT}/final_log.txt", "w") as fh: fh.write("\n".join(log) + "\n")
