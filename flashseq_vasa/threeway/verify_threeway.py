#!/usr/bin/env python3
"""Re-derive and ASSERT every number quoted in THREEWAY_COMPOSITION.md.

Project rule: verify interactive numbers with a committed script before they
enter a report. Nothing here is transcribed -- every constant is asserted
against the TSV that produced it, and the two cross-track constants
(annotation_crossplate_gap.tsv) are asserted against the OTHER track's output,
so a change in either track breaks this script rather than silently
invalidating the report.

Usage: verify_threeway.py <RESDIR>   -- exits non-zero on any mismatch.
"""
import sys
import numpy as np
import pandas as pd

R = sys.argv[1]
STRUCTURAL = ["MiscRna", "snRNA", "snoRNA", "scaRNA", "ribozyme"]
ok, fail = [], []


def chk(name, got, want, tol=5e-4):
    good = (np.isnan(got) and np.isnan(want)) if isinstance(want, float) and np.isnan(want) \
        else abs(got - want) <= tol
    (ok if good else fail).append(f"{'PASS' if good else 'FAIL'} {name}: got {got!r} want {want!r}")


comp = pd.read_csv(f"{R}/composition_threeway.tsv", sep="\t", index_col=0)
den = pd.read_csv(f"{R}/threeway_denominators.tsv", sep="\t", index_col=0)
st = pd.read_csv(f"{R}/threeway_structural.tsv", sep="\t", index_col=0)
al = pd.read_csv(f"{R}/threeway_allocation_rules.tsv", sep="\t", index_col=0)
rvm = pd.read_csv(f"{R}/threeway_reads_vs_molecules.tsv", sep="\t", index_col=0)
geo = pd.read_csv(f"{R}/threeway_readlength_geometry.tsv", sep="\t")
pu = pd.read_csv(f"{R}/threeway_structural_per_unit.tsv", sep="\t")
rel = pd.read_csv(f"{R}/threeway_release_control.tsv", sep="\t", index_col=0)
disc = pd.read_csv(f"{R}/threeway_discordant_combos.tsv", sep="\t")
cpg = pd.read_csv(f"{R}/annotation_crossplate_gap.tsv", sep="\t")
churn = pd.read_csv(f"{R}/annotation_biotype_churn_ledger.tsv", sep="\t", index_col=0)

# --- 1. headline composition -------------------------------------------------
s = st["structural_pct_of_nonrRNA"]
chk("structural published", s["published_vasa"], 2.2474)
chk("structural own", s["own_vasa"], 20.5906)
chk("structural FS native", s["flashseq_native"], 0.2592)
chk("structural FS vasalen", s["flashseq_vasalen"], 0.2585)
chk("fold own/FSnative", s["own_vasa"] / s["flashseq_native"], 79.44, 0.02)
chk("fold published/FSnative", s["published_vasa"] / s["flashseq_native"], 8.67, 0.02)
chk("fold own/published", s["own_vasa"] / s["published_vasa"], 9.16, 0.01)
chk("readlen effect pp", s["flashseq_vasalen"] - s["flashseq_native"], -0.0007, 1e-4)

# structural total must equal the sum of its five members in composition_threeway
for k in ("published_vasa", "own_vasa", "flashseq_native", "flashseq_vasalen"):
    chk(f"structural sum consistency {k}", float(comp.loc[STRUCTURAL, k].sum()), float(s[k]), 1e-3)

# --- 2. denominators and units ----------------------------------------------
chk("n_units published", float(den.loc["published_vasa", "n_units"]), 173.0, 0)
chk("n_units own", float(den.loc["own_vasa", "n_units"]), 12.0, 0)
chk("n_units FS", float(den.loc["flashseq_native", "n_units"]), 10.0, 0)
chk("denom published", float(den.loc["published_vasa", "denominator_nonrRNA_reads"]), 42050713.0, 1)
chk("denom own", float(den.loc["own_vasa", "denominator_nonrRNA_reads"]), 52823542.0, 1)
chk("denom FS native", float(den.loc["flashseq_native", "denominator_nonrRNA_reads"]), 238184504.0, 1)
chk("species-filtered reads published", float(den.loc["published_vasa", "reads_dropped_species_filter"]),
    1700818.0, 1)
chk("pct combos own", float(den.loc["own_vasa", "pct_reads_in_combination_keys"]), 20.5741, 1e-3)
chk("pct combos published", float(den.loc["published_vasa", "pct_reads_in_combination_keys"]), 3.7031, 1e-3)

# --- 3. ProteinCoding / lncRNA -------------------------------------------
chk("PC published", float(comp.loc["ProteinCoding", "published_vasa"]), 92.2232)
chk("PC own", float(comp.loc["ProteinCoding", "own_vasa"]), 64.1192)
chk("PC FS native", float(comp.loc["ProteinCoding", "flashseq_native"]), 89.6426)
chk("lncRNA own", float(comp.loc["lncRNA", "own_vasa"]), 13.9149)
chk("lncRNA published", float(comp.loc["lncRNA", "published_vasa"]), 2.4443)
chk("snRNA own", float(comp.loc["snRNA", "own_vasa"]), 7.0327)
chk("snRNA published", float(comp.loc["snRNA", "published_vasa"]), 0.1625)
chk("snRNA FS native", float(comp.loc["snRNA", "flashseq_native"]), 0.000364, 1e-5)

# --- 4. allocation invariance -----------------------------------------------
cols = ["structural_pct_last", "structural_pct_unan", "structural_pct_frac"]
chk("allocation rules identical (max spread)",
    float((al[cols].max(axis=1) - al[cols].min(axis=1)).max()), 0.0, 1e-9)
chk("discordant combos w/ structural, published",
    float(disc.set_index("dataset").loc["published_vasa", "discordant_reads_containing_structural"]), 0.0, 0)
chk("discordant combos w/ structural, own",
    float(disc.set_index("dataset").loc["own_vasa", "discordant_reads_containing_structural"]), 0.0, 0)

# --- 5. duplication ---------------------------------------------------------
sr = {k: float(rvm.loc[STRUCTURAL, f"{k}_pct_reads"].sum()) for k in ("published_vasa", "own_vasa")}
sm = {k: float(rvm.loc[STRUCTURAL, f"{k}_pct_molecules"].sum()) for k in ("published_vasa", "own_vasa")}
chk("structural gap reads pp", sr["own_vasa"] - sr["published_vasa"], 18.3432, 1e-3)
chk("structural gap molecules pp", sm["own_vasa"] - sm["published_vasa"], 16.7988, 1e-3)
chk("fold in molecule space", sm["own_vasa"] / sm["published_vasa"], 8.599, 0.01)
chk("dedup moves gap pp", (sr["own_vasa"] - sr["published_vasa"]) - (sm["own_vasa"] - sm["published_vasa"]),
    1.5444, 1e-3)

# --- 6. geometry ------------------------------------------------------------
g = geo.groupby("dataset")[["span_median", "pct_jS_IN_structural", "pct_rows_structural"]].median()
chk("median span own", float(g.loc["own_vasa", "span_median"]), 130.0, 0)
chk("median span published", float(g.loc["published_vasa", "span_median"]), 74.0, 0)
chk("jS:IN structural own", float(g.loc["own_vasa", "pct_jS_IN_structural"]), 82.972, 1e-2)
chk("jS:IN structural published", float(g.loc["published_vasa", "pct_jS_IN_structural"]), 82.165, 1e-2)
chk("BED structural row share own", float(g.loc["own_vasa", "pct_rows_structural"]), 4.835, 1e-2)
chk("BED structural row share published", float(g.loc["published_vasa", "pct_rows_structural"]), 0.515, 1e-2)
chk("BED-level fold own/published",
    float(g.loc["own_vasa", "pct_rows_structural"] / g.loc["published_vasa", "pct_rows_structural"]),
    9.384, 0.01)
# the direction claim: own reads must be LONGER, so geometry works against the gap
assert g.loc["own_vasa", "span_median"] > g.loc["published_vasa", "span_median"], \
    "direction claim broken: own reads are NOT longer than published"
ok.append("PASS own reads longer than published (geometry opposes the gap)")

# --- 7. pooling / overlap ---------------------------------------------------
lo = float(pu[pu.dataset == "own_vasa"].structural_pct.min())
hi = float(pu[pu.dataset == "published_vasa"].structural_pct.max())
chk("own-plate minimum", lo, 13.8326, 1e-3)
chk("published maximum", hi, 11.6269, 1e-3)
chk("published cells at/above own minimum",
    float((pu[pu.dataset == "published_vasa"].structural_pct >= lo).sum()), 0.0, 0)
chk("n published cells", float((pu.dataset == "published_vasa").sum()), 173.0, 0)

# --- 8. release control, this track ----------------------------------------
pub_e99 = float(rel.loc[STRUCTURAL, "pub_shared_E99"].sum())
pub_e116 = float(rel.loc[STRUCTURAL, "pub_shared_E116"].sum())
own_e116 = float(rel.loc[STRUCTURAL, "own_shared_E116"].sum())
chk("published structural on shared rows, E99 labels", pub_e99, 0.5850, 1e-3)
chk("published structural on shared rows, E116 labels", pub_e116, 0.5850, 1e-3)
chk("relabel effect on structural (pp)", abs(pub_e116 - pub_e99), 0.0, 1e-6)
chk("own structural on shared rows", own_e116, 3.4726, 1e-3)
chk("restricted fold own/published", own_e116 / pub_e116, 5.936, 0.01)
for b in STRUCTURAL:
    chk(f"churn: {b} lost to other class", float(churn.loc[b, "n_lost_to_other_class"]), 0.0, 0)
    chk(f"churn: {b} gained from other class", float(churn.loc[b, "n_gained_from_other_class"]), 0.0, 0)

# --- 9. cross-track: the annotation control's own number -------------------
snc = cpg[(cpg.level == "coarse") & (cpg["class"] == "sncRNA")].iloc[0]
chk("annot track sncRNA published", float(snc.published_all), 0.9666, 1e-3)
chk("annot track sncRNA own", float(snc.own_all), 6.3704, 1e-3)
chk("annot track sncRNA gap all pp", float(snc.gap_all_pp), 5.4037, 1e-3)
chk("annot track sncRNA gap shared pp", float(snc.gap_shared_pp), 2.7515, 1e-3)
chk("annot track gap closed pp", float(snc.gap_closed_pp), 2.6522, 1e-3)
chk("annot track pct of gap closed", 100 * float(snc.gap_closed_pp) / float(snc.gap_all_pp), 49.08, 0.05)
chk("annot track raw fold", float(snc.own_all / snc.published_all), 6.59, 0.01)
chk("annot track restricted fold", float(snc.own_shared / snc.published_shared), 5.51, 0.01)

# --- 10. both tracks must agree the residual is 5-6x ----------------------
r1 = float(snc.own_shared / snc.published_shared)
r2 = own_e116 / pub_e116
assert 5.0 <= r1 <= 6.0 and 5.0 <= r2 <= 6.0, f"residual folds not both in [5,6]: {r1}, {r2}"
ok.append(f"PASS both tracks residual in 5-6x ({r1:.2f}x, {r2:.2f}x)")

print("\n".join(ok))
if fail:
    print("\n".join(fail))
    sys.exit(f"{len(fail)} of {len(ok) + len(fail)} checks FAILED")
print(f"\nall {len(ok)} checks passed")
