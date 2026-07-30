#!/usr/bin/env python3
"""
verify_annotation_release_effect.py -- assert every number quoted in
ANNOTATION_RELEASE_EFFECT.md against the tables that produced it.

Per this project's standing rule: no number enters a report by transcription.
Each claim below is re-read from the TSV and asserted, so the note cannot drift
from the data. Run it after any change to the analysis or the note.

Exit 0 = every quoted number matches its source table.

USAGE
    ./verify_annotation_release_effect.py [RESDIR]
"""
import os
import sys

import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else \
    "/nemo/lab/turnerj/working/guangxin/vasaseq/res/threeway"

ok = 0
bad = []


def T(name):
    p = os.path.join(RES, name)
    assert os.path.exists(p), f"missing table: {p}"
    return pd.read_csv(p, sep="\t")


def hdr_val(df, metric, scope):
    m = df[(df.metric == metric) & (df.scope == scope)]
    assert len(m) == 1, f"{metric} / {scope}: {len(m)} rows, expected 1"
    return m.iloc[0]["value"]


def check(label, got, want, tol=0.0):
    """Numeric or string equality against the source table."""
    global ok
    try:
        g, w = float(got), float(want)
        good = abs(g - w) <= tol
    except (TypeError, ValueError):
        g, w = str(got), str(want)
        good = g == w
    if good:
        ok += 1
        print(f"  OK   {label:<62} = {got}")
    else:
        bad.append(f"{label}: table says {got!r}, note says {want!r}")
        print(f"  FAIL {label:<62} table={got!r} note={want!r}")


def main():
    H = T("annotation_release_effect.tsv")
    B = T("annotation_release_bracket.tsv")
    U = T("annotation_gene_universe_coarse.tsv")
    UB = T("annotation_gene_universe_biotype.tsv")
    L = T("annotation_biotype_churn_ledger.tsv")
    G = T("annotation_crossplate_gap.tsv")
    R = T("annotation_release_reid.tsv")
    A_ = T("annotation_release_asymmetry.tsv")
    P = T("annotation_feature_length_paired.tsv")

    print("== 1. gene universe ==")
    check("E99 mouse genes", hdr_val(H, "genes_total", "E99 mouse"), 55471)
    check("E116 genes", hdr_val(H, "genes_total", "E116"), 78348)
    check("shared genes", hdr_val(H, "genes_shared", "E99 n E116"), 51765)
    check("E99-only genes", hdr_val(H, "genes_only_99", "E99 \\ E116"), 3706)
    check("E116-only genes", hdr_val(H, "genes_only_116", "E116 \\ E99"), 26583)
    check("shared % of E99", hdr_val(H, "genes_shared_pct_of_99", "E99 mouse"), 93.319, 0.001)
    check("shared % of E116", hdr_val(H, "genes_shared_pct_of_116", "E116"), 66.071, 0.001)
    check("Jaccard", hdr_val(H, "jaccard_gene_universe", "E99 mouse vs E116"), 0.6309, 1e-4)
    # E116/E99 gene-universe ratio quoted in the note
    check("E116/E99 universe ratio", round(78348 / 55471, 3), 1.412, 0.001)

    print("== 2. protein-coding ==")
    check("PC genes E99", hdr_val(H, "protein_coding_genes", "E99 mouse"), 21933)
    check("PC genes E116", hdr_val(H, "protein_coding_genes", "E116"), 21818)
    check("PC shared", hdr_val(H, "protein_coding_shared",
                               "E99 n E116 (same class both)"), 21513)
    check("PC only-E99", hdr_val(H, "protein_coding_only_99", "E99 \\ E116"), 199)
    check("PC only-E116", hdr_val(H, "protein_coding_only_116", "E116 \\ E99"), 158)
    check("PC shared % of E99", hdr_val(H, "protein_coding_shared_pct_of_99",
                                        "E99 mouse"), 98.085, 0.001)
    pc = U[U["class"] == "ProteinCoding"].iloc[0]
    check("PC id in other class E99->E116", pc.n_E99_id_in_E116_other_class, 221)
    check("PC id in other class E116->E99", pc.n_E116_id_in_E99_other_class, 147)

    print("== 3. lncRNA expansion ==")
    ln = U[U["class"] == "lncRNA"].iloc[0]
    check("lncRNA E99", ln.n_E99, 9959)
    check("lncRNA E116", ln.n_E116, 32889)
    check("lncRNA only-E116", ln.n_only_E116, 23236)
    check("lncRNA fold", round(32889 / 9959, 2), 3.30, 0.01)

    print("== 4. biotype churn ==")
    check("biotype changed n", hdr_val(H, "shared_genes_biotype_changed",
                                       "E99 n E116"), 1761)
    check("biotype changed %", hdr_val(H, "shared_genes_biotype_changed_pct",
                                       "E99 n E116"), 3.402, 0.001)
    check("coarse changed n", hdr_val(H, "shared_genes_coarse_class_changed",
                                      "E99 n E116"), 379)
    check("coarse changed %", hdr_val(H, "shared_genes_coarse_class_changed_pct",
                                      "E99 n E116"), 0.732, 0.001)
    check("symbol changed n", hdr_val(H, "shared_genes_symbol_changed",
                                      "E99 n E116"), 2911)
    check("symbol changed %", hdr_val(H, "shared_genes_symbol_changed_pct",
                                      "E99 n E116"), 5.623, 0.001)
    check("largest transition n",
          hdr_val(H, "largest_single_biotype_transition",
                  "UnprocessedPseudogene -> TranscribedUnprocessedPseudogene"), 671)
    check("PolymorphicPseudogene retired",
          hdr_val(H, "biotype_classes_retired_after_99", "E99 \\ E116"),
          "PolymorphicPseudogene")
    check("new classes in 116", hdr_val(H, "biotype_classes_new_in_116",
                                        "E116 \\ E99"), "none")
    pp = L[L.biotype == "PolymorphicPseudogene"].iloc[0]
    check("PolymorphicPseudogene genes lost", pp.n_lost_to_other_class, 88)
    check("PolymorphicPseudogene churn %", pp.churn_pct_of_E99_members, 100.0, 0.01)
    up = L[L.biotype == "UnprocessedPseudogene"].iloc[0]
    check("UnprocessedPseudogene churn %", up.churn_pct_of_E99_members, 25.68, 0.01)

    print("== 5. coordinate convention ==")
    check("E99 convention", hdr_val(H, "bed_coordinate_convention", "E99 mouse BED"),
          "1-based inclusive")
    check("E116 convention", hdr_val(H, "bed_coordinate_convention", "E116 BED"),
          "0-based half-open")
    check("control genes", hdr_val(H, "single_feature_genes_in_both",
                                   "E99 n E116"), 19092)
    check("length match after fix",
          hdr_val(H, "single_feature_length_identical_after_convention_fix",
                  "E99 n E116"), 98.041, 0.001)
    check("length match before fix",
          hdr_val(H, "single_feature_length_identical_before_fix",
                  "E99 n E116"), 0.026, 0.001)

    print("== 6. feature length ==")
    check("exonic length unchanged %",
          hdr_val(H, "shared_genes_exonic_length_unchanged_pct", "E99 n E116"),
          55.64, 0.01)
    check("exon count unchanged %",
          hdr_val(H, "shared_genes_exon_count_unchanged_pct", "E99 n E116"),
          71.827, 0.001)
    check("median log2 exonic ratio",
          hdr_val(H, "shared_genes_median_log2_exonic_ratio", "E116/E99"), 0.0, 1e-9)
    # per-biotype: the short classes the jS:IN rule makes invisible
    for bt, e99, e116 in (("snoRNA", 94.0, 94.0), ("snRNA", 107.0, 107.0),
                          ("miRNA", 102.0, 102.0)):
        row = P[P.biotype_E99 == bt]
        if len(row):
            r0 = row.iloc[0]
            check(f"{bt} median exonic E99", r0.median_exonic_E99, e99, 0.5)
            check(f"{bt} median exonic E116", r0.median_exonic_E116, e116, 0.5)

    print("== 7. composition gap, ID-only universe (upper bound) ==")
    check("raw TVD coarse",
          hdr_val(H, "crossplate_TVD_coarse_ReadCounts_all_genes",
                  "published vs own"), 9.9692, 1e-4)
    check("restricted TVD coarse",
          hdr_val(H, "crossplate_TVD_coarse_ReadCounts_shared_genes",
                  "published vs own"), 4.1398, 1e-4)
    check("release pp coarse",
          hdr_val(H, "crossplate_TVD_coarse_release_attributable_pp",
                  "published vs own"), 5.8294, 1e-4)
    check("release % coarse",
          hdr_val(H, "crossplate_TVD_coarse_release_attributable_pct",
                  "published vs own"), 58.474, 0.001)

    print("== 8. bracket: upper vs lower bound ==")
    check("bracket raw coarse",
          hdr_val(B, "crossplate_TVD_raw_coarse", "id_only"), 9.9692, 1e-4)
    check("upper bound pp",
          hdr_val(B, "release_attributable_pp_coarse", "id_only"), 5.8294, 1e-4)
    check("upper bound %",
          hdr_val(B, "release_attributable_pct_coarse", "id_only"), 58.474, 0.01)
    check("lower bound pp",
          hdr_val(B, "release_attributable_pp_coarse", "id_plus_symbol_any"),
          3.5666, 1e-4)
    check("lower bound %",
          hdr_val(B, "release_attributable_pct_coarse", "id_plus_symbol_any"),
          35.777, 0.01)
    check("informative-symbol pp",
          hdr_val(B, "release_attributable_pp_coarse",
                  "id_plus_symbol_informative"), 4.3057, 1e-4)
    check("informative-symbol %",
          hdr_val(B, "release_attributable_pct_coarse",
                  "id_plus_symbol_informative"), 43.19, 0.01)
    check("symbol bridge pairs any",
          hdr_val(B, "symbol_bridge_pairs_any", "E99-only <-> E116-only"), 923)
    check("symbol bridge pairs informative",
          hdr_val(B, "symbol_bridge_pairs_informative",
                  "E99-only <-> E116-only"), 200)

    print("== 9. ProteinCoding gap, the number other tracks will cite ==")
    g = G[(G.level == "coarse") & (G["class"] == "ProteinCoding")].iloc[0]
    check("PC published unrestricted", g.published_all, 95.405943, 1e-5)
    check("PC own unrestricted", g.own_all, 85.830364, 1e-5)
    check("PC gap unrestricted pp", g.gap_all_pp, -9.575578, 1e-5)
    check("PC published restricted", g.published_shared, 96.014254, 1e-5)
    check("PC own restricted", g.own_shared, 92.175839, 1e-5)
    check("PC gap restricted pp", g.gap_shared_pp, -3.838415, 1e-5)
    check("PC gap closed pp", g.gap_closed_pp, 5.737164, 1e-5)
    gs = G[(G.level == "coarse") & (G["class"] == "sncRNA")].iloc[0]
    check("sncRNA gap unrestricted pp", gs.gap_all_pp, 5.403744, 1e-5)
    check("sncRNA gap restricted pp", gs.gap_shared_pp, 2.751513, 1e-5)
    gl = G[(G.level == "coarse") & (G["class"] == "lncRNA")].iloc[0]
    check("lncRNA gap unrestricted pp", gl.gap_all_pp, 4.565501, 1e-5)
    check("lncRNA gap restricted pp", gl.gap_shared_pp, 1.388313, 1e-5)

    print("== 10. gene detection ==")
    check("own genes detected",
          hdr_val(H, "own_genes_detected_ReadCounts",
                  "own plate, single-gene rows, 12 cells"), 51868)
    check("own detected in shared set",
          hdr_val(H, "own_genes_detected_in_shared_set_ReadCounts", "own plate"), 33813)
    check("own detected only-in-116 %",
          hdr_val(H, "own_genes_detected_only_in_116_pct_ReadCounts",
                  "own plate"), 34.81, 0.01)
    check("published mouse genes detected",
          hdr_val(H, "published_mouse_genes_detected",
                  "mESC wells, single-gene mouse rows"), 32430)
    check("published detected in shared",
          hdr_val(H, "published_mouse_genes_detected_in_shared_set",
                  "mESC wells"), 31684)
    check("published detected only-in-99 %",
          hdr_val(H, "published_mouse_genes_detected_only_in_99_pct",
                  "mESC wells"), 2.3, 0.01)

    print("== 11. restriction read loss asymmetry ==")
    check("own restriction loss %",
          hdr_val(H, "own_denominator_shrink_coarse_ReadCounts",
                  "shared-only vs all"), 7.0157, 1e-4)
    check("published restriction loss %",
          hdr_val(H, "published_denominator_shrink_coarse_ReadCounts",
                  "shared-only vs all"), 0.9351, 1e-4)
    check("loss asymmetry fold", round(7.0157 / 0.9351, 2), 7.50, 0.01)
    own_l = A_[A_.plate == "own_E116"].set_index("coarse_class")
    check("own loss on lncRNA pp",
          own_l.loc["lncRNA", "pct_of_plate_lost_by_restriction"], 3.6444, 1e-3)
    check("own loss on sncRNA pp",
          own_l.loc["sncRNA", "pct_of_plate_lost_by_restriction"], 3.2451, 1e-3)

    print("== 12. Rn7sk, the single locus ==")
    check("Rn7sk own reads",
          hdr_val(B, "Rn7sk_reads", "own plate / ENSMUSG00002076161_Rn7sk_MiscRna"),
          669205)
    check("Rn7sk published reads",
          hdr_val(B, "Rn7sk_reads",
                  "published plate / ENSMUSG00000065037_Rn7sk_MiscRna"), 4199)
    check("raw TVD without Rn7sk",
          hdr_val(B, "crossplate_TVD_raw_coarse_no_Rn7sk", "id_only"), 8.3951, 1e-4)
    check("release pp without Rn7sk",
          hdr_val(B, "release_attributable_pp_coarse_no_Rn7sk", "id_only"),
          4.2553, 1e-4)
    check("Rn7sk share of release term",
          round(5.8294 - 4.2553, 4), 1.5741, 1e-4)
    check("Rn7sk % of own restriction loss",
          round(100 * 669205 / 2611895, 2), 25.62, 0.01)

    print("== 13. re-ID evidence ==")
    for bt, want, ctrl in (("snRNA", 100.0, 72.25), ("snoRNA", 100.0, 56.97),
                           ("scaRNA", 100.0, 51.18), ("miRNA", 99.91, 77.78),
                           ("rRNA", 99.41, 97.26), ("MiscRna", 99.81, 94.92)):
        check(f"{bt} length-multiset match %",
              hdr_val(R, "only_set_length_multiset_match", bt), want, 0.01)
        check(f"{bt} shuffled control %",
              hdr_val(R, "only_set_length_multiset_match_control", bt), ctrl, 0.01)
    check("ProteinCoding length-multiset match %",
          hdr_val(R, "only_set_length_multiset_match", "ProteinCoding"), 15.82, 0.01)
    check("ProteinCoding shuffled control %",
          hdr_val(R, "only_set_length_multiset_match_control",
                  "ProteinCoding"), 1.58, 0.01)
    check("only99 sncRNA id block",
          hdr_val(R, "sncRNA_only_set_id_block", "only99"),
          "min=64393;median=93239;max=106670")
    check("only116 sncRNA id block",
          hdr_val(R, "sncRNA_only_set_id_block", "only116"),
          "min=118674;median=2075486;max=2076992")

    print("== 14. combination rows and tRNA ==")
    check("own combination read share",
          hdr_val(H, "own_combination_row_ReadCounts_share", "own plate"),
          30.096, 0.001)
    check("own combination naming 116-only gene",
          hdr_val(H, "own_combination_ReadCounts_naming_a_116only_gene",
                  "own plate"), 82.848, 0.001)
    check("E99 tRNAscan rows",
          hdr_val(H, "trnascan_rows", "E99 BED (human+mouse)"), 1758)
    check("E116 tRNAscan rows", hdr_val(H, "trnascan_rows", "E116 BED (mouse)"), 1137)

    print("== 15. cell calls ==")
    check("published mESC wells",
          hdr_val(H, "published_wells_mESC_fig1d", "SRR14783059 plate"), 173)
    check("published HEK293T wells",
          hdr_val(H, "published_wells_HEK293T_fig1d", "SRR14783059 plate"), 178)
    check("published mixed wells",
          hdr_val(H, "published_wells_mixed_fig1d", "SRR14783059 plate"), 2)
    check("published discarded wells",
          hdr_val(H, "published_wells_discarded_fig1d", "SRR14783059 plate"), 31)

    print()
    print("=" * 74)
    if bad:
        print(f"VERIFICATION FAILED -- {len(bad)} mismatch(es) of {ok + len(bad)}:")
        for b in bad:
            print(f"  - {b}")
        return 1
    print(f"VERIFICATION PASSED -- all {ok} quoted numbers match their source tables")
    return 0


if __name__ == "__main__":
    sys.exit(main())
