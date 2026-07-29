#!/bin/bash
###############################################################################
# shortspecies_from_bed.sh -- measure the short-biotype recovery DIRECTLY from
# the step-5 BEDs, without waiting for steps 6-7.
#
# WHY THIS IS A REAL MEASUREMENT AND NOT A SHORTCUT
# -------------------------------------------------
# Step 6 keeps a read on a non-spliceable biotype ONLY when its jS tag is IN
# (countTables_2pickle_cellsSpliced.py: it drops any row whose Biotype is not in
# biotypeWsplicing unless jSs == 'IN'). tRNA, miRNA, snoRNA, snRNA, MiscRna and
# rRNA are all outside biotypeWsplicing. So for those classes the question
# "does this arm detect them at all" is decided entirely at the BED, by whether
# any row for that biotype carries jS:IN -- which is what this script counts.
#
# That makes this an independent check on the step-7 tables rather than a
# substitute for them: step 7's numbers additionally pass through gene assignment
# (best nM, jS:IN priority, non-spliceable-biotype priority, tie-breaking) and
# reduceGeneName, so its per-gene totals will be SMALLER than the row counts
# here. What must agree is the qualitative answer -- zero jS:IN rows for a
# biotype cannot become a non-empty table for it, and a large jS:IN count cannot
# become an empty one.
#
# The BED columns, from deal_with_singlemappers.sh's final print:
#   1 chr  2 readstart  3 readend  4 readname  5 readstrand  6 refname
#   7 info(;CG:..;nM:..;jS:..)  8 refend-refstart  9 cov
# The biotype is the second-to-last '_'-separated field of refname (the last is
# the exon/intron label), except for tRNA rows, whose refname is rewritten by
# step 6; here the raw form is matched on the substring 'tRNA', exactly as step 6
# does before it rewrites.
#
# Usage: shortspecies_from_bed.sh <scratch_root> <out.tsv> <libs...>
###############################################################################
set -eo pipefail

ROOT=${1:?scratch root}
OUT=${2:?output tsv}
shift 2
LIBS="$*"
[ -n "$LIBS" ] || { echo "give at least one library" >&2; exit 1; }

printf 'arm\tlibrary\tbiotype\trows_total\trows_jsIN\tdistinct_features_jsIN\n' > "$OUT"

for arm in native vasalen; do
    S="$ROOT/$arm/cells"
    [ -d "$S" ] || { echo "no $S -- skipping $arm"; continue; }
    for lib in $LIBS; do
        bed="$S/${lib}_cbc_noumi_E99_Aligned.out.singlemappers_genes.bed.gz"
        [ -s "$bed" ] || { echo "  no BED for $arm/$lib" >&2; continue; }
        # One pass per BED, all biotypes at once. Two counters per class: every
        # row that overlaps such a feature, and the subset with jS:IN -- the
        # difference is precisely what the containment rule discards.
        zcat "$bed" | awk -v arm="$arm" -v lib="$lib" -v OFS='\t' '
        {
            ref=$6; info=$7
            # biotype = second-to-last _-field of refname; tRNA is matched as a
            # substring because its refname carries the isotype form.
            n=split(ref, p, "_")
            bt=(n>=2 ? p[n-1] : "NA")
            if (ref ~ /tRNA/) bt="tRNA"
            if (bt=="tRNA" || bt=="miRNA" || bt=="snoRNA" || bt=="snRNA" ||
                bt=="MiscRna" || bt=="scaRNA" || bt=="rRNA" || bt=="ribozyme" ||
                bt=="MtTrna" || bt=="MtRrna") {
                tot[bt]++
                js=info; sub(/.*;jS:/,"",js)
                if (js=="IN") { inn[bt]++; feat[bt";"ref]=1 }
            }
        }
        END {
            for (b in tot) {
                nf=0
                for (k in feat) { split(k,q,";"); if (q[1]==b) nf++ }
                print arm, lib, b, tot[b], inn[b]+0, nf
            }
        }' >> "$OUT"
        echo "  done $arm/$lib" >&2
    done
done

echo "=== pooled per arm x biotype ==="
awk -F'\t' 'NR>1{t[$1"\t"$3]+=$4; i[$1"\t"$3]+=$5}
    END{printf "%-9s %-10s %14s %14s %9s\n","arm","biotype","rows_total","rows_jsIN","%IN";
        for(k in t){split(k,p,"\t"); printf "%-9s %-10s %14d %14d %8.2f%%\n",p[1],p[2],t[k],i[k],100*i[k]/t[k]}}' "$OUT" | sort
echo "written: $OUT"
