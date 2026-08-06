#!/bin/bash
###############################################################################
# 04_step5_biotype.sh -- what step 5 assigned each read to, single vs multi,
# over all four datasets.
#
# Same tally as demo_analyze/scripts/12_step5_biotype.sh, which is validated on
# three datasets; the changes are that the job list comes from
# scripts/datasets.sh and that FLASH-seq is now in it.
#
# WHAT IT WAS WRITTEN FOR, AND WHAT IT IS FOR NOW
# -----------------------------------------------
# On three datasets it answered "why does the own library multimap 2.4x the
# published plate at matched read length" -- not residual rDNA (0.91% of
# multimapper reads even touch an rRNA gene) but small RNA: snRNA 15x and miRNA
# 2.6x the plate's rate after standardising to its read-length distribution.
# See demo_analyze/README.md, "What the multimappers are (step 5)".
#
# FLASH-seq is the fourth dataset and the first NON-VASA protocol in it, so the
# question it answers here is different: is the small-RNA excess a VASA
# property or a property of this library. FLASH-seq reaches STAR through
# 00_fs_ribo.sh, so its rRNA is removed by the same script and the same
# reference -- without that stage the comparison would be measuring the pipeline.
#
# NOTHING IS RE-RUN. Read-only tally of
#   *_Aligned.out.singlemappers_genes.bed.gz          (STAR NH:i:1)
#   *_Aligned.out.nsorted.multimappers_genes.bed.gz   (STAR NH:i:>=2)
#
# WHAT ONE ROW OF THOSE BEDS IS
# -----------------------------
# One read x one alignment locus x one overlapping annotation feature, so a
# read appears many times: 1.4 rows/read in the singlemapper BED, 3.7 in the
# multimapper one (cell 002). Every count here is therefore taken PER READ,
# grouping on col 4 (the read name). Both files are already grouped by read --
# the multimapper BED by construction (`sort -k4`, hence `nsorted`), the
# singlemapper BED because it is position-sorted and a unique read's rows all
# share one position. Verified on cell 002: 0 non-contiguous read names in
# either file. The tally hard-errors if that ever stops holding.
#
# THE NAME FIELD, AND ITS ONE EXCEPTION
# -------------------------------------
# col 6 is GENEID_SYMBOL_BIOTYPE_LABEL (gtf2bed_vasa.py rule 5; no symbol
# contains '_', so the 4-way split is safe) -- e.g.
#   ENSMUSG00000051951_Xkr4_ProteinCoding_intron
# tRNA rows are the exception: they carry the GtRNAdb locus name alone and no
# underscore at all -- e.g. `1.tRNA1006-GluCTC`. They are counted as biotype
# and label "tRNA", which is what they are. 4,525 of 2,229,555 rows on cell 002.
#
# THREE COUNTS PER READ, AND WHY IT IS NOT ONE
# --------------------------------------------
# A multimapper's loci usually do NOT agree on a biotype -- 84.8% of cell 002's
# multimapper reads touch more than one. So a single "the biotype of this read"
# column would put most of the library in a MIXED bucket and answer nothing.
# The three metrics split that apart:
#
#   any     each read counted ONCE PER DISTINCT BIOTYPE it touches. Columns sum
#           to more than the read count; this is the one that answers "do the
#           multimappers touch rRNA at all".
#   ngenes  distinct gene ids the read was assigned to (capped at 20, which is
#           also STAR's --outFilterMultimapNmax) -- how ambiguous, not what.
#   tot     per-read bookkeeping: reads in the BED, rows, whether the read's
#           biotypes agreed, and its exon/intron call.
#
# DENOMINATORS
# ------------
# Not every mapped read reaches a gene. `bed_reads / star_reads` is the share
# that did, against that cell's own Log.final.txt (uniquely-mapped for single,
# multiple-loci for multi). Reads are lost between STAR and the BED for reasons
# that live in deal_with_*mappers.sh: no overlapping feature, antisense to it
# (`stranded=y` keeps readstrand==refstrand only), or spanning past both ends
# of the feature (`jS:OUT`, dropped). `too many loci` reads are absent from the
# BAM as alignments, so they are in NEITHER file and in no denominator here.
#
# THE PLATE IS SCORED ON MOUSE ANNOTATION ONLY -- see MOUSECONTIG below.
#
# WHAT THIS CANNOT CONTROL FOR, AND IT MATTERS
# --------------------------------------------
# Three of the four are annotated against GRCm39 Ensembl 116 (both VASA own runs
# and FLASH-seq -- pipeline_fs.sh points at the same IntronExonTrna.v2.bed). The
# published plate is the exception: the mixed GRCh38+GRCm38 E99 BED.
# Counted per distinct gene id on the mouse side of each, the releases agree on
# almost every biotype -- miRNA 2206 vs 2207, snoRNA 1507 vs 1507, snRNA 1381
# vs 1385, MiscRna 562 vs 562, rRNA 354 vs 354, ProteinCoding 21818 vs 21933,
# ProcessedPseudogene 9312 vs 10003 -- with ONE large exception: lncRNA, 32889
# in E116 against 9959 in E99, a 3.3x difference. So any lncRNA gap involving THE
# PLATE is confounded by annotation content and cannot be read as biology; the
# small-RNA and pseudogene classes can, and the own-vs-FLASH-seq comparison is
# free of the release effect entirely because both sit on E116.
#
# Output: tables/cross/step5_biotype_per_cell.tsv   dataset cell class biotype reads
#         tables/cross/step5_genes_per_read.tsv     dataset cell class ngenes  reads
#         tables/cross/step5_assign_totals.tsv      one row per dataset/unit/class
#         tables/cross/step5_biotype_by_length.tsv  biotype x read length, pooled
#
# COST
#   On the three VASA datasets alone (410 BEDs / 7.5 GB) this took 6 min at -P 8
#   and peaked at 5.8 GB. FLASH-seq adds 20 BEDs but they are far bigger -- one
#   library is ~24M reads against a VASA cell's ~1M -- so the peak is set by the
#   per-read-name hash of ONE FLASH-seq BED, not by the count of them.
#   MEASURED HERE (job 51370952): 852 BEDs / 17.12 GB -- own130 2.61, own75 3.01,
#   plate 4.42, FLASH-seq 7.08 -- in 18m16s, peak 19.7 GB against the 64G request.
#   The guess above was right about the mechanism: FLASH-seq is 20 of the 852
#   files and 41% of the bytes, and the peak tripled (5.8 -> 19.7 GB) while the
#   file count doubled. 64G stays; 32G would be uncomfortably close.
#
#   sbatch -c 16 --mem=64G -t 240 --wrap="scripts/04_step5_biotype.sh"
###############################################################################
set -euo pipefail
ROOT=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
cd "$ROOT"

source scripts/datasets.sh          # the four datasets and their per-unit stems

# The published plate is mapped and annotated against the MIXED human+mouse
# reference, so a mouse well's reads can pick up human annotation -- measured on
# three wells, 12.6-15.8% of multimapper BED rows and ~20% of multimapper reads
# touch a GRCh38_ contig. Comparing that against three mouse-only libraries would
# make the annotation, not the biology, part of the answer. Plate rows are
# therefore restricted to mouse contigs, and what that discards is written out as
# rows_offspecies / reads_offspecies rather than hidden. Set MOUSECONTIG='' to
# score the plate against both species. The other three datasets are mouse-only
# references and are never filtered.
MOUSECONTIG=${MOUSECONTIG-^GRCm38_}

NPAR=${NPAR:-8}          # BEDs decompressed at once
NGCAP=${NGCAP:-20}       # ngenes is capped here; STAR's own limit is 20 loci

BIO=$OUTROOT/tables/cross/step5_biotype_per_cell.tsv
NGN=$OUTROOT/tables/cross/step5_genes_per_read.tsv
TOT=$OUTROOT/tables/cross/step5_assign_totals.tsv
BLN=$OUTROOT/tables/cross/step5_biotype_by_length.tsv
MLD=$OUTROOT/tables/cross/mapped_length_dist.tsv   # scripts/03, only for the cross-check

TMP=$(mktemp -d "${TMPDIR:-/tmp}/s5bt.XXXXXX")
trap 'rm -rf "$TMP"' EXIT


# --- one BED ---------------------------------------------------------------
# Emits long form: dataset <TAB> cell <TAB> class <TAB> metric <TAB> key <TAB> value

tally() {
    key=$1; cell=$2; cls=$3; bed=$4; out=$5
    keep=''
    if [ "$key" = plate ]; then keep=$MOUSECONTIG; fi   # not `&&`: the child runs under set -e
    zcat "$bed" | awk -v ds="$(ds_label "$key")" -v cell="$cell" -v cls="$cls" \
                      -v cap="$NGCAP" -v src="$bed" -v keep="$keep" '
      # Read length out of the CIGAR that deal_with_*mappers.sh parked in col 7
      # as `CG:109M21S;nM:0;jS:IN`. M/I/S/=/X consume the read, D/N/H/P do not,
      # and STAR never hard-clips -- so this is the full length entering STAR,
      # not the aligned span. Every locus of one read agrees on it (checked on
      # cell 002: 0 of 602,844 multimapper reads disagreed), which is why the
      # per-read value can be taken off any row; max() is belt and braces.
      function clen(c,   i, ch, num, s) {
          sub(/;.*/, "", c); sub(/^CG:/, "", c); s = 0; num = ""
          for (i = 1; i <= length(c); i++) { ch = substr(c, i, 1)
              if (ch ~ /[0-9]/) num = num ch
              else { if (ch ~ /[MIS=X]/) s += num + 0; num = "" } }
          return s
      }
      function flush(   b, g, nb, ng, nl, only) {
          if (!open) return
          open = 0
          if (!kept) { t["reads_offspecies"]++; return }   # every row filtered out
          reads++
          nb = 0; for (b in bt) { nb++; any[b]++; bl[b SUBSEP rlen]++ }
          ng = 0; for (g in gid) ng++
          nl = 0; for (b in lab) { nl++; only = b }
          t[nb == 1 ? "reads_1biotype" : "reads_mixed"]++
          t[nl == 1 ? only : "mixed"]++
          ngn[ng > cap ? cap : ng]++
          len[rlen]++
          delete bt; delete gid; delete lab
      }
      { if ($4 != prev) {
            flush()
            # Grouping is by adjacency, so a read name coming back after an
            # intervening one would be silently counted twice. Refuse instead.
            if ($4 in seen) { printf "read %s not contiguous in %s\n", $4, src > "/dev/stderr"; exit 3 }
            seen[$4] = 1; prev = $4; open = 1; kept = 0; rlen = 0 }
        if (keep != "" && $1 !~ keep) { skipped++; next }
        rows++; kept++
        L = clen($7); if (L > rlen) rlen = L
        n = split($6, f, "_")
        if (n == 1) { bt["tRNA"] = 1; gid[$6]  = 1; lab["tRNA"] = 1 }   # GtRNAdb locus name
        else        { bt[f[3]]   = 1; gid[f[1]] = 1; lab[f[4]]  = 1 } }
      END { flush()
            for (b in any) printf "%s\t%s\t%s\tany\t%s\t%d\n",    ds, cell, cls, b, any[b]
            for (g in ngn) printf "%s\t%s\t%s\tngenes\t%s\t%d\n", ds, cell, cls, g, ngn[g]
            for (k in len) printf "%s\t%s\t%s\tlen\t%s\t%d\n",    ds, cell, cls, k, len[k]
            for (k in bl)  { split(k, a, SUBSEP)
                             printf "%s\t%s\t%s\tbtlen\t%s|%s\t%d\n", ds, cell, cls, a[1], a[2], bl[k] }
            for (k in t)   printf "%s\t%s\t%s\ttot\t%s\t%d\n",    ds, cell, cls, k, t[k]
            printf "%s\t%s\t%s\ttot\tbed_reads\t%d\n",    ds, cell, cls, reads
            printf "%s\t%s\t%s\ttot\tbed_rows\t%d\n",     ds, cell, cls, rows
            printf "%s\t%s\t%s\ttot\trows_offspecies\t%d\n", ds, cell, cls, skipped }' > "$out"
}
export -f tally ds_label
export NGCAP MOUSECONTIG

# --- job list: key unit class bed outfile ----------------------------------
# Space-separated, and that is safe here because every field is a KEY, a unit id
# or a path -- none contains whitespace. The dataset LABEL (which does) is
# resolved inside the worker via ds_label, never carried through xargs.

: > "$TMP/jobs"
i=0
add() {   # add <key> <unit> <class> <bed>
    [ -s "$4" ] || { echo "missing $4" >&2; exit 1; }
    i=$((i+1))
    printf '%s %s %s %s %s\n' "$1" "$2" "$3" "$4" "$TMP/t.$i" >> "$TMP/jobs"
}
for key in $DS_KEYS; do
    while IFS=$'\t' read -r unit map_stem ribo_stem; do
        add "$key" "$unit" single "${map_stem}_E99_Aligned.out.singlemappers_genes.bed.gz"
        add "$key" "$unit" multi  "${map_stem}_E99_Aligned.out.nsorted.multimappers_genes.bed.gz"
    done < <(ds_units "$key")
done

echo "scanning $(wc -l < "$TMP/jobs") BEDs, $NPAR at a time (plate source=$PSOURCE)..."
# `set -eo pipefail` in the child so a zcat that dies takes the shard with it
# rather than leaving awk to exit 0 over an empty stream.
xargs -P "$NPAR" -n 5 bash -c 'set -eo pipefail; tally "$0" "$1" "$2" "$3" "$4"' < "$TMP/jobs"

cat "$TMP"/t.* > "$TMP/all"

# --- split into the three tables -------------------------------------------

{ printf 'dataset\tcell\tclass\tbiotype\treads\n'
  awk -F'\t' '$4=="any" {print $1"\t"$2"\t"$3"\t"$5"\t"$6}' "$TMP/all" \
    | sort -t$'\t' -k1,1 -k2,2 -k3,3 -k5,5nr
} > "$BIO"
echo "wrote $BIO ($(( $(wc -l < "$BIO") - 1 )) rows)"

{ printf 'dataset\tcell\tclass\tngenes\treads\n'
  awk -F'\t' '$4=="ngenes" {print $1"\t"$2"\t"$3"\t"$5"\t"$6}' "$TMP/all" \
    | sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4n
} > "$NGN"
echo "wrote $NGN ($(( $(wc -l < "$NGN") - 1 )) rows)"

# Biotype x read length. POOLED OVER CELLS, unlike the three tables above: kept
# per cell it is 1.2M rows, and no question this table answers uses the cell
# dimension. `reads_at_length` is that class's own read count at that length, so
# a row is self-contained -- pct = reads / reads_at_length.
{ printf 'dataset\tclass\tlength\tbiotype\treads\treads_at_length\n'
  awk -F'\t' '$4=="len"   { n[$1 SUBSEP $3 SUBSEP $5] += $6 }
              $4=="btlen" { split($5, a, "|"); b[$1 SUBSEP $3 SUBSEP a[2] SUBSEP a[1]] += $6 }
    END { for (i in b) { split(i, a, SUBSEP)
            printf "%s\t%s\t%s\t%s\t%d\t%d\n", a[1], a[2], a[3], a[4], b[i],
                   n[a[1] SUBSEP a[2] SUBSEP a[3]] } }' "$TMP/all" \
    | sort -t$'\t' -k1,1 -k2,2 -k3,3n -k4,4
} > "$BLN"
echo "wrote $BLN ($(( $(wc -l < "$BLN") - 1 )) rows)"

# star_reads: the class's own denominator out of that cell's Log.final.txt.
# The logs are streamed one at a time -- 384-cell runs are what made ARG_MAX a
# real bug in this folder once.
awk '{print $1"\t"$2"\t"$3"\t"$4}' "$TMP/jobs" | sort -u | while IFS=$'\t' read -r key cell cls bed; do
    log=${bed%%_Aligned.out.*}_Log.final.txt
    [ -s "$log" ] || { echo "missing $log" >&2; exit 1; }
    pat='Uniquely mapped reads number'
    [ "$cls" = single ] || pat='Number of reads mapped to multiple loci'
    printf '%s\t%s\t%s\t%s\n' "$(ds_label "$key")" "$cell" "$cls" \
        "$(awk -F'\t' -v p="$pat" 'index($1,p){gsub(/[^0-9]/,"",$2); print $2; exit}' "$log")"
done > "$TMP/star"

{ printf 'dataset\tcell\tclass\tstar_reads\tbed_reads\tbed_rows\treads_1biotype\treads_mixed\texon\tintron\ttRNA\tmixed_label\treads_offspecies\trows_offspecies\n'
  awk -F'\t' 'NR==FNR { s[$1 SUBSEP $2 SUBSEP $3] = $4; next }
              $4=="tot" { v[$1 SUBSEP $2 SUBSEP $3, $5] = $6; k[$1 SUBSEP $2 SUBSEP $3] = 1 }
    END { for (i in k) { split(i, a, SUBSEP)
            printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", a[1], a[2], a[3],
              s[i], v[i,"bed_reads"], v[i,"bed_rows"], v[i,"reads_1biotype"], v[i,"reads_mixed"],
              v[i,"exon"], v[i,"intron"], v[i,"tRNA"], v[i,"mixed"],
              v[i,"reads_offspecies"], v[i,"rows_offspecies"] } }' "$TMP/star" "$TMP/all" \
    | sort -t$'\t' -k1,1 -k2,2 -k3,3
} > "$TOT"
echo "wrote $TOT ($(( $(wc -l < "$TOT") - 1 )) rows)"

# --- sanity: a BED cannot hold reads STAR never placed in that class --------

awk -F'\t' 'NR>1 && $5 + $13 > $4 { printf "  %s cell %s %s: bed_reads %d + offspecies %d > star_reads %d\n", $1,$2,$3,$5,$13,$4; bad++ }
            END { exit (bad>0) }' "$TOT" \
  || { echo "step-5 BED has more reads than STAR reported -- not trusting $TOT" >&2; exit 1; }
echo "every cell's BED read count is within its STAR class total."

# --- sanity: MOUSECONTIG must not have discarded a whole dataset ------------
# The plate filter is a string match on the contig name, so it silently matches
# NOTHING the moment the plate is re-mapped to a reference whose contigs are not
# prefixed. That is not hypothetical: the E116 remap (2026-08-06) put the plate
# on mouse-only GRCm39, whose contigs are plain `1`, `2`, ..., and the default
# `^GRCm38_` threw away all 173 wells. Every table was still written, with zeros,
# and the failure only surfaced 17 minutes later as `division by zero` inside an
# awk summary -- which names neither the dataset nor the cause.
awk -F'\t' 'NR>1 { t[$1 SUBSEP $3] += $4; s[$1 SUBSEP $3] += $5 }
  END { for (i in t) if (t[i] > 0 && s[i] == 0) { split(i,a,SUBSEP)
          printf "  %s / %s: STAR placed %d reads, none reached a gene\n", a[1], a[2], t[i]; bad++ }
        exit (bad>0) }' "$TOT" \
  || { echo "a whole dataset x class has zero BED reads. If that is the plate," >&2
       echo "MOUSECONTIG=${MOUSECONTIG:-<unset>} matched no contig -- re-run with" >&2
       echo "MOUSECONTIG='' when the plate is on a mouse-only reference." >&2; exit 1; }

# --- sanity: the CIGAR-derived lengths against step 10's independent tally ---
# scripts/10 reads length(SEQ) straight out of the BAM; this script rebuilds it
# from the CIGAR two stages downstream. They must agree in the only way they
# can: at every length, the BED (which loses reads that reached no gene, and on
# the plate reads with only human annotation) is a SUBSET of the BAM. A
# CIGAR-parsing slip -- counting D/N, or dropping S -- shifts the whole
# distribution and breaks this immediately.

if [ -s "$MLD" ]; then
    awk -F'\t' -v s5="$BLN" '
        NR>1 { cat = ($3=="unique") ? "single" : (($3=="multi") ? "multi" : "")
               if (cat != "") bam[$1 SUBSEP cat SUBSEP $4] += $5 }
        END { while ((getline l < s5) > 0) { split(l, f, "\t")
                 if (f[1] == "dataset") continue
                 k = f[1] SUBSEP f[2] SUBSEP f[3]
                 if (k in done) continue; done[k] = 1
                 if (f[6] > bam[k]) { printf "  %s %s %s nt: bed %d > bam %d\n", f[1], f[2], f[3], f[6], bam[k]; bad++ } }
              exit (bad > 0) }' "$MLD" \
      || { echo "step-5 read lengths do not fit inside step-4 tally -- not trusting $BLN" >&2; exit 1; }
    echo "read lengths fit inside scripts/10's BAM tally at every length."
else
    echo "note: $MLD missing -- length cross-check skipped (run scripts/10 first)"
fi

# --- what it says ----------------------------------------------------------

echo
echo "reads reaching a gene, % of that STAR class:"
awk -F'\t' 'NR>1 { s[$1 SUBSEP $3] += $5; t[$1 SUBSEP $3] += $4; m[$1 SUBSEP $3] += $8
                   r[$1 SUBSEP $3] += $6; e[$1 SUBSEP $3] += $9; o[$1 SUBSEP $3] += $13 }
  END { printf "%-24s %-7s %12s %12s %7s %8s %9s %7s %7s\n", "dataset","class","star_reads","bed_reads","pct","rows/read","%mixedBT","%exon","%offsp"
        for (i in s) { split(i,a,SUBSEP)
          # Guarded: the check above already refuses to get here with s[i]==0,
          # but a summary print must not be the thing that decides a run failed.
          printf "%-24s %-7s %12d %12d %6.1f%% %8.2f %8.1f%% %6.1f%% %6.1f%%\n", a[1], a[2], t[i], s[i],
                 t[i] ? 100*s[i]/t[i] : 0, s[i] ? r[i]/s[i] : 0, s[i] ? 100*m[i]/s[i] : 0,
                 s[i] ? 100*e[i]/s[i] : 0, t[i] ? 100*o[i]/t[i] : 0 } }' "$TOT" | (read -r h; echo "$h"; sort)

echo
echo "biotypes touched, % of that class's BED reads (top 12):"
awk -F'\t' 'NR==FNR && FNR>1 { n[$1 SUBSEP $3] += $5; next }
            FNR>1 { b[$1 SUBSEP $3 SUBSEP $4] += $5 }
  END { for (i in b) { split(i,a,SUBSEP)
          printf "%s\t%s\t%s\t%.2f\n", a[1], a[2], a[3], 100*b[i]/n[a[1] SUBSEP a[2]] } }' "$TOT" "$BIO" \
  | sort -t$'\t' -k1,1 -k2,2 -k4,4nr \
  | awk -F'\t' '{ if ($1"\t"$2 != last) { c=0; last=$1"\t"$2; printf "\n%s / %s\n", $1, $2 }
                  if (++c <= 12) printf "  %-32s %6.2f%%\n", $3, $4 }'

echo
echo "the same, AT MATCHED READ LENGTH -- % of multimapper reads of that length:"
# The dataset labels are passed in from ds_label rather than written out again:
# a second copy of them is exactly what scripts/datasets.sh exists to prevent,
# and a typo in one would print a row of "-" that reads as missing data.
DSLABELS=$(for k in $DS_KEYS; do printf '%s|' "$(ds_label "$k")"; done)
awk -F'\t' -v want="snRNA,miRNA,snoRNA,MiscRna,ProcessedPseudogene" -v dsl="${DSLABELS%|}" '
  BEGIN { split(want, w, ","); for (i in w) keep[w[i]] = 1
          split("25,40,55,70,75", L, ",")
          nd = split(dsl, D, "|") }
  NR>1 && $2=="multi" && ($4 in keep) && $6 >= 1000 { p[$4 SUBSEP $1 SUBSEP $3] = 100*$5/$6 }
  END { for (b in keep) { printf "\n  %s\n", b
          printf "  %-26s", "% at length"; for (i=1; i<=5; i++) printf "%8s", L[i] " nt"; printf "\n"
          for (j=1; j<=nd; j++) { printf "  %-26s", D[j]
            for (i=1; i<=5; i++) { k = b SUBSEP D[j] SUBSEP L[i]
              if (k in p) printf "%7.1f%%", p[k]; else printf "%8s", "-" }
            printf "\n" } } }' "$BLN"
