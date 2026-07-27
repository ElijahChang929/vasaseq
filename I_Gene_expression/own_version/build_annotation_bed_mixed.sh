#!/bin/bash
###############################################################################
# build_annotation_bed_mixed.sh -- the human+mouse ("mixed") annotation BED that
# step 5 intersects against, rebuilt from primary sources WITH cytoplasmic tRNA.
#
# THIS SCRIPT LIVES IN THE GIT REPO ON PURPOSE. The file it replaces was built
# by ${VASA_REFS}/mixed/build/gtf_to_homemade_bed.py, which sits untracked next
# to its own outputs -- the exact arrangement that already lost
# build_mouse_reference.sh. Keep reference builders in git.
#
# Sibling of build_annotation_bed.sh (mouse GRCm39 / Ensembl 116), deliberately
# NOT merged into it, for the same reasons build_rrna_reference_mixed.sh is a
# separate script: two species, a different Ensembl release, species-prefixed
# contigs, and two GtRNAdb sources instead of one.
#
# ===========================================================================
# WHAT WAS WRONG WITH THE FILE THIS REPLACES
# ===========================================================================
# Human_Mouse_ensembl99.homemade_IntronExonTrna.bed is 1,165,949 rows and
# contains ZERO tRNA rows despite the name. Consequence, measured on the
# authors' own species-mixing control (vasaplate_out_rrnav2, 384 cells):
#
#     vasaplate_out_rrnav2_mapStats.log :  Total reads assigned to tRNA: 0
#     vasaplate_out_rrnav2_tRNA.ReadCounts.tsv : header only, zero rows
#
# while the published table for the same library (GEO GSM5369495) carries
# 1,804 tRNA rows. Of its 77,207 simple (non-combination) gene rows our re-run
# already reproduces 72,613; of the 4,594 we miss, 1,130 are tRNA. tRNA is the
# single largest structural gap and the only one with a systematic cause.
#
# ===========================================================================
# TWO DECISIONS, BOTH SETTLED BY MEASUREMENT AGAINST THE PUBLISHED TABLE
# ===========================================================================
# The brief was to match what the original authors did rather than to improve
# on it. Both open questions turned out to be empirically decidable.
#
# (1) DO tRNA NAMES CARRY A SPECIES TAG?  No -- and they must not.
#
#     Human and mouse both have chr1..chr19, so a bare tRNAscan id such as
#     `13.tRNA1495-ValAAC` is in principle ambiguous in a two-species
#     reference. Measured, converting both GtRNAdb sets through the name_map:
#
#         hg38 loci                619
#         mm10 loci               1139
#         names colliding            1   <- one, out of 1757
#
#     and the published table's tRNA rows carry no species marker of any kind
#     (`6.tRNA60.IleAAT`, `14.tRNA3.ProTGG`). So the authors used one flat
#     namespace and accepted the collision. We reproduce that. The BED's
#     CONTIG column is still species-prefixed -- it has to be, to match the
#     STAR index -- but the NAME field is bare, exactly as published.
#
# (2) WHICH GtRNAdb RELEASE?  Not one that still exists.
#
#     Of the 1,130 simple tRNA names in the published table, only 503 appear
#     in the current hg38+mm10 union. The tell is the isotype label: 110
#     published rows say `Undet`, the retired tRNAscan-SE spelling, and the
#     current GtRNAdb writes `Und` -- zero published rows use `Und`, zero
#     current rows use `Undet`. Normalising that raises the match only to 511,
#     because GtRNAdb has also re-run tRNAscan-SE and RENUMBERED the loci.
#
#     LOCUS-LEVEL IDS THEREFORE CANNOT BE REPRODUCED, and chasing them is
#     wasted effort. What survives is the biology: 59 of the 62 published
#     isotype+anticodon classes (95%) are present in the current sets -- and
#     isotype is exactly what countTables_fromPickle.py collapses tRNA to
#     (`t.rsplit('.')[-1]` -> `ValAAC`). Compare at that level, not by id.
#
# ===========================================================================
# COORDINATES: asis, BECAUSE THAT IS WHAT THE AUTHORS USED
# ===========================================================================
# The mouse builder defaults to BED_COORD=fix because GTF is 1-based inclusive
# and BED is 0-based half-open, so copying GTF numbers verbatim puts every
# start 1 bp too high and costs reads on short non-splicing features -- +8.3%
# on real mouse cell 002 when corrected.
#
# Here the authors' convention is knowable, so it wins. Test: if they had used
# true 0-based while our current run uses asis, short non-splicing biotypes
# would sit ~8% BELOW protein-coding in a per-gene ours/published ratio.
# Measured on 72,613 shared genes of vasaplate_out_rrnav2_total.UFICounts.tsv
# against GSM5369495:
#
#     ProteinCoding, pub>=100    n=22523   median log2(ours/pub) = -0.0248
#     non-splicing,  pub>=100    n=  213   median log2(ours/pub) = -0.0087
#     offset (non-splicing - protein-coding)      = +0.0161 log2 = +1.12%
#
# +1.12%, not -8%, and in the wrong direction for a coordinate mismatch. The
# authors used asis. BED_COORD=asis is therefore the DEFAULT here, the inverse
# of the mouse script's default, and the difference is deliberate.
#
# BED_COORD=fix is still available and is the better file for one's own
# science -- it just is not the file that reproduces this paper.
#
# ===========================================================================
# A QUIRK THAT IS REPRODUCED ON PURPOSE, NOT FIXED
# ===========================================================================
# countTables_2pickle_cellsSpliced.py:93 diverts any gene whose name contains
# the literal substring `tRNA` into the tRNA table. Exactly one non-tRNA row in
# this BED matches, via the `l`+`tRNA` in "vaultRNA":
#
#     ENSG00000270123_VTRNA2.1_vaultRNA
#
# It is absent from the published table AND from our own rrnav2 table, so the
# authors' pipeline did the same thing ours does. Left alone. The human
# VTRNA1.x genes do not match -- their biotype is MiscRna and `VTRNA` is
# uppercase, and the test is case-sensitive.
#
# ===========================================================================
# OUTPUTS (in $OUTDIR)
# ===========================================================================
#   Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed
#   Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.provenance.txt
#
# The file this replaces is left untouched so the two can be compared.
# Idempotent: re-running overwrites its own outputs and reuses any download.
#
# Usage:  ./build_annotation_bed_mixed.sh              # asis + tRNA + validate
#         BED_COORD=fix ./build_annotation_bed_mixed.sh
#         WITH_TRNA=no  ./build_annotation_bed_mixed.sh
#         VALIDATE=no   ./build_annotation_bed_mixed.sh
###############################################################################
set -euo pipefail

# ---------------------------------------------------------------------------
# EDIT ME
# ---------------------------------------------------------------------------
VASA_REFS="${VASA_REFS:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq}"
OUTDIR="${OUTDIR:-${VASA_REFS}/mixed}"
BUILD="${BUILD:-${OUTDIR}/build/annotation}"

# The species-prefixed combined GTF the STAR index was built from. Using it
# rather than the two per-species GTFs is what makes the contig prefix
# automatic for gene rows -- and guarantees the BED agrees with the index,
# since it is literally the same file.
GTF="${GTF:-${OUTDIR}/build/combined.gtf}"

# Authority for "does this contig exist". The mouse builder checks tRNA contigs
# against the gene rows; that is too strict here, because mm10 puts 3 tRNA loci
# on GL456213.1 and GL456367.1, which are in the genome and in the STAR index
# but carry no Ensembl gene. Reads can map there, so the tRNAs are real.
CHRNAMES="${CHRNAMES:-${OUTDIR}/star_index_74/chrName.txt}"

OUT_BED="${OUT_BED:-${OUTDIR}/Human_Mouse_ensembl99.homemade_IntronExonTrna.v2.bed}"

# The file being replaced. Used ONLY by the validation step as the
# reproduction target; the build does not depend on it existing.
V1_BED="${V1_BED:-${OUTDIR}/Human_Mouse_ensembl99.homemade_IntronExonTrna.bed}"

# GtRNAdb, per species: <assembly-dir> <file-prefix> <contig-prefix>
# Hsapi38 = hg38 = GRCh38 and Mmusc10 = mm10 = GRCm38, both matching Ensembl 99.
# Mmusc39/mm39 is GRCm39 and is the WRONG assembly here -- coordinates would
# silently disagree with the BAMs.
GTRNADB_SPECS="${GTRNADB_SPECS:-Hsapi38:hg38:GRCh38_ Mmusc10:mm10:GRCm38_}"
GTRNADB_BASE="${GTRNADB_BASE:-https://gtrnadb.ucsc.edu/genomes/eukaryota}"

WITH_TRNA="${WITH_TRNA:-yes}"   # yes | no
TRNA_SET="${TRNA_SET:-all}"     # all (incl. low-confidence tRX-*) | high
BED_COORD="${BED_COORD:-asis}"  # asis (the authors' convention) | fix
VALIDATE="${VALIDATE:-yes}"     # yes | no
# ---------------------------------------------------------------------------

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GTF2BED="${HERE}/gtf2bed_vasa.py"
PY="${PY:-/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/python}"

say() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

case "$BED_COORD" in fix|asis) ;; *) die "BED_COORD must be fix or asis" ;; esac
case "$TRNA_SET"  in all|high) ;; *) die "TRNA_SET must be all or high"  ;; esac
case "$WITH_TRNA" in yes|no)   ;; *) die "WITH_TRNA must be yes or no"   ;; esac
[ -x "$PY" ]      || die "python not found: $PY"
[ -s "$GTF2BED" ] || die "gtf2bed_vasa.py not found next to this script"
[ -s "$GTF" ]     || die "combined GTF not found: $GTF"
[ -s "$CHRNAMES" ] || die "STAR chrName.txt not found: $CHRNAMES"
mkdir -p "$BUILD"

###############################################################################
say "=== STEP 1: combined GTF -> exon/intron rows (coord=$BED_COORD) ==="
###############################################################################
# gtf2bed_vasa.py is used UNMODIFIED and shared with the mouse builder. The
# contig prefix needs no special handling because $GTF already carries it.
say "GTF: $(awk -F'\t' '$3=="gene"' "$GTF" | wc -l) gene features"
GENEBED="${BUILD}/genes_mixed_${BED_COORD}.bed"
"$PY" "$GTF2BED" "$GTF" "$GENEBED" --coord "$BED_COORD" || die "gtf2bed_vasa.py failed"
say "gene rows: $(wc -l < "$GENEBED")"

###############################################################################
say "=== STEP 2: validate the builder by reproducing the existing BED ==="
###############################################################################
# The file being replaced has no tracked builder, so this test is what
# establishes that its rules really are its rules. Compared as a SORTED SET:
# row order is irrelevant to bedtools intersect and the two writers differ in
# it (gtf2bed_vasa.py emits in GTF gene order, gtf_to_homemade_bed.py sorts).
if [ "$VALIDATE" = "yes" ] && [ -s "$V1_BED" ]; then
    REPRO="${BUILD}/repro_mixed_v1.bed"
    if [ "$BED_COORD" = "asis" ]; then cp "$GENEBED" "$REPRO"
    else "$PY" "$GTF2BED" "$GTF" "$REPRO" --coord asis || die "repro build failed"; fi
    sort "$V1_BED" > "${BUILD}/v1.sorted"
    sort "$REPRO"  > "${BUILD}/repro.sorted"
    if cmp -s "${BUILD}/v1.sorted" "${BUILD}/repro.sorted"; then
        say "    ok   reproduces the existing BED exactly ($(wc -l < "${BUILD}/v1.sorted") rows, sorted-set identical)"
    else
        echo "    FAIL does not reproduce the existing BED:"
        echo "      only in existing: $(comm -23 "${BUILD}/v1.sorted" "${BUILD}/repro.sorted" | wc -l)"
        echo "      only in rebuild : $(comm -13 "${BUILD}/v1.sorted" "${BUILD}/repro.sorted" | wc -l)"
        comm -3 "${BUILD}/v1.sorted" "${BUILD}/repro.sorted" | head -6
        die "builder does not match the file it replaces -- resolve before writing"
    fi
elif [ "$VALIDATE" = "yes" ]; then
    echo "    SKIP existing BED not found ($V1_BED) -- reproduction test not run"
fi

###############################################################################
say "=== STEP 3: tRNA rows, per species ==="
###############################################################################
#   chrom  start  end  strand  NAME  genelen  genestart  geneend
# For a tRNA the feature IS the gene, so genelen=end-start and genestart/geneend
# repeat start/end.
#
# NAME must be BARE -- no _Biotype_Label suffix, unlike every other row -- must
# contain the literal string "tRNA", and must contain no "_". All three are
# forced by countTables_2pickle_cellsSpliced.py:93, which does
#     Gene.replace('-','.') + '_tRNA'  if 'tRNA' in Gene
#     Biotype = Gene.rsplit('_')[-2]
# so `13.tRNA1495-ValAAC` becomes `13.tRNA1495.ValAAC_tRNA` -> Label=tRNA,
# Biotype=13.tRNA1495.ValAAC, which is the published value exactly. GtRNAdb's
# low-confidence `tRX-*` symbols contain no "tRNA" and no "_", and
# rsplit('_')[-2] on them raises IndexError -- pasting GtRNAdb's BED in
# unconverted would CRASH step 6, not degrade it. Hence the name_map rebuild.
#
# Contig names UCSC -> Ensembl -> species prefix:
#   chr13 -> 13 -> GRCh38_13 ; chrX/chrY -> X/Y ; chrM -> MT
#   hg38 chr1_KI270713v1_random -> KI270713.1   (v<N> -> .<N>)
#   mm10 chr1_GL456210_random   -> GL456210.1   (mm10 carries no version, so .1)
TRNABED="${BUILD}/trna_rows_mixed.bed"
: > "$TRNABED"

if [ "$WITH_TRNA" = "yes" ]; then
    for spec in $GTRNADB_SPECS; do
        IFS=: read -r gdir gpfx cpfx <<< "$spec"
        TARBALL="${BUILD}/${gpfx}-tRNAs.tar.gz"
        SRC="${BUILD}/gtrnadb_${gpfx}"
        BEDIN="${SRC}/${gpfx}-tRNAs.bed"
        NAMEMAP="${SRC}/${gpfx}-tRNAs_name_map.txt"

        if [ ! -s "$TARBALL" ]; then
            say "downloading ${GTRNADB_BASE}/${gdir}/${gpfx}-tRNAs.tar.gz"
            wget -q -O "${TARBALL}.tmp" "${GTRNADB_BASE}/${gdir}/${gpfx}-tRNAs.tar.gz" \
                || die "GtRNAdb download failed for ${gpfx}"
            mv "${TARBALL}.tmp" "$TARBALL"
        fi
        mkdir -p "$SRC"; tar xzf "$TARBALL" -C "$SRC"
        [ -s "$BEDIN" ]   || die "no BED in tarball: $BEDIN"
        [ -s "$NAMEMAP" ] || die "no name_map in tarball: $NAMEMAP"

        n_before=$(wc -l < "$TRNABED")
        awk -v set="$TRNA_SET" -v coord="$BED_COORD" -v cpfx="$cpfx" 'BEGIN{FS=OFS="\t"}
            NR==FNR { if (FNR>1) scan[$2]=$1; next }
            {
                gid=$4
                if (set=="high" && gid !~ /^tRNA-/) next
                sid=scan[gid]
                if (sid=="") { print "no name_map entry for " gid > "/dev/stderr"; bad++; next }
                i=index(sid, ".trna")
                if (i==0) { print "unparsable tRNAscan id " sid > "/dev/stderr"; bad++; next }
                uchr=substr(sid,1,i-1); num=substr(sid,i+5)
                if (uchr ~ /^chr[0-9]+$/ || uchr=="chrX" || uchr=="chrY") { chrom=substr(uchr,4) }
                else if (uchr=="chrM")                                    { chrom="MT" }
                else if (match(uchr, /(GL|JH|KI)[0-9]+(v[0-9]+)?/)) {
                    tok=substr(uchr,RSTART,RLENGTH)
                    if (tok ~ /v[0-9]+$/) sub(/v/,".",tok); else tok=tok".1"
                    chrom=tok
                }
                else { print "unmapped contig " uchr > "/dev/stderr"; bad++; next }
                n=split(gid,p,"-")
                if (n<3) { print "unparsable GtRNAdb id " gid > "/dev/stderr"; bad++; next }
                # NAME stays species-free (the authors chose one flat namespace);
                # only the CONTIG is prefixed, so it matches the STAR index.
                name = chrom ".tRNA" num "-" p[2] p[3]
                start=$2; end=$3                    # GtRNAdb BED is true 0-based
                if (coord=="asis") start=start+1     # match the authors 1-bp-high convention
                print cpfx chrom, start, end, $6, name, end-start, start, end
            }
            END { if (bad>0) { print "FAILED: " bad " rows unconverted" > "/dev/stderr"; exit 1 } }
        ' "$NAMEMAP" "$BEDIN" >> "$TRNABED" || die "tRNA conversion failed for ${gpfx}"
        say "    ${gpfx}: $(( $(wc -l < "$TRNABED") - n_before )) tRNA rows (prefix ${cpfx})"
    done
    say "tRNA rows: $(wc -l < "$TRNABED") total (TRNA_SET=$TRNA_SET)"

    # The one collision the header documents. Report it rather than hide it:
    # two loci sharing a name are merged by step 6 into a single count row.
    ncol=$(cut -f5 "$TRNABED" | sort | uniq -d | wc -l)
    say "    names shared between the two species: ${ncol}"
else
    say "tRNA rows: 0 (WITH_TRNA=no)"
fi

###############################################################################
say "=== STEP 4: validate the output ==="
###############################################################################
fail=0
chk() { if [ "$2" -eq 0 ]; then echo "    ok   $1"; else echo "    FAIL $1 ($2)"; fail=$((fail+1)); fi; }

ALL="${BUILD}/all_rows_mixed.bed"
cat "$GENEBED" "$TRNABED" > "$ALL"

chk "8 columns on every row"        "$(awk -F'\t' 'NF!=8' "$ALL" | wc -l)"
chk "strand is + or -"              "$(awk -F'\t' '$4!="+" && $4!="-"' "$ALL" | wc -l)"
chk "genelen == geneend-genestart"  "$(awk -F'\t' '$6!=($8-$7)' "$ALL" | wc -l)"
chk "feature within its gene span"  "$(awk -F'\t' '$2<$7 || $3>$8' "$ALL" | wc -l)"
chk "name has no whitespace"        "$(awk -F'\t' '$5 ~ /[ \t]/' "$ALL" | wc -l)"
chk "every contig is in the STAR index" \
    "$(cut -f1 "$ALL" | sort -u | comm -23 - <(sort -u "$CHRNAMES") | wc -l)"

# Zero-length rows: 1 bp GTF features (start==end) that the asis convention
# turns into empty intervals bedtools can never match. Expected under asis --
# they are part of what "reproduce the authors" means -- impossible under fix.
ndeg=$(awk -F'\t' '$2<0 || $3<=$2' "$ALL" | wc -l)
if [ "$BED_COORD" = "fix" ]; then
    chk "start>=0 and end>start" "$ndeg"
elif [ "$ndeg" -eq 0 ]; then
    echo "    ok   start>=0 and end>start"
else
    echo "    WARN $ndeg zero-length rows (expected under BED_COORD=asis; bedtools"
    echo "         cannot match them -- BED_COORD=fix makes them 1 bp intervals)"
fi

if [ "$WITH_TRNA" = "yes" ]; then
    chk "tRNA: every name contains 'tRNA'" "$(awk -F'\t' '$5 !~ /tRNA/' "$TRNABED" | wc -l)"
    chk "tRNA: no '_' in any name"         "$(awk -F'\t' '$5 ~ /_/'     "$TRNABED" | wc -l)"

    # The only check that proves the file cannot crash step 6: replay its three
    # parsing lines on every tRNA row.
    if "$PY" - "$TRNABED" <<'PY'
import sys
bad = n = 0
for line in open(sys.argv[1]):
    name = line.rstrip('\n').split('\t')[4]; n += 1
    g = name.replace('-', '.') + '_tRNA' if 'tRNA' in name else name
    try:
        biotype = g.rsplit('_')[-2]; label = g.rsplit('_')[-1]
    except IndexError:
        bad += 1; continue
    if label != 'tRNA' or '_' in biotype:
        bad += 1
print(f"    {'ok  ' if bad == 0 else 'FAIL'} step-6 name parse on all {n} tRNA rows"
      + (f" ({bad} bad)" if bad else ""))
sys.exit(1 if bad else 0)
PY
    then :; else fail=$((fail+1)); fi
fi

[ "$fail" -eq 0 ] || die "$fail validation check(s) failed -- not writing $OUT_BED"

###############################################################################
say "=== STEP 5: write ==="
###############################################################################
sort -k1,1 -k2,2n "$ALL" > "${OUT_BED}.tmp" && mv "${OUT_BED}.tmp" "$OUT_BED"
ngene=$(wc -l < "$GENEBED"); ntrna=$(wc -l < "$TRNABED"); nout=$(wc -l < "$OUT_BED")
[ "$nout" -eq "$((ngene + ntrna))" ] || die "row count mismatch: $ngene + $ntrna != $nout"
say "wrote $OUT_BED ($ngene gene rows + $ntrna tRNA rows = $nout)"

###############################################################################
say "=== STEP 6: provenance ==="
###############################################################################
PROV="${OUT_BED%.bed}.provenance.txt"
{
    echo "$(basename "$OUT_BED")"
    echo "built    : $(date -Is)"
    echo "by       : build_annotation_bed_mixed.sh + gtf2bed_vasa.py"
    echo "           (both tracked in code/I_Gene_expression/own_version/)"
    echo "host     : $(hostname)"
    echo
    echo "BUILT FROM PRIMARY SOURCES -- nothing inherited from an existing BED."
    echo
    echo "gene rows"
    echo "  gtf        : ${GTF}"
    echo "               (Ensembl 99, GRCh38 + GRCm38, contigs already prefixed"
    echo "                GRCh38_ / GRCm38_; the same file the STAR index was built from)"
    echo "  converter  : gtf2bed_vasa.py --coord ${BED_COORD}"
    echo "  rows       : ${ngene}"
    echo
    echo "tRNA rows"
    if [ "$WITH_TRNA" = "yes" ]; then
        for spec in $GTRNADB_SPECS; do
            IFS=: read -r gdir gpfx cpfx <<< "$spec"
            echo "  ${gpfx}       : ${GTRNADB_BASE}/${gdir}/${gpfx}-tRNAs.tar.gz  -> contigs ${cpfx}*"
        done
        echo "  set        : ${TRNA_SET}"
        echo "  rows       : ${ntrna}"
        echo "  naming     : bare tRNAscan-SE ids (<chrom>.tRNA<n>-<AA><Anticodon>),"
        echo "               NO species tag -- matches the published tables, which"
        echo "               carry no species marker. Contig column IS prefixed."
        echo "  collisions : $(cut -f5 "$TRNABED" | sort | uniq -d | wc -l) name(s) shared between the two species"
    else
        echo "  (none -- WITH_TRNA=no)"
    fi
    echo
    echo "coordinates : ${BED_COORD}"
    if [ "$BED_COORD" = "asis" ]; then
        echo "  GTF 1-based numbers emitted verbatim into a 0-based BED, so every"
        echo "  start is 1 bp high. This is the AUTHORS' convention, established by"
        echo "  measurement, not assumed: on 72,613 genes shared with GSM5369495 the"
        echo "  offset between short non-splicing biotypes and protein-coding is"
        echo "  +1.12%, where an asis/fix mismatch would show about -8%."
        echo "  BED_COORD=fix produces the corrected file."
    else
        echo "  True 0-based half-open. Correct, but NOT what the authors used --"
        echo "  use BED_COORD=asis to reproduce the published tables."
    fi
    echo
    echo "total rows  : ${nout}"
    echo "replaces    : $(basename "$V1_BED") (${V1_BED}), which has 0 tRNA rows"
    echo "consumed by : deal_with_singlemappers.sh / deal_with_multimappers.sh (step 5)"
    echo
    echo "md5"
    md5sum "$OUT_BED" | sed 's/^/  /'
    [ -s "$V1_BED" ] && md5sum "$V1_BED" | sed 's/^/  /'
} > "$PROV"
say "wrote $PROV"

say "=== DONE ==="
echo "REF_BED=${OUT_BED}"
