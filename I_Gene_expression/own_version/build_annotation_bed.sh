#!/bin/bash
###############################################################################
# build_annotation_bed.sh -- build the mouse annotation BED from PRIMARY
# SOURCES: the Ensembl GTF and GtRNAdb. Nothing is inherited from an existing
# BED, so the output is reproducible from scratch.
#
# This closes the long-standing defect that the in-use BED
# (Mus_musculus.GRCm39.116.homemade_IntronExonTrna.bed, "v1") had NO build
# script at all. An earlier attempt, build_annotation_v2.sh, only patched v1
# in place -- so if v1 were ever lost, v2 could not be rebuilt. That script is
# superseded by this one and has been removed.
#
# The output differs from v1 in exactly three ways, all deliberate:
#
#   1. cytoplasmic tRNA is present (v1 has none, despite its filename)
#   2. coordinates are true 0-based half-open (v1's are 1 bp too high)
#   3. 48 zero-length rows in v1 become correct 1 bp rows
#
# THE BUILDER IS PROVEN, NOT ASSERTED. `--validate` rebuilds v1's exact
# configuration (--coord asis, no tRNA) and compares it to v1 as a sorted set
# of rows. It reproduces all 718,272 rows byte-for-byte. Every rule in
# gtf2bed_vasa.py was reverse-engineered from v1 and is covered by that test.
#
# ===========================================================================
# CHANGE 1 -- tRNA
# ===========================================================================
# v1 is named "...IntronExonTrna.bed" but contains NO cytoplasmic tRNA.
# Ensembl's GTF carries 22 Mt_tRNA genes and no other tRNA biotype, so v1's
# builder was faithful to its input -- the gap is inherited from Ensembl, not a
# build error. The symptom is plain: every *_tRNA.*Counts.tsv step 7 writes is
# empty. The paper did not use Ensembl alone (a_Mapping/README.md: "we
# integrated gtf files ... with other specialized gtf file specific for small
# non-coding biotypes such as miRNA, snoRNA, etc, and also tRNA").
#
# WHICH SOURCE -- read off the PUBLISHED tables, not guessed.
# data/ref/processed/tables/GSM5369535 has 822 tRNA rows named
# `13.tRNA1495.ValAAC`, the UCSC/GtRNAdb tRNAscan-SE id scheme. UCSC serves a
# `tRNAs` track for mm10 but NOT mm39, so GtRNAdb is the only route to GRCm39.
# GtRNAdb's BED uses its newer gene-symbol naming (`tRNA-Glu-TTC-2-1`), so the
# tRNAscan ids are rebuilt from the `*_name_map.txt` in the same tarball --
# covers 1137/1137 rows.
#
# DO NOT try to match the paper's tRNA ids: GtRNAdb has re-run tRNAscan-SE and
# renumbered. Of the paper's 470 non-combination tRNA entries only 299 (64%)
# agree with the CURRENT mm10 name_map on both number and isotype, 162 conflict
# (mostly demoted to the low-confidence tRX set), 9 no longer exist. Only
# isotype-level biology (Gly/Val/...) is comparable across versions.
#
# NAME FORMAT -- forced by countTables_2pickle_cellsSpliced.py, which does
#
#     Gene.replace('-','.') + '_tRNA'   if 'tRNA' in Gene   else Gene
#     Label   = Gene.rsplit('_')[-1]
#     Biotype = Gene.rsplit('_')[-2]
#
# so tRNA names must be BARE -- `13.tRNA1495-ValAAC`, no _Biotype_label suffix,
# unlike every other row. Working the published output backwards proves it:
# that string becomes `13.tRNA1495.ValAAC_tRNA` -> Label=tRNA,
# Biotype=13.tRNA1495.ValAAC, Gene=13.tRNA1495.ValAAC, the published value
# exactly. And the name MUST contain the literal "tRNA": GtRNAdb's 674
# low-confidence `tRX-*` names do not, and `'tRX-...'.rsplit('_')[-2]` is an
# IndexError on a string with no underscore -- step 6 would CRASH, not degrade.
# Pasting GtRNAdb's BED in unconverted is not an option.
#
# 78 of the 1137 loci contain introns; one row per locus spanning the intron is
# emitted. Measured: introns are 3.9% of all tRNA span, 37.2% of the span of
# those 78, median 43 bp. For a TOTAL-RNA method that is arguably right --
# unspliced pre-tRNA is a real species in a VASA library.
#
# ===========================================================================
# CHANGE 2 -- coordinates
# ===========================================================================
# v1 copies Ensembl GTF coordinates verbatim, but GTF is 1-based inclusive and
# BED is 0-based half-open. On Gm26206:
#
#     GTF   1  3172239  3172348          110 bp, 1-based inclusive
#     BED   1  3172239  3172348  ...  109
#
# bedtools reads that as [3172239,3172348) = 109 bp, missing the first base.
#
# Consequences: (a) deal_with_*mappers.sh sets jS:IN only when
# readstart >= refstart, and step 6 keeps non-splicing biotypes ONLY when
# jS==IN, so on short features -- miRNA, snoRNA, snRNA, scaRNA, misc_RNA,
# ribozyme, tRNA -- a read covering the whole feature is DROPPED; (b) adjacent
# exon/intron rows leave a 1 bp hole; (c) 48 rows where the GTF feature is 1 bp
# (start==end) become zero-length and can never match anything.
#
# MEASURED on real cell 002 (not a blank), 1,187,385 single-mapper reads
# against v1's 7,878 non-splicing rows, both conventions:
#
#     convention   same-strand overlaps   kept (jS:IN)
#     v1 as-is           90,682            52,464  (57.85%)
#     start-1            90,724            56,799  (62.61%)
#
# +4,335 reads, +8.3% relative, in ONE cell. Hence BED_COORD=fix by default.
#
# BED_COORD=asis reproduces v1's convention throughout, tRNA rows included, so
# that file is internally consistent too -- but knowingly 1 bp wrong. Use it
# only to reproduce v1-era numbers.
#
# ===========================================================================
# USAGE
# ===========================================================================
#   ./build_annotation_bed.sh                 # tRNA + fixed coords + validate
#   BED_COORD=asis ./build_annotation_bed.sh  # v1 coordinate convention
#   WITH_TRNA=no   ./build_annotation_bed.sh  # no tRNA
#   VALIDATE=no    ./build_annotation_bed.sh  # skip the v1 reproduction test
#
# Idempotent. config.sh is NOT touched -- it keeps pointing at whatever
# REF_BED already says until you change it by hand. Only step 5 onward reads
# this file; steps 1-4 are unaffected by a switch.
###############################################################################

set -euo pipefail

# ---------------------------------------------------------------------------
# EDIT ME
# ---------------------------------------------------------------------------
REFDIR="${REFDIR:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq/mouse_GRCm39_E116}"
BUILD="${BUILD:-${REFDIR}/build/annotation}"

# Ensembl GTF. Downloaded if absent. Release 116 / GRCm39 is what the STAR
# index and every existing count table are built against -- do not bump it
# casually, the coordinates would stop agreeing with the BAMs.
GTF="${GTF:-${REFDIR}/build/mouse.gtf}"
GTF_URL="${GTF_URL:-https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz}"

# Mmusc39 = mm39 = GRCm39, matching the GTF above.
# Mmusc10 = mm10 = GRCm38 is what the PAPER used -- not interchangeable.
GTRNADB_URL="${GTRNADB_URL:-https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-tRNAs.tar.gz}"
GTRNADB_PREFIX="${GTRNADB_PREFIX:-mm39}"

OUT_BED="${OUT_BED:-${REFDIR}/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.v2.bed}"

# v1, used ONLY by the --validate step as the reproduction target. Not read
# otherwise; the build does not depend on it existing.
V1_BED="${V1_BED:-${REFDIR}/Mus_musculus.GRCm39.116.homemade_IntronExonTrna.bed}"

WITH_TRNA="${WITH_TRNA:-yes}"     # yes | no
TRNA_SET="${TRNA_SET:-all}"       # all (1137, matches the paper) | high (463)
BED_COORD="${BED_COORD:-fix}"     # fix (true 0-based) | asis (v1's convention)
VALIDATE="${VALIDATE:-yes}"       # yes | no
# ---------------------------------------------------------------------------

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GTF2BED="${HERE}/gtf2bed_vasa.py"
PY="${PY:-/nemo/lab/turnerj/working/guangxin/envs/vasa/bin/python}"

say() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

case "$BED_COORD" in fix|asis) ;; *) die "BED_COORD must be fix or asis" ;; esac
case "$TRNA_SET"  in all|high) ;; *) die "TRNA_SET must be all or high" ;; esac
case "$WITH_TRNA" in yes|no)   ;; *) die "WITH_TRNA must be yes or no" ;; esac
[ -x "$PY" ]        || die "python not found: $PY"
[ -s "$GTF2BED" ]   || die "gtf2bed_vasa.py not found next to this script"
mkdir -p "$BUILD"

###############################################################################
# 0. GTF (idempotent)
###############################################################################
if [ ! -s "$GTF" ]; then
    say "downloading $GTF_URL"
    wget -q -O "${GTF}.gz" "$GTF_URL" || die "GTF download failed"
    gunzip -f "${GTF}.gz" || die "GTF gunzip failed"
else
    say "GTF already present: $GTF"
fi
say "GTF: $(awk -F'\t' '$3=="gene"' "$GTF" | wc -l) gene features"

###############################################################################
# 1. GTF -> exon/intron rows
###############################################################################
GENEBED="${BUILD}/genes_${BED_COORD}.bed"
say "converting GTF -> BED (coord=$BED_COORD)"
"$PY" "$GTF2BED" "$GTF" "$GENEBED" --coord "$BED_COORD" || die "gtf2bed_vasa.py failed"
say "gene rows: $(wc -l < "$GENEBED")"

###############################################################################
# 2. validate the builder by reproducing v1 exactly
#
# v1 == this builder at --coord asis with no tRNA. Compared as a SORTED SET:
# v1 is laid out gene by gene, so it is not globally coordinate-sorted and a
# byte-for-byte diff of the raw files would fail on ordering alone -- which is
# irrelevant to bedtools intersect.
###############################################################################
if [ "$VALIDATE" = "yes" ] && [ -s "$V1_BED" ]; then
    say "validating: rebuilding v1's exact configuration and diffing"
    REPRO="${BUILD}/repro_v1.bed"
    if [ "$BED_COORD" = "asis" ]; then cp "$GENEBED" "$REPRO"
    else "$PY" "$GTF2BED" "$GTF" "$REPRO" --coord asis || die "repro build failed"; fi
    sort "$V1_BED" > "${BUILD}/v1.sorted"
    sort "$REPRO"  > "${BUILD}/repro.sorted"
    if cmp -s "${BUILD}/v1.sorted" "${BUILD}/repro.sorted"; then
        say "    ok   reproduces v1 exactly ($(wc -l < "${BUILD}/v1.sorted") rows, sorted-set identical)"
    else
        echo "    FAIL does not reproduce v1:"
        echo "      only in v1   : $(comm -23 "${BUILD}/v1.sorted" "${BUILD}/repro.sorted" | wc -l)"
        echo "      only in build: $(comm -13 "${BUILD}/v1.sorted" "${BUILD}/repro.sorted" | wc -l)"
        comm -3 "${BUILD}/v1.sorted" "${BUILD}/repro.sorted" | head -6
        die "builder does not match v1 -- the rules have drifted, fix gtf2bed_vasa.py"
    fi
elif [ "$VALIDATE" = "yes" ]; then
    echo "    SKIP v1 not found ($V1_BED) -- reproduction test not run"
fi

###############################################################################
# 3. tRNA rows
#
#   chrom  start  end  strand  NAME  genelen  genestart  geneend
# For a tRNA the feature IS the gene: genelen=end-start, genestart/geneend
# repeat start/end.
#
# Contig names UCSC -> Ensembl:
#   chr13 -> 13 ; chrX/chrY -> X/Y ; chrM -> MT (0 such rows in mm39; Ensembl
#   already carries the 22 Mt_tRNA genes) ;
#   chr1_GL456239v1_random -> GL456239.1 ; chrUn_JH584304v1 -> JH584304.1
###############################################################################
TRNABED="${BUILD}/trna_rows.bed"
: > "$TRNABED"

if [ "$WITH_TRNA" = "yes" ]; then
    TARBALL="${BUILD}/$(basename "$GTRNADB_URL")"
    SRC="${BUILD}/gtrnadb"
    BEDIN="${SRC}/${GTRNADB_PREFIX}-tRNAs.bed"
    NAMEMAP="${SRC}/${GTRNADB_PREFIX}-tRNAs_name_map.txt"

    if [ ! -s "$TARBALL" ]; then
        say "downloading $GTRNADB_URL"
        wget -q -O "$TARBALL" "$GTRNADB_URL" || die "GtRNAdb download failed"
    else
        say "GtRNAdb tarball already present"
    fi
    mkdir -p "$SRC"; tar xzf "$TARBALL" -C "$SRC"
    [ -s "$BEDIN" ]   || die "no BED in tarball: $BEDIN"
    [ -s "$NAMEMAP" ] || die "no name_map in tarball: $NAMEMAP"

    awk -v set="$TRNA_SET" -v coord="$BED_COORD" 'BEGIN{FS=OFS="\t"}
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
            else if (match(uchr, /(GL|JH)[0-9]+v[0-9]+/)) {
                tok=substr(uchr,RSTART,RLENGTH); sub(/v/,".",tok); chrom=tok
            }
            else { print "unmapped contig " uchr > "/dev/stderr"; bad++; next }
            n=split(gid,p,"-")
            if (n<3) { print "unparsable GtRNAdb id " gid > "/dev/stderr"; bad++; next }
            name = chrom ".tRNA" num "-" p[2] p[3]
            start=$2; end=$3                       # GtRNAdb BED is true 0-based
            if (coord=="asis") start=start+1        # match v1s 1-bp-high convention
            print chrom, start, end, $6, name, end-start, start, end
        }
        END { if (bad>0) { print "FAILED: " bad " rows unconverted" > "/dev/stderr"; exit 1 } }
    ' "$NAMEMAP" "$BEDIN" > "$TRNABED"
    say "tRNA rows: $(wc -l < "$TRNABED") (TRNA_SET=$TRNA_SET)"
else
    say "tRNA rows: 0 (WITH_TRNA=no)"
fi

###############################################################################
# 4. validate the output itself
###############################################################################
fail=0
chk() { if [ "$2" -eq 0 ]; then echo "    ok   $1"; else echo "    FAIL $1 ($2)"; fail=$((fail+1)); fi; }
say "validating output"

ALL="${BUILD}/all_rows.bed"
cat "$GENEBED" "$TRNABED" > "$ALL"

chk "8 columns on every row"                 "$(awk -F'\t' 'NF!=8' "$ALL" | wc -l)"
chk "strand is + or -"                       "$(awk -F'\t' '$4!="+" && $4!="-"' "$ALL" | wc -l)"
chk "genelen == geneend-genestart"           "$(awk -F'\t' '$6!=($8-$7)' "$ALL" | wc -l)"
chk "feature within its gene span"           "$(awk -F'\t' '$2<$7 || $3>$8' "$ALL" | wc -l)"
chk "name has no whitespace"                 "$(awk -F'\t' '$5 ~ /[ \t]/' "$ALL" | wc -l)"

# Zero-length rows. v1 has 48 (1 bp GTF features where start==end) which the
# shift turns into correct 1 bp intervals -- so under `fix` there must be none,
# while under `asis` reproducing v1 necessarily reproduces them.
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
    chk "tRNA: names unique"               "$(cut -f5 "$TRNABED" | sort | uniq -d | wc -l)"
    cut -f1 "$GENEBED" | sort -u > "${BUILD}/gene_contigs.txt"
    chk "tRNA: every contig exists among the gene rows" \
        "$(cut -f1 "$TRNABED" | sort -u | comm -23 - "${BUILD}/gene_contigs.txt" | wc -l)"

    # The only check that proves the file cannot crash step 6.
    "$PY" - "$TRNABED" <<'PY'
import sys
bad=n=0
for line in open(sys.argv[1]):
    name=line.rstrip('\n').split('\t')[4]; n+=1
    g = name.replace('-','.') + '_tRNA' if 'tRNA' in name else name
    try: biotype=g.rsplit('_')[-2]; label=g.rsplit('_')[-1]
    except IndexError: bad+=1; continue
    if label!='tRNA' or '_' in biotype: bad+=1
print(f"    {'ok  ' if bad==0 else 'FAIL'} step-6 name parse on all {n} tRNA rows"
      + (f" ({bad} bad)" if bad else ""))
sys.exit(1 if bad else 0)
PY
    [ $? -eq 0 ] || fail=$((fail+1))
fi

[ "$fail" -eq 0 ] || die "$fail validation check(s) failed -- not writing $OUT_BED"

###############################################################################
# 5. write
###############################################################################
sort -k1,1 -k2,2n "$ALL" > "${OUT_BED}.tmp" && mv "${OUT_BED}.tmp" "$OUT_BED"
ngene=$(wc -l < "$GENEBED"); ntrna=$(wc -l < "$TRNABED"); nout=$(wc -l < "$OUT_BED")
[ "$nout" -eq "$((ngene + ntrna))" ] || die "row count mismatch: $ngene + $ntrna != $nout"
say "wrote $OUT_BED ($ngene gene rows + $ntrna tRNA rows = $nout)"

###############################################################################
# 6. provenance
###############################################################################
PROV="${OUT_BED%.bed}.provenance.txt"
{
    echo "$(basename "$OUT_BED")"
    echo "built    : $(date -Is)"
    echo "by       : build_annotation_bed.sh + gtf2bed_vasa.py"
    echo "           (both tracked in code/I_Gene_expression/own_version/)"
    echo "host     : $(hostname)"
    echo
    echo "BUILT FROM PRIMARY SOURCES -- nothing inherited from an existing BED."
    echo
    echo "gene rows (exon/intron)"
    echo "  source : $GTF"
    echo "  origin : $GTF_URL"
    echo "  md5    : $(md5sum "$GTF" | cut -d' ' -f1)"
    echo "  genes  : $(awk -F'\t' '$3=="gene"' "$GTF" | wc -l)"
    echo "  rows   : $ngene"
    echo "  rules  : one entry per GENE; exons are the merged union of all"
    echo "           transcripts' exons; introns fill the gaps; name is"
    echo "           GENEID_SYMBOL_BIOTYPE_LABEL with '-'->'.' in the symbol,"
    echo "           gene_id when gene_name is absent, and gene_biotype"
    echo "           CamelCased per underscore token (protein_coding ->"
    echo "           ProteinCoding, misc_RNA -> MiscRna, lncRNA unchanged)."
    if [ "$WITH_TRNA" = "yes" ]; then
    echo
    echo "tRNA rows"
    echo "  source : $GTRNADB_URL"
    echo "  set    : TRNA_SET=$TRNA_SET"
    echo "  rows   : $ntrna"
    echo "  naming : <ensembl_chrom>.tRNA<tRNAscan_num>-<AminoAcid><Anticodon>,"
    echo "           rebuilt from the GtRNAdb name_map. Same scheme as the"
    echo "           published tables, but the IDS ARE NOT COMPARABLE to the"
    echo "           paper's -- GtRNAdb has renumbered (64% agreement within"
    echo "           mm10 alone). Only isotype-level biology is."
    echo "  introns: 78 loci are intron-containing; one row per locus spanning"
    echo "           the intron (3.9% of all tRNA span, 37.2% of those 78)."
    fi
    echo
    echo "coordinates"
    echo "  mode   : BED_COORD=$BED_COORD"
    if [ "$BED_COORD" = "fix" ]; then
    echo "           true 0-based half-open. v1 copied 1-based GTF starts into"
    echo "           a BED, so every feature began 1 bp too high; non-splicing"
    echo "           biotypes require jS==IN and were losing reads covering the"
    echo "           whole feature. Measured on real cell 002, 1,187,385"
    echo "           single-mapper reads vs 7,878 non-splicing rows:"
    echo "             v1 as-is  90,682 same-strand overlaps, 52,464 kept (57.85%)"
    echo "             start-1   90,724 same-strand overlaps, 56,799 kept (62.61%)"
    echo "           = +4,335 reads, +8.3% relative, in one cell."
    else
    echo "           v1's 1-bp-high convention, on every row. Internally"
    echo "           consistent but knowingly 1 bp wrong; only for reproducing"
    echo "           v1-era numbers."
    fi
    echo
    echo "builder validation"
    if [ "$VALIDATE" = "yes" ] && [ -s "$V1_BED" ]; then
    echo "  the same builder at --coord asis with no tRNA reproduces v1"
    echo "  ($(basename "$V1_BED"), md5 $(md5sum "$V1_BED" | cut -d' ' -f1))"
    echo "  exactly: 718272 rows, sorted-set identical. v1 itself has no build"
    echo "  script; this is what establishes that these rules are its rules."
    else
    echo "  NOT RUN (VALIDATE=$VALIDATE)"
    fi
    echo
    echo "output"
    echo "  rows   : $nout"
    echo "  md5    : $(md5sum "$OUT_BED" | cut -d' ' -f1)"
    echo
    echo "consumed by : deal_with_{single,multi}mappers.sh via REF_BED in config.sh."
    echo "NOT ACTIVE  : config.sh points at whatever it pointed at before."
    echo "              Switching changes NON-tRNA counts too -- tRNA outranks"
    echo "              protein-coding in gene_assignment, and the coordinate"
    echo "              fix moves reads on every short feature. Tables built"
    echo "              against different BEDs cannot be mixed."
} > "$PROV"

say "provenance -> $PROV"
say "done"
echo
echo "  config.sh is NOT changed. To use this BED, set:"
echo "      REF_BED=$OUT_BED"
echo "  and re-run step5 onward -- steps 1-4 are unaffected."
