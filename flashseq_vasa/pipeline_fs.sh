#!/bin/bash
###############################################################################
# pipeline_fs.sh -- push the ten FLASH-seq libraries through the VASA
# count-table path, so both protocols are quantified by the SAME code.
#
# HOW TO USE
#   1. Edit config.sh (or override FSV_* on the command line)
#   2. ./pipeline_fs.sh check                  <- verifies every path and tool
#   3. ./pipeline_fs.sh lut                    <- measure VASA's read lengths (once)
#      FSV_ARM=native  ./pipeline_fs.sh prep map assign
#      FSV_LIBS=<one>  FSV_ARM=native ./pipeline_fs.sh pickle1   <- x10, one per job
#      FSV_ARM=native  ./pipeline_fs.sh pickle_merge tables recon
#   or ./pipeline_fs.sh status                 <- how many files exist per stage
#
#   FSV_LIBS='ZHA8833A9' ./pipeline_fs.sh prep map    <- one library, for testing
#
# THE STAGES (each one's output is the next one's input)
#   lut          VASA's step-3 fastqs -> vasa_starinput_len_lut.txt (once, both arms)
#   prep         delivered R1         -> cells/<LIB>_cbc_noumi_R1.fq.gz
#   map          that                 -> cells/<LIB>_cbc_noumi_E99_Aligned.out.bam
#   assign       that                 -> cells/*_genes.bed.gz       (reads -> genes)
#   pickle1      ONE library's beds   -> pickle/<LIB>/<LIB>dict.pickle   (step 6)
#   pickle_merge all of those         -> <SAMPLE>.pickle.gz
#   tables       pickle               -> <SAMPLE>_*.tsv             (step 7)
#   recon        everything above     -> reconciliation.tsv + report
#
# pickle1 is per library because step 6 is the memory-hungry stage and its cost
# is linear in BED size -- see the stage's own header for the measurements and for
# why splitting it computes the same object rather than an approximation.
#
# WHY THIS IS A FORK AND NOT A PATCH
#   Nothing in code/I_Gene_expression/a_Mapping/ is edited -- it is published
#   code, and this repo's rule is that it gets explanatory comments, not logic
#   changes. Steps 6 and 7 are called from here exactly as they ship, with
#   different arguments (protocol=smartseq_noUMI, filt_unigenes=n,
#   CELLID_FROM=f). Step 5 uses own_version's two forks, i.e. the SAME scripts
#   that produced the VASA tables. The four things that genuinely had to be
#   written -- the R1-only single-end prep, the adapter trim, the vasalen hard
#   trim, and the reconciliation -- live here.
#
#   Every argument choice is justified in config.sh's header and in
#   NOUMI_PATH.md; none of them is a preference.
#
# RE-RUN SAFETY
#   Four upstream scripts end with a bare `gzip <file>`, which REFUSES to
#   replace an existing <file>.gz -- it prints "already exists; not overwritten",
#   exits 2, and none of the four checks the status. So re-running a stage over a
#   directory that still holds the previous run's output leaves the FRESH data
#   uncompressed beside the STALE .gz, and the next stage reads the stale one.
#   That cost a whole run on the VASA side on 2026-07-27. rm_stale is called
#   before every such stage, per library, and every stage additionally verifies
#   its outputs are NEWER than their inputs. Counting files is not proof.
###############################################################################

set -u -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# `sbatch pipeline_fs.sh` copies the script into a spool dir and runs the copy,
# so BASH_SOURCE points there and config.sh is absent; fall back to the
# directory sbatch was invoked from.
if [ ! -f "${HERE}/config.sh" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.sh" ]; then
    HERE="${SLURM_SUBMIT_DIR}"
fi
[ -f "${HERE}/config.sh" ] || { echo "ERROR: cannot find config.sh next to pipeline_fs.sh (looked in ${HERE})" >&2; exit 1; }
source "${HERE}/config.sh"

say()  { echo "[$(date '+%H:%M:%S')] $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }
rule() { echo "-------------------------------------------------------------"; }
rm_stale() { [ $# -gt 0 ] && rm -f "$@"; return 0; }

# Per-library stems. The literal '_cbc' is LOAD-BEARING: step 6 with
# CELLID_FROM=f does cellfile[:cellfile.index('_cbc')] and raises
# `ValueError: substring not found` without it. FLASH-seq has no cell barcode --
# the token is in the filename purely to satisfy that contract, so a filename is
# transformed rather than any read. '_noumi' is there to say so out loud.
stem()   { echo "${1}_cbc_noumi"; }
prep_fq(){ echo "${FSV_CELLDIR}/$(stem "$1")_R1.fq.gz"; }
bam()    { echo "${FSV_CELLDIR}/$(stem "$1")_E99_Aligned.out.bam"; }
sglbed() { echo "${FSV_CELLDIR}/$(stem "$1")_E99_Aligned.out.singlemappers_genes.bed.gz"; }
mltbed() { echo "${FSV_CELLDIR}/$(stem "$1")_E99_Aligned.out.nsorted.multimappers_genes.bed.gz"; }

# The delivered R1 for a library. The S-number differs per library and is not
# derivable, so it is globbed -- but exactly one match is required, because two
# matches would silently mean a lane or a re-run got mixed in.
#
# `find -L`, NOT bare `find`. $FSV_FASTQ is a SYMLINK (the delivery directory is
# .../guangxin.zhang/RN26038/<run>/fastq -> /nemo/stp/sequencing/inputs/...), and
# find does not follow a symlink at its start point unless told to: bare
# `find "$FSV_FASTQ" -maxdepth 1 -name ...` stats the link itself, descends into
# nothing, prints nothing, and EXITS 0. That is the worst shape of failure --
# it looks exactly like "the files are not there". Measured on cn079: bare find
# returned nothing while a shell glob resolved the file, which is why the login
# node's `ls` probe disagreed with the job's check. -L makes find traverse it.
find_r1() {
    local lib=$1 hits n
    hits=$(find -L "$FSV_FASTQ" -maxdepth 1 -name "${lib}_S*_R1_001.fastq.gz" | sort)
    n=$(printf '%s\n' "$hits" | grep -c . || true)
    [ "$n" -eq 1 ] || die "$lib: expected exactly 1 delivered R1, found $n"
    printf '%s\n' "$hits"
}

libs() { echo $FSV_LIBS; }
nlibs() { libs | wc -w; }

###############################################################################
# check
###############################################################################
stage_check() {
    local bad=0
    rule; say "checking configuration"; rule
    echo "arm        : ${FSV_ARM}"
    echo "sample     : ${FSV_SAMPLE}"
    echo "libraries  : $(nlibs) -- $(libs)"
    echo "outdir     : ${FSV_OUTDIR}"
    echo "tables to  : ${FSV_TABLES}"
    echo "protocol   : ${FSV_PROTOCOL}  stranded=${FSV_STRANDED}  cellid=${FSV_CELLID_FROM}  filt_unigenes=${FSV_FILT_UNIGENES}"
    echo

    case "$FSV_ARM" in
        native|vasalen) ;;
        *) echo "  PROBLEM FSV_ARM='${FSV_ARM}' -- must be 'native' or 'vasalen'"; bad=1 ;;
    esac
    # These four are what make the run a method-matched comparison rather than a
    # second method. Wrong values do not crash -- they quietly halve or reshape
    # every figure -- so they are checked here, loudly.
    [ "$FSV_STRANDED"      = "n" ] || { echo "  PROBLEM stranded='${FSV_STRANDED}' -- FLASH-seq is unstranded; y would halve every biotype figure"; bad=1; }
    [ "$FSV_PROTOCOL"      = "smartseq_noUMI" ] || { echo "  PROBLEM protocol='${FSV_PROTOCOL}' -- must be smartseq_noUMI (vasa raises KeyError 'RX' on these read names)"; bad=1; }
    [ "$FSV_CELLID_FROM"   = "f" ] || { echo "  PROBLEM cellid_from='${FSV_CELLID_FROM}' -- 'r' needs ;SM:<lib> injected into every read name"; bad=1; }
    [ "$FSV_FILT_UNIGENES" = "n" ] || { echo "  PROBLEM filt_unigenes='${FSV_FILT_UNIGENES}' -- at y the threshold is max(5,1% of cols) = 5 of 10 = 50%"; bad=1; }
    # The '_cbc' contract, tested rather than trusted.
    case "$(stem ZHA8833A9)" in
        *_cbc*) echo "  OK   per-library stem contains '_cbc' ($(stem ZHA8833A9))" ;;
        *) echo "  PROBLEM stem '$(stem ZHA8833A9)' has no '_cbc' -- step 6 mode f will raise ValueError"; bad=1 ;;
    esac

    # -- input fastqs --
    local lib r1
    for lib in $(libs); do
        if r1=$(find_r1 "$lib" 2>/dev/null); then
            printf "  OK   %-11s %s\n" "$lib" "$(basename "$r1")"
        else
            echo "  MISS $lib -- no unique ${lib}_S*_R1_001.fastq.gz under $FSV_FASTQ"; bad=1
        fi
    done

    # -- scratch --
    if mkdir -p "$FSV_CELLDIR" "$FSV_LOGDIR" 2>/dev/null && [ -w "$FSV_OUTDIR" ]; then
        echo "  OK   $FSV_OUTDIR"
    else
        echo "  MISS $FSV_OUTDIR (not creatable/writable)"; bad=1
    fi
    mkdir -p "$FSV_RES" "$FSV_TABLES" 2>/dev/null || true
    [ -w "$FSV_RES" ]    && echo "  OK   $FSV_RES"    || { echo "  MISS $FSV_RES";    bad=1; }
    [ -w "$FSV_TABLES" ] && echo "  OK   $FSV_TABLES" || { echo "  MISS $FSV_TABLES"; bad=1; }

    # -- references --
    [ -s "$FSV_REF_BED" ] && echo "  OK   $FSV_REF_BED" || { echo "  MISS $FSV_REF_BED"; bad=1; }
    # Test for SA and Genome, not the directory: genomeGenerate creates the
    # directory up front, so a failed or still-building index would otherwise
    # pass check and blow up an hour into mapping.
    if [ -s "${FSV_STAR_INDEX}/SA" ] && [ -s "${FSV_STAR_INDEX}/Genome" ]; then
        echo "  OK   $FSV_STAR_INDEX"
        local ovh; ovh=$(awk '$1=="sjdbOverhang"{print $2; exit}' "${FSV_STAR_INDEX}/genomeParameters.txt" 2>/dev/null || true)
        # 151 nt reads need sjdbOverhang >= 150. Larger is harmless (STAR stores
        # flanking sequence it never uses); smaller silently costs
        # junction-spanning sensitivity, which is the one failure worth catching.
        case "$ovh" in
            ''|*[!0-9]*) echo "       (could not read sjdbOverhang)" ;;
            *) if [ "$ovh" -ge 150 ]; then echo "       sjdbOverhang=$ovh (needs >= 150 for 151 nt reads) -- OK"
               else echo "  PROBLEM sjdbOverhang=$ovh < 150"; bad=1; fi ;;
        esac
    else
        echo "  MISS ${FSV_STAR_INDEX}/SA -- index absent or still building; there is NO build script for it"; bad=1
    fi

    # -- scripts --
    local s
    for s in "$FSV_ASSIGN_SINGLE_SH" "$FSV_ASSIGN_MULTI_SH"; do
        [ -x "$s" ] && echo "  OK   $s" || { echo "  MISS/NOT-EXEC $s"; bad=1; }
    done
    for s in countTables_2pickle_cellsSpliced.py countTables_fromPickle.py; do
        [ -x "${FSV_VASA_SCRIPTS}/$s" ] && echo "  OK   ${FSV_VASA_SCRIPTS}/$s" \
            || { echo "  MISS/NOT-EXEC ${FSV_VASA_SCRIPTS}/$s"; bad=1; }
    done
    for s in "$FSV_TRIM_VASALEN" "$FSV_MEASURE_LEN" "$FSV_BUILD_TABLES" "$FSV_RECON" "$FSV_MERGE_PICKLES"; do
        [ -e "$s" ] && echo "  OK   $s" || { echo "  MISS $s"; bad=1; }
    done

    # -- tools --
    local t
    for t in "${FSV_P2STAR}/STAR" "${FSV_P2SAMTOOLS}/samtools" "${FSV_P2BEDTOOLS}/bedtools" "$FSV_CUTADAPT" "$FSV_PYTHON"; do
        [ -x "$t" ] && echo "  OK   $t" || { echo "  MISS $t"; bad=1; }
    done
    # pandas MUST be 2.x -- countTables_fromPickle.py calls DataFrame.applymap,
    # removed in pandas 3.0. This is the check that catches a silently upgraded env.
    if "$FSV_PYTHON" -c 'import pandas,sys; sys.exit(0 if hasattr(pandas.DataFrame,"applymap") else 1)' 2>/dev/null; then
        echo "  OK   pandas $("$FSV_PYTHON" -c 'import pandas;print(pandas.__version__)') has DataFrame.applymap"
    else
        echo "  PROBLEM pandas in $FSV_CONDA_ENV has no DataFrame.applymap -- step 7 will fail"; bad=1
    fi

    # -- the vasalen arm's LUT --
    if [ "$FSV_ARM" = vasalen ]; then
        [ -s "$FSV_VASALEN_LUT" ] && echo "  OK   $FSV_VASALEN_LUT ($(wc -l < "$FSV_VASALEN_LUT") entries)" \
            || { echo "  MISS $FSV_VASALEN_LUT -- run './pipeline_fs.sh lut' first"; bad=1; }
    fi

    # -- the saved TrimGalore reports the adapter check compares against --
    local nrep; nrep=$(find "$FSV_TG_REPORTS" -name '*_trimmed_1.fastq.gz_trimming_report.txt' 2>/dev/null | wc -l)
    echo "  $([ "$nrep" -ge 1 ] && echo OK || echo WARN) $nrep TrimGalore R1 reports under $FSV_TG_REPORTS"

    rule
    [ $bad -eq 0 ] && say "ALL CHECKS PASSED" || die "fix the MISS/PROBLEM lines above"
}

###############################################################################
# lut -- measure VASA's read-length distribution. Once; both arms read it.
###############################################################################
stage_lut() {
    say "lut: measuring VASA's STAR-input read-length distribution"
    eval "$FSV_CONDA_ACTIVATE"
    mkdir -p "$FSV_RES"
    "$FSV_PYTHON" "$FSV_MEASURE_LEN" \
        "${FSV_VASA_OUT}/cells" "$FSV_VASA_SAMPLE" \
        "$FSV_VASA_REALCELLS" "$FSV_VASA_BLANKS" \
        "$FSV_VASALEN_LUT" "${FSV_RES}/vasa_readlen_report.txt" \
        || die "measure_vasa_readlen.py failed"
    say "lut done: $FSV_VASALEN_LUT"
}

###############################################################################
# prep -- delivered R1 -> the fastq STAR will read
#
# native : adapter+quality trim only, natural length
# vasalen: the same, then hard-trimmed per read to a draw from VASA's LUT
#
# R1 ONLY. Paired-end cannot reach step 6 (bamtobed's /1,/2 suffix lands on the
# end of the nM value and step 6 int()s it -- 100% of 1,033,528 PE rows
# unparseable, measured). R1-only also keeps the unit of measurement identical to
# VASA's one-biological-read-per-fragment, and is the same choice
# code/flashseq/05_rrna_bwa.sh made for the rRNA leg.
###############################################################################
do_prep() {
    local lib=$1 r1 out log
    r1=$(find_r1 "$lib")
    out=$(prep_fq "$lib")
    log="${FSV_LOGDIR}/prep_${lib}.log"
    rm_stale "$out" "${out%.gz}"

    local trimmed="$out"
    [ "$FSV_ARM" = vasalen ] && trimmed="${FSV_CELLDIR}/$(stem "$lib")_R1.native.fq.gz"
    rm_stale "$trimmed"

    # The cutadapt call TrimGalore itself issued in the nf-core run, taken from
    # results/trimgalore/*_trimming_report.txt. One recipe, reused -- not a
    # second trimming method. Single-end here (TrimGalore ran --paired), which is
    # the one deliberate difference and is stated in config.sh.
    "$FSV_CUTADAPT" -j 4 -e "$FSV_TRIM_ERR" -q "$FSV_TRIM_Q" -O "$FSV_TRIM_OVERLAP" \
        -a "$FSV_TRIM_ADAPTER" -m "$FSV_TRIM_MINLEN" \
        -o "$trimmed" "$r1" > "${FSV_LOGDIR}/cutadapt_${lib}.log" 2>&1 \
        || { echo "  FAILED cutadapt: $lib"; return 1; }

    if [ "$FSV_ARM" = vasalen ]; then
        "$FSV_TRIM_VASALEN" "$trimmed" "$FSV_VASALEN_LUT" "$out" \
            "$FSV_VASALEN_SEED" "${FSV_LOGDIR}/vasalen_${lib}.tsv" \
            > "$log" 2>&1 || { echo "  FAILED vasalen trim: $lib"; return 1; }
        # The native intermediate is a full-size fastq on a filesystem at 83%.
        # It is reproducible from the delivered R1 by the cutadapt call above, so
        # it is not kept.
        rm -f "$trimmed"
    else
        : > "$log"
    fi

    # Record counts at every point, so `recon` never has to re-derive them.
    local n_raw n_out
    # cutadapt prints "Total reads processed:  38,710,566" -- thousands
    # separators and all. Take the whole comma-formatted number off the END of
    # the line and strip the commas from it; do NOT put ',' in the field
    # separator, which would split 38,710,566 into three fields and silently
    # record 38 as the chain-start count. recon consumes this value rather than
    # re-deriving it, so a truncated number would corrupt the reconciliation.
    n_raw=$(awk '/^Total reads processed:/{n=$NF; gsub(/,/,"",n); print n; exit}' "${FSV_LOGDIR}/cutadapt_${lib}.log")
    n_out=$(( $(zcat "$out" | wc -l) / 4 ))
    printf 'library\t%s\narm\t%s\nfastq_records\t%s\nprep_records\t%s\n' \
        "$lib" "$FSV_ARM" "${n_raw:-NA}" "$n_out" > "${FSV_LOGDIR}/prep_${lib}.tsv"
}

stage_prep() {
    say "prep: $(nlibs) libraries, arm=${FSV_ARM}, ${FSV_NCORES} at a time"
    mkdir -p "$FSV_CELLDIR" "$FSV_LOGDIR"
    if [ "$FSV_ARM" = vasalen ]; then
        [ -s "$FSV_VASALEN_LUT" ] || die "no LUT at $FSV_VASALEN_LUT -- run './pipeline_fs.sh lut'"
    fi
    export -f do_prep say die rm_stale stem prep_fq find_r1
    export FSV_ARM FSV_CELLDIR FSV_LOGDIR FSV_FASTQ FSV_CUTADAPT FSV_TRIM_ADAPTER
    export FSV_TRIM_Q FSV_TRIM_OVERLAP FSV_TRIM_ERR FSV_TRIM_MINLEN
    export FSV_TRIM_VASALEN FSV_VASALEN_LUT FSV_VASALEN_SEED FSV_LIBS
    # Workers run in a fresh `bash -c`, which inherits neither functions nor
    # unexported variables -- forgetting one of these cost a whole step-5 launch
    # on the VASA side ("rm_stale: command not found", and the ASSIGN_*_SH empty
    # besides). Every function and variable do_prep touches is above.
    libs | tr ' ' '\n' | xargs -P "$FSV_NCORES" -I{} bash -c 'do_prep {}' _

    # Counting files is not proof. Verify per library that the output exists and
    # is NEWER than the delivered R1, and fail loudly.
    local bad=0 lib out r1
    for lib in $(libs); do
        out=$(prep_fq "$lib"); r1=$(find_r1 "$lib")
        if [ ! -s "$out" ]; then echo "  BAD prep: $lib -- missing $(basename "$out")"; bad=$((bad+1))
        elif [ "$out" -ot "$r1" ]; then echo "  BAD prep: $lib -- output OLDER than the delivered R1"; bad=$((bad+1)); fi
    done
    [ "$bad" -eq 0 ] || die "prep FAILED for $bad libraries -- do NOT run map"

    # Did we reproduce TrimGalore's adapter rate? Both run over the whole
    # library, so they should agree to well under a point. A mismatch means the
    # reproduction is wrong and every downstream number would be wrong with it.
    # Checked on the native arm only: the vasalen arm's cutadapt call is the same
    # one, and its hard trim happens after.
    if [ "$FSV_ARM" = native ]; then
        say "prep: adapter-rate reproduction check against the saved TrimGalore reports"
        local rep ref our
        for lib in $(libs); do
            rep="${FSV_TG_REPORTS}/${lib}_trimmed_1.fastq.gz_trimming_report.txt"
            [ -s "$rep" ] || { echo "  SKIP $lib -- no saved report"; continue; }
            ref=$(awk -F'[()%]' '/^Reads with adapters:/{print $2; exit}' "$rep")
            our=$(awk -F'[()%]' '/^Reads with adapters:/{print $2; exit}' "${FSV_LOGDIR}/cutadapt_${lib}.log")
            [ -n "$our" ] || our=$(awk -F'[()%]' '/Read 1 with adapter:/{print $2; exit}' "${FSV_LOGDIR}/cutadapt_${lib}.log")
            awk -v a="$ref" -v b="$our" -v t="$FSV_TRIM_TOL" -v l="$lib" \
                'BEGIN{d=a-b; if(d<0)d=-d; printf "  %-11s TrimGalore %s%%  ours %s%%  diff %.2f pt %s\n", l, a, b, d, (d<=t?"OK":"MISMATCH"); exit !(a!="" && b!="" && d<=t)}' \
                || die "$lib: cutadapt reproduction disagrees with the TrimGalore report beyond ${FSV_TRIM_TOL} pt"
        done
    fi
    say "prep done: $(find "$FSV_CELLDIR" -name '*_cbc_noumi_R1.fq.gz' | wc -l) fastqs"
}

###############################################################################
# map -- STAR, single-end, same index and same flag set as the VASA side
#
# --genomeLoad is deliberately left at its default NoSharedMemory, as upstream
# map_star.sh does. It costs an index load per library but cannot leave 27 GB
# stranded in a node's shm if the job dies. Libraries are mapped SERIALLY with
# STAR threaded internally, because two concurrent STARs would each load the
# index.
###############################################################################
stage_map() {
    say "map: STAR over $(nlibs) libraries, ${FSV_STAR_THREADS} threads each, arm=${FSV_ARM}"
    eval "$FSV_ML_STAR"
    export TMPDIR="${FSV_OUTDIR}/tmp"; mkdir -p "$TMPDIR"
    local i=0 n lib fq pfx
    n=$(nlibs)
    for lib in $(libs); do
        i=$((i+1))
        fq=$(prep_fq "$lib")
        pfx="${FSV_CELLDIR}/$(stem "$lib")_E99_"
        [ -s "$fq" ] || { echo "  SKIP (no prep fastq): $lib"; continue; }
        say "  mapping $i/$n  $lib"
        # Flag set copied from upstream map_star.sh / the A9 dry run. --outSAMattributes
        # All is required: deal_with_*.sh reads NH and nM off the SAM text.
        "${FSV_P2STAR}/STAR" --runThreadN "$FSV_STAR_THREADS" --genomeDir "$FSV_STAR_INDEX" \
            --readFilesIn "$fq" --readFilesCommand zcat \
            --outFilterMultimapNmax 20 --outSAMunmapped Within \
            --outSAMtype BAM Unsorted --outSAMattributes All \
            --outFileNamePrefix "$pfx" > /dev/null \
            || echo "  FAILED star: $lib"
        rm -rf "${pfx}_STARtmp" "${pfx}Log.progress.out"
        mv -f "${pfx}Log.final.out" "${pfx}Log.final.txt" 2>/dev/null
        mv -f "${pfx}Log.out"       "${pfx}Log.txt"       2>/dev/null
        grep -E 'input reads|Uniquely mapped reads number|multiple loci|too short' \
            "${pfx}Log.final.txt" | sed 's/^/      /'
    done
    local bad=0
    for lib in $(libs); do
        local b; b=$(bam "$lib")
        if [ ! -s "$b" ]; then echo "  BAD map: $lib -- missing BAM"; bad=$((bad+1))
        elif [ "$b" -ot "$(prep_fq "$lib")" ]; then echo "  BAD map: $lib -- BAM OLDER than its fastq"; bad=$((bad+1)); fi
    done
    [ "$bad" -eq 0 ] || die "map FAILED for $bad libraries -- do NOT run assign"
    say "map done: $(find "$FSV_CELLDIR" -name '*_E99_Aligned.out.bam' | wc -l) BAMs"
}

###############################################################################
# assign -- reads to genes, own_version's two forks, stranded=n
#
# Both forks, not just the single one: step 6's per-library glob is
# cell + '*_genes.bed.gz', which matches the multimapper BED too, so multimappers
# are part of the measurement. own_version rather than upstream on BOTH so that
# the NH:i:10..19 rule is the same on the two sides of the comparison.
###############################################################################
do_assign() {
    local lib=$1 b s m stem_path
    b=$(bam "$lib")
    [ -s "$b" ] || { echo "  SKIP (no bam): $lib"; return 0; }
    s=$(sglbed "$lib"); m=$(mltbed "$lib")
    # Both scripts end with a bare `gzip` -- clear their targets. This one
    # matters most: step 6 GLOBS *_genes.bed.gz, so a stale file here is silently
    # folded into the count tables.
    rm_stale "$s" "${s%.gz}" "$m" "${m%.gz}"
    "$FSV_ASSIGN_SINGLE_SH" "$b" "$FSV_REF_BED" "$FSV_STRANDED" \
        "$FSV_P2SAMTOOLS" "$FSV_P2BEDTOOLS" > "${FSV_LOGDIR}/assign_single_${lib}.log" 2>&1 \
        || echo "  FAILED single: $lib"
    "$FSV_ASSIGN_MULTI_SH" "$b" "$FSV_REF_BED" "$FSV_STRANDED" \
        "$FSV_P2SAMTOOLS" "$FSV_P2BEDTOOLS" > "${FSV_LOGDIR}/assign_multi_${lib}.log" 2>&1 \
        || echo "  FAILED multi: $lib"
    # BED row counts, recorded now so `recon` does not re-decompress 10 GB.
    printf 'library\t%s\nsingle_bed_rows\t%s\nmulti_bed_rows\t%s\n' "$lib" \
        "$( [ -s "$s" ] && zcat "$s" | wc -l || echo 0 )" \
        "$( [ -s "$m" ] && zcat "$m" | wc -l || echo 0 )" \
        > "${FSV_LOGDIR}/assign_${lib}.tsv"
}

stage_assign() {
    say "assign: $(nlibs) libraries, stranded=${FSV_STRANDED}, ${FSV_NCORES} at a time"
    eval "$FSV_ML_BED"
    # deal_with_multimappers.sh ends with `sort -k4`, which spills to $TMPDIR on
    # large inputs. Left at the default it lands on node-local /tmp.
    export TMPDIR="${FSV_OUTDIR}/tmp"; mkdir -p "$TMPDIR"
    export -f do_assign say die rm_stale stem bam sglbed mltbed
    export FSV_CELLDIR FSV_LOGDIR FSV_REF_BED FSV_STRANDED FSV_P2SAMTOOLS FSV_P2BEDTOOLS
    export FSV_ASSIGN_SINGLE_SH FSV_ASSIGN_MULTI_SH FSV_LIBS
    libs | tr ' ' '\n' | xargs -P "$FSV_NCORES" -I{} bash -c 'do_assign {}' _

    local bad=0 lib
    for lib in $(libs); do
        local b; b=$(bam "$lib")
        for out in "$(sglbed "$lib")" "$(mltbed "$lib")"; do
            if [ ! -s "$out" ]; then echo "  BAD assign: $lib -- missing $(basename "$out")"; bad=$((bad+1))
            elif [ "$out" -ot "$b" ]; then echo "  BAD assign: $lib -- $(basename "$out") OLDER than its BAM"; bad=$((bad+1)); fi
        done
    done
    [ "$bad" -eq 0 ] || die "assign FAILED: $bad missing/stale outputs -- do NOT run pickle"

    # THE nM CONTRACT, on the real BEDs. Paired-end input produces 'nM:0/2',
    # which step 6 int()s and dies on. We map R1 only so this must be clean --
    # checked here rather than discovered an hour into step 6, which is where it
    # would otherwise surface. Costs seconds; the 200 sampled rows are a
    # contract check, not a statistic.
    say "assign: nM-field contract check (step 6 int()s this field)"
    for lib in $(libs); do
        local nbad
        nbad=$(zcat "$(sglbed "$lib")" | awk 'NR%997==1' | head -200 \
               | awk -F'\t' '{n=$7; sub(/.*;nM:/,"",n); sub(/;jS:.*/,"",n); if (n !~ /^[0-9]+$/) c++} END{print c+0}')
        [ "$nbad" -eq 0 ] || die "$lib: $nbad of 200 sampled BED rows have an unparseable nM field -- paired-end contamination"
    done
    say "assign done: $(find "$FSV_CELLDIR" -name '*singlemappers_genes.bed.gz' | wc -l) single, $(find "$FSV_CELLDIR" -name '*multimappers_genes.bed.gz' | wc -l) multi"
}

###############################################################################
# pickle -- step 6, upstream script, run ONE LIBRARY AT A TIME, then merged
#
# WHY THIS IS SPLIT, and why the split is exact rather than an approximation
# ------------------------------------------------------------------------
# Upstream step 6 runs once over the whole cells/ folder with an 8-wide pool.
# That shape does not fit ten FLASH-seq libraries. Sizing from the VASA side,
# where the same script was measured at 69 MB of *_genes.bed.gz -> 28m36s /
# 3.24 GB and 254 MB -> 1h47m27s / 12.22 GB (45.6 GB of RSS per GB of BED,
# scaling exponent 1.02 on both time and memory, i.e. linear): a full FLASH-seq
# library yields roughly 870 MB of BED (640 MB single + 230 MB multi,
# extrapolated from the stride-64 dry run at 32 bytes per uniquely mapped read),
# so one library costs ~5-6 h and ~40 GB. Eight of those concurrently in one
# process is ~250 GB in a single node, and a failure anywhere loses the whole
# run.
#
# get_cellDict(cell) is per-library independent -- it globs cell + '*_genes.bed.gz'
# and shares no state with any other cell -- and the parent does nothing but
# collect: `gcnt[cell] = cnt` for each result, then `pd.DataFrame(gcnt)`. So
# running the script once per library and merging the dicts computes the SAME
# object, not an approximation of it. merge_pickles.py does only what the
# parent's tail does, and asserts the properties that make it equivalent
# (disjoint column keys, no cross-library gene collision handling needed because
# DataFrame construction aligns on the union of row keys). Column ORDER cannot
# matter: step 7's first act is `cntdf = cntdf[sorted(cntdf.columns)]`.
#
# The cell id must come out the same as it would have from a single run, so each
# per-library working directory contains a folder literally named 'cells' holding
# symlinks to that library's two BEDs, and step 6 is invoked from its parent with
# the relative name 'cells'. That makes the id 'cells/<LIB>' for every library,
# exactly as a combined run would have produced.
#
#   ./pipeline_fs.sh pickle1              # every library, serially (rarely wanted)
#   FSV_LIBS=ZHA8833A9 ./pipeline_fs.sh pickle1   # one -- this is the array unit
#   ./pipeline_fs.sh pickle_merge         # combine into <SAMPLE>.pickle.gz
###############################################################################
pickle_dir() { echo "${FSV_OUTDIR}/pickle/$1"; }

do_pickle1() {
    local lib=$1 d s m
    d=$(pickle_dir "$lib")
    s=$(sglbed "$lib"); m=$(mltbed "$lib")
    [ -s "$s" ] || die "$lib: no singlemapper BED -- run assign first"
    rm -rf "$d"; mkdir -p "$d/cells"
    # Symlinks, not copies: 870 MB per library and the BEDs are read-only input.
    ln -sf "$s" "$d/cells/$(basename "$s")"
    [ -s "$m" ] && ln -sf "$m" "$d/cells/$(basename "$m")"
    say "pickle1 $lib: $(du -Lsh "$d/cells" | cut -f1) of BED"
    # Step 6 ends with a bare `gzip` on <output>.pickle, which refuses to
    # clobber; the rm -rf above already guarantees a clean directory.
    ( cd "$d" && "${FSV_VASA_SCRIPTS}/countTables_2pickle_cellsSpliced.py" \
        "cells" "$lib" "$FSV_PROTOCOL" "$FSV_CELLID_FROM" ) || die "$lib: step 6 failed"
    [ -s "$d/${lib}dict.pickle" ] || die "$lib: step 6 wrote no dict.pickle"
    # du, not `ls -lh | awk '{print $5}'` -- the group here is "domain users",
    # so the space in it shifts every ls field and $5 prints "users", not a size.
    say "pickle1 $lib done: $(du -Lh "$d/${lib}dict.pickle" | cut -f1)"
}

stage_pickle1() {
    eval "$FSV_CONDA_ACTIVATE"
    local lib
    for lib in $(libs); do do_pickle1 "$lib"; done
}

stage_pickle_merge() {
    say "pickle_merge: combining $(nlibs) per-library dicts into ${FSV_SAMPLE}.pickle.gz"
    eval "$FSV_CONDA_ACTIVATE"
    rm_stale "${FSV_OUTDIR}/${FSV_SAMPLE}.pickle.gz" "${FSV_OUTDIR}/${FSV_SAMPLE}.pickle" \
             "${FSV_OUTDIR}/${FSV_SAMPLE}dict.pickle"
    local args=()
    local lib
    for lib in $(libs); do
        local p; p="$(pickle_dir "$lib")/${lib}dict.pickle"
        [ -s "$p" ] || die "$lib: missing $p -- run pickle1 for it first"
        args+=("$p")
    done
    "$FSV_PYTHON" "$FSV_MERGE_PICKLES" "${FSV_OUTDIR}/${FSV_SAMPLE}" "${args[@]}" \
        || die "merge_pickles.py failed"
    [ -s "${FSV_OUTDIR}/${FSV_SAMPLE}.pickle.gz" ] || die "merge wrote no pickle.gz"
    say "pickle_merge done: ${FSV_OUTDIR}/${FSV_SAMPLE}.pickle.gz"
}

###############################################################################
# tables -- step 7, upstream, protocol=smartseq_noUMI, filt_unigenes=n
###############################################################################
stage_tables() {
    say "tables: step 7, protocol=${FSV_PROTOCOL}, filt_unigenes=${FSV_FILT_UNIGENES}"
    eval "$FSV_CONDA_ACTIVATE"
    ( cd "$FSV_OUTDIR" && "${FSV_VASA_SCRIPTS}/countTables_fromPickle.py" \
        "${FSV_SAMPLE}.pickle.gz" "$FSV_SAMPLE" "$FSV_PROTOCOL" "$FSV_FILT_UNIGENES" ) \
        || die "step 7 failed"
    local nl; nl=$(wc -l < "${FSV_OUTDIR}/${FSV_SAMPLE}_mapStats.log")
    say "  mapStats.log: $nl lines (a complete VASA log is 21 -- CLAUDE.md's 22 is wrong)"
    sed 's/^/      /' "${FSV_OUTDIR}/${FSV_SAMPLE}_mapStats.log"
    say "  tables written: $(find "$FSV_OUTDIR" -maxdepth 1 -name "${FSV_SAMPLE}*.tsv" | wc -l)"
}

###############################################################################
# recon -- the reconciliation chain and the analysis tables
###############################################################################
stage_recon() {
    say "recon: reconciliation + analysis tables, arm=${FSV_ARM}"
    eval "$FSV_CONDA_ACTIVATE"
    "$FSV_PYTHON" "$FSV_RECON" "$FSV_OUTDIR" "$FSV_SAMPLE" "$FSV_CELLDIR" \
        "$FSV_LOGDIR" "$FSV_ARM" "${FSV_RES}/reconciliation_${FSV_ARM}.tsv" \
        "${FSV_RES}/recon_report_${FSV_ARM}.txt" || die "reconcile.py failed"
    "$FSV_PYTHON" "$FSV_BUILD_TABLES" "$FSV_OUTDIR" "$FSV_SAMPLE" "$FSV_ARM" \
        "$FSV_TABLES" "${FSV_RES}/analysis_filter_${FSV_ARM}.txt" || die "build_analysis_tables.py failed"
    say "recon done"
}

###############################################################################
# status
###############################################################################
stage_status() {
    rule; say "status -- arm=${FSV_ARM} sample=${FSV_SAMPLE}"; rule
    printf "  %-34s %s\n" "libraries configured" "$(nlibs)"
    printf "  %-34s %s\n" "LUT"    "$([ -s "$FSV_VASALEN_LUT" ] && wc -l < "$FSV_VASALEN_LUT" || echo no)"
    printf "  %-34s %s\n" "prep fastq"  "$(find "$FSV_CELLDIR" -name '*_cbc_noumi_R1.fq.gz' 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "map BAM"     "$(find "$FSV_CELLDIR" -name '*_E99_Aligned.out.bam' 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "assign single bed" "$(find "$FSV_CELLDIR" -name '*singlemappers_genes.bed.gz' 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "assign multi bed"  "$(find "$FSV_CELLDIR" -name '*multimappers_genes.bed.gz' 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "pickle1 per-library dicts" "$(find "${FSV_OUTDIR}/pickle" -name '*dict.pickle' 2>/dev/null | wc -l)"
    printf "  %-34s %s\n" "merged pickle" "$([ -s "${FSV_OUTDIR}/${FSV_SAMPLE}.pickle.gz" ] && echo yes || echo no)"
    printf "  %-34s %s\n" "step7 tables" "$(find "$FSV_OUTDIR" -maxdepth 1 -name "${FSV_SAMPLE}*.tsv" 2>/dev/null | wc -l)"
    rule
}

usage() { sed -n '2,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

mkdir -p "$FSV_LOGDIR" 2>/dev/null || true
[ $# -gt 0 ] || usage
for cmd in "$@"; do
    case "$cmd" in
        check)  stage_check ;;
        lut)    stage_lut    2>&1 | tee "${FSV_LOGDIR}/lut.log" ;;
        prep)   stage_prep   2>&1 | tee "${FSV_LOGDIR}/prep.log" ;;
        map)    stage_map    2>&1 | tee "${FSV_LOGDIR}/map.log" ;;
        assign) stage_assign 2>&1 | tee "${FSV_LOGDIR}/assign.log" ;;
        pickle1)      stage_pickle1      2>&1 | tee -a "${FSV_LOGDIR}/pickle1.log" ;;
        pickle_merge) stage_pickle_merge 2>&1 | tee "${FSV_LOGDIR}/pickle_merge.log" ;;
        tables) stage_tables 2>&1 | tee "${FSV_LOGDIR}/tables.log" ;;
        recon)  stage_recon  2>&1 | tee "${FSV_LOGDIR}/recon.log" ;;
        status) stage_status ;;
        help|-h|--help) usage ;;
        *) echo "unknown stage: $cmd"; echo; usage ;;
    esac
done
