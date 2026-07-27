#!/bin/bash
# Shared paths for the FLASH-seq QC scripts. Source this, do not run it.
#
# Every other file in code/flashseq/ reads its paths from here, so this is the
# one file to edit if anything moves. Mirrors the arrangement in
# code/I_Gene_expression/own_version/config.sh.

# --- where the nf-core/rnaseq run lives -------------------------------------
# The run itself: data/flashseq/nfcore_rnaseq_all.sh, nf-core/rnaseq 3.22.2,
# --aligner star_rsem, GRCm39 + Ensembl release-116.
export FS_ROOT="${FS_ROOT:-/nemo/lab/turnerj/working/guangxin/vasaseq}"
export FS_DATA="${FS_DATA:-$FS_ROOT/data/flashseq}"
export FS_RESULTS="${FS_RESULTS:-$FS_DATA/results}"
export FS_MQC="${FS_MQC:-$FS_RESULTS/multiqc/star_rsem/multiqc_report_data}"
export FS_CODE="${FS_CODE:-$FS_ROOT/code/flashseq}"
export FS_OUT="${FS_OUT:-$FS_ROOT/res/flashseq}"

# --- raw FASTQs and the STP's own QC ----------------------------------------
# Read-only delivery directory. fastq_screen output is the STP's, not ours, and
# is the independent check on the mouse/human split reported in README.md.
export FS_FASTQ="${FS_FASTQ:-/nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/RN26038/20260325_LH00442_0237_B23GT7GLT3/fastq}"
export FS_SCREEN="${FS_SCREEN:-/nemo/lab/turnerj/inputs/genomics-stp/guangxin.zhang/RN26038/20260325_LH00442_0237_B23GT7GLT3/fastq_screen}"

# --- references -------------------------------------------------------------
# Same root and same override variable as submit_vasaplate_map.sh, so the
# FLASH-seq and VASA sides resolve references identically.
export VASA_REFS="${VASA_REFS:-/nemo/lab/turnerj/working/guangxin/reference/vasaseq}"
export FS_REF_GENOMES="${FS_REF_GENOMES:-/nemo/lab/turnerj/working/guangxin/reference/genomes}"

# The rRNA reference that actually contains the 47S unit. The Ensembl-only v1
# is deliberately NOT used -- see README.md "rRNA is not 0.8%".
export FS_RRNA_FA="${FS_RRNA_FA:-$VASA_REFS/mouse_GRCm39_E116/unique_rRNA_mouse.v2.fa}"
export FS_MOUSE_FA="${FS_MOUSE_FA:-$FS_REF_GENOMES/mus_musculus/GRCm39/genome.fa}"
export FS_MOUSE_GTF="${FS_MOUSE_GTF:-$FS_REF_GENOMES/mus_musculus/GRCm39/annotation/release-116/gtf/Mus_musculus.GRCm39.116.gtf}"
# Human genome/GTF come from the mixed-reference build tree -- they were already
# downloaded for the VASA-plate barnyard reference, so nothing new is fetched.
export FS_HUMAN_FA="${FS_HUMAN_FA:-$VASA_REFS/mixed/build/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"
export FS_HUMAN_GTF="${FS_HUMAN_GTF:-$VASA_REFS/mixed/build/Homo_sapiens.GRCh38.99.gtf.gz}"
export FS_ERCC_FA="${FS_ERCC_FA:-/nemo/lab/turnerj/working/guangxin/reference/ercc/ERCC92.fa}"
export FS_RRNA_INTERVALS="${FS_RRNA_INTERVALS:-$(dirname "$FS_RRNA_FA")/rrna_intervals.tsv}"

# --- the VASA side, borrowed by 05_rrna_bwa.sh ------------------------------
# 05 does NOT reimplement rRNA detection. It calls the same two scripts the VASA
# run called, with the same reference, so the two percentages are the same
# measurement rather than two methods that happen to be about rRNA. Anything
# reimplemented here would silently become a second method.
export FS_VASA_SCRIPTS="${FS_VASA_SCRIPTS:-$FS_ROOT/code/I_Gene_expression/a_Mapping}"
export FS_VASA_OWN="${FS_VASA_OWN:-$FS_ROOT/code/I_Gene_expression/own_version}"

# --- cluster tools, for the one stage that needs binaries -------------------
# Same module builds and same absolute-path convention as own_version/config.sh.
# Trim_Galore is deliberately NOT loaded: it is a foss-2018b build and dragging
# its libstdc++ in alongside BWA/SAMtools breaks them (see repo CLAUDE.md). 05
# therefore reproduces TrimGalore's cutadapt call directly -- see its header.
export FS_EBROOT="${FS_EBROOT:-/camp/apps/eb/software}"
export FS_P2BWA="${FS_P2BWA:-$FS_EBROOT/BWA/0.7.17-GCC-10.3.0/bin}"
export FS_P2SAMTOOLS="${FS_P2SAMTOOLS:-$FS_EBROOT/SAMtools/1.11-GCC-10.2.0/bin}"

# riboread-selection.py needs pysam, which the Anaconda base below does NOT
# have. That one script therefore runs under the VASA pipeline's env, exactly as
# it does on the VASA side. Nothing is installed into it here -- 05 only reads.
export FS_VASA_ENV="${FS_VASA_ENV:-/nemo/lab/turnerj/working/guangxin/envs/vasa}"
export FS_CUTADAPT="${FS_CUTADAPT:-$FS_VASA_ENV/bin/cutadapt}"
export FS_ML_RIBO="${FS_ML_RIBO:-source /usr/share/lmod/lmod/init/bash; export MODULEPATH=$FS_EBROOT/../modules/all; module load BWA/0.7.17-GCC-10.3.0 SAMtools/1.11-GCC-10.2.0; source $FS_EBROOT/Anaconda3/2024.10-1/etc/profile.d/conda.sh; conda activate $FS_VASA_ENV}"

# --- python -----------------------------------------------------------------
# The cluster's Anaconda3 base interpreter, used directly. Nothing is built or
# installed: it already carries numpy 1.26, pandas 2.2, matplotlib 3.9,
# scipy 1.13, and the full nbformat/nbclient/nbconvert/ipykernel set, which is
# everything these scripts and the notebook need. Verified 2026-07-27.
#
# Three environments were tried before settling here, and the reasons they were
# rejected are worth keeping so nobody repeats them:
#   - envs/vasa (the VASA pipeline's env) has pandas but no matplotlib. Adding
#     to it was rejected: it is in live use by the VASA mapping run, and a conda
#     solve could move numpy/pandas underneath it.
#   - A venv on the matplotlib/3.7.2-gfbf-2023a module works only with a
#     PYTHONPATH prepend hack: EasyBuild modules put their packages on
#     PYTHONPATH, which outranks venv site-packages, and that module's
#     Python-bundle-PyPI 2023.06 ships a typing_extensions too old for
#     jupyter_client (TypedDict 'extra_items' -> TypeError on import).
#   - A fresh conda env: the shared Anaconda's classic solver ran >20 min
#     without producing an environment.
#
# PYTHONPATH must be CLEARED, not extended -- if any EasyBuild module is loaded
# in the calling shell its packages shadow Anaconda's and the older
# typing_extensions comes back. fs_python does that for you.
export FS_CONDA="${FS_CONDA:-/camp/apps/eb/software/Anaconda3/2024.10-1}"
export FS_PY="${FS_PY:-$FS_CONDA/bin/python}"

# Run python with a clean module environment. Use this rather than "$FS_PY".
fs_python() {
    env -u PYTHONPATH "$FS_PY" "$@"
}

# Verify the interpreter has what the scripts need. Returns non-zero if not.
fs_check_python() {
    fs_python - <<'PYEOF'
import sys
need = ["numpy", "pandas", "matplotlib", "scipy",
        "nbformat", "nbclient", "nbconvert", "ipykernel"]
missing = []
for m in need:
    try:
        __import__(m)
    except ImportError:
        missing.append(m)
if missing:
    sys.exit("FATAL: $FS_PY is missing " + ", ".join(missing))
print(f"python {sys.version.split()[0]} at {sys.executable}: all modules present")
PYEOF
}

# --- R, for flashseq_qc.Rmd -------------------------------------------------
# The R report reads the SAME TSVs in res/flashseq/ that the Python notebook
# does. It recomputes nothing, so the two reports cannot disagree on a number.
#
# Interpreter: the shared conda env envs/r4.3 (R 4.3.3, 305 packages incl.
# Seurat 5.3.0 and pandoc 3.8.3), used READ-ONLY. New packages go to FS_R_LIB
# instead, for the same reason envs/vasa is left alone above -- and with a
# sharper warning next door: envs/sct_R (R 4.5.2) is already broken, 280 package
# directories of which 14 load, because conda clobbered their DESCRIPTION files
# into DESCRIPTION.c~. Do not install into a shared env.
export FS_R_ENV="${FS_R_ENV:-/nemo/lab/turnerj/working/guangxin/envs/r4.3}"
export FS_R_LIB="${FS_R_LIB:-/nemo/lab/turnerj/working/guangxin/envs/Rlib_flashseq_4.3}"
export FS_R_NCPUS="${FS_R_NCPUS:-4}"

# Run R with FS_R_LIB ahead of the env's own library. Use this, not a bare R.
fs_R() {
    ( . "$FS_CONDA/etc/profile.d/conda.sh"
      conda activate "$FS_R_ENV"
      R_LIBS_USER="$FS_R_LIB" R "$@" )
}
fs_Rscript() {
    ( . "$FS_CONDA/etc/profile.d/conda.sh"
      conda activate "$FS_R_ENV"
      R_LIBS_USER="$FS_R_LIB" Rscript "$@" )
}

# --- sampling ---------------------------------------------------------------
# Reads per library for the two FASTQ-streaming screens. 400k R1 reads takes a
# few minutes per library and puts the sampling error on a 4% rate at ~0.03%,
# which is far below the effects being measured (rRNA 3.3-6.3%, CALB1 0-15%).
export FS_SUBSAMPLE="${FS_SUBSAMPLE:-400000}"

# TAKING THE FIRST N READS IS BIASED -- USE FS_STRIDE INSTEAD.
# Measured on ZHA8833A1 (2026-07-27) with the identical cutadapt call, reading
# the adapter-containing rate off the head of the file:
#
#     first 20k   58.5%      first 400k  57.7%      first 2M  56.4%
#     every 64th  55.2%      whole library (saved TrimGalore report)  55.1%
#
# The head of a fastq is one end of the flowcell, and adapter content -- hence
# insert length, hence anything that depends on it -- drifts along it. Head
# sampling overstated the rate by 2.6 points at 400k; a uniform every-Nth-read
# sample reproduces the whole-library figure, and reproduces its quality-trimmed
# (0.7%) and bases-written (89.2%) rates too. So sample with a stride.
#
# 64 gives ~456k reads from a ~29M-read library: binomial s.e. 0.03% on a 5%
# rate, against a 4x effect. Denser costs almost nothing in bwa (5 s per 456k
# reads) but the decompression pass -- 57 s per file, the real cost -- is paid
# whatever the stride.
export FS_STRIDE="${FS_STRIDE:-64}"
export FS_RRNA_BWA_DIR="${FS_RRNA_BWA_DIR:-$FS_OUT/rrna_bwa}"

mkdir -p "$FS_OUT"
