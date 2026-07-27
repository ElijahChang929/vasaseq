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

# --- sampling ---------------------------------------------------------------
# Reads per library for the two FASTQ-streaming screens. 400k R1 reads takes a
# few minutes per library and puts the sampling error on a 4% rate at ~0.03%,
# which is far below the effects being measured (rRNA 3.3-6.3%, CALB1 0-15%).
export FS_SUBSAMPLE="${FS_SUBSAMPLE:-400000}"

mkdir -p "$FS_OUT"
