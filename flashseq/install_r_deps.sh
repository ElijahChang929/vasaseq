#!/bin/bash
###############################################################################
# install_r_deps.sh -- put scplotter and plotthis where flashseq_qc.Rmd can
# find them, WITHOUT touching a shared environment.
#
# WHY A SEPARATE LIBRARY
# ----------------------
# The R that runs the report is the shared conda env envs/r4.3 (R 4.3.3, 305
# packages, Seurat 5.3.0, pandoc 3.8.3). Nothing is installed INTO it: it is
# shared, and envs/sct_R next door is a live demonstration of what happens when
# a conda env gets clobbered -- 280 package directories on disk, 14 of them
# loadable, because their DESCRIPTION files were renamed DESCRIPTION.c~. That
# env is unusable and this script deliberately does not try to repair it.
#
# So: r4.3 provides the interpreter and its 305 packages read-only, and the
# handful of new ones land in FS_R_LIB, which is prepended to .libPaths().
#
# WHAT IS INSTALLED AND WHY
#   plotthis    CRAN. scplotter's own plotting engine, same author. This is what
#               actually draws the report's figures -- see the Rmd's header for
#               why scplotter's own functions do not apply to per-library QC
#               tables.
#   scplotter   NOT ON CRAN -- github.com/pwwang/scplotter only, so remotes is
#               needed too. Used for the one panel where it genuinely applies.
#   the rest    scplotter's declared Imports that r4.3 lacks.
#
# BINARIES, NOT SOURCE
# `R CMD config CC` in this env reports x86_64-conda-linux-gnu-cc, a compiler
# that is not installed -- so a source build of anything with C code fails
# before it starts. Posit Package Manager serves prebuilt binaries for RHEL 8
# (this host is el8), which sidesteps compilation entirely. R only gets offered
# them if it sends a platform-specific User-Agent, hence the HTTPUserAgent
# option below. CC/CXX are also pointed at the system gcc so that any package
# PPM has no binary for can still fall back to source.
#
#   bash code/flashseq/install_r_deps.sh
#
# Idempotent: already-installed packages are skipped. Ends by loading every
# package and printing its version, so a partially failed install cannot look
# like a successful one.
###############################################################################
set -eu

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$HERE/config.sh"

mkdir -p "$FS_R_LIB"

# System gcc, for anything PPM cannot supply as a binary. GCC 8.5 (el8).
export CC=/usr/bin/gcc CXX=/usr/bin/g++ FC=/usr/bin/gfortran F77=/usr/bin/gfortran
export MAKEFLAGS="-j${FS_R_NCPUS:-4}"

# shellcheck source=/dev/null
source "$FS_CONDA/etc/profile.d/conda.sh"
conda activate "$FS_R_ENV"

R_LIB="$FS_R_LIB" NCPUS="${FS_R_NCPUS:-4}" R --no-save --quiet <<'RSCRIPT'
lib <- Sys.getenv("R_LIB")
.libPaths(c(lib, .libPaths()))

# PPM serves binaries only when R identifies its platform in the User-Agent.
options(
    repos = c(PPM = "https://packagemanager.posit.co/cran/__linux__/centos8/latest",
              CRAN = "https://cloud.r-project.org"),
    HTTPUserAgent = sprintf(
        "R/%s R (%s)", getRversion(),
        paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])),
    Ncpus = as.integer(Sys.getenv("NCPUS", "4"))
)
cat("library:", lib, "\n")
cat("R:", R.version.string, "\n\n")

need <- function(p) !requireNamespace(p, quietly = TRUE)

# CRAN first -- scplotter's Imports that r4.3 does not already carry.
from_cran <- c("remotes", "plotthis", "circlize", "ggnewscale", "scRepertoire")
todo <- Filter(need, from_cran)
if (length(todo)) {
    cat("installing from CRAN/PPM:", paste(todo, collapse = ", "), "\n")
    install.packages(todo, lib = lib)
} else {
    cat("CRAN packages already present\n")
}

# gglogger and scplotter are GitHub-only.
if (need("gglogger")) {
    cat("installing gglogger from GitHub\n")
    remotes::install_github("pwwang/gglogger", lib = lib, upgrade = "never")
}
if (need("scplotter")) {
    cat("installing scplotter from GitHub\n")
    remotes::install_github("pwwang/scplotter", lib = lib, upgrade = "never")
}

# Loading is the test. install.packages() warns rather than errors when a
# package fails to install, so "no error" proves nothing on its own.
cat("\n--- verifying by loading ---\n")
ok <- TRUE
for (p in c("plotthis", "scplotter", "ggplot2", "dplyr", "tidyr", "rmarkdown", "knitr")) {
    v <- tryCatch({
        suppressPackageStartupMessages(library(p, character.only = TRUE))
        as.character(packageVersion(p))
    }, error = function(e) {ok <<- FALSE; paste("FAILED:", conditionMessage(e))})
    cat(sprintf("  %-12s %s\n", p, v))
}
if (!ok) quit(status = 1)
cat("\nall packages load.\n")
RSCRIPT

echo "OK -- FS_R_LIB=$FS_R_LIB"
