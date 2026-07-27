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
# THE CURL SHIM -- WITHOUT THIS, scplotter CANNOT BE INSTALLED
# The r4.3 env ships its own curl (8.18.0, conda build, OpenSSL 3.6.1) at
# $FS_R_ENV/bin/curl, and conda activation puts it ahead of /usr/bin/curl. That
# build CANNOT REACH Bioconductor's archive hosts: bioconductor.org 302-redirects
# /packages/3.18/... to an Open Storage Network bucket (mghp.osn.xsede.org), and
# the conda curl times out on it after 300 s while /usr/bin/curl fetches it in 5.
# R inherits the problem twice over -- its linked libcurl is the same build, and
# download.file(method="curl") finds the conda binary on PATH.
#
# Measured 2026-07-27: available.packages() on Bioconductor 3.18 returns 0
# packages as shipped, and 2216 with a /usr/bin/curl shim ahead on PATH. That is
# the whole difference between scplotter installing and not, because scplotter
# Imports scRepertoire, which is Bioconductor-only. The earlier "dependency 'gsl'
# is not available" and "download of package 'quantreg' failed" were downstream
# of the same broken network, not separate problems: the GSL C library and
# gsl-config 2.7 are already inside the env.
#
# So: a one-symlink shim directory is created and prepended to PATH, and R is
# told to download via that curl rather than its own libcurl. Nothing in the
# shared env is touched.
#
#   bash code/flashseq/install_r_deps.sh
#
# Idempotent: already-installed packages are skipped. Ends by loading every
# package and printing its version, so a partially failed install cannot look
# like a successful one -- install.packages() only WARNS on failure, which is
# exactly how the first attempt at this looked successful while installing
# nothing.
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

# The shim, created AFTER activation so it wins against the env's own curl.
# See "THE CURL SHIM" above -- this is load-bearing, not a tidy-up.
SHIM="$FS_R_LIB/.bin"
mkdir -p "$SHIM"
ln -sf /usr/bin/curl "$SHIM/curl"
export PATH="$SHIM:$PATH"
[ "$(command -v curl)" = "$SHIM/curl" ] || { echo "FATAL: shim not on PATH"; exit 1; }

R_LIB="$FS_R_LIB" NCPUS="${FS_R_NCPUS:-4}" R --no-save --quiet <<'RSCRIPT'
lib <- Sys.getenv("R_LIB")
.libPaths(c(lib, .libPaths()))

# PPM serves binaries only when R identifies its platform in the User-Agent.
BIOC <- "3.18"   # the release that goes with R 4.3; 3.17 is archived the same way
options(
    repos = c(PPM  = "https://packagemanager.posit.co/cran/__linux__/centos8/latest",
              CRAN = "https://cloud.r-project.org",
              BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", BIOC),
              BioCann  = sprintf("https://bioconductor.org/packages/%s/data/annotation", BIOC),
              BioCexp  = sprintf("https://bioconductor.org/packages/%s/data/experiment", BIOC)),
    HTTPUserAgent = sprintf(
        "R/%s R (%s)", getRversion(),
        paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])),
    # Both of these matter: the method sends downloads through the shimmed
    # /usr/bin/curl instead of the env's libcurl, and the archive redirect is
    # slow enough that the 60 s default trips on the larger tarballs.
    download.file.method = "curl",
    download.file.extra = "-L -f -s",
    timeout = 600,
    Ncpus = as.integer(Sys.getenv("NCPUS", "4"))
)
cat("library:", lib, "\n")
cat("R:", R.version.string, "\n")
cat("curl:", Sys.which("curl"), "\n")

# Prove the shim works before spending twenty minutes finding out it does not.
nbioc <- tryCatch(nrow(available.packages(repos = getOption("repos")[["BioCsoft"]])),
                  error = function(e) 0)
cat("Bioconductor", BIOC, "packages visible:", nbioc, "\n")
if (nbioc < 100) {
    stop("Bioconductor is unreachable (", nbioc, " packages). The curl shim is ",
         "the fix for this -- see the script header. scplotter Imports ",
         "scRepertoire, which is Bioconductor-only, so this is fatal.")
}
cat("\n")

need <- function(p) !requireNamespace(p, quietly = TRUE)

# gsl must be pinned and must come first. Current CRAN gsl (2.1-9) declares
# R (>= 4.5.0) and this env is R 4.3.3, so install.packages() refuses it -- and
# the refusal surfaces four levels away, as "lazy loading failed for powerTCR",
# because the chain is scplotter -> scRepertoire -> powerTCR -> gsl. 2.1-8 is
# the last version without that floor. It builds from source against the GSL C
# library already inside the env (gsl-config 2.7), using the system gcc set
# above.
if (need("gsl")) {
    cat("installing gsl 2.1-8 from the CRAN archive (2.1-9 needs R >= 4.5)\n")
    install.packages(
        "https://cran.r-project.org/src/contrib/Archive/gsl/gsl_2.1-8.tar.gz",
        repos = NULL, type = "source", lib = lib)
}

# scplotter's Imports that r4.3 does not already carry, from CRAN/PPM.
from_repos <- c("remotes", "plotthis", "circlize", "ggnewscale")
todo <- Filter(need, from_repos)
if (length(todo)) {
    cat("installing:", paste(todo, collapse = ", "), "\n")
    install.packages(todo, lib = lib)
} else {
    cat("repository packages already present\n")
}

# scRepertoire MUST BE >= 2, and Bioconductor 3.18 -- the release paired with
# R 4.3 -- ships 1.12.0. Installing that one succeeds and then scplotter fails
# at byte-compile with "object 'clonalAbundance' is not exported by
# 'namespace:scRepertoire'", because the v1 -> v2 rewrite renamed the API. So
# the version is what is checked, not merely the presence, and 2.x comes from
# GitHub. Its own dependency powerTCR is what needs the pinned gsl above.
if (!requireNamespace("scRepertoire", quietly = TRUE) ||
    packageVersion("scRepertoire") < "2.0.0") {
    have <- if (requireNamespace("scRepertoire", quietly = TRUE))
                as.character(packageVersion("scRepertoire")) else "none"
    cat("scRepertoire", have, "-> installing v2.2.1 from GitHub (scplotter needs >= 2)\n")
    if (need("remotes")) install.packages("remotes", lib = lib)
    remotes::install_github("BorchLab/scRepertoire@v2.2.1", lib = lib,
                            upgrade = "never")
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

cat("\ninstalls done; verifying in a FRESH session\n")
RSCRIPT

###############################################################################
# Verification, deliberately in a separate R process.
#
# It cannot share the installing session: that one calls requireNamespace() to
# check scRepertoire's version, which LOADS scRepertoire 1.12.0's namespace and
# keeps it for the rest of the session. scplotter then byte-compiles against the
# stale namespace and reports "object 'clonalAbundance' is not exported" even
# though 2.2.1 is on disk by then. A new process reads what was actually
# installed.
#
# Loading is the test at all because install.packages() only WARNS when an
# install fails -- the first run of this script "succeeded" while installing
# nothing.
###############################################################################
R_LIB="$FS_R_LIB" R --no-save --quiet <<'RSCRIPT'
.libPaths(c(Sys.getenv("R_LIB"), .libPaths()))
cat("\n--- verifying by loading ---\n")
ok <- TRUE
for (p in c("plotthis", "scplotter", "scRepertoire", "ggplot2", "dplyr",
            "tidyr", "rmarkdown", "knitr")) {
    v <- tryCatch({
        suppressPackageStartupMessages(library(p, character.only = TRUE))
        as.character(packageVersion(p))
    }, error = function(e) {ok <<- FALSE; paste("FAILED:", conditionMessage(e))})
    cat(sprintf("  %-13s %s\n", p, v))
}
if (!ok) quit(status = 1)
cat("\nall packages load.\n")
RSCRIPT

echo "OK -- FS_R_LIB=$FS_R_LIB"
