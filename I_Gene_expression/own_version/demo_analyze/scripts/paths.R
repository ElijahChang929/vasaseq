# Where everything lives. Sourced by every R script here so the layout is
# defined once: scripts/ holds code, tables/ holds TSVs (one directory per
# dataset plus cross/), figures/ holds PNG+PDF grouped by pipeline step.
#
# ROOT is derived from this file's own location, so the scripts work from any
# working directory and moving the tree needs no edits.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- normalizePath(file.path(dirname(.self), ".."))
TAB  <- file.path(ROOT, "tables")
FIG  <- file.path(ROOT, "figures")
OWNV <- normalizePath(file.path(ROOT, ".."))       # own_version/
FIGDIR <- c(reads = "01_reads", trim = "02_trim", rrna = "03_rrna",
            mapping = "04_mapping", assign = "05_assign", cross = "cross")
for (d in FIGDIR) dir.create(file.path(FIG, d), showWarnings = FALSE, recursive = TRUE)
