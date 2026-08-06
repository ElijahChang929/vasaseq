# Where everything lives. Sourced by every R script here, so the layout is
# defined once. ROOT is derived from this file's own location, so the scripts
# work from any working directory and moving the tree needs no edits.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- normalizePath(file.path(dirname(.self), ".."))
TAB  <- file.path(ROOT, "tables")
FIG  <- file.path(ROOT, "figures")

VASA <- "/nemo/lab/turnerj/working/guangxin/vasaseq"
# demo_analyze already resolves the VASA design xlsx into a sample sheet and the
# plate's species calls into groups. Those are read, not re-derived -- a second
# implementation of the barcode-to-well join is exactly the kind of thing that
# drifts silently, and demo_analyze's version hard-errors on a failed join.
DEMO <- file.path(VASA, "code/I_Gene_expression/own_version/demo_analyze/tables")
FSMETA <- file.path(VASA, "data/flashseq_vasa/FS_native_library_metadata.tsv")

# The four datasets, in the order every figure shows them. These strings must
# match scripts/datasets.sh's ds_label exactly -- the scanners write them into
# the TSVs and the plots join on them.
DATASETS <- c("VASA own, 130 nt", "VASA own, 75 nt",
              "VASA published, mouse", "FLASH-seq")

# Numbers are kept ALIGNED WITH demo_analyze's, so 03_rrna here means the same
# thing as 03_rrna there. 02 is deliberately absent: demo_analyze's 02_trim holds
# the trim-funnel and read-class figures, and those have no four-way counterpart
# (FLASH-seq trims in one cutadapt pass, VASA in three -- see README). A gap in
# the sequence is the honest form; an empty 02_length directory just looks like
# a figure group that failed to render.
FIGDIR <- c(reads = "01_reads", rrna = "03_rrna",
            mapping = "04_mapping", assign = "05_assign",
            coverage = "06_coverage")   # gene body coverage; no demo_analyze twin
for (d in FIGDIR) dir.create(file.path(FIG, d), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(TAB, "cross"), showWarnings = FALSE, recursive = TRUE)

theme_demo <- function(base = 14) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}
save_fig <- function(p, name, width, height, sub) {
  d <- file.path(FIG, FIGDIR[[sub]])
  ggplot2::ggsave(file.path(d, paste0(name, ".png")), p, width = width, height = height,
                  dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(d, paste0(name, ".pdf")), p, width = width, height = height)
  cat("wrote figures/", FIGDIR[[sub]], "/", name, ".{png,pdf}\n", sep = "")
}
