#!/usr/bin/env Rscript
# 384-well plate layout for the VASA-seq experimental design.
#
# Design: one barcode per ROW. Every well within a row carries the same
# barcode, so the 16 rows (A-P) of a 384 plate map one-to-one onto the 16
# cell barcodes of bc_PM26037_6nt.tsv, and the 24 columns are replicates of
# that barcode.
#
# Run with the dedicated env:
#   /nemo/lab/turnerj/working/guangxin/envs/ggplate/bin/Rscript plate_design.R

suppressPackageStartupMessages({
  library(ggplate)
  library(ggplot2)
})

outdir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))[1]))

# --- barcodes -------------------------------------------------------------
# Sequences read off bc_PM26037_6nt.tsv (16 x 6 nt). Kept here as data so the
# figure script does not depend on a path outside this directory; the order is
# the file's order, i.e. barcode 01..16.
barcodes <- c(
  "ACTCGA", "AGACTC", "AGCTAG", "AGCTCA",
  "AGCTTC", "CAGATC", "CATGAG", "CATGCA",
  "CATGTC", "GTCTAG", "GTGACA", "GTTGCA",
  "TCACAG", "TGCAGA", "TGTCAC", "TGTCGA"
)

rows <- LETTERS[1:16]   # A-P
cols <- 1:24

plate <- expand.grid(row = rows, col = cols, stringsAsFactors = FALSE)
plate$position <- paste0(plate$row, plate$col)
plate$barcode_n <- match(plate$row, rows)
plate$barcode_seq <- barcodes[plate$barcode_n]
# Character, not factor: plate_plot() calls min()/max() on the value column,
# which errors on a factor. Zero-padding keeps the alphabetical order that
# ggplate uses identical to the numeric order 1..16.
plate$barcode <- sprintf("barcode %02d", plate$barcode_n)

# Only column 1 is used: 16 wells, one per barcode. The rest of the plate is
# left empty (NA), so the figure keeps the 384 frame while showing what was
# actually loaded.
plate$barcode[plate$col != 1] <- NA

# --- samples --------------------------------------------------------------
# The 16 loaded wells are four samples of four barcodes each, top to bottom.
sample_groups <- data.frame(
  sample    = c("XY", "XO", "EpiLCs", "NTC"),
  first_bc  = c(1, 5,  9, 13),
  last_bc   = c(4, 8, 12, 16),
  stringsAsFactors = FALSE
)
# Ordered factor levels = plate order, so the legend/table read A->P.
sample_of_bc <- rep(sample_groups$sample, times = with(sample_groups, last_bc - first_bc + 1))
plate$sample <- sample_of_bc[plate$barcode_n]
plate$sample[is.na(plate$barcode)] <- NA

plate <- plate[order(plate$barcode_n, plate$col), ]

write.table(plate[, c("position", "row", "col", "sample", "barcode", "barcode_seq")],
            file.path(outdir, "plate_layout_384.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- plot -----------------------------------------------------------------
# 16 discrete levels: ggplot2's own default discrete palette, expanded to 16.
pal <- scales::hue_pal()(length(rows))

p <- plate_plot(
  data        = plate,
  position    = position,
  value       = barcode,
  plate_size  = 384,
  plate_type  = "round",
  colour      = pal,
  title       = "",
  show_legend = FALSE,
  scale       = 1
)

# Sample brackets to the right of column 1. plate_plot()'s panel is a plain
# continuous grid: x = column number (1..24), y = row counted from the bottom,
# so row A (barcode 1) sits at y = 16 and row P at y = 1.
brk <- transform(sample_groups,
                 ytop    = length(rows) + 1 - first_bc,
                 ybottom = length(rows) + 1 - last_bc)
brk$ymid <- (brk$ytop + brk$ybottom) / 2

p <- p +
  geom_segment(data = brk,
               aes(x = 1.7, xend = 1.7, y = ytop + 0.35, yend = ybottom - 0.35),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = brk,
            aes(x = 1.95, y = ymid, label = sample),
            inherit.aes = FALSE, hjust = 0, size = 4)

ggsave(file.path(outdir, "plate_layout_384.pdf"), p, width = 8, height = 5.5)
ggsave(file.path(outdir, "plate_layout_384.png"), p, width = 8, height = 5.5, dpi = 300)

cat("wrote:\n  plate_layout_384.tsv\n  plate_layout_384.pdf\n  plate_layout_384.png\n")
