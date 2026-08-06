#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Gene-detection saturation, four datasets, log-x.
#
# Input:  tables/cross/saturation.tsv, saturation_qc.tsv   (scripts/09)
# Output: figures/07_saturation/saturation.{png,pdf}
#         figures/07_saturation/saturation_75k.{png,pdf}
#
# WHY LOG X AND WHY THE 75,000 RULE IS DRAWN
# ------------------------------------------
# Depth spans 5e3 to 1e6 across these libraries, so a linear axis would put
# every point that matters in the leftmost centimetre. The dashed vertical is
# 75,000 reads per cell -- the depth every sensitivity number in the paper is
# quoted at, and the depth at which the paper states the curve has NOT
# flattened. Drawing it is what stops a reader treating the right-hand end of
# one curve as that method's sensitivity.
#
# EACH POINT USES ONLY THE CELLS THAT REACH ITS DEPTH
# ---------------------------------------------------
# The paper's own rule ("only cells that were sequenced to at least 75,000
# reads were used"). So the cell set changes along the x axis, and the point
# labels carry n. A curve that ends early has run out of cells deep enough,
# which is information, not a gap -- the alternative, extrapolating shallow
# cells, would invent the very sensitivity the figure is measuring.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(readr); library(dplyr); library(ggplot2)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

f <- file.path(TAB, "cross", "saturation.tsv")
if (!file.exists(f)) stop("run scripts/09_saturation.py first -- no ", f)
sat <- read_tsv(f, show_col_types = FALSE) %>%
  mutate(dataset = factor(dataset, levels = DATASETS))
qc  <- read_tsv(file.path(TAB, "cross", "saturation_qc.tsv"), show_col_types = FALSE)
if (any(is.na(sat$dataset))) stop("a dataset label does not match paths.R's DATASETS")

p <- ggplot(sat, aes(depth, mean_genes, colour = dataset)) +
  geom_vline(xintercept = 75000, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_x_log10(breaks = c(5e3, 1e4, 2.5e4, 7.5e4, 2e5, 1e6),
                labels = c("5k", "10k", "25k", "75k", "200k", "1M")) +
  labs(title = "Gene detection vs sequencing depth",
       subtitle = "dashed = 75,000 reads/cell, the paper's quoted depth",
       x = "reads assigned to a gene, per cell (log)",
       y = "genes detected per cell", colour = NULL) +
  theme_demo() + theme(legend.position = "bottom")
save_fig(p, "saturation", 9, 6, "saturation")

# The single-depth comparison, which is the only fair "sensitivity" bar chart:
# same depth for every dataset, cells that cannot reach it excluded.
at <- sat %>% filter(depth == 75000)
if (nrow(at) > 0) {
  q <- ggplot(at, aes(dataset, mean_genes, fill = dataset)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_errorbar(aes(ymin = pmax(0, mean_genes - sd_genes), ymax = mean_genes + sd_genes),
                  width = 0.2) +
    geom_text(aes(label = sprintf("%.0f\n(n=%d)", mean_genes, cells)), vjust = -0.3, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Genes per cell at 75,000 assigned reads",
         subtitle = "error bars = s.d. across cells deep enough to reach it",
         x = NULL, y = "genes detected per cell") +
    theme_demo() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_fig(q, "saturation_75k", 8, 6, "saturation")
} else {
  cat("no dataset has >=3 cells at 75,000 reads -- skipped the bar figure\n")
}

cat("\ncells available at each depth:\n")
sat %>% select(dataset, depth, cells, mean_genes) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = c(cells, mean_genes)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nper-dataset depth range:\n")
as.data.frame(qc) %>% print(row.names = FALSE)
