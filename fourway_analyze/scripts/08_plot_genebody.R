#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Gene body coverage, 5' -> 3', four datasets on one axis.
#
# Input:  tables/cross/genebody_coverage.tsv  (scripts/07)
#         tables/cross/genebody_qc.tsv
# Output: figures/06_coverage/genebody_coverage.{png,pdf}
#         figures/06_coverage/genebody_intronic.{png,pdf}
#
# THE Y AXIS IS A SHARE, NOT A DEPTH
# ----------------------------------
# Each gene's coverage is normalised to sum to 1 across the 100 bins before the
# genes are averaged, so the curve is "what fraction of this gene's reads land
# here" and a flat line at 0.01 is perfectly uniform. Without that per-gene
# normalisation the curve is whatever the few most-expressed genes happen to do,
# which differs between protocols and would be read as positional bias.
#
# All four curves are built from THE SAME genes -- the intersection that passes
# the depth floor in every dataset (scripts/07). Different gene sets would put
# gene composition on the same axis as positional bias.
#
# The companion intronic-fraction figure is not decoration: gene body coverage
# is exon-only by definition, so each dataset silently drops a different share
# of its reads to get onto this axis. Showing what was dropped next to the curve
# is what stops the main figure from being misread.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(readr); library(dplyr); library(ggplot2)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

f_cov <- file.path(TAB, "cross", "genebody_coverage.tsv")
f_qc  <- file.path(TAB, "cross", "genebody_qc.tsv")
if (!file.exists(f_cov)) stop("run scripts/07_genebody_coverage.py first -- no ", f_cov)

cov <- read_tsv(f_cov, show_col_types = FALSE) %>%
  mutate(dataset = factor(dataset, levels = DATASETS))
qc  <- read_tsv(f_qc,  show_col_types = FALSE) %>%
  mutate(dataset = factor(dataset, levels = DATASETS))
if (any(is.na(cov$dataset))) stop("a dataset label does not match paths.R's DATASETS")

ngene <- unique(cov$genes)
stopifnot(length(ngene) == 1)

p <- ggplot(cov, aes(bin, coverage, colour = dataset)) +
  geom_hline(yintercept = 1 / 100, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("5'", "25%", "50%", "75%", "3'")) +
  labs(title = "Gene body coverage",
       subtitle = paste0(ngene, " shared protein-coding genes  |  dashed = uniform"),
       x = "position along mature transcript", y = "fraction of gene's reads per percentile",
       colour = NULL) +
  theme_demo() + theme(legend.position = "bottom")
save_fig(p, "genebody_coverage", 9, 6, "coverage")

q <- ggplot(qc, aes(dataset, pct_intronic, fill = dataset)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct_intronic)), vjust = -0.4, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Reads excluded from the coverage curve",
       subtitle = "intronic reads have no mature-transcript coordinate",
       x = NULL, y = "% of gene-assigned reads that are intronic") +
  theme_demo() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(q, "genebody_intronic", 8, 6, "coverage")

cat("\ncoverage summary (mean of first/last 10 percentiles):\n")
cov %>% group_by(dataset) %>%
  summarise(`5p_10pct` = mean(coverage[bin <= 10]),
            mid        = mean(coverage[bin > 45 & bin <= 55]),
            `3p_10pct` = mean(coverage[bin > 90]),
            .groups = "drop") %>%
  mutate(`3p/5p` = round(`3p_10pct` / `5p_10pct`, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nintronic fraction:\n")
qc %>% select(dataset, exonic_reads, intronic_reads, pct_intronic) %>%
  as.data.frame() %>% print(row.names = FALSE)
