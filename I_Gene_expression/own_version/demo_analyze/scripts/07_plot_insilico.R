#!/usr/bin/env Rscript
# What the in-silico rRNA depletion (step 3) removed, as a share of the reads
# entering it. Read straight from each run's own step3_report / depletion table
# -- nothing is recomputed.
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

d <- read_tsv(file.path(TAB, "cross", "insilico_depletion.tsv"), show_col_types = FALSE) %>%
  mutate(dataset = factor(dataset, levels = dataset))
p <- ggplot(d, aes(dataset, pct, fill = dataset)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", pct)), vjust = -0.4, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Removed by in-silico depletion (% of reads in)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(FIG, FIGDIR[["rrna"]], "insilico_depletion.png"), p, width = 7.5, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(FIG, FIGDIR[["rrna"]], "insilico_depletion.pdf"), p, width = 7.5, height = 4.8)
cat("wrote figures/03_rrna/insilico_depletion.{png,pdf}\n")
