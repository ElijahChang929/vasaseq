#!/usr/bin/env Rscript
# One number per dataset: probe-targeted rRNA the RNase H reaction left behind,
# as a share of reads entering step 3. Scored on the 47S contig only -- the one
# contig both rRNA references name identically -- so the two sides are exactly
# comparable. See probe_qc.sh and probe_reference/README.md.
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

d <- read_tsv(file.path(TAB, "cross", "probe_qc.tsv"), show_col_types = FALSE) %>%
  mutate(dataset = factor(dataset, levels = dataset))
p <- ggplot(d, aes(dataset, pct, fill = dataset)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.2f%%", pct)), vjust = -0.4, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Probe-targeted rRNA remaining (% of reads)", title = NULL) +
  theme_classic(base_size = 14) + theme(legend.position = "none")
ggsave(file.path(FIG, FIGDIR[["rrna"]], "probe_qc.png"), p, width = 6, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(FIG, FIGDIR[["rrna"]], "probe_qc.pdf"), p, width = 6, height = 4.5)
cat("wrote figures/03_rrna/probe_qc.{png,pdf}\n")
