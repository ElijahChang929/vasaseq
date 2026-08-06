#!/usr/bin/env Rscript
# What STAR (step 4) does with the reads in-silico depletion passed through.
# Read straight from each cell's Log.final.txt -- nothing recomputed, nothing re-run.
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr); library(tidyr); library(scales)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

FATES <- c(unique = "Uniquely mapped", multi = "Multi-mapped",
           toomany = "Too many loci", unmapped = "Unmapped")
d <- read_tsv(file.path(TAB, "cross", "step4_mapping.tsv"), show_col_types = FALSE) %>%
  select(dataset, all_of(names(FATES))) %>%
  pivot_longer(-dataset, names_to = "fate", values_to = "reads") %>%
  group_by(dataset) %>% mutate(pct = 100 * reads / sum(reads)) %>% ungroup() %>%
  mutate(dataset = factor(dataset, levels = unique(dataset)),
         fate = factor(unname(FATES[fate]), levels = unname(FATES)))
p <- ggplot(d, aes(dataset, reads, fill = fate)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.6,
           colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct >= 4, sprintf("%.1f%%", pct), "")),
            position = position_fill(vjust = 0.5, reverse = TRUE), size = 4) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "Share of reads entering STAR", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 15, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1))
ggsave(file.path(FIG, FIGDIR[["mapping"]], "step4_mapping.png"), p, width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIG, FIGDIR[["mapping"]], "step4_mapping.pdf"), p, width = 8, height = 5)
cat("wrote figures/04_mapping/step4_mapping.{png,pdf}\n")
