#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Read length by STAR outcome (step 4), three libraries.
#
# step4_mapping.png says the own library multimaps 28.6% (130 nt) / 33.0%
# (75 nt) against the published plate's 13.3%. These two figures ask why.
#
#   A  step4_length_by_fate         the length distribution of each class,
#                                   one panel per library
#   B  step4_multi_rate_by_length   of the reads STAR PLACED at a given length,
#                                   what share it placed at more than one locus
#
# B is the one that separates the two explanations. If the three curves lie on
# top of each other, the answer is "our reads are shorter". If the own library
# sits above the plate at the SAME length, length is not the cause and the next
# question is what the multimappers are.
#
# Note the length axis is not neutral between the libraries: own130 reads run to
# 130 nt, the two 75 nt libraries stop at 75, and the published plate is mapped
# against a human+mouse MIXED index -- a doubled reference, which should raise
# its multimapping, not lower it.
#
# `toomany` (STAR's uT:A:3, >20 loci) is folded into `multi` HERE, not in the
# tally: a read at 20+ loci is a multimapper, but the split is 0.4-0.9% and
# worth keeping on disk. Percentages therefore match step4_mapping.png's
# Multi-mapped + Too many loci, not its Multi-mapped alone.
#
# Input:  tables/cross/mapped_length_dist.tsv  (scripts/10_mapped_length_dist.sh)
# Output: figures/04_mapping/step4_length_by_fate.{png,pdf}
#         figures/04_mapping/step4_multi_rate_by_length.{png,pdf}
#         tables/cross/step4_length_stats.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr); library(tidyr)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

SRC <- file.path(TAB, "cross", "mapped_length_dist.tsv")
if (!file.exists(SRC)) {
  cat("  (skipping step-4 length figures -- run scripts/10_mapped_length_dist.sh first)\n")
  quit(save = "no")
}

DATASETS <- c("own library, 130 nt", "own library, 75 nt", "published, mouse wells")
FATES    <- c(unique = "Uniquely mapped", multi = "Multi-mapped", unmapped = "Unmapped")

# A length is dropped from figure B below this many placed reads in that
# library. Without it the far tail is one cell's handful of reads and the curve
# swings between 0 and 100% on counts of three.
MIN_PLACED <- 1000

d <- read_tsv(SRC, show_col_types = FALSE) %>%
  mutate(category = ifelse(category == "toomany", "multi", category)) %>%
  group_by(dataset, category, length) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%      # pool cells
  mutate(dataset  = factor(dataset, levels = DATASETS),
         category = factor(unname(FATES[category]), levels = unname(FATES)))
stopifnot(!any(is.na(d$dataset)), !any(is.na(d$category)))

# --- A. length distribution per class --------------------------------------
# Normalised WITHIN dataset x class: the classes differ 5-fold in size and the
# question is the shape of each, not how big it is (step4_mapping.png already
# answers that).
#
# CUMULATIVE, and that is not a stylistic choice. Every one of these libraries
# has a huge point mass at its full read length -- 41-53% of reads sit on the
# single 130 nt / 75 nt bar -- so on a raw per-length curve that one spike owns
# the whole y axis and the entire body of the distribution is pressed flat onto
# zero. All three classes then look identical, which is the opposite of the
# finding. A log y axis rescues the body but turns the 15-70 nt range into
# noise. The cumulative curve keeps both: the body is the rise, the spike is the
# final vertical jump, and "which class is shorter" is just which curve is left.

da <- d %>%
  group_by(dataset, category) %>% arrange(length, .by_group = TRUE) %>%
  mutate(cum = 100 * cumsum(reads) / sum(reads)) %>% ungroup()

p <- ggplot(da, aes(length, cum, colour = category)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~dataset, nrow = 1) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 130)) +
  labs(x = "Read length entering STAR (nt)",
       y = "% of that class's reads at or below", colour = NULL,
       title = "Read length by STAR outcome") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))
ggsave(file.path(FIG, FIGDIR[["mapping"]], "step4_length_by_fate.png"), p,
       width = 11, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIG, FIGDIR[["mapping"]], "step4_length_by_fate.pdf"), p, width = 11, height = 5)
cat("wrote figures/04_mapping/step4_length_by_fate.{png,pdf}\n")

# --- B. multimapping rate as a function of length --------------------------
# Denominator is PLACED reads (unique + multi), not all reads at that length.
# Including the unmapped would fold mappability into a curve that is meant to
# be about ambiguity alone.

db <- d %>%
  select(dataset, category, length, reads) %>%
  pivot_wider(names_from = category, values_from = reads, values_fill = 0) %>%
  transmute(dataset, length,
            placed     = `Uniquely mapped` + `Multi-mapped`,
            multi_rate = 100 * `Multi-mapped` / placed) %>%
  filter(placed >= MIN_PLACED)

p <- ggplot(db, aes(length, multi_rate, colour = dataset)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 130)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Read length entering STAR (nt)",
       y = "% of placed reads mapped to >1 locus", colour = NULL,
       title = "Multimapping rate against read length") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))
ggsave(file.path(FIG, FIGDIR[["mapping"]], "step4_multi_rate_by_length.png"), p,
       width = 8.5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIG, FIGDIR[["mapping"]], "step4_multi_rate_by_length.pdf"), p, width = 8.5, height = 5)
cat("wrote figures/04_mapping/step4_multi_rate_by_length.{png,pdf}\n")

# --- the numbers behind the curves -----------------------------------------
# So that nothing quoted in README.md is read off a plotted line.

wq <- function(len, w, p) {                # weighted quantile of a length tally
  o <- order(len); len <- len[o]; w <- w[o]
  len[which(cumsum(w) >= p * sum(w))[1]]
}
st <- d %>%
  group_by(dataset, category) %>%
  summarise(n = sum(reads), mean = sum(length * reads) / sum(reads),
            median = wq(length, reads, 0.5),
            q1 = wq(length, reads, 0.25), q3 = wq(length, reads, 0.75),
            .groups = "drop") %>%
  mutate(mean = round(mean, 1))
write_tsv(st, file.path(TAB, "cross", "step4_length_stats.tsv"))
cat("wrote tables/cross/step4_length_stats.tsv\n")
print(as.data.frame(st), row.names = FALSE)

cat("\nmultimapping rate at matched length (% of placed reads):\n")
db %>% filter(length %in% c(25, 40, 55, 70, 75)) %>%
  select(dataset, length, multi_rate) %>%
  mutate(multi_rate = round(multi_rate, 1)) %>%
  pivot_wider(names_from = length, values_from = multi_rate) %>%
  as.data.frame() %>% print(row.names = FALSE)
