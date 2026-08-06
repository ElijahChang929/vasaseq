#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# What step 5 assigned the reads to -- single-mappers vs multi-mappers.
#
#   A  step5_gene_assignment    of the reads STAR placed, the share step 5 put
#                               on a gene, by class and library
#   B  step5_biotype_by_class   the biotypes those reads touch, own vs published
#   C  step5_genes_per_read     how many distinct genes one multimapper reaches
#   D  step5_biotype_by_length  B again, but at matched read length
#
# B is the one the previous step asked for. step4_multi_rate_by_length.png
# showed the own library multimaps ~2.4x the plate AT THE SAME LENGTH, and the
# standing hypothesis in README.md was residual rDNA -- a repeat array
# multimaps by construction. B tests it: if it is right, the own library's
# multimappers must sit on rRNA (and the rDNA-adjacent classes) far more than
# the plate's.
#
# B IS NOT A COMPOSITION, AND IS NOT STACKED. A multimapper's loci usually do
# not agree on a biotype -- most multimapper reads touch more than one -- so a
# read is counted once per DISTINCT biotype it touches and the bars sum past
# 100%. Stacking them would assert a partition that does not exist. Dodged bars
# with the library on the fill make the own-vs-published comparison the thing
# being read, which is the question.
#
# D IS THE CONTROL ON B, and it is the same trick step 4 used. B's small-RNA gap
# has an innocent explanation available: snRNA/miRNA/snoRNA genes are short, a
# short read fits inside one and a long read does not (step 5 drops a read that
# spans past both ends of a feature, `jS:OUT`), and our library's reads are
# shorter. So "more small RNA" could just be "more short reads" restated. D holds
# length fixed and asks the question again. If the three curves lie on top of
# each other, B was measuring read length. If ours sits above the plate at every
# length, B is measuring the library.
#
# Input:  tables/cross/step5_biotype_per_cell.tsv   (scripts/12_step5_biotype.sh)
#         tables/cross/step5_genes_per_read.tsv
#         tables/cross/step5_assign_totals.tsv
#         tables/cross/step5_biotype_by_length.tsv
# Output: figures/05_assign/step5_gene_assignment.{png,pdf}
#         figures/05_assign/step5_biotype_by_class.{png,pdf}
#         figures/05_assign/step5_genes_per_read.{png,pdf}
#         figures/05_assign/step5_biotype_by_length.{png,pdf}
#         tables/cross/step5_biotype_summary.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr); library(tidyr)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

TOT <- file.path(TAB, "cross", "step5_assign_totals.tsv")
BIO <- file.path(TAB, "cross", "step5_biotype_per_cell.tsv")
NGN <- file.path(TAB, "cross", "step5_genes_per_read.tsv")
BLN <- file.path(TAB, "cross", "step5_biotype_by_length.tsv")
if (!all(file.exists(TOT, BIO, NGN))) {
  cat("  (skipping step-5 figures -- run scripts/12_step5_biotype.sh first)\n")
  quit(save = "no")
}

DATASETS <- c("own library, 130 nt", "own library, 75 nt", "published, mouse wells")
CLASSES  <- c(single = "Uniquely mapped", multi = "Multi-mapped")
NTOP     <- 12          # biotypes shown in B, ranked by their multimapper share

fac <- function(d) mutate(d, dataset = factor(dataset, levels = DATASETS),
                             class   = factor(unname(CLASSES[class]), levels = unname(CLASSES)))
save_fig <- function(p, name, width, height) {
  d <- file.path(FIG, FIGDIR[["assign"]])
  ggsave(file.path(d, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(d, paste0(name, ".pdf")), p, width = width, height = height)
  cat("wrote figures/05_assign/", name, ".{png,pdf}\n", sep = "")
}
theme_demo <- function() theme_classic(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))

tot <- fac(read_tsv(TOT, show_col_types = FALSE))
stopifnot(!any(is.na(tot$dataset)), !any(is.na(tot$class)))

# --- A. did the read reach a gene ------------------------------------------
# Denominator is that class's own STAR total, per cell, summed over cells --
# so this is reads-weighted, not a mean of per-cell rates.

a <- tot %>%
  group_by(dataset, class) %>%
  summarise(star = sum(star_reads), bed = sum(bed_reads), .groups = "drop") %>%
  mutate(pct = 100 * bed / star)

p <- ggplot(a, aes(class, pct, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3.9) +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% of that class's STAR reads assigned to a gene",
       fill = NULL, title = "Reads reaching a gene (step 5)") +
  theme_demo() + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))
save_fig(p, "step5_gene_assignment", 8.5, 5.2)

# --- B. what they are ------------------------------------------------------

bio <- fac(read_tsv(BIO, show_col_types = FALSE)) %>%
  group_by(dataset, class, biotype) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  left_join(select(a, dataset, class, bed), by = c("dataset", "class")) %>%
  mutate(pct = 100 * reads / bed)

write_tsv(arrange(bio, dataset, class, desc(pct)),
          file.path(TAB, "cross", "step5_biotype_summary.tsv"))
cat("wrote tables/cross/step5_biotype_summary.tsv\n")

# Ranked on the multimappers, because that is the class the question is about;
# the same biotypes are then shown for the singlemappers so the two panels are
# directly comparable rather than each having its own top 12.
top <- bio %>% filter(class == unname(CLASSES["multi"])) %>%
  group_by(biotype) %>% summarise(m = max(pct), .groups = "drop") %>%
  slice_max(m, n = NTOP) %>% pull(biotype)

pb <- bio %>% filter(biotype %in% top) %>%
  mutate(biotype = factor(biotype, levels = rev(top)))

p <- ggplot(pb, aes(biotype, pct, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           colour = "white", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~class, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "% of that class's assigned reads touching the biotype",
       fill = NULL, title = "Biotypes the assigned reads touch (step 5)") +
  theme_demo() + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))
save_fig(p, "step5_biotype_by_class", 11, 6.5)

# --- C. how ambiguous ------------------------------------------------------
# Distinct GENES, not loci: two alignments inside one gene are one gene here,
# which is the thing that makes a read hard to count. Capped at 20 by the
# tally, so the last bar is "20 or more".

ng <- fac(read_tsv(NGN, show_col_types = FALSE)) %>%
  filter(class == unname(CLASSES["multi"])) %>%
  group_by(dataset, ngenes) %>% summarise(reads = sum(reads), .groups = "drop") %>%
  group_by(dataset) %>% mutate(pct = 100 * reads / sum(reads)) %>% ungroup()

# x is made discrete before plotting. On a continuous axis position_dodge()
# spreads the three bars of a group over most of an integer step, so a group
# straddles two ticks and the eye reads every bar off the wrong count. The
# distribution really is this jagged -- 3,5,7,9 above their even neighbours,
# and a spike at 12 -- which is exactly why the x reading has to be unambiguous.
lv <- as.character(1:max(ng$ngenes))
lv[length(lv)] <- paste0(max(ng$ngenes), "+")               # the tally's cap
ngp <- mutate(ng, ngenes = factor(ngenes, levels = seq_along(lv), labels = lv))

p <- ggplot(ngp, aes(ngenes, pct, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Distinct genes the read was assigned to",
       y = "% of multimapper reads", fill = NULL,
       title = "How many genes one multimapper reaches") +
  theme_demo() + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))
save_fig(p, "step5_genes_per_read", 8.5, 5.2)

# --- D. the same, at matched read length -----------------------------------
# Multimappers only: that is the population step 4 left a question about.
#
# Denominator is reads OF THAT LENGTH in that library's multimapper BED, so
# every point is a within-library, within-length share -- the read-length
# distributions differ ~2x between the libraries and must not leak in.
#
# Free y scales, deliberately: snRNA runs to ~40% and snoRNA to ~6%, and on a
# shared axis the three small panels flatten to a line at zero. The comparison
# being made is between the three curves inside one panel, never across panels.

MIN_AT_LENGTH <- 1000   # lengths thinner than this are one cell's noise

SMALL <- c("snRNA", "miRNA", "snoRNA", "MiscRna", "ProcessedPseudogene")
if (file.exists(BLN)) {
  bl <- fac(read_tsv(BLN, show_col_types = FALSE)) %>%
    filter(class == unname(CLASSES["multi"]), biotype %in% SMALL,
           reads_at_length >= MIN_AT_LENGTH) %>%
    mutate(pct = 100 * reads / reads_at_length,
           biotype = factor(biotype, levels = SMALL))

  p <- ggplot(bl, aes(length, pct, colour = dataset)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~biotype, nrow = 2, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 130)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Read length entering STAR (nt)",
         y = "% of multimapper reads of that length touching the biotype",
         colour = NULL, title = "Biotype against read length (multimappers)") +
    theme_demo() + theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 1))
  save_fig(p, "step5_biotype_by_length", 11, 7)

  cat("\nat matched length, % of multimapper reads of that length:\n")
  bl %>% filter(length %in% c(25, 40, 55, 70, 75)) %>%
    mutate(pct = round(pct, 1)) %>%
    select(biotype, dataset, length, pct) %>%
    pivot_wider(names_from = length, values_from = pct) %>%
    arrange(biotype, dataset) %>% as.data.frame() %>% print(row.names = FALSE)

  # --- the one number: direct standardisation to the plate's read lengths ---
  # Reading five lengths off a curve is not an answer. This reweights each own
  # library's per-length rates by the PLATE's read-length distribution -- "what
  # our number would be if our reads were as long as theirs" -- over the lengths
  # both libraries actually have (>= MIN_AT_LENGTH in each). Anything the
  # standardised number keeps is not read length.
  #
  # The plate stops at 75 nt, so this window is 15-75 nt and own130's 40% of
  # reads above 75 nt are outside it. That is the price of a matched comparison,
  # not a defect: above 75 nt there is nothing to match against.
  ref <- bl %>% filter(dataset == DATASETS[3]) %>%
    distinct(length, w = reads_at_length)
  std <- bl %>% inner_join(ref, by = "length") %>%
    group_by(biotype, dataset) %>%
    summarise(crude = 100 * sum(reads) / sum(reads_at_length),
              standardised = sum(w * pct) / sum(w),
              lengths = n(), .groups = "drop") %>%
    mutate(across(c(crude, standardised), ~round(.x, 1)))

  # The plate standardised to its own length distribution must equal its crude
  # rate. It does, to the decimal -- that is the formula checking itself.
  cat("\nstandardised to the plate's read-length distribution (15-75 nt window):\n")
  std %>%
    left_join(std %>% filter(dataset == DATASETS[3]) %>% select(biotype, ref = standardised),
              by = "biotype") %>%
    mutate(vs_plate = round(standardised / ref, 1)) %>%
    select(biotype, dataset, crude, standardised, vs_plate) %>%
    arrange(biotype, dataset) %>% as.data.frame() %>% print(row.names = FALSE)
  write_tsv(std, file.path(TAB, "cross", "step5_biotype_length_standardised.tsv"))
  cat("wrote tables/cross/step5_biotype_length_standardised.tsv\n")
} else {
  cat("  (skipping figure D -- rerun scripts/12_step5_biotype.sh for the length table)\n")
}

# --- the numbers behind the figures ----------------------------------------

cat("\nreads reaching a gene:\n")
a %>% mutate(pct = round(pct, 1)) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nbiotypes touched, % of that class's assigned reads:\n")
bio %>% filter(biotype %in% top) %>%
  mutate(pct = round(pct, 2), biotype = factor(biotype, levels = top)) %>%
  select(dataset, class, biotype, pct) %>%
  pivot_wider(names_from = dataset, values_from = pct) %>%
  arrange(class, biotype) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\nmultimapper reads on >1 gene, and median genes per read:\n")
ng %>% arrange(dataset, ngenes) %>% group_by(dataset) %>%   # cumsum below needs the order
  summarise(pct_1gene = round(100 * sum(reads[ngenes == 1]) / sum(reads), 1),
            median_genes = ngenes[which(cumsum(reads) >= 0.5 * sum(reads))[1]],
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
