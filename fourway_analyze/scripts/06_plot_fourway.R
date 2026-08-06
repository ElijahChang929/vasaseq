#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Every cross-dataset figure, four datasets, one script.
#
# WHY ONE SCRIPT. demo_analyze learned this the hard way: its two datasets began
# with a plot script each, and restyling one silently left the other behind. One
# theme_demo(), one save_fig(), one dataset order.
#
# WHAT IS HERE AND WHAT IS DELIBERATELY NOT
# -----------------------------------------
# Only figures that all four datasets can honestly appear in. The VASA-internal
# ones stay in demo_analyze and are not reproduced here:
#
#   read classes (5-way)      needs the barcode read-through recomputed from a
#                             cell barcode inside the read. FLASH-seq has none.
#   pass-2 adapter attribution  VASA trims in three passes with a measured
#                             read-through adapter and a 5' polyT; FLASH-seq is
#                             one cutadapt call with the Nextera mosaic end.
#                             There is no shared axis to put them on.
#   UFI per unit              FLASH-seq has no UMI, so UFICounts degenerates to
#                             a detection mask. Reads are the only unit all four
#                             protocols measure.
#
# Naming these rather than quietly dropping them is the point: a four-way figure
# with one dataset missing looks like a pipeline failure.
#
# Input:  tables/units.tsv                            (scripts/05)
#         tables/cross/insilico_depletion*.tsv        (scripts/01)
#         tables/cross/probe_qc.tsv                   (scripts/02)
#         tables/cross/mapped_length_dist.tsv, step4_mapping.tsv   (scripts/03)
#         tables/cross/step5_*.tsv                    (scripts/04)
# Output: figures/{01_reads,02_length,03_rrna,04_mapping,05_assign}/*
#         tables/cross/step5_biotype_summary.tsv
#         tables/cross/step5_biotype_length_standardised.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(dplyr); library(tidyr); library(scales)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

fac <- function(d) mutate(d, dataset = factor(dataset, levels = DATASETS))
have <- function(p) { ok <- file.exists(p); if (!ok) cat("  (skipping -- no ", basename(p), ")\n", sep = ""); ok }
SEG <- list(colour = "white", linewidth = 0.4)

units <- fac(read_tsv(file.path(TAB, "units.tsv"), show_col_types = FALSE))

# =========================================================================
# 01  depth per unit
# =========================================================================
# Reads entering step 3, i.e. post-trim. See 05_build_units.R for why that is
# the only point the four share. Free x scale per facet: 16 barcodes, 173 wells
# and 10 libraries cannot share one axis.
p <- ggplot(units, aes(reorder(unit, -reads_in), reads_in / 1e6, fill = group)) +
  geom_col() +
  facet_wrap(~dataset, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Reads entering step 3 (millions)", fill = NULL,
       title = "Depth per unit") +
  theme_demo() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
save_fig(p, "depth_per_unit", 13, 5.5, "reads")

# Per group: the same numbers, but as the medians a reader would quote.
g <- units %>% group_by(dataset, group) %>%
  summarise(units = n(), total = sum(reads_in), median = median(reads_in), .groups = "drop")
p <- ggplot(g, aes(group, median / 1e6, fill = dataset)) +
  geom_col(colour = SEG$colour, linewidth = SEG$linewidth) +
  facet_grid(~dataset, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Median reads per unit (millions)", fill = NULL,
       title = "Depth per unit, by group") +
  theme_demo() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
save_fig(p, "depth_per_group", 12, 5.5, "reads")
write_tsv(g, file.path(TAB, "cross", "depth_per_group.tsv"))

# =========================================================================
# 03  in-silico rRNA depletion
# =========================================================================
# THE STRAND FLAG IS NOT SHARED. The three VASA runs ran riboread-selection.py
# with stranded=y and FLASH-seq with n, because FLASH-seq is unstranded (49-50%
# forward) and VASA is not (76% forward). `y` costs VASA ~24% of its ribosomal
# reads and would cost FLASH-seq ~50% of its. The bars are therefore each under
# the flag that is right for that protocol, which is stated on the figure rather
# than left for the reader to discover.
DEPA <- file.path(TAB, "cross", "insilico_depletion.tsv")
if (have(DEPA)) {
  d <- fac(read_tsv(DEPA, show_col_types = FALSE))
  p <- ggplot(d, aes(dataset, pct, fill = dataset)) +
    geom_col(width = 0.62, colour = SEG$colour, linewidth = SEG$linewidth) +
    geom_text(aes(label = sprintf("%.2f%%", pct)), vjust = -0.4, size = 3.9) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "% of post-trim reads called ribosomal",
         title = "In-silico rRNA depletion (step 3)",
         caption = "VASA rows measured with stranded=y, FLASH-seq with n -- each protocol's correct flag") +
    theme_demo() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
  save_fig(p, "insilico_depletion", 8.5, 5.2, "rrna")

  pu <- fac(read_tsv(file.path(TAB, "cross", "insilico_depletion_per_unit.tsv"),
                     show_col_types = FALSE)) %>%
    left_join(select(units, dataset, unit, group), by = c("dataset", "unit"))
  p <- ggplot(pu, aes(dataset, pct, colour = dataset)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
    geom_jitter(width = 0.14, height = 0, size = 1.6, alpha = 0.7) +
    labs(x = NULL, y = "% of that unit's reads called ribosomal",
         title = "In-silico rRNA depletion, per unit") +
    theme_demo() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
  save_fig(p, "insilico_depletion_per_unit", 8.5, 5.2, "rrna")
}

PQC <- file.path(TAB, "cross", "probe_qc.tsv")
if (have(PQC)) {
  d <- fac(read_tsv(PQC, show_col_types = FALSE))
  p <- ggplot(d, aes(dataset, pct, fill = dataset)) +
    geom_col(width = 0.62, colour = SEG$colour, linewidth = SEG$linewidth) +
    geom_text(aes(label = sprintf("%.2f%%", pct)), vjust = -0.4, size = 3.9) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "% of reads entering step 3 on a probe target",
         title = "Probe-target rRNA residual, 47S only") +
    theme_demo() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
  save_fig(p, "probe_qc", 8.5, 5.2, "rrna")
}

# =========================================================================
# 04  what STAR did
# =========================================================================
FATES <- c(unique = "Uniquely mapped", multi = "Multi-mapped",
           toomany = "Too many loci", unmapped = "Unmapped")
S4 <- file.path(TAB, "cross", "step4_mapping.tsv")
if (have(S4)) {
  d <- fac(read_tsv(S4, show_col_types = FALSE)) %>%
    select(dataset, all_of(names(FATES))) %>%
    pivot_longer(-dataset, names_to = "fate", values_to = "reads") %>%
    group_by(dataset) %>% mutate(pct = 100 * reads / sum(reads)) %>% ungroup() %>%
    mutate(fate = factor(unname(FATES[fate]), levels = unname(FATES)))
  p <- ggplot(d, aes(dataset, reads, fill = fate)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.6,
             colour = SEG$colour, linewidth = SEG$linewidth) +
    geom_text(aes(label = ifelse(pct >= 4, sprintf("%.1f%%", pct), "")),
              position = position_fill(vjust = 0.5, reverse = TRUE), size = 3.9) +
    scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0.02, 0.04))) +
    labs(x = NULL, y = "Share of reads entering STAR", fill = NULL,
         title = "What STAR did with the reads (step 4)") +
    theme_demo() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 15, hjust = 1)) +
    guides(fill = guide_legend(nrow = 1))
  save_fig(p, "step4_mapping", 9, 5.5, "mapping")
}

MLD <- file.path(TAB, "cross", "mapped_length_dist.tsv")
if (have(MLD)) {
  FT <- c(unique = "Uniquely mapped", multi = "Multi-mapped", unmapped = "Unmapped")
  d <- fac(read_tsv(MLD, show_col_types = FALSE)) %>%
    mutate(category = ifelse(category == "toomany", "multi", category)) %>%
    group_by(dataset, category, length) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    mutate(category = factor(unname(FT[category]), levels = unname(FT)))

  # CUMULATIVE, not a per-length curve: every library puts 40-55% of its reads
  # on the single full-length bar, so a raw histogram gives that one spike the
  # whole y axis and every class looks identical. The cumulative form keeps the
  # body (the rise) and the spike (the final jump) in one picture.
  da <- d %>% group_by(dataset, category) %>% arrange(length, .by_group = TRUE) %>%
    mutate(cum = 100 * cumsum(reads) / sum(reads)) %>% ungroup()
  p <- ggplot(da, aes(length, cum, colour = category)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~dataset, nrow = 1) +
    labs(x = "Read length entering STAR (nt)", y = "% of that class at or below",
         colour = NULL, title = "Read length by STAR outcome") +
    theme_demo() + theme(legend.position = "bottom")
  save_fig(p, "step4_length_by_fate", 14, 5, "mapping")

  # Denominator is PLACED reads (unique + multi), never all reads at that
  # length: including the unmapped folds mappability into a curve that is about
  # ambiguity alone. Lengths under 1,000 placed reads are dropped.
  db <- d %>% select(dataset, category, length, reads) %>%
    pivot_wider(names_from = category, values_from = reads, values_fill = 0) %>%
    transmute(dataset, length,
              placed = `Uniquely mapped` + `Multi-mapped`,
              multi_rate = 100 * `Multi-mapped` / placed) %>%
    filter(placed >= 1000)
  p <- ggplot(db, aes(length, multi_rate, colour = dataset)) +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Read length entering STAR (nt)",
         y = "% of placed reads mapped to >1 locus", colour = NULL,
         title = "Multimapping rate against read length") +
    theme_demo() + theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2))
  save_fig(p, "step4_multi_rate_by_length", 9, 5.5, "mapping")

  wq <- function(len, w, p) { o <- order(len); len <- len[o]; w <- w[o]
                              len[which(cumsum(w) >= p * sum(w))[1]] }
  st <- d %>% group_by(dataset, category) %>%
    summarise(n = sum(reads), median = wq(length, reads, 0.5),
              q1 = wq(length, reads, 0.25), q3 = wq(length, reads, 0.75), .groups = "drop")
  write_tsv(st, file.path(TAB, "cross", "step4_length_stats.tsv"))
  cat("\nread length by STAR outcome:\n"); print(as.data.frame(st), row.names = FALSE)
  cat("\nmultimapping rate at matched length (% of placed reads):\n")
  db %>% filter(length %in% c(25, 40, 55, 70, 75)) %>%
    mutate(multi_rate = round(multi_rate, 1)) %>%
    pivot_wider(names_from = length, values_from = multi_rate) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

# =========================================================================
# 05  what step 5 assigned them to
# =========================================================================
CLASSES <- c(single = "Uniquely mapped", multi = "Multi-mapped")
TOT <- file.path(TAB, "cross", "step5_assign_totals.tsv")
BIO <- file.path(TAB, "cross", "step5_biotype_per_cell.tsv")
NGN <- file.path(TAB, "cross", "step5_genes_per_read.tsv")
BLN <- file.path(TAB, "cross", "step5_biotype_by_length.tsv")
facc <- function(d) fac(d) %>% mutate(class = factor(unname(CLASSES[class]),
                                                     levels = unname(CLASSES)))
if (have(TOT)) {
  tot <- facc(read_tsv(TOT, show_col_types = FALSE))
  a <- tot %>% group_by(dataset, class) %>%
    summarise(star = sum(star_reads), bed = sum(bed_reads), .groups = "drop") %>%
    mutate(pct = 100 * bed / star)
  p <- ggplot(a, aes(class, pct, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72,
             colour = SEG$colour, linewidth = SEG$linewidth) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(width = 0.78), vjust = -0.4, size = 3.4) +
    scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "% of that class's STAR reads assigned to a gene",
         fill = NULL, title = "Reads reaching a gene (step 5)") +
    theme_demo() + theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))
  save_fig(p, "step5_gene_assignment", 9, 5.5, "assign")

  # B is NOT a composition and is NOT stacked: a multimapper's loci usually
  # disagree on a biotype, so a read is counted once per DISTINCT biotype it
  # touches and the bars sum past 100%. Stacking would assert a partition that
  # does not exist.
  bio <- facc(read_tsv(BIO, show_col_types = FALSE)) %>%
    group_by(dataset, class, biotype) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    left_join(select(a, dataset, class, bed), by = c("dataset", "class")) %>%
    mutate(pct = 100 * reads / bed)
  write_tsv(arrange(bio, dataset, class, desc(pct)),
            file.path(TAB, "cross", "step5_biotype_summary.tsv"))
  top <- bio %>% filter(class == unname(CLASSES["multi"])) %>%
    group_by(biotype) %>% summarise(m = max(pct), .groups = "drop") %>%
    slice_max(m, n = 12) %>% pull(biotype)
  p <- ggplot(filter(bio, biotype %in% top) %>%
                mutate(biotype = factor(biotype, levels = rev(top))),
              aes(biotype, pct, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.74,
             colour = SEG$colour, linewidth = 0.25) +
    coord_flip() + facet_wrap(~class, nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "% of that class's assigned reads touching the biotype",
         fill = NULL, title = "Biotypes the assigned reads touch (step 5)") +
    theme_demo() + theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))
  save_fig(p, "step5_biotype_by_class", 12, 7, "assign")

  ng <- facc(read_tsv(NGN, show_col_types = FALSE)) %>%
    filter(class == unname(CLASSES["multi"])) %>%
    group_by(dataset, ngenes) %>% summarise(reads = sum(reads), .groups = "drop") %>%
    group_by(dataset) %>% mutate(pct = 100 * reads / sum(reads)) %>% ungroup()
  # x made discrete before plotting: on a continuous axis position_dodge spreads
  # a group over most of an integer step, so a group straddles two ticks and
  # every bar is read off the wrong count.
  lv <- as.character(1:max(ng$ngenes)); lv[length(lv)] <- paste0(max(ng$ngenes), "+")
  p <- ggplot(mutate(ng, ngenes = factor(ngenes, levels = seq_along(lv), labels = lv)),
              aes(ngenes, pct, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.88), width = 0.82) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Distinct genes the read was assigned to", y = "% of multimapper reads",
         fill = NULL, title = "How many genes one multimapper reaches") +
    theme_demo() + theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))
  save_fig(p, "step5_genes_per_read", 9.5, 5.5, "assign")

  cat("\nreads reaching a gene:\n")
  a %>% mutate(pct = round(pct, 1)) %>% as.data.frame() %>% print(row.names = FALSE)
  cat("\nbiotypes touched, % of that class's assigned reads:\n")
  bio %>% filter(biotype %in% top) %>%
    mutate(pct = round(pct, 2), biotype = factor(biotype, levels = top)) %>%
    select(dataset, class, biotype, pct) %>%
    pivot_wider(names_from = dataset, values_from = pct) %>%
    arrange(class, biotype) %>% as.data.frame() %>% print(row.names = FALSE)
}

# --- the length control, four ways -----------------------------------------
# The same trick step 4 used. If the biotype gaps above are just "some libraries
# have shorter reads", holding length fixed makes them close.
if (have(BLN)) {
  SMALL <- c("snRNA", "miRNA", "snoRNA", "MiscRna", "ProcessedPseudogene")
  bl <- facc(read_tsv(BLN, show_col_types = FALSE)) %>%
    filter(class == unname(CLASSES["multi"]), biotype %in% SMALL,
           reads_at_length >= 1000) %>%
    mutate(pct = 100 * reads / reads_at_length,
           biotype = factor(biotype, levels = SMALL))
  p <- ggplot(bl, aes(length, pct, colour = dataset)) +
    geom_line(linewidth = 0.85) +
    facet_wrap(~biotype, nrow = 2, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Read length entering STAR (nt)",
         y = "% of multimapper reads of that length touching the biotype",
         colour = NULL, title = "Biotype against read length (multimappers)") +
    theme_demo() + theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2))
  save_fig(p, "step5_biotype_by_length", 13, 7.5, "assign")

  # Direct standardisation onto the PUBLISHED PLATE's read-length distribution
  # -- it is the shortest-read dataset, so it is the window all four overlap in,
  # and it is the one demo_analyze already standardised to, which keeps the two
  # folders' numbers comparable. A dataset standardised to its own distribution
  # must reproduce its crude rate; the plate row does, which is the formula
  # checking itself.
  ref <- bl %>% filter(dataset == DATASETS[3]) %>% distinct(length, w = reads_at_length)
  std <- bl %>% inner_join(ref, by = "length") %>%
    group_by(biotype, dataset) %>%
    summarise(crude = 100 * sum(reads) / sum(reads_at_length),
              standardised = sum(w * pct) / sum(w), lengths = n(), .groups = "drop") %>%
    mutate(across(c(crude, standardised), ~round(.x, 1)))
  write_tsv(std, file.path(TAB, "cross", "step5_biotype_length_standardised.tsv"))
  cat("\nstandardised to the published plate's read-length distribution:\n")
  std %>% left_join(std %>% filter(dataset == DATASETS[3]) %>%
                      select(biotype, ref = standardised), by = "biotype") %>%
    mutate(vs_plate = round(standardised / ref, 1)) %>%
    select(biotype, dataset, crude, standardised, vs_plate) %>%
    arrange(biotype, dataset) %>% as.data.frame() %>% print(row.names = FALSE)
}
cat("\ndone.\n")
