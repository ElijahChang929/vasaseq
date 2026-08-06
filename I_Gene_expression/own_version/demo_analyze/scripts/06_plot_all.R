#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_all.R
#
# Every figure for both datasets, from one script, in one style.
#
#   own      PM26037 / ZHA9292A1   16 barcodes, grouped by genotype
#   plate    SRR14783059           384 barcodes, grouped by species call
#
# Reads the tables written by ./build_tables.R and vasaplate/build_tables.R.
# Writes individual PNG+PDF into figures/ -- nothing is a composite, so any one
# of them can go into a slide on its own:
#
#   figures/own_reads_per_cell.*        reads per barcode, by genotype
#   figures/own_reads_per_group.*       mean +- SEM per genotype
#   figures/plate_reads_per_cell.*      reads per barcode, by species
#   figures/plate_reads_per_group.*     mean +- SEM per species
#   figures/both_reads_distribution.*   the two datasets side by side
#   figures/both_ufi_distribution.*     ditto, UFIs
#   figures/own_trim_loss_per_cell.*    where step 2 loses reads, per barcode
#   figures/own_trim_loss_per_genotype.* the same summed per genotype, labelled
#
# The two datasets used to have their own plot scripts, which meant restyling
# one silently left the other behind. One script, one theme_demo(), one save():
# a style change lands everywhere at once.
#
# NOTE the plate's species split uses the **Fig. 1d UFI-fraction** doublet rule.
# The paper's Methods gene-fraction rule gives a ~6x different answer on this
# library. It is set by RULE in vasaplate/build_tables.R and explained in
# vasaplate/README.md -- the figures do not carry it, so name it yourself
# whenever one of these leaves the folder.
#
# Run:
#   /nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/Rscript plot_all.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(forcats)
  library(ggplot2)
  library(scales)
})

source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))
OUT <- file.path(TAB, "own130")

OWN_LEVELS   <- c("XY", "XO", "EpiLCs", "Blank control")
# "Mixed (doublet)" is deliberately absent -- 2 wells, 0.17% of reads, dropped in
# vasaplate/build_tables.R (see DROP there for why and for the denominator caveat).
PLATE_LEVELS <- c("Human (HEK293T)", "Mouse (mESC)", "Below UFI cutoff")

# --- the one place style lives ---------------------------------------------

# base_size 14, not ggplot's 11: these are read at slide and print size, where
# the default axis text is too small. The in-panel geom_text sizes below are
# scaled to match -- ggplot's `size` is in mm and does NOT follow base_size, so
# raising one without the other leaves the labels looking shrunken.
BASE <- 14
LAB  <- 3.9        # in-panel value labels (was 3.0 at base_size 11)

theme_demo <- function() {
  theme_classic(base_size = BASE) +
    theme(plot.title = element_text(hjust = 0.5))
}

# `sub` picks the figures/ subdirectory: reads / trim / rrna / mapping / cross.
save_fig <- function(p, name, width, height, sub = "trim") {
  d <- file.path(FIG, FIGDIR[[sub]])
  ggsave(file.path(d, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(d, paste0(name, ".pdf")), p, width = width, height = height)
  cat("  figures/", FIGDIR[[sub]], "/", name, ".{png,pdf}\n", sep = "")
}

# --- load, and harmonise the two tables to the same column names -----------

own <- read_tsv(file.path(OUT, "reads_per_cell.tsv"), show_col_types = FALSE) %>%
  transmute(dataset = "own library",
            group   = factor(genotype, levels = OWN_LEVELS),
            replicate,
            xlabel  = paste0("seq", seq_no, "\n", cell),
            reads   = demultiplexed,
            ufi_total) %>%
  arrange(group, replicate)

plate <- read_tsv(file.path(TAB, "plate", "reads_per_cell.tsv"), show_col_types = FALSE) %>%
  transmute(dataset = "published VASA-plate",
            group   = factor(group, levels = PLATE_LEVELS),
            replicate = well_n,
            xlabel  = well,
            reads, ufi_total)

both <- bind_rows(own, plate) %>%
  mutate(dataset = factor(dataset, levels = c("own library", "published VASA-plate")))

# --- shared plot builders --------------------------------------------------

# mean +- SEM with every barcode as a point. Same code for 4 replicates and for
# 178 -- only the point size changes.
p_group <- function(d, ptsize = 2) {
  s <- d %>%
    group_by(group) %>%
    summarise(mean_reads = mean(reads),
              sem_reads  = sd(reads) / sqrt(n()), .groups = "drop")
  ggplot(s, aes(group, mean_reads / 1e6, fill = group)) +
    geom_col() +
    geom_errorbar(aes(ymin = (mean_reads - sem_reads) / 1e6,
                      ymax = (mean_reads + sem_reads) / 1e6), width = 0.2) +
    geom_point(data = d, aes(group, reads / 1e6), inherit.aes = FALSE,
               position = position_jitter(width = 0.15, height = 0, seed = 1),
               size = ptsize, alpha = if (ptsize < 1.5) 0.5 else 1) +
    labs(x = NULL, y = "Reads (millions)", fill = NULL) +
    theme_demo()
}

# box + every point, log scale, one panel per dataset. Used for reads and UFIs.
p_distribution <- function(d, value, ylab, title) {
  v <- rlang::ensym(value)
  lab <- d %>% group_by(dataset, group) %>%
    summarise(n = n(), y = max(!!v), .groups = "drop")
  ggplot(d, aes(group, !!v, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.55) +
    geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.6) +
    geom_text(data = lab, aes(group, y, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -1, size = LAB) +
    facet_wrap(~ dataset, scales = "free_x") +
    scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0.06, 0.10))) +
    labs(x = NULL, y = ylab, title = title) +
    theme_demo() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

cat("wrote:\n")

# --- own: reads per barcode ------------------------------------------------
# 16 bars, so every bar is labelled with its design well and cell id.

p <- ggplot(own, aes(fct_inorder(xlabel), reads / 1e6, fill = group)) +
  geom_col() +
  labs(x = NULL, y = "Reads (millions)", fill = NULL,
       title = "Reads per cell barcode") +
  theme_demo()
save_fig(p, "own_reads_per_cell", 9, 4, sub = "reads")

save_fig(p_group(own), "own_reads_per_group", sub = "reads", 6, 4)

# --- plate: reads per barcode ----------------------------------------------
# 384 bars cannot be labelled, so they are sorted by reads inside each group and
# the group carries its n in the strip. Equal panel widths, not
# space = "free_x": the n = 31 group otherwise collapses to a sliver.

plate_sorted <- plate %>%
  group_by(group) %>%
  arrange(desc(reads), .by_group = TRUE) %>%
  mutate(rank = row_number(), n = n()) %>%
  ungroup() %>%
  mutate(panel = factor(paste0(group, "\nn = ", n),
                        levels = unique(paste0(group, "\nn = ", n)[order(group)])))

p <- ggplot(plate_sorted, aes(rank, reads / 1e6, fill = group)) +
  geom_col() +
  facet_grid(~ panel, scales = "free_x") +
  labs(x = "Barcode, sorted by reads within group",
       y = "Reads (millions)", fill = NULL,
       title = "Reads per cell barcode") +
  theme_demo() +
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p, "plate_reads_per_cell", 10, 4, sub = "reads")

save_fig(p_group(plate, ptsize = 0.8), "plate_reads_per_group", sub = "reads", 7, 4)
# --- own-style figures, for both the full-length run and the 75 nt one -------
# `own_composition(dir, prefix)` draws the whole own set from one directory of
# tables. demo_analyze/ holds the 130 nt run; demo_analyze/own75/ holds the same
# library truncated to 75 nt of biological sequence to match the plate. Calling
# it twice is what keeps the two comparable -- one code path, one style.
SEG <- list(colour = "white", linewidth = 0.4)   # the border between segments

stacked_abs <- function(d, xvar, fillvar, ylab, title, legend_rows = 1) {
  ggplot(d, aes({{ xvar }}, reads / 1e6, fill = {{ fillvar }})) +
    geom_col(position = position_stack(reverse = TRUE),
             colour = SEG$colour, linewidth = SEG$linewidth) +
    labs(x = NULL, y = ylab, fill = NULL, title = title) +
    theme_demo() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = legend_rows))
}

stacked_pct <- function(d, fillvar, totals, title, min_label = 4, legend_rows = 1) {
  ggplot(d, aes(genotype, reads, fill = {{ fillvar }})) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.62,
             colour = SEG$colour, linewidth = SEG$linewidth) +
    geom_text(aes(label = ifelse(pct >= min_label, sprintf("%.1fM\n%.1f%%", reads / 1e6, pct), "")),
              position = position_fill(vjust = 0.5, reverse = TRUE),
              size = LAB, lineheight = 0.95) +
    geom_text(data = totals, aes(genotype, 1.05, label = sprintf("%.1fM", tot / 1e6)),
              inherit.aes = FALSE, size = LAB) +
    scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0.02, 0.09))) +
    labs(x = NULL, y = "Share of demultiplexed reads", fill = NULL, title = title) +
    theme_demo() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = legend_rows))
}

own_composition <- function(dir, prefix) {
  OUT <- dir            # the helpers below all read tables relative to OUT
  sf <- function(p, n, w, h) save_fig(p, paste0(prefix, n), w, h)

  # --- own: the composition figures -------------------------------------------
  # One idiom for all three: the CATEGORY gets the fill (distinct hues), and a 2px
  # white border separates neighbouring segments. Genotype is read off the x axis
  # or the facet strip instead of the fill.
  #
  # This replaced encoding the category on alpha over the genotype hue. That works
  # at three levels and fails at five -- the bottom of the ramp (0.29 / 0.26 /
  # 0.14) is not separable, which is exactly where the small classes live. Hatching
  # was considered instead and rejected: ggpattern is not installed, and a texture
  # inside a 2%-tall segment reads as noise at print size.
  #
  # Genotype keeps its own palette in every figure where genotype IS the variable
  # (the reads and UFI figures) -- it is only displaced where a second categorical
  # dimension has to be shown at the same time.
  #
  # Pass 0 (the barcode anchor) is absent from these on purpose: it has no length
  # filter, so it removes bases and never a whole read. Its base-level cost is in
  # trim_loss_per_cell.tsv as `anchor_bases`.


  FATES <- c("Kept", "Lost at pass 2 (cutadapt -m)", "Lost at pass 1 (TrimGalore --length)")

  tl_cell <- read_tsv(file.path(OUT, "trim_loss_per_cell_long.tsv"), show_col_types = FALSE) %>%
    mutate(genotype = factor(genotype, levels = OWN_LEVELS),
           fate     = factor(sub(" 20\\)$", ")", fate), levels = FATES),
           label    = paste0("seq", seq_no, "\n", cell)) %>%
    arrange(genotype, replicate) %>%
    mutate(label = fct_inorder(label))
  stopifnot(!any(is.na(tl_cell$fate)), !any(is.na(tl_cell$genotype)))

  p <- stacked_abs(tl_cell, label, fate, "Reads (millions)",
                   "Read loss through step 2, per cell barcode") +
    facet_grid(~ genotype, scales = "free_x")
  sf(p, "trim_loss_per_cell", 11, 5)

  tl_geno <- read_tsv(file.path(OUT, "trim_loss_per_genotype_long.tsv"), show_col_types = FALSE) %>%
    mutate(genotype = factor(genotype, levels = OWN_LEVELS),
           fate     = factor(sub(" 20\\)$", ")", fate), levels = FATES))
  tl_tot <- tl_geno %>% distinct(genotype, tot = demultiplexed)
  sf(stacked_pct(tl_geno, fate, tl_tot, "Read loss through step 2, by genotype"),
           "trim_loss_per_genotype", 7.5, 5.4)

  # --- the five read classes ---
  CLASSES <- c("1 read-through, kept", "5 no read-through, kept", "2 insert too short",
               "3 pass 1 trimmed it short", "4 pass 2 trimmed it short")

  cls_f <- file.path(OUT, "read_classes_per_cell.tsv")
  if (file.exists(cls_f)) {
    rc_cell <- read_tsv(cls_f, show_col_types = FALSE) %>%
      mutate(genotype = factor(genotype, levels = OWN_LEVELS),
             class    = factor(class, levels = CLASSES),
             label    = paste0("seq", seq_no, "\n", cell)) %>%
      arrange(genotype, replicate) %>%
      mutate(label = fct_inorder(label))
    stopifnot(!any(is.na(rc_cell$class)), !any(is.na(rc_cell$genotype)))

    p <- stacked_abs(rc_cell, label, class, "Reads (millions)",
                     "Read classes, per cell barcode", legend_rows = 2) +
      facet_grid(~ genotype, scales = "free_x")
    sf(p, "read_classes_per_cell", 11, 5.6)

    rc_geno <- read_tsv(file.path(OUT, "read_classes_per_genotype.tsv"), show_col_types = FALSE) %>%
      mutate(genotype = factor(genotype, levels = OWN_LEVELS),
             class    = factor(class, levels = CLASSES))
    rc_tot <- rc_geno %>% group_by(genotype) %>% summarise(tot = sum(reads), .groups = "drop")
    sf(stacked_pct(rc_geno, class, rc_tot, "Read classes, by genotype", legend_rows = 2),
             "read_classes_per_genotype", 8, 6)
  } else {
    cat("  (skipping read-class figures -- run own_version/diagnostics/classify_reads.py first)\n")
  }

  # --- what cut the reads pass 2 dropped ---
  ADAPTERS <- c("polyA", "polyT5", "rt", "polyG", "no_adapter")

  att_f <- file.path(OUT, "pass2_attribution_per_genotype.tsv")
  if (file.exists(att_f)) {
    att <- read_tsv(att_f, show_col_types = FALSE) %>%
      mutate(genotype = factor(genotype, levels = OWN_LEVELS),
             adapter  = factor(adapter, levels = ADAPTERS))
    stopifnot(!any(is.na(att$adapter)), !any(is.na(att$genotype)))
    att_tot <- att %>% group_by(genotype) %>% summarise(tot = sum(reads), .groups = "drop")

    p <- stacked_pct(att, adapter, att_tot, "What cut the reads pass 2 dropped",
                     min_label = 5) +
      labs(y = "Share of reads dropped at pass 2")
    sf(p, "pass2_cause_per_genotype", 8, 5.6)
  } else {
    cat("  (skipping pass-2 cause figure -- run own_version/diagnostics/attribute_pass2_loss.py first)\n")
  }

  # --- length of what survives step 2 ---
  # Normalised within genotype: the depths differ ~20x, so on raw counts the
  # blanks are a flat line and only EpiLCs has a visible shape. The question is
  # the shape.
  len_f <- file.path(OUT, "read_length_per_genotype.tsv")
  if (file.exists(len_f)) {
    lg <- read_tsv(len_f, show_col_types = FALSE) %>%
      mutate(genotype = factor(genotype, levels = OWN_LEVELS))
    st <- read_tsv(file.path(OUT, "read_length_stats.tsv"), show_col_types = FALSE) %>%
      mutate(genotype = factor(genotype, levels = OWN_LEVELS))

    ML <- st$minlen[1]                    # the -m floor the run actually used
    PRE <- isTRUE(st$prefilter[1])
    # Shade what the floor throws away. On the pre-filter tally that region is
    # real data; on the post-filter one it is empty by construction, so say which
    # is being shown rather than letting an empty band look like a finding.
    sub <- if (PRE) sprintf("shaded: < %d nt, discarded by cutadapt -m %d", ML, ML)
           else "POST-filter tally -- the discarded tail is missing"

    p <- ggplot(lg, aes(length, pct, colour = genotype)) +
      annotate("rect", xmin = -Inf, xmax = ML, ymin = -Inf, ymax = Inf,
               fill = "grey80", alpha = 0.45) +
      geom_line(linewidth = 0.7) +
      geom_vline(xintercept = ML, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      annotate("text", x = ML + 2, y = Inf, label = sub,
               hjust = 0, vjust = 1.6, size = 3.6, colour = "grey30") +
      scale_x_continuous(breaks = c(0, ML, seq(20, 130, 20))) +
      labs(x = "Read length once trimming can remove nothing more (nt)",
           y = "% of that genotype's reads entering pass 2",
           colour = NULL, title = "Fragment length after trimming") +
      theme_demo() +
      theme(legend.position = "bottom")
    sf(p, "read_length", 9.5, 5)

    # how much of each genotype the floor removes, which the curve only implies
    p <- ggplot(st, aes(genotype, pct_discarded, fill = genotype)) +
      geom_col(width = 0.62) +
      geom_text(aes(label = sprintf("%.1f%%", pct_discarded)), vjust = -0.4, size = LAB) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = NULL, y = sprintf("%% of reads shorter than %d nt", ML), fill = NULL,
           title = sprintf("Discarded as too short (cutadapt -m %d)", ML)) +
      theme_demo() + theme(legend.position = "none")
    sf(p, "read_length_discarded", 6.5, 4.5)

    p <- ggplot(st, aes(genotype, median, fill = genotype)) +
      geom_col(width = 0.62) +
      geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.2) +
      geom_text(aes(y = q3, label = sprintf("%d nt\nIQR %d-%d", median, q1, q3)),
                vjust = -0.4, size = LAB, lineheight = 0.95) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(x = NULL, y = "Read length after trimming (nt)", fill = NULL,
           title = "Median read length after trimming (bar = median, whisker = IQR)") +
      theme_demo() + theme(legend.position = "none")
    sf(p, "read_length_median", 6.5, 4.5)
  } else {
    cat("  (skipping length figures -- run own_version/diagnostics/read_length_dist.sh first)\n")
  }
}

# 130 nt, the run as sequenced
own_composition(OUT, "own_")

# 75 nt: the same library truncated to the plate's read length, so that any
# remaining difference against the plate is not read length. Absent until the
# out75/ pipeline has finished, hence the guard rather than an error.
OUT75 <- file.path(TAB, "own75")
if (file.exists(file.path(OUT75, "trim_loss_per_cell_long.tsv"))) {
  own_composition(OUT75, "own75_")
} else {
  cat("  (skipping own75 figures -- out75 pipeline not built yet)\n")
}

# --- plate: the same two, under the UPSTREAM pipeline ------------------------
# a_Mapping/trim.sh, i.e. what the paper ran: no pass 0 at all, and pass 2 is
# four 6-mer homopolymers at -n 1. The five-class figure has no counterpart --
# classes 1 and 5 split on reading through into the cell barcode, which upstream
# never looks for, and which barely happens here anyway (0.3-0.5% of these 75 nt
# reads, against ~25% of the own library's 130 nt ones).
#
# stacked_pct() is written against a column called `genotype`, so the plate's
# `group` is renamed rather than the function being duplicated.

plate_pct <- function(d, fillvar, totals, title, ylab, min_label = 4, rows = 1) {
  stacked_pct(rename(d, genotype = group), {{ fillvar }},
              rename(totals, genotype = group), title,
              min_label = min_label, legend_rows = rows) +
    labs(y = ylab)
}

tlp_f <- file.path(TAB, "plate", "trim_loss_per_group_long.tsv")
if (file.exists(tlp_f)) {
  tlp <- read_tsv(tlp_f, show_col_types = FALSE) %>%
    mutate(group = factor(group, levels = PLATE_LEVELS),
           fate  = factor(fate, levels = unique(fate)))
  tlp_tot <- tlp %>% distinct(group, tot = demultiplexed)
  save_fig(plate_pct(tlp, fate, tlp_tot,
                     "Read loss through trimming, by species (upstream pipeline)",
                     "Share of demultiplexed reads"),
           "plate_trim_loss_per_group", 8, 5.4)
} else {
  cat("  (skipping plate trim-loss figure -- vasaplate/build_tables.R found no counts)\n")
}

lenp_f <- file.path(TAB, "plate", "read_length_per_group.tsv")
if (file.exists(lenp_f)) {
  lgp <- read_tsv(lenp_f, show_col_types = FALSE) %>%
    mutate(group = factor(group, levels = PLATE_LEVELS))
  stp <- read_tsv(file.path(TAB, "plate", "read_length_stats.tsv"),
                  show_col_types = FALSE) %>%
    mutate(group = factor(group, levels = PLATE_LEVELS))
  MLp <- stp$minlen[1]

  p <- ggplot(lgp, aes(length, pct, colour = group)) +
    annotate("rect", xmin = -Inf, xmax = MLp, ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.45) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = MLp, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    annotate("text", x = MLp + 2, y = Inf,
             label = sprintf("shaded: < %d nt, discarded by cutadapt -m %d", MLp, MLp),
             hjust = 0, vjust = 1.6, size = 3.6, colour = "grey30") +
    scale_x_continuous(breaks = c(0, MLp, seq(20, 80, 20))) +
    labs(x = "Read length once trimming can remove nothing more (nt)",
         y = "% of that group's reads entering pass 2", colour = NULL,
         title = "Fragment length after trimming (upstream pipeline)") +
    theme_demo() + theme(legend.position = "bottom")
  save_fig(p, "plate_read_length", 9.5, 5)

  p <- ggplot(stp, aes(group, pct_discarded, fill = group)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = sprintf("%.1f%%", pct_discarded)), vjust = -0.4, size = LAB) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = sprintf("%% of reads shorter than %d nt", MLp), fill = NULL,
         title = sprintf("Discarded as too short (cutadapt -m %d)", MLp)) +
    theme_demo() + theme(legend.position = "none")
  save_fig(p, "plate_read_length_discarded", 6.5, 4.5)
} else {
  cat("  (skipping plate length figures -- no read_length_per_group.tsv)\n")
}

attp_f <- file.path(TAB, "plate", "pass2_attribution_per_group.tsv")
if (file.exists(attp_f)) {
  LEGACY_ADAPTERS <- c("polyA1", "polyT1", "polyG1", "polyC1", "no_adapter")
  attp <- read_tsv(attp_f, show_col_types = FALSE) %>%
    mutate(group   = factor(group, levels = PLATE_LEVELS),
           adapter = factor(adapter, levels = LEGACY_ADAPTERS))
  stopifnot(!any(is.na(attp$adapter)))
  attp_tot <- attp %>% group_by(group) %>% summarise(tot = sum(reads), .groups = "drop")
  save_fig(plate_pct(attp, adapter, attp_tot,
                     "What cut the reads pass 2 dropped (upstream pipeline)",
                     "Share of reads dropped at pass 2", min_label = 5),
           "plate_pass2_cause_per_group", 8, 5.6)
} else {
  cat("  (skipping plate pass-2 cause figure -- no attribution table)\n")
}

# --- both: the comparison figures ------------------------------------------

# --- all three fragment-length curves on one axis ---------------------------
# The three separate figures cannot be compared by eye: their y-axes top out at
# 25%, 27% and 44%. Pooled to one curve per dataset (cells summed, then
# normalised) so the shapes sit on the same scale.

len_sources <- list(
  c(f = file.path(OUT, "read_length_per_genotype.tsv"),          lab = "own library, 130 nt"),
  c(f = file.path(TAB, "own75", "read_length_per_genotype.tsv"), lab = "own library, 75 nt"),
  c(f = file.path(TAB, "plate", "read_length_per_group.tsv"),lab = "published VASA-plate, 75 nt")
)
have <- vapply(len_sources, function(s) file.exists(s[["f"]]), logical(1))
if (all(have)) {
  pooled <- bind_rows(lapply(len_sources, function(s)
    read_tsv(s[["f"]], show_col_types = FALSE) %>%
      group_by(length) %>% summarise(reads = sum(reads), .groups = "drop") %>%
      mutate(dataset = s[["lab"]], pct = 100 * reads / sum(reads)))) %>%
    mutate(dataset = factor(dataset, levels = vapply(len_sources, `[[`, "", "lab")))

  ML <- 15
  p <- ggplot(pooled, aes(length, pct, colour = dataset)) +
    annotate("rect", xmin = -Inf, xmax = ML, ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.45) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = ML, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    annotate("text", x = ML + 2, y = Inf,
             label = sprintf("shaded: < %d nt, discarded by cutadapt -m %d", ML, ML),
             hjust = 0, vjust = 1.6, size = 3.6, colour = "grey30") +
    scale_x_continuous(breaks = c(0, ML, seq(25, 130, 25))) +
    labs(x = "Read length once trimming can remove nothing more (nt)",
         y = "% of that library's reads entering pass 2", colour = NULL,
         title = "Fragment length after trimming") +
    theme_demo() +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 1))
  save_fig(p, "both_read_length", 10, 5.5, sub = "cross")

  # --- the same three distributions, binned ---------------------------------
  # The per-nt curve above hides two things. (1) Each curve is a distribution
  # summing to 100%, but with 131 points at 0-20% each you cannot read any
  # share off it. (2) Below 15 nt the 130 nt and 75 nt curves of the *same*
  # library sit on top of each other (6.7% vs 7.2% at length 0), so the one
  # drawn last hides the other -- and it reads as if the 130 nt library had no
  # short reads at all, when in fact 47.0% of its reads are < 15 nt (vs 47.8%
  # and 10.0%). Binning fixes both: eight bars per library, summing to 100%.
  #
  # Bins are absolute nt, except the last: "untrimmed" is the library's own
  # full read length (130 or 75), pulled out of its bin because a read that
  # trimming never touched is a different thing from a merely long one, and
  # because that bar is the one directly comparable across the three.
  # 15 = cutadapt -m; 30 = below it bwa mem reports nothing (-T 30, see
  # own_version/README.md), so the first two bars are both "lost", differently.
  LEN_BINS <- c("0-14", "15-29", "30-49", "50-69",
                "70-89", "90-109", "110-129", "untrimmed")
  binned <- pooled %>%
    group_by(dataset) %>% mutate(full = max(length)) %>% ungroup() %>%
    mutate(bin = ifelse(length == full, "untrimmed",
                        as.character(cut(length, c(-1, 14, 29, 49, 69, 89, 109, 129),
                                         labels = LEN_BINS[1:7]))),
           bin = factor(bin, levels = LEN_BINS)) %>%
    group_by(dataset, full, bin, .drop = FALSE) %>%
    summarise(pct = sum(pct), .groups = "drop") %>%
    # a 75 nt library has no 90-109 nt read to have; that is not 0%, it is not
    # a measurable bin, so drop the bar rather than draw an empty one.
    mutate(lo = c(0, 15, 30, 50, 70, 90, 110, NA)[as.integer(bin)]) %>%
    filter(is.na(lo) | lo < full) %>% select(-lo, -full)
  stopifnot(all(abs(tapply(binned$pct, binned$dataset, sum) - 100) < 1e-6))

  p <- ggplot(binned, aes(bin, pct, fill = dataset)) +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.45) +
    # position_dodge, not dodge2: where a bin has only the 130 nt library in it
    # dodge2 would centre that lone bar in the bin, dodge keeps it in the first
    # slot -- so a bar's horizontal position means the same thing in every bin.
    geom_col(position = position_dodge(width = 0.8, preserve = "single"),
             width = 0.72) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              position = position_dodge(width = 0.8, preserve = "single"),
              vjust = -0.4, size = LAB - 0.9) +
    annotate("text", x = 1.6, y = Inf,
             label = sprintf("shaded: < %d nt, discarded by cutadapt -m %d", ML, ML),
             hjust = 0, vjust = 1.6, size = 3.6, colour = "grey30") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Read length once trimming can remove nothing more (nt)",
         y = "% of that library's reads entering pass 2", fill = NULL,
         title = "Fragment length after trimming, binned",
         subtitle = "each library's eight bars sum to 100%") +
    theme_demo() +
    theme(plot.subtitle = element_text(hjust = 0.5, colour = "grey30"),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1))
  save_fig(p, "both_read_length_binned", 10, 5.5, sub = "cross")
  write_tsv(binned, file.path(TAB, "cross", "read_length_binned.tsv"))
  cat("  tables/cross/read_length_binned.tsv\n")
} else {
  cat("  (skipping combined length figure -- missing",
      paste(vapply(len_sources[!have], `[[`, "", "lab"), collapse = ", "), ")\n")
}

save_fig(p_distribution(both, reads, "Reads per cell barcode (log scale)",
                        "Read distribution per cell barcode"),
         "both_reads_distribution", 10, 5, sub = "cross")

save_fig(p_distribution(both, ufi_total, "UFIs per cell barcode (log scale)",
                        "UFI distribution per cell barcode"),
         "both_ufi_distribution", 10, 5, sub = "cross")

# --- summary table across both datasets ------------------------------------

summ <- both %>%
  group_by(dataset, group) %>%
  summarise(n = n(),
            median_reads = median(reads), mean_reads = mean(reads),
            median_ufi   = median(ufi_total), mean_ufi = mean(ufi_total),
            .groups = "drop")
write_tsv(summ, file.path(OUT, "summary_both_datasets.tsv"))
cat("  summary_both_datasets.tsv\n")

print(as.data.frame(summ %>%
  transmute(dataset, group, n,
            median_reads = comma(median_reads), median_ufi = comma(median_ufi))),
  row.names = FALSE)
