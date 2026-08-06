#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_tables.R  (VASA-plate, published library)
#
# The same measurement as ../build_tables.R, on the published plate library
# SRR14783059 / GSM5369495 (vasaplate-HEK293T-mESC): reads per cell barcode,
# grouped.
#
# There is no genotype design here -- this library is a species-mixing control,
# 384 CEL-seq2 barcodes of HEK293T (human) and mESC (mouse). The grouping that
# plays the role "genotype" played in PM26037 is therefore the **species call**,
# taken from res/vasaplate/per_cell.tsv (source `ours_v3`, the anchor run):
#
#   human / mouse   the two cell types actually plated
#   mixed           heterotypic doublet
#   discarded       below the 7,500-UFI cutoff -- empty or near-empty wells,
#                   the closest thing this library has to a blank control
#
# IMPORTANT: the paper states two doublet rules that disagree ~6x on this data
# (Fig. 1d thresholds UFI fraction, Methods thresholds gene fraction: 0.85% vs
# 4.82%). This uses **call_fig1d**, and the figures say so. Never quote a
# barnyard number from here without naming the rule.
#
# Counts come from demux_read_counts.tsv (count_demux_reads.sh -- exact
# `wc -l / 4` over the step-1 per-cell FASTQs).
#
# Run:
#   /nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/Rscript build_tables.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(scales)
  library(tidyr)
})

source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))
OUT <- file.path(TAB, "plate")

PER_CELL <- "/nemo/lab/turnerj/working/guangxin/vasaseq/res/vasaplate/per_cell.tsv"
FQDIR    <- "/nemo/lab/turnerj/working/guangxin/vasaseq/data/ref/fastq_vasaplate"
SOURCE   <- "ours_v3"          # the anchor run; see vasaplate_check/BENCHMARK_published_plate.md
RULE     <- "call_fig1d"       # UFI-fraction rule, the paper's Fig. 1d

GROUPS <- c(human     = "Human (HEK293T)",
            mouse     = "Mouse (mESC)",
            mixed     = "Mixed (doublet)",
            discarded = "Below UFI cutoff")

# Wells excluded from every table and figure downstream. `mixed` is 2 of 384
# wells and 0.17% of the reads -- too few to say anything, and at that size its
# bars carry "0.0M" labels that read as an error rather than as a number.
#
# NOTE this makes the plate tables cover 382 wells, not 384, and every
# percentage below is a share of those 382. It is NOT the same denominator the
# barnyard numbers in vasaplate_check/ use, which keep all 384 -- the doublet
# rate is the whole point there. Do not carry a percentage across.
DROP <- "mixed"

# --- 1. reads per barcode --------------------------------------------------

demux_f <- file.path(OUT, "demux_read_counts.tsv")
if (!file.exists(demux_f)) stop("missing ", demux_f, " -- run ./count_demux_reads.sh first")

demux <- read_tsv(demux_f, col_types = "cd") %>%
  mutate(well = str_pad(well, 3, pad = "0"))

# --- 2. species call per barcode -------------------------------------------

calls <- read_tsv(PER_CELL, col_types = cols(well = col_character(), .default = col_guess())) %>%
  filter(source == SOURCE) %>%
  mutate(well = str_pad(well, 3, pad = "0")) %>%
  select(well, ufi_total, frac_ufi_human, call = all_of(RULE))

per_cell <- demux %>%
  left_join(calls, by = "well") %>%
  mutate(group = factor(unname(GROUPS[call]), levels = unname(GROUPS)),
         well_n = as.integer(well)) %>%
  arrange(well_n)

# A barcode with reads but no call means the two tables disagree about which
# wells exist -- fail rather than plot a grey bar.
if (any(is.na(per_cell$group))) {
  stop("no species call for well(s): ",
       paste(per_cell$well[is.na(per_cell$group)], collapse = ", "))
}

# Drop the excluded calls AFTER the NA check above, so a well with no call at
# all still errors rather than being silently swept out with them. Percentages
# are computed after the drop, i.e. over what remains.
if (length(DROP)) {
  gone <- per_cell %>% filter(call %in% DROP)
  if (nrow(gone)) {
    message(sprintf("dropping %d well(s) called %s: %s reads (%.2f%% of the library)",
                    nrow(gone), paste(DROP, collapse = "/"),
                    format(sum(gone$reads), big.mark = ","),
                    100 * sum(gone$reads) / sum(per_cell$reads)))
  }
  per_cell <- per_cell %>% filter(!call %in% DROP) %>%
    mutate(group = droplevels(group))
}

per_cell <- per_cell %>%
  mutate(pct_of_library = 100 * reads / sum(reads)) %>%
  select(well, well_n, group, call, reads, pct_of_library, ufi_total, frac_ufi_human)

per_group <- per_cell %>%
  group_by(group) %>%
  summarise(n_wells        = n(),
            total_reads    = sum(reads),
            mean_reads     = mean(reads),
            sd_reads       = sd(reads),
            sem_reads      = sd(reads) / sqrt(n()),
            median_reads   = median(reads),
            min_reads      = min(reads),
            max_reads      = max(reads),
            pct_of_library = 100 * sum(reads) / sum(per_cell$reads),
            .groups = "drop")

write_tsv(per_cell,  file.path(OUT, "reads_per_cell.tsv"))
write_tsv(per_group, file.path(OUT, "reads_per_group.tsv"))

# --- 3. trimming, the upstream way -----------------------------------------
#
# This library was trimmed by a_Mapping/trim.sh, the published pipeline, which
# is NOT what own_version runs. Two differences shape everything below:
#
#   * there is no pass 0. Upstream never anchors on the cell barcode, so the
#     own-library split into "read through into its own barcode" and "did not"
#     has no counterpart here. It would also be pointless: the plate reads are
#     75 nt, and only 0.3-0.5% of them reach their own barcode at all, against
#     ~25% in the 130 nt own library. Read-through is an own-library problem.
#   * pass 2 is four 6-mer homopolymers (polyG1/polyC1/polyT1/polyA1) at -n 1,
#     not the measured adapter + 20-mers + 5' poly-T + --poly-a.
#
# So the plate gets the two figures that do transfer -- where reads are lost,
# and what cut the ones pass 2 dropped -- and not the five-class one.

FATES <- c(kept       = "Kept",
           lost_pass2 = "Lost at pass 2 (cutadapt -m 15)",
           lost_pass1 = "Lost at pass 1 (TrimGalore)")

stage_f <- file.path(FQDIR, "trim_stage_counts.tsv")
att_f   <- file.path(FQDIR, "pass2_adapter_attribution.tsv")

if (file.exists(stage_f)) {
  stage <- read_tsv(stage_f, col_types = "cddd") %>%
    mutate(well = str_pad(well, 3, pad = "0")) %>%
    inner_join(select(per_cell, well, group), by = "well") %>%
    mutate(lost_pass1 = demultiplexed - after_pass1,
           lost_pass2 = after_pass1 - after_pass2,
           kept       = after_pass2)
  stopifnot(all(with(stage, kept + lost_pass1 + lost_pass2 == demultiplexed)))

  trim_group <- stage %>%
    group_by(group) %>%
    summarise(demultiplexed = sum(demultiplexed),
              across(all_of(names(FATES)), sum), .groups = "drop")

  trim_group_long <- trim_group %>%
    pivot_longer(all_of(names(FATES)), names_to = "fate", values_to = "reads") %>%
    mutate(fate = factor(unname(FATES[fate]), levels = unname(FATES)),
           pct  = 100 * reads / demultiplexed)

  write_tsv(stage,           file.path(OUT, "trim_loss_per_cell.tsv"))
  write_tsv(trim_group_long, file.path(OUT, "trim_loss_per_group_long.tsv"))

  cat("\n== trimming loss per group (upstream pipeline) ==\n")
  # percentages first: formatting `demultiplexed` with comma() inside the same
  # transmute turns it into a string, and every later division silently breaks
  print(as.data.frame(trim_group %>%
    mutate(p1 = 100 * lost_pass1 / demultiplexed,
           p2 = 100 * lost_pass2 / demultiplexed,
           pk = 100 * kept       / demultiplexed) %>%
    transmute(group, demultiplexed = comma(demultiplexed),
              lost_pass1 = sprintf("%s (%.1f%%)", comma(lost_pass1), p1),
              lost_pass2 = sprintf("%s (%.1f%%)", comma(lost_pass2), p2),
              kept       = sprintf("%s (%.1f%%)", comma(kept), pk))),
    row.names = FALSE)
} else {
  message("no ", stage_f, " -- trim-loss tables skipped")
}

# --- fragment length, before the -m filter ----------------------------------
# From own_version/diagnostics/read_length_dist.sh --prefilter --mode legacy, i.e. re-derived with
# the FOUR 6-mer homopolymers the published pipeline used, not own_version's set.
MINLEN <- 15L
len_f <- file.path(FQDIR, "read_length_prefilter.tsv")
if (file.exists(len_f)) {
  len <- read_tsv(len_f, col_types = "cdd") %>%
    mutate(well = str_pad(cell, 3, pad = "0")) %>%
    inner_join(select(per_cell, well, group), by = "well")

  len_group <- len %>%
    group_by(group, length) %>%
    summarise(reads = sum(reads), .groups = "drop_last") %>%
    mutate(pct = 100 * reads / sum(reads)) %>%
    ungroup()

  q_at <- function(len, n, p) { o <- order(len); c <- cumsum(n[o]); sort(len)[which(c >= sum(n)*p)[1]] }
  len_stats <- len %>%
    group_by(group) %>%
    summarise(n = sum(reads),
              q1 = q_at(length, reads, 0.25),
              median = q_at(length, reads, 0.50),
              q3 = q_at(length, reads, 0.75),
              pct_discarded = 100 * sum(reads[length < MINLEN]) / sum(reads),
              .groups = "drop") %>%
    mutate(prefilter = TRUE, minlen = MINLEN)

  write_tsv(len_group, file.path(OUT, "read_length_per_group.tsv"))
  write_tsv(len_stats, file.path(OUT, "read_length_stats.tsv"))

  cat("\n== fragment length after trimming (upstream pipeline) ==\n")
  print(as.data.frame(len_stats %>%
    transmute(group, reads = comma(n), Q1 = q1, median, Q3 = q3,
              `discarded <15nt` = sprintf("%.1f%%", pct_discarded))), row.names = FALSE)
} else {
  message("no ", len_f, " -- length tables skipped")
}

if (file.exists(att_f)) {
  att <- read_tsv(att_f, col_types = "ccddd") %>%
    mutate(well = str_pad(cell, 3, pad = "0")) %>%
    inner_join(select(per_cell, well, group), by = "well") %>%
    group_by(group, adapter) %>%
    summarise(reads = sum(reads), bases_removed = sum(bases_removed),
              .groups = "drop_last") %>%
    mutate(pct = 100 * reads / sum(reads)) %>%
    ungroup()
  write_tsv(att, file.path(OUT, "pass2_attribution_per_group.tsv"))
} else {
  message("no ", att_f, " -- attribution table skipped")
}

# --- 3. console summary ----------------------------------------------------

cat("\nlibrary : SRR14783059 / GSM5369495 (vasaplate-HEK293T-mESC)\n")
cat("source  : ", SOURCE, "   doublet rule: ", RULE, " (UFI fraction, paper Fig. 1d)\n", sep = "")
cat("total   : ", comma(sum(per_cell$reads)), " reads over ",
    nrow(per_cell), " barcodes\n\n", sep = "")
print(as.data.frame(per_group %>%
  transmute(group, n_wells,
            total = comma(total_reads), mean = comma(round(mean_reads)),
            sd = comma(round(sd_reads)), median = comma(median_reads),
            min = comma(min_reads), max = comma(max_reads),
            pct = sprintf("%.2f%%", pct_of_library))),
  row.names = FALSE)
cat("\nwrote: reads_per_cell.tsv, reads_per_group.tsv\n")
