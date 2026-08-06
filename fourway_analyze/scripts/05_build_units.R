#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# One row per unit across the four datasets: what it is, what group it is in,
# and how deep it is.
#
# WHAT A "UNIT" IS, AND WHY THE WORD IS NOT "CELL"
# ------------------------------------------------
# The four datasets do not count the same object. VASA's units are barcoded
# cells inside one library; FLASH-seq's are ten separate libraries in an input
# titration, each of which is bulk RNA, not a cell. Calling them all "cells"
# would make the reads-per-unit figure read as a depth-per-cell comparison,
# which it is not. They are units, and the figures say which is which.
#
# DEPTH IS MEASURED AT STEP 3's INPUT, ON EVERY SIDE
# --------------------------------------------------
# `reads_in` = reads entering rRNA depletion = post-trim reads. That is the
# first point in the pipeline where all four datasets have run identically:
# before it VASA has been demultiplexed off a barcode read and FLASH-seq has
# not, so "demultiplexed reads" has no FLASH-seq counterpart at all. Raw and
# post-trim counts are carried alongside where they exist, and the trim funnel
# is reported per protocol rather than pooled -- VASA trims in three passes and
# FLASH-seq in one, so a shared "pass 1 / pass 2" axis would be a fiction.
#
# GROUPS
#   VASA own       genotype x replicate, from demo_analyze's sample sheet
#                  (which resolves the design xlsx and hard-errors on a bad join)
#   VASA published species call, mouse wells only
#   FLASH-seq      input amount, 30 ng -> 30 pg, two replicates each
#
# A8 IS KEPT. FS_native_library_metadata.tsv marks ZHA8833A8 "exclude" for 18.3%
# human CALB1 contamination and A7 "caveat" at 3.6%. Both are processed and
# reported identically to the rest and flagged in the `qc` column -- the same
# rule this project already applies to VASA's four blank cells. Dropping a unit
# because it is odd is how an artefact stops being visible.
#
# Input:  tables/cross/insilico_depletion_per_unit.tsv  (scripts/01)
#         demo_analyze's sample sheets, FS_native_library_metadata.tsv,
#         the FLASH-seq cutadapt logs
# Output: tables/units.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(stringr)})
source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))

DEP <- file.path(TAB, "cross", "insilico_depletion_per_unit.tsv")
if (!file.exists(DEP)) stop("run scripts/01_insilico_depletion.sh first -- no ", DEP)
dep <- read_tsv(DEP, show_col_types = FALSE)

# --- VASA own: genotype from demo_analyze's resolved sheet -------------------
own <- read_tsv(file.path(DEMO, "own130", "sample_sheet.tsv"), show_col_types = FALSE) %>%
  transmute(unit = cell, group = genotype, replicate = as.character(replicate), qc = "")

# --- VASA published: species call, mouse wells only --------------------------
plate <- read_tsv(file.path(DEMO, "plate", "reads_per_cell.tsv"), show_col_types = FALSE) %>%
  filter(call == "mouse") %>%
  transmute(unit = well, group = "Mouse (mESC)", replicate = NA_character_, qc = "")

# --- FLASH-seq: input titration ---------------------------------------------
fs <- read_tsv(FSMETA, show_col_types = FALSE) %>%
  transmute(unit = library, group = input_amount, replicate = as.character(replicate),
            qc = ifelse(qc_verdict == "ok", "", qc_verdict), input_pg)
# Ordered by input, descending, so the group axis reads as a titration and not
# alphabetically ("1.5 ng" < "3 ng" < "30 ng" as strings).
FS_ORDER <- fs %>% distinct(group, input_pg) %>% arrange(desc(input_pg)) %>% pull(group)

groups <- bind_rows(
  mutate(own,   dataset = DATASETS[1]),
  mutate(own,   dataset = DATASETS[2]),   # own75 is the same 16 barcodes
  mutate(plate, dataset = DATASETS[3]),
  mutate(select(fs, -input_pg), dataset = DATASETS[4]))

# --- the FLASH-seq trim funnel, from its own cutadapt logs -------------------
# Its trim is ONE cutadapt call (Nextera adapter, q20, -m 20), not VASA's three
# passes, so this is reported as its own two numbers and never merged into
# VASA's pass-1/pass-2 table.
fslog <- file.path(VASA, "data/flashseq_vasa/run/native/logs")
fs_trim <- if (dir.exists(fslog)) {
  lapply(fs$unit, function(u) {
    f <- file.path(fslog, paste0("cutadapt_", u, ".log"))
    if (!file.exists(f)) return(NULL)
    txt <- readLines(f, warn = FALSE)
    num <- function(pat) {
      l <- grep(pat, txt, value = TRUE)[1]
      if (is.na(l)) return(NA_real_)
      as.numeric(gsub(",", "", str_extract(l, "[0-9,]+(?= *(bp|\\())|[0-9,]+$")))
    }
    tibble(dataset = DATASETS[4], unit = u,
           reads_raw = num("^Total reads processed:"),
           reads_trimmed = num("^Reads written \\(passing filters\\):"))
  }) %>% bind_rows()
} else tibble(dataset = character(), unit = character(),
              reads_raw = numeric(), reads_trimmed = numeric())

units <- dep %>%
  left_join(groups, by = c("dataset", "unit")) %>%
  left_join(fs_trim, by = c("dataset", "unit")) %>%
  mutate(dataset = factor(dataset, levels = DATASETS),
         group = factor(group, levels = c(
           "XY", "XO", "EpiLCs", "Blank control", "Mouse (mESC)", FS_ORDER))) %>%
  arrange(dataset, group, unit)

# A unit with no group is a join failure, and a join failure here just looks
# like an empty facet. Refuse instead.
bad <- units %>% filter(is.na(group))
if (nrow(bad)) {
  print(as.data.frame(bad))
  stop(nrow(bad), " units have no group -- the sample-sheet join failed")
}
if (any(is.na(units$dataset))) stop("a dataset label does not match paths.R's DATASETS")

write_tsv(units, file.path(TAB, "units.tsv"))
cat("wrote tables/units.tsv\n")
units %>% count(dataset, group, name = "units") %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nflagged units:\n")
units %>% filter(qc != "") %>% select(dataset, unit, group, qc) %>%
  as.data.frame() %>% print(row.names = FALSE)
