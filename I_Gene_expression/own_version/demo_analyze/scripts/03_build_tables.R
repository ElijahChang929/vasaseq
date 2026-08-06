#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_tables.R
#
# Builds every table in this folder. Draws nothing -- plot_all.R does that,
# reading only the TSVs written here and in vasaplate/.
#
# Reads per cell barcode for the PM26037 / ZHA9292A1 VASA library, grouped by
# the experimental design (4 genotypes x 4 replicates).
#
# It builds the sample sheet itself, from the two files that define it:
#
#   "vasa cbc 6bp.xlsx"     the ordered primer design (seq1..seq48). The 6 nt
#                           cell barcode sits immediately after the 6 N UMI:
#                           45 nt 5' handle + NNNNNN + 6 nt CBC + polyT.
#   bc_PM26037_6nt.tsv      the barcodes the pipeline actually detected in the
#                           data, and the cell id (001..016) it gave each one.
#
# Joining the two on the barcode sequence is what maps a design well to a
# pipeline cell id -- nothing here is hard-coded to a cell number.
#
# The genotype block assignment (GENOTYPE_ORDER below) is the one thing that
# comes from the user, not from a file: the 16 detected design wells, taken in
# ascending seq order, are 4 consecutive replicates of each genotype.
#
# Counts come from demux_read_counts.tsv (written by count_demux_reads.sh --
# an exact `wc -l / 4` over the step-1 per-cell FASTQs). The later pipeline
# stages are parsed out of the run's own step2/step3 reports and STAR logs, so
# no number in the output was typed by hand.
#
# Run with the r4.3 env:
#   /nemo/lab/turnerj/working/guangxin/envs/r4.3/bin/Rscript build_tables.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
  library(jsonlite)
})

# --- config ----------------------------------------------------------------

source(file.path(dirname(sub("^--file=", "",
       grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "paths.R"))
OWN <- OWNV                       # own_version/, where the xlsx and barcode file live
# Both overridable, so the same script serves the full-length run and the 75 nt
# truncation without a second copy:
#
#   DEMO_RUN  the pipeline output tree to read   (out/ or out75/)
#   DEMO_OUT  where the tables go                (demo_analyze/ or demo_analyze/own75/)
#
# The sample sheet, the barcode file and the design xlsx are shared -- truncating
# the biological read changes nothing about which barcode a read belongs to.
OUT  <- Sys.getenv("DEMO_OUT", unset = file.path(TAB, "own130"))
RUN  <- Sys.getenv("DEMO_RUN",
                   unset = "/nemo/lab/turnerj/working/guangxin/vasaseq/data/PM26037/out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

XLSX <- file.path(OWN, "vasa cbc 6bp.xlsx")
BCTSV <- file.path(OWN, "bc_PM26037_6nt.tsv")

# Layout of the ordered primer: 45 nt handle, 6 nt UMI, then the 6 nt barcode.
CBC_START <- 52L
CBC_END   <- 57L

# Design wells in ascending seq order -> genotype, 4 replicates each.
GENOTYPE_ORDER <- c("XY", "XO", "EpiLCs", "Blank control")
N_REPLICATE    <- 4L

# --- 1. sample sheet: design -> barcode -> pipeline cell id -----------------

design <- read_excel(XLSX) %>%
  transmute(
    seq_id  = str_remove(str_trim(Name), "^seq"),
    # some cells carry a UTF-8 BOM from Excel; strip it before slicing
    seqfull = str_remove(Sequence, "^﻿")
  ) %>%
  mutate(
    seq_no  = as.integer(seq_id),
    barcode = str_sub(seqfull, CBC_START, CBC_END)
  ) %>%
  arrange(seq_no)

stopifnot(all(str_detect(design$barcode, "^[ACGT]{6}$")))

detected <- read_tsv(BCTSV, col_names = c("barcode", "cell"),
                     col_types = "cc", trim_ws = TRUE) %>%
  # the barcode file writes cell ids 2-wide (01..16); every pipeline output
  # writes them 3-wide (001..016). Pad here so all later joins line up.
  mutate(cell = str_pad(cell, 3, pad = "0"))

sheet <- design %>%
  left_join(detected, by = "barcode")

# The design has one well the pipeline never saw; say so out loud rather than
# letting it disappear into a join.
missing <- sheet %>% filter(is.na(cell))
if (nrow(missing) > 0) {
  message("design wells with no detected barcode (dropped): ",
          paste0("seq", missing$seq_no, " (", missing$barcode, ")", collapse = ", "))
}
orphan <- setdiff(detected$barcode, design$barcode)
if (length(orphan) > 0) {
  stop("detected barcodes absent from the design: ", paste(orphan, collapse = ", "))
}

sheet <- sheet %>% filter(!is.na(cell)) %>% arrange(seq_no)

if (nrow(sheet) != length(GENOTYPE_ORDER) * N_REPLICATE) {
  stop(sprintf("expected %d detected wells, got %d -- the genotype blocks below assume %d x %d",
               length(GENOTYPE_ORDER) * N_REPLICATE, nrow(sheet),
               length(GENOTYPE_ORDER), N_REPLICATE))
}

sheet <- sheet %>%
  mutate(
    genotype  = factor(rep(GENOTYPE_ORDER, each = N_REPLICATE), levels = GENOTYPE_ORDER),
    replicate = rep(seq_len(N_REPLICATE), times = length(GENOTYPE_ORDER)),
    label     = sprintf("seq%d\ncell %s", seq_no, cell)
  ) %>%
  select(genotype, replicate, seq_no, barcode, cell, label)

# --- 2. read counts --------------------------------------------------------

# 2a. exact reads per barcode, out of step 1
demux_f <- file.path(OUT, "demux_read_counts.tsv")
if (!file.exists(demux_f)) {
  stop("missing ", demux_f, " -- run ./count_demux_reads.sh first")
}
demux <- read_tsv(demux_f, col_types = "cd") %>%
  transmute(cell = str_pad(cell, 3, pad = "0"), demultiplexed = reads)

# 2b. surviving each later stage, from the run's own reports
#     step2_report.txt: `in`, anchor, anc%, exact%, ->tg, kept, end%
#     step3_report.txt (TABLE 1): `in`, ribo, ribo%, kept, aln, mem, both, aln_len
parse_report <- function(path, ncol_expect, colnames_out, stop_at_table2 = FALSE) {
  ln <- read_lines(path)
  if (stop_at_table2) ln <- ln[seq_len(which(str_detect(ln, "^TABLE 2"))[1] - 1)]
  ln <- ln[str_detect(ln, "^\\s+\\d{3}\\s")]           # data rows only, not ALL
  parts <- str_split(str_trim(ln), "\\s+")
  stopifnot(all(lengths(parts) == ncol_expect))
  m <- do.call(rbind, parts)
  colnames(m) <- colnames_out
  as_tibble(m) %>%
    mutate(across(-cell, ~ as.numeric(str_remove_all(.x, "[,%]"))))
}

step2 <- parse_report(file.path(RUN, "logs/step2_report.txt"), 8,
                      c("cell", "in2", "anchor", "anc_pct", "exact_pct",
                        "to_tg", "trimmed", "end_pct"))

# Steps 3-7 can legitimately be absent: step 2 is re-run on its own whenever the
# trimming changes, and the downstream logs are cleared with it. Everything that
# depends only on step 2 must still build, so these degrade to NA columns with a
# message rather than aborting the whole script.
missing_downstream <- character()
step3_f <- file.path(RUN, "logs/step3_report.txt")
step3 <- if (file.exists(step3_f)) {
  parse_report(step3_f, 9,
               c("cell", "in3", "ribo", "ribo_pct", "non_rrna",
                 "aln", "mem", "both", "aln_len"), stop_at_table2 = TRUE)
} else {
  missing_downstream <- c(missing_downstream, "step3_report.txt (ribo, non_rrna)")
  tibble(cell = character(), ribo = numeric(), non_rrna = numeric())
}

# 2c. uniquely mapped, from each cell's STAR log (absent until step 4 has run)
if (length(Sys.glob(file.path(RUN, "cells", "*_E99_Log.final.txt"))) == 0) {
  missing_downstream <- c(missing_downstream, "STAR *_Log.final.txt (uniquely_mapped)")
}
star <- lapply(sheet$cell, function(cl) {
  f <- Sys.glob(file.path(RUN, "cells", sprintf("*_%s_*_E99_Log.final.txt", cl)))
  if (length(f) != 1) return(tibble(cell = cl, uniquely_mapped = NA_real_))
  ln <- read_lines(f)
  grab <- function(pat) as.numeric(str_trim(str_split_fixed(
    ln[str_detect(ln, pat)][1], "\\|", 2)[, 2]))
  tibble(cell = cl, uniquely_mapped = grab("Uniquely mapped reads number"))
}) %>% bind_rows()

# 2d. reads assigned to genes, as column sums of the step-7 read-count table
genes_f <- file.path(RUN, "ZHA9292A1_total.ReadCounts.tsv")
assigned <- if (file.exists(genes_f)) {
  # the table's header line starts with an empty field (the gene-name column),
  # so name it explicitly rather than letting the reader invent one
  hdr <- str_split(read_lines(genes_f, n_max = 1), "\t")[[1]]
  hdr[1] <- "gene"
  tab <- read_tsv(genes_f, skip = 1, col_names = hdr,
                  col_types = paste0("c", strrep("d", length(hdr) - 1L)))
  num <- tab[, -1, drop = FALSE]
  tibble(cell = names(num), assigned_to_genes = unname(colSums(num, na.rm = TRUE)))
} else {
  tibble(cell = character(), assigned_to_genes = numeric())
}

# 2e. UFIs per cell, as column sums of the step-7 UFI table.
#
# Restricted to rows carrying an Ensembl gene id, which drops tRNA and the few
# other non-gene rows (323 of 305,399). That is deliberate: it is the same
# definition the VASA-plate side uses -- vasaplate_check/vp_common.py sums
# species-assignable rows only (`species_of`), so `ufi_total` there is
# human-gene rows + mouse-gene rows. Matching it keeps the two datasets
# comparable in demo_analyze/vasaplate.
ufi_f <- file.path(RUN, "ZHA9292A1_total.UFICounts.tsv")
ufis <- if (file.exists(ufi_f)) {
  hdr <- str_split(read_lines(ufi_f, n_max = 1), "\t")[[1]]
  hdr[1] <- "gene"
  tab <- read_tsv(ufi_f, skip = 1, col_names = hdr,
                  col_types = paste0("c", strrep("d", length(hdr) - 1L)))
  tab <- tab[str_detect(tab$gene, "ENSMUSG"), , drop = FALSE]
  num <- tab[, -1, drop = FALSE]
  tibble(cell = names(num), ufi_total = unname(colSums(num, na.rm = TRUE)))
} else {
  tibble(cell = character(), ufi_total = numeric())
}

# 2f. step 2 stage by stage -- where reads are lost inside trimming.
#
# Reads step 2's own per-cell outputs, so it works whenever step 2 has run no
# matter what state steps 3-7 are in:
#   <cell>_cbc_bcanchor.log    pass 0: reads in, bases cut at the anchor
#   <cell>_cbc_cutadapt.json   pass 2: reads in (= TrimGalore's output) and out
#
# Pass 0 removes BASES, never whole reads -- it has no length filter at all -- so
# the read count is flat across it, and every read step 2 loses is lost either at
# pass 1 (TrimGalore --length 20 + Phred 20) or at pass 2 (cutadapt -m 20).
# That is why the per-stage figure shows pass 0 as no loss; it is a fact about
# the design, not a missing number.
trim_stages <- lapply(sheet$cell, function(cl) {
  bl <- Sys.glob(file.path(RUN, "cells", sprintf("*_%s_cbc_bcanchor.log", cl)))
  jf <- Sys.glob(file.path(RUN, "cells", sprintf("*_%s_cbc_cutadapt.json", cl)))
  if (length(bl) != 1 || length(jf) != 1) {
    return(tibble(cell = cl, after_pass0 = NA_real_, anchor_bases = NA_real_,
                  after_pass1 = NA_real_, after_pass2 = NA_real_))
  }
  t <- read_lines(bl)
  j <- fromJSON(jf)
  tibble(cell         = cl,
         after_pass0  = as.numeric(str_match(t[str_detect(t, "^reads ")], "^reads (\\d+)")[, 2]),
         anchor_bases = as.numeric(str_match(t[str_detect(t, "^bases removed")],
                                             "(\\d+)$")[, 2]),
         after_pass1  = j$read_counts$input,
         after_pass2  = j$read_counts$output)
}) %>% bind_rows()

if (any(is.na(trim_stages$after_pass2))) {
  missing_downstream <- c(missing_downstream,
                          "some *_cbc_cutadapt.json (trim stages) -- re-run step2")
}

per_cell <- sheet %>%
  left_join(trim_stages, by = "cell") %>%
  left_join(demux, by = "cell") %>%
  left_join(select(step2, cell, trimmed), by = "cell") %>%
  left_join(select(step3, cell, ribo, non_rrna), by = "cell") %>%
  left_join(star, by = "cell") %>%
  left_join(assigned, by = "cell") %>%
  left_join(ufis, by = "cell") %>%
  mutate(pct_of_library = 100 * demultiplexed / sum(demultiplexed))

# A failed join shows up as silent NAs and an empty-looking plot, so make the
# join itself the thing that fails.
if (any(is.na(per_cell$demultiplexed))) {
  stop("no read count for cell(s): ",
       paste(per_cell$cell[is.na(per_cell$demultiplexed)], collapse = ", "))
}

# Cross-check the exact count against the run's own bookkeeping. step2's `in`
# is read out of the barcode-anchor log; if the two disagree, something moved.
chk <- per_cell %>% left_join(select(step2, cell, in2), by = "cell") %>%
  mutate(d = abs(demultiplexed - in2) / demultiplexed)
if (any(chk$d > 0.01, na.rm = TRUE)) {
  message("NOTE: exact FASTQ count and step2 `in` differ by >1% for cell(s): ",
          paste(chk$cell[chk$d > 0.01], collapse = ", "))
}

per_genotype <- per_cell %>%
  group_by(genotype) %>%
  summarise(
    n_cells        = n(),
    total_reads    = sum(demultiplexed),
    mean_reads     = mean(demultiplexed),
    sd_reads       = sd(demultiplexed),
    sem_reads      = sd(demultiplexed) / sqrt(n()),
    median_reads   = median(demultiplexed),
    min_reads      = min(demultiplexed),
    max_reads      = max(demultiplexed),
    pct_of_library = 100 * sum(demultiplexed) / sum(per_cell$demultiplexed),
    mean_ufi       = mean(ufi_total),
    median_ufi     = median(ufi_total),
    .groups = "drop"
  )

# `label` carries a newline for the plot axis; flatten it so the TSVs stay
# one-record-per-line.
flat <- function(d) mutate(d, across(any_of("label"), ~ str_replace_all(.x, "\n", " ")))
write_tsv(flat(sheet),    file.path(OUT, "sample_sheet.tsv"))
write_tsv(flat(per_cell), file.path(OUT, "reads_per_cell.tsv"))
write_tsv(per_genotype,   file.path(OUT, "reads_per_genotype.tsv"))

# --- 3. long form: reads surviving each pipeline stage ----------------------
# Note this is NOT a monotone funnel: assigned_to_genes exceeds uniquely_mapped
# because the pipeline keeps multimappers (--outFilterMultimapNmax 20) and
# rescues them at assignment, so a read can reach a gene without mapping uniquely.

stages <- c(demultiplexed     = "Demultiplexed",
            trimmed           = "After trimming",
            non_rrna          = "After rRNA depletion",
            uniquely_mapped   = "Uniquely mapped",
            assigned_to_genes = "Assigned to genes")
present <- intersect(names(stages), names(per_cell))
present <- present[sapply(present, function(k) !all(is.na(per_cell[[k]])))]

stage_geno <- per_cell %>%
  select(genotype, all_of(present)) %>%
  pivot_longer(-genotype, names_to = "stage", values_to = "reads") %>%
  group_by(genotype, stage) %>%
  summarise(total_reads = sum(reads), .groups = "drop") %>%
  mutate(stage = factor(stage, levels = present, labels = unname(stages[present])))

write_tsv(stage_geno, file.path(OUT, "reads_by_stage_genotype.tsv"))

# --- 3b. trimming loss, long form ------------------------------------------
# One row per cell per fate. `lost_pass1` + `lost_pass2` + `kept` = demultiplexed
# by construction, so the stacked figure always totals the input.

TRIM_FATES <- c(kept       = "Kept",
                lost_pass2 = "Lost at pass 2 (cutadapt -m 20)",
                lost_pass1 = "Lost at pass 1 (TrimGalore --length 20)")

# No `label` column anywhere in these: the axis label carries a newline
# ("seq1\n002") and a newline inside a TSV field splits the record in two, which
# reads back as twice the rows with everything after it shifted. plot_all.R
# builds the label from seq_no and cell instead.
trim_cell <- per_cell %>%
  transmute(genotype, replicate, cell, seq_no,
            demultiplexed = after_pass0,
            lost_pass1  = after_pass0 - after_pass1,
            lost_pass2  = after_pass1 - after_pass2,
            kept        = after_pass2,
            pct_lost    = 100 * (after_pass0 - after_pass2) / after_pass0,
            pct_kept    = 100 * after_pass2 / after_pass0,
            anchor_bases)

stopifnot(all(with(trim_cell,
                   abs(kept + lost_pass1 + lost_pass2 - demultiplexed) < 1e-6)))

trim_long <- trim_cell %>%
  select(genotype, replicate, cell, seq_no, demultiplexed,
         all_of(names(TRIM_FATES))) %>%
  pivot_longer(all_of(names(TRIM_FATES)), names_to = "fate", values_to = "reads") %>%
  mutate(fate = factor(unname(TRIM_FATES[fate]), levels = unname(TRIM_FATES)),
         pct  = 100 * reads / demultiplexed)

trim_geno <- trim_cell %>%
  group_by(genotype) %>%
  summarise(n_cells       = n(),
            demultiplexed = sum(demultiplexed),
            lost_pass1    = sum(lost_pass1),
            lost_pass2    = sum(lost_pass2),
            kept          = sum(kept),
            .groups = "drop") %>%
  mutate(lost_total   = lost_pass1 + lost_pass2,
         pct_pass1    = 100 * lost_pass1 / demultiplexed,
         pct_pass2    = 100 * lost_pass2 / demultiplexed,
         pct_lost     = 100 * lost_total / demultiplexed,
         pct_kept     = 100 * kept / demultiplexed)

trim_geno_long <- trim_geno %>%
  select(genotype, demultiplexed, all_of(names(TRIM_FATES))) %>%
  pivot_longer(all_of(names(TRIM_FATES)), names_to = "fate", values_to = "reads") %>%
  mutate(fate = factor(unname(TRIM_FATES[fate]), levels = unname(TRIM_FATES)),
         pct  = 100 * reads / demultiplexed)

# --- 3c. read classes -------------------------------------------------------
# From own_version/diagnostics/classify_reads.py: every read in exactly one of five boxes, crossing
# "did it read through into its barcode" with "where did it end up". Present
# only once that has been run over the current step 2 output.
cls_f <- file.path(RUN, "read_classes.tsv")
if (file.exists(cls_f)) {
  cls <- read_tsv(cls_f, col_types = "ccdd") %>%
    mutate(cell = str_pad(cell, 3, pad = "0")) %>%
    inner_join(select(sheet, genotype, replicate, cell, seq_no), by = "cell")

  # classes 1 and 5 are the two kept boxes; they must reproduce step 2's `kept`
  chk_kept <- cls %>% filter(str_starts(class, "1|5")) %>% summarise(s = sum(reads)) %>% pull(s)
  if (abs(chk_kept - sum(per_cell$after_pass2)) > 0) {
    warning("read classes 1+5 (", comma(chk_kept), ") disagree with step 2 kept (",
            comma(sum(per_cell$after_pass2)), ")")
  }

  cls_geno <- cls %>%
    group_by(genotype, class) %>%
    summarise(reads = sum(reads), .groups = "drop_last") %>%
    mutate(pct = 100 * reads / sum(reads)) %>%
    ungroup()

  write_tsv(cls,      file.path(OUT, "read_classes_per_cell.tsv"))
  write_tsv(cls_geno, file.path(OUT, "read_classes_per_genotype.tsv"))
} else {
  missing_downstream <- c(missing_downstream,
                          "read_classes.tsv -- run own_version/diagnostics/classify_reads.py")
}

# --- 3d. what cut the reads pass 2 dropped ----------------------------------
# From own_version/diagnostics/attribute_pass2_loss.py: for every read in class 4, which adapter
# removed the most bases from it. Answers "cutadapt -m dropped it" with a cause.
att_f <- file.path(RUN, "pass2_adapter_attribution.tsv")
if (file.exists(att_f)) {
  att <- read_tsv(att_f, col_types = "ccddd") %>%
    mutate(cell = str_pad(cell, 3, pad = "0")) %>%
    inner_join(select(sheet, genotype, replicate, cell, seq_no), by = "cell")

  att_geno <- att %>%
    group_by(genotype, adapter) %>%
    summarise(reads = sum(reads), bases_removed = sum(bases_removed),
              .groups = "drop_last") %>%
    mutate(pct = 100 * reads / sum(reads)) %>%
    ungroup()

  write_tsv(att,      file.path(OUT, "pass2_attribution_per_cell.tsv"))
  write_tsv(att_geno, file.path(OUT, "pass2_attribution_per_genotype.tsv"))
} else {
  missing_downstream <- c(missing_downstream,
                          "pass2_adapter_attribution.tsv -- run own_version/diagnostics/attribute_pass2_loss.py")
}

# --- 3e. length of what survives step 2 -------------------------------------
# From own_version/diagnostics/read_length_dist.sh. Normalised within genotype, because the depths
# differ ~20x and the question is the SHAPE of the distribution, not its height.
# The floor that decides "too short". Read from the environment so the table
# cannot drift from what the run actually used (config.sh sets it to 15).
TRIM_MINLEN <- as.integer(Sys.getenv("TRIM_MINLEN", unset = "15"))

# Prefer the PRE-filter tally: it measures every read at its final trimmed
# length, including the ones -m discards. The post-filter file has had its whole
# left tail removed by the length floor, so the part that matters -- how much of
# the library ends up too short, and by how much -- is simply absent from it.
len_f <- file.path(RUN, "read_length_prefilter.tsv")
prefilter <- file.exists(len_f)
if (!prefilter) len_f <- file.path(RUN, "read_length_dist.tsv")
if (file.exists(len_f)) {
  if (!prefilter) {
    message("read_length: only the post-filter tally is present; the <", TRIM_MINLEN,
            " nt tail is missing. Re-run own_version/diagnostics/read_length_dist.sh --prefilter")
  }
  len <- read_tsv(len_f, col_types = "cdd") %>%
    mutate(cell = str_pad(cell, 3, pad = "0")) %>%
    inner_join(select(sheet, genotype, replicate, cell), by = "cell")

  len_geno <- len %>%
    group_by(genotype, length) %>%
    summarise(reads = sum(reads), .groups = "drop_last") %>%
    mutate(pct = 100 * reads / sum(reads)) %>%
    ungroup()

  # the quartiles are what the figure cannot show precisely
  len_stats <- len %>%
    group_by(genotype) %>%
    summarise(n = sum(reads),
              median = { c <- cumsum(reads[order(length)]); l <- sort(length)
                         l[which(c >= sum(reads) / 2)[1]] },
              q1 = { c <- cumsum(reads[order(length)]); l <- sort(length)
                     l[which(c >= sum(reads) * 0.25)[1]] },
              q3 = { c <- cumsum(reads[order(length)]); l <- sort(length)
                     l[which(c >= sum(reads) * 0.75)[1]] },
              pct_under_30 = 100 * sum(reads[length < 30]) / sum(reads),
              pct_discarded = 100 * sum(reads[length < TRIM_MINLEN]) / sum(reads),
              .groups = "drop") %>%
    mutate(prefilter = prefilter, minlen = TRIM_MINLEN)

  write_tsv(len_geno,  file.path(OUT, "read_length_per_genotype.tsv"))
  write_tsv(len_stats, file.path(OUT, "read_length_stats.tsv"))

  cat("\n== length once trimming can remove nothing more ==",
      if (prefilter) "" else "  (POST-filter: left tail missing)", "\n", sep = "")
  print(as.data.frame(len_stats %>%
    transmute(genotype, reads = comma(n), Q1 = q1, median, Q3 = q3,
              `<30nt` = sprintf("%.1f%%", pct_under_30),
              `discarded <15nt` = sprintf("%.1f%%", pct_discarded))), row.names = FALSE)
} else {
  missing_downstream <- c(missing_downstream,
                          "read_length_dist.tsv -- run own_version/diagnostics/read_length_dist.sh")
}

write_tsv(trim_cell,      file.path(OUT, "trim_loss_per_cell.tsv"))
write_tsv(trim_long,      file.path(OUT, "trim_loss_per_cell_long.tsv"))
write_tsv(trim_geno,      file.path(OUT, "trim_loss_per_genotype.tsv"))
write_tsv(trim_geno_long, file.path(OUT, "trim_loss_per_genotype_long.tsv"))

# --- 4. console summary ----------------------------------------------------

cat("\n== sample sheet ==\n")
print(as.data.frame(sheet), row.names = FALSE)
cat("\n== reads per cell ==\n")
print(as.data.frame(per_cell %>%
  transmute(genotype, replicate, seq = paste0("seq", seq_no), cell, barcode,
            reads = comma(demultiplexed), pct = sprintf("%.2f%%", pct_of_library))),
  row.names = FALSE)
cat("\n== reads per genotype ==\n")
print(as.data.frame(per_genotype %>%
  transmute(genotype, n_cells,
            total = comma(total_reads), mean = comma(round(mean_reads)),
            sd = comma(round(sd_reads)), min = comma(min_reads), max = comma(max_reads),
            pct = sprintf("%.2f%%", pct_of_library))),
  row.names = FALSE)
cat("\n== trimming loss per genotype ==\n")
print(as.data.frame(trim_geno %>%
  transmute(genotype, n_cells,
            demultiplexed = comma(demultiplexed),
            lost_pass1 = sprintf("%s (%.1f%%)", comma(lost_pass1), pct_pass1),
            lost_pass2 = sprintf("%s (%.1f%%)", comma(lost_pass2), pct_pass2),
            kept       = sprintf("%s (%.1f%%)", comma(kept), pct_kept))),
  row.names = FALSE)

if (length(missing_downstream)) {
  cat("\nNOTE: built without these (steps 3-7 not re-run yet) --",
      "affected columns are NA:\n", sep = " ")
  for (m in missing_downstream) cat("  -", m, "\n")
}

cat("\nwrote: sample_sheet.tsv, reads_per_cell.tsv, reads_per_genotype.tsv,\n",
    "       reads_by_stage_genotype.tsv,\n",
    "       trim_loss_per_cell.tsv, trim_loss_per_cell_long.tsv,\n",
    "       trim_loss_per_genotype.tsv, trim_loss_per_genotype_long.tsv\n", sep = "")
