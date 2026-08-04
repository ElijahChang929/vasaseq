# Experimental design — 384-well plate layout

The plate figure for the user's own VASA library (`PM26037` / `ZHA9292A1`).

## The design being shown

The library uses **16 cell barcodes** (`../I_Gene_expression/own_version/bc_PM26037_6nt.tsv`,
6 nt each, numbered 01–16). A 384-well plate is 16 rows × 24 columns, so the
barcodes map onto the plate **one per row**: every well in a row would carry the
same barcode, and the 24 columns of a row are replicates of it.

Only **column 1** was actually loaded — 16 wells, A1–P1, one cell per barcode.
The figure therefore fills column 1 and leaves columns 2–24 as empty wells, so
the barcode-per-row scheme and what was used are both visible in one picture.

## Files

| file | what it is |
|---|---|
| `plate_design.R` | the whole thing: builds the layout, writes the table, draws the plate |
| `plate_layout_384.png` / `.pdf` | the figure (8 × 5.5 in; PDF is the one to put in a paper/slide) |
| `plate_layout_384.tsv` | the layout as a table — `position`, `row`, `col`, `barcode`, `barcode_seq`. `barcode` is `NA` for columns 2–24, matching the figure |

## Environment

Drawing uses [`ggplate`](https://cran.r-project.org/package=ggplate) 0.3.1,
which is not in any of the pre-existing R envs here. A dedicated conda env was
built for it (R 4.4 + ggplate + ggplot2 + dplyr, all from conda-forge):

```bash
conda create -y -p /nemo/lab/turnerj/working/guangxin/envs/ggplate \
  -c conda-forge r-base=4.4 r-ggplate=0.3.1 r-ggplot2 r-dplyr
```

## Running it

```bash
cd "code/experimental design"
/nemo/lab/turnerj/working/guangxin/envs/ggplate/bin/Rscript plate_design.R
```

Re-runnable; it overwrites its three outputs and writes nothing else.

## Notes on the code (things that bite)

- **The value column must be `character`, not `factor`.** `plate_plot()` calls
  `min()`/`max()` on it, which errors on a factor (`'min' not meaningful for
  factors`). The barcode labels are zero-padded (`barcode 01`…`barcode 16`) so
  the alphabetical order ggplate uses is the numeric order.
- **Colours are `scales::hue_pal()(16)`** — ggplot2's own default discrete
  palette, expanded to 16. ggplate's built-in discrete palettes do not reach 16
  levels and it errors if you have more levels than colours.
- Empty wells come from setting `barcode` to `NA`; ggplate's `remove_na = TRUE`
  (the default) draws those as outlines rather than filling them with
  `na_fill`. To drop them entirely, filter the data frame before plotting.
- The barcode sequences are inlined in the script rather than read from
  `bc_PM26037_6nt.tsv`, so the figure does not depend on a path outside this
  folder. **If that barcode file ever changes, update the vector here too.**

## Changing the layout

Everything is driven by `plate$barcode_n` (which barcode a well gets) and the
`plate$barcode[plate$col != 1] <- NA` line (which wells are used). For barcodes
cycling *along* each row instead of down the rows, that is a one-line change to
`barcode_n`.
