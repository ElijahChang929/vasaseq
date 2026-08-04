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

The 16 loaded wells are four samples, four barcodes each:

| wells | barcodes | sample |
|---|---|---|
| A1–D1 | 01–04 | XY |
| E1–H1 | 05–08 | XO |
| I1–L1 | 09–12 | EpiLCs |
| M1–P1 | 13–16 | NTC (no-template control) |

Worth knowing when reading the mapping output: `../I_Gene_expression/own_version/README.md`
reports **cells 001, 014, 015 and 016** as 5–30× lower than the rest and
behaving as blanks. Three of those (014–016) are NTC wells, as expected — but
**001 is an XY sample**, not a control.

## Files

| file | what it is |
|---|---|
| `plate_design.R` | the whole thing: builds the layout, writes the table, draws the plate |
| `plate_layout_384.png` / `.pdf` | the figure (8 × 5.5 in; PDF is the one to put in a paper/slide) |
| `plate_layout_384.tsv` | the layout as a table — `position`, `row`, `col`, `sample`, `barcode`, `barcode_seq`. `sample` and `barcode` are `NA` for columns 2–24, matching the figure |

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

The sample brackets come from the `sample_groups` data frame — edit
`first_bc`/`last_bc` there and both the figure and the TSV follow. They are
drawn as plain `geom_segment` + `geom_text` on top of the ggplate object:
`plate_plot()`'s panel is an ordinary continuous grid where **x is the column
number (1–24) and y counts rows from the bottom**, so row A is `y = 16` and row
P is `y = 1`.
