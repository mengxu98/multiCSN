# Plot pseudotime states

Plot pseudotime states

## Usage

``` r
plot_states(
  object,
  pseudotime_column = NULL,
  group_column = NULL,
  palette = NULL,
  palette_name = "Chinese"
)
```

## Arguments

- object:

  A `Seurat` object.

- pseudotime_column:

  Optional pseudotime column. If missing, the function reuses the value
  stored in the active dynamic network or in
  `DensityPoints_<pseudotime_column>`.

- group_column:

  Optional grouping column for density curves.

- palette:

  Optional colors.

- palette_name:

  Palette name used by `.multicsn_palette_colors()`.

## Value

A data.frame of kept state windows.
