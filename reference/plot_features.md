# Plot top network features per state

Plot top network features per state

## Usage

``` r
plot_features(
  object = NULL,
  data = NULL,
  type = c("regulators", "targets"),
  network = NULL,
  celltypes = NULL,
  method = c("page_rank", "degree_distribution"),
  top_n = 10,
  value_column = NULL,
  weight_column = "weight",
  view = c("ranking", "dynamics"),
  scale_value = c("zscore", "raw"),
  x_text_angle = NULL,
  cell_size = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  return_grob = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- data:

  Optional precomputed ranked-feature table used instead of `object`.

- type:

  Feature type, `"regulators"` or `"targets"`.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- celltypes:

  Optional subset of states.

- method:

  Gene-rank method used when `type = "regulators"`.

- top_n:

  Number of features to display per state.

- value_column:

  Optional score column to plot.

- weight_column:

  Weight column used when `type = "targets"`.

- view:

  View mode, `"ranking"` for per-state rankings or `"dynamics"` for
  per-state feature-score heatmaps (`type = "regulators"` only).

- scale_value:

  Score scaling for `view = "dynamics"`: `"zscore"` to standardize each
  feature across states, or `"raw"` for raw scores.

- x_text_angle:

  Column-label angle for dynamics heatmaps. If `NULL`, choose
  automatically from the number and length of state labels.

- cell_size:

  Side length of each square heatmap cell in mm or an absolute grid
  unit. If `NULL`, choose automatically from the matrix dimensions.

- cluster_rows:

  Logical, whether to cluster rows in dynamics heatmap. Default is
  `FALSE`.

- cluster_columns:

  Logical, whether to cluster columns in dynamics heatmap. Default is
  `FALSE`.

- return_grob:

  Logical, whether to return a grid `grob` instead of a
  `ggplot`/`patchwork` object. Default is `FALSE`.

## Value

A `ggplot` or `patchwork` object (or a grid `grob` if
`return_grob = TRUE`).
