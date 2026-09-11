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
  scale_value = c("zscore", "raw")
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

## Value

A `ggplot` object.
