# Infer dynamic state networks from density-partitioned pseudotime states

Infer dynamic state networks from density-partitioned pseudotime states

## Usage

``` r
state_dynamic(
  object,
  pseudotime_column,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  penalty = "L0",
  r_squared_threshold = 0,
  regulators = NULL,
  targets = NULL,
  cores = 1,
  verbose = TRUE,
  dynamic_features = NULL,
  group_column = NULL,
  posterior_threshold = 0.8,
  n_grid = 512,
  n_boot = 50,
  bandwidth = "nrd0",
  overwrite_density_points = FALSE,
  min_cells = NULL,
  ...
)
```

## Arguments

- object:

  A `Seurat` object carrying CSNObject state.

- pseudotime_column:

  Metadata column containing pseudotime values.

- method:

  Network inference method.

- penalty:

  Penalty passed to
  [`inferCSN()`](https://mengxu98.github.io/multiCSN/reference/inferCSN.md)
  for `Network`.

- r_squared_threshold:

  R-squared threshold.

- regulators:

  Optional regulator set.

- targets:

  Optional target set.

- cores:

  Number of cores.

- verbose:

  Logical.

- dynamic_features:

  Optional dynamic feature filtering config.

- group_column:

  Optional grouping column reused by
  [`density_points()`](https://mengxu98.github.io/multiCSN/reference/density_points.md).

- posterior_threshold:

  Maximum pairwise posterior dominance allowed inside a transition
  window; lower values make transition windows narrower.

- n_grid:

  Number of grid points used for kernel density evaluation.

- n_boot:

  Number of bootstrap resamples used to assess boundary stability.

- bandwidth:

  Bandwidth passed to
  [`stats::density()`](https://rdrr.io/r/stats/density.html).

- overwrite_density_points:

  Logical; when `TRUE`, recompute and overwrite a cached density
  partition.

- min_cells:

  Ignored. Retained for backwards compatibility; transition states are
  determined automatically.

- ...:

  Additional parameters forwarded to
  [`inferCSN()`](https://mengxu98.github.io/multiCSN/reference/inferCSN.md)
  on `Network` objects.

## Value

A `Seurat` object with a dynamic network stored in CSNObject state.
