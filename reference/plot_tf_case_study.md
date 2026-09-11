# Plot single-TF case study across ordered states

Plot single-TF case study across ordered states

## Usage

``` r
plot_tf_case_study(
  object,
  tf,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  method = c("page_rank", "degree_distribution"),
  weight_column = "weight",
  top_targets = 15,
  ranked_features = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- tf:

  A transcription factor name.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- pseudotime_column:

  Optional pseudotime column.

- assay:

  Assay used to fetch expression.

- method:

  Gene-rank method.

- weight_column:

  Weight column used when exporting networks.

- top_targets:

  Number of top targets to display in the target heatmap.

- ranked_features:

  Optional precomputed regulator ranking table; when `NULL`,
  `rank_features(type = "regulators")` is called.

## Value

A patchwork/ggplot object.
