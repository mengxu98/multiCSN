# Plot manuscript Figure 3 key-TF summary

Plot manuscript Figure 3 key-TF summary

## Usage

``` r
plot_key_tf_summary_figure(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  tfs = NULL,
  method = c("page_rank", "degree_distribution"),
  weight_column = "weight",
  top_tfs_heatmap = 25,
  n_case_tfs = 6,
  feature_scores = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- pseudotime_column:

  Optional pseudotime column.

- assay:

  Assay used to fetch expression.

- tfs:

  Optional TF vector. If `NULL`, select automatically.

- method:

  Gene-rank method.

- weight_column:

  Weight column used when exporting networks.

- top_tfs_heatmap:

  Number of TFs shown in the global TF-dynamics heatmap.

- n_case_tfs:

  Number of TFs selected automatically when `tfs=NULL`.

- feature_scores:

  Optional integrated key-TF score table; when `NULL`, recomputed with
  `prioritize_features(strategy = "integrated")`.

## Value

A patchwork/ggplot object.
