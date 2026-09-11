# Plot manuscript Figure 2 rewiring overview

Plot manuscript Figure 2 rewiring overview

## Usage

``` r
plot_rewiring_overview_figure(
  object,
  network = DefaultNetwork(object),
  state_ids = NULL,
  top_tfs = 20,
  metric = c("rewiring_score", "regulon_jaccard"),
  weight_column = "weight"
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- state_ids:

  Optional subset of states.

- top_tfs:

  Number of rewired TFs kept per transition in the heatmap.

- metric:

  Rewiring metric used for the heatmap, `"rewiring_score"` or
  `"regulon_jaccard"`.

- weight_column:

  Weight column used when exporting networks.

## Value

A patchwork/ggplot object.
