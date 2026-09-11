# Plot cell type rewiring heatmap

Plot cell type rewiring heatmap

## Usage

``` r
plot_celltype_specificity_heatmap(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  top_tfs = 20,
  metric = c("rewiring_score", "regulon_jaccard"),
  weight_column = "weight"
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Static network name. Defaults to `DefaultNetwork(object)`.

- celltypes:

  Optional subset of cell types.

- top_tfs:

  Number of TFs to display.

- metric:

  Heatmap value. Supported values include `"rewiring_score"`,
  `"weighted_regulon_jaccard"`, `"binary_rewiring_score"`, and
  `"regulon_jaccard"`.

- weight_column:

  Weight column used when exporting networks.

## Value

A `ggplot` object.
