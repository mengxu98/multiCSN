# Get network graph

Get network graph

Get network graph

Compute network graph embedding using UMAP.

## Usage

``` r
get_network_graph(object, ...)

get_network_graph(object, ...)

# S4 method for class 'Seurat'
get_network_graph(
  object,
  network = NULL,
  celltypes = NULL,
  graph_name = "module_graph",
  rna_assay = "RNA",
  rna_layer = "data",
  umap_method = c("weighted", "coef", "corr", "none"),
  features = NULL,
  seed = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  CSNObject or Network

- ...:

  Additional arguments for
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html).

- network:

  Name of the network to use.

- celltypes:

  Celltypes to plot. If `NULL`, all celltypes are plotted.

- graph_name:

  Name of the graph.

- rna_assay:

  Name of the RNA assay.

- rna_layer:

  Name of the RNA slot to use.

- umap_method:

  Method to compute edge weights for UMAP: \* `'weighted'` - Correlation
  weighted by GRN coefficient. \* `'corr'` - Only correlation. \*
  `'coef'` - Only GRN coefficient. \* `'none'` - Don't compute UMAP and
  create the graph directly from modules.

- features:

  Features to use to create the graph. If `NULL` uses all features in
  the network.

- seed:

  Random seed for UMAP computation

- verbose:

  Print messages.

## Value

A CSNObject object.
