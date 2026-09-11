# Plot network graph

Plot network graph

Plot network graph

Plot network graph.

## Usage

``` r
plot_network_graph(object, ...)

plot_network_graph(object, ...)

# S4 method for class 'Seurat'
plot_network_graph(
  object,
  network = NULL,
  celltypes = NULL,
  graph = "module_graph",
  layout = "umap",
  edge_width = 0.2,
  edge_color = c(`-1` = "darkgrey", `1` = "orange"),
  node_color = pals::magma(100),
  node_size = c(1, 5),
  text_size = 10,
  color_nodes = TRUE,
  label_nodes = TRUE,
  color_edges = TRUE,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...
)
```

## Arguments

- object:

  CSNObject or Network

- ...:

  Additional arguments

- network:

  Name of the network to use.

- celltypes:

  Celltypes to plot. If `NULL`, all celltypes are plotted.

- graph:

  Name of the graph.

- layout:

  Layout for the graph. Can be 'umap' or any force-directed layout
  implemented in
  [`ggraph`](https://ggraph.data-imaginist.com/reference/ggraph.html)

- edge_width:

  Edge width.

- edge_color:

  Edge color.

- node_color:

  Node color or color gradient.

- node_size:

  Node size range.

- text_size:

  Font size for labels.

- color_nodes:

  Logical, Whether to color nodes by centrality.

- label_nodes:

  Logical, Whether to label nodes with gene name.

- color_edges:

  Logical, Whether to color edges by direction.

- combine:

  Logical. If `TRUE` (default), return a combined patchwork plot; if
  `FALSE`, return a named list of individual plots.

- nrow, ncol:

  Integer or `NULL`. Number of rows/columns for
  [`patchwork::wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html).

- byrow:

  Logical. If `TRUE`, fill plots by row; passed to `wrap_plots`.

## Value

A combined ggplot (combine=TRUE) or list of ggplots (combine=FALSE).
