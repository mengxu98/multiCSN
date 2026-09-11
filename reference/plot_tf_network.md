# Plot TF network

Plot TF network

Plot TF network

Plot sub-network centered around one TF.

## Usage

``` r
plot_tf_network(object, ...)

plot_tf_network(object, ...)

# S4 method for class 'Seurat'
plot_tf_network(
  object,
  tfs = NULL,
  network = NULL,
  celltypes = NULL,
  graph = "module_graph",
  circular = TRUE,
  edge_width = 0.2,
  edge_color = c(`-1` = "darkgrey", `1` = "orange"),
  node_size = 3,
  text_size = 10,
  label_nodes = c("tfs", "all", "none"),
  color_edges = TRUE,
  ...
)
```

## Arguments

- object:

  CSNObject or Network

- ...:

  Additional arguments

- tfs:

  The transcription factor to center around.

- network:

  Name of the network to use.

- celltypes:

  Celltypes to plot. If `NULL`, all celltypes are plotted.

- graph:

  Name of the graph.

- circular:

  Logical. Layout tree in circular layout.

- edge_width:

  Edge width.

- edge_color:

  Edge color.

- node_size:

  Node size.

- text_size:

  Font size for labels.

- label_nodes:

  String, indicating what to label. \* `'tfs'` - Label all TFs. \*
  `'all'` - Label all genes. \* `'none'` - Label nothing (except the
  root TF).

- color_edges:

  Logical, whether to color edges by direction.

## Value

A CSNObject object.
