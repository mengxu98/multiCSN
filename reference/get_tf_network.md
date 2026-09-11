# Get TF network

Get TF network

Get TF network

Get sub-network centered around one TF.

## Usage

``` r
get_tf_network(object, ...)

get_tf_network(object, ...)

# S4 method for class 'Seurat'
get_tf_network(
  object,
  celltypes = NULL,
  tfs = NULL,
  features = NULL,
  network = NULL,
  graph = "module_graph",
  order = 3,
  keep_all_edges = FALSE,
  verbose = TRUE,
  cores = 1,
  ...
)
```

## Arguments

- object:

  CSNObject or Network

- ...:

  Additional arguments

- celltypes:

  Celltypes to plot. If `NULL`, all celltypes are plotted.

- tfs:

  The transcription factors to center around.

- features:

  Features to use. If `NULL` uses all features in the graph.

- network:

  Name of the network to use.

- graph:

  Name of the graph.

- order:

  Integer indicating the maximal order of the graph.

- keep_all_edges:

  Logical, whether to maintain all edges to each leaf or prune to the
  strongest overall connection.

- verbose:

  Logical. Whether to print messages.

- cores:

  Logical. Whether to parallelize the computation with
  [`foreach`](https://rdrr.io/pkg/foreach/man/foreach.html).

## Value

A CSNObject object.
