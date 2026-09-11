# Get network graph

Get network graph

Get network graph

## Usage

``` r
NetworkGraph(object, ...)

NetworkGraph(object, ...)

# S4 method for class 'Seurat'
NetworkGraph(
  object,
  network = DefaultNetwork(object),
  graph = "module_graph",
  celltypes = NULL,
  ...
)

# S4 method for class 'Network'
NetworkGraph(object, graph = "module_graph", ...)
```

## Arguments

- object:

  The input data, a csn object.

- ...:

  Arguments for other methods

- network:

  network

- graph:

  graph

- celltypes:

  cell types to analyze, NULL for all cell types
