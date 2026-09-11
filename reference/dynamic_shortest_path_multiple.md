# Function to return shortest path from multiple TFs to multiple targets in a dynamic network

Function to return shortest path from multiple TFs to multiple targets
in a dynamic network

## Usage

``` r
dynamic_shortest_path_multiple(object, ...)

# S4 method for class 'list'
dynamic_shortest_path_multiple(
  object,
  regulators,
  targets,
  weight_column = "weight"
)

# S4 method for class 'Network'
dynamic_shortest_path_multiple(
  object,
  regulators,
  targets,
  weight_column = "weight"
)

# S4 method for class 'Seurat'
dynamic_shortest_path_multiple(
  object,
  network = NULL,
  celltypes = NULL,
  regulators,
  targets,
  weight_column = "weight",
  weight_cutoff = NULL,
  corr_column = "mean_corr"
)

# S4 method for class 'CSNObject'
dynamic_shortest_path_multiple(object, ...)
```

## Arguments

- object:

  A CSNObject with a dynamic network (states in
  object@networks\[\[network\]\]).

- ...:

  Arguments for other methods.

- regulators:

  Starting genes/TFs.

- targets:

  End genes.

- weight_column:

  Column name for weights.

- network:

  Name of the dynamic network.

- celltypes:

  States/cell types to include; NULL for all.

- weight_cutoff:

  Optional cutoff passed to export_csn.

- corr_column:

  Column name for correlation; if absent, fall back to weight.

## Value

dataframe with shortest path, distance, normalized distance, and action
