# Summarize per-state network architecture

Summarize per-state network architecture

## Usage

``` r
summarize_state_networks(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  celltypes = NULL,
  weight_cutoff = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- pseudotime_column:

  Optional pseudotime column.

- celltypes:

  Optional subset of states.

- weight_cutoff:

  Optional absolute-weight cutoff for exported edges.

## Value

A data.frame summarizing cells, edges, regulators, and targets per
state.
