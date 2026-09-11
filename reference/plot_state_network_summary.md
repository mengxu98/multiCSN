# Plot state-level network summary metrics

Plot state-level network summary metrics

## Usage

``` r
plot_state_network_summary(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  celltypes = NULL,
  metrics = c("n_cells", "n_edges", "n_regulators", "n_targets"),
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

- metrics:

  Metrics to display.

- weight_cutoff:

  Optional absolute-weight cutoff for exported edges.

## Value

A `ggplot` object.
