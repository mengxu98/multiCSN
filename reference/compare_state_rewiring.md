# Compare adjacent state rewiring

Compare adjacent state rewiring

## Usage

``` r
compare_state_rewiring(
  object,
  network = DefaultNetwork(object),
  state_pairs = NULL,
  weight_column = "weight"
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- state_pairs:

  Optional two-column data.frame with columns `state_from` and
  `state_to`. If `NULL`, adjacent ordered states are compared.

- weight_column:

  Weight column used when exporting networks.

## Value

A list with pair summaries and per-TF rewiring tables.
