# Summarize transition-specific TF rewiring and target gain/loss

Summarize transition-specific TF rewiring and target gain/loss

## Usage

``` r
summarize_transition_analysis(
  object,
  network = DefaultNetwork(object),
  state_pairs = NULL,
  weight_column = "weight",
  top_tfs = 10,
  top_targets = 10
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

- top_tfs:

  Number of rewired TFs kept per transition.

- top_targets:

  Number of gain/loss targets kept per TF and transition.

## Value

A list with transition summary, TF summary, and target summary tables.
