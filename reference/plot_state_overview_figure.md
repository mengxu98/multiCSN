# Plot paper-style ordered-state overview

Plot paper-style ordered-state overview

## Usage

``` r
plot_state_overview_figure(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  state_ids = NULL,
  top_tfs_network = 10,
  top_targets_per_tf = 12,
  max_targets = 60,
  top_tfs_bar = 8,
  palette = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- pseudotime_column:

  Optional pseudotime column.

- assay:

  Assay used to summarize state-level expression in subnetworks.

- state_ids:

  Optional subset of states. Defaults to kept inferable states.

- top_tfs_network:

  Number of TFs shown in each state subnetwork.

- top_targets_per_tf:

  Number of top targets retained per TF in each subnetwork.

- max_targets:

  Maximum union of targets shown per state subnetwork.

- top_tfs_bar:

  Number of TFs shown in the bottom target-count summary.

- palette:

  Optional named palette for the density panel.

## Value

A patchwork/ggplot object.
