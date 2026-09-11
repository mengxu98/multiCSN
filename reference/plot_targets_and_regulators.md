# Updated plot of top regulators given targets in dynamic networks based on a weight column. Top regulators computed for each epoch, but maintained in plot across epochs if present in epoch subnetwork.

Updated plot of top regulators given targets in dynamic networks based
on a weight column. Top regulators computed for each epoch, but
maintained in plot across epochs if present in epoch subnetwork.

## Usage

``` r
plot_targets_and_regulators(
  network,
  targets,
  epochs = NULL,
  weight_column = "weight",
  gene_ranks = NULL,
  regulators_num = 5,
  network_order = NULL,
  fixed_layout = TRUE,
  declutter = TRUE,
  show_expression = TRUE,
  node_size = 5,
  label_size = 2,
  title_size = 10,
  legend_size = 8
)
```

## Arguments

- network:

  the result of running epochGRN

- targets:

  targets

- epochs:

  result of running assign_epochs

- weight_column:

  column name containing reconstruction weights to use

- gene_ranks:

  gene_ranks

- regulators_num:

  number of top regulators to plot

- network_order:

  which epochs or transitions to plot

- fixed_layout:

  whether or not to fix node positions across epoch networks

- declutter:

  if TRUE, will only label nodes with active interactions in given
  network

- show_expression:

  if TRUE, size and shade of node indicates mean expression in a given
  epoch.

- node_size:

  node_size

- label_size:

  label_size

- title_size:

  title_size

- legend_size:

  legend_size

## Value

plot
